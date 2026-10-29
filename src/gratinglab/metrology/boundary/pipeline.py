"""
Building a boundary profile, separated from writing one.

The computation used to live inside ``workflows.run_boundary_profile_export``,
fused with the file writing. Separating the two lets the GUI preview a profile,
and recompute it as controls move, without touching the disk, while the CLI
mode writes exactly what the panel showed.

The flattening here differs from the profile flattening used for blaze-angle
work. ``flatten_endpoints`` forces the ends of the trace to zero so one groove
tiles seamlessly into the next; profile flattening removes a fitted background
so a facet can be measured. The two serve different purposes, and substituting
one for the other would produce a profile with a step at the period boundary.
That defect makes a PCGrate efficiency curve wrong while the file still looks
correct.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.processing import find_groove_positions, raw_data
from ..core.region_mask import (GROOVE_REACH_PERIODS, bridge_nan,
                                describe_mask, groove_spacings_px, has_mask,
                                reject_grooves_near_nan,
                                tip_correct_around_nan)
from .average import (average_grooves, flatten_endpoints, groove_residuals,
                      normalize_profile, profile_metrics, select_grooves)

__all__ = ["BoundaryProfile", "build_boundary_profile"]


@dataclass(frozen=True)
class BoundaryProfile:
    """One averaged, normalised groove, and how it was arrived at."""

    #: Normalised for export: x over [0, 1] across one period, y as a fraction
    #: of the period. These two are what lands in the .ggp file.
    x_norm: np.ndarray
    y_norm: np.ndarray

    #: The averaged groove before normalisation, in µm and nm, with the spread
    #: across the grooves that went into it. For the +/-1 sigma band.
    x_avg_um: np.ndarray
    y_avg_nm: np.ndarray
    y_std_nm: np.ndarray

    metrics: dict
    period_nm: float
    n_grooves: int          # detected, after edge exclusion
    n_used: int             # averaged; the rest were too near an edge
    n_edge_rejected: int
    n_mask_rejected: int = 0  # dropped for reaching into a masked column

    #: The full flattened trace and its detected groove centres, so a caller can
    #: show where the averaged groove came from.
    profile_x_um: np.ndarray = field(repr=False, default=None)
    profile_y_nm: np.ndarray = field(repr=False, default=None)
    groove_centers: np.ndarray = field(repr=False, default=None)

    @property
    def summary(self) -> str:
        """One-line description, for a status bar."""
        return (f"{self.n_used} of {self.n_grooves} grooves averaged, "
                f"period {self.period_nm:.2f} nm, "
                f"depth {self.metrics['groove_depth']:.4f} of period")

    @property
    def groove_residual_rms_nm(self):
        """RMS residual of the grooves about their average, in nm.

        None when fewer than two grooves were averaged. What it measures, and
        why `to_problem` does not use it unless asked to, is in that method's
        docstring and in docs/theory/metrology.md section 4.
        """
        return self.metrics.get('groove_residual_rms_nm')

    def to_problem(self, *, coating=None, substrate=None, roughness=0.0,
                   blaze_angle=None):
        """This measured groove as an efficiency problem.

        This is the in-process handoff to the rest of gratinglab, and the
        reason the two halves are one package. Writing a .ggp and reading it
        back preserves the shape but loses the scale: the file format is
        normalised and carries no period, so the period has to be supplied
        again from elsewhere. `benchmarks/corpus.toml` serves that purpose, as
        a sidecar of periods recovered by hand after the fact. Here
        `period_nm` is the period this scan measured, from the spacing of the
        grooves it detected.

        The period therefore enters the result physically. `y_norm` is a
        fraction of the period, so the solver's phase term is
        `(2*pi/lambda) * y_norm * period * ...` and an error in the period
        scales the groove depth with it.

        Parameters:
            coating, substrate: material names, resolved by
                `gratinglab.materials` at solve time.
            roughness: RMS surface roughness in nm, the Nevot-Croce sigma.
                0.0 means "not supplied", not "smooth". It stays 0.0 unless
                the caller passes a number, including the groove residual
                this class measures; the note below says why that is not
                the default.
            blaze_angle: fitted facet angle in degrees, if the same scan was
                also run through `core.analysis.extract_blaze_angle`. Without
                it the scalar solver evaluates reflection on the mean surface,
                which for a sawtooth at grazing incidence is off by the whole
                blaze angle.

        A note on roughness: `Problem.roughness` means high-spatial-frequency
        surface roughness, the quantity with which Nevot-Croce damps a Fresnel
        coefficient. `y_std_nm` is the groove-to-groove form spread and
        `metrics['rms_slope']` is a slope statistic, and neither is that.

        `groove_residual_rms_nm` is the nearest thing this class measures: the
        RMS of each groove about the average once registered to it, which is
        the part of the profile that does not repeat with the period. It is
        still not passed by default, for four reasons, each of which makes it
        a different number from the sigma Nevot-Croce wants:

        - It is measured on the row-averaged profile. Structure that varies
          along a groove averages away over the scan lines, and for isotropic
          roughness with a correlation length well under the scan height that
          is nearly all of it. The number keeps only what is coherent along
          the groove for the height of the scan, so it is a lower bound on the
          surface's roughness in its band.
        - Its band is not the sigma's band. It sees nothing finer than the
          pixel pitch or the tip, a few nanometres, while the structure that
          damps a soft X-ray reflection extends down to wavelengths near the
          X-ray wavelength itself. At the other end it includes groove-scale
          form error, whose correlation length is about a period, where
          Nevot-Croce assumes a short one.
        - Groove-to-groove form error removes light from the orders by
          scattering it between them, a Debye-Waller factor on each order's
          own momentum transfer, not a damping of the Fresnel coefficient at
          the facet. Folding it into the facet reflectivity applies the wrong
          factor to the right effect.
        - It carries instrument terms: AFM noise divided by the square root of
          the number of scan lines, and a floor where the surface has corners
          sharper than a pixel (0.21 nm on the synthetic sawtooth).

        A caller who accepts those can pass it, and the name then appears in
        their own code, which is the point:
        `profile.to_problem(roughness=profile.groove_residual_rms_nm)`.
        """
        # Imported here rather than at module scope: `boundary` is reachable
        # from the GUI preview, which has no reason to pull the solver stack in
        # just to draw a curve.
        from ...problem import Problem
        from ...profiles import FromProfileData

        return Problem(
            period=self.period_nm,
            # x_norm advances along +d-hat, as the .ggp written from it
            # does; the profile parameter runs against it.
            profile=FromProfileData.along_dispersion(
                self.x_norm, self.y_norm, blaze_angle=blaze_angle
            ),
            coating=coating,
            substrate=substrate,
            roughness=roughness,
            # Measured, and free to carry. The scalar solver does
            # not apply an interference factor to order efficiencies, so this
            # is metadata rather than an input -- but it is the real count, and
            # a later method that wants it should not have to guess.
            n_grooves=self.n_used,
        )


def build_boundary_profile(data, scan_x_size, settings, orientation=None) -> BoundaryProfile:
    """
    Average the grooves of one scan into a normalised boundary profile.

    Parameters:
        data: 2-D height array in metres, already image-flattened by the caller
        scan_x_size: scan width in microns
        settings: AnalysisSettings; uses the period estimate, groove detection
            and edge exclusion shared with the blaze-angle path, plus the ggp_*
            fields for the export itself.
        orientation: the OrientationReport for this scan, when the caller
            measured one (core.orientation.apply_groove_orientation). The
            correction itself is the caller's, made on the whole image before
            this is called, because a swath passed in here is too short to
            measure an angle from. This only records it in the metrics, so the
            sidecar says whether the period and shape are along the scan line
            or along the groove normal.

    Returns:
        BoundaryProfile

    Raises:
        ValueError if no grooves are detected, or none survive the minimum
        half-width, since there is nothing to average in either case.
    """
    # Tip correction happens here, on the 2-D array, because this function is
    # the one seam both boundary callers (CLI export and GUI preview) pass
    # through -- the next line averages the rows away. The correction settings
    # and the certain fraction land in the metrics: a corrected depth and an
    # uncorrected one are different measurements and the sidecar must say
    # which one it is describing.
    #
    # The region mask, scar repair and orientation are the callers', all before
    # this point: the mask before levelling, which reads whole lines. It is not
    # re-applied here, because the orientation correction shifts and crops the
    # columns its pixel indices refer to. Excluded pixels arrive as NaN.
    n_masked = int(np.isnan(data).sum())
    data, tip = tip_correct_around_nan(data, scan_x_size, settings)

    raw_x, measured_y = raw_data(data, scan_x_size)
    raw_y = bridge_nan(measured_y)

    # Endpoints to zero, then remove the residual tilt. See the module docstring
    # for why this is not the profile flattening used elsewhere.
    flat_y = flatten_endpoints(raw_x, raw_y)
    # The tilt is fitted to measured columns only: a bridged gap is a straight
    # line, and fitting it would tilt every groove by however much it weighs.
    measured = np.isfinite(measured_y)
    flat_y = flat_y - np.polyval(
        np.polyfit(raw_x[measured], flat_y[measured], 1), raw_x)

    scan_width_nm = scan_x_size * 1000
    period_nm = scan_width_nm / max(2, int(scan_width_nm / settings.period_est))

    groove_centers, n_edge_rejected = find_groove_positions(
        raw_x, flat_y, period_nm,
        prominence_factor=settings.prominence_factor,
        distance_factor=settings.distance_factor,
        edge_exclusion=settings.edge_exclusion_periods,
        return_n_edge_rejected=True)
    groove_centers, n_mask_rejected = reject_grooves_near_nan(
        groove_centers, measured_y,
        GROOVE_REACH_PERIODS * period_nm / ((raw_x[1] - raw_x[0]) * 1000))

    if len(groove_centers) == 0:
        raise ValueError(
            "no grooves detected - check that the estimated period matches "
            "this grating, or lower the detection prominence"
            + (", or shrink the region mask" if n_mask_rejected else ""))

    # Measured period, once there is more than one groove to measure between.
    # A spacing across a masked gap is not a period (groove_spacings_px).
    spacings = groove_spacings_px(groove_centers, measured_y)
    if np.isfinite(spacings).any():
        period_nm = float(np.nanmean(spacings) * (raw_x[1] - raw_x[0]) * 1000)

    x_avg, y_avg, y_std, n_used = average_grooves(
        raw_x, flat_y, groove_centers, period_nm,
        margin=0.0, n_points=settings.ggp_n_points,
        min_half_width=settings.ggp_min_half_width)

    x_norm, y_norm, _edge_height = normalize_profile(
        x_avg, y_avg, period_nm,
        apply_smoothing=settings.ggp_apply_smoothing,
        smoothing_window=settings.ggp_smoothing_window)

    metrics = profile_metrics(x_norm, y_norm, period_nm, n_used)
    # The same grooves average_grooves used, chosen by the same rules.
    residual = groove_residuals(
        select_grooves(raw_x, flat_y, groove_centers, period_nm, margin=0.0,
                       min_half_width=settings.ggp_min_half_width),
        x_avg, y_avg)
    metrics['groove_residual_rms_nm'] = residual['rms_nm']
    metrics['groove_residual_per_groove_nm'] = residual['per_groove_nm'].tolist()
    metrics['groove_residual_max_shift_px'] = (
        float(np.max(np.abs(residual['shift_px'])))
        if residual['rms_nm'] is not None else None)
    if orientation is not None:
        metrics.update(orientation.as_metrics())
    metrics['tip_correction'] = settings.tip_correction
    if tip is not None:
        metrics['tip_radius_nm'] = tip.radius_nm
        metrics['tip_half_angle_deg'] = tip.half_angle_deg
        if not tip.symmetric:
            metrics['tip_half_angle_right_deg'] = tip.half_angle_right_deg
        metrics['tip_certain_fraction'] = tip.certain_fraction
    if n_masked and has_mask(settings):
        metrics['region_mask'] = describe_mask(settings)
        metrics['masked_pixels'] = n_masked
        metrics['mask_rejected_grooves'] = n_mask_rejected

    return BoundaryProfile(
        x_norm=x_norm,
        y_norm=y_norm,
        x_avg_um=x_avg,
        y_avg_nm=y_avg,
        y_std_nm=y_std,
        metrics=metrics,
        period_nm=period_nm,
        n_grooves=len(groove_centers),
        n_used=n_used,
        n_edge_rejected=n_edge_rejected,
        n_mask_rejected=n_mask_rejected,
        profile_x_um=raw_x,
        profile_y_nm=flat_y,
        groove_centers=groove_centers,
    )
