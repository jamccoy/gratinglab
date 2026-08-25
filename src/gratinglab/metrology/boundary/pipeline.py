"""
Building a boundary profile, separated from writing one.

The computation used to live inside ``workflows.run_boundary_profile_export``,
fused with the file writing. Splitting it lets the GUI preview a profile - and
recompute it as controls move - without touching the disk, while the CLI mode
writes exactly what the panel showed.

The flattening here is deliberately *not* the profile flattening used for
blaze-angle work. ``flatten_endpoints`` forces the ends of the trace to zero so
one groove tiles seamlessly into the next; profile flattening removes a fitted
background so a facet can be measured. Different jobs, and swapping one for the
other would produce a profile with a step at the period boundary - which is
exactly the defect that makes a PCGrate efficiency curve wrong while the file
still looks fine.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.processing import find_groove_positions, raw_data
from ..core.tip import apply_tip_correction
from .average import (average_grooves, flatten_endpoints, normalize_profile,
                      profile_metrics)

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
    n_used: int             # actually averaged; the rest were too near an edge
    n_edge_rejected: int

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

    def to_problem(self, *, coating=None, substrate=None, roughness=0.0,
                   blaze_angle=None):
        """This measured groove as an efficiency problem.

        The in-process handoff to the rest of gratinglab, and the reason the
        two halves are one package. Writing a .ggp and reading it back gets
        the same *shape* but loses the scale: the file format is normalised
        and carries no period, so the period has to be supplied again from
        somewhere -- which is what `benchmarks/corpus.toml` is, a sidecar of
        periods recovered by hand after the fact. Here `period_nm` is the one
        this scan measured, from the spacing of the grooves it detected.

        That matters beyond bookkeeping. `y_norm` is a fraction of the period,
        so the solver's phase term is `(2*pi/lambda) * y_norm * period * ...`
        and an error in the period scales the groove depth with it. The period
        is not a label on the result; it is inside the physics.

        Parameters:
            coating, substrate: material names, resolved by
                `gratinglab.materials` at solve time.
            roughness: RMS surface roughness in nm. See the note below --
                this class does not measure it, and 0.0 means "not supplied",
                not "smooth".
            blaze_angle: fitted facet angle in degrees, if the same scan was
                also run through `core.analysis.extract_blaze_angle`. Without
                it the scalar solver evaluates reflection on the mean surface,
                which for a sawtooth at grazing incidence is off by the whole
                blaze angle.

        **On roughness.** `Problem.roughness` means high-spatial-frequency
        surface roughness, the quantity Nevot-Croce damps a Fresnel
        coefficient with. This class computes no such number. `y_std_nm` is
        the groove-to-groove *form* spread and `metrics['rms_slope']` is a
        slope statistic; neither is surface microroughness, and passing either
        one in here would produce a confident, wrong efficiency rather than an
        obviously missing one. The honest quantity -- the RMS residual of each
        groove about this average -- is not computed yet.
        """
        # Imported here rather than at module scope: `boundary` is reachable
        # from the GUI preview, which has no reason to pull the solver stack in
        # just to draw a curve.
        from ...problem import Problem
        from ...profiles import FromProfileData

        return Problem(
            period=self.period_nm,
            profile=FromProfileData(
                t=tuple(float(v) for v in self.x_norm),
                y=tuple(float(v) for v in self.y_norm),
                blaze_angle=blaze_angle,
            ),
            coating=coating,
            substrate=substrate,
            roughness=roughness,
            # Measured, and free to carry. The scalar solver deliberately does
            # not apply an interference factor to order efficiencies, so this
            # is metadata rather than an input -- but it is the real count, and
            # a later method that wants it should not have to guess.
            n_grooves=self.n_used,
        )


def build_boundary_profile(data, scan_x_size, settings) -> BoundaryProfile:
    """
    Average the grooves of one scan into a normalised boundary profile.

    Parameters:
        data: 2-D height array in metres, already image-flattened by the caller
        scan_x_size: scan width in microns
        settings: AnalysisSettings; uses the period estimate, groove detection
            and edge exclusion shared with the blaze-angle path, plus the ggp_*
            fields for the export itself.

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
    data, tip = apply_tip_correction(data, scan_x_size, settings)

    raw_x, raw_y = raw_data(data, scan_x_size)

    # Endpoints to zero, then remove the residual tilt. See the module docstring
    # for why this is not the profile flattening used elsewhere.
    flat_y = flatten_endpoints(raw_x, raw_y)
    flat_y = flat_y - np.polyval(np.polyfit(raw_x, flat_y, 1), raw_x)

    scan_width_nm = scan_x_size * 1000
    period_nm = scan_width_nm / max(2, int(scan_width_nm / settings.period_est))

    groove_centers, n_edge_rejected = find_groove_positions(
        raw_x, flat_y, period_nm,
        prominence_factor=settings.prominence_factor,
        distance_factor=settings.distance_factor,
        edge_exclusion=settings.edge_exclusion_periods,
        return_n_edge_rejected=True)

    if len(groove_centers) == 0:
        raise ValueError(
            "no grooves detected - check that the estimated period matches "
            "this grating, or lower the detection prominence")

    # Measured period, once there is more than one groove to measure between.
    if len(groove_centers) > 1:
        period_nm = float(np.mean(
            np.diff(groove_centers) * (raw_x[1] - raw_x[0]) * 1000))

    x_avg, y_avg, y_std, n_used = average_grooves(
        raw_x, flat_y, groove_centers, period_nm,
        margin=0.0, n_points=settings.ggp_n_points,
        min_half_width=settings.ggp_min_half_width)

    x_norm, y_norm, _edge_height = normalize_profile(
        x_avg, y_avg, period_nm,
        apply_smoothing=settings.ggp_apply_smoothing,
        smoothing_window=settings.ggp_smoothing_window)

    metrics = profile_metrics(x_norm, y_norm, period_nm, n_used)
    metrics['tip_correction'] = settings.tip_correction
    if tip is not None:
        metrics['tip_radius_nm'] = tip.radius_nm
        metrics['tip_half_angle_deg'] = tip.half_angle_deg
        metrics['tip_certain_fraction'] = tip.certain_fraction

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
        profile_x_um=raw_x,
        profile_y_nm=flat_y,
        groove_centers=groove_centers,
    )
