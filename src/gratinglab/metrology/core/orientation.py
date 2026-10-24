"""
Groove orientation: measuring, and undoing, a sample rotated on the stage.

When the grooves are not perpendicular to the fast scan axis, by the degree or
two a sample is easily rotated when it is mounted, every scan line crosses them
obliquely. Two things then go wrong, and nothing downstream can see either:

1. **Each row is stretched.** A row crosses one period in ``d / cos(theta)`` of
   travel, so the measured period is long by ``1 / cos(theta)`` and every facet
   slope is shallow by ``cos(theta)``: ``tan(blaze_measured) =
   tan(blaze) * cos(theta)``.
2. **Averaging rows smears the groove.** The grooves sit at a different x in
   every row, drifting by ``tan(theta)`` per unit of slow-axis travel, so a row
   group, and much more so the whole-image average the boundary profile is
   built from, averages a groove with displaced copies of itself.

The first is small (0.14% of the period at 3 degrees) and the second is not:
on a 2 µm square scan at 1 degree the grooves drift 35 nm from the top row to
the bottom one, 11% of a 315 nm period.

**Measuring it.** The image is cut into horizontal bands and the grooves are
found in each band with the same `find_groove_positions` the analysis uses,
refined to a fraction of a pixel by a parabola through the minimum. Positions
are linked band to band into one track per groove, and one common slope
``dx/dy`` is fitted to all of the tracks at once, each with its own intercept
(a pooled within-track regression). ``theta = atan(dx/dy)``, and its standard
error comes from the scatter of the positions about that fit. A groove that
cannot be linked starts a track of its own, and a track seen in one band only
carries no slope information and is dropped from the fit.

This needs the slow-axis pixel size, since ``dx/dy`` is a ratio of lengths and
the image only gives pixels. It comes from the file (`slow_axis_size_um`), and
when the file does not say, the angle is not measured: an assumed square scan
would turn a 2:1 scan's rotation into a wrong number with no sign that it was
assumed.

Sign: ``theta`` is positive when the grooves move toward larger x as the row
index increases, that is, down the array as it is stored.

**Correcting it.** Each row is shifted along x, by linear interpolation, so the
grooves line up from row to row, and the fast-axis pixel is then taken as
``dx * cos(theta)`` instead of ``dx``. For a grating, whose surface does not
change along the grooves, this is exact: a row of the rotated surface is
``f(cos(theta) * (x + y * tan(theta)))``, so the shift removes the ``y`` term
and the pixel change removes the ``cos(theta)``. Heights are untouched, so
slopes come out right on their own.

That was chosen over rotating the 2-D image (``scipy.ndimage.rotate``) for
three reasons. A rotation interpolates in two dimensions where this
interpolates in one, along the row, so a facet is resampled along its own
slope and nowhere else. A rotation leaves four invalid corners that must be
cropped from both axes, where this only loses the columns the shift slides off
the ends of a row, at most ``tan(theta)`` of the scan height. And a rotation
needs square pixels, so a scan with different fast and slow pixel sizes would
have to be resampled first, which this never does. What it gives up is
generality: it is right for a surface that is invariant along one direction,
which is what a grating is and what this package measures.

The shift is centred on the middle row, so the columns lost are split between
the two ends. Cropped columns are reported.

The interpolation has a cost, and it is measured rather than assumed. A
corner sharper than a pixel is blunted by it, by about one pixel over one
period: on an ideal knife-edge sawtooth (the test suite's) the corrected
whole-image groove reads 3% shallow at 33 pixels per period, 1.5% at 66 and
under 0.5% at 120. Left uncorrected, the same grooves lose 9% of their depth
to smearing at 1 degree and 32% at 3 degrees. A real groove is rounded by the
tip over several nanometres, so it has no sub-pixel corner to blunt.

**When.** `AnalysisSettings.groove_orientation` chooses:

``auto`` (the default)
    Measure, and correct when ``|theta|`` is at least
    ``groove_orientation_threshold_deg`` (0.2 degrees) *and* at least three
    standard errors from zero. Below 0.2 degrees the whole-image drift on a
    square scan is under 0.4% of the scan width, smaller than the tip's own
    effect on the groove, and correcting a rotation the data cannot tell from
    zero only adds interpolation.
``always``
    Correct whatever was measured. Refuses, by raising, when the angle cannot
    be measured, because the request was explicit.
``detect``
    Measure and report; change nothing.
``off``
    Do nothing, and report that nothing was done.

Every result carries a `OrientationReport`, applied or not, so a reader can see
the angle and why it was or was not corrected.
"""
from __future__ import annotations

import contextlib
import io
import math
from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np

from .processing import find_groove_positions

__all__ = [
    "VALID_GROOVE_ORIENTATIONS",
    "DEFAULT_THRESHOLD_DEG",
    "OrientationMeasurement",
    "OrientationReport",
    "measure_groove_orientation",
    "align_grooves",
    "apply_groove_orientation",
]

VALID_GROOVE_ORIENTATIONS = ('auto', 'always', 'detect', 'off')

#: Below this the correction is not worth its interpolation; see the module
#: docstring for the arithmetic.
DEFAULT_THRESHOLD_DEG = 0.2

#: ``auto`` also needs the angle this many standard errors from zero.
SIGNIFICANCE = 3.0

#: Bands the image is cut into for measuring. Each must keep at least
#: ``MIN_ROWS_PER_BAND`` rows so the band average still suppresses line noise.
MAX_BANDS = 32
MIN_ROWS_PER_BAND = 3
MIN_BANDS = 3

#: A groove links to a track ending within this fraction of a period.
LINK_FRACTION = 0.35

#: The correction is refused when it would keep less than this fraction of the
#: columns; that is a rotation of tens of degrees, not a mounting error.
MIN_KEPT_FRACTION = 0.5


@dataclass(frozen=True)
class OrientationMeasurement:
    """The fitted groove angle and what it was fitted from."""

    angle_deg: float
    stderr_deg: float
    n_tracks: int        # grooves followed through two or more bands
    n_points: int        # positions in the fit
    n_bands: int
    pixel_x_nm: float
    pixel_y_nm: float


@dataclass(frozen=True)
class OrientationReport:
    """What the orientation step measured and did, for results and sidecars."""

    mode: str
    applied: bool
    reason: str
    angle_deg: Optional[float] = None
    stderr_deg: Optional[float] = None
    n_tracks: int = 0
    n_points: int = 0
    scan_y_size_um: Optional[float] = None
    columns_cropped: int = 0
    scan_x_size_in_um: Optional[float] = None
    scan_x_size_out_um: Optional[float] = None

    @property
    def measured(self) -> bool:
        return self.angle_deg is not None

    def summary(self) -> str:
        """One line for a log, a results panel or a metrics sidecar."""
        if not self.measured:
            return f"groove orientation not measured: {self.reason}"
        text = (f"groove rotation {self.angle_deg:+.3f} ± "
                f"{self.stderr_deg:.3f} deg "
                f"({self.n_tracks} grooves, {self.n_points} positions)")
        if self.applied:
            text += (f", corrected: {self.columns_cropped} columns cropped, "
                     f"width {self.scan_x_size_in_um:.4f} -> "
                     f"{self.scan_x_size_out_um:.4f} µm along the groove normal")
        else:
            text += f", not corrected: {self.reason}"
        return text

    def as_dict(self) -> dict:
        """Plain values, for a results dictionary or a JSON sidecar."""
        return asdict(self)

    def as_metrics(self) -> dict:
        """The keys a boundary profile's metrics carry."""
        return {
            'groove_orientation': self.mode,
            'groove_rotation_deg': self.angle_deg,
            'groove_rotation_stderr_deg': self.stderr_deg,
            'groove_rotation_corrected': self.applied,
            'groove_rotation_note': self.reason,
        }


def _subpixel(profile, index):
    """Vertex of the parabola through a minimum and its two neighbours."""
    if index <= 0 or index >= len(profile) - 1:
        return float(index)
    left, mid, right = profile[index - 1], profile[index], profile[index + 1]
    curvature = left - 2.0 * mid + right
    if curvature <= 0:
        return float(index)
    offset = 0.5 * (left - right) / curvature
    return index + float(np.clip(offset, -0.5, 0.5))


def _bands(n_rows):
    n_bands = min(MAX_BANDS, n_rows // MIN_ROWS_PER_BAND)
    if n_bands < MIN_BANDS:
        return []
    edges = np.linspace(0, n_rows, n_bands + 1).round().astype(int)
    return [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:]) if b > a]


def _link(band_positions, period_nm):
    """Positions per band -> tracks of (y_nm, x_nm), one per groove."""
    tracks = []
    for y_nm, xs in band_positions:
        ends = np.array([t[-1][1] for t in tracks]) if tracks else np.empty(0)
        claimed = {}
        for x in xs:
            if len(ends):
                nearest = int(np.argmin(np.abs(ends - x)))
                distance = abs(ends[nearest] - x)
                if distance <= LINK_FRACTION * period_nm:
                    held = claimed.get(nearest)
                    if held is None or distance < held[0]:
                        claimed[nearest] = (distance, x)
                    continue
            tracks.append([(y_nm, x)])
        for index, (_distance, x) in claimed.items():
            tracks[index].append((y_nm, x))
    return [t for t in tracks if len(t) >= 2]


def _pooled_slope(tracks):
    """Common dx/dy across tracks, each with its own intercept."""
    dy, dx, groups = [], [], []
    for number, track in enumerate(tracks):
        ys = np.array([p[0] for p in track])
        xs = np.array([p[1] for p in track])
        dy.append(ys - ys.mean())
        dx.append(xs - xs.mean())
        groups.append(np.full(len(ys), number))
    dy, dx, groups = np.concatenate(dy), np.concatenate(dx), np.concatenate(groups)
    sxx = float(np.sum(dy * dy))
    if sxx <= 0:
        return None
    slope = float(np.sum(dy * dx) / sxx)
    residual = dx - slope * dy
    dof = len(dy) - len(tracks) - 1
    if dof < 1:
        return None
    sigma2 = float(np.sum(residual ** 2) / dof)
    return slope, math.sqrt(sigma2 / sxx), residual, groups, dof


def measure_groove_orientation(data, scan_x_size_um, scan_y_size_um, period_nm,
                               prominence_factor=0.01, distance_factor=0.3):
    """
    Measure how far the grooves are rotated from the slow scan axis.

    Parameters:
        data: 2-D height array, rows are scan lines, ideally image-flattened
        scan_x_size_um, scan_y_size_um: the physical size of each axis
        period_nm: the period estimate, as the analysis uses it
        prominence_factor, distance_factor: groove detection, as the analysis

    Returns:
        OrientationMeasurement, or a string saying why there is none.
    """
    data = np.asarray(data, dtype=float)
    n_rows, n_cols = data.shape
    if n_cols < 3:
        return "the image has fewer than three columns"
    pixel_x = scan_x_size_um * 1000.0 / (n_cols - 1)
    pixel_y = scan_y_size_um * 1000.0 / max(n_rows - 1, 1)
    bands = _bands(n_rows)
    if not bands:
        return (f"{n_rows} rows make fewer than {MIN_BANDS} bands of "
                f"{MIN_ROWS_PER_BAND}")

    x_um = np.arange(n_cols) * pixel_x / 1000.0
    band_positions = []
    for first, end in bands:
        profile = np.nanmean(data[first:end, :], axis=0) * 1e9
        if not np.all(np.isfinite(profile)):
            continue
        profile = profile - np.polyval(np.polyfit(x_um, profile, 1), x_um)
        with contextlib.redirect_stdout(io.StringIO()):
            centers = find_groove_positions(
                x_um, profile, period_nm, prominence_factor=prominence_factor,
                distance_factor=distance_factor, edge_exclusion=0.5)
        if len(centers):
            xs = [_subpixel(profile, int(c)) * pixel_x for c in centers]
            band_positions.append((0.5 * (first + end - 1) * pixel_y, xs))

    if len(band_positions) < MIN_BANDS:
        return (f"grooves were found in {len(band_positions)} of {len(bands)} "
                "bands, too few to follow")

    tracks = _link(band_positions, period_nm)
    if not tracks:
        return "no groove could be followed from one band to the next"
    fit = _pooled_slope(tracks)
    if fit is None:
        return "too few linked positions for a slope and its uncertainty"

    # One pass of outlier rejection: a groove linked to its neighbour in one
    # band puts a single point a period away, which a least-squares slope
    # would follow.
    slope, stderr, residual, groups, _dof = fit
    scale = 1.4826 * float(np.median(np.abs(residual)))
    if scale > 0:
        keep = np.abs(residual) <= 4.0 * scale
        if not keep.all():
            flat = [p for t in tracks for p in t]
            survivors = {}
            for point, group, kept in zip(flat, groups, keep):
                if kept:
                    survivors.setdefault(int(group), []).append(point)
            trimmed = [t for t in survivors.values() if len(t) >= 2]
            refit = _pooled_slope(trimmed) if trimmed else None
            if refit is not None:
                tracks = trimmed
                slope, stderr = refit[0], refit[1]

    angle = math.degrees(math.atan(slope))
    angle_err = math.degrees(stderr / (1.0 + slope * slope))
    return OrientationMeasurement(
        angle_deg=angle, stderr_deg=angle_err, n_tracks=len(tracks),
        n_points=sum(len(t) for t in tracks), n_bands=len(band_positions),
        pixel_x_nm=pixel_x, pixel_y_nm=pixel_y)


def align_grooves(data, scan_x_size_um, scan_y_size_um, angle_deg):
    """
    Undo a groove rotation: shift each row, then shrink the fast-axis pixel.

    Returns:
        (data, scan_x_size_um, columns_cropped). The width returned follows the
        loaders' convention, ``pixel = width / (columns - 1)``, with the pixel
        multiplied by ``cos(theta)``: it is the width along the groove normal.

    Raises:
        ValueError when the shift would leave under half the columns.
    """
    data = np.asarray(data, dtype=float)
    n_rows, n_cols = data.shape
    pixel_x = scan_x_size_um * 1000.0 / (n_cols - 1)
    pixel_y = scan_y_size_um * 1000.0 / max(n_rows - 1, 1)
    theta = math.radians(angle_deg)

    rows = np.arange(n_rows) - 0.5 * (n_rows - 1)
    shifts = math.tan(theta) * rows * pixel_y / pixel_x     # in columns
    first = int(math.ceil(max(0.0, -shifts.min()) - 1e-9))
    last = int(math.floor(n_cols - 1 - max(0.0, shifts.max()) + 1e-9))
    kept = last - first + 1
    if kept < MIN_KEPT_FRACTION * n_cols or kept < 3:
        raise ValueError(
            f"a {angle_deg:.2f} deg correction would keep {kept} of {n_cols} "
            "columns; that is not a mounting error this can undo")

    columns = np.arange(n_cols, dtype=float)
    wanted = np.arange(first, last + 1, dtype=float)
    aligned = np.empty((n_rows, kept))
    for r in range(n_rows):
        aligned[r] = np.interp(wanted + shifts[r], columns, data[r])

    width = pixel_x * math.cos(theta) * (kept - 1) / 1000.0
    return aligned, width, n_cols - kept


def _period_estimate(scan_x_size_um, period_est):
    """The analysis's own rounding of the period estimate to the scan width."""
    width = scan_x_size_um * 1000.0
    return width / max(2, int(width / period_est))


def apply_groove_orientation(data, scan_x_size_um, scan_y_size_um, settings):
    """
    Measure the groove rotation and correct it as `settings` asks.

    The one entry point the pipeline calls, after image flattening and before
    tip correction. Tip correction comes after deliberately: an axisymmetric
    tip on a surface that is invariant along the grooves reduces to a 1-D
    dilation along the groove normal, which is the profile this produces and
    not the oblique one the scan recorded.

    Parameters:
        scan_y_size_um: the slow-axis size from the file, or None when the
            file does not record it.

    Returns:
        (data, scan_x_size_um, OrientationReport). The data and width are
        unchanged unless the correction was applied.
    """
    mode = getattr(settings, 'groove_orientation', 'auto')
    threshold = getattr(settings, 'groove_orientation_threshold_deg',
                        DEFAULT_THRESHOLD_DEG)
    if mode not in VALID_GROOVE_ORIENTATIONS:
        raise ValueError(f"groove_orientation must be one of "
                         f"{', '.join(VALID_GROOVE_ORIENTATIONS)}, not {mode!r}")

    if mode == 'off':
        return data, scan_x_size_um, OrientationReport(
            mode=mode, applied=False, reason="switched off",
            scan_y_size_um=scan_y_size_um)

    if scan_y_size_um is None or not scan_y_size_um > 0:
        reason = ("the file does not record the slow-axis scan size, and a "
                  "square scan is not assumed")
        if mode == 'always':
            raise ValueError(f"groove orientation correction was asked for, but {reason}")
        return data, scan_x_size_um, OrientationReport(
            mode=mode, applied=False, reason=reason)

    measured = measure_groove_orientation(
        data, scan_x_size_um, scan_y_size_um,
        _period_estimate(scan_x_size_um, settings.period_est),
        prominence_factor=settings.prominence_factor,
        distance_factor=settings.distance_factor)
    if isinstance(measured, str):
        if mode == 'always':
            raise ValueError(f"groove orientation could not be measured: {measured}")
        return data, scan_x_size_um, OrientationReport(
            mode=mode, applied=False, reason=measured,
            scan_y_size_um=scan_y_size_um)

    base = dict(mode=mode, angle_deg=measured.angle_deg,
                stderr_deg=measured.stderr_deg, n_tracks=measured.n_tracks,
                n_points=measured.n_points, scan_y_size_um=scan_y_size_um)
    angle = abs(measured.angle_deg)

    if mode == 'detect':
        return data, scan_x_size_um, OrientationReport(
            applied=False, reason="detect only", **base)
    if mode == 'auto':
        if angle < threshold:
            return data, scan_x_size_um, OrientationReport(
                applied=False, reason=f"below the {threshold:g} deg threshold",
                **base)
        if angle < SIGNIFICANCE * measured.stderr_deg:
            return data, scan_x_size_um, OrientationReport(
                applied=False,
                reason=f"within {SIGNIFICANCE:g} standard errors of zero",
                **base)

    try:
        aligned, width, cropped = align_grooves(
            data, scan_x_size_um, scan_y_size_um, measured.angle_deg)
    except ValueError as exc:
        if mode == 'always':
            raise
        return data, scan_x_size_um, OrientationReport(
            applied=False, reason=str(exc), **base)

    reason = "always" if mode == 'always' else (
        f"at least {threshold:g} deg and {SIGNIFICANCE:g} standard errors")
    return aligned, width, OrientationReport(
        applied=True, reason=reason, columns_cropped=cropped,
        scan_x_size_in_um=scan_x_size_um, scan_x_size_out_um=width, **base)
