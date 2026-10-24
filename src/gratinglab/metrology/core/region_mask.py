"""
Region masking: leaving part of a scan out of the analysis.

A particle on the surface, a scratch, a patch where the tip picked something up:
Gwyddion's answer is to draw a mask and crop or exclude it before export. Here
the mask is a *setting*, so a batch or CLI run reproduces exactly what was drawn
in the window, and it is recorded beside the result.

Three shapes, all in pixel indices of the image as loaded (row 0 is the first
scan line in the file, column 0 the left edge), half-open like a Python slice:

- ``mask_rect``: one rectangle ``(row0, row1, col0, col1)``, or None;
- ``mask_rows``: whole scan-line ranges ``((row0, row1), ...)``;
- ``mask_cols``: whole column ranges ``((col0, col1), ...)``.

Pixels rather than micrometres because the slow axis has no reliable physical
size in every format, and because the mask belongs to one scan: it is meant to
reproduce that scan's analysis, not to transfer between scans of another
resolution. Ranges that run past the image are clipped.

How excluded pixels leave the analysis
--------------------------------------
They become NaN **first**, straight after loading and before image flattening,
because levelling reads the whole image: ``align_rows`` takes each line's
median and ``plane`` fits every pixel, so a particle left in would shift every
line it sits on before any later mask could hide it. Both methods skip NaN.

Tip erosion cannot take a NaN (its relief is ``data.max() - data.min()``), so
``tip_correct_around_nan`` hands it a copy with each gap bridged along its line
and puts the NaN back afterwards. Pixels next to a masked particle are then
eroded against the bridge rather than against the particle, which is what
excluding the particle means.

From there:

- A scan line that is entirely NaN is skipped by row grouping, so masked rows
  contribute to no row group and the groups are formed from what is left.
- Within a row group, each column is averaged over the lines that are not
  masked there. A column masked across a whole group is NaN in that group's
  profile, and any groove whose analysis window touches such a column is
  discarded (``reject_grooves_near_nan``), so masked columns contribute no
  grooves.

The profile still has to be flattened and searched for grooves across the gap,
so NaN columns are bridged by linear interpolation for those two steps only
(``bridge_nan``). The bridged values are never measured: any groove that would
read them is the one rejected above.
"""
from __future__ import annotations

import warnings

import numpy as np

__all__ = ["region_mask", "apply_region_mask", "describe_mask", "has_mask",
           "bridge_nan", "reject_grooves_near_nan", "nan_mean_rows",
           "pixels_from_drag", "tip_correct_around_nan",
           "GROOVE_REACH_PERIODS"]

#: How far either side of its centre a groove's analysis reads, in periods.
#: ``extract_blaze_angle`` takes half a period plus a 20% margin (0.6), and
#: the boundary average half a period; 0.65 leaves room for a local period a
#: little longer than the estimate.
GROOVE_REACH_PERIODS = 0.65


def tip_correct_around_nan(data, scan_x_size, settings):
    """
    ``apply_tip_correction`` for an image that may hold NaN.

    The gaps are bridged along each line for the erosion and restored after,
    so NaN in means the same NaN out. With nothing masked, this is exactly
    ``apply_tip_correction``.
    """
    from .tip import apply_tip_correction

    gaps = ~np.isfinite(data)
    if not gaps.any() or settings.tip_correction == 'none':
        return apply_tip_correction(data, scan_x_size, settings)
    filled = np.array(data, dtype=float)
    columns = np.arange(filled.shape[1])
    for row in np.flatnonzero(gaps.any(axis=1)):
        known = ~gaps[row]
        if known.any():
            filled[row, ~known] = np.interp(columns[~known], columns[known],
                                            filled[row, known])
        else:
            filled[row] = 0.0  # a whole masked line; restored to NaN below
    corrected, correction = apply_tip_correction(filled, scan_x_size, settings)
    corrected = np.array(corrected, dtype=float)
    corrected[gaps] = np.nan
    return corrected, correction


def has_mask(settings) -> bool:
    return bool(getattr(settings, 'mask_rect', None)
                or getattr(settings, 'mask_rows', ())
                or getattr(settings, 'mask_cols', ()))


def region_mask(shape, settings) -> np.ndarray:
    """Boolean array of ``shape``, True where the settings exclude a pixel."""
    n_rows, n_cols = shape
    mask = np.zeros(shape, dtype=bool)
    rect = getattr(settings, 'mask_rect', None)
    if rect:
        r0, r1, c0, c1 = (int(v) for v in rect)
        mask[max(r0, 0):max(r1, 0), max(c0, 0):max(c1, 0)] = True
    for r0, r1 in getattr(settings, 'mask_rows', ()) or ():
        mask[max(int(r0), 0):max(int(r1), 0), :] = True
    for c0, c1 in getattr(settings, 'mask_cols', ()) or ():
        mask[:, max(int(c0), 0):max(int(c1), 0)] = True
    return mask


def apply_region_mask(data, settings):
    """
    ``data`` with the excluded pixels set to NaN, and how many there were.

    Returns ``(data, n_excluded)``. With no mask set the input array is
    returned unchanged, so an unmasked run is bit-identical to before.
    """
    if not has_mask(settings):
        return data, 0
    data = np.array(data, dtype=float)
    mask = region_mask(data.shape, settings)
    data[mask] = np.nan
    return data, int(mask.sum())


def describe_mask(settings) -> str:
    """One line for logs, results files and the window. Empty when unmasked."""
    parts = []
    rect = getattr(settings, 'mask_rect', None)
    if rect:
        r0, r1, c0, c1 = rect
        parts.append(f"rect rows {r0}-{r1} cols {c0}-{c1}")
    rows = getattr(settings, 'mask_rows', ()) or ()
    if rows:
        parts.append("rows " + ", ".join(f"{a}-{b}" for a, b in rows))
    cols = getattr(settings, 'mask_cols', ()) or ()
    if cols:
        parts.append("cols " + ", ".join(f"{a}-{b}" for a, b in cols))
    return "; ".join(parts)


def pixels_from_drag(x0, y0, x1, y1, shape, scan_x_size):
    """
    A drag on the displayed image, as ``(row0, row1, col0, col1)`` pixels.

    The image is drawn with ``extent=[0, scan_x_size, n_rows, 0]``: x in
    micrometres across ``n_cols`` equal pixels, y in scan lines with line 0 at
    the top, so line ``r`` covers ``r <= y < r + 1``. Every pixel the drag
    touches is included. Kept here, without Qt, so the mapping is testable.
    """
    n_rows, n_cols = shape
    x_lo, x_hi = sorted((float(x0), float(x1)))
    y_lo, y_hi = sorted((float(y0), float(y1)))
    col0 = int(np.clip(np.floor(x_lo / scan_x_size * n_cols), 0, n_cols - 1))
    col1 = int(np.clip(np.floor(x_hi / scan_x_size * n_cols), 0, n_cols - 1)) + 1
    row0 = int(np.clip(np.floor(y_lo), 0, n_rows - 1))
    row1 = int(np.clip(np.floor(y_hi), 0, n_rows - 1)) + 1
    return row0, row1, col0, col1


def nan_mean_rows(data) -> np.ndarray:
    """
    Column means over the lines that are not NaN.

    Uses plain ``np.mean`` when nothing is NaN, so a scan without masks or
    dropped lines averages exactly as it always did and the stored baselines
    still match. A column NaN on every line stays NaN.
    """
    if np.isfinite(data).all():
        return np.mean(data, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)  # all-NaN columns
        return np.nanmean(data, axis=0)


def bridge_nan(y) -> np.ndarray:
    """``y`` with NaN runs replaced by straight lines between their neighbours."""
    y = np.asarray(y, dtype=float)
    bad = ~np.isfinite(y)
    if not bad.any():
        return y
    good = ~bad
    if not good.any():
        raise ValueError("the whole profile is masked; nothing left to analyse")
    idx = np.arange(len(y))
    out = y.copy()
    out[bad] = np.interp(idx[bad], idx[good], y[good])
    return out


def reject_grooves_near_nan(centers, profile, half_width_px):
    """
    Drop groove centres whose window of +/- ``half_width_px`` touches a NaN.

    Returns ``(kept_centres, n_rejected)``. ``profile`` is the averaged
    profile before ``bridge_nan``, so its NaN columns are exactly the masked
    ones.
    """
    centers = np.asarray(centers, dtype=int)
    bad = ~np.isfinite(np.asarray(profile, dtype=float))
    if not bad.any() or centers.size == 0:
        return centers, 0
    n = len(bad)
    half = int(np.ceil(half_width_px))
    keep = np.array([not bad[max(c - half, 0):min(c + half + 1, n)].any()
                     for c in centers], dtype=bool)
    return centers[keep], int((~keep).sum())
