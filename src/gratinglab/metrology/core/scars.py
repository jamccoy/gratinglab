"""
Scar repair: finding scan lines that disagree with their neighbours, and fixing them.

A scar is a short run of one or a few scan lines that sits above (or below)
the lines on both sides of it: the tip skipping over a contaminant, the
feedback ringing after a step, a line where the Z piezo briefly lost track.
Gwyddion calls this "Mark scars" and "Remove scars", and the detection here
follows the same idea, so a user moving across from Gwyddion meets the same
behaviour with the same three knobs.

**Why this is a stage of its own and not an entry in IMAGE_FLATTEN_METHODS.**
That registry is a choice of one: a scan is levelled by ``none``, ``plane`` or
``align_rows``, never two at once. Scar repair is orthogonal to levelling and
wants to run *after* it, because an un-levelled row offset is exactly the
signature a scar detector looks for, and every row would be marked. Making it
a registry entry would force a choice between levelling and repair, or a
combinatorial set of entries (``align_rows+scars``, ``plane+scars``), which is
the thing a registry exists to avoid. So it has its own setting,
``scar_repair``, and runs between image flattening and tip correction.

Before tip correction rather than after, because erosion spreads each pixel
sideways by the tip's reach, so a scar left in place would be widened into its
neighbours before anything could find it.

Detection
---------
The vertical (row-to-row) difference of the image is noise plus whatever the
grooves contribute, and for a grating whose grooves run along the slow axis
that is mostly noise. Its scale ``sigma`` is taken robustly (1.4826 x MAD), so
the scars themselves and any steep facet cannot inflate it.

A band of ``w`` adjacent lines (``1 <= w <= max_width``) is a scar candidate at
a pixel where every line of the band lies above both the line just above the
band and the line just below it by more than ``threshold * sigma``, or below
both by as much. Candidates are grown along the line while the same test holds
at half that threshold (hysteresis, as in Gwyddion), and a run is kept only if
it is at least ``min_length`` pixels long. The length rule is what keeps noise
out: a 4-sigma excursion on both sides is rare, and sixteen of them in a row
along a line essentially never happens by chance.

A local maximum in the *vertical* direction is required, not just a large
vertical difference, which is why a groove that runs at a slight angle to the
slow axis (a monotonic shift from line to line) is not mistaken for a scar.

Repair
------
``interpolate`` replaces each marked pixel by linear interpolation, down its
own column, between the nearest unmarked pixels above and below. Down the
column is deliberate: grooves run that way, so the neighbouring lines are the
best estimate of what the tip should have seen.

``drop`` leaves the values alone and removes every line that holds a marked
pixel from the analysis, by setting it to NaN. Row grouping skips NaN lines
(``core.processing.raw_data_multi_group``), so a dropped line contributes to
no row group. This trades data for not inventing any.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["VALID_SCAR_REPAIRS", "DEFAULT_SCAR_REPAIR", "ScarReport",
           "mark_scars", "repair_scars", "repair_scars_for"]

VALID_SCAR_REPAIRS = ('none', 'interpolate', 'drop')
DEFAULT_SCAR_REPAIR = 'none'

#: MAD to standard deviation, for Gaussian noise.
_MAD_TO_SIGMA = 1.4826


@dataclass(frozen=True)
class ScarReport:
    """What scar repair found and did, for the log and the results file."""

    method: str
    n_lines: int        # scan lines holding at least one scar pixel
    n_pixels: int       # scar pixels in total
    lines: tuple        # the indices of those lines, ascending

    @property
    def summary(self) -> str:
        if self.method == 'none':
            return "scar repair off"
        if not self.n_lines:
            return f"{self.method}: no scars found"
        verb = "interpolated" if self.method == 'interpolate' else "dropped"
        return (f"{self.method}: {self.n_lines} line(s) {verb}, "
                f"{self.n_pixels} scar pixel(s)")


def _runs(flags):
    """(start, stop) of each run of True in a 1-D boolean array."""
    padded = np.concatenate(([False], flags, [False]))
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    return edges[0::2], edges[1::2]


def _robust_sigma(data):
    """Scale of the row-to-row differences, immune to the scars themselves."""
    diff = np.diff(data, axis=0)
    diff = diff[np.isfinite(diff)]
    if diff.size == 0:
        return 0.0
    return float(_MAD_TO_SIGMA * np.median(np.abs(diff - np.median(diff))))


def mark_scars(data, threshold: float = 4.0, min_length: int = 16,
               max_width: int = 4) -> np.ndarray:
    """
    Boolean mask of scar pixels, the same shape as ``data``.

    Parameters:
        data: 2-D height array, already levelled (see the module docstring)
        threshold: how far a scar must stand off both neighbours, in units of
            the robust row-to-row noise
        min_length: shortest run along a line that counts, in pixels
        max_width: thickest scar, in lines

    Returns an all-False mask when the noise scale is zero (a perfectly
    synthetic image), since no threshold can then be stated.
    """
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"expected a 2-D image, got shape {data.shape}")
    n_rows, n_cols = data.shape
    marks = np.zeros(data.shape, dtype=bool)

    sigma = _robust_sigma(data)
    if not sigma > 0:
        return marks
    high = threshold * sigma
    low = 0.5 * high

    # Per pixel, whether some band containing it stands off its neighbours,
    # at either threshold. Collected over every width first and only then
    # tested for length: two overlapping scars on adjacent lines change band
    # width part-way along, and a length test per width would cut the run at
    # each change and lose the short pieces between.
    strong = np.zeros(data.shape, dtype=bool)
    weak = np.zeros(data.shape, dtype=bool)
    for width in range(1, max_width + 1):
        if n_rows < width + 2:
            break
        above = data[:n_rows - width - 1]          # the line before each band
        below = data[width + 1:]                   # the line after each band
        # Every line of the band compared against both neighbours; the band's
        # weakest margin is the one that decides.
        band = np.stack([data[1 + j:n_rows - width + j] for j in range(width)])
        up = np.min(band, axis=0) - np.maximum(above, below)
        down = np.minimum(above, below) - np.max(band, axis=0)
        margin = np.maximum(up, down)
        # A NaN line (dropped earlier, or masked) cannot vote either way.
        margin = np.where(np.isfinite(margin), margin, -np.inf)
        for j in range(width):
            strong[1 + j:n_rows - width + j] |= margin > high
            weak[1 + j:n_rows - width + j] |= margin > low

    # The length rule applies to the strong core, and only a core that passes
    # it is then grown along the line at the low threshold. Testing the grown
    # run instead lets a ten-pixel spike seed a weak run the length of the
    # whole line, wherever a line sits a couple of sigma off its neighbours.
    for row in np.flatnonzero(strong.any(axis=1)):
        seeds = np.zeros(n_cols, dtype=bool)
        for a, b in zip(*_runs(strong[row])):
            if b - a >= min_length:
                seeds[a:b] = True
        if not seeds.any():
            continue
        for a, b in zip(*_runs(weak[row])):
            if seeds[a:b].any():
                marks[row, a:b] = True

    return marks


def repair_scars(data, method: str = 'interpolate', threshold: float = 4.0,
                 min_length: int = 16, max_width: int = 4):
    """
    Mark scars and repair them by ``method``.

    Returns ``(data, ScarReport)``. The input is never modified in place.
    """
    if method not in VALID_SCAR_REPAIRS:
        raise ValueError(f"unknown scar repair {method!r}. "
                         f"Available: {', '.join(VALID_SCAR_REPAIRS)}")
    data = np.array(data, dtype=float)
    if method == 'none':
        return data, ScarReport('none', 0, 0, ())

    marks = mark_scars(data, threshold=threshold, min_length=min_length,
                       max_width=max_width)
    lines = tuple(int(r) for r in np.flatnonzero(marks.any(axis=1)))
    report = ScarReport(method, len(lines), int(marks.sum()), lines)
    if not lines:
        return data, report

    if method == 'drop':
        data[list(lines), :] = np.nan
        return data, report

    rows = np.arange(data.shape[0])
    for col in np.flatnonzero(marks.any(axis=0)):
        bad = marks[:, col]
        good = ~bad & np.isfinite(data[:, col])
        if good.sum() < 2:
            continue  # nothing to interpolate from; leave the column alone
        data[bad, col] = np.interp(rows[bad], rows[good], data[good, col])
    return data, report


def repair_scars_for(data, settings):
    """Pipeline entry point: ``repair_scars`` with the settings' parameters."""
    return repair_scars(data,
                        method=getattr(settings, 'scar_repair', 'none'),
                        threshold=getattr(settings, 'scar_threshold', 4.0),
                        min_length=getattr(settings, 'scar_min_length', 16),
                        max_width=getattr(settings, 'scar_max_width', 4))
