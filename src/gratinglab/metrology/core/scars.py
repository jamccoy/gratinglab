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

#: How much stronger a neighbouring scar of the other sign must be before a
#: candidate is taken for its shadow. See ``_suppress_shadows``.
_SHADOW_RATIO = 2.0


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


def _level_by_steps(data):
    """
    Each line shifted so the median step from the line before it is zero.

    The detector's own levelling, used for marking only. ``align_rows`` takes
    each line's median height, and a grating's heights spread over its whole
    depth, so a scar along part of a line moves that median and leaves the rest
    of the line offset by tens of nanometres. On a rotated grating the lines
    next to such a line then look like scars too: 19 lines were marked for 3.
    The median *step* between neighbouring lines is a different statistic,
    immune to a scar covering under half the line, so after this only the scar
    stands off. Lines that are NaN throughout (masked, dropped) are stepped
    over, to the next line that has data.
    """
    out = np.array(data, dtype=float)
    previous = None
    for row in range(out.shape[0]):
        if not np.isfinite(out[row]).any():
            continue
        if previous is not None:
            # Against the previous line as already levelled, so the shifts
            # chain without being added up twice.
            step = out[row] - previous
            step = step[np.isfinite(step)]
            if step.size:
                out[row] -= float(np.median(step))
        previous = out[row]
    return out


def _relevel_scarred_lines(data, marks, min_length):
    """
    Remove a scarred line's offset against its clean neighbours, in place.

    Measured on the line's unmarked pixels against the straight line between
    the nearest wholly clean lines above and below. This is the offset the
    levelling put there (see ``_level_by_steps``); interpolating only the
    marked pixels would otherwise leave the rest of the line sitting tens of
    nanometres off.
    """
    clean = np.flatnonzero(~marks.any(axis=1) & np.isfinite(data).any(axis=1))
    for row in np.flatnonzero(marks.any(axis=1)):
        above = clean[clean < row]
        below = clean[clean > row]
        if above.size and below.size:
            a, b = above[-1], below[0]
            ref = data[a] + (data[b] - data[a]) * (row - a) / (b - a)
        elif above.size or below.size:
            ref = data[above[-1] if above.size else below[0]]
        else:
            continue
        cols = ~marks[row] & np.isfinite(data[row]) & np.isfinite(ref)
        if cols.sum() >= min_length:
            data[row] -= np.median(data[row, cols] - ref[cols])


def _suppress_shadows(best_up, best_down, low):
    """
    Drop candidates that only exist because a stronger scar is next to them.

    A band compared against a scar is judged against a neighbour that is not
    the surface. On a grating at an angle to the slow axis, neighbouring lines
    differ by several nanometres at each steep facet, so the lines on the far
    side of a raised scar sit below both of their neighbours there, and look
    like a lowered scar 21 pixels long at every groove: on a 3 degree grating
    with three scarred lines, eight more were marked. So, column by column, a
    vertical run of candidates of one sign that touches a run of the other
    sign reaching at least ``_SHADOW_RATIO`` times its margin is set aside.

    Twice, and not merely larger, because of the sandwich: two raised lines
    with an untouched one between them are, locally, the same picture as one
    lowered line, with the same margin, and only the rest of the lines can say
    which. Left alone, the length rule says it. A shadow is far weaker than
    the scar casting it (3 to 12 nm against 60 in the case above). A real
    raised scar beside a real lowered one under half its height would lose the
    lowered one; ringing of that shape is rare, and the stronger half is the
    one that matters.

    Modifies both arrays in place.
    """
    up = best_up > low
    down = best_down > low
    sign = np.where(up & (best_up >= best_down), 1,
                    np.where(down, -1, 0))
    size = np.where(sign > 0, best_up, np.where(sign < 0, best_down, 0.0))
    n_rows = sign.shape[0]
    for col in np.flatnonzero((sign != 0).any(axis=0)):
        s, m = sign[:, col], size[:, col]
        # Runs of one sign down this column, with the largest margin in each.
        runs = []
        r = 0
        while r < n_rows:
            if s[r] == 0:
                r += 1
                continue
            start = r
            while r < n_rows and s[r] == s[start]:
                r += 1
            runs.append((start, r, s[start], m[start:r].max()))
        for i, (start, stop, sgn, peak) in enumerate(runs):
            for j in (i - 1, i + 1):
                if not 0 <= j < len(runs):
                    continue
                o_start, o_stop, o_sgn, o_peak = runs[j]
                touching = o_stop == start or o_start == stop
                if touching and o_sgn == -sgn and o_peak > _SHADOW_RATIO * peak:
                    best_up[start:stop, col] = -np.inf
                    best_down[start:stop, col] = -np.inf
                    break


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

    # Judged on lines levelled against each other, not as they arrive.
    data = _level_by_steps(data)
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
    # The best margin each pixel reaches above and below, over every band
    # width it is part of. A NaN line (dropped earlier, or masked) cannot
    # vote either way.
    best_up = np.full(data.shape, -np.inf)
    best_down = np.full(data.shape, -np.inf)
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
        up = np.where(np.isfinite(up), up, -np.inf)
        down = np.where(np.isfinite(down), down, -np.inf)
        for j in range(width):
            rows = slice(1 + j, n_rows - width + j)
            np.maximum(best_up[rows], up, out=best_up[rows])
            np.maximum(best_down[rows], down, out=best_down[rows])

    _suppress_shadows(best_up, best_down, low)
    margin = np.maximum(best_up, best_down)
    strong = margin > high
    weak = margin > low

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

    _relevel_scarred_lines(data, marks, min_length)
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
