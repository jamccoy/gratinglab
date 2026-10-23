"""Groove shape statistics across swaths: how much the shape varies down a scan.

The boundary profile averages every row of a scan into one trace and every
groove along it into one shape. That answers "what does the groove look like"
but not "how much does it vary". This module answers the second question by
building the same averaged groove separately for each swath (the horizontal
row bands the blaze-angle analysis calls row groups, with exactly its banding)
and comparing them.

Two spreads are produced, and the result of interest is the comparison
between them:

``between_rms_nm``
    How much the averaged groove differs from swath to swath: variation down
    the image, along the length of the grooves.
``within_rms_nm``
    How much individual grooves differ from their swath's average: variation
    along the scan line, groove to groove.

A between-swath spread well above the within-swath one says the shape changes
along the grooves (a ruling or etch gradient, or a drifting scan); the reverse
says the grooves differ from each other more than any one groove changes along
its length.

Each swath's groove is aligned to the others by a circular shift before it is
compared, so a groove that is the same shape but displaced counts as the same
shape. Normalisation already rolls every groove so its minimum sits at the
boundary, so the alignment only removes the remaining offset. The
lateral position is measured separately, from where each swath finds its
grooves, and reported as ``shift_nm``: a steady drift down the image is a fact
about the scan (the grooves rotated relative to the fast axis, or thermal
drift) and should not be mistaken for a shape change.

Swaths are not independent samples: they re-measure the same grooves, so a
spread here describes variation along the grooves, not the uncertainty of a
mean. The blaze-angle row groups carry the same caveat and correct for it;
see the ``icc-correction`` wiki page. The ``swath-statistics`` page is the
reader's guide to everything here.
"""

from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core.processing import raw_data_multi_group
from .pipeline import BoundaryProfile, build_boundary_profile

__all__ = ["SwathStatistics", "swath_statistics"]

#: Per-swath scalars reported, from each swath's profile metrics.
METRICS = ("groove_depth", "peak_to_valley", "max_angle_deg", "rms_slope", "max_curvature")


@dataclass(frozen=True)
class SwathStatistics:
    """The averaged groove of each swath, aligned, and how much they differ."""

    #: The common axis, fractions of one period, and the whole-scan profile the
    #: swaths were aligned to.
    t: np.ndarray
    whole: BoundaryProfile
    #: One row per swath that produced a profile, aligned to one another.
    profiles: np.ndarray
    #: (first row, end row) of each of those swaths, in scan rows.
    row_ranges: list
    #: Lateral offset of each swath's groove from the whole-scan one, nm.
    shift_nm: np.ndarray
    #: Per-swath scalars and counts, one dict per swath.
    per_swath: list
    #: Swaths that produced no profile, as ((first, end), reason).
    skipped: list = field(default_factory=list)

    @property
    def n_swaths(self) -> int:
        return len(self.profiles)

    @property
    def mean(self) -> np.ndarray:
        return self.profiles.mean(axis=0)

    @property
    def std(self) -> np.ndarray:
        """Pointwise spread across swaths, fractions of the period."""
        return self.profiles.std(axis=0, ddof=1) if self.n_swaths > 1 else np.zeros_like(self.t)

    @property
    def between_rms_nm(self) -> float:
        """Swath-to-swath shape spread, RMS over the period, nm."""
        return float(np.sqrt(np.mean(self.std ** 2)) * self.whole.period_nm)

    @property
    def within_rms_nm(self) -> float:
        """Groove-to-groove spread within a swath, RMS over the period, averaged over swaths, nm."""
        values = [s["within_rms_nm"] for s in self.per_swath]
        return float(np.mean(values)) if values else float("nan")

    @property
    def whole_vs_swaths_rms_nm(self) -> float:
        """How far the whole-scan groove sits from the swaths' mean, RMS, nm.

        Large when the grooves drift sideways down the image: averaging every
        row then blurs the shape, so the whole-scan groove, the one exported,
        is rounder than any swath. Compare with ``shift_nm``.
        """
        mean = self.mean[:-1]
        whole, _ = _align(self.whole.y_norm[:-1], mean)
        return float(np.sqrt(np.mean((whole - mean) ** 2)) * self.whole.period_nm)

    def metric_summary(self) -> dict:
        """mean, std, min and max of each per-swath scalar."""
        out = {}
        for name in METRICS + ("period_nm",):
            values = np.array([s[name] for s in self.per_swath], dtype=float)
            out[name] = {
                "mean": float(values.mean()),
                "std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
                "min": float(values.min()),
                "max": float(values.max()),
            }
        return out

    def summary(self) -> str:
        """A few lines for a panel or a log."""
        m = self.metric_summary()
        lines = [
            f"{self.n_swaths} swaths" + (f" ({len(self.skipped)} skipped)" if self.skipped else ""),
            f"Shape spread between swaths : {self.between_rms_nm:.2f} nm RMS",
            f"Groove-to-groove within one : {self.within_rms_nm:.2f} nm RMS",
            f"Whole scan vs swath mean    : {self.whole_vs_swaths_rms_nm:.2f} nm RMS",
            f"Depth        : {m['groove_depth']['mean']:.4f} +/- {m['groove_depth']['std']:.4f} of period",
            f"Max sidewall : {m['max_angle_deg']['mean']:.2f} +/- {m['max_angle_deg']['std']:.2f} deg",
            f"Period       : {m['period_nm']['mean']:.2f} +/- {m['period_nm']['std']:.2f} nm",
            f"Lateral drift: {np.ptp(self.shift_nm):.1f} nm across the swaths",
        ]
        return "\n".join(lines)

    def write_csv(self, path) -> Path:
        """One row per swath: rows, shift, counts and scalars, under # comment lines."""
        path = Path(path)
        columns = ["first_row", "end_row", "shift_nm", "n_grooves_used", "period_nm",
                   "within_rms_nm", *METRICS]
        lines = [",".join(columns),
                 f"# between-swath shape spread: {self.between_rms_nm!r} nm RMS",
                 f"# within-swath groove-to-groove spread: {self.within_rms_nm!r} nm RMS",
                 f"# whole-scan groove vs swath mean: {self.whole_vs_swaths_rms_nm!r} nm RMS",
                 "# swaths re-measure the same grooves, so these describe variation, not independent samples"]
        for (first, end), shift, s in zip(self.row_ranges, self.shift_nm, self.per_swath):
            values = [first, end, float(shift), s["n_used"], s["period_nm"], s["within_rms_nm"],
                      *(s[name] for name in METRICS)]
            lines.append(",".join(repr(v) if isinstance(v, float) else str(v) for v in values))
        path.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
        return path


def swath_statistics(data, scan_x_size, settings, n_swaths=None) -> SwathStatistics:
    """Build the averaged groove for each swath of ``data`` and compare them.

    Parameters:
        data: 2-D height array in metres, image-flattened by the caller, the
            same input ``build_boundary_profile`` takes.
        scan_x_size: scan width in microns.
        settings: AnalysisSettings, applied to each swath exactly as to the
            whole scan (tip correction, edge exclusion, the ggp_* fields).
        n_swaths: how many bands; defaults to ``settings.n_row_groups``. Bands
            are at least 3 rows, so a short scan gets fewer, as with the
            blaze-angle row groups whose banding this reuses.

    Raises:
        ValueError if the whole scan gives no profile, or fewer than two swaths
        do, since there is then nothing to compare.
    """
    n = settings.n_row_groups if n_swaths is None else int(n_swaths)
    whole = build_boundary_profile(data, scan_x_size, settings)
    reference = whole.y_norm[:-1]  # closed profile: the last point repeats the first

    with contextlib.redirect_stdout(io.StringIO()):  # it narrates its banding
        _, _, info = raw_data_multi_group(data, scan_x_size, n_groups=n)

    whole_phase = _groove_phase(whole)
    profiles, ranges, phases, per_swath, skipped = [], [], [], [], []
    for first, end in info["group_ranges"]:
        try:
            swath = build_boundary_profile(data[first:end, :], scan_x_size, settings)
        except ValueError as exc:
            skipped.append(((first, end), str(exc)))
            continue
        if len(swath.y_norm) != len(whole.y_norm):
            skipped.append(((first, end), "profile length differs from the whole scan's"))
            continue
        aligned, _ = _align(swath.y_norm[:-1], reference)
        profiles.append(np.append(aligned, aligned[0]))
        ranges.append((first, end))
        phases.append(_groove_phase(swath) - whole_phase)
        per_swath.append({
            **{name: float(swath.metrics[name]) for name in METRICS},
            "period_nm": float(swath.period_nm),
            "n_used": int(swath.n_used),
            "within_rms_nm": float(np.sqrt(np.mean(np.asarray(swath.y_std_nm) ** 2))),
        })

    if len(profiles) < 2:
        raise ValueError(
            f"only {len(profiles)} of {len(info['group_ranges'])} swaths produced a "
            "groove profile, so there is nothing to compare; use fewer swaths")

    # The whole-scan groove is only a first reference: if the grooves drift
    # down the image, averaging every row blurs it, so it is rounder than any
    # swath. Re-align to the swaths' own mean, then roll the set so that mean
    # has its minimum at the boundary, the convention normalisation uses.
    closed = np.array(profiles)
    open_ = closed[:, :-1]
    for _ in range(2):
        target = open_.mean(axis=0)
        open_ = np.array([_align(y, target)[0] for y in open_])
    open_ = np.roll(open_, -int(np.argmin(open_.mean(axis=0))), axis=1)
    profiles = np.hstack([open_, open_[:, :1]])

    return SwathStatistics(
        t=whole.x_norm,
        whole=whole,
        profiles=profiles,
        row_ranges=ranges,
        # Unwrapped down the scan, so a drift past half a period reads as a
        # trend rather than jumping back.
        shift_nm=np.unwrap(np.array(phases)) / (2 * np.pi) * whole.period_nm,
        per_swath=per_swath,
        skipped=skipped,
    )


def _align(y: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, int]:
    """Circularly shift ``y`` onto ``reference``; return it and the signed shift.

    The shift minimising the squared difference is the one maximising the
    circular cross-correlation, found by FFT. The shift is signed, so a drift
    down the image reads as a trend instead of wrapping.
    """
    a = y - y.mean()
    b = reference - reference.mean()
    correlation = np.fft.irfft(np.fft.rfft(b) * np.conj(np.fft.rfft(a)), n=len(a))
    k = int(np.argmax(correlation))
    if k > len(a) // 2:
        k -= len(a)
    return np.roll(y, k), k


def _groove_phase(profile: BoundaryProfile) -> float:
    """Where the grooves sit along the scan line, as a phase of the period.

    The circular mean of the groove centres, so it needs no matching of one
    swath's grooves to another's and does not care which grooves the edge
    exclusion dropped.
    """
    x_nm = np.asarray(profile.profile_x_um)[np.asarray(profile.groove_centers)] * 1000.0
    return float(np.angle(np.mean(np.exp(2j * np.pi * x_nm / profile.period_nm))))
