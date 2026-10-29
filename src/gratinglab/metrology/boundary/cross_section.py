"""
From an SEM cross-section to a groove profile.

A cleaved cross-section shows the grating side-on: vacuum above, the grooves
in silhouette, solid material below. This module turns such an image, read by
:func:`~gratinglab.metrology.io.sem.read_sem`, into the two things the rest of
the package consumes:

* the averaged groove, normalised to the period, as a
  :class:`~gratinglab.profiles.FromProfileData`;
* a fitted :class:`~gratinglab.profiles.Trapezoidal` -- depth, crest and base
  widths, sidewall angle, corner radius -- with uncertainties and a statement
  of how well it fits.

The method, in the order it runs, is in ``docs/theory/metrology.md``,
section 6. Briefly:

1. **Edge.** Smooth lightly, set an iso-level halfway between the vacuum
   level (the top rows) and the material level (the bottom rows), and in each
   column take the *first* crossing of that level coming down from the vacuum
   side, interpolated linearly between the two pixels that straddle it. The
   first crossing is what makes a bright edge harmless: secondary electrons
   escaping near a surface make a rim brighter than the bulk, and the outer
   flank of that rim is the surface while its inner flank is not. A crossing
   must be followed by a short run of material, so an isolated charging speck
   in the vacuum is not taken for the surface. Where neighbouring columns
   cross more than a pixel apart -- a sidewall -- the rows between them are
   crossed horizontally too, so a steep wall is located to a fraction of a
   pixel across as well as along.
2. **Grooves.** The period comes from the strongest spatial frequency of the
   traced height, then from a straight-line fit to the ridge centres (the
   midpoint of each ridge's two half-height crossings), which also gives each
   groove's own centre and the scatter of their placement. Any tilt of the
   cleave edge is removed by fitting a line to the per-groove mean heights,
   which, unlike a least-squares line through the whole trace, is unbiased by
   the groove shape.
3. **Fit.** Every complete groove is folded into one period about its own
   centre, and a :class:`~gratinglab.profiles.Trapezoidal` is fitted to all the
   folded points at once by *orthogonal* distance, so a steep wall counts as
   much as a flat. The same fit run groove by groove gives the uncertainties:
   the groove-to-groove standard error, which carries the real variation of the
   etch and not only the pixel noise.
4. **Average.** The averaged groove is the fitted shape plus the mean of the
   orthogonal residuals, binned along its arc length. It therefore keeps
   whatever the trapezoid does not describe -- a bowed floor, unequal crest
   and floor radii -- while the noise averages down.

Three things are inputs rather than guesses, because the file cannot say:

``foreshortened``
    which image axis a tilted stage compresses. :attr:`SemImage.foreshortening`
    is reported by the reader and not applied; here it must be applied, and a
    cross-section from a tilted stage with the correction off is refused until
    the caller names the axis.
``material_below``
    which side of the cleave edge is solid.
``image_x_along_dispersion``
    which way the image's columns run against the grating's dispersion
    direction. Irrelevant to a symmetric groove, and it decides the sign of
    ``t`` (``docs/conventions.md``, "Direction of the profile parameter").

What it cannot do is listed in the theory page: undercut walls (the trace is
a height function), a protective cap deposited over the grooves (it traces
whatever the outermost surface is), and an image that holds the grating's top
surface as well as the cleave face.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import ndimage, optimize
from scipy.spatial import cKDTree

from ...profiles import FromProfileData, Trapezoidal
from ..io.sem import SemImage

__all__ = ["CrossSection", "TrapezoidFit", "trace_cross_section"]

Foreshortened = Literal["vertical", "horizontal", "neither"]

#: Fraction of the region's rows read as pure vacuum (top) and pure material
#: (bottom) to set the two levels the edge is placed between.
_BAND_FRACTION = 0.06

#: Smallest material-minus-vacuum contrast accepted, in units of the vacuum
#: noise. Below it the iso-level sits in the noise and the trace is noise.
_MIN_CONTRAST_SIGMA = 4.0

#: Hysteresis for finding a ridge: a crossing of half height counts only once
#: the trace has gone from below this fraction of the depth to above one minus
#: it. Stops noise on a wall registering as two walls.
_HYSTERESIS = 0.25

#: Points on the model curve, per period, for the orthogonal distances.
_MODEL_POINTS = 900

#: Most passes of the averaging; it stops sooner once a pass moves nothing
#: by more than a hundredth of a pixel.
_AVERAGE_PASSES = 6

#: Parameter bounds for the fit, as (depth/P, crest width, base-width share,
#: radius share, horizontal offset, vertical offset). See `_trapezoid`.
_LOWER = np.array([1e-3, 0.02, 0.0, 0.0, -0.25, -np.inf])
_UPPER = np.array([3.0, 0.97, 0.97, 0.95, 0.25, np.inf])


# ── Results ─────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class TrapezoidFit:
    """A :class:`Trapezoidal` fitted to a traced cross-section.

    Widths are fractions of the period, as :class:`Trapezoidal` holds them;
    the trench widths a micrograph is usually read for are one minus the ridge
    widths. Every ``*_err`` is one standard error: the groove-to-groove spread
    divided by the square root of the number of grooves, unless only one groove
    was complete, in which case it is the least-squares covariance and says so
    in :attr:`uncertainty`.
    """

    profile: Trapezoidal
    period_nm: float
    depth_nm: float
    depth_nm_err: float
    ridge_top_width: float
    ridge_top_width_err: float
    ridge_base_width: float
    ridge_base_width_err: float
    sidewall_angle_deg: float
    sidewall_angle_deg_err: float
    corner_radius_nm: float
    corner_radius_nm_err: float
    #: Orthogonal RMS distance of every traced point from the fitted shape.
    #: Pixel noise and groove-to-groove variation, plus any misfit.
    rms_residual_nm: float
    #: RMS of the averaged groove about the fitted shape. The noise has
    #: averaged down here, so what remains is mostly shape the trapezoid
    #: does not describe.
    mean_profile_rms_nm: float
    #: Largest excursion of the averaged groove from the fitted shape.
    mean_profile_max_nm: float
    n_grooves: int
    n_points: int
    uncertainty: Literal["groove-to-groove", "covariance"]
    #: The same five quantities fitted groove by groove, for a histogram or a
    #: drift plot: keys ``depth_nm``, ``ridge_top_width``, ``ridge_base_width``,
    #: ``sidewall_angle_deg``, ``corner_radius_nm``.
    per_groove: dict[str, NDArray[np.float64]]

    @property
    def trench_width_mouth(self) -> float:
        """Width of the trench at the crest, as a fraction of the period."""
        return 1.0 - self.ridge_top_width

    @property
    def trench_width_floor(self) -> float:
        """Width of the trench where the walls meet the floor."""
        return 1.0 - self.ridge_base_width

    @property
    def summary(self) -> str:
        return (
            f"depth {self.depth_nm:.1f} +/- {self.depth_nm_err:.1f} nm, "
            f"trench {self.trench_width_mouth:.2f} -> {self.trench_width_floor:.2f} "
            f"of the period, sidewall {self.sidewall_angle_deg:.1f} +/- "
            f"{self.sidewall_angle_deg_err:.1f} deg, corner radius "
            f"{self.corner_radius_nm:.1f} +/- {self.corner_radius_nm_err:.1f} nm; "
            f"{self.n_grooves} grooves, residual {self.rms_residual_nm:.2f} nm rms, "
            f"averaged groove {self.mean_profile_rms_nm:.2f} nm rms from the trapezoid"
        )


@dataclass(frozen=True)
class CrossSection:
    """A groove boundary traced off a cleave edge, averaged and fitted."""

    #: The traced boundary across the whole region, nm, levelled. ``x_nm`` runs
    #: along the image columns and ``height_nm`` upward, from an arbitrary zero.
    x_nm: NDArray[np.float64]
    height_nm: NDArray[np.float64]
    #: The period the profile is normalised to.
    period_nm: float
    #: Where it came from: ``"measured"`` from the ridge lattice, or
    #: ``"given"`` when the caller supplied it.
    period_source: Literal["measured", "given"]
    #: Period from the ridge lattice, and its standard error. Reported even
    #: when a given period is the one used, so the two can be compared.
    measured_period_nm: float
    measured_period_err_nm: float
    #: Centre of every complete ridge, nm, along ``x_nm``.
    groove_centres_nm: NDArray[np.float64]
    #: RMS departure of those centres from a perfect lattice.
    placement_rms_nm: float
    #: The averaged groove, normalised: ``u`` along the image columns over one
    #: period with the ridge centred on 0.5, ``y`` in periods with its minimum
    #: at zero. :meth:`to_profile` turns it into ``t``.
    u: NDArray[np.float64]
    y: NDArray[np.float64]
    #: Standard error of each averaged point, normal to the groove, in periods.
    normal_err: NDArray[np.float64]
    fit: TrapezoidFit
    #: The iso-level the edge was placed at, and the two levels it lies
    #: between, in detector units of the image as given.
    level: float
    vacuum_level: float
    material_level: float
    #: Lengths per pixel actually used, after any foreshortening, nm.
    pixel_x_nm: float
    pixel_y_nm: float
    #: Which axis was stretched by ``1 / foreshortening``, and by how much.
    foreshortened: Foreshortened
    foreshortening: float
    image_x_along_dispersion: bool

    @property
    def depth_nm(self) -> float:
        """Peak-to-valley depth of the averaged groove."""
        return float(np.ptp(self.y) * self.period_nm)

    def to_profile(self) -> FromProfileData:
        """The averaged groove in the package's conventions.

        ``t`` runs against the dispersion direction. With the default
        ``image_x_along_dispersion=True`` the image's columns are taken to run
        along it, as the AFM pipeline's ``x_norm`` does, and the conversion is
        :meth:`FromProfileData.along_dispersion`.
        """
        u = np.clip(self.u, 0.0, 1.0)
        if self.image_x_along_dispersion:
            return FromProfileData.along_dispersion(u, self.y)
        return FromProfileData(t=tuple(u), y=tuple(self.y))

    @property
    def summary(self) -> str:
        return (
            f"period {self.period_nm:.1f} nm ({self.period_source}), "
            f"{self.fit.n_grooves} grooves; {self.fit.summary}"
        )


# ── The edge ────────────────────────────────────────────────────────────────


def _levels(
    smooth: NDArray[np.float64], level: float | None, level_fraction: float
) -> tuple[float, float, float, bool]:
    """Vacuum level, material level, iso-level and whether contrast is inverted.

    Works on the image with material at the bottom. Returns the levels in the
    image's own units, and ``True`` when material is the darker of the two, in
    which case the caller negates the image so material is always brighter.
    """
    rows = smooth.shape[0]
    band = max(3, int(round(rows * _BAND_FRACTION)))
    top, bottom = smooth[:band], smooth[-band:]
    vacuum, material = float(np.median(top)), float(np.median(bottom))
    inverted = material < vacuum
    noise = float(np.std(top))
    if level is None:
        if abs(material - vacuum) < _MIN_CONTRAST_SIGMA * max(noise, 1e-12):
            raise ValueError(
                f"material ({material:.1f}) and vacuum ({vacuum:.1f}) differ by "
                f"less than {_MIN_CONTRAST_SIGMA:g} times the vacuum noise "
                f"({noise:.2f}), so a level between them would trace noise. "
                "Crop the region so its top rows are vacuum and its bottom rows "
                "are solid material, or pass `level` explicitly"
            )
        if not 0.0 < level_fraction < 1.0:
            raise ValueError(f"level_fraction must lie in (0, 1), got {level_fraction}")
        level = vacuum + level_fraction * (material - vacuum)
    elif not min(vacuum, material) < level < max(vacuum, material):
        raise ValueError(
            f"level {level} does not lie between the vacuum ({vacuum:.1f}) and "
            f"material ({material:.1f}) levels of the region"
        )
    return vacuum, material, float(level), inverted


def _trace(
    work: NDArray[np.float64], level: float, material: float, run: int
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """The boundary as an ordered polyline, in pixel-centre coordinates.

    ``work`` has material at the bottom and brighter than vacuum. Returns
    ``(x, y)`` of the polyline, with ``y`` the row (increasing downward), and
    the per-column crossing row (NaN where a column has none).
    """
    rows, cols = work.shape
    above = work >= level
    # A crossing must be followed by `run` rows of material.
    ok = above.copy()
    for k in range(1, run):
        ok[: rows - k] &= above[k:]
    first = np.argmax(ok, axis=0)
    valid = ok.any(axis=0) & (first > 0)

    col_y = np.full(cols, np.nan)
    j = np.flatnonzero(valid)
    i = first[j]
    a, b = work[i - 1, j], work[i, j]
    frac = np.clip((level - a) / np.where(b > a, b - a, 1.0), 0.0, 1.0)
    col_y[j] = i - 1 + frac

    # Below its crossing every pixel is material: the bulk may be darker than
    # a bright rim, and a row crossing at a wall needs the material side on the
    # right side of the level. Only a pixel that reads below the level is
    # raised, and only to the material level. Raising every pixel below the
    # crossing would also lift the partly covered pixel at the wall itself,
    # and that alone moved each wall 0.17 nm into the vacuum on a clean image.
    below = np.arange(rows)[:, None] >= np.where(valid, first, rows)[None, :]
    filled = np.where(below & (work < level), material, work)

    xs: list[float] = []
    ys: list[float] = []
    for c in range(cols):
        if not valid[c]:
            continue
        xs.append(float(c))
        ys.append(float(col_y[c]))
        if c + 1 >= cols or not valid[c + 1]:
            continue
        y0, y1 = col_y[c], col_y[c + 1]
        lo, hi = min(y0, y1), max(y0, y1)
        between = np.arange(int(np.floor(lo)) + 1, int(np.ceil(hi)))
        between = between[(between > lo) & (between < hi)]
        if between.size == 0:
            continue
        if y1 < y0:
            between = between[::-1]  # rising wall: walk upward
        va, vb = filled[between, c], filled[between, c + 1]
        step = vb - va
        f = np.clip(
            np.divide(level - va, step, out=np.full_like(va, 0.5), where=step != 0),
            0.0, 1.0,
        )
        xs.extend((c + f).tolist())
        ys.extend(between.astype(np.float64).tolist())
    return np.asarray(xs), np.asarray(ys), col_y


# ── The grooves ─────────────────────────────────────────────────────────────


def _fft_period(x: NDArray[np.float64], h: NDArray[np.float64]) -> float:
    """Period of a uniformly sampled trace from its strongest frequency."""
    dx = x[1] - x[0]
    n = len(h)
    detrended = h - np.polyval(np.polyfit(x, h, 1), x)
    pad = 16 * n
    spectrum = np.abs(np.fft.rfft((detrended - detrended.mean()) * np.hanning(n), pad))
    freq = np.fft.rfftfreq(pad, d=dx)
    usable = freq >= 2.0 / (n * dx)  # at least two periods in the region
    if not usable.any():
        raise ValueError("the region is too narrow to hold two periods")
    k = int(np.flatnonzero(usable)[np.argmax(spectrum[usable])])
    if 0 < k < len(spectrum) - 1:
        a, b, c = spectrum[k - 1], spectrum[k], spectrum[k + 1]
        denom = a - 2 * b + c
        shift = 0.5 * (a - c) / denom if denom != 0 else 0.0
    else:
        shift = 0.0
    return float(1.0 / ((k + shift) * (freq[1] - freq[0])))


def _ridge_centres(
    x: NDArray[np.float64], h: NDArray[np.float64], floor: float, crest: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Centre and half-height width of every ridge the trace crosses fully."""
    depth = crest - floor
    lo, hi = floor + _HYSTERESIS * depth, crest - _HYSTERESIS * depth
    mid = 0.5 * (floor + crest)
    state: str | None = None
    last: float | None = None
    events: list[tuple[int, float]] = []
    for k in range(len(h) - 1):
        a, b = h[k], h[k + 1]
        if (a - mid) * (b - mid) <= 0 and a != b:
            last = x[k] + (mid - a) / (b - a) * (x[k + 1] - x[k])
        if state != "high" and b >= hi:
            if state == "low" and last is not None:
                events.append((+1, last))
            state, last = "high", None
        elif state != "low" and b <= lo:
            if state == "high" and last is not None:
                events.append((-1, last))
            state, last = "low", None
    centres, widths = [], []
    for (s0, x0), (s1, x1) in zip(events[:-1], events[1:]):
        if s0 == +1 and s1 == -1:
            centres.append(0.5 * (x0 + x1))
            widths.append(x1 - x0)
    return np.asarray(centres), np.asarray(widths)


def _lattice(
    centres: NDArray[np.float64], period: float
) -> tuple[float, float, float, NDArray[np.float64]]:
    """Fit a lattice to ridge centres: period, its error, placement RMS."""
    n = np.round((centres - centres[0]) / period).astype(int)
    if len(np.unique(n)) != len(n):
        raise ValueError(
            "two ridges fall in one period of the lattice; the period estimate "
            "is wrong for this image. Pass period_nm"
        )
    if len(centres) < 2:
        return period, float("nan"), 0.0, centres
    coef, cov = np.polyfit(n, centres, 1, cov=len(centres) > 3)
    if len(centres) > 3:
        err = float(np.sqrt(cov[0, 0]))
    else:
        err = float("nan")
    resid = centres - np.polyval(coef, n)
    rms = float(np.sqrt(np.mean(resid**2)))
    return float(coef[0]), err, rms, centres


def _window_mean(
    xc: NDArray[np.float64], hc: NDArray[np.float64], centre: float, period: float
) -> float:
    grid = np.linspace(centre - period / 2, centre + period / 2, 257)
    return float(np.mean(np.interp(grid, xc, hc)[:-1]))


# ── The fit ─────────────────────────────────────────────────────────────────


def _trapezoid(p: NDArray[np.float64]) -> Trapezoidal:
    """The profile for fit parameters ``p``.

    Parametrised so every point inside the bounds is a valid groove: the crest
    width ``w_t``, then the base width as a share ``s`` of the room left,
    ``w_b = w_t + s (1 - w_t)``, so the floor never closes, and the radius as a
    share of the largest that fits. ``s = 0`` is a vertical wall.
    """
    depth, top, share, rho = (float(v) for v in p[:4])
    base = top + share * (1.0 - top)
    run = 0.5 * (base - top)
    angle = 90.0 if run <= 1e-12 else float(np.degrees(np.arctan2(depth, run)))
    sharp = Trapezoidal(depth_fraction=depth, duty_cycle=top, sidewall_angle=angle)
    radius = rho * sharp.max_corner_radius * (1.0 - 1e-9)
    if radius <= 0.0:
        return sharp
    return Trapezoidal(
        depth_fraction=depth, duty_cycle=top, sidewall_angle=angle,
        corner_radius=radius,
    )


def _model_curve(
    p: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The fitted curve over three periods, offset, as ordered points."""
    curve = _trapezoid(p).boundary(_MODEL_POINTS)
    t = np.concatenate([curve.t - 1.0, curve.t, curve.t + 1.0, [2.0 + curve.t[0]]])
    y = np.concatenate([curve.y, curve.y, curve.y, [curve.y[0]]])
    return t + p[4], y + p[5]


def _project(
    cx: NDArray[np.float64], cy: NDArray[np.float64],
    px: NDArray[np.float64], py: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Signed orthogonal distance of points to a polyline, and the arc position.

    Positive on the vacuum side. The polyline advances in ``+u`` with material
    below and to its right, so the vacuum is always on its left.
    """
    seg = np.diff(np.column_stack([cx, cy]), axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    arc = np.concatenate([[0.0], np.cumsum(seglen)])
    _, nearest = cKDTree(np.column_stack([cx, cy])).query(np.column_stack([px, py]))

    best_d = np.full(px.shape, np.inf)
    best_signed = np.zeros_like(px)
    best_s = np.zeros_like(px)
    for offset in (-1, 0):
        k = np.clip(nearest + offset, 0, len(seg) - 1)
        ax, ay = cx[k], cy[k]
        dx, dy = seg[k, 0], seg[k, 1]
        L2 = np.maximum(seglen[k] ** 2, 1e-30)
        tau = np.clip(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0)
        qx, qy = ax + tau * dx, ay + tau * dy
        rx, ry = px - qx, py - qy
        d = np.hypot(rx, ry)
        cross = dx * ry - dy * rx
        better = d < best_d
        best_d = np.where(better, d, best_d)
        best_signed = np.where(better, np.where(cross >= 0, d, -d), best_signed)
        best_s = np.where(better, arc[k] + tau * seglen[k], best_s)
    return best_signed, best_s


def _residuals(p, u, y):
    cx, cy = _model_curve(p)
    return _project(cx, cy, u, y)[0]


def _fit(p0, u, y):
    p0 = np.clip(p0, _LOWER + 1e-9, np.where(np.isfinite(_UPPER), _UPPER - 1e-9, _UPPER))
    return optimize.least_squares(
        _residuals, p0, args=(u, y), bounds=(_LOWER, _UPPER),
        diff_step=1e-3, x_scale="jac", max_nfev=400,
    )


def _derived(p: NDArray[np.float64], period: float) -> NDArray[np.float64]:
    """(depth nm, crest width, base width, sidewall angle, radius nm)."""
    shape = _trapezoid(p)
    return np.array([
        shape.depth_fraction * period,
        shape.ridge_top_width,
        shape.ridge_base_width,
        shape.sidewall_angle,
        shape.corner_radius * period,
    ])


_DERIVED_KEYS = (
    "depth_nm", "ridge_top_width", "ridge_base_width",
    "sidewall_angle_deg", "corner_radius_nm",
)


def _tile(cu, cy):
    """A closed one-period curve repeated over three periods, for projection."""
    body_u, body_y = cu[:-1], cy[:-1]
    return (
        np.concatenate([body_u - 1.0, body_u, body_u + 1.0, [cu[-1] + 1.0]]),
        np.concatenate([body_y, body_y, body_y, [cy[-1]]]),
    )


def _binned_average(
    cu: NDArray[np.float64], cy: NDArray[np.float64],
    u: NDArray[np.float64], y: NDArray[np.float64], n_bins: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64],
           NDArray[np.float64]]:
    """Move a reference curve by the mean residual of the points about it.

    ``(cu, cy)`` is closed over one period, last point one period after the
    first. The points are projected onto it orthogonally, binned by the arc
    position of their foot, and each bin's centre is moved along the normal
    by the bin's mean signed distance. Returns the moved bin centres, the
    move, and the standard error of each bin's mean (NaN for a bin of one).
    """
    seg = np.hypot(np.diff(cu), np.diff(cy))
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    total = arc[-1]
    signed, s = _project(*_tile(cu, cy), u, y)
    s = np.mod(s - total, total)  # the first tile is one period long

    edges = np.linspace(0.0, total, n_bins + 1)
    which = np.clip(np.searchsorted(edges, s, side="right") - 1, 0, n_bins - 1)
    count = np.bincount(which, minlength=n_bins)
    sum_r = np.bincount(which, weights=signed, minlength=n_bins)
    sum_r2 = np.bincount(which, weights=signed**2, minlength=n_bins)
    centre = 0.5 * (edges[:-1] + edges[1:])
    filled = count > 0
    mean_r = np.zeros(n_bins)
    mean_r[filled] = sum_r[filled] / count[filled]
    mean_r = np.interp(centre, centre[filled], mean_r[filled], period=total)
    err = np.full(n_bins, np.nan)
    many = count > 1
    var = (sum_r2[many] - count[many] * mean_r[many] ** 2) / (count[many] - 1)
    err[many] = np.sqrt(np.maximum(var, 0.0) / count[many])

    bu, by = np.interp(centre, arc, cu), np.interp(centre, arc, cy)
    # Tangent by central differences, wrapping round the period.
    wu = np.concatenate([[bu[-1] - 1.0], bu, [bu[0] + 1.0]])
    wy = np.concatenate([[by[-1]], by, [by[0]]])
    tu, ty = wu[2:] - wu[:-2], wy[2:] - wy[:-2]
    norm = np.maximum(np.hypot(tu, ty), 1e-30)
    nu, ny = -ty / norm, tu / norm  # left of travel: the vacuum side
    return bu + mean_r * nu, by + mean_r * ny, mean_r, err


# ── The whole step ──────────────────────────────────────────────────────────


def trace_cross_section(
    image: SemImage | ArrayLike,
    *,
    pixel_size_nm: float | None = None,
    foreshortened: Foreshortened | None = None,
    foreshortening: float | None = None,
    region: tuple[int, int, int, int] | None = None,
    material_below: bool = True,
    period_nm: float | None = None,
    smoothing_px: float = 0.5,
    level: float | None = None,
    level_fraction: float = 0.5,
    image_x_along_dispersion: bool = True,
) -> CrossSection:
    """Trace the groove boundary along a cleave edge and fit a trapezoid to it.

    Parameters
    ----------
    image
        A :class:`SemImage` from :func:`read_sem`, or a 2-D array of detector
        levels, in which case ``pixel_size_nm`` is required.
    pixel_size_nm
        Overrides the image's own pixel size. Required for an array.
    foreshortened
        Which image axis a tilted stage compressed: ``"vertical"`` (heights),
        ``"horizontal"`` (along the grooves' period) or ``"neither"``. Required
        whenever the foreshortening factor is below one, because the file does
        not record which; the factor is then undone on that axis alone.
    foreshortening
        Overrides :attr:`SemImage.foreshortening`, the ``cos`` of the stage
        tilt. Defaults to the image's own, or 1 for an array.
    region
        ``(row0, row1, col0, col1)`` to crop to, as array slices. The region's
        top rows must be vacuum and its bottom rows solid material, because
        the two levels are read from them.
    material_below
        Whether the solid is below the cleave edge in the image.
    period_nm
        A period known better than a few dozen grooves can give it, from
        diffraction for instance. Used to normalise; the image's own lattice
        period is still measured and reported beside it.
    smoothing_px
        Gaussian smoothing before the edge is located, in pixels. Zero to
        disable. It does not move a straight edge, but it rounds every corner
        exactly as the beam does and spreads a bright rim into the vacuum, so
        it is kept small: at 1 px the fitted radius of a 6 nm fillet read 3 nm
        high, at 0.5 px 1.6 nm, while the scatter it saves is negligible.
    level, level_fraction
        The iso-level, in the image's detector units, or its place between the
        vacuum and material levels (0.5, halfway, by default).
    image_x_along_dispersion
        Whether the image's columns advance along the dispersion direction.
        Decides the sign of ``t`` in :meth:`CrossSection.to_profile`.

    Raises
    ------
    ValueError
        When the region holds too little contrast, fewer than one complete
        groove, or a tilted-stage image whose foreshortened axis is unnamed.
    """
    if isinstance(image, SemImage):
        data = np.asarray(image.data, dtype=np.float64)
        px = pixel_size_nm if pixel_size_nm is not None else image.pixel_size_nm
        factor = image.foreshortening if foreshortening is None else foreshortening
    else:
        data = np.asarray(image, dtype=np.float64)
        if pixel_size_nm is None:
            raise ValueError("pixel_size_nm is required when image is an array")
        px = pixel_size_nm
        factor = 1.0 if foreshortening is None else foreshortening
    if data.ndim != 2:
        raise ValueError(f"expected a 2-D image, got shape {data.shape}")
    if not px or px <= 0:
        raise ValueError(f"pixel_size_nm must be positive, got {px}")
    if not 0.0 < factor <= 1.0:
        raise ValueError(f"foreshortening must lie in (0, 1], got {factor}")

    if factor < 1.0 - 1e-9 and foreshortened is None:
        raise ValueError(
            f"this image is foreshortened by {factor:.4f} (a tilted stage with "
            "the tilt correction off), and the file does not say which image "
            "axis lies across the tilt. Pass foreshortened='vertical' if the "
            "heights are compressed, 'horizontal' if the period is, or "
            "'neither' if the cleave face was square to the beam"
        )
    if foreshortened not in (None, "vertical", "horizontal", "neither"):
        raise ValueError(f"foreshortened must be vertical, horizontal or neither, "
                         f"got {foreshortened!r}")
    axis: Foreshortened = foreshortened or "neither"
    sx = px / factor if axis == "horizontal" else px
    sy = px / factor if axis == "vertical" else px

    if region is not None:
        r0, r1, c0, c1 = region
        data = data[r0:r1, c0:c1]
    if not material_below:
        data = data[::-1]
    if min(data.shape) < 16:
        raise ValueError(f"region {data.shape} is too small to trace")

    smooth = ndimage.gaussian_filter(data, smoothing_px, mode="nearest") \
        if smoothing_px > 0 else data
    vacuum, material, iso, inverted = _levels(smooth, level, level_fraction)
    sign = -1.0 if inverted else 1.0
    run = max(2, int(np.ceil(2 * smoothing_px)))
    xp, yp, col_y = _trace(sign * smooth, sign * iso, sign * material, run)

    good = np.isfinite(col_y)
    if good.sum() < 0.5 * len(col_y):
        raise ValueError(
            f"only {good.sum()} of {len(col_y)} columns cross from vacuum into "
            "material; check that the region's top is vacuum, its bottom is "
            "solid, and material_below is right"
        )

    x_nm = xp * sx
    h_nm = -yp * sy
    xc = np.flatnonzero(good) * sx
    hc = -col_y[good] * sy

    # Level crudely, find the period, then level properly on whole grooves.
    tilt = np.polyfit(xc, hc, 1)
    period0 = period_nm if period_nm is not None else _fft_period(
        np.arange(len(col_y)) * sx,
        np.interp(np.arange(len(col_y)) * sx, xc, hc),
    )
    centres = np.array([])
    period_meas = period0
    for _ in range(3):
        hcl = hc - np.polyval(tilt, xc)
        floor, crest = np.percentile(hcl, [10, 90])
        centres, _ = _ridge_centres(x_nm, h_nm - np.polyval(tilt, x_nm), floor, crest)
        if centres.size == 0:
            break
        period_meas, period_err, placement, centres = _lattice(centres, period0)
        if period_nm is None:
            period0 = period_meas
        whole = centres[
            (centres - period0 / 2 >= xc[0]) & (centres + period0 / 2 <= xc[-1])
        ]
        if whole.size < 2:
            break
        means = [_window_mean(xc, hcl, c, period0) for c in whole]
        extra = np.polyfit(whole, means, 1)
        tilt = tilt + np.array([extra[0], 0.0])
    if centres.size == 0:
        raise ValueError("no complete ridge was found in the region")
    period = period_nm if period_nm is not None else period_meas

    h_nm = h_nm - np.polyval(tilt, x_nm)
    hcl = hc - np.polyval(tilt, xc)
    floor, crest = np.percentile(hcl, [10, 90])
    centres, widths = _ridge_centres(x_nm, h_nm, floor, crest)
    if centres.size == 0:
        raise ValueError("no complete ridge was found in the region")
    period_meas, period_err, placement, centres = _lattice(centres, period0)
    if period_nm is None:
        period = period_meas

    # Fold each complete groove into one period about its own centre.
    gap = 1.5 * sx
    grooves = []
    for c in centres:
        lo, hi = c - period / 2, c + period / 2
        if lo < xc[0] or hi > xc[-1]:
            continue
        inside = (xc >= lo) & (xc <= hi)
        if np.any(np.diff(xc[inside]) > gap):
            continue  # a column without a crossing breaks this groove
        sel = (x_nm >= lo) & (x_nm < hi)
        grooves.append(((x_nm[sel] - c) / period + 0.5, h_nm[sel] / period, c))
    if not grooves:
        raise ValueError(
            "no groove lies wholly inside the region; widen it to hold at "
            "least one full period around a ridge"
        )

    u = np.concatenate([g[0] for g in grooves])
    y = np.concatenate([g[1] for g in grooves])

    mid_width = float(np.median(widths)) / period if widths.size else 0.5
    depth0 = (crest - floor) / period
    run0 = depth0 / np.tan(np.radians(75.0))
    top0 = float(np.clip(mid_width - run0, 0.05, 0.9))
    base0 = float(np.clip(mid_width + run0, top0 + 1e-3, 0.95))
    p0 = np.array([depth0, top0, (base0 - top0) / (1 - top0), 0.2, 0.0, floor / period])
    pooled = _fit(p0, u, y)
    p = pooled.x

    derived = _derived(p, period)
    per = np.array([_derived(_fit(p, gu, gy).x, period) for gu, gy, _ in grooves]) \
        if len(grooves) > 1 else derived[None, :]

    if len(grooves) > 1:
        err = per.std(axis=0, ddof=1) / np.sqrt(len(grooves))
        method: Literal["groove-to-groove", "covariance"] = "groove-to-groove"
    else:
        J = pooled.jac
        dof = max(len(u) - len(p), 1)
        s2 = float(np.sum(pooled.fun**2) / dof)
        cov = np.linalg.pinv(J.T @ J) * s2
        step = 1e-5 * np.maximum(np.abs(p), 1e-3)
        G = np.empty((len(derived), len(p)))
        for k in range(len(p)):
            q = p.copy()
            q[k] = np.clip(q[k] + step[k], _LOWER[k], _UPPER[k])
            G[:, k] = (_derived(q, period) - derived) / max(q[k] - p[k], 1e-15)
        err = np.sqrt(np.maximum(np.diag(G @ cov @ G.T), 0.0))
        method = "covariance"

    # Lengths in nm also carry the period's own uncertainty.
    rel_p = 0.0 if period_nm is not None or not np.isfinite(period_err) \
        else period_err / period
    for k in (0, 4):
        err[k] = float(np.hypot(err[k], derived[k] * rel_p))

    # The averaged groove: start from the fitted shape, bin the orthogonal
    # residuals along it, move it by their mean, and repeat on the result.
    # One pass is exact only where the trapezoid already fits; a groove it
    # does not describe needs the reference to follow the data, and three
    # more passes bring an asymmetric groove's steep wall back to its width.
    shape = _trapezoid(p)
    curve = shape.boundary(_MODEL_POINTS)
    model_u = np.concatenate([curve.t, [curve.t[0] + 1.0]])
    model_y = np.concatenate([curve.y, [curve.y[0]]])
    fu, fy = np.mod(u - p[4], 1.0), y - p[5]
    ref_u, ref_y = model_u, model_y
    arc_total = float(np.sum(np.hypot(np.diff(model_u), np.diff(model_y))))
    n_bins = max(64, int(round(arc_total * period / min(sx, sy))))
    for _ in range(_AVERAGE_PASSES):
        u_avg, y_avg, step, normal_err = _binned_average(ref_u, ref_y, fu, fy, n_bins)
        ref_u = np.concatenate([u_avg, [u_avg[0] + 1.0]])
        ref_y = np.concatenate([y_avg, [y_avg[0]]])
        if np.max(np.abs(step)) * period < 0.01 * min(sx, sy):
            break

    # How far the averaged groove sits from the fitted shape: the misfit,
    # with the noise averaged down.
    misfit, _ = _project(*_tile(model_u, model_y), u_avg, y_avg)

    # The trace is a height function and cannot see an undercut, so a point
    # that steps back in u can only be noise in the average on a steep wall.
    # Left in, it would make FromProfileData read the groove as undercut and
    # refuse height(); clamped to the point before it, it would make a
    # vertical jump. It is dropped.
    keep = np.concatenate([[True], u_avg[1:] > np.maximum.accumulate(u_avg)[:-1]])
    u_avg, y_avg, normal_err = u_avg[keep], y_avg[keep], normal_err[keep]
    y_avg = y_avg - y_avg.min()

    fit = TrapezoidFit(
        profile=shape,
        period_nm=float(period),
        depth_nm=float(derived[0]), depth_nm_err=float(err[0]),
        ridge_top_width=float(derived[1]), ridge_top_width_err=float(err[1]),
        ridge_base_width=float(derived[2]), ridge_base_width_err=float(err[2]),
        sidewall_angle_deg=float(derived[3]), sidewall_angle_deg_err=float(err[3]),
        corner_radius_nm=float(derived[4]), corner_radius_nm_err=float(err[4]),
        rms_residual_nm=float(np.sqrt(np.mean(pooled.fun**2)) * period),
        mean_profile_rms_nm=float(np.sqrt(np.mean(misfit**2)) * period),
        mean_profile_max_nm=float(np.max(np.abs(misfit)) * period),
        n_grooves=len(grooves),
        n_points=len(u),
        uncertainty=method,
        per_groove={key: per[:, k].copy() for k, key in enumerate(_DERIVED_KEYS)},
    )

    lattice_err = period_err if centres.size > 3 else float("nan")
    level_out = iso
    return CrossSection(
        x_nm=x_nm,
        height_nm=h_nm,
        period_nm=float(period),
        period_source="given" if period_nm is not None else "measured",
        measured_period_nm=float(period_meas),
        measured_period_err_nm=float(lattice_err),
        groove_centres_nm=np.array([g[2] for g in grooves]),
        placement_rms_nm=float(placement),
        u=u_avg,
        y=y_avg,
        normal_err=normal_err,
        fit=fit,
        level=float(level_out),
        vacuum_level=vacuum,
        material_level=material,
        pixel_x_nm=float(sx),
        pixel_y_nm=float(sy),
        foreshortened=axis,
        foreshortening=float(factor),
        image_x_along_dispersion=image_x_along_dispersion,
    )
