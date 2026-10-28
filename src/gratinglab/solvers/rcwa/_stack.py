"""The layered structure RCWA solves, and how a profile becomes one.

A :class:`Stack` is the whole problem along the normal: a superstrate, any
homogeneous layers above the grating region, the grating's slabs, any
homogeneous layers below it, and a substrate. Every slab has the same shape,
a background index plus segments of other indices along ``x``, so a
homogeneous film, a lamellar ridge and a staircase step with a conformal
coating on it are one type. That is what keeps a layer stack from being a
retrofit (``docs/theory/rcwa.md`` section 6).

Positions along ``x`` are fractions of the period, ``xi = x / period``.
Indices are complex refractive indices at one wavelength; the backend builds
a stack per wavelength, since a dispersive material changes only the indices.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...profiles import Blazed, Lamellar, _runs_to_intervals

__all__ = ["Slab", "Stack", "Groove", "groove_geometry", "fourier_matrix"]


@dataclass(frozen=True, slots=True)
class Slab:
    """One slab between two planes normal to the grating.

    ``segments`` are ``(xi_start, xi_end, index)``, non-overlapping, where
    the slab holds ``index`` instead of ``background``. An interval may
    start below 0 or end above 1; the structure is periodic, so it means the
    same thing wrapped. No segments is a homogeneous layer.
    """

    thickness: float  # nm
    background: complex
    segments: tuple[tuple[float, float, complex], ...] = ()

    @property
    def homogeneous(self) -> bool:
        return not self.segments or all(
            index == self.background for _, _, index in self.segments
        )


@dataclass(frozen=True, slots=True)
class Stack:
    """Superstrate, slabs top to bottom, substrate.

    ``above`` and ``below`` are the homogeneous layers on either side of the
    grating region. Nothing in the solve distinguishes them from the grating
    slabs; they are kept apart so a caller can say which is which and a
    reader of a stack can see the grating.
    """

    superstrate: complex
    grating: tuple[Slab, ...]
    substrate: complex
    above: tuple[Slab, ...] = ()
    below: tuple[Slab, ...] = ()

    @property
    def slabs(self) -> tuple[Slab, ...]:
        return self.above + self.grating + self.below


@dataclass(frozen=True, slots=True)
class Groove:
    """The grating region's geometry, independent of the materials.

    ``slabs`` are ``(thickness_nm, segments)`` top to bottom, where each
    segment is ``(xi_start, xi_end, key)``: key ``0`` is the material the
    grating is cut into and key ``k`` the ``k``-th conformal layer. Anything
    no segment covers is the superstrate. ``exact`` is False for a staircase.
    """

    slabs: tuple[tuple[float, tuple[tuple[float, float, int], ...]], ...]
    exact: bool

    def stack(self, indices, *, superstrate: complex = 1.0) -> Stack:
        """Superstrate above, the slabs, and the grating's material below.

        ``indices[k]`` is the index of key ``k``: the grating's own material,
        then each layer. The substrate is the grating's material, which is the
        mapping a :class:`~gratinglab.problem.Problem` makes (``rcwa.md``
        section 1).
        """
        grating = tuple(
            Slab(
                thickness=thickness,
                background=superstrate,
                segments=tuple((a, b, complex(indices[k])) for a, b, k in segments),
            )
            for thickness, segments in self.slabs
        )
        return Stack(superstrate=superstrate, grating=grating, substrate=indices[0])


#: Bisection steps refining a staircase wall from its sampled position. The
#: bracket starts at a few samples wide (~1e-4 of the period), so 40 halvings
#: put the wall within ~1e-16, which is below anything a truncation can see.
_BISECTIONS = 40

#: Points on which a smooth profile's surface is sampled as a polyline before
#: it is offset by a layer. Its chord error on a sinusoid of depth 0.3 period
#: is ~2e-7 of the period.
_POLYLINE_POINTS = 4096


def _refine(height, guesses, level: float, step: float):
    """Where ``height(xi)`` crosses ``level`` near each guess, by bisection.

    Vectorised over the guesses. A guess whose bracket holds no sign change
    (a surface touching the level without crossing it) is kept as it is.
    """
    guesses = np.asarray(guesses, dtype=np.float64)
    low, high = guesses - 2.0 * step, guesses + 2.0 * step
    f_low = height(low) - level >= 0.0
    f_high = height(high) - level >= 0.0
    crossing = f_low != f_high
    for _ in range(_BISECTIONS):
        middle = 0.5 * (low + high)
        same = (height(middle) - level >= 0.0) == f_low
        low = np.where(same, middle, low)
        high = np.where(same, high, middle)
    return np.where(crossing, 0.5 * (low + high), guesses)


def _surface(profile, period: float):
    """One period of the surface as a polyline along ``+x``, in nm.

    Exact vertices for the polygonal profiles, so a vertical wall is a
    vertical segment; a dense sampling of ``height`` otherwise.
    """
    if isinstance(profile, Lamellar):
        w, d = profile.duty_cycle, profile.depth
        t = np.array([0.0, w, w, 1.0, 1.0])
        y = np.array([d, d, 0.0, 0.0, d])
    elif isinstance(profile, Blazed):
        d = profile.depth
        if profile.has_vertical_facet:
            t, y = np.array([0.0, 1.0, 1.0]), np.array([0.0, d, 0.0])
        else:
            t, y = np.array([0.0, profile.apex, 1.0]), np.array([0.0, d, 0.0])
    else:
        t = np.linspace(0.0, 1.0, _POLYLINE_POINTS + 1)
        y = np.asarray(profile.height(t), dtype=np.float64)
    # xi = 1 - t (conventions section 3); reversed so x ascends.
    return ((1.0 - t) * period)[::-1], (y * period)[::-1]


def _offset_height(x_surface, y_surface, period: float, thickness: float):
    """``H(xi)``: the top of every point within ``thickness`` of the surface.

    That is the true offset (``ConformalLayer``). For each segment, the
    highest point of its stadium of radius ``thickness`` above ``x``, which is
    concave along the segment and so maximised in closed form, then the
    largest over segments. Periodic copies cover ``thickness`` beyond a
    period.
    """
    copies = int(np.ceil(thickness / period)) + 1
    shifts = np.arange(-copies, copies + 1) * period
    xs = (x_surface[None, :] + shifts[:, None]).ravel()
    ys = np.tile(y_surface, len(shifts))
    xa, xb, ya, yb = xs[:-1], xs[1:], ys[:-1], ys[1:]
    lo, hi = np.minimum(xa, xb), np.maximum(xa, xb)
    run = xb - xa
    vertical = np.abs(run) <= 1e-12 * period
    slope = np.where(vertical, 0.0, (yb - ya) / np.where(vertical, 1.0, run))
    lean = thickness * slope / np.sqrt(1.0 + slope**2)
    top = np.maximum(ya, yb)

    def height(xi):
        xi = np.atleast_1d(np.asarray(xi, dtype=np.float64))
        out = np.full(xi.shape, -np.inf)
        for chunk in range(0, len(xi), 256):
            x = xi[chunk : chunk + 256, None] * period
            left = np.maximum(lo, x - thickness)
            right = np.minimum(hi, x + thickness)
            best = np.clip(x + lean, left, right)
            along = np.where(vertical, top, ya + slope * (best - xa))
            reach = np.where(vertical, lo, best) - x
            value = along + np.sqrt(np.maximum(thickness**2 - reach**2, 0.0))
            value = np.where(
                (left <= right) & (np.abs(reach) <= thickness), value, -np.inf
            )
            out[chunk : chunk + 256] = value.max(axis=1)
        return out

    return height


def _wrap(intervals):
    """Intervals folded into ``[0, 1)``, sorted and merged."""
    pieces = []
    for a, b in intervals:
        if b - a >= 1.0:
            return [(0.0, 1.0)]
        a0 = a % 1.0
        b0 = a0 + (b - a)
        if b0 <= 1.0:
            pieces.append((a0, b0))
        else:
            pieces += [(a0, 1.0), (0.0, b0 - 1.0)]
    pieces.sort()
    merged: list[tuple[float, float]] = []
    for a, b in pieces:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def _subtract(outer, inner):
    """``outer`` minus ``inner``, both wrapped interval lists."""
    out = []
    for a, b in outer:
        start = a
        for c, d in inner:
            if d <= start or c >= b:
                continue
            if c > start:
                out.append((start, c))
            start = max(start, d)
        if start < b:
            out.append((start, b))
    return [(a, b) for a, b in out if b - a > 1e-15]


def _inside(height, level: float, samples: int):
    """The spans of ``xi`` where ``level <= height(xi)``, walls refined."""
    xi = np.arange(samples) / samples
    spans = _runs_to_intervals(height(xi) >= level, xi)
    partial = [span for span in spans if span[1] - span[0] < 1.0]
    if not partial:
        return list(spans)
    ends = _refine(height, [edge for span in partial for edge in span], level, 1.0 / samples)
    refined = iter(zip(ends[0::2], ends[1::2]))
    return [span if span[1] - span[0] >= 1.0 else next(refined) for span in spans]


def _layered(profile, period: float, slices: int, thicknesses) -> Groove:
    """The staircase of a grating with conformal layers on it.

    Surface ``k`` is the profile offset by the first ``k`` thicknesses, and
    ``k = 0`` is the profile itself; a point is key ``k`` when it lies under
    surface ``k`` and over surface ``k - 1``. Slab boundaries are placed at
    every height where a surface has a flat (the profile's top and bottom,
    and each layer's offsets of them), so a lamellar groove's flats are exact
    and only the rounded corners are staircased.
    """
    x_surface, y_surface = _surface(profile, period)
    cumulative = np.cumsum(thicknesses)

    def bare(xi):
        return period * np.asarray(profile.height(1.0 - xi), dtype=np.float64)

    heights = [bare] + [
        _offset_height(x_surface, y_surface, period, float(c)) for c in cumulative
    ]
    bottom, depth = float(y_surface.min()), float(y_surface.max())
    keys = sorted(
        {bottom, depth}
        | {bottom + float(c) for c in cumulative}
        | {depth + float(c) for c in cumulative}
    )
    total = keys[-1] - keys[0]
    edges: list[float] = []
    for low, high in zip(keys[:-1], keys[1:]):
        count = max(1, int(round(slices * (high - low) / total)))
        edges.extend(np.linspace(low, high, count + 1)[:-1].tolist())
    edges.append(keys[-1])

    samples = max(2048, 64 * slices)
    slabs = []
    for low, high in reversed(list(zip(edges[:-1], edges[1:]))):
        level = 0.5 * (low + high)
        segments = []
        below: list[tuple[float, float]] = []
        for key, height in enumerate(heights):
            covered = _wrap(_inside(height, level, samples))
            segments += [(a, b, key) for a, b in _subtract(covered, below)]
            below = covered
        slabs.append((high - low, tuple(segments)))
    return Groove(slabs=tuple(slabs), exact=False)


def groove_geometry(profile, period: float, slices: int, thicknesses=()) -> Groove:
    """The grating region of ``profile``, as slabs top to bottom.

    Bare, a :class:`~gratinglab.profiles.Lamellar` groove is one exact slab,
    whatever ``slices`` says, and anything else a staircase of ``slices``
    equal slabs by the mid-height rule of ``Profile.slice_layers``, with each
    wall moved from its sampled position to the profile's own crossing.
    ``thicknesses`` adds conformal layers (:func:`_layered`), and the region
    then reaches the outermost layer's top.

    The profile's ``t`` runs along ``-d`` (``docs/conventions.md`` section 3),
    so ``t`` becomes ``xi = 1 - t``, the placement the integral solver uses.
    """
    if len(thicknesses):
        return _layered(profile, period, slices, thicknesses)
    if isinstance(profile, Lamellar):
        segments = ((1.0 - profile.duty_cycle, 1.0, 0),)
        return Groove(slabs=((profile.depth * period, segments),), exact=True)

    def height(t):
        return np.asarray(profile.height(t), dtype=np.float64)

    layers = profile.slice_layers(slices)
    samples = max(2048, 64 * slices)  # the sampling slice_layers used
    slabs = []
    for layer in reversed(layers):  # slice_layers is bottom to top
        level = 0.5 * (layer.y_lower + layer.y_upper)
        segments = []
        for start, end in layer.intervals:
            if end - start >= 1.0:
                segments.append((0.0, 1.0, 0))
                continue
            start, end = _refine(height, [start, end], level, 1.0 / samples)
            segments.append((1.0 - float(end), 1.0 - float(start), 0))
        slabs.append((layer.thickness * period, tuple(segments)))
    return Groove(slabs=tuple(slabs), exact=False)


def _coefficients(slab: Slab, power: int, harmonics: int) -> NDArray[np.complex128]:
    """Fourier coefficients ``c_n``, ``n = -2M..2M``, of ``eps(xi)**power``.

    ``eps(xi) = sum_n c_n exp(2 pi i n xi)``, so ``c_n`` is the period average
    of ``eps**power * exp(-2 pi i n xi)``: exact for a piecewise-constant slab.
    """
    n = np.arange(-2 * harmonics, 2 * harmonics + 1)
    base = complex(slab.background) ** 2
    background = base**power
    coefficients = np.where(n == 0, background, 0.0).astype(np.complex128)
    nonzero = n != 0
    for start, end, index in slab.segments:
        value = (complex(index) ** 2) ** power - background
        if value == 0.0:
            continue
        integral = np.empty(len(n), dtype=np.complex128)
        integral[~nonzero] = end - start
        k = -2j * np.pi * n[nonzero]
        integral[nonzero] = (np.exp(k * end) - np.exp(k * start)) / k
        coefficients += value * integral
    return coefficients


def fourier_matrix(slab: Slab, harmonics: int, *, inverse: bool) -> NDArray[np.complex128]:
    """The ``(2M+1)``-square Toeplitz matrix for a product ``eps E``.

    ``inverse=False`` is Laurent's rule, ``[[eps]]``. ``inverse=True`` is the
    inverse rule, ``[[1/eps]]^-1``, which is the correct one for a field
    component that jumps where ``eps`` does (Li 1996; ``rcwa.md`` section 3).
    """
    size = 2 * harmonics + 1
    coefficients = _coefficients(slab, -1 if inverse else 1, harmonics)
    rows = np.arange(size)
    toeplitz = coefficients[(rows[:, None] - rows[None, :]) + 2 * harmonics]
    return np.linalg.inv(toeplitz) if inverse else toeplitz
