"""Groove profiles, in the four representations the solver methods need.

Profiles are normalised. Both the coordinate ``t`` and the returned height
are expressed in units of the period, so a profile is a pure shape and the
period is stored once, on :class:`~gratinglab.problem.Problem`. A solver that needs
physical height multiplies by ``problem.period``. This also matches PCGrate's
``.ggp`` format, which stores ``Period: 1`` and normalised points.

The four methods require different geometry from the same profile:

=================  ==========================================================
Method             Representation
=================  ==========================================================
Scalar             :meth:`Profile.height` -- the groove phase function
RCWA               :meth:`Profile.slice_layers` -- piecewise-constant slabs
C-method           :meth:`Profile.slope` -- smooth dy/dt; no vertical facets
Integral           :meth:`Profile.boundary` -- parametrised curve with
                   outward normals and arc-length quadrature weights
=================  ==========================================================

All four were designed together deliberately, because retrofitting normals and
arc lengths onto a shipped API once the integral method arrived would have
required changes to every profile type.

A profile that cannot supply a representation raises an error. ``Blazed`` with
a vertical anti-blaze facet has no finite slope there, so :meth:`Blazed.slope`
refuses instead of returning a large number that would silently corrupt a
C-method solve.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pydantic import BaseModel, ConfigDict, Field, model_validator

__all__ = [
    "Profile",
    "Layer",
    "BoundaryCurve",
    "ProfileRepresentationError",
    "Blazed",
    "Lamellar",
    "Sinusoidal",
    "FromProfileData",
]


class ProfileRepresentationError(NotImplementedError):
    """A profile cannot be expressed in the representation a solver asked for.

    The error is raised in place of an approximation, because a C-method solve
    given a silently smoothed vertical facet produces a plausible but wrong
    answer.
    """


@dataclass(frozen=True, slots=True)
class Layer:
    """One piecewise-constant slab of an RCWA discretisation.

    ``intervals`` lists the ``(t_start, t_end)`` spans, normalised to the
    period, where grating material occupies this slab. The complement is the
    incident medium.
    """

    y_lower: float
    y_upper: float
    intervals: tuple[tuple[float, float], ...]

    @property
    def thickness(self) -> float:
        return self.y_upper - self.y_lower

    @property
    def fill_factor(self) -> float:
        """Fraction of the period occupied by grating material."""
        return sum(end - start for start, end in self.intervals)


@dataclass(frozen=True, slots=True)
class BoundaryCurve:
    """A parametrised boundary over one period, for integral-equation solvers.

    Attributes
    ----------
    t, y
        Points along the profile, normalised to the period.
    nx, ny
        Outward unit normal, pointing into the incident medium. ``ny > 0``
        wherever the boundary is a height function; a vertical facet has a
        horizontal normal (``ny = 0``), and an undercut portion can have
        ``ny < 0``. Never inward.
    ds
        Arc-length weight for each point, for quadrature. Sums to the total
        profile arc length over one period. By default points are spaced at
        equal arc length, which the integral method's rectangular quadrature
        rule needs for its spectral accuracy on a smooth curve.
    graded
        True when the points come from a corner-graded mesh
        (``boundary(n, grading=p)`` on a profile with corners). The points
        are then equally spaced in a mesh parameter, not in arc length, and
        ``ds`` holds the exact quadrature weight ``s'(tau) / n`` of each
        point, which the integral method must use node by node.
    polygon
        On a graded mesh, the open vertex list ``(t, y)`` of the polygon the
        points lie on, over one period and normalised like ``t`` and ``y``.
        The integral method integrates its corner kernels over these straight
        pieces in closed form. ``None`` on the equal-arc mesh.
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]
    nx: NDArray[np.float64]
    ny: NDArray[np.float64]
    ds: NDArray[np.float64]
    graded: bool = False
    polygon: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None

    @property
    def arc_length(self) -> float:
        return float(self.ds.sum())


@runtime_checkable
class Profile(Protocol):
    """The interface every groove profile provides."""

    @property
    def depth(self) -> float:
        """Peak-to-valley height, in units of the period."""

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        """Profile height at normalised position ``t``, periodic in ``t``."""

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        """dy/dt. Raises :class:`ProfileRepresentationError` if not finite."""

    def slice_layers(self, n: int) -> tuple[Layer, ...]:
        """Piecewise-constant slabs, bottom to top."""

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Parametrised boundary with outward normals and arc-length weights.

        ``grading = p >= 2`` clusters the points at the profile's corners
        with the sigmoidal transformation of order ``p`` (see
        :func:`_polyline_boundary`); ``0`` keeps equal arc length.
        """

    def is_single_valued(self) -> bool:
        """False for undercut profiles, which RCWA tolerates and C-method cannot."""


class _BaseProfile(BaseModel):
    """Shared machinery. Concrete profiles supply ``depth`` and ``height``."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    def height(self, t: ArrayLike) -> NDArray[np.float64]:  # pragma: no cover
        raise NotImplementedError

    @property
    def depth(self) -> float:  # pragma: no cover
        raise NotImplementedError

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        """Numerical dy/dt by central difference.

        Subclasses with an analytic derivative should override; subclasses with
        a discontinuity must raise :class:`ProfileRepresentationError`.
        """
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        h = 1e-7
        return (self.height(t + h) - self.height(t - h)) / (2 * h)

    def is_single_valued(self) -> bool:
        """True unless a subclass models an undercut."""
        return True

    def slice_layers(self, n: int) -> tuple[Layer, ...]:
        """Horizontal slabs, bottom to top.

        Each slab records where ``height(t)`` reaches at least the slab's
        mid-height, which is the standard staircase approximation RCWA uses.
        Accuracy improves with ``n``; the convergence harness sweeps it.
        """
        if n < 1:
            raise ValueError(f"need at least one layer, got {n}")
        if not self.is_single_valued():
            raise ProfileRepresentationError(
                "cannot slice an undercut profile into horizontal layers by "
                "height alone; the boundary is not a function of t"
            )

        edges = np.linspace(0.0, self.depth, n + 1)
        # Sample finely so interval edges are located to better than a layer.
        samples = max(2048, 64 * n)
        t = np.linspace(0.0, 1.0, samples, endpoint=False)
        y = self.height(t)

        layers = []
        for lower, upper in zip(edges[:-1], edges[1:]):
            inside = y >= 0.5 * (lower + upper)
            layers.append(
                Layer(
                    y_lower=float(lower),
                    y_upper=float(upper),
                    intervals=_runs_to_intervals(inside, t),
                )
            )
        return tuple(layers)

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Boundary resampled at equal arc length, with outward normals.

        Samples ``height`` finely and resamples the polyline at equal arc
        length, which is the node spacing the integral method's rectangular
        quadrature rule needs. Profiles whose geometry is an exact polygon
        (:class:`Blazed`, :class:`Lamellar`) override with their vertices so
        facet corners and vertical segments carry no sampling error. A
        finely sampled smooth height has no corners, so ``grading`` leaves
        it on the equal-arc mesh.
        """
        if n < 3:
            raise ValueError(f"need at least 3 boundary points, got {n}")
        fine = max(4096, 8 * n)
        t = np.linspace(0.0, 1.0, fine, endpoint=False)
        return _polyline_boundary(t, self.height(t), n, grading=grading)


#: A polygon vertex turning the boundary by more than this is a corner, at
#: which a graded mesh clusters its points. Exact polygons (``Blazed``,
#: ``Lamellar``) turn by tens of degrees at every vertex; a finely sampled
#: smooth height turns by a fraction of a degree, so it has none.
CORNER_TURN_DEGREES = 10.0

#: Fewest points a graded mesh gives one side between two corners. Below
#: this the sigmoid cannot both cluster at each end and cover the middle.
_MIN_POINTS_PER_SIDE = 4


def _grading_map(
    u: NDArray[np.float64], order: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Kress's sigmoidal transformation of ``[0, 1]`` onto itself, and its slope.

    ``w(u) = v^p / (v^p + (1 - v)^p)`` with the cubic
    ``v = (1/2 - 1/p) xi^3 + xi/p + 1/2`` of Kress (1990), ``xi = 2u - 1``,
    which keeps the middle of a side at twice the equal-arc spacing rather
    than ``p`` times it.
    Every derivative of ``w`` below order ``p`` vanishes at both ends, so a
    corner density that is singular in arc length becomes, multiplied by the
    slope, smooth in ``u``, and the rectangular rule recovers high order.
    The map is odd about ``u = 1/2`` (``w(1 - u) = 1 - w(u)``), which the
    caller uses to place points near the far corner without cancellation.
    """
    p = float(order)
    xi = 2.0 * u - 1.0
    c = 1.0 / p - 0.5
    v = -c * xi**3 + xi / p + 0.5
    dv = 2.0 * (-3.0 * c * xi**2 + 1.0 / p)
    vp, wp = v**p, (1.0 - v) ** p
    w = vp / (vp + wp)
    dw = p * (v * (1.0 - v)) ** (p - 1.0) / (vp + wp) ** 2 * dv
    return w, dw


def _graded_polyline(
    tx: NDArray[np.float64],
    ty: NDArray[np.float64],
    n: int,
    order: int,
) -> BoundaryCurve | None:
    """A corner-graded mesh on the closed polyline ``(tx, ty)``.

    Returns ``None`` when the polyline has no corner, so the caller keeps the
    equal-arc mesh. Otherwise each side between consecutive corners receives
    points in proportion to its length, placed at the midpoints of an equal
    split of the mesh parameter and mapped through :func:`_grading_map`, so
    no point sits on a corner and every point lies on a straight piece with
    an exact normal. The weights are the exact ``ds/dtau / n``.
    """
    # Drop repeated vertices: a zero-length edge has no direction.
    keep = np.hypot(np.diff(tx), np.diff(ty)) > 0.0
    vx = tx[:-1][keep]
    vy = ty[:-1][keep]
    m = len(vx)
    if m < 2:
        return None
    ex = np.append(vx[1:], vx[0] + 1.0) - vx
    ey = np.append(vy[1:], vy[0]) - vy
    prev_x, prev_y = np.roll(ex, 1), np.roll(ey, 1)
    turn = np.degrees(
        np.abs(np.arctan2(prev_x * ey - prev_y * ex, prev_x * ex + prev_y * ey))
    )
    corners = np.flatnonzero(turn > CORNER_TURN_DEGREES)
    if corners.size == 0:
        return None
    if n < _MIN_POINTS_PER_SIDE * corners.size:
        raise ValueError(
            f"a graded mesh needs at least {_MIN_POINTS_PER_SIDE} points on "
            f"each of the profile's {corners.size} sides between corners, "
            f"so n >= {_MIN_POINTS_PER_SIDE * corners.size}; got {n}"
        )

    # Two periods of vertices, so a side may run through the period wrap.
    big_x = np.concatenate((vx, vx + 1.0, [vx[0] + 2.0]))
    big_y = np.concatenate((vy, vy, [vy[0]]))
    edge = np.hypot(np.diff(big_x), np.diff(big_y))
    cumulative = np.concatenate(([0.0], np.cumsum(edge)))
    total = float(cumulative[m])

    starts = corners
    ends = np.append(corners[1:], corners[0] + m)
    lengths = cumulative[ends] - cumulative[starts]

    # Largest-remainder allocation, then the per-side floor.
    share = n * lengths / total
    counts = np.maximum(np.floor(share).astype(int), _MIN_POINTS_PER_SIDE)
    while counts.sum() < n:
        counts[int(np.argmax(share - counts))] += 1
    while counts.sum() > n:
        spare = np.where(counts > _MIN_POINTS_PER_SIDE, counts - share, -np.inf)
        counts[int(np.argmax(spare))] -= 1

    def locate(first: int, last: int, offset, forward: bool):
        """Points at ``offset`` along the side, measured from its nearer end.

        Measuring from the nearer corner keeps a point a distance 1e-9 from
        its corner exactly that far, rather than the difference of two
        nearly equal arc positions.
        """
        index = np.arange(first, last + 1)
        if not forward:
            index = index[::-1]
        px, py = big_x[index], big_y[index]
        steps = np.hypot(np.diff(px), np.diff(py))
        along = np.concatenate(([0.0], np.cumsum(steps)))
        piece = np.clip(
            np.searchsorted(along, offset, side="right") - 1, 0, len(steps) - 1
        )
        fraction = (offset - along[piece]) / steps[piece]
        dx = px[piece + 1] - px[piece]
        dy = py[piece + 1] - py[piece]
        sign = 1.0 if forward else -1.0
        return (
            px[piece] + fraction * dx,
            py[piece] + fraction * dy,
            sign * dx / steps[piece],
            sign * dy / steps[piece],
        )

    xs, ys, txs, tys, ds, arc = [], [], [], [], [], []
    for first, last, length, count in zip(starts, ends, lengths, counts):
        u = (np.arange(count) + 0.5) / count
        w, dw = _grading_map(u, order)
        near = u <= 0.5
        far_w, _ = _grading_map(1.0 - u[~near], order)
        px = np.empty(count)
        py = np.empty(count)
        qx = np.empty(count)
        qy = np.empty(count)
        px[near], py[near], qx[near], qy[near] = locate(
            first, last, length * w[near], True
        )
        px[~near], py[~near], qx[~near], qy[~near] = locate(
            first, last, length * far_w, False
        )
        xs.append(px)
        ys.append(py)
        txs.append(qx)
        tys.append(qy)
        ds.append(length * dw / count)
        arc.append(cumulative[first] + length * w)

    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    tx_all = np.concatenate(txs)
    ty_all = np.concatenate(tys)
    ds_all = np.concatenate(ds)
    arc_all = np.concatenate(arc)

    # Points past one period belong to the start of this one: shift them
    # back and rotate them to the front, so the order runs from x[0] like
    # the equal-arc mesh and only the last point wraps.
    over = arc_all >= total
    x_all[over] -= 1.0
    shift = int(over.sum())
    order_index = np.roll(np.arange(n), shift)
    # Tangent (tx, ty) rotated +90 deg is the outward normal (see below).
    return BoundaryCurve(
        t=x_all[order_index],
        y=y_all[order_index],
        nx=-ty_all[order_index],
        ny=tx_all[order_index],
        ds=ds_all[order_index],
        graded=True,
        polygon=(vx.copy(), vy.copy()),
    )


def _polyline_boundary(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    n: int,
    *,
    grading: int = 0,
) -> BoundaryCurve:
    """Boundary samples of a polyline spanning one period.

    ``x`` and ``y`` are the open polygon vertices, normalised to the period;
    the curve is closed by continuing to ``(x[0] + 1, y[0])``, so a trailing
    vertical segment (a lamellar sidewall, an ideal sawtooth's anti-blaze
    facet) is represented as real geometry instead of a steep approximation. Handles
    boundaries that double back in ``x`` (undercut), which is why every
    ``boundary()`` implementation funnels through here.

    With ``grading = 0`` the samples are at equal arc length. With
    ``grading = p >= 2`` and at least one corner (a vertex turning by more
    than :data:`CORNER_TURN_DEGREES`) they are graded towards the corners by
    Kress's sigmoidal transformation of order ``p`` (M&P section 4.6.6;
    Kress 1990), which restores high-order convergence where the field is
    singular. With no corner the equal-arc mesh is returned unchanged.
    """
    if n < 3:
        raise ValueError(f"need at least 3 boundary points, got {n}")
    if grading < 0 or grading == 1:
        raise ValueError(
            f"grading must be 0 (equal arc length) or at least 2, got "
            f"{grading}; order 1 is the identity map and clusters nothing"
        )

    # Close the polygon: the curve continues into the next period.
    tx = np.append(x, x[0] + 1.0)
    ty = np.append(y, y[0])

    segment = np.hypot(np.diff(tx), np.diff(ty))
    cumulative = np.concatenate(([0.0], np.cumsum(segment)))
    total = cumulative[-1]
    if total <= 0:
        raise ValueError("profile has zero arc length")

    if grading:
        graded = _graded_polyline(tx, ty, n, grading)
        if graded is not None:
            return graded

    stations = np.linspace(0.0, total, n, endpoint=False)
    xs = np.interp(stations, cumulative, tx)
    ys = np.interp(stations, cumulative, ty)

    # Periodic central differences; x wraps by exactly one period.
    dx = np.roll(xs, -1) - np.roll(xs, 1)
    dx[0] += 1.0
    dx[-1] += 1.0
    dy = np.roll(ys, -1) - np.roll(ys, 1)
    norm = np.hypot(dx, dy)

    # Tangent (dx, dy) rotated +90 deg is (-dy, dx), which points into the
    # incident medium for a curve traversed in the +t direction.
    forward = np.hypot(np.roll(xs, -1) - xs, np.roll(ys, -1) - ys)
    forward[-1] = np.hypot(xs[0] + 1.0 - xs[-1], ys[0] - ys[-1])
    backward = np.roll(forward, 1)

    return BoundaryCurve(
        t=xs,
        y=ys,
        nx=-dy / norm,
        ny=dx / norm,
        ds=0.5 * (forward + backward),
    )


def _runs_to_intervals(
    mask: NDArray[np.bool_], t: NDArray[np.float64]
) -> tuple[tuple[float, float], ...]:
    """Contiguous True runs of ``mask`` as ``(start, end)`` spans in ``t``.

    Runs that wrap the period boundary are merged, since the profile is
    periodic and a split run would misrepresent a single feature as two.
    """
    if not mask.any():
        return ()
    if mask.all():
        return ((0.0, 1.0),)

    step = 1.0 / len(t)
    edges = np.flatnonzero(np.diff(mask.astype(np.int8)))
    starts = [i + 1 for i in edges if not mask[i]]
    ends = [i + 1 for i in edges if mask[i]]

    if mask[0]:
        starts.insert(0, 0)
    if mask[-1]:
        ends.append(len(mask))

    spans = [(t[s], t[e - 1] + step) for s, e in zip(starts, ends)]

    # Merge a run touching t=0 with one touching t=1: they are one feature.
    if len(spans) > 1 and spans[0][0] == 0.0 and spans[-1][1] >= 1.0:
        spans = [(spans[-1][0] - 1.0, spans[0][1])] + spans[1:-1]

    return tuple((float(a), float(b)) for a, b in spans)


class Blazed(_BaseProfile):
    """Sawtooth groove, the workhorse of blazed reflection gratings.

    ``depth = 1 / (cot(blaze) + cot(antiblaze))`` in units of the period, which
    reduces to ``tan(blaze)`` for the ideal sawtooth with a vertical anti-blaze
    facet (``docs/conventions.md`` 9).

    >>> round(Blazed(blaze_angle=30.0, antiblaze_angle=80.0).depth, 4)
    0.5241
    """

    blaze_angle: float = Field(gt=0.0, lt=90.0, description="Active facet angle, deg")
    antiblaze_angle: float = Field(
        default=90.0, gt=0.0, le=90.0, description="Opposite facet angle, deg"
    )

    @property
    def _runs(self) -> tuple[float, float]:
        """Horizontal extent of each facet, normalised. Sums to 1."""
        cot_b = 1.0 / np.tan(np.radians(self.blaze_angle))
        cot_a = 1.0 / np.tan(np.radians(self.antiblaze_angle))
        total = cot_b + cot_a
        return cot_b / total, cot_a / total

    @property
    def depth(self) -> float:
        cot_b = 1.0 / np.tan(np.radians(self.blaze_angle))
        cot_a = 1.0 / np.tan(np.radians(self.antiblaze_angle))
        return 1.0 / (cot_b + cot_a)

    @property
    def apex(self) -> float:
        """Normalised position of the facet crest."""
        return self._runs[0]

    @property
    def has_vertical_facet(self) -> bool:
        return bool(np.isclose(self.antiblaze_angle, 90.0))

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        apex, depth = self.apex, self.depth
        rising = depth * t / apex
        if self.has_vertical_facet:
            return np.where(t < 1.0, rising, 0.0)
        falling = depth * (1.0 - t) / (1.0 - apex)
        return np.where(t <= apex, rising, falling)

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        if self.has_vertical_facet:
            raise ProfileRepresentationError(
                f"{type(self).__name__} with antiblaze_angle=90 has a vertical "
                "facet and no finite slope there. The C-method cannot represent "
                "it; use a slightly shallower antiblaze_angle, or a method that "
                "handles discontinuities (RCWA, integral)."
            )
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        up = np.tan(np.radians(self.blaze_angle))
        down = -np.tan(np.radians(self.antiblaze_angle))
        return np.where(t <= self.apex, up, down)

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Exact sawtooth polygon.

        A vertical anti-blaze facet is the closing segment of the polygon, a
        piece of geometry that a height function cannot express.
        """
        if self.has_vertical_facet:
            x, y = np.array([0.0, 1.0]), np.array([0.0, self.depth])
        else:
            x, y = np.array([0.0, self.apex]), np.array([0.0, self.depth])
        return _polyline_boundary(x, y, n, grading=grading)


class Lamellar(_BaseProfile):
    """Rectangular (binary, laminar) groove."""

    depth_fraction: float = Field(gt=0.0, description="Depth / period")
    duty_cycle: float = Field(
        default=0.5, gt=0.0, lt=1.0, description="Fraction of the period at full height"
    )

    @property
    def depth(self) -> float:
        return self.depth_fraction

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        return np.where(t < self.duty_cycle, self.depth, 0.0)

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        raise ProfileRepresentationError(
            "Lamellar has vertical sidewalls and no finite slope. The C-method "
            "cannot represent it; use RCWA or the integral method."
        )

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Exact rectangular polygon; both sidewalls are genuine vertical segments."""
        d, w = self.depth, self.duty_cycle
        x = np.array([0.0, w, w, 1.0])
        y = np.array([d, d, 0.0, 0.0])
        return _polyline_boundary(x, y, n, grading=grading)


class Sinusoidal(_BaseProfile):
    """Holographic (sinusoidal) groove, smooth and therefore suited to the C-method."""

    depth_fraction: float = Field(gt=0.0, description="Peak-to-valley / period")

    @property
    def depth(self) -> float:
        return self.depth_fraction

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.asarray(t, dtype=np.float64)
        return 0.5 * self.depth * (1.0 - np.cos(2.0 * np.pi * t))

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.asarray(t, dtype=np.float64)
        return np.pi * self.depth * np.sin(2.0 * np.pi * t)


class FromProfileData(_BaseProfile):
    """A measured profile, from AFM or from PCGrate's ``.ggp`` polygonal border.

    Points are normalised: ``t`` spans one period on ``[0, 1]`` and ``y`` is in
    units of the period. Use :meth:`from_measurements` to normalise raw nm data.

    The profile is stored as given and interpolated linearly, matching how
    PCGrate treats a polygonal border, so a scalar calculation and a PCGrate
    run can be driven from the identical geometry.

    Undercut profiles are representable. If ``t`` doubles back, the boundary
    is not a function of position and :meth:`height` and
    :meth:`~_BaseProfile.slice_layers` refuse, but :meth:`boundary` still
    works, because an integral-equation solver parametrises the curve instead
    of assuming a height function. This asymmetry reflects the physics of the
    methods and is not a limitation of this class: RCWA and scalar theory
    cannot represent an undercut, and the integral method can.

    ``blaze_angle`` is optional and has a different meaning here from the one
    it has on :class:`Blazed`. There it is a definition: the facet has that
    angle because it was specified so. Here it is a measurement, obtained by
    fitting a line to a real facet, and it carries an uncertainty. Both feed
    ``solvers.scalar._reflecting_graze`` identically, and the provenance string
    records which kind was used, because a reader deciding whether to trust an
    absolute efficiency needs to know.

    Leaving it ``None`` has consequences. The reflection then falls back to
    the mean surface, a zero-tilt facet, which for a sawtooth measured at
    grazing incidence is the wrong plane by the whole blaze angle. The default
    is ``None`` because most measured profiles arrive without a fit, and not
    because ``None`` is the better answer.
    """

    t: tuple[float, ...]
    y: tuple[float, ...]
    blaze_angle: float | None = None

    @model_validator(mode="after")
    def _check_points(self) -> "FromProfileData":
        if len(self.t) != len(self.y):
            raise ValueError(f"t has {len(self.t)} points, y has {len(self.y)}")
        if len(self.t) < 3:
            raise ValueError(f"need at least 3 points, got {len(self.t)}")
        t = np.asarray(self.t)
        if t.min() < 0.0 or t.max() > 1.0:
            raise ValueError(f"t must lie in [0, 1], got [{t.min()}, {t.max()}]")
        if self.blaze_angle is not None and not 0.0 < self.blaze_angle < 90.0:
            # The same open interval `Blazed` enforces. A fit that lands outside
            # it is not a shallow facet, it is a fit that found the wrong facet
            # or the wrong sign convention, and passing it on as geometry would
            # produce a confident, wrong reflectivity.
            raise ValueError(
                f"blaze_angle must lie in (0, 90) degrees, got {self.blaze_angle}"
            )
        return self

    @classmethod
    def from_measurements(
        cls,
        x_nm: ArrayLike,
        y_nm: ArrayLike,
        period_nm: float,
        *,
        blaze_angle: float | None = None,
    ) -> "FromProfileData":
        """Normalise a measured profile in nm, referencing its minimum to zero.

        Pass ``blaze_angle`` when the same measurement also produced a facet
        fit; see the class docstring for why it is worth carrying.
        """
        x = np.asarray(x_nm, dtype=np.float64)
        y = np.asarray(y_nm, dtype=np.float64)
        if period_nm <= 0:
            raise ValueError(f"period_nm must be positive, got {period_nm}")
        return cls(
            t=tuple((x - x.min()) / period_nm),
            y=tuple((y - y.min()) / period_nm),
            blaze_angle=blaze_angle,
        )

    @property
    def depth(self) -> float:
        return float(np.ptp(np.asarray(self.y)))

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        if not self.is_single_valued():
            raise ProfileRepresentationError(
                "this profile is undercut -- its boundary doubles back in t, so "
                "height is not a function of position. Use boundary() with an "
                "integral-equation solver; scalar theory and RCWA cannot "
                "represent an undercut."
            )
        query = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        return np.interp(query, np.asarray(self.t), np.asarray(self.y), period=1.0)

    def is_single_valued(self) -> bool:
        """A measured boundary that doubles back in ``t`` is undercut."""
        return bool((np.diff(np.asarray(self.t)) >= 0).all())

    @property
    def apex(self) -> float:
        """Normalised position of the highest measured point."""
        return float(np.asarray(self.t)[int(np.argmax(np.asarray(self.y)))])

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Resample the stored polygon at equal arc length, or graded.

        Works from the stored points directly instead of through ``height``,
        which allows it to handle an undercut boundary.
        """
        return _polyline_boundary(
            np.asarray(self.t, dtype=np.float64),
            np.asarray(self.y, dtype=np.float64),
            n,
            grading=grading,
        )
