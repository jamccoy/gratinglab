"""Polygons of (t, height) points: finding their corners, smoothing between
them, and meshing the boundary.

Functions on arrays, with no state. The profile classes call these for their
boundaries, the integral solver (``_boundary.py``, ``_nystrom.py``) for measured
polygons, and RCWA's ``_stack.py`` for ``_runs_to_intervals``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._types import BoundaryCurve


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


def _vertex_turns(
    tx: NDArray[np.float64], ty: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """The distinct vertices of the closed polyline ``(tx, ty)`` and the
    angle in degrees the boundary turns through at each.

    ``(tx, ty)`` ends on the first vertex one period on. Repeated vertices
    are dropped, because a zero-length edge has no direction.
    """
    keep = np.hypot(np.diff(tx), np.diff(ty)) > 0.0
    vx = tx[:-1][keep]
    vy = ty[:-1][keep]
    if len(vx) < 2:
        return vx, vy, np.zeros(len(vx))
    ex = np.append(vx[1:], vx[0] + 1.0) - vx
    ey = np.append(vy[1:], vy[0]) - vy
    prev_x, prev_y = np.roll(ex, 1), np.roll(ey, 1)
    turn = np.degrees(
        np.abs(np.arctan2(prev_x * ey - prev_y * ex, prev_x * ex + prev_y * ey))
    )
    return vx, vy, turn


def polyline_corner_count(x: ArrayLike, y: ArrayLike) -> int:
    """How many vertices of the open vertex list ``(x, y)`` of a polygon
    spanning one period turn by more than :data:`CORNER_TURN_DEGREES`,
    closing through ``(x[0] + 1, y[0])`` as :func:`_polyline_boundary` does."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    _, _, turn = _vertex_turns(np.append(x, x[0] + 1.0), np.append(y, y[0]))
    return int((turn > CORNER_TURN_DEGREES).sum())


def polyline_has_corners(x: ArrayLike, y: ArrayLike) -> bool:
    """Whether :func:`polyline_corner_count` finds any corner."""
    return polyline_corner_count(x, y) > 0


def smoothed_polyline(
    x: ArrayLike,
    y: ArrayLike,
    sigma: float,
    samples: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The closed polygon ``(x, y)`` low-passed along its own arc length.

    ``x`` and ``y`` are an open vertex list normalised to the period, closed
    through ``(x[0] + 1, y[0])``. The polygon is resampled at ``samples``
    equal arc-length stations, and the periodic parts of both coordinates,
    ``x(s) - s / L`` and ``y(s)``, are convolved with a Gaussian of standard
    deviation ``sigma`` (in units of the period, measured along the arc) by
    multiplying their spectra by ``exp(-2 pi^2 sigma^2 f^2)``. A Gaussian
    rather than a sharp spectral cut because truncating the spectrum of a
    polygon rings, and the ringing is curvature the integral method would
    have to resolve. Parametrising by arc length rather than by ``x`` keeps
    an undercut profile representable.

    Returns the open list of ``samples`` smoothed points, ready for
    :func:`_polyline_boundary`.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if not sigma > 0.0:
        raise ValueError(f"sigma must be positive, got {sigma}")
    tx = np.append(x, x[0] + 1.0)
    ty = np.append(y, y[0])
    cumulative = np.concatenate(
        ([0.0], np.cumsum(np.hypot(np.diff(tx), np.diff(ty))))
    )
    total = float(cumulative[-1])
    if total <= 0.0:
        raise ValueError("profile has zero arc length")
    s = np.arange(samples) * (total / samples)
    drift = s / total
    xs = np.interp(s, cumulative, tx) - drift
    ys = np.interp(s, cumulative, ty)
    frequency = np.fft.rfftfreq(samples, d=total / samples)
    gain = np.exp(-2.0 * np.pi**2 * sigma**2 * frequency**2)
    xs = np.fft.irfft(np.fft.rfft(xs) * gain, n=samples) + drift
    ys = np.fft.irfft(np.fft.rfft(ys) * gain, n=samples)
    return xs, ys


#: A vertex turning by no more than this, in degrees, is collinear with its
#: neighbours to rounding: a side made only of such vertices is straight and
#: :func:`smoothed_polygon_sides` leaves it exactly as given.
_STRAIGHT_TURN_DEGREES = 1e-6


def polyline_rounded_turn(
    x: ArrayLike, y: ArrayLike, window: float
) -> float:
    """The most turning, in degrees, that vertices each below
    :data:`CORNER_TURN_DEGREES` carry within any stretch ``window`` long
    (in units of the period, along the arc) of the polygon ``(x, y)``.

    Large where a corner has been rounded into many slight bends, which is
    what averaging and resampling do to a measured valley: the turning is
    all there, but no single vertex crosses the threshold the graded mesh
    looks for.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    vx, vy, turn = _vertex_turns(np.append(x, x[0] + 1.0), np.append(y, y[0]))
    slight = np.where(turn <= CORNER_TURN_DEGREES, turn, 0.0)
    edge = np.hypot(np.diff(np.append(vx, vx[0] + 1.0)), np.diff(np.append(vy, vy[0])))
    at = np.concatenate(([0.0], np.cumsum(edge)))
    length, at = at[-1], at[:-1]
    # Two laps, so a window may run through the period wrap.
    laps = np.concatenate((at, at + length))
    running = np.concatenate(([0.0], np.cumsum(np.concatenate((slight, slight)))))
    reach = np.searchsorted(laps, at + window)
    return float((running[reach] - running[np.arange(len(at))]).max())


def polyline_has_bends(x: ArrayLike, y: ArrayLike) -> bool:
    """Whether the polygon ``(x, y)`` has a vertex that turns, but by no more
    than :data:`CORNER_TURN_DEGREES`: the slight corners a measured profile is
    made of, which no graded mesh resolves."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    _, _, turn = _vertex_turns(np.append(x, x[0] + 1.0), np.append(y, y[0]))
    slight = (turn > _STRAIGHT_TURN_DEGREES) & (turn <= CORNER_TURN_DEGREES)
    return bool(slight.any())


def smoothed_polygon_sides(
    x: ArrayLike,
    y: ArrayLike,
    sigma: float,
    samples: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.intp]]:
    """The polygon ``(x, y)`` with each side between its corners low-passed,
    and the corners kept.

    The corners are the vertices turning by more than
    :data:`CORNER_TURN_DEGREES`, found on the polygon as given. Each side
    between two of them is resampled at equal arc length and its departure
    from the chord joining its corners is convolved with the Gaussian of
    :func:`smoothed_polyline`. The departure is extended oddly about both
    ends before filtering, so it stays zero at the corners, which do not
    move, and the side keeps its slope into each corner. A side with no
    bend is returned as given. ``samples`` is spread over the sides in
    proportion to their length, with a few per side at least.

    Returns the open vertex list, starting at the polygon's first vertex as
    :func:`_polyline_boundary` expects, and the indices of the corners in
    it, for :func:`_graded_polyline`. With no corner this is
    :func:`smoothed_polyline`'s job, and it raises.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if not sigma > 0.0:
        raise ValueError(f"sigma must be positive, got {sigma}")
    vx, vy, turn = _vertex_turns(np.append(x, x[0] + 1.0), np.append(y, y[0]))
    m = len(vx)
    corners = np.flatnonzero(turn > CORNER_TURN_DEGREES)
    if corners.size == 0:
        raise ValueError("the polygon has no corner; use smoothed_polyline")

    # Two periods of vertices, so a side may run through the period wrap.
    big_x = np.concatenate((vx, vx + 1.0, [vx[0] + 2.0]))
    big_y = np.concatenate((vy, vy, [vy[0]]))
    big_turn = np.concatenate((turn, turn, [turn[0]]))
    edge = np.hypot(np.diff(big_x), np.diff(big_y))
    cumulative = np.concatenate(([0.0], np.cumsum(edge)))
    total = float(cumulative[m])

    pieces_x, pieces_y, starts, positions = [], [], [], []
    count = 0
    for first, last in zip(corners, np.append(corners[1:], corners[0] + m)):
        side_x, side_y = big_x[first : last + 1], big_y[first : last + 1]
        starts.append(count)
        if not (big_turn[first + 1 : last] > _STRAIGHT_TURN_DEGREES).any():
            pieces_x.append(side_x[:-1])
            pieces_y.append(side_y[:-1])
            positions.append(cumulative[first:last])
            count += len(side_x) - 1
            continue
        along = cumulative[first : last + 1] - cumulative[first]
        length = float(along[-1])
        stations = max(8, int(np.ceil(samples * length / total)))
        s = np.arange(stations + 1) * (length / stations)
        fraction = s / length
        smoothed = []
        for coordinate in (side_x, side_y):
            chord = coordinate[0] + (coordinate[-1] - coordinate[0]) * fraction
            departure = np.interp(s, along, coordinate) - chord
            odd = np.concatenate((departure[:-1], -departure[:0:-1]))
            frequency = np.fft.rfftfreq(2 * stations, d=length / stations)
            gain = np.exp(-2.0 * np.pi**2 * sigma**2 * frequency**2)
            filtered = np.fft.irfft(np.fft.rfft(odd) * gain, n=2 * stations)
            smoothed.append(chord + filtered[: stations + 1])
        pieces_x.append(np.concatenate(([side_x[0]], smoothed[0][1:-1])))
        pieces_y.append(np.concatenate(([side_y[0]], smoothed[1][1:-1])))
        positions.append(cumulative[first] + s[:-1])
        count += stations

    dense_x = np.concatenate(pieces_x)
    dense_y = np.concatenate(pieces_y)
    corner_at = np.array(starts, dtype=np.intp)
    # The sides start at the first corner; put the stretch past one period,
    # measured along the polygon as given, back at the front, so the list
    # starts where the polygon did.
    over = np.concatenate(positions) >= total
    dense_x[over] -= 1.0
    shift = int(over.sum())
    order = np.roll(np.arange(len(dense_x)), shift)
    corner_at = (corner_at + shift) % len(dense_x)
    return dense_x[order], dense_y[order], np.sort(corner_at)


def _graded_polyline(
    tx: NDArray[np.float64],
    ty: NDArray[np.float64],
    n: int,
    order: int,
    corners: NDArray[np.intp] | None = None,
) -> BoundaryCurve | None:
    """A corner-graded mesh on the closed polyline ``(tx, ty)``.

    Returns ``None`` when the polyline has no corner, so the caller keeps the
    equal-arc mesh. Otherwise each side between consecutive corners receives
    points in proportion to its length, placed at the midpoints of an equal
    split of the mesh parameter and mapped through :func:`_grading_map`, so
    no point sits on a corner and every point lies on a straight piece with
    an exact normal. The weights are the exact ``ds/dtau / n``.

    ``corners``, indices into the distinct vertices, names the corners
    instead of detecting them, which a smoothed polygon needs
    (:func:`smoothed_polygon_sides`): smoothing a side can move a corner's
    turn across the threshold, and its own dense vertices must not count.
    """
    vx, vy, turn = _vertex_turns(tx, ty)
    m = len(vx)
    if m < 2:
        return None
    if corners is None:
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
    corners: NDArray[np.intp] | None = None,
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
    ``corners`` names the corners instead (:func:`_graded_polyline`).
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
        graded = _graded_polyline(tx, ty, n, grading, corners)
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
