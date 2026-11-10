""":class:`Trapezoidal`, a binary groove with optional rounded corners, and
``_Piece``, one section of its outline.

A frozen pydantic model, built by callers and read by the solvers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pydantic import Field, model_validator

from ._base import _BaseProfile
from ._polyline import _polyline_boundary
from ._types import BoundaryCurve, ProfileRepresentationError


#: Floor on the radicand of an arc's derivative. The derivative is finite
#: everywhere on a fillet of a trapezoid, because the arc never reaches the
#: point where its tangent turns vertical unless the sidewall is itself
#: vertical -- and that case refuses before it gets here.
_TINY = 1e-300


@dataclass(frozen=True, slots=True)
class _Piece:
    """One x-interval of a trapezoidal groove: a straight run or a fillet arc.

    Straight pieces are ``y = y0 + gradient * (x - x0)``. Arc pieces are
    ``y = cy + bulge * sqrt(r**2 - (x - cx)**2)`` with ``bulge`` +1 where the
    centre lies below the surface (a convex crest corner) and -1 where it lies
    above (a concave floor corner).
    """

    x0: float
    x1: float
    y0: float = 0.0
    gradient: float = 0.0
    cx: float = 0.0
    cy: float = 0.0
    r: float = 0.0
    bulge: float = 0.0

    @property
    def is_arc(self) -> bool:
        return self.bulge != 0.0

    def height(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        if not self.is_arc:
            return self.y0 + self.gradient * (x - self.x0)
        # Clipped at zero: a sample exactly on the tangent point can land a
        # rounding error outside the disc, and a NaN there would propagate
        # silently through an efficiency calculation.
        dx = x - self.cx
        return self.cy + self.bulge * np.sqrt(np.maximum(self.r**2 - dx**2, 0.0))

    def slope(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        if not self.is_arc:
            return np.full_like(x, self.gradient)
        dx = x - self.cx
        return -self.bulge * dx / np.sqrt(np.maximum(self.r**2 - dx**2, _TINY))


class Trapezoidal(_BaseProfile):
    """Lamellar groove with sloped sidewalls and rounded corners.

    The shape an etched grating has, as opposed to the rectangle
    :class:`Lamellar` describes. An etch leaves the trench wider at its mouth
    than at its floor, and the crest and the floor each meet the sidewall on a
    radius rather than at a point.

    ``sidewall_angle = 90`` and ``corner_radius = 0`` recover :class:`Lamellar`
    exactly, up to a translation: the ridge here is centred on ``t = 0.5`` and
    the floor straddles the period boundary, which puts all four corners inside
    one period so each carries a whole fillet. Translating a profile shifts the
    phase of every order equally and leaves the efficiencies unchanged.

    Undercut sidewalls are deliberately out of range: ``sidewall_angle`` stops
    at 90 degrees, and a re-entrant groove is a measured shape rather than a
    designed one, so it belongs to :class:`FromProfileData`, whose boundary may
    double back.

    >>> t = Trapezoidal(depth_fraction=0.43, duty_cycle=0.36, sidewall_angle=76.0)
    >>> round(t.ridge_base_width, 4)
    0.5744
    >>> round(float(t.height(0.5)), 4)
    0.43
    """

    depth_fraction: float = Field(gt=0.0, description="Depth / period")
    duty_cycle: float = Field(
        default=0.5,
        gt=0.0,
        lt=1.0,
        description="Fraction of the period at full height, measured at the crest",
    )
    sidewall_angle: float = Field(
        default=90.0,
        gt=0.0,
        le=90.0,
        description="Sidewall angle from the horizontal, deg; 90 is vertical",
    )
    corner_radius: float = Field(
        default=0.0,
        ge=0.0,
        description="Fillet radius at all four corners, in units of the period",
    )

    @model_validator(mode="after")
    def _check_fits_in_one_period(self) -> "Trapezoidal":
        if self.floor_width <= 0.0:
            raise ValueError(
                f"the sidewalls consume the whole period: duty_cycle "
                f"{self.duty_cycle} plus two runs of {self.wall_run:.4g} leaves "
                f"no floor. Reduce the depth or the duty cycle, steepen the "
                f"sidewall, or use Blazed for a groove with no floor at all"
            )
        largest = self.max_corner_radius
        if self.corner_radius > largest * (1.0 + 1e-12):
            raise ValueError(
                f"corner_radius {self.corner_radius} is too large for this "
                f"groove; adjacent fillets would overlap. The largest that fits "
                f"is {largest:.6g}"
            )
        return self

    @property
    def depth(self) -> float:
        return self.depth_fraction

    @property
    def wall_run(self) -> float:
        """Horizontal distance a sidewall covers, normalised. Zero when vertical."""
        if self.has_vertical_sidewall:
            return 0.0
        return float(self.depth / np.tan(np.radians(self.sidewall_angle)))

    @property
    def wall_length(self) -> float:
        """Slant length of one sidewall, normalised."""
        return float(self.depth / np.sin(np.radians(self.sidewall_angle)))

    @property
    def has_vertical_sidewall(self) -> bool:
        return bool(np.isclose(self.sidewall_angle, 90.0))

    @property
    def ridge_top_width(self) -> float:
        """Width of the crest, normalised. The same quantity as ``duty_cycle``."""
        return self.duty_cycle

    @property
    def ridge_base_width(self) -> float:
        """Width of the ridge where it meets the floor, normalised.

        With :attr:`ridge_top_width` this is the pair a cross-section measures
        directly, one at each end of the sidewall.
        """
        return self.duty_cycle + 2.0 * self.wall_run

    @property
    def floor_width(self) -> float:
        """Width of the trench floor, normalised."""
        return 1.0 - self.ridge_base_width

    @property
    def max_corner_radius(self) -> float:
        """Largest fillet the three flats and the sidewall leave room for.

        Every corner of a trapezoid turns by the sidewall angle, so each
        consumes the same tangent length ``R * tan(phi / 2)`` from the sides
        meeting there. The binding side is whichever of the floor, the crest
        and the sidewall is shortest; each is shared by two fillets, so each
        may give up half its length.
        """
        shortest = min(self.floor_width, self.duty_cycle, self.wall_length)
        half_angle = np.radians(self.sidewall_angle) / 2.0
        return float(0.5 * shortest / np.tan(half_angle))

    def _corners(self) -> tuple[tuple[float, float], ...]:
        """The four vertices of the unrounded trapezoid, left to right."""
        run, half_floor = self.wall_run, 0.5 * self.floor_width
        d = self.depth
        return (
            (half_floor, 0.0),
            (half_floor + run, d),
            (half_floor + run + self.duty_cycle, d),
            (1.0 - half_floor, 0.0),
        )

    def _fillet(
        self, before: tuple[float, float], at: tuple[float, float],
        after: tuple[float, float],
    ) -> tuple[_Piece, tuple[float, float], tuple[float, float]]:
        """The arc replacing vertex ``at``, and its two tangent points.

        Standard construction: the fillet touches both sides at the tangent
        length ``R / tan(theta)`` from the vertex, where ``2 theta`` is the
        angle between them, and its centre sits ``R / sin(theta)`` from the
        vertex along the bisector.
        """
        px, py = at
        u1 = np.array(before, dtype=np.float64) - np.array(at, dtype=np.float64)
        u2 = np.array(after, dtype=np.float64) - np.array(at, dtype=np.float64)
        u1 /= np.hypot(*u1)
        u2 /= np.hypot(*u2)

        half_angle = 0.5 * np.arccos(np.clip(float(u1 @ u2), -1.0, 1.0))
        radius = self.corner_radius
        tangent = radius / np.tan(half_angle)

        bisector = u1 + u2
        bisector /= np.hypot(*bisector)
        cx, cy = np.array([px, py]) + bisector * (radius / np.sin(half_angle))

        start = (px + u1[0] * tangent, py + u1[1] * tangent)
        end = (px + u2[0] * tangent, py + u2[1] * tangent)
        bulge = -1.0 if cy > py else 1.0
        arc = _Piece(
            x0=start[0], x1=end[0], cx=float(cx), cy=float(cy),
            r=radius, bulge=bulge,
        )
        return arc, start, end

    def _pieces(self) -> tuple[_Piece, ...]:
        """The profile over ``[0, 1]`` as consecutive x-intervals, left to right.

        A vertical sidewall contributes a zero-width interval, which the
        lookup in :meth:`height` never selects for an interior point. The
        height genuinely jumps there, exactly as it does for
        :class:`Lamellar`.
        """
        corners = self._corners()

        if self.corner_radius == 0.0:
            nodes = ((0.0, 0.0), *corners, (1.0, 0.0))
            return tuple(
                _Piece(
                    x0=a[0], x1=b[0], y0=a[1],
                    gradient=0.0 if b[0] == a[0] else (b[1] - a[1]) / (b[0] - a[0]),
                )
                for a, b in zip(nodes[:-1], nodes[1:])
            )

        # Each corner's neighbours are the corners either side, with the floor
        # continuing into the previous and next periods at the two ends.
        neighbours = (
            ((corners[3][0] - 1.0, 0.0), corners[1]),
            (corners[0], corners[2]),
            (corners[1], corners[3]),
            (corners[2], (corners[0][0] + 1.0, 0.0)),
        )

        arcs, starts, ends = [], [], []
        for corner, (before, after) in zip(corners, neighbours):
            arc, start, end = self._fillet(before, corner, after)
            arcs.append(arc)
            starts.append(start)
            ends.append(end)

        pieces = [_Piece(x0=0.0, x1=starts[0][0], y0=0.0)]
        for i, arc in enumerate(arcs):
            pieces.append(arc)
            if i < 3:
                a, b = ends[i], starts[i + 1]
                pieces.append(
                    _Piece(
                        x0=a[0], x1=b[0], y0=a[1],
                        gradient=0.0 if b[0] == a[0] else (b[1] - a[1]) / (b[0] - a[0]),
                    )
                )
        pieces.append(_Piece(x0=ends[3][0], x1=1.0, y0=0.0))
        return tuple(pieces)

    def _evaluate(self, t: ArrayLike, derivative: bool) -> NDArray[np.float64]:
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        pieces = self._pieces()
        edges = np.array([p.x1 for p in pieces])
        # side="right" makes the profile right-continuous at a discontinuity,
        # so a vertical sidewall takes the height of the interval beginning
        # there. That is the convention `Lamellar.height` already follows, and
        # the two agree point for point because of it.
        index = np.clip(np.searchsorted(edges, t, side="right"), 0, len(pieces) - 1)
        out = np.empty_like(t)
        for i, piece in enumerate(pieces):
            here = index == i
            if here.any():
                x = t[here]
                out[here] = piece.slope(x) if derivative else piece.height(x)
        return out

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        return self._evaluate(t, derivative=False)

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        """dy/dt, analytic on the flats, the sidewalls and the fillets alike.

        A vertical sidewall has no finite slope and refuses, as
        :class:`Lamellar` does. Rounding the corners does not rescue it: the
        wall between the two fillets is still vertical.
        """
        if self.has_vertical_sidewall:
            raise ProfileRepresentationError(
                "Trapezoidal with sidewall_angle=90 has vertical sidewalls and "
                "no finite slope there. The C-method cannot represent it; use a "
                "slightly shallower sidewall_angle, or a method that handles "
                "discontinuities (RCWA, integral)."
            )
        return self._evaluate(t, derivative=True)

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Exact polygon when unrounded; the fillets sampled finely when not.

        Without rounding the four corners are genuine vertices, which
        ``grading`` clusters points at. With rounding there are no corners and
        no field singularity to resolve, so a graded mesh is neither offered
        nor needed, and :func:`_polyline_boundary` returns the equal-arc mesh.
        """
        corners = self._corners()
        if self.corner_radius == 0.0:
            x = np.array([0.0, *(c[0] for c in corners)])
            y = np.array([0.0, *(c[1] for c in corners)])
            return _polyline_boundary(x, y, n, grading=grading)

        xs: list[float] = [0.0]
        ys: list[float] = [0.0]
        for piece in self._pieces():
            if not piece.is_arc:
                # A run the fillets consume entirely (at max_corner_radius)
                # has zero width, and its end is the previous arc's end again
                # to within rounding. Appending it leaves a segment of ~1e-17
                # pointing anywhere, which the graded mesh took for a corner
                # and clustered two nodes onto (measured 2026-10-29: the
                # integral solver refused every boundary_points at 72 deg).
                if piece.x1 > piece.x0:
                    xs.append(piece.x1)
                    ys.append(float(piece.height(np.array([piece.x1]))[0]))
                continue
            # One point per degree of sweep, and never fewer than eight: the
            # chord sagitta is then below 1e-5 of the radius, well under the
            # resampling the integral method applies afterwards.
            sweep = max(8, int(np.ceil(self.sidewall_angle)))
            a0 = np.arctan2(
                float(piece.height(np.array([piece.x0]))[0]) - piece.cy,
                piece.x0 - piece.cx,
            )
            a1 = np.arctan2(
                float(piece.height(np.array([piece.x1]))[0]) - piece.cy,
                piece.x1 - piece.cx,
            )
            angles = np.linspace(a0, a1, sweep + 1)[1:]
            xs.extend(piece.cx + piece.r * np.cos(angles))
            ys.extend(piece.cy + piece.r * np.sin(angles))
        # The closing segment of the polygon runs from the last point to
        # (x[0] + 1, y[0]), which is the rest of the floor.
        return _polyline_boundary(
            np.array(xs[:-1]), np.array(ys[:-1]), n, grading=grading
        )
