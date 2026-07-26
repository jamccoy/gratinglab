"""Groove profiles, in the four representations the solver methods need.

**Profiles are normalised.** Both the coordinate ``t`` and the returned height
are expressed in units of the period, so a profile is a pure *shape* and the
period lives once, on :class:`~gratinglab.problem.Problem`. A solver that needs
physical height multiplies by ``problem.period``. This also matches PCGrate's
``.ggp`` format, which stores ``Period: 1`` and normalised points.

The four methods want genuinely different geometry from the same profile:

=================  ==========================================================
Method             Representation
=================  ==========================================================
Scalar             :meth:`Profile.height` -- the groove phase function
RCWA               :meth:`Profile.slice_layers` -- piecewise-constant slabs
C-method           :meth:`Profile.slope` -- smooth dy/dt; no vertical facets
Integral           :meth:`Profile.boundary` -- parametrised curve with
                   outward normals and arc-length quadrature weights
=================  ==========================================================

Designing all four now is deliberate. Retrofitting normals and arc lengths onto
a shipped API once the integral method arrives would churn every profile type.

**A profile that cannot supply a representation raises.** ``Blazed`` with a
vertical anti-blaze facet has no finite slope there, so :meth:`Blazed.slope`
refuses rather than returning a large number that would quietly poison a
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

    Raised rather than approximated: a C-method solve fed a silently smoothed
    vertical facet produces a plausible, wrong answer.
    """


@dataclass(frozen=True, slots=True)
class Layer:
    """One piecewise-constant slab of an RCWA discretisation.

    ``intervals`` lists the ``(t_start, t_end)`` spans, normalised to the
    period, where *grating material* occupies this slab. The complement is the
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
        Outward unit normal, pointing into the incident medium (``ny > 0``).
    ds
        Arc-length weight for each point, for quadrature. Sums to the total
        profile arc length over one period.
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]
    nx: NDArray[np.float64]
    ny: NDArray[np.float64]
    ds: NDArray[np.float64]

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

    def boundary(self, n: int) -> BoundaryCurve:
        """Parametrised boundary with outward normals and arc-length weights."""

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

    def boundary(self, n: int) -> BoundaryCurve:
        """Uniformly sampled boundary with outward normals and arc-length weights.

        Normals point into the incident medium (``ny > 0``). Weights use the
        trapezoidal arc length around each point, so ``ds.sum()`` is the profile
        length over one period.
        """
        if n < 3:
            raise ValueError(f"need at least 3 boundary points, got {n}")

        t = np.linspace(0.0, 1.0, n, endpoint=False)
        y = self.height(t)

        # Periodic central differences for the tangent.
        dt = 1.0 / n
        dy = (np.roll(y, -1) - np.roll(y, 1)) / (2 * dt)
        norm = np.hypot(1.0, dy)

        # Tangent (1, dy)/norm rotated -90 degrees gives outward (dy, -1);
        # flip so the normal points away from the material, into +y.
        nx = -dy / norm
        ny = np.ones_like(dy) / norm

        # Trapezoidal arc length: half the chord to each neighbour.
        forward = np.hypot(dt, np.roll(y, -1) - y)
        backward = np.hypot(dt, y - np.roll(y, 1))
        ds = 0.5 * (forward + backward)

        return BoundaryCurve(t=t, y=y, nx=nx, ny=ny, ds=ds)


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


class Sinusoidal(_BaseProfile):
    """Holographic (sinusoidal) groove -- smooth, so C-method friendly."""

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
    """A measured profile -- AFM, or PCGrate's ``.ggp`` polygonal border.

    Points are normalised: ``t`` spans one period on ``[0, 1]`` and ``y`` is in
    units of the period. Use :meth:`from_measurements` to normalise raw nm data.

    The profile is stored as given and interpolated linearly, matching how
    PCGrate treats a polygonal border, so a scalar calculation and a PCGrate
    run can be driven from the identical geometry.

    **Undercut profiles are representable.** If ``t`` doubles back, the boundary
    is not a function of position and :meth:`height` and
    :meth:`~_BaseProfile.slice_layers` refuse -- but :meth:`boundary` still
    works, because an integral-equation solver parametrises the curve rather
    than assuming a height function. That asymmetry is real physics, not a
    limitation of this class: RCWA and scalar theory genuinely cannot represent
    an undercut, and the integral method genuinely can.
    """

    t: tuple[float, ...]
    y: tuple[float, ...]

    @model_validator(mode="after")
    def _check_points(self) -> "FromProfileData":
        if len(self.t) != len(self.y):
            raise ValueError(f"t has {len(self.t)} points, y has {len(self.y)}")
        if len(self.t) < 3:
            raise ValueError(f"need at least 3 points, got {len(self.t)}")
        t = np.asarray(self.t)
        if t.min() < 0.0 or t.max() > 1.0:
            raise ValueError(f"t must lie in [0, 1], got [{t.min()}, {t.max()}]")
        return self

    @classmethod
    def from_measurements(
        cls, x_nm: ArrayLike, y_nm: ArrayLike, period_nm: float
    ) -> "FromProfileData":
        """Normalise a measured profile in nm, referencing its minimum to zero."""
        x = np.asarray(x_nm, dtype=np.float64)
        y = np.asarray(y_nm, dtype=np.float64)
        if period_nm <= 0:
            raise ValueError(f"period_nm must be positive, got {period_nm}")
        return cls(
            t=tuple((x - x.min()) / period_nm),
            y=tuple((y - y.min()) / period_nm),
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

    def boundary(self, n: int) -> BoundaryCurve:
        """Resample the stored polygon at equal arc length.

        Overrides the base implementation, which assumes a height function.
        Working from the stored points directly is what lets this handle an
        undercut boundary.
        """
        if n < 3:
            raise ValueError(f"need at least 3 boundary points, got {n}")

        # Close the polygon: the curve continues into the next period.
        tx = np.append(np.asarray(self.t, dtype=np.float64), self.t[0] + 1.0)
        ty = np.append(np.asarray(self.y, dtype=np.float64), self.y[0])

        segment = np.hypot(np.diff(tx), np.diff(ty))
        cumulative = np.concatenate(([0.0], np.cumsum(segment)))
        total = cumulative[-1]
        if total <= 0:
            raise ValueError("profile has zero arc length")

        stations = np.linspace(0.0, total, n, endpoint=False)
        x = np.interp(stations, cumulative, tx)
        y = np.interp(stations, cumulative, ty)

        # Periodic central differences; x wraps by exactly one period.
        dx = np.roll(x, -1) - np.roll(x, 1)
        dx[0] += 1.0
        dx[-1] += 1.0
        dy = np.roll(y, -1) - np.roll(y, 1)
        norm = np.hypot(dx, dy)

        # Tangent (dx, dy) rotated +90 deg is (-dy, dx), which points into the
        # incident medium for a curve traversed in the +t direction.
        forward = np.hypot(np.roll(x, -1) - x, np.roll(y, -1) - y)
        forward[-1] = np.hypot(x[0] + 1.0 - x[-1], y[0] - y[-1])
        backward = np.roll(forward, 1)

        return BoundaryCurve(
            t=x,
            y=y,
            nx=-dy / norm,
            ny=dx / norm,
            ds=0.5 * (forward + backward),
        )
