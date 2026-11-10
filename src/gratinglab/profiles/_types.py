"""What a groove profile is: the :class:`Profile` protocol, and the slabs and
boundary curves it hands to the solvers.

No state. Everything else in this package imports from here, and this module
imports nothing from the package.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike, NDArray


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
