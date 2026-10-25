"""The boundary curve in physical units, computed in a single place.

This mirrors ``Problem.height_nm``: the normalised profile is converted to nm
in exactly one place, so no other integral-method module multiplies by the
period.

Every quadrature in the integral method is written in a mesh parameter
``tau`` on ``[0, 1)`` whose nodes are equally spaced, with node ``j``
carrying the weight ``w_j = s'(tau_j) / P``. On the default equal-arc mesh
``tau`` is normalised arc length, so ``s' = L`` and every weight is the
spacing ``L / P``; on a corner-graded mesh the weights vary and are exact.
The operators are written once, in that form, and reduce to the equal-arc
formulas term by term.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...problem import Problem

__all__ = ["PhysicalBoundary", "physical_boundary"]


@dataclass(frozen=True, slots=True)
class PhysicalBoundary:
    """Boundary samples of one groove period, in nm.

    ``weights`` is ``None`` on the equal-arc mesh, where every node carries
    ``arc_length / len(x)`` (see :meth:`gratinglab.profiles.Profile.boundary`),
    and holds the per-node quadrature weights of a graded mesh otherwise.
    ``polygon`` is a graded mesh's open vertex list ``(x, y)`` in nm, on
    whose straight pieces every node lies, and ``None`` otherwise.
    """

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    nx: NDArray[np.float64]
    ny: NDArray[np.float64]
    arc_length: float
    weights: NDArray[np.float64] | None = None
    polygon: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None

    @property
    def graded(self) -> bool:
        return self.weights is not None

    @property
    def spacing(self) -> float:
        """The mean node spacing, which is every node's on the equal-arc mesh."""
        return self.arc_length / len(self.x)

    @property
    def quadrature(self) -> NDArray[np.float64]:
        """The weight ``w_j = s'(tau_j) / P`` of each node, in nm."""
        if self.weights is None:
            return np.full(len(self.x), self.spacing)
        return self.weights

    @property
    def speed(self) -> NDArray[np.float64]:
        """``ds/dtau`` at each node, in nm: ``L`` everywhere when equal-arc."""
        if self.weights is None:
            return np.full(len(self.x), self.arc_length)
        return self.weights * len(self.x)

    @property
    def resolution_length(self) -> float:
        """``P`` times the largest node spacing, in nm.

        The arc length on the equal-arc mesh. The resolution floor asks the
        coarsest spacing, not the mean, to resolve the field, and a graded
        mesh's coarsest spacing is about twice its mean.
        """
        if self.weights is None:
            return self.arc_length
        return float(self.weights.max()) * len(self.x)


def physical_boundary(
    problem: Problem, n: int, *, grading: int = 0
) -> PhysicalBoundary:
    """Sample ``problem.profile`` at ``n`` boundary points, in nm.

    ``grading`` is passed to the profile only when non-zero, so a profile
    written against the older ``boundary(n)`` signature still serves the
    equal-arc mesh.

    The profile parameter runs against the dispersion direction,
    ``t-hat = -d-hat`` (``docs/conventions.md``, "Direction of the profile
    parameter"), so the curve is reflected onto ``x = (1 - t) P``. The nodes
    are reversed with it, which keeps them advancing in ``+x`` as every
    operator here assumes, and ``nx`` changes sign. Through 0.1.12 this placed
    ``x = t P``, so every asymmetric groove was solved as its mirror image:
    a blazed grating sent its blaze to the wrong side of the normal while the
    scalar solver, on the same problem, sent it to the documented one.
    """
    curve = (
        problem.profile.boundary(n, grading=grading)
        if grading
        else problem.profile.boundary(n)
    )
    period = problem.period
    graded = bool(getattr(curve, "graded", False))
    polygon = getattr(curve, "polygon", None) if graded else None
    return PhysicalBoundary(
        x=(1.0 - np.asarray(curve.t, dtype=np.float64)[::-1]) * period,
        y=np.asarray(curve.y, dtype=np.float64)[::-1] * period,
        nx=-np.asarray(curve.nx, dtype=np.float64)[::-1],
        ny=np.asarray(curve.ny, dtype=np.float64)[::-1],
        arc_length=float(curve.arc_length) * period,
        weights=(
            np.asarray(curve.ds, dtype=np.float64)[::-1] * period if graded else None
        ),
        polygon=None if polygon is None else _mirrored_polygon(polygon, period),
    )


def _mirrored_polygon(
    polygon: tuple[NDArray[np.float64], NDArray[np.float64]], period: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The open vertex list ``(t, y)`` reflected onto ``x = (1 - t) P``.

    The list closes by continuing to ``(t[0] + 1, y[0])``, so reversing the
    traversal brings that closing segment to the front: the reflected list
    starts at the image of ``t[0] + 1``, which is ``-t[0]``, then runs
    ``1 - t`` back down to vertex 1, and closes on ``1 - t[0]``. Reversing
    the list alone would leave the first wall outside the three periods the
    corner subtraction integrates, offsetting every node's correction.
    """
    t = np.asarray(polygon[0], dtype=np.float64)
    y = np.asarray(polygon[1], dtype=np.float64)
    x = np.concatenate(([-t[0]], 1.0 - t[:0:-1]))
    return x * period, np.concatenate(([y[0]], y[:0:-1])) * period
