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
        x=np.asarray(curve.t, dtype=np.float64) * period,
        y=np.asarray(curve.y, dtype=np.float64) * period,
        nx=np.asarray(curve.nx, dtype=np.float64),
        ny=np.asarray(curve.ny, dtype=np.float64),
        arc_length=float(curve.arc_length) * period,
        weights=(
            np.asarray(curve.ds, dtype=np.float64) * period if graded else None
        ),
        polygon=(
            None
            if polygon is None
            else tuple(np.asarray(v, dtype=np.float64) * period for v in polygon)
        ),
    )
