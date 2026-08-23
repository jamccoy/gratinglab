"""The boundary curve in physical units, once.

Mirrors ``Problem.height_nm``: the normalised profile is converted to nm in
exactly one place, so no other integral-method module multiplies by the
period.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...problem import Problem

__all__ = ["PhysicalBoundary", "physical_boundary"]


@dataclass(frozen=True, slots=True)
class PhysicalBoundary:
    """Equal-arc-length boundary samples of one groove period, in nm.

    ``spacing`` is the quadrature weight of the rectangular rule,
    ``arc_length / len(x)`` -- the nodes are equally spaced in arc length by
    construction (see :meth:`gratinglab.profiles.Profile.boundary`).
    """

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    nx: NDArray[np.float64]
    ny: NDArray[np.float64]
    arc_length: float

    @property
    def spacing(self) -> float:
        return self.arc_length / len(self.x)


def physical_boundary(problem: Problem, n: int) -> PhysicalBoundary:
    """Sample ``problem.profile`` at ``n`` boundary points, in nm."""
    curve = problem.profile.boundary(n)
    period = problem.period
    return PhysicalBoundary(
        x=np.asarray(curve.t, dtype=np.float64) * period,
        y=np.asarray(curve.y, dtype=np.float64) * period,
        nx=np.asarray(curve.nx, dtype=np.float64),
        ny=np.asarray(curve.ny, dtype=np.float64),
        arc_length=float(curve.arc_length) * period,
    )
