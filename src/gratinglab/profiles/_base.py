"""The machinery every concrete profile inherits: a numerical slope, slab
slicing and a meshed boundary, all built from a subclass's ``height``.

Profiles are frozen pydantic models, so there is no mutable state here.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pydantic import BaseModel, ConfigDict

from ._polyline import _polyline_boundary, _runs_to_intervals
from ._types import BoundaryCurve, Layer, ProfileRepresentationError


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
