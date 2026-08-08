"""The scalar solver's own form: quadrature points, and their Nyquist guard.

Split out of :mod:`gratinglab.gui.state` when the window grew a tab per
solver. Geometry (period, profile, mount, wavelengths) is shared across every
solver a tab might run and lives in :class:`gratinglab.gui.state.FormState`;
this module holds the one thing that was never geometry to begin with --
``quadrature_points`` is scalar's own accuracy knob, and validating it needs
the *result* of parsing geometry (the propagating-order range) without being
part of parsing geometry itself. Keeping it here, rather than back inside
``state.build()``, means a second solver's own knob and guard get their own
module the same way, instead of `build()` growing a branch per solver.

Nothing here imports a toolkit; it is tested the same way `state.py` is,
without a window.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..geometry import order_range
from ..illumination import Illumination
from ..problem import Problem
from .state import FieldError, FormErrors, _number

__all__ = ["ScalarOptionsState", "build_options"]


@dataclass(frozen=True, slots=True)
class ScalarOptionsState:
    """Scalar's own form fields, exactly as typed."""

    quadrature_points: str = "2048"

    def with_field(self, name: str, value: Any) -> "ScalarOptionsState":
        """Return a copy with one field changed."""
        return replace(self, **{name: value})


def build_options(
    problem: Problem,
    illumination: Illumination,
    wavelengths: NDArray[np.float64],
    form: ScalarOptionsState,
) -> dict[str, Any]:
    """Parse scalar's options and check they can resolve the geometry's orders.

    Raises
    ------
    FormErrors
        Same contract as :func:`gratinglab.gui.state.build`: every problem
        found, not just the first.

    The Nyquist check needs a solved geometry (the propagating-order range,
    from the *shortest* wavelength in the scan, where the most orders
    propagate) to say anything, which is exactly why it cannot live inside
    geometry parsing -- geometry alone does not know which solver, or which
    solver option, is about to consume it.
    """
    errors: list[FieldError] = []
    quadrature = _number(
        form.quadrature_points, "quadrature_points", errors, minimum=16, integer=True
    )
    if errors:
        raise FormErrors(tuple(errors))
    assert quadrature is not None

    orders = order_range(
        float(wavelengths.min()), problem.period,
        illumination.sin_alpha, illumination.sin_gamma,
    )
    needed = 2 * int(np.abs(orders).max()) + 1
    if quadrature < needed:
        raise FormErrors(
            (
                FieldError(
                    "quadrature_points",
                    f"needs at least {needed} to resolve order "
                    f"{int(np.abs(orders).max())} at {wavelengths.min():g} nm",
                ),
            )
        )

    return {"quadrature_points": int(quadrature)}
