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

__all__ = [
    "REFLECTIVITY_MODELS",
    "ROUGHNESS_MODELS",
    "VISIBILITY_MODES",
    "ScalarOptionsState",
    "build_options",
]


#: Selectable values, in the order a form should offer them. The first is the
#: solver's own default in each case, so a freshly opened window and a bare
#: ``scalar.solve`` call agree.
REFLECTIVITY_MODELS = ("local", "average", "facet")
ROUGHNESS_MODELS = ("nevot-croce", "debye-waller", "none")
VISIBILITY_MODES = ("facet-normal", "horizon")


@dataclass(frozen=True, slots=True)
class ScalarOptionsState:
    """Scalar's own form fields, exactly as typed."""

    quadrature_points: str = "2048"
    reflectivity_model: str = REFLECTIVITY_MODELS[0]
    roughness_model: str = ROUGHNESS_MODELS[0]
    visibility: str = VISIBILITY_MODES[0]

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
    # Validated here rather than left to the solver so a bad value surfaces on
    # the field that produced it, alongside every other form error, instead of
    # as a ValueError from a worker thread.
    for name, allowed in (
        ("reflectivity_model", REFLECTIVITY_MODELS),
        ("roughness_model", ROUGHNESS_MODELS),
        ("visibility", VISIBILITY_MODES),
    ):
        if getattr(form, name) not in allowed:
            errors.append(
                FieldError(name, f"must be one of {', '.join(allowed)}")
            )
    # The one combination the solver refuses: with a coating, "facet" applies
    # one reflectivity to the whole groove and carries no per-point masks for
    # a horizon to narrow. Uncoated it is fine -- the model is inert and the
    # masks run at unit amplitude.
    if (
        form.reflectivity_model == "facet"
        and form.visibility == "horizon"
        and problem.coating is not None
    ):
        errors.append(
            FieldError(
                "visibility",
                "'horizon' cannot act under reflectivity_model='facet' with "
                "a coating: that model has no per-point masks for a horizon "
                "to narrow. Use 'local' or 'average', or keep 'facet-normal'.",
            )
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

    return {
        "quadrature_points": int(quadrature),
        "reflectivity_model": form.reflectivity_model,
        "roughness_model": form.roughness_model,
        "visibility": form.visibility,
    }
