"""The C-method's own form: boundary condition, truncation, and their guards.

The third options module, written to the pattern of ``integral_options.py``:
its own file, its own frozen state and its own ``build_options`` raising
:class:`~gratinglab.gui.state.FormErrors`, and no toolkit import.

The form has two settings, because the solver has two. ``conductivity`` is
named and offered as the integral solver offers it, so one geometry reaches
either method with the same choice: ``"perfect"`` (relative efficiency, no
material) or ``"tabulated"`` (one interface into the coating, read per
wavelength, with absolute efficiency, absorption and, into a lossless
material, the transmitted orders). ``truncation`` is the accuracy knob: ``N``
harmonics either side of the centre of the propagating orders.

The guards mirror the solver's refusals so an error lands on the field that
caused it instead of arriving as a ValueError from the worker thread. The
truncation floors are not re-derived here: they are
:func:`gratinglab.solvers.cmethod.truncation_floor` (the orders) and
:func:`~gratinglab.solvers.cmethod.wall_floor` (a steep wall), the functions
the solver itself applies, so the form and the solver cannot disagree.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..illumination import Illumination
from ..problem import Problem
from ..solvers.base import UnsupportedConfiguration
from ..solvers.cmethod import _MAX_TRUNCATION, truncation_floor, wall_floor
from .integral_options import _tabulated_index
from .state import FieldError, FormErrors, _number

__all__ = [
    "CONDUCTIVITY_MODES",
    "CONDUCTIVITY_NOTES",
    "CMethodOptionsState",
    "build_options",
]

#: Selectable boundary conditions, the solver's default first.
CONDUCTIVITY_MODES = ("perfect", "tabulated")

#: What each mode means, shown beneath the selector, for the reason
#: ``integral_options.CONDUCTIVITY_NOTES`` gives: a claim about the result
#: must not be hardcoded beside a control that can falsify it.
CONDUCTIVITY_NOTES = {
    "perfect": (
        "Perfect conductor in any mount: efficiencies are relative to a "
        "perfect reflector and sum to 1. A coating named in the geometry "
        "panel is not consulted (the provenance panel says so)."
    ),
    "tabulated": (
        "One interface into the coating named in the geometry panel, read "
        "as a complex index at every wavelength. Efficiencies are absolute "
        "and absorption is recorded; into a lossless material the "
        "transmitted orders are recorded too, and R + T + A = 1 is a check."
    ),
}


@dataclass(frozen=True, slots=True)
class CMethodOptionsState:
    """The C-method's own form fields, exactly as typed."""

    truncation: str = "32"
    conductivity: str = CONDUCTIVITY_MODES[0]

    def with_field(self, name: str, value: Any) -> "CMethodOptionsState":
        """Return a copy with one field changed."""
        return replace(self, **{name: value})


def build_options(
    problem: Problem,
    illumination: Illumination,
    wavelengths: NDArray[np.float64],
    form: CMethodOptionsState,
) -> dict[str, Any]:
    """Parse the C-method options and check the truncation holds every order.

    Raises
    ------
    FormErrors
        Same contract as :func:`gratinglab.gui.state.build`.
    """
    errors: list[FieldError] = []
    truncation = _number(
        form.truncation, "truncation", errors, minimum=1, integer=True
    )
    if truncation is not None and truncation > _MAX_TRUNCATION:
        errors.append(
            FieldError(
                "truncation",
                f"at most {_MAX_TRUNCATION}: the eigenproblem is dense and its "
                "cost grows as N^3, so a sweep that has not settled by then is "
                "the answer",
            )
        )
    if form.conductivity not in CONDUCTIVITY_MODES:
        errors.append(
            FieldError(
                "conductivity",
                f"must be one of {', '.join(CONDUCTIVITY_MODES)}",
            )
        )
    if errors:
        raise FormErrors(tuple(errors))
    assert truncation is not None

    indices = None
    if form.conductivity == "tabulated":
        _name, indices, error = _tabulated_index(problem, wavelengths)
        if error is not None:
            raise FormErrors((error,))

    floor = truncation_floor(problem, illumination, wavelengths, indices)
    if truncation < floor.needed:
        raise FormErrors((FieldError("truncation", floor.message),))
    # A steep wall's floor goes on the field too, but only while a truncation
    # the field accepts can meet it. Past the ceiling, or on a profile the
    # method cannot describe at all, the solver refuses the geometry with
    # UnsupportedConfiguration, and the tab names that refusal instead.
    try:
        wall = wall_floor(problem, illumination, wavelengths, indices)
    except UnsupportedConfiguration:
        wall = None
    if wall is not None and truncation < wall.needed <= _MAX_TRUNCATION:
        raise FormErrors((FieldError("truncation", wall.message),))

    return {"truncation": int(truncation), "conductivity": form.conductivity}
