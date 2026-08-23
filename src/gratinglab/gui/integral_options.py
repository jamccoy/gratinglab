"""The integral solver's own form: boundary points, and their sampling guard.

The second solver's options module, made the way ``scalar_options.py`` said a
second one should be: its own file, its own frozen state, its own
``build_options`` raising :class:`~gratinglab.gui.state.FormErrors`, no
toolkit import anywhere.

The guard mirrors the solver's own refusal (``solvers/integral``): the
boundary field oscillates at the *reduced* wavelength ``lambda / sin(gamma)``,
and a mesh below ~6 nodes per oscillation cannot resolve it. Validating here
too means the error lands on the field that caused it, alongside every other
form error, instead of surfacing as a ValueError from the worker thread.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..illumination import Illumination
from ..problem import Problem
from ..solvers.integral import _POINTS_PER_WAVELENGTH
from .state import FieldError, FormErrors, _number

__all__ = ["IntegralOptionsState", "build_options"]

#: Boundary sampling used only to estimate the arc length for the guard.
#: Coarse is fine -- arc length converges long before the solve does.
_ARC_SAMPLES = 512


@dataclass(frozen=True, slots=True)
class IntegralOptionsState:
    """Integral's own form fields, exactly as typed."""

    boundary_points: str = "400"

    def with_field(self, name: str, value: Any) -> "IntegralOptionsState":
        """Return a copy with one field changed."""
        return replace(self, **{name: value})


def build_options(
    problem: Problem,
    illumination: Illumination,
    wavelengths: NDArray[np.float64],
    form: IntegralOptionsState,
) -> dict[str, Any]:
    """Parse the integral options and check the mesh can resolve the field.

    Raises
    ------
    FormErrors
        Same contract as :func:`gratinglab.gui.state.build`.
    """
    errors: list[FieldError] = []
    points = _number(
        form.boundary_points, "boundary_points", errors, minimum=32, integer=True
    )
    if errors:
        raise FormErrors(tuple(errors))
    assert points is not None

    arc_length = (
        problem.profile.boundary(_ARC_SAMPLES).arc_length * problem.period
    )
    reduced = float(wavelengths.min()) / illumination.sin_gamma
    needed = int(np.ceil(_POINTS_PER_WAVELENGTH * arc_length / reduced))
    if points < needed:
        raise FormErrors(
            (
                FieldError(
                    "boundary_points",
                    f"needs at least {needed} to give "
                    f"{_POINTS_PER_WAVELENGTH:g} nodes per transverse "
                    f"wavelength ({reduced:.4g} nm) along the "
                    f"{arc_length:.4g} nm boundary",
                ),
            )
        )

    return {"boundary_points": int(points)}
