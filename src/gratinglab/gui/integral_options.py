"""The integral solver's own form: boundary condition, points, and their guards.

The second solver's options module, made the way ``scalar_options.py`` said a
second one should be: its own file, its own frozen state, its own
``build_options`` raising :class:`~gratinglab.gui.state.FormErrors`, no
toolkit import anywhere.

Two knobs now, because the solver has two. ``conductivity`` chooses the
boundary condition -- ``"perfect"`` (relative efficiency, no material) or
``"tabulated"`` (the finite-conductivity system of Goray & Schmidt 2010,
absolute efficiency, absorption recorded). The mesh guard follows that
choice: the boundary field oscillates at the *reduced* wavelength
``lambda / sin(gamma)`` on the vacuum side, and under ``"tabulated"`` also at
``lambda / |sqrt(n^2 - cos^2 gamma)|` inside the material, which for anything
but a soft-X-ray index is the shorter of the two. Both floors mirror the
solver's own refusals (``solvers/integral``); validating here too means the
error lands on the field that caused it, alongside every other form error,
instead of surfacing as a ValueError from the worker thread.
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

__all__ = [
    "CONDUCTIVITY_MODES",
    "CONDUCTIVITY_NOTES",
    "IntegralOptionsState",
    "build_options",
]

#: Boundary sampling used only to estimate the arc length for the guards.
#: Coarse is fine -- arc length converges long before the solve does.
_ARC_SAMPLES = 512

#: Selectable boundary conditions, in the order a form should offer them. The
#: first is the solver's own default, so a freshly opened window and a bare
#: ``integral.solve`` call agree.
CONDUCTIVITY_MODES = ("perfect", "tabulated")

#: What each mode means, in the words the tab shows beneath the selector. Kept
#: here rather than in the widget for the reason `provenance.py` gives: what a
#: panel *claims about a result* is a correctness question, and the old static
#: label went stale the moment the solver grew a second boundary condition.
CONDUCTIVITY_NOTES = {
    "perfect": (
        "Perfectly conducting boundary: efficiencies are relative to a "
        "perfect reflector and sum to 1. A coating named in the geometry "
        "panel is not consulted (the provenance panel says so)."
    ),
    "tabulated": (
        "Finite conductivity (Goray & Schmidt 2010): the coating named in "
        "the geometry panel is read as a complex index at every wavelength. "
        "Efficiencies are absolute, absorption is recorded, and R + A = 1. "
        "Costs roughly twice as much per wavelength, on a finer mesh."
    ),
}


@dataclass(frozen=True, slots=True)
class IntegralOptionsState:
    """Integral's own form fields, exactly as typed."""

    boundary_points: str = "400"
    conductivity: str = CONDUCTIVITY_MODES[0]

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
    # Validated here rather than left to the solver so a bad value surfaces on
    # the field that produced it rather than as an UnsupportedConfiguration
    # from a worker thread.
    if form.conductivity not in CONDUCTIVITY_MODES:
        errors.append(
            FieldError(
                "conductivity",
                f"must be one of {', '.join(CONDUCTIVITY_MODES)}",
            )
        )
    if errors:
        raise FormErrors(tuple(errors))
    assert points is not None

    arc_length = (
        problem.profile.boundary(_ARC_SAMPLES).arc_length * problem.period
    )
    reduced = float(wavelengths.min()) / illumination.sin_gamma
    needed = int(np.ceil(_POINTS_PER_WAVELENGTH * arc_length / reduced))
    message = (
        f"needs at least {needed} to give "
        f"{_POINTS_PER_WAVELENGTH:g} nodes per transverse "
        f"wavelength ({reduced:.4g} nm) along the "
        f"{arc_length:.4g} nm boundary"
    )

    if form.conductivity == "tabulated":
        # Sequential prerequisites, not parallel checks: there is no index to
        # range-check without a material, and no metal-side floor without an
        # index. So these short-circuit rather than accumulating.
        name, indices, error = _tabulated_index(problem, wavelengths)
        if error is not None:
            raise FormErrors((error,))
        assert name is not None and indices is not None
        metal_needed, metal_message = _metal_floor(
            name, indices, illumination, wavelengths, arc_length
        )
        # One error on one field: whichever side of the boundary binds harder.
        # Reporting both would put two messages on `boundary_points` that a
        # single number answers.
        if metal_needed > needed:
            needed, message = metal_needed, metal_message

    if points < needed:
        raise FormErrors((FieldError("boundary_points", message),))

    return {"boundary_points": int(points), "conductivity": form.conductivity}


def _tabulated_index(
    problem: Problem, wavelengths: NDArray[np.float64]
) -> "tuple[str | None, NDArray[np.complex128] | None, FieldError | None]":
    """Resolve the material the tabulated boundary condition reads.

    Mirrors ``solvers.integral._resolve_index``, but returns a
    :class:`~gratinglab.gui.state.FieldError` instead of raising, so a missing
    or unusable material lands on a form field rather than arriving as a
    traceback from the worker thread.

    ``lookup`` is called without an :class:`UnknownMaterial` guard on purpose:
    :func:`gratinglab.gui.state.build` has already rejected any coating name
    that is not in the library, so that branch is unreachable from the window
    and a handler for it would be untested dead code.
    """
    from ..materials import lookup

    name = problem.coating or problem.substrate
    if name is None:
        return None, None, FieldError(
            "coating",
            'conductivity="tabulated" reads a complex index at every '
            "wavelength; name a coating in the geometry panel, or set "
            'conductivity to "perfect" (which needs no material).',
        )

    table = lookup(name)
    if not table.covers(wavelengths).all():
        low, high = table.range_nm
        return name, None, FieldError(
            "wavelength_start",
            f"{name} is tabulated over {low:.4g}-{high:.4g} nm; this scan "
            f"runs {wavelengths.min():.4g}-{wavelengths.max():.4g} nm. "
            "Extrapolating an optical constant is refused, so shorten the "
            "scan or install a table that covers it.",
        )

    indices = np.asarray(table.n(wavelengths), dtype=np.complex128)
    if (indices.imag < 0.0).any():
        return name, None, FieldError(
            "coating",
            f"{name} tabulates a negative absorption (gain) inside this "
            "scan; the transmission problem loses uniqueness there (Goray & "
            "Schmidt 2010, section 2.B).",
        )
    return name, indices, None


def _metal_floor(
    name: str,
    indices: NDArray[np.complex128],
    illumination: Illumination,
    wavelengths: NDArray[np.float64],
    arc_length: float,
) -> tuple[int, str]:
    """The mesh floor the material side imposes, and what to say about it.

    Same formula as the solver's own guard, so the two numbers agree: the
    field inside a conductive material oscillates and decays on the scale
    ``lambda / |sqrt(n^2 - cos^2 gamma)|``. For the shipped soft-X-ray tables
    ``|n|`` is within a hair of 1 and this barely moves; for a visible-light
    index of a few it is the binding constraint by a wide margin.
    """
    cos_gamma = float(np.sqrt(1.0 - illumination.sin_gamma**2))
    factors = np.abs(np.sqrt(indices**2 - cos_gamma**2))
    reduced = float((wavelengths / factors).min())
    needed = int(np.ceil(_POINTS_PER_WAVELENGTH * arc_length / reduced))
    return needed, (
        f"needs at least {needed} with coating {name}: the field inside a "
        f"conductive material oscillates and decays on the scale lambda/|n| "
        f"({reduced:.4g} nm here), and the mesh must resolve it along the "
        f'{arc_length:.4g} nm boundary. Raise it, or use conductivity="perfect".'
    )
