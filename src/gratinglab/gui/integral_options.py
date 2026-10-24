"""The integral solver's own form: boundary condition, points, and their guards.

This is the options module for the second solver, written as
``scalar_options.py`` recommended for a second one: it has its own file, its
own frozen state and its own ``build_options`` raising
:class:`~gratinglab.gui.state.FormErrors`, and it imports no toolkit.

The form has three settings, because the solver has three. ``conductivity``
chooses the boundary condition, either ``"perfect"`` (relative efficiency, no
material) or ``"tabulated"`` (the finite-conductivity system of Goray &
Schmidt 2010, with absolute efficiency and absorption recorded). The mesh
guard depends on that choice. The boundary field oscillates at the reduced
wavelength ``lambda / sin(gamma)`` on the vacuum side, and under
``"tabulated"`` also at
``lambda / |sqrt(n^2 - cos^2 gamma)|`` inside the material, which is the
shorter of the two for any index other than a soft-X-ray one. Both floors
mirror the solver's own refusals (``solvers/integral``). Validating here as
well places the error on the field that caused it, alongside every other form
error, instead of letting it surface as a ValueError from the worker thread.

``corner_grading`` chooses the mesh: equal arc length, or graded towards the
profile's corners at an order from 2 to 6. A graded mesh changes the floor,
because its coarsest spacing (mid-side) is about twice the mean, so the
guard is computed on the mesh the solver will actually build, and the
solver's own check for nodes too close to separate is applied here too.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..illumination import Illumination
from ..problem import Problem
from ..solvers.integral import (
    _MAX_GRADING,
    _POINTS_PER_WAVELENGTH,
    _coincident_nodes,
)
from ..solvers.integral._boundary import physical_boundary
from .state import FieldError, FormErrors, _number

__all__ = [
    "CONDUCTIVITY_MODES",
    "CONDUCTIVITY_NOTES",
    "GRADING_CHOICES",
    "GRADING_NOTES",
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


#: The corner-mesh choices as the form shows them: ``off`` is the equal-arc
#: mesh, the solver's default, and a number is the grading order.
GRADING_CHOICES = ("off",) + tuple(str(p) for p in range(2, _MAX_GRADING + 1))

#: What the choice means, shown beneath the selector, for the reason
#: ``CONDUCTIVITY_NOTES`` gives.
GRADING_NOTES = {
    "off": (
        "Equal-arc mesh. On a profile with corners TM converges only at "
        "first order (with finite conductivity, both polarizations do); "
        "the provenance panel's energy-balance deviation shows how far."
    ),
    "graded": (
        "Nodes clustered at the profile's corners (order {p}), where the "
        "field is singular: cornered profiles converge far faster, and a "
        "profile with no corner is unaffected. Needs about twice the "
        "boundary points for the same resolution floor. 4 is recommended."
    ),
}


@dataclass(frozen=True, slots=True)
class IntegralOptionsState:
    """Integral's own form fields, exactly as typed."""

    boundary_points: str = "400"
    conductivity: str = CONDUCTIVITY_MODES[0]
    corner_grading: str = GRADING_CHOICES[0]

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
    if form.corner_grading not in GRADING_CHOICES:
        errors.append(
            FieldError(
                "corner_grading",
                f"must be one of {', '.join(GRADING_CHOICES)}",
            )
        )
    if errors:
        raise FormErrors(tuple(errors))
    assert points is not None
    grading = 0 if form.corner_grading == "off" else int(form.corner_grading)

    arc_length = (
        problem.profile.boundary(_ARC_SAMPLES).arc_length * problem.period
    )
    # The floor is set by the coarsest spacing. On the equal-arc mesh that is
    # the mean, arc_length / P; a graded mesh is built here exactly as the
    # solver will build it, since its coarsest spacing depends on P.
    resolution_length = arc_length
    if grading:
        try:
            mesh = physical_boundary(problem, int(points), grading=grading)
        except ValueError as refused:  # too few points for the corners
            raise FormErrors((FieldError("boundary_points", str(refused)),)) from None
        coincident = _coincident_nodes(mesh, problem.period)
        if coincident is not None:
            raise FormErrors((FieldError(
                "corner_grading",
                f"order {grading} at {int(points)} points puts two nodes "
                f"{coincident[0]:.3g} nm apart, closer than double precision "
                "can separate on this groove; choose a lower order (4 is "
                "recommended) or fewer points.",
            ),))
        if mesh.graded:
            resolution_length = mesh.resolution_length
    reduced = float(wavelengths.min()) / illumination.sin_gamma
    needed = int(np.ceil(_POINTS_PER_WAVELENGTH * resolution_length / reduced))
    message = (
        f"needs at least {needed} to give "
        f"{_POINTS_PER_WAVELENGTH:g} nodes per transverse "
        f"wavelength ({reduced:.4g} nm) along the "
        f"{arc_length:.4g} nm boundary"
        + (", where the graded mesh's coarsest spacing sets the floor"
           if resolution_length != arc_length else "")
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
            name, indices, illumination, wavelengths, resolution_length
        )
        # One error on one field: whichever side of the boundary binds harder.
        # Reporting both would put two messages on `boundary_points` that a
        # single number answers.
        if metal_needed > needed:
            needed, message = metal_needed, metal_message

    if points < needed:
        raise FormErrors((FieldError("boundary_points", message),))

    return {
        "boundary_points": int(points),
        "conductivity": form.conductivity,
        "corner_grading": grading,
    }


def _tabulated_index(
    problem: Problem, wavelengths: NDArray[np.float64]
) -> "tuple[str | None, NDArray[np.complex128] | None, FieldError | None]":
    """Resolve the material the tabulated boundary condition reads.

    Mirrors ``solvers.integral._resolve_index``, but returns a
    :class:`~gratinglab.gui.state.FieldError` instead of raising, so a missing
    or unusable material is reported on a form field instead of as a
    traceback from the worker thread.

    ``lookup`` is called without an :class:`UnknownMaterial` guard because
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

    It uses the same formula as the solver's own guard, so the two numbers
    agree. The field inside a conductive material oscillates and decays on the
    scale ``lambda / |sqrt(n^2 - cos^2 gamma)|``. For the shipped soft-X-ray
    tables ``|n|`` is very close to 1 and this floor changes little; for a
    visible-light index of a few it is by far the binding constraint.
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
