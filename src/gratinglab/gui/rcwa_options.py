"""The RCWA solver's own form: truncation, slices, and their guards.

Written as ``integral_options.py`` was: its own file, its own frozen state and
its own ``build_options`` raising :class:`~gratinglab.gui.state.FormErrors`,
importing no toolkit.

Two settings, because the solver has two. ``truncation`` is the accuracy knob
the convergence sweep varies. ``slices`` is the staircase's slab count for a
profile that is not lamellar; left blank it follows the truncation, which is
the solver's own default, so the one knob refines both.

The guards mirror the solver's refusals, so a bad setting is reported on the
field that caused it rather than as a ValueError from the worker thread: the
grating needs a material (RCWA has no perfectly conducting mode), the
material's table must cover the scan, it must not tabulate gain, and the
truncation must hold every order the scan needs. That last one calls the
solver's own :func:`~gratinglab.solvers.rcwa.required_truncation`, so the form
and the solver cannot disagree about it.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..illumination import Illumination
from ..problem import Problem
from ..profiles import Lamellar
from ..solvers.rcwa import _MIN_SLICES, required_truncation
from .state import FieldError, FormErrors, _number

__all__ = ["RcwaOptionsState", "build_options", "profile_note"]


@dataclass(frozen=True, slots=True)
class RcwaOptionsState:
    """RCWA's own form fields, exactly as typed. Blank ``slices`` means auto."""

    truncation: str = "20"
    slices: str = ""

    def with_field(self, name: str, value: Any) -> "RcwaOptionsState":
        """Return a copy with one field changed."""
        return replace(self, **{name: value})


def profile_note(problem: Problem | None, slices: str, truncation: str) -> str:
    """What the solver will do with this profile, in the words the tab shows.

    A claim about the result, so it is computed from the problem rather than
    written once beside a control that can falsify it (``integral_options``
    makes the same point about its boundary condition).
    """
    if problem is None:
        return ""
    if isinstance(problem.profile, Lamellar) and not problem.layers:
        return (
            "Lamellar: one exact slab, so slices are not used and the only "
            "approximation is the truncation."
        )
    try:
        count = int(slices) if slices.strip() else max(
            _MIN_SLICES, 2 * int(truncation)
        )
        how = f"{count} slabs"
    except ValueError:
        how = "a staircase of slabs"
    return (
        f"Not lamellar: the groove is replaced by {how}. In TM on a metal a "
        "staircase converges slowly (its corners carry field singularities "
        "the true profile lacks), and the provenance panel says so; sweep "
        "slices as well as truncation before trusting a TM result."
    )


def build_options(
    problem: Problem,
    illumination: Illumination,
    wavelengths: NDArray[np.float64],
    form: RcwaOptionsState,
) -> dict[str, Any]:
    """Parse the RCWA options and check them as the solver will.

    Raises
    ------
    FormErrors
        Same contract as :func:`gratinglab.gui.state.build`.
    """
    errors: list[FieldError] = []
    truncation = _number(form.truncation, "truncation", errors, minimum=1, integer=True)
    slices = None
    if form.slices.strip():
        slices = _number(form.slices, "slices", errors, minimum=1, integer=True)
    if errors:
        raise FormErrors(tuple(errors))
    assert truncation is not None

    # Sequential prerequisites, as in integral_options: no window check
    # without an index, and no index without a material.
    name, indices, error = _material_index(problem, wavelengths)
    if error is not None:
        raise FormErrors((error,))

    needed, reason = required_truncation(
        problem.period,
        illumination.sin_alpha,
        illumination.sin_gamma,
        wavelengths,
        indices,
        name,
    )
    if truncation < needed:
        raise FormErrors((FieldError(
            "truncation",
            f"needs at least {needed}: {reason}, and every order the scan "
            f"reports must sit inside the harmonic window -M..M.",
        ),))

    options: dict[str, Any] = {"truncation": int(truncation)}
    if slices is not None:
        options["slices"] = int(slices)
    return options


def _material_index(
    problem: Problem, wavelengths: NDArray[np.float64]
) -> "tuple[str | None, NDArray[np.complex128] | None, FieldError | None]":
    """The grating's material and its index over the scan, or the field at fault.

    Mirrors ``solvers.rcwa._resolve_index``. ``lookup`` needs no
    ``UnknownMaterial`` guard, because :func:`gratinglab.gui.state.build` has
    already rejected a coating that is not in the library.
    """
    from ..materials import lookup

    name = problem.coating or problem.substrate
    if name is None:
        return None, None, FieldError(
            "coating",
            "RCWA expands the grating's permittivity, so it needs a material: "
            "name a coating in the geometry panel. (A perfect conductor is the "
            'integral solver\'s conductivity="perfect".)',
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
            f"{name} tabulates a negative absorption (gain) inside this scan; "
            "a passive grating cannot be built from it.",
        )
    return name, indices, None
