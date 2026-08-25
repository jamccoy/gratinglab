"""Integral's own solver options: boundary condition, points, their guards.

Headless, same discipline as ``test_gui_scalar_options.py`` -- the second
solver's options module, exercised the way the first one's docstring promised
a second one would be.
"""

import re

import pytest

from gratinglab.gui.integral_options import (
    CONDUCTIVITY_MODES,
    CONDUCTIVITY_NOTES,
    IntegralOptionsState,
    build_options,
)
from gratinglab.gui.state import FormErrors, FormState, build
from gratinglab.solvers import integral


class TestSamplingGuard:
    def test_catches_a_mesh_too_coarse_for_the_reduced_wavelength(self):
        """The default geometry is soft X-ray off-plane; 32 points cannot
        resolve the transverse wavelength and the error must land on the
        field, with the number needed in the message."""
        parsed = build(FormState())
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(boundary_points="32"),
            )
        errors = excinfo.value.errors
        assert errors[0].field == "boundary_points"
        assert "at least" in errors[0].message

    def test_a_valid_form_never_trips_the_solvers_own_guard(self):
        """If build_options() accepts it, solve() must not reject it."""
        parsed = build(
            FormState(
                period="600", profile_kind="Sinusoidal", depth_fraction="0.3",
                mount="Classical", alpha="25",
                wavelength_start="500", wavelength_stop="600",
                wavelength_count="2",
            )
        )
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            IntegralOptionsState(boundary_points="64"),
        )
        integral.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )

    def test_a_non_numeric_count_is_a_field_error_not_a_traceback(self):
        parsed = build(FormState())
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(boundary_points="abc"),
            )
        assert excinfo.value.errors[0].field == "boundary_points"


class TestOptions:
    def test_passes_solver_options_through(self):
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            IntegralOptionsState(boundary_points="400"),
        )
        assert options == {"boundary_points": 400, "conductivity": "perfect"}

    def test_the_default_is_the_solvers_own(self):
        import inspect

        signature = inspect.signature(integral.solve)
        assert (
            int(IntegralOptionsState().boundary_points)
            == signature.parameters["boundary_points"].default
        )

    def test_the_conductivity_default_is_the_solvers_own(self):
        import inspect

        signature = inspect.signature(integral.solve)
        assert (
            IntegralOptionsState().conductivity
            == signature.parameters["conductivity"].default
        )

    def test_the_offered_modes_are_exactly_the_solvers_own(self):
        """A mode the form offers that the solver does not implement would be
        an UnsupportedConfiguration from a worker thread; one it implements
        but the form omits is unreachable. Read them off the signature so the
        tuple cannot drift."""
        import typing

        hints = typing.get_type_hints(integral.solve)
        assert set(CONDUCTIVITY_MODES) == set(typing.get_args(hints["conductivity"]))

    def test_every_mode_has_a_note(self):
        assert set(CONDUCTIVITY_NOTES) == set(CONDUCTIVITY_MODES)

    def test_with_field_does_not_mutate(self):
        original = IntegralOptionsState()
        changed = original.with_field("boundary_points", "800")
        assert original.boundary_points == "400"
        assert changed.boundary_points == "800"


class TestTabulatedConductivity:
    """The finite-conductivity mode's refusals, lifted onto form fields.

    Every one of these is a refusal the solver would make anyway. The point of
    repeating it here is *where it lands*: on the field that caused it, next to
    every other form error, rather than as a traceback from the worker thread
    after the user has waited for a solve to start.
    """

    def test_an_unknown_mode_is_a_field_error(self):
        parsed = build(FormState())
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(conductivity="leontovich"),
            )
        assert excinfo.value.errors[0].field == "conductivity"

    def test_without_a_coating_the_error_names_the_coating_field(self):
        """There is no rigorous answer without an index, and the form says so
        on the field that would supply one -- naming the way out of it as
        well, since "perfect" needs no material."""
        parsed = build(FormState())
        assert parsed.problem.coating is None
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(conductivity="tabulated"),
            )
        error = excinfo.value.errors[0]
        assert error.field == "coating"
        assert "perfect" in error.message

    def test_a_scan_outside_the_table_names_the_range(self):
        """Extrapolating an optical constant is refused, not silently done."""
        parsed = build(
            FormState(
                coating="Au", wavelength_start="400", wavelength_stop="700",
                wavelength_count="3",
            )
        )
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(conductivity="tabulated"),
            )
        error = excinfo.value.errors[0]
        assert error.field == "wavelength_start"
        assert "Au is tabulated over" in error.message

    def test_the_metal_side_mesh_floor_matches_the_solvers_own(self):
        """The guard exists to pre-empt the solver's refusal, so the number it
        demands has to be the number the solver demands -- one more point than
        it asks for must actually solve, and one fewer must be refused by
        both."""
        parsed = build(FormState(coating="Au"))
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(
                    boundary_points="32", conductivity="tabulated"
                ),
            )
        error = excinfo.value.errors[0]
        assert error.field == "boundary_points"
        needed = int(re.search(r"at least (\d+)", error.message).group(1))
        assert "Au" in error.message

        # One below the stated floor: both refuse.
        with pytest.raises(FormErrors):
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                IntegralOptionsState(
                    boundary_points=str(needed - 1), conductivity="tabulated"
                ),
            )
        with pytest.raises(ValueError):
            integral.solve(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                boundary_points=needed - 1, conductivity="tabulated",
            )
        # At the floor the form accepts, and the solver's own guard agrees --
        # checked on a single wavelength so the assertion is about the guard
        # rather than about waiting for a scan.
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            IntegralOptionsState(
                boundary_points=str(needed), conductivity="tabulated"
            ),
        )
        assert options["boundary_points"] == needed

    def test_a_valid_form_never_trips_the_solvers_own_guard(self):
        """The same property the perfect boundary carries, for the coupled
        system: if build_options() accepts it, solve() must not reject it.

        Two wavelengths and the coarsest accepted mesh, because this is a
        statement about the guards agreeing, not about the physics."""
        parsed = build(
            FormState(
                coating="Au", wavelength_start="4.0", wavelength_stop="4.5",
                wavelength_count="2",
            )
        )
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            IntegralOptionsState(boundary_points="128", conductivity="tabulated"),
        )
        scan = integral.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        # The mode really ran: only the tabulated condition records absorption.
        assert scan.absorption is not None
        assert scan.provenance.notes["material"] == "Au"


class TestPurity:
    def test_imports_no_toolkit(self):
        import ast
        import inspect

        import gratinglab.gui.integral_options as module

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(module))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6", "matplotlib"} & imported
