"""Integral's own solver options: boundary points, and their sampling guard.

Headless, same discipline as ``test_gui_scalar_options.py`` -- the second
solver's options module, exercised the way the first one's docstring promised
a second one would be.
"""

import pytest

from gratinglab.gui.integral_options import IntegralOptionsState, build_options
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
        assert options == {"boundary_points": 400}

    def test_the_default_is_the_solvers_own(self):
        import inspect

        signature = inspect.signature(integral.solve)
        assert (
            int(IntegralOptionsState().boundary_points)
            == signature.parameters["boundary_points"].default
        )

    def test_with_field_does_not_mutate(self):
        original = IntegralOptionsState()
        changed = original.with_field("boundary_points", "800")
        assert original.boundary_points == "400"
        assert changed.boundary_points == "800"


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
