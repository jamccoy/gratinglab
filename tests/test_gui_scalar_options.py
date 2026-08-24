"""Scalar's own solver options: quadrature points, and their Nyquist guard.

Headless, same discipline as `test_gui_state.py` -- this is exactly the logic
that used to live inside `state.build()` before a second solver's own options
proved it did not belong there. See `gratinglab.gui.scalar_options`.
"""

import pytest

from gratinglab.gui.scalar_options import ScalarOptionsState, build_options
from gratinglab.gui.state import FormErrors, FormState, build
from gratinglab.solvers import scalar


class TestNyquistGuard:
    def test_catches_insufficient_quadrature_before_the_solver_does(self):
        """Turns a solver traceback into an actionable field error."""
        parsed = build(
            FormState(
                period="1400", wavelength_start="20", wavelength_stop="30",
                mount="Classical", alpha="10",
            )
        )
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                ScalarOptionsState(quadrature_points="16"),
            )
        errors = excinfo.value.errors
        assert any(e.field == "quadrature_points" for e in errors)
        assert "at least" in errors[0].message

    def test_a_valid_form_never_trips_the_solvers_own_nyquist_check(self):
        """If build_options() accepts it, solve() must not reject it."""
        parsed = build(
            FormState(
                period="1400", wavelength_start="400", wavelength_stop="700",
                mount="Classical", alpha="10",
            )
        )
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(quadrature_points="4096"),
        )
        scalar.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )

    def test_a_non_numeric_quadrature_is_a_field_error_not_a_traceback(self):
        parsed = build(FormState())
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                ScalarOptionsState(quadrature_points="abc"),
            )
        assert excinfo.value.errors[0].field == "quadrature_points"


class TestOptions:
    def test_passes_solver_options_through(self):
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(quadrature_points="4096"),
        )
        assert options == {
            "quadrature_points": 4096,
            "reflectivity_model": "local",
            "roughness_model": "nevot-croce",
            "visibility": "facet-normal",
        }

    @pytest.mark.parametrize(
        "field,value",
        [
            ("reflectivity_model", "mean-surface"),
            ("roughness_model", "gaussian"),
            ("visibility", "shadow-map"),
        ],
    )
    def test_a_model_the_solver_does_not_know_is_a_form_error(self, field, value):
        """Caught on the field that produced it, alongside every other form
        error, rather than surfacing later as a ValueError from the worker."""
        parsed = build(FormState())
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                ScalarOptionsState(**{field: value}),
            )
        assert excinfo.value.errors[0].field == field

    def test_facet_plus_horizon_with_a_coating_is_a_form_error(self):
        """The one combination the solver refuses -- caught on the field that
        produced it, same as an unknown model, rather than surfacing as a
        ValueError from the worker."""
        parsed = build(FormState(coating="Au"))
        with pytest.raises(FormErrors) as excinfo:
            build_options(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                ScalarOptionsState(
                    reflectivity_model="facet", visibility="horizon"
                ),
            )
        error = excinfo.value.errors[0]
        assert error.field == "visibility"
        assert "facet" in error.message

    def test_facet_plus_horizon_without_a_coating_is_accepted(self):
        """Uncoated, the facet model is inert and the masks run at unit
        amplitude -- if build_options() accepts it, solve() must not reject."""
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(reflectivity_model="facet", visibility="horizon"),
        )
        scalar.solve(parsed.problem, parsed.illumination, [2.4], **options)

    def test_horizon_with_a_coating_and_a_resolved_model_reaches_the_solver(self):
        """The guard refuses only the facet model: coated horizon runs under
        the default per-point reflectivity."""
        parsed = build(FormState(coating="Au"))
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(visibility="horizon"),
        )
        scalar.solve(parsed.problem, parsed.illumination, [2.4], **options)

    def test_the_defaults_are_the_solvers_own(self):
        """A freshly opened window and a bare `scalar.solve` must agree, or the
        GUI quietly answers a different question from the API."""
        import inspect

        signature = inspect.signature(scalar.solve)
        defaults = ScalarOptionsState()
        for name in ("reflectivity_model", "roughness_model", "visibility"):
            assert getattr(defaults, name) == signature.parameters[name].default

    def test_options_are_exactly_what_the_solver_accepts(self):
        """A renamed solver keyword would break the GUI silently otherwise."""
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(),
        )
        scalar.solve(parsed.problem, parsed.illumination, [2.4], **options)

    def test_with_field_does_not_mutate(self):
        original = ScalarOptionsState()
        changed = original.with_field("quadrature_points", "8192")
        assert original.quadrature_points == "2048"
        assert changed.quadrature_points == "8192"


class TestPurity:
    def test_imports_no_toolkit(self):
        import ast
        import inspect

        import gratinglab.gui.scalar_options as module

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(module))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6", "matplotlib"} & imported
