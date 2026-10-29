"""The C-method's own solver options: boundary condition, truncation, guards.

Headless, the discipline of ``test_gui_integral_options.py``: the third
solver's options module, tested without a toolkit.
"""

import inspect

import numpy as np
import pytest

from gratinglab.gui.cmethod_options import (
    CONDUCTIVITY_MODES,
    CONDUCTIVITY_NOTES,
    CMethodOptionsState,
    build_options,
)
from gratinglab.gui.state import FormErrors, FormState, build
from gratinglab.solvers.cmethod import _MAX_TRUNCATION, cmethod, truncation_floor, wall_floor

#: A scan short enough to solve in a test, on a smooth profile.
SMOOTH = dict(
    period="500", profile_kind="Sinusoidal", depth_fraction="0.2",
    mount="Classical", alpha="10",
    wavelength_start="450", wavelength_stop="600", wavelength_count="2",
)


def _options(form=FormState(), **fields):
    parsed = build(form)
    return build_options(
        parsed.problem, parsed.illumination, parsed.wavelengths,
        CMethodOptionsState(**fields),
    )


def _field_errors(form=FormState(), **fields):
    with pytest.raises(FormErrors) as excinfo:
        _options(form, **fields)
    return excinfo.value.errors


class TestOptions:
    def test_passes_solver_options_through(self):
        assert _options() == {"truncation": 32, "conductivity": "perfect"}

    def test_the_defaults_are_the_solvers_own(self):
        """A freshly opened tab and a bare `cmethod.solve` call agree."""
        parameters = inspect.signature(cmethod.solve).parameters
        defaults = CMethodOptionsState()
        assert int(defaults.truncation) == parameters["truncation"].default
        assert defaults.conductivity == parameters["conductivity"].default

    def test_the_modes_are_the_integral_solvers_names(self):
        """One geometry and one choice reach either rigorous method."""
        from gratinglab.gui import integral_options

        assert CONDUCTIVITY_MODES == integral_options.CONDUCTIVITY_MODES

    def test_every_mode_has_a_note(self):
        assert set(CONDUCTIVITY_NOTES) == set(CONDUCTIVITY_MODES)

    def test_an_unknown_mode_is_a_field_error(self):
        errors = _field_errors(conductivity="lossy")
        assert errors[0].field == "conductivity"


class TestTruncationRange:
    @pytest.mark.parametrize("raw", ["0", "-3", "abc", "2.5", ""])
    def test_not_a_positive_integer_is_a_field_error(self, raw):
        assert _field_errors(truncation=raw)[0].field == "truncation"

    def test_above_the_solvers_ceiling_is_a_field_error(self):
        errors = _field_errors(truncation=str(_MAX_TRUNCATION + 1))
        assert errors[0].field == "truncation"
        assert str(_MAX_TRUNCATION) in errors[0].message

    def test_the_ceiling_itself_is_accepted(self):
        assert _options(truncation=str(_MAX_TRUNCATION))["truncation"] == _MAX_TRUNCATION

    def test_both_bad_fields_are_reported_together(self):
        errors = _field_errors(truncation="abc", conductivity="lossy")
        assert {e.field for e in errors} == {"truncation", "conductivity"}


class TestTruncationFloor:
    def test_too_few_harmonics_for_the_orders_lands_on_the_field(self):
        """The default soft X-ray off-plane scan carries orders -4 to 11."""
        parsed = build(FormState())
        floor = truncation_floor(parsed.problem, parsed.illumination, parsed.wavelengths)
        errors = _field_errors(truncation=str(floor.needed - 1))
        assert errors[0].field == "truncation"
        assert f"need at least {floor.needed}" in errors[0].message

    def test_the_floor_is_accepted(self):
        """On the default scan the 70.5 deg facet binds, not the orders."""
        parsed = build(FormState())
        floor = max(
            truncation_floor(parsed.problem, parsed.illumination, parsed.wavelengths).needed,
            wall_floor(parsed.problem, parsed.illumination, parsed.wavelengths).needed,
        )
        assert _options(truncation=str(floor))["truncation"] == floor

    def test_the_form_and_the_solver_agree_at_the_floor(self):
        """The form uses the solver's own function, so what the form accepts
        the solver solves, and one below it the solver refuses too."""
        parsed = build(FormState(**SMOOTH))
        floor = truncation_floor(parsed.problem, parsed.illumination, parsed.wavelengths)
        options = _options(FormState(**SMOOTH), truncation=str(floor.needed))
        cmethod.solve(parsed.problem, parsed.illumination, parsed.wavelengths, **options)
        with pytest.raises(ValueError, match=f"need at least {floor.needed}"):
            cmethod.solve(
                parsed.problem, parsed.illumination, parsed.wavelengths,
                truncation=floor.needed - 1,
            )

    def test_a_lossless_material_raises_the_floor_by_its_transmitted_orders(self):
        """Into a high-index dielectric (n = 3) more orders transmit than
        reflect, and the form names them."""
        parsed = build(FormState(**SMOOTH))
        glass = np.full(len(parsed.wavelengths), 3.0 + 0j)
        vacuum = truncation_floor(parsed.problem, parsed.illumination, parsed.wavelengths)
        into = truncation_floor(
            parsed.problem, parsed.illumination, parsed.wavelengths, glass
        )
        assert into.needed > vacuum.needed
        assert "transmitted" in into.message


class TestSteepWallFloor:
    """Issue #7: the solver's wall floor, on the field while it can be met."""

    def test_too_few_harmonics_for_a_steep_wall_lands_on_the_field(self):
        """The default scan starts at 1 nm, where its 70.5 deg facet needs 29
        harmonics; 24 was 2.5e-2 from the integral solver there."""
        parsed = build(FormState())
        wall = wall_floor(parsed.problem, parsed.illumination, parsed.wavelengths)
        assert wall.needed == 29
        errors = _field_errors(truncation=str(wall.needed - 1))
        assert errors[0].field == "truncation"
        assert f"needs about {wall.needed} harmonics" in errors[0].message

    def test_a_wall_beyond_the_ceiling_is_left_to_the_solver(self):
        """No value of the field could meet it, so the form does not ask for
        one: the solver refuses the geometry, and the tab says why."""
        form = FormState(antiblaze_angle="89.0", wavelength_count="2")
        parsed = build(form)
        assert wall_floor(
            parsed.problem, parsed.illumination, parsed.wavelengths
        ).needed > _MAX_TRUNCATION
        options = _options(form, truncation="64")
        from gratinglab.solvers import UnsupportedConfiguration

        with pytest.raises(UnsupportedConfiguration, match="too steep"):
            cmethod.solve(parsed.problem, parsed.illumination, parsed.wavelengths, **options)


class TestTabulated:
    def test_without_a_coating_is_refused_on_the_coating_field(self):
        errors = _field_errors(conductivity="tabulated")
        assert errors[0].field == "coating"

    def test_a_scan_outside_the_table_is_refused_on_the_wavelengths(self):
        form = FormState(coating="Au", wavelength_start="0.001", wavelength_stop="0.002",
                         wavelength_count="2")
        errors = _field_errors(form, conductivity="tabulated")
        assert errors[0].field == "wavelength_start"

    def test_a_covered_coating_passes(self):
        form = FormState(coating="Au", wavelength_start="4.0", wavelength_stop="4.5",
                         wavelength_count="2")
        assert _options(form, conductivity="tabulated", truncation="16") == {
            "truncation": 16, "conductivity": "tabulated",
        }
