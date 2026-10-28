"""RCWA's own solver options: truncation, slices, and their guards.

Headless, the discipline of ``test_gui_integral_options.py``: the third
solver's options module, which imports no toolkit and is tested without one.
The material is the bundled Au table, which is what the window offers.
"""

import pytest

from gratinglab.gui.rcwa_options import RcwaOptionsState, build_options, profile_note
from gratinglab.gui.state import FormErrors, FormState, build
from gratinglab.solvers import rcwa

#: Soft X-ray, off-plane, lamellar, two wavelengths: small enough to solve in a
#: test, and inside the Au table's range.
FORM = dict(
    profile_kind="Lamellar", coating="Au",
    wavelength_start="4.0", wavelength_stop="4.5", wavelength_count="2",
)


def parsed(**overrides):
    return build(FormState(**{**FORM, **overrides}))


def options(form: RcwaOptionsState, **overrides):
    p = parsed(**overrides)
    return build_options(p.problem, p.illumination, p.wavelengths, form)


class TestTheGuards:
    def test_a_valid_form_never_trips_the_solvers_own_guard(self):
        """If build_options() accepts it, solve() must not reject it."""
        p = parsed()
        chosen = build_options(
            p.problem, p.illumination, p.wavelengths, RcwaOptionsState(truncation="10")
        )
        rcwa.solve(p.problem, p.illumination, p.wavelengths, **chosen)

    def test_a_truncation_short_of_an_order_lands_on_the_field(self):
        """In-plane at 4 nm on a 315 nm period, orders out to about 78
        propagate; the message names the number needed and the order that
        sets it, from the same function the solver refuses with."""
        with pytest.raises(FormErrors) as caught:
            options(RcwaOptionsState(truncation="10"), mount="Classical", alpha="10")
        (error,) = caught.value.errors
        assert error.field == "truncation"
        assert "needs at least" in error.message and "propagates" in error.message
        needed = int(error.message.split("needs at least ")[1].split(":")[0])
        options(RcwaOptionsState(truncation=str(needed)), mount="Classical", alpha="10")
        with pytest.raises(ValueError, match=f"need at least {needed}"):
            p = parsed(mount="Classical", alpha="10")
            rcwa.solve(p.problem, p.illumination, p.wavelengths, truncation=needed - 1)

    def test_no_material_is_refused_on_the_coating_field(self):
        with pytest.raises(FormErrors) as caught:
            options(RcwaOptionsState(), coating="")
        (error,) = caught.value.errors
        assert error.field == "coating" and "needs a material" in error.message

    def test_a_scan_outside_the_table_is_refused_on_the_wavelength(self):
        with pytest.raises(FormErrors) as caught:
            options(RcwaOptionsState(), wavelength_start="0.1", wavelength_stop="0.2")
        assert caught.value.errors[0].field == "wavelength_start"

    @pytest.mark.parametrize("field", ["truncation", "slices"])
    @pytest.mark.parametrize("value", ["0", "-3", "2.5", "many"])
    def test_bad_numbers_are_field_errors(self, field, value):
        with pytest.raises(FormErrors) as caught:
            options(RcwaOptionsState(**{field: value}))
        assert {e.field for e in caught.value.errors} == {field}


class TestWhatItPasses:
    def test_blank_slices_is_left_to_the_solver(self):
        assert options(RcwaOptionsState(truncation="10")) == {"truncation": 10}

    def test_a_number_pins_them(self):
        chosen = options(RcwaOptionsState(truncation="10", slices="24"))
        assert chosen == {"truncation": 10, "slices": 24}


class TestTheProfileNote:
    def test_lamellar_is_exact(self):
        assert "exact" in profile_note(parsed().problem, "", "20")

    def test_anything_else_is_a_staircase_with_its_count(self):
        problem = parsed(profile_kind="Blazed").problem
        assert "40 slabs" in profile_note(problem, "", "20")  # max(8, 2M)
        assert "12 slabs" in profile_note(problem, "12", "20")
        assert "TM" in profile_note(problem, "", "20")

    def test_a_half_typed_number_does_not_break_it(self):
        problem = parsed(profile_kind="Blazed").problem
        assert "staircase" in profile_note(problem, "", "lots")
