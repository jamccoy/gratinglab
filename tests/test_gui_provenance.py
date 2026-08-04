r"""What the provenance panel says.

Headless, and that is the point: which line carries which colour is a
correctness question -- it decides whether a user reads a valid run as broken
-- so it is tested here rather than by walking tag ranges in a live window.
"""

import numpy as np
import pytest

from gratinglab.checks import check_energy_balance
from gratinglab.gui.provenance import (
    ALARM_TAGS,
    TAG_COLORS,
    Line,
    about_text,
    alarming_text,
    error_lines,
    provenance_lines,
    solving_lines,
    to_html,
)
from gratinglab.gui.state import FieldError, FormState, build
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.result import EfficiencyScan, Provenance
from gratinglab.solvers import scalar

UNPOL = "unpolarized"


@pytest.fixture(scope="module")
def default_run():
    """The scan a user gets by opening the app and pressing nothing."""
    parsed = build(FormState())
    scan = scalar.solve(
        parsed.problem, parsed.illumination, parsed.wavelengths, **parsed.options
    )
    return scan, check_energy_balance(scan), parsed.lambda_over_period


def make_scan(*, converged=None, warnings=(), notes=None, totals=(0.5, 0.6)):
    """A minimal scan with exactly the provenance under test."""
    values = np.asarray(totals, dtype=float)
    return EfficiencyScan(
        wavelengths=np.arange(1.0, len(values) + 1.0),
        orders=np.array([0, 1]),
        efficiency=np.column_stack([values * 0.6, values * 0.4]),
        propagating=np.ones((len(values), 2), dtype=bool),
        provenance=Provenance(
            "scalar",
            version="test",
            truncation=2048,
            converged=converged,
            wall_time_s=0.071,
            warnings=tuple(warnings),
            notes=dict(notes or {}),
        ),
    )


def lines_for(scan, lambda_over_period=0.0159, **kwargs):
    return provenance_lines(
        scan, check_energy_balance(scan), lambda_over_period, **kwargs
    )


class TestNothingCorrectLooksWrong:
    """The regression this module was extracted to protect.

    A run in its default, correct mode used to render as a list of problems:
    `converged is None` and a coating-free relative normalization were tagged
    identically to a real validity breach.
    """

    def test_not_yet_checked_is_not_an_alarm(self):
        lines = lines_for(make_scan(converged=None))
        assert "not yet checked" in "".join(l.text for l in lines)
        assert "not yet checked" not in alarming_text(lines)

    def test_relative_normalization_is_not_an_alarm(self):
        lines = lines_for(make_scan(notes={"normalization": "relative"}))
        assert "normalization: " in "".join(l.text for l in lines)
        assert "coating" not in alarming_text(lines)

    def test_a_demonstrated_failure_to_converge_IS_an_alarm(self):
        """The other half of the rule -- the check must still be able to fail,
        or dimming everything would have been the same fix."""
        lines = lines_for(make_scan(converged=False))
        assert "NO" in alarming_text(lines)

    def test_convergence_success_reads_as_success(self):
        lines = lines_for(make_scan(converged=True))
        assert any(l.tag == "ok" and l.text.strip() == "yes" for l in lines)


class TestWarnings:
    def test_every_solver_warning_is_shown_verbatim(self):
        """Never summarised, never truncated: the warnings are why the panel
        exists."""
        message = "λ/period is 0.6, past the scalar validity bound"
        lines = lines_for(make_scan(warnings=[message]))
        assert message in alarming_text(lines)

    def test_warnings_keep_their_order(self):
        lines = lines_for(make_scan(warnings=["first", "second", "third"]))
        text = "".join(l.text for l in lines)
        assert text.index("first") < text.index("second") < text.index("third")

    def test_the_scalar_energy_warning_reaches_the_panel(self, default_run):
        """End to end on the real default case: the solver raises it, the
        panel shows it, and it is styled as a warning."""
        scan, energy, ratio = default_run
        lines = provenance_lines(scan, energy, ratio)
        assert "summed efficiency" in alarming_text(lines)


class TestEnergyLine:
    def test_a_deficit_reads_as_a_pass(self):
        lines = lines_for(make_scan(totals=(0.4, 0.7)))
        energy = next(l for l in lines if "Σ ∈" in l.text)
        assert energy.tag == "ok"
        assert "EXCEEDS UNITY" not in energy.text

    def test_an_excess_reads_as_a_failure(self):
        lines = lines_for(make_scan(totals=(0.9, 1.5)))
        energy = next(l for l in lines if "Σ ∈" in l.text)
        assert energy.tag == "bad"
        assert "EXCEEDS UNITY" in energy.text

    def test_the_range_is_reported_not_a_single_number(self):
        """Σ varies across the scan; one number would hide that."""
        lines = lines_for(make_scan(totals=(0.4, 0.9)))
        assert "Σ ∈ [0.4000, 0.9000]" in "".join(l.text for l in lines)


class TestHeader:
    def test_names_the_method_version_and_cost(self, default_run):
        scan, energy, ratio = default_run
        header = provenance_lines(scan, energy, ratio)[0].text
        assert scan.provenance.method in header
        assert "quadrature pts" in header
        assert "ms" in header
        assert "λ/period" in header

    def test_survives_a_solver_that_recorded_no_wall_time(self):
        """Imported reference data has no timing. Formatting None with a
        multiplication would be a crash in the one panel meant to explain
        things."""
        reference = make_scan()
        imported = EfficiencyScan(
            wavelengths=reference.wavelengths,
            orders=reference.orders,
            efficiency=reference.efficiency,
            propagating=reference.propagating,
            provenance=Provenance("integral", version="PCGrate-SX 6.7.1"),
        )
        header = lines_for(imported)[0].text
        assert "time not recorded" in header


class TestCancellation:
    def test_says_what_is_actually_still_happening(self):
        """Cancel abandons the result; it cannot stop the CPU. Saying
        'cancelled' alone would overstate it."""
        lines = lines_for(make_scan(), cancelled=True)
        text = "".join(l.text for l in lines)
        assert "previous result" in text
        assert "still finishing" in text

    def test_is_not_styled_as_an_error(self):
        lines = lines_for(make_scan(), cancelled=True)
        assert "cancelled" not in alarming_text(lines)

    def test_absent_unless_asked_for(self):
        assert "cancelled" not in "".join(l.text for l in lines_for(make_scan()))


class TestErrors:
    def test_counts_and_names_every_bad_field(self):
        lines = error_lines(
            (FieldError("period", "must be positive"), FieldError("alpha", "too big"))
        )
        text = "".join(l.text for l in lines)
        assert "2 field(s) need attention" in text
        assert "period: must be positive" in text
        assert "alpha: too big" in text

    def test_all_of_it_is_an_alarm(self):
        """Unlike the provenance case, a rejected form genuinely is a
        problem, and every line should look like one."""
        lines = error_lines((FieldError("period", "bad"),))
        assert all(line.tag in ALARM_TAGS for line in lines)


class TestSolvingLines:
    def test_states_what_is_running(self):
        text = "".join(l.text for l in solving_lines("scalar", 200))
        assert "scalar" in text and "200" in text

    def test_is_not_an_alarm(self):
        assert alarming_text(solving_lines("scalar", 200)) == ""


class TestHtml:
    def test_every_tag_has_a_colour(self):
        """A typo in a tag name would otherwise render as unstyled text with
        no failure anywhere."""
        used = {"warn", "bad", "ok", "dim"}
        assert used <= set(TAG_COLORS)

    def test_tagged_lines_carry_their_colour(self):
        html = to_html([Line("boom\n", "bad"), Line("plain\n")])
        assert TAG_COLORS["bad"] in html
        assert "plain" in html

    def test_untagged_lines_get_no_span(self):
        assert "<span" not in to_html([Line("plain\n")])

    def test_markup_in_a_warning_cannot_become_markup(self):
        """Solver warnings are prose written elsewhere; one containing a
        bracket must not silently vanish into a tag."""
        html = to_html([Line("λ < 2 nm & shallow\n", "warn")])
        assert "&lt;" in html and "&amp;" in html
        assert "<2" not in html

    def test_newlines_survive_as_line_breaks(self):
        assert to_html([Line("a\nb\n")]).count("<br>") == 2

    def test_renders_the_real_default_run(self, default_run):
        scan, energy, ratio = default_run
        html = to_html(provenance_lines(scan, energy, ratio))
        assert "summed efficiency" in html
        assert TAG_COLORS["warn"] in html


class TestAbout:
    def test_reports_the_real_version_and_licence(self):
        from gratinglab import __version__

        text = about_text(__version__)
        assert __version__ in text
        assert "BSD-3-Clause" in text
        assert "GratingLab" in text


class TestPurity:
    def test_imports_no_toolkit_and_no_plotting(self):
        import ast
        import inspect

        import gratinglab.gui.provenance as module

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(module))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6", "matplotlib"} & imported


class TestAgainstARealSolve:
    def test_an_off_plane_run_produces_a_complete_panel(self):
        """Every section present, on the geometry the app opens with."""
        scan = scalar.solve(
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)),
            Illumination(alpha_deg=25.0, gamma_deg=1.5, polarization=UNPOL),
            np.linspace(1.0, 5.0, 20),
            quadrature_points=2048,
        )
        text = "".join(
            l.text for l in provenance_lines(scan, check_energy_balance(scan), 0.0159)
        )
        for expected in ("scalar", "convergence:", "normalization:", "energy balance:"):
            assert expected in text
