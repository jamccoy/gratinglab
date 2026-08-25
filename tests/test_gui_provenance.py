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
from gratinglab.gui.scalar_options import ScalarOptionsState, build_options
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
    options = build_options(
        parsed.problem, parsed.illumination, parsed.wavelengths, ScalarOptionsState()
    )
    scan = scalar.solve(
        parsed.problem, parsed.illumination, parsed.wavelengths, **options
    )
    return scan, check_energy_balance(scan), parsed.lambda_over_period


def make_scan(
    *,
    converged=None,
    warnings=(),
    notes=None,
    totals=(0.5, 0.6),
    method="scalar",
    truncation=2048,
    absorption=None,
):
    """A minimal scan with exactly the provenance under test."""
    values = np.asarray(totals, dtype=float)
    return EfficiencyScan(
        wavelengths=np.arange(1.0, len(values) + 1.0),
        orders=np.array([0, 1]),
        efficiency=np.column_stack([values * 0.6, values * 0.4]),
        propagating=np.ones((len(values), 2), dtype=bool),
        provenance=Provenance(
            method,
            version="test",
            truncation=truncation,
            converged=converged,
            wall_time_s=0.071,
            warnings=tuple(warnings),
            notes=dict(notes or {}),
        ),
        absorption=None if absorption is None else np.asarray(absorption, dtype=float),
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

    def test_not_checked_is_not_an_alarm(self):
        """A plain solve does not sweep -- that costs several solves and is
        the caller's call, not a hidden tax on pressing Solve. So this is the
        ordinary state, and it points at the way out rather than sitting
        there as a bare negative."""
        lines = lines_for(make_scan(converged=None))
        shown = "".join(line.text for line in lines)
        assert "not checked" in shown
        assert "gratinglab.convergence" in shown
        assert "not checked" not in alarming_text(lines)

    def test_relative_normalization_is_not_an_alarm(self):
        lines = lines_for(make_scan(notes={"normalization": "relative"}))
        assert "normalization: " in "".join(line.text for line in lines)
        assert "coating" not in alarming_text(lines)

    def test_a_resolved_but_unapplied_coating_is_named_not_denied(self):
        """"relative (Au ...)" is a real state: the material is known and its
        reflectivity has not been applied. Wording taken from `normalization`
        alone would report "no coating" for it, which is false."""
        lines = lines_for(
            make_scan(
                notes={"normalization": "relative", "coating": "Au (CXRO export)"}
            )
        )
        shown = "".join(line.text for line in lines)
        assert "relative (Au (CXRO export))" in shown
        assert "no coating" not in shown

    def test_and_a_bare_run_still_says_so(self):
        """Non-vacuity: the "no coating" wording has to survive for the case it
        was written for."""
        shown = "".join(
            line.text for line in lines_for(make_scan(notes={"normalization": "relative"}))
        )
        assert "relative (no coating)" in shown

    def test_a_verdict_carries_the_evidence_behind_it(self):
        """A bare "yes" is a claim. `converged_at` is what a reader can check
        -- and the actionable half, since it is usually cheaper than the value
        it took to prove it."""
        study = {"knob": "quadrature_points", "converged_at": 4096,
                 "values": [256, 512, 1024, 2048, 4096, 8192, 16384]}
        shown = "".join(
            line.text for line in lines_for(
                make_scan(converged=True, notes={"convergence": study})
            )
        )
        assert "quadrature_points=4096 is enough" in shown
        assert "swept to 16384" in shown

    def test_and_a_verdict_without_one_says_less_rather_than_inventing_it(self):
        """Non-vacuity for the test above, and the honest degradation: an
        imported scan can carry `converged` with no study behind it."""
        shown = "".join(line.text for line in lines_for(make_scan(converged=True)))
        assert "convergence: yes" in shown
        assert "is enough" not in shown

    def test_a_demonstrated_failure_to_converge_IS_an_alarm(self):
        """The other half of the rule -- the check must still be able to fail,
        or dimming everything would have been the same fix."""
        lines = lines_for(make_scan(converged=False))
        assert "NO" in alarming_text(lines)

    def test_convergence_success_reads_as_success(self):
        lines = lines_for(make_scan(converged=True))
        assert any(line.tag == "ok" and line.text.strip() == "yes" for line in lines)


class TestWarnings:
    def test_every_solver_warning_is_shown_verbatim(self):
        """Never summarised, never truncated: the warnings are why the panel
        exists."""
        message = "λ/period is 0.6, past the scalar validity bound"
        lines = lines_for(make_scan(warnings=[message]))
        assert message in alarming_text(lines)

    def test_warnings_keep_their_order(self):
        lines = lines_for(make_scan(warnings=["first", "second", "third"]))
        text = "".join(line.text for line in lines)
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
        energy = next(line for line in lines if "Σ ∈" in line.text)
        assert energy.tag == "ok"
        assert "EXCEEDS UNITY" not in energy.text

    def test_an_excess_reads_as_a_failure(self):
        lines = lines_for(make_scan(totals=(0.9, 1.5)))
        energy = next(line for line in lines if "Σ ∈" in line.text)
        assert energy.tag == "bad"
        assert "EXCEEDS UNITY" in energy.text

    def test_the_range_is_reported_not_a_single_number(self):
        """Σ varies across the scan; one number would hide that."""
        lines = lines_for(make_scan(totals=(0.4, 0.9)))
        assert "Σ ∈ [0.4000, 0.9000]" in "".join(line.text for line in lines)


class TestHeader:
    def test_names_the_method_version_and_cost(self, default_run):
        scan, energy, ratio = default_run
        header = provenance_lines(scan, energy, ratio)[0].text
        assert scan.provenance.method in header
        assert "quadrature points" in header
        assert "ms" in header
        assert "λ/period" in header

    def test_the_knob_is_named_per_method(self):
        """The header said "quadrature pts" whatever produced the number,
        which was wrong for the integral solver's boundary points from the day
        that tab existed. Read off the solver's declared accuracy knob so a
        backend added later names its own."""
        header = lines_for(make_scan(method="integral", truncation=400))[0].text
        assert "400 boundary points" in header
        assert "quadrature" not in header

    def test_an_unregistered_method_still_renders(self):
        """Imported reference data carries a foreign method string. The panel
        that exists to explain things must not be the thing that crashes."""
        header = lines_for(make_scan(method="pcgrate-sx", truncation=91))[0].text
        assert "91 truncation" in header

    def test_a_scan_with_no_truncation_says_nothing_rather_than_None(self):
        header = lines_for(make_scan(truncation=None))[0].text
        assert "None" not in header
        assert "ms" in header

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
        text = "".join(line.text for line in lines)
        assert "previous result" in text
        assert "the calculation stopped" in text

    def test_is_not_styled_as_an_error(self):
        lines = lines_for(make_scan(), cancelled=True)
        assert "cancelled" not in alarming_text(lines)

    def test_absent_unless_asked_for(self):
        assert "cancelled" not in "".join(line.text for line in lines_for(make_scan()))


class TestErrors:
    def test_counts_and_names_every_bad_field(self):
        lines = error_lines(
            (FieldError("period", "must be positive"), FieldError("alpha", "too big"))
        )
        text = "".join(line.text for line in lines)
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
        text = "".join(line.text for line in solving_lines("scalar", 200))
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


class TestAbsorption:
    """What a finite-conductivity scan adds, and what it must not look like.

    `check_energy_balance` folds absorption into the total when a scan carries
    one, so the number beside "energy balance" is R + A there while the Σ curve
    on the plot is R alone. Two quantities under one symbol is the same class
    of mistake as tagging a correct default as a problem.
    """

    def test_the_absorbed_fraction_is_reported(self):
        shown = "".join(
            line.text
            for line in lines_for(
                make_scan(
                    method="integral",
                    totals=(0.7, 0.75),
                    absorption=(0.3, 0.25),
                    notes={"normalization": "absolute", "material": "Au"},
                )
            )
        )
        assert "absorption: A ∈ [0.2500, 0.3000]" in shown

    def test_absorption_is_not_an_alarm(self):
        """Power going into the metal is the answer, not a warning about it."""
        lines = lines_for(
            make_scan(method="integral", totals=(0.7,), absorption=(0.3,))
        )
        assert "absorption" not in alarming_text(lines)

    def test_the_conserved_quantity_is_named_when_it_includes_absorption(self):
        shown = "".join(
            line.text
            for line in lines_for(
                make_scan(method="integral", totals=(0.7, 0.75), absorption=(0.3, 0.25))
            )
        )
        assert "Σ + A ∈ [1.0000, 1.0000]" in shown

    def test_a_scan_without_absorption_still_says_plain_sigma(self):
        """Non-vacuity: every scalar run, and every perfect-conductivity one."""
        shown = "".join(line.text for line in lines_for(make_scan(totals=(0.4, 0.9))))
        assert "Σ ∈ [0.4000, 0.9000]" in shown
        assert "absorption" not in shown

    def test_the_material_a_boundary_condition_read_is_named(self):
        """The integral solver records `material`, not `coating`. Reading only
        the latter reported "absolute (no coating)" for a run whose whole
        result came from the material."""
        shown = "".join(
            line.text
            for line in lines_for(
                make_scan(
                    method="integral",
                    totals=(0.7,),
                    absorption=(0.3,),
                    notes={"normalization": "absolute", "material": "Au"},
                )
            )
        )
        assert "absolute (Au)" in shown


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
            line.text for line in provenance_lines(scan, check_energy_balance(scan), 0.0159)
        )
        for expected in ("scalar", "convergence:", "normalization:", "energy balance:"):
            assert expected in text
