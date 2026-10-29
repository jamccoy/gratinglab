"""The Compare tab: several solvers on one geometry, refusals named, cancellable.

Two halves, as the solver tabs' tests are split: the drawing and the panel's
wording are decided in ``gui/comparison.py`` without Qt and are tested
directly; the tab's wiring (which solvers run, with whose settings, what a
refusal looks like, Cancel) is tested through a real window. The physics
belongs to each solver's own tests and the agreement between them to
``test_cross_method.py``; what is checked here is that the tab shows the
scans the solvers produced, unaltered.
"""

import numpy as np
import pytest

from gratinglab.compare import align, sweep
from gratinglab.gui.comparison import (
    VERIFICATION_NOTE,
    ComparisonOutcome,
    comparison_lines,
    draw_comparison,
    line_style,
    pick_reference,
)
from gratinglab.gui.provenance import alarming_text
from gratinglab.gui.state import FormState, build
from gratinglab.gui.theme import DARK, LIGHT, order_colors

SOLVE_TIMEOUT_MS = 60_000

#: Smooth, classical, two wavelengths: the C-method and scalar solve it in
#: well under a second, so a comparison test is quick.
SMOOTH = dict(
    period="500", profile_kind="Sinusoidal", depth_fraction="0.2",
    mount="Classical", alpha="10",
    wavelength_start="450", wavelength_stop="600", wavelength_count="2",
)

#: The same, but lamellar: vertical sidewalls, which the C-method refuses.
LAMELLAR = dict(SMOOTH, profile_kind="Lamellar", duty_cycle="0.5")

SOLVERS = ("cmethod", "integral", "rcwa", "scalar")


def _outcome(methods=("scalar", "cmethod"), *, refused=None, failed=None, form=SMOOTH):
    parsed = build(FormState(**form))
    scans = sweep(parsed.problem, parsed.illumination, parsed.wavelengths, list(methods))
    requested = tuple(methods) + tuple(refused or ()) + tuple(failed or ())
    return ComparisonOutcome(
        requested=requested,
        comparison=align(scans),
        refused=dict(refused or {}),
        failed=dict(failed or {}),
    )


# -- without Qt --------------------------------------------------------------


class TestLineStyles:
    def test_each_solver_keeps_its_style_whatever_else_is_ticked(self):
        styles = [line_style(n, SOLVERS) for n in SOLVERS]
        assert len(set(styles)) == len(SOLVERS)
        assert line_style("scalar", SOLVERS) == styles[SOLVERS.index("scalar")]

    def test_an_unknown_solver_still_gets_a_style(self):
        assert line_style("pcgrate", SOLVERS) in ("-", "--", "-.", ":")


class TestPickReference:
    def test_the_chosen_one_when_it_solved(self):
        assert pick_reference("scalar", ("cmethod", "scalar")) == "scalar"

    def test_a_rigorous_one_when_the_chosen_one_did_not_solve(self):
        """The scalar solver is approximate, so it is never the fallback
        while a rigorous solver is available."""
        assert pick_reference("integral", ("scalar", "cmethod")) == "cmethod"

    def test_the_first_when_none_is_rigorous(self):
        assert pick_reference("integral", ("scalar",)) == "scalar"

    def test_none_when_nothing_solved(self):
        assert pick_reference("integral", ()) is None


class TestThePanel:
    def test_a_refusal_is_named_with_its_reason_and_not_as_a_fault(self):
        """A solver declining a problem it cannot model is correct behaviour,
        so it is a warning, never the red of a failure."""
        outcome = _outcome(refused={"rcwa": "needs a coating"})
        lines = comparison_lines(outcome, "cmethod")
        text = "".join(line.text for line in lines)
        assert "RCWA: refused: needs a coating" in text
        refusal = next(line for line in lines if "needs a coating" in line.text)
        assert refusal.tag == "warn"

    def test_a_failure_is_named_and_is_a_fault(self):
        outcome = _outcome(failed={"integral": "ValueError: contrived"})
        lines = comparison_lines(outcome, "cmethod")
        failure = next(line for line in lines if "contrived" in line.text)
        assert failure.tag == "bad"
        assert "Integral: " in "".join(line.text for line in lines)

    def test_every_requested_solver_is_accounted_for(self):
        outcome = _outcome(refused={"rcwa": "r"}, failed={"integral": "f"})
        text = "".join(line.text for line in comparison_lines(outcome, "cmethod"))
        for title in ("C-method:", "Scalar:", "RCWA:", "Integral:"):
            assert title in text

    def test_it_states_the_disagreement_against_the_reference(self):
        outcome = _outcome()
        summary = outcome.comparison.summary("scalar", "cmethod")
        text = "".join(line.text for line in comparison_lines(outcome, "cmethod"))
        assert "reference: C-method" in text
        assert f"max |Δ| {summary['max_abs_difference']:.3g}" in text

    def test_a_substituted_reference_is_said_out_loud(self):
        outcome = _outcome()
        text = "".join(line.text for line in comparison_lines(outcome, "integral"))
        assert "Integral did not solve" in text
        assert "differences are from C-method" in text

    def test_a_clean_run_has_nothing_alarming(self):
        """Non-vacuity for the two tests above: an ordinary comparison where
        every solver solved must not read as a problem."""
        lines = comparison_lines(_outcome(), "cmethod")
        assert alarming_text(lines) == ""

    def test_agreement_is_never_presented_as_verification(self):
        """However closely the solvers agree, the panel says that agreement
        is not verification, and says it quietly: it is a fact about every
        method here, not a fault in this run."""
        lines = comparison_lines(_outcome(), "cmethod")
        note = next(line for line in lines if line.text == VERIFICATION_NOTE)
        assert note.tag == "dim"
        assert "not verification" in note.text

    def test_nothing_solved_says_so(self):
        outcome = ComparisonOutcome(
            requested=("cmethod",), comparison=None, refused={"cmethod": "no"}
        )
        text = "".join(line.text for line in comparison_lines(outcome, "cmethod"))
        assert "no solver produced a result" in text


class TestTheDrawing:
    @pytest.fixture
    def axes(self):
        from matplotlib.figure import Figure

        figure = Figure()
        return figure.add_subplot(2, 1, 1), figure.add_subplot(2, 1, 2)

    def test_orders_take_the_themes_order_colours(self, axes):
        outcome = _outcome()
        comparison = outcome.comparison
        colors = draw_comparison(
            *axes, comparison, "cmethod", set(int(m) for m in comparison.orders),
            LIGHT, solvers=SOLVERS,
        )
        assert colors == order_colors((int(m) for m in comparison.orders), LIGHT)

    def test_the_curves_are_the_scans_unaltered(self, axes):
        """Upper panel: one curve per solver and order. Lower: solver minus
        reference, for every solver but the reference."""
        comparison = _outcome().comparison
        visible = {int(m) for m in comparison.orders}
        draw_comparison(*axes, comparison, "cmethod", visible, LIGHT, solvers=SOLVERS)
        upper, lower = axes
        assert len(upper.get_lines()) == 2 * len(visible) + 1  # + the unity line
        expected = comparison.difference("scalar", "cmethod")
        # Every curve is drawn against the wavelengths; the zero line is not.
        drawn = [line.get_ydata() for line in lower.get_lines()
                 if np.array_equal(line.get_xdata(), comparison.wavelengths)]
        assert len(drawn) == len(visible)
        for column in range(len(comparison.orders)):
            assert any(np.array_equal(y, expected[:, column]) for y in drawn)

    def test_log_difference_drops_exact_zeros_rather_than_inventing_a_floor(self, axes):
        comparison = _outcome().comparison
        visible = {int(m) for m in comparison.orders}
        draw_comparison(
            *axes, comparison, "cmethod", visible, LIGHT,
            solvers=SOLVERS, log_difference=True,
        )
        lower = axes[1]
        assert lower.get_yscale() == "log"
        for line in lower.get_lines():
            values = np.asarray(line.get_ydata(), dtype=float)
            assert not (values <= 0).any()

    def test_hidden_orders_are_not_drawn(self, axes):
        comparison = _outcome().comparison
        draw_comparison(*axes, comparison, "cmethod", set(), LIGHT, solvers=SOLVERS)
        assert len(axes[0].get_lines()) == 1  # the unity line alone


# -- through the window ------------------------------------------------------


pytest.importorskip("PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"')


@pytest.fixture
def win(qtbot):
    from gratinglab.gui.qt.main_window import MainWindow

    window = MainWindow()
    qtbot.addWidget(window)
    with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
        pass
    return window


@pytest.fixture
def tab(win):
    return win.compare_tab


def set_geometry(win, **fields):
    for field, value in fields.items():
        widget = win.geometry._fields[field]
        if hasattr(widget, "setCurrentIndex"):
            index = widget.findData(value)
            if index < 0:
                index = widget.findText(value)
            assert index >= 0, f"{field} offers no {value!r}"
            widget.setCurrentIndex(index)
        else:
            widget.setText(value)


def compare(qtbot, win, *names):
    win.compare_tab.set_selected(names)
    with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
        win.compare()


class TestItAppears:
    def test_it_is_the_last_tab_and_not_a_solver(self, win, tab):
        bar = win._tab_widget
        assert bar.widget(bar.count() - 1) is tab
        assert bar.tabText(bar.count() - 1) == "Compare"
        assert tab not in win.tabs.values()

    def test_it_offers_every_solver_tab(self, win, tab):
        from gratinglab.gui.docs import short_title

        assert tuple(tab._checks) == tuple(win.tabs)
        for name, box in tab._checks.items():
            assert box.text() == short_title(name)

    def test_it_does_not_solve_on_opening(self, tab):
        """Every rigorous solver on the default 200-wavelength scan takes
        minutes; that is a deliberate click, not a side effect of opening."""
        assert tab._outcome is None

    def test_nothing_ticked_means_nothing_to_press(self, tab):
        tab.set_selected(())
        assert not tab._compare_button.isEnabled()
        tab.set_selected(("scalar",))
        assert tab._compare_button.isEnabled()


class TestItCompares:
    def test_the_plotted_scans_are_a_direct_sweep(self, qtbot, win, tab):
        """Not a second source of truth: what the tab holds is what
        `compare.sweep` and `compare.align` return for the same problem with
        each tab's own options."""
        set_geometry(win, **SMOOTH)
        win.tabs["cmethod"]._fields["truncation"].setText("12")
        compare(qtbot, win, "scalar", "cmethod")

        parsed = build(FormState(**SMOOTH))
        expected = align(sweep(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ["cmethod", "scalar"], options={"cmethod": {"truncation": 12}},
        ))
        got = tab._outcome.comparison
        assert got.methods == ("cmethod", "scalar")
        assert np.array_equal(got.efficiency, expected.efficiency)

    def test_each_solver_uses_its_own_tabs_settings(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        win.tabs["cmethod"]._fields["truncation"].setText("9")
        compare(qtbot, win, "cmethod")
        scan = tab._outcome.comparison.scans[0]
        assert scan.provenance.truncation == 9

    def test_the_reference_changes_without_solving_again(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        compare(qtbot, win, "scalar", "cmethod")
        solves = []
        win._worker.compared.connect(lambda *a: solves.append(a))
        tab.set_reference("scalar")
        assert tab._difference_axes.get_title(loc="left") == "Difference from Scalar"
        tab.set_reference("cmethod")
        assert tab._difference_axes.get_title(loc="left") == "Difference from C-method"
        qtbot.wait(50)
        assert solves == []

    def test_an_unticked_reference_falls_back_and_says_so(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        tab.set_reference("integral")
        compare(qtbot, win, "scalar", "cmethod")
        assert tab._difference_axes.get_title(loc="left") == "Difference from C-method"
        assert "Integral did not solve" in tab._panel.toPlainText()

    def test_the_order_list_shows_every_order_with_its_colour(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        compare(qtbot, win, "scalar", "cmethod")
        comparison = tab._outcome.comparison
        assert tab._order_list.count() == len(comparison.orders)
        assert tab._order_colors == order_colors(
            (int(m) for m in comparison.orders), LIGHT
        )

    def test_unchecking_an_order_removes_its_curves(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        compare(qtbot, win, "scalar", "cmethod")
        from PySide6.QtCore import Qt

        tab._show_all_orders()
        before = len(tab._axes.get_lines())
        tab._order_list.item(0).setCheckState(Qt.CheckState.Unchecked)
        assert len(tab._axes.get_lines()) == before - 2  # one per solver

    def test_enter_in_the_geometry_compares_when_this_tab_is_showing(
        self, qtbot, win, tab
    ):
        set_geometry(win, **SMOOTH)
        tab.set_selected(("scalar",))
        win._tab_widget.setCurrentWidget(tab)
        with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
            win._solve_active_tab()
        assert tab._outcome.comparison.methods == ("scalar",)


class TestRefusals:
    def test_a_solver_that_cannot_model_the_profile_is_named_and_the_rest_run(
        self, qtbot, win, tab
    ):
        """The C-method has no finite slope on a lamellar sidewall and says
        so; that is reported by name, and the scalar answer still arrives."""
        set_geometry(win, **LAMELLAR)
        compare(qtbot, win, "scalar", "cmethod")
        outcome = tab._outcome
        assert outcome.solved == ("scalar",)
        assert "cmethod" in outcome.refused
        assert "sidewalls" in outcome.refused["cmethod"]
        shown = tab._panel.toPlainText()
        assert "C-method: refused:" in shown
        assert "sidewalls" in shown
        assert "Scalar: solved" in shown

    def test_a_form_its_own_tab_rejects_is_a_refusal_not_a_stop(self, qtbot, win, tab):
        """Too few harmonics on the C-method tab: the field shows the error
        there, the comparison names it, and the scalar solver still runs."""
        set_geometry(win, **SMOOTH)
        win.tabs["cmethod"]._fields["truncation"].setText("0")
        compare(qtbot, win, "scalar", "cmethod")
        assert tab._outcome.solved == ("scalar",)
        assert "on the C-method tab" in tab._outcome.refused["cmethod"]
        assert "need attention" in win.tabs["cmethod"]._provenance.toPlainText()

    def test_layers_are_refused_by_name_for_the_solvers_that_cannot(
        self, qtbot, win, tab
    ):
        set_geometry(win, **SMOOTH, coating="Au", layers="Au 5")
        set_geometry(win, wavelength_start="4.0", wavelength_stop="4.5")
        compare(qtbot, win, "scalar", "cmethod")
        assert tab._outcome.comparison is None
        assert set(tab._outcome.refused) == {"scalar", "cmethod"}
        assert "no solver produced a result" in tab._panel.toPlainText()

    def test_a_geometry_error_stops_the_whole_comparison(self, win, tab):
        set_geometry(win, period="-1")
        tab.set_selected(("scalar",))
        win.compare()
        assert not win._running
        assert "need attention" in tab._panel.toPlainText()


class TestRunning:
    def test_every_solve_button_is_disabled_while_it_runs(self, win, tab):
        tab.set_selected(("scalar",))
        win.compare()
        assert all(not t._solve_button.isEnabled() for t in win.tabs.values())
        assert not tab._compare_button.isEnabled()
        assert tab._cancel_button.isEnabled()
        win.cancel()
        assert all(t._solve_button.isEnabled() for t in win.tabs.values())
        assert tab._compare_button.isEnabled()

    def test_a_single_solve_disables_compare_too(self, win, tab):
        win.solve("scalar")
        assert not tab._compare_button.isEnabled()
        win.cancel()
        assert tab._compare_button.isEnabled()

    def test_cancel_stops_it_mid_run_and_is_not_a_failure(self, qtbot, win, tab):
        """The default scan is 200 wavelengths, about 0.1 s each for the
        C-method, so a Cancel after the first progress report lands inside
        the first solver's block and the second never starts."""
        reached, events = [], []
        win._worker.progress.connect(lambda t, d, n: reached.append((d, n)))
        win._worker.cancelled.connect(lambda t, m: events.append(("cancelled", m)))
        win._worker.failed.connect(lambda t, m, msg: events.append(("failed", m)))
        win._worker.compared.connect(lambda t, o: events.append(("compared", None)))

        tab.set_selected(("cmethod", "scalar"))
        win.compare()
        qtbot.waitUntil(lambda: any(d > 0 for d, _ in reached), timeout=SOLVE_TIMEOUT_MS)
        win.cancel()
        qtbot.waitUntil(lambda: bool(events), timeout=SOLVE_TIMEOUT_MS)

        assert events == [("cancelled", "compare")]
        done, total = max(reached)
        assert total == 400, "one bar over both solvers' wavelengths"
        assert done < 200, f"reached {done} of {total} despite Cancel"
        assert tab._outcome is None
        assert "stopped" in tab._panel.toPlainText()
        assert tab._compare_button.isEnabled()

    def test_a_cancelled_comparison_keeps_the_previous_one(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        compare(qtbot, win, "scalar")
        previous = tab._outcome
        win.compare()
        win.cancel()
        assert tab._outcome is previous
        assert "showing the previous one" in tab._panel.toPlainText()

    def test_a_stale_outcome_is_dropped_when_it_arrives(self, win, tab):
        tab.set_selected(("scalar",))
        win.compare()
        stale = win._tokens.current
        win.cancel()
        win._on_compared(stale, object())  # would crash if accepted
        assert tab._outcome is None

    def test_a_failure_of_the_whole_run_is_a_message_not_a_crash(self, win, tab):
        win._tokens.invalidate()
        win._active_name = tab.name
        win._set_running(True)
        win._on_failed(win._tokens.current, "compare", "ValueError: contrived")
        assert "comparison failed" in tab._panel.toPlainText()
        assert tab._compare_button.isEnabled()

    def test_progress_reaches_this_tabs_bar(self, win, tab):
        tab.set_selected(("scalar",))
        win.compare()
        win._on_progress(win._tokens.current, 3, 10)
        assert tab._progress.maximum() == 10
        assert tab._progress.value() == 3
        win.cancel()


class TestTheme:
    def test_a_theme_change_restyles_the_drawn_comparison(self, qtbot, win, tab):
        from matplotlib.colors import to_hex

        set_geometry(win, **SMOOTH)
        compare(qtbot, win, "scalar", "cmethod")
        win.set_theme("dark")
        try:
            assert to_hex(tab._axes.get_facecolor()) == DARK.axes
            assert to_hex(tab._difference_axes.get_facecolor()) == DARK.axes
            assert tab._order_colors == order_colors(
                (int(m) for m in tab._outcome.comparison.orders), DARK
            )
        finally:
            win.set_theme("light")
