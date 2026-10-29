"""The C-method's tab: it appears, solves, warns, cancels and sweeps.

Wiring only, as ``test_gui_qt.py``'s ``TestTheIntegralTab`` is: the guards
belong to ``test_gui_cmethod_options.py`` and the physics to
``test_cmethod.py``. Kept in its own file so a later solver's tab tests do
not collide with these in one shared module.
"""

import numpy as np
import pytest

pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"'
)

from gratinglab.gui.cmethod_options import (  # noqa: E402
    CONDUCTIVITY_MODES,
    CMethodOptionsState,
    build_options,
)
from gratinglab.gui.state import FormState, build  # noqa: E402
from gratinglab.solvers.cmethod import cmethod  # noqa: E402

SOLVE_TIMEOUT_MS = 30_000

#: Smooth, classical, two wavelengths: solves in well under a second and
#: plateaus inside the sweep's ladder.
SMOOTH = dict(
    period="500", profile_kind="Sinusoidal", depth_fraction="0.2",
    mount="Classical", alpha="10",
    wavelength_start="450", wavelength_stop="600", wavelength_count="2",
)

#: Short soft X-ray scan with a coating, for the tabulated wiring.
TABULATED = dict(
    coating="Au", wavelength_start="4.0", wavelength_stop="4.5",
    wavelength_count="2",
)


@pytest.fixture
def win(qtbot):
    """A constructed window with its first (scalar) solve already drawn."""
    from gratinglab.gui.qt.main_window import MainWindow

    window = MainWindow()
    qtbot.addWidget(window)
    with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
        pass
    return window


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


def resolve(qtbot, win):
    with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
        win.solve("cmethod")


@pytest.fixture
def tab(win):
    return win.tabs["cmethod"]


class TestItAppears:
    def test_the_tab_is_registered_under_its_solver_name(self, win):
        from gratinglab.gui.qt.cmethod_tab import CMethodTab

        assert isinstance(win.tabs["cmethod"], CMethodTab)

    def test_its_title_is_the_help_menus(self, win):
        """The Help menu lists the theory page as "C-method (Chandezon)"; the
        tab is the same solver, labelled with the short form of that name and
        carrying the full one as its tooltip, as every solver tab does."""
        from gratinglab.gui.docs import display_title, short_title

        bar = win._tab_widget
        index = next(i for i in range(bar.count()) if bar.widget(i) is win.tabs["cmethod"])
        assert display_title("cmethod") == "C-method (Chandezon)"
        assert bar.tabText(index) == short_title("cmethod") == "C-method"
        assert bar.tabToolTip(index) == display_title("cmethod")

    def test_the_window_still_opens_on_scalar(self, win):
        """A rigorous method is a deliberate click away, not the opening solve."""
        assert win._tab_widget.currentWidget() is win.tabs["scalar"]
        assert win.tabs["cmethod"]._scan is None


class TestTheControls:
    def test_the_conductivity_selector_offers_the_solvers_modes(self, tab):
        combo = tab._fields["conductivity"]
        assert tuple(combo.itemText(i) for i in range(combo.count())) == CONDUCTIVITY_MODES
        assert combo.currentText() == "perfect"

    def test_the_truncation_starts_at_the_solvers_default(self, tab):
        assert tab._fields["truncation"].text() == CMethodOptionsState().truncation

    def test_the_note_follows_the_selector(self, tab):
        before = tab._conductivity_note.text()
        assert "Perfect conductor" in before
        tab._fields["conductivity"].setCurrentText("tabulated")
        after = tab._conductivity_note.text()
        assert after != before
        assert "absolute" in after

    def test_both_fields_reach_build_options(self, tab):
        tab._fields["conductivity"].setCurrentText("tabulated")
        tab._fields["truncation"].setText("16")
        parsed = build(FormState(**TABULATED))
        assert tab.build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths
        ) == {"truncation": 16, "conductivity": "tabulated"}


class TestItSolves:
    def test_a_solve_equals_a_direct_one(self, qtbot, win, tab):
        """Not a second source of truth: the tab passes the scan through."""
        set_geometry(win, **SMOOTH)
        tab._fields["truncation"].setText("12")
        resolve(qtbot, win)

        parsed = build(FormState(**SMOOTH))
        expected = cmethod.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, truncation=12
        )
        assert np.array_equal(tab._scan.efficiency, expected.efficiency)
        assert "12 truncation" in tab._provenance.toPlainText()

    def test_a_tabulated_solve_equals_a_direct_one(self, qtbot, win, tab):
        set_geometry(win, **TABULATED)
        tab._fields["truncation"].setText("16")
        tab._fields["conductivity"].setCurrentText("tabulated")
        resolve(qtbot, win)

        parsed = build(FormState(**TABULATED))
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            CMethodOptionsState(truncation="16", conductivity="tabulated"),
        )
        expected = cmethod.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(tab._scan.efficiency, expected.efficiency)
        assert np.array_equal(tab._scan.absorption, expected.absorption)
        shown = tab._provenance.toPlainText()
        assert "absolute (Au)" in shown
        assert "Σ + A" in shown

    def test_the_solvers_warnings_reach_the_panel(self, qtbot, win, tab):
        """The default blazed profile has corners, and at 32 harmonics its
        70.5 deg facet is within a factor of two of the 29 it needs at 1 nm:
        both of the solver's warnings are printed as warnings, word for
        word."""
        set_geometry(win, wavelength_count="3")
        tab._fields["truncation"].setText("32")
        resolve(qtbot, win)

        warnings = tab._scan.provenance.warnings
        assert any("corners" in w for w in warnings)
        assert any("steep wall" in w for w in warnings)
        shown = tab._provenance.toPlainText()
        for warning in warnings:
            assert f"⚠ {warning}" in shown
        assert any(line.tag == "warn" for line in tab._lines)

    def test_a_smooth_profile_carries_no_corner_warning(self, qtbot, win, tab):
        """Non-vacuity for the test above."""
        set_geometry(win, **SMOOTH)
        resolve(qtbot, win)
        assert not any("corners" in w for w in tab._scan.provenance.warnings)
        assert "corners" not in tab._provenance.toPlainText()


class TestRefusedOnTheForm:
    def test_too_few_harmonics_never_reaches_the_worker(self, win, tab):
        tab._fields["truncation"].setText("2")
        win.solve("cmethod")
        assert not win._running
        shown = tab._provenance.toPlainText()
        assert "field(s) need attention" in shown
        assert "truncation" in shown

    def test_a_wall_too_steep_is_refused_by_name(self, qtbot, win, tab):
        """Issue #7: an 89 deg anti-blaze facet at 1 nm needs more harmonics
        than the solver accepts. The form cannot fix that, so the solver's
        refusal reaches the panel with its reason."""
        set_geometry(win, antiblaze_angle="89.0", wavelength_count="2")
        win.solve("cmethod")
        qtbot.waitUntil(lambda: not win._running, timeout=SOLVE_TIMEOUT_MS)
        shown = tab._provenance.toPlainText()
        assert "UnsupportedConfiguration" in shown
        assert "too steep" in shown
        assert "89 deg" in shown

    def test_tabulated_without_a_coating_names_the_coating(self, win, tab):
        tab._fields["conductivity"].setCurrentText("tabulated")
        win.solve("cmethod")
        shown = tab._provenance.toPlainText()
        assert "coating" in shown
        assert "field(s) need attention" in shown


class TestItCancels:
    def test_the_solver_reports_progress_so_cancel_is_offered(self, win, tab):
        assert cmethod.capabilities.reports_progress
        win.solve("cmethod")
        assert tab._cancel_button.isEnabled()
        assert not tab._solve_button.isEnabled()
        assert not tab._fields["truncation"].isEnabled()
        win.cancel()

    def test_cancel_stops_it_mid_scan_and_is_not_a_failure(self, qtbot, win, tab):
        """The default scan is 200 wavelengths at about 0.1 s each, so a
        Cancel after the first progress report lands well inside it."""
        reached, events = [], []
        win._worker.progress.connect(lambda t, d, n: reached.append(d))
        win._worker.cancelled.connect(lambda t, m: events.append("cancelled"))
        win._worker.failed.connect(lambda t, m, msg: events.append("failed"))

        win.solve("cmethod")
        qtbot.waitUntil(lambda: any(d > 0 for d in reached), timeout=SOLVE_TIMEOUT_MS)
        win.cancel()
        qtbot.waitUntil(lambda: bool(events), timeout=SOLVE_TIMEOUT_MS)

        assert events == ["cancelled"]
        assert max(reached) < 100, f"reached {max(reached)} of 200 despite Cancel"
        assert tab._scan is None
        assert tab._fields["truncation"].isEnabled()


class TestItSweeps:
    def _sweep(self, qtbot, win, tab):
        with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
            tab._converge_button.click()

    def test_the_sweep_plateaus_and_says_so(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        self._sweep(qtbot, win, tab)
        shown = tab._provenance.toPlainText()
        assert "convergence: yes" in shown
        assert "truncation=4 is enough" in shown

    def test_the_knob_ends_at_the_coarsest_defensible_setting(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        field = tab._fields["truncation"]
        assert field.text() == "32"
        self._sweep(qtbot, win, tab)
        assert field.text() == "4"
        assert field.isEnabled()

    def test_the_ladder_stops_at_the_solvers_ceiling(self, qtbot, win, tab):
        """Without the ceiling the default ladder runs on to 2048, and a rung
        the solver refuses after others solved fails the whole sweep. The bar
        counts the rungs, so its range is the ladder's length: 4 to 256."""
        set_geometry(win, **SMOOTH)
        self._sweep(qtbot, win, tab)
        assert tab._progress.maximum() == 7

    def test_cancelling_a_sweep_keeps_the_previous_result(self, qtbot, win, tab):
        set_geometry(win, **SMOOTH)
        resolve(qtbot, win)
        previous = tab._scan
        win.solve("cmethod", sweep=True)
        win.cancel()
        assert tab._scan is previous
