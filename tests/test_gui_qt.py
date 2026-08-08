"""The Qt widget layer.

Thin on purpose, like its Tk predecessor was: the logic worth testing lives in
``gui/state.py``, ``gui/provenance.py``, ``gui/richtext.py`` and
``gui/orders.py``, all of which are pure and tested headlessly. What is checked
here is the wiring -- and above all that the window is not a second source of
truth: whatever it plots must equal what the core produces for the same inputs.

Unlike the Tk suite this replaces, these tests need no display. ``conftest.py``
sets ``QT_QPA_PLATFORM=offscreen``, under which a window always constructs, so
they run on every CI job rather than only on macOS.
"""

import numpy as np
import pytest

pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"'
)

from gratinglab.gui.scalar_options import ScalarOptionsState, build_options  # noqa: E402
from gratinglab.gui.state import FormState, build  # noqa: E402
from gratinglab.solvers import scalar  # noqa: E402

#: Generous, because it covers a whole solve plus a redraw on a loaded CI
#: runner -- but finite, because a hung worker should fail rather than hang the
#: suite. The default scalar solve takes about 70 ms.
SOLVE_TIMEOUT_MS = 10_000


@pytest.fixture
def win(qtbot):
    """A constructed window with its first solve already drawn."""
    from gratinglab.gui.qt.main_window import MainWindow

    window = MainWindow()
    qtbot.addWidget(window)  # closes and deletes it at teardown
    with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
        pass  # the solve the constructor started
    return window


def resolve(qtbot, window, name="scalar"):
    """Press Solve on the named tab and wait for the result to land."""
    with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
        window.solve(name)


class TestToolchain:
    """That embedding matplotlib in Qt works at all.

    Deliberately the first thing built. Everything downstream assumes a
    ``FigureCanvasQTAgg`` can be parented into a Qt window on a headless
    machine, and discovering otherwise five milestones later would be
    expensive. The toolbar is here because it is the one piece of matplotlib's
    Qt backend that reaches for icons and a style, which is where a missing
    system library shows up first.
    """

    def test_a_canvas_and_toolbar_embed_in_a_window(self, qtbot):
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure
        from PySide6.QtWidgets import QMainWindow, QVBoxLayout, QWidget

        window = QMainWindow()
        qtbot.addWidget(window)

        figure = Figure(figsize=(4, 3), layout="constrained")
        figure.add_subplot().plot([0, 1], [0, 1])
        canvas = FigureCanvasQTAgg(figure)
        toolbar = NavigationToolbar2QT(canvas, window)

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        window.setCentralWidget(central)

        canvas.draw()
        assert canvas.width() > 0
        assert toolbar.actions()

    def test_embedding_does_not_hijack_the_global_backend(self, qtbot):
        """Importing the Qt backend must not switch matplotlib out from under
        the rest of the suite.

        The Tk implementation called ``matplotlib.use("TkAgg")`` from inside a
        widget method, which mutated global state for every test that ran
        afterwards. Explicit canvas embedding needs no such call, and this
        pins that it stays true.
        """
        import matplotlib

        before = matplotlib.get_backend()
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure

        canvas = FigureCanvasQTAgg(Figure())
        canvas.draw()
        assert matplotlib.get_backend() == before


class TestNotASecondSourceOfTruth:
    """The single most important property of this layer.

    It survived both the toolkit change and the move to tabs unchanged, which
    is the point of having written it this way: if a tab ever computed
    anything itself, these would drift.
    """

    def test_plotted_values_equal_a_direct_solve(self, win):
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(),
        )
        expected = scalar.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        scalar_tab = win.tabs["scalar"]
        assert np.array_equal(scalar_tab._scan.efficiency, expected.efficiency)
        assert np.array_equal(scalar_tab._scan.orders, expected.orders)

    def test_changing_a_field_changes_the_result_consistently(self, qtbot, win):
        win.geometry._fields["blaze_angle"].setText("12.0")
        resolve(qtbot, win)

        parsed = build(FormState(blaze_angle="12.0"))
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(),
        )
        expected = scalar.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(win.tabs["scalar"]._scan.efficiency, expected.efficiency)


class TestBehaviour:
    def test_solves_on_construction(self, win):
        """Opening the app shows something, not an empty window."""
        scan = win.tabs["scalar"]._scan
        assert scan is not None
        assert len(scan) == 200

    def test_invalid_input_reports_errors_without_crashing(self, win):
        win.geometry._fields["period"].setText("not a number")
        win.solve("scalar")  # rejected synchronously; never reaches the worker
        text = win.tabs["scalar"]._provenance.toPlainText()
        assert "period" in text
        assert "need attention" in text

    def test_a_failed_solve_leaves_the_previous_result_intact(self, win):
        good = win.tabs["scalar"]._scan
        win.geometry._fields["period"].setText("")
        win.solve("scalar")
        assert win.tabs["scalar"]._scan is good

    def test_the_panel_shows_what_provenance_decided(self, win):
        """All this layer owes: paint the lines it was given. What each line
        says is checked in tests/test_gui_provenance.py."""
        from gratinglab.gui.provenance import provenance_lines

        scalar_tab = win.tabs["scalar"]
        expected = provenance_lines(
            scalar_tab._scan, scalar_tab._energy, scalar_tab._lambda_over_period
        )
        shown = scalar_tab._provenance.toPlainText()
        for line in expected:
            assert line.text.strip() in shown

    def test_reads_every_form_field(self, win):
        form = win.geometry.read_form()
        assert isinstance(form, FormState)
        assert form.period == "315.15"

    def test_a_geometry_field_missing_from_formstate_stops_the_window_opening(
        self, qtbot
    ):
        """`GeometryPanel._fields`' keys must exactly match FormState's field
        names -- `read_form` builds one by keyword. Previously a mismatch was
        a TypeError at solve time, visible only to whoever pressed Solve."""
        from gratinglab.gui.qt import geometry_panel as module

        original = module.GeometryPanel._build

        def drop_a_field(self):
            original(self)
            del self._fields["period"]

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(module.GeometryPanel, "_build", drop_a_field)
            with pytest.raises(AssertionError, match="period"):
                module.GeometryPanel()

    def test_a_scalar_option_field_missing_stops_the_tab_opening(self, qtbot):
        """The other half of the same guard, now that scalar's own options
        are checked inside ScalarTab rather than on MainWindow."""
        from gratinglab.gui.qt import scalar_tab as module

        original = module.ScalarTab._build_options_group

        def drop_a_field(self):
            widget = original(self)
            del self._fields["quadrature_points"]
            return widget

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(module.ScalarTab, "_build_options_group", drop_a_field)
            with pytest.raises(AssertionError, match="quadrature_points"):
                module.ScalarTab()

    def test_mount_change_relabels_the_angle_fields(self, win):
        win.geometry._fields["mount"].setCurrentText("Classical")
        assert "α" in win.geometry._angle_labels["alpha"].text()

    def test_classical_hides_the_second_angle(self, win):
        """Checked with isHidden, not isVisible: under the offscreen platform
        an unshown window reports every child invisible, which would make this
        pass without testing anything."""
        win.geometry._fields["mount"].setCurrentText("Classical")
        assert win.geometry._fields["gamma"].isHidden()

    def test_off_plane_shows_the_second_angle(self, win):
        win.geometry._fields["mount"].setCurrentText("Off-plane")
        assert not win.geometry._fields["gamma"].isHidden()

    def test_profile_change_hides_irrelevant_fields(self, win):
        """duty_cycle is meaningless for a sinusoid. The mapping itself is
        tested headlessly against PROFILE_FIELDS."""
        win.geometry._fields["profile_kind"].setCurrentText("Sinusoidal")
        assert win.geometry._profile_rows["duty_cycle"].isHidden()
        assert not win.geometry._profile_rows["depth_fraction"].isHidden()

    def test_from_file_reveals_the_load_button(self, win):
        win.geometry._fields["profile_kind"].setCurrentText("From file")
        assert not win.geometry._load_button.isHidden()
        assert win.geometry._profile_rows["blaze_angle"].isHidden()


class TestBackgroundSolve:
    def test_the_solve_does_not_run_on_the_ui_thread(self, qtbot, win):
        """`moveToThread` only redirects signal delivery, so calling the
        worker's method directly would quietly solve on the UI thread and
        freeze the window -- the exact failure this design exists to avoid."""
        import threading

        seen = {}

        from gratinglab.gui.qt import worker as worker_module

        original_run = worker_module.SolveWorker.run

        def record(self, token, method, geometry, options):
            seen["thread"] = threading.current_thread().ident
            return original_run(self, token, method, geometry, options)

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(worker_module.SolveWorker, "run", record)
            resolve(qtbot, win)

        assert "thread" in seen, "the worker slot never ran"
        assert seen["thread"] != threading.current_thread().ident

    def test_solve_is_disabled_while_one_is_running(self, win):
        win.solve("scalar")
        scalar_tab = win.tabs["scalar"]
        assert not scalar_tab._solve_button.isEnabled()
        assert scalar_tab._cancel_button.isEnabled()

    def test_a_second_solve_while_running_is_ignored(self, win):
        """One window, one result. Superseding would stack CPU-burning
        threads; queueing would surprise the user with a stale-then-fresh
        double redraw."""
        win.solve("scalar")
        token = win._token
        win.solve("scalar")
        assert win._token == token

    def test_cancel_keeps_the_previous_result_on_screen(self, win):
        scalar_tab = win.tabs["scalar"]
        previous = scalar_tab._scan
        win.solve("scalar")
        win.cancel()
        assert scalar_tab._scan is previous
        assert scalar_tab._solve_button.isEnabled()

    def test_cancel_says_the_calculation_is_still_finishing(self, win):
        """Cancel abandons a result; it cannot stop the CPU. Claiming
        otherwise would be the panel's first lie."""
        win.solve("scalar")
        win.cancel()
        assert "still finishing" in win.tabs["scalar"]._provenance.toPlainText()

    def test_a_cancelled_result_is_dropped_when_it_arrives(self, win):
        from gratinglab.gui.qt.worker import SolveResult

        win.solve("scalar")
        stale = win._token
        win.cancel()
        # would crash if accepted
        win._on_solved(stale, "scalar", SolveResult(None, None))
        assert win.tabs["scalar"]._scan is not None

    def test_a_solver_exception_becomes_a_message_not_a_crash(self, win):
        """A Python exception escaping a slot under PySide6 can abort the
        process rather than surfacing."""
        win._token += 1
        win._active_name = "scalar"
        win._set_running(True)
        win._on_failed(win._token, "scalar", "ValueError: contrived")
        scalar_tab = win.tabs["scalar"]
        assert "solve failed" in scalar_tab._provenance.toPlainText()
        assert scalar_tab._solve_button.isEnabled()

    def test_no_progress_bar_for_a_fast_solve(self, win):
        """The scalar solve takes about 70 ms; flashing a bar on and off
        within a frame reads as a glitch."""
        assert win.tabs["scalar"]._progress.isHidden()

    def test_every_tabs_solve_button_disables_together(self, win):
        """Cross-tab coordination: only one solve can be in flight
        window-wide, so no tab's Solve button may stay clickable while
        another tab's solve runs.

        With only Setup (no Solve button, added at M11-D) and Scalar
        registered, there is no *second* Solve button to prove this against
        yet -- this asserts the mechanism (`_set_running` loops every tab),
        not a second tab's behaviour. Worth revisiting once a second solver
        tab exists.
        """
        win.solve("scalar")
        assert all(not t._solve_button.isEnabled() for t in win.tabs.values())
        win.cancel()
        assert all(t._solve_button.isEnabled() for t in win.tabs.values())

    def test_closing_joins_the_worker_thread(self, qtbot):
        """Otherwise Qt reports 'QThread: Destroyed while thread is still
        running' and aborts -- most visibly under offscreen CI, where windows
        are built and torn down back to back.

        Constructed without qtbot.addWidget so that closing it here is the
        only close, and the assertion is about closeEvent rather than about
        teardown order.
        """
        from gratinglab.gui.qt.main_window import MainWindow

        window = MainWindow()
        assert window._thread.isRunning()
        window.close()
        assert not window._thread.isRunning()

    def test_the_progress_timer_does_not_outlive_the_tab_that_armed_it(self, qtbot):
        """Solve arms a 150 ms timer holding the tab's progress bar. If the
        window goes away inside that delay the timer used to fire on a deleted
        QProgressBar, and an exception escaping a slot under PySide6 can abort
        the process rather than surfacing.

        Caught by macOS CI, where the stale timer landed in the *next* test's
        setup; the race is timing-dependent and does not reproduce reliably
        here, so this exercises the path rather than pinning the crash. The
        fix is the context-object overload of `QTimer.singleShot`, which Qt
        cancels when the context dies -- verified directly: a timer whose
        context is destroyed does not fire.
        """
        from gratinglab.gui.qt.main_window import _PROGRESS_DELAY_MS, MainWindow

        window = MainWindow()
        qtbot.addWidget(window)
        window.solve("scalar")
        window.close()
        window.deleteLater()
        # pytest-qt fails the test if any slot raises during this wait.
        qtbot.wait(_PROGRESS_DELAY_MS * 3)


class TestOrderPanel:
    """Wiring for `gui/orders.py` inside `ScalarTab`.

    The rule itself -- which orders default to visible, how a re-solve
    carries a selection forward -- is already covered headlessly in
    `tests/test_gui_orders.py`. This is the wiring: does toggling a checkbox
    actually redraw, and does it do so without re-solving.
    """

    def test_every_order_from_the_default_scan_is_listed(self, win):
        scalar_tab = win.tabs["scalar"]
        assert scalar_tab._order_list.count() == len(scalar_tab._scan.orders)

    def test_the_count_label_matches_the_checked_items(self, win):
        scalar_tab = win.tabs["scalar"]
        shown = len(scalar_tab._visible_orders)
        assert scalar_tab._order_count_label.text() == f"{shown} of {scalar_tab._order_list.count()} shown"

    def test_unchecking_an_order_removes_its_curve_without_a_resolve(self, win):
        from PySide6.QtCore import Qt as QtCoreQt

        from gratinglab.gui.qt.scalar_tab import _ORDER_ROLE

        scalar_tab = win.tabs["scalar"]
        token_before = win._token
        item = scalar_tab._order_list.item(0)
        order = item.data(_ORDER_ROLE)
        assert order in scalar_tab._visible_orders

        item.setCheckState(QtCoreQt.CheckState.Unchecked)

        assert win._token == token_before, "toggling an order must not trigger a solve"
        assert order not in scalar_tab._visible_orders

    def test_the_none_button_hides_every_curve(self, win):
        scalar_tab = win.tabs["scalar"]
        scalar_tab._show_no_orders()
        assert scalar_tab._visible_orders == frozenset()
        assert scalar_tab._order_count_label.text().startswith("0 of")

    def test_the_all_button_shows_every_order_that_existed_last_solve(self, win):
        scalar_tab = win.tabs["scalar"]
        scalar_tab._show_no_orders()
        scalar_tab._show_all_orders()
        assert scalar_tab._visible_orders == scalar_tab._previous_orders

    def test_the_default_button_reapplies_the_default_rule(self, win):
        from gratinglab.gui.orders import default_visible, summarize

        scalar_tab = win.tabs["scalar"]
        scalar_tab._show_no_orders()
        scalar_tab._show_default_orders()
        assert scalar_tab._visible_orders == default_visible(summarize(scalar_tab._scan))

    def test_a_deselected_order_survives_a_resolve(self, qtbot, win):
        """The behaviour that makes the panel feel right: nudge a geometry
        field, re-solve, and an order you unchecked stays unchecked."""
        scalar_tab = win.tabs["scalar"]
        deselected = next(iter(scalar_tab._previous_orders))
        scalar_tab._visible_orders = scalar_tab._visible_orders - {deselected}
        scalar_tab._resync_order_checkboxes()

        win.geometry._fields["blaze_angle"].setText("12.0")
        resolve(qtbot, win)

        assert deselected not in win.tabs["scalar"]._visible_orders


class TestLayoutFloors:
    """Sizing regressions, which is the class of bug that has shipped here
    twice unnoticed because nothing asserted on a pixel."""

    def test_no_scroll_area_promises_less_width_than_its_content(self, win):
        """A minimumWidth below the inner widget's own minimumSizeHint is not
        a tight budget -- it is a constraint Qt cannot satisfy, and it answers
        with a horizontal scrollbar and a clipped control.

        Would have failed before M13-B: ScalarTab declared 240 against 276 px
        of content. Needs no laid-out geometry, so it is deterministic.
        """
        from PySide6.QtWidgets import QScrollArea

        areas = [a for a in win.findChildren(QScrollArea) if a.widget() is not None]
        assert areas, "vacuous unless the window actually has scroll areas"

        for area in areas:
            needed = (
                area.widget().minimumSizeHint().width()
                + area.verticalScrollBar().sizeHint().width()
                + 2 * area.frameWidth()
            )
            assert area.minimumWidth() >= needed, (
                f"{type(area.widget()).__name__} promises "
                f"{area.minimumWidth()} px but needs {needed}"
            )


class TestGeometryDock:
    """The geometry inputs, in a panel that can be got out of the way."""

    def test_the_geometry_panel_lives_in_a_dock(self, win):
        from PySide6.QtWidgets import QDockWidget

        assert isinstance(win.geometry_dock, QDockWidget)
        assert win.geometry_dock.widget().widget() is win.geometry

    def test_it_docks_left_and_refuses_top_and_bottom(self, win):
        """A tall column of form rows docked along the top becomes a
        1180-px-wide strip of three fields: allowed by Qt, useless here."""
        from PySide6.QtCore import Qt

        areas = win.geometry_dock.allowedAreas()
        assert areas & Qt.DockWidgetArea.LeftDockWidgetArea
        assert areas & Qt.DockWidgetArea.RightDockWidgetArea
        assert not (areas & Qt.DockWidgetArea.TopDockWidgetArea)
        assert not (areas & Qt.DockWidgetArea.BottomDockWidgetArea)

    def test_it_is_closable(self, win):
        from PySide6.QtWidgets import QDockWidget

        features = win.geometry_dock.features()
        assert features & QDockWidget.DockWidgetFeature.DockWidgetClosable

    def test_closing_it_leaves_a_way_back(self, win):
        """A closable control with no route to reopen it is a trap."""
        toggle = win._view_menu.actions()[0]
        win.geometry_dock.hide()
        assert win.geometry_dock.isHidden()

        toggle.trigger()
        assert not win.geometry_dock.isHidden()

    def test_the_view_menu_holds_exactly_the_docks_own_toggle(self, win):
        """The dock's own action, so its checked state and its label cannot
        drift from the dock -- which a hand-rolled checkable action would, the
        moment someone clicked the dock's own close button."""
        actions = win._view_menu.actions()
        assert len(actions) == 1
        assert actions[0] is win.geometry_dock.toggleViewAction()

    def test_the_toggle_is_reachable_by_keyboard(self, win):
        from PySide6.QtGui import QKeySequence

        assert win._view_menu.actions()[0].shortcut() == QKeySequence("Ctrl+G")

    def test_the_dock_is_outside_the_central_widget(self, win):
        """What lets closing it hand the reclaimed width to the tabs. (The
        central widget still carries the profile plot above the tabs until
        M13-C moves it into the geometry tab.)"""
        assert not win.centralWidget().isAncestorOf(win.geometry)


class TestGeometryTab:
    """The tab and its wiring. What the *drawing* contains is decided and
    tested in `tests/test_gui_diagram.py`, without a window."""

    def test_it_sits_between_setup_and_the_solver_tabs(self, win):
        titles = [win._tab_widget.tabText(i) for i in range(win._tab_widget.count())]
        assert titles[0] == "Setup"
        assert titles[1] == "Grating Geometry"
        assert "Scalar" in titles[2]

    def test_it_is_not_in_the_solver_tab_registry(self, win):
        """It implements none of the solve/cancel contract, so it must never
        be somewhere `solve(name)` could reach."""
        assert win.geometry_tab not in win.tabs.values()
        assert not hasattr(win.geometry_tab, "name")

    def test_it_draws_before_any_solve_completes(self, qtbot):
        """Geometry needs no solver, so it should be right immediately rather
        than showing an empty panel until a result exists."""
        from gratinglab.gui.qt.main_window import MainWindow

        window = MainWindow()
        qtbot.addWidget(window)
        assert window.geometry_tab._diagram is not None

    def test_what_it_draws_equals_what_the_pure_module_returns(self, win):
        """TestNotASecondSourceOfTruth, extended to the diagram. If the tab
        ever computed an angle itself, this drifts."""
        from gratinglab.gui import diagram

        tab = win.geometry_tab
        expected = diagram.build(
            tab._parsed.problem, tab._parsed.illumination, tab._diagram.wavelength
        )
        assert [m.order for m in tab._diagram.orders] == [m.order for m in expected.orders]
        assert [m.sin_beta for m in tab._diagram.orders] == [
            m.sin_beta for m in expected.orders
        ]
        assert [a.tag for a in tab._diagram.arrows] == [a.tag for a in expected.arrows]

    def test_the_slider_spans_the_scan_and_starts_in_the_middle(self, win):
        tab = win.geometry_tab
        assert tab._slider.maximum() == len(tab._parsed.wavelengths) - 1
        assert tab._slider.value() == len(tab._parsed.wavelengths) // 2

    def test_moving_the_slider_redraws_and_never_solves(self, win):
        """The geometry is already known; nothing here depends on a solver."""
        tab = win.geometry_tab
        before = tab._diagram.wavelength
        token = win._token

        tab._slider.setValue(tab._slider.value() + 20)

        assert tab._diagram.wavelength != before
        assert win._token == token, "moving the slider must not trigger a solve"

    def test_the_blaze_button_jumps_to_the_blaze_wavelength(self, win):
        """At the shipped defaults lambda_b(m=2) = 4.05 nm, inside the 1-5 nm
        scan."""
        tab = win.geometry_tab
        tab._blaze_button.click()
        assert 4.0 < tab._diagram.wavelength < 4.1

    def test_the_blaze_button_is_disabled_with_a_reason_when_it_cannot_apply(
        self, qtbot, win
    ):
        """A dead control that explains itself, not a dead control. An enabled
        one that silently did nothing is the mistake M9 named."""
        win.geometry._fields["profile_kind"].setCurrentText("Sinusoidal")
        win._refresh_geometry_tab()

        tab = win.geometry_tab
        assert not tab._blaze_button.isEnabled()
        assert "no blaze angle" in tab._blaze_button.toolTip()

    def test_the_blaze_button_is_enabled_for_a_blazed_profile(self, win):
        """Non-vacuity for the test above."""
        assert win.geometry_tab._blaze_button.isEnabled()
        assert win.geometry_tab._blaze_button.toolTip() == ""

    def test_the_profile_plot_lives_here_now(self, win):
        """It was above the tabs, taking 876x345 -- more pixels than the whole
        geometry canvas got. A groove's shape is geometry, so it belongs in
        this tab."""
        assert win.geometry_tab.isAncestorOf(win.geometry_tab.profile_panel)
        assert not hasattr(win, "profile_panel")

    def test_the_hero_canvas_is_larger_than_the_side_panels(self, qtbot, win):
        """The point of the restructure. Needs real geometry, so the window is
        shown and laid out."""
        win.resize(1180, 820)
        win.show()
        qtbot.waitExposed(win)
        win._tab_widget.setCurrentWidget(win.geometry_tab)
        qtbot.wait(50)

        tab = win.geometry_tab
        hero = tab._scene3d._canvas.size()
        side = tab._canvas.size()
        assert hero.width() > side.width()
        assert hero.width() * hero.height() > 2 * side.width() * side.height()

    def test_the_hero_is_a_real_3d_axes(self, win):
        assert win.geometry_tab._scene3d._axes.name == "3d"

    def test_the_3d_scene_matches_what_the_pure_module_returns(self, win):
        """TestNotASecondSourceOfTruth, extended to 3D. If the panel ever
        computed a direction itself, this drifts."""
        from gratinglab.gui import diagram3d

        tab = win.geometry_tab
        expected = diagram3d.build_scene(
            tab._parsed.problem, tab._parsed.illumination, tab._scene.wavelength
        )
        drawn = {r.order: r.direction for r in tab._scene.rays if r.order is not None}
        wanted = {r.order: r.direction for r in expected.rays if r.order is not None}
        assert drawn == wanted

    def test_the_3d_view_offers_no_pan_or_zoom(self, win):
        """A NavigationToolbar2QT zoom rescales an Axes3D's data limits, which
        breaks the cube box aspect the true angles depend on. A control that
        can make the picture lie must not be offered."""
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT

        assert win.geometry_tab._scene3d.findChildren(NavigationToolbar2QT) == []

    def test_but_the_2d_figures_do_keep_theirs(self, win):
        """Non-vacuity for the test above: pan and zoom cannot distort an
        angle the 2D panels are not claiming."""
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT

        assert win.geometry_tab.profile_panel.findChildren(NavigationToolbar2QT)

    def test_the_3d_box_is_a_cube(self, win):
        """Without equal box aspect every angle in the view is a lie."""
        tab = win.geometry_tab
        tab._scene3d.show_scene(tab._scene)
        aspect = tab._scene3d._axes.get_box_aspect()
        assert aspect[0] == pytest.approx(aspect[1]) == pytest.approx(aspect[2])

    def test_the_captions_reach_the_panel(self, win):
        text = win.geometry_tab._captions.toPlainText()
        assert "groove axis" in text
        assert "evanescent" in text


class TestRimFocusAndPresets:
    """M13-G. All on the shared `win` fixture -- test_gui_qt.py is already 76%
    of suite wall time, so no test here builds its own window."""

    def _panel(self, win):
        """The tab defers `_paint` to `showEvent` and the shared fixture never
        shows it, so paint once on demand rather than making the tab current
        (which would cost every test in this class a full relayout)."""
        panel = win.geometry_tab._scene3d
        if not panel._artists["fine"]:
            panel.show_scene(win.geometry_tab._scene)
        return panel

    def test_the_preset_round_trips_back_to_the_physical_direction(self, win):
        """`_view_angles` inverts `diagram3d.view_direction` after the display
        permutation. Setting the preset and reading the camera back recovers
        the vector the pure module declared."""
        from gratinglab.gui import diagram3d
        from gratinglab.gui.qt.scene3d_panel import _display

        panel = self._panel(win)
        panel._preset_box.setCurrentText(diagram3d.RIM_PRESET)
        looking = diagram3d.view_direction(panel._axes.elev, panel._axes.azim)
        wanted = np.asarray(_display(diagram3d.PRESET_VIEWS[diagram3d.RIM_PRESET]))
        assert looking == pytest.approx(wanted, abs=1e-9)

    def test_and_that_direction_is_the_groove_axis(self, win):
        """Non-vacuity: the round-trip above would pass on any preset, so pin
        which one this is."""
        from gratinglab.gui import diagram3d

        assert diagram3d.PRESET_VIEWS[diagram3d.RIM_PRESET] == pytest.approx(
            diagram3d.G_HAT, abs=1e-12
        )

    def test_choosing_it_also_focuses_on_the_rim(self, win):
        """The pairing is the entire point; separating them would make the one
        useful combination a two-step discovery."""
        from gratinglab.gui import diagram3d

        panel = self._panel(win)
        panel._focus_button.setChecked(False)
        panel._preset_box.setCurrentText(diagram3d.RIM_PRESET)
        assert panel._focus_button.isChecked()
        assert panel._focus == "rim"

    def test_another_preset_does_not(self, win):
        """Non-vacuity for the test above."""
        panel = self._panel(win)
        panel._focus_button.setChecked(False)
        panel._preset_box.setCurrentText("Oblique")
        assert not panel._focus_button.isChecked()

    def test_focusing_shrinks_the_frame_by_the_stated_magnification(self, win):
        panel = self._panel(win)
        scene = win.geometry_tab._scene

        panel._focus_button.setChecked(False)
        wide = panel._axes.get_xlim3d()
        panel._focus_button.setChecked(True)
        close = panel._axes.get_xlim3d()

        assert (wide[1] - wide[0]) / (close[1] - close[0]) == pytest.approx(
            scene.rim_magnification, rel=1e-9
        )

    def test_and_says_so_on_the_canvas(self, win):
        panel = self._panel(win)
        panel._focus_button.setChecked(True)
        assert "27" in panel._magnification.text()

    def test_dragging_swaps_in_the_coarse_artists(self, win):
        """Both levels of detail are built once by `build_scene`; the drag
        flips visibility and resamples nothing."""
        from matplotlib.backend_bases import MouseButton, MouseEvent

        panel = self._panel(win)
        panel._focus_button.setChecked(False)
        assert panel._lod == "fine"

        press = MouseEvent(
            "button_press_event", panel._canvas, 200, 200, MouseButton.LEFT
        )
        press.inaxes = panel._axes
        panel._on_press(press)
        assert panel._lod == "coarse"
        assert any(a.get_visible() for a in panel._artists["coarse"])
        assert not any(a.get_visible() for a in panel._artists["fine"])

        panel._on_release(None)
        assert panel._lod == "fine"
        assert any(a.get_visible() for a in panel._artists["fine"])

    def test_the_two_artist_groups_are_different_objects(self, win):
        """Companion: without this the swap test would pass on an
        implementation that swapped nothing."""
        panel = self._panel(win)
        fine, coarse = panel._artists["fine"], panel._artists["coarse"]
        assert fine and coarse
        assert not (set(map(id, fine)) & set(map(id, coarse)))

    def test_rim_focus_is_offered_at_grazing_incidence(self, win):
        panel = self._panel(win)
        assert panel._focus_button.isEnabled()
        assert panel._focus_button.toolTip() == ""

    def test_and_refused_with_a_reason_in_plane(self, win):
        """The `blaze_jump` contract: a disabled control explains itself. At
        gamma = 90 the rim already is the whole scene."""
        from gratinglab.gui import diagram3d
        from gratinglab.illumination import Illumination

        tab = win.geometry_tab
        panel = self._panel(win)
        scene = diagram3d.build_scene(
            tab._parsed.problem,
            Illumination.classical(alpha=25.0, polarization="unpolarized"),
            600.0,
        )
        panel.show_scene(scene)
        assert not panel._focus_button.isEnabled()
        assert panel._focus_button.toolTip()

        panel.show_scene(tab._scene)  # leave the fixture as we found it


class TestSetupTab:
    """The stub, and the guard that keeps it honest."""

    def test_it_is_the_leftmost_tab(self, win):
        assert win._tab_widget.tabText(0) == "Setup"

    def test_a_solver_tab_is_what_is_actually_selected_on_open(self, win):
        """Regression, and a real bug caught by this milestone: adding Setup
        at index 0 made it the *current* tab, so `_solve_active_tab` found no
        `name` on it and the construction-time solve silently never ran --
        every test then sat out its full waitSignal timeout.

        Setup is leftmost because it is easy to find there; it is not what
        should be showing when the window opens.
        """
        current = win._tab_widget.currentWidget()
        assert getattr(current, "name", None) == "scalar"

    def test_it_explains_what_will_live_there(self, win):
        text = win._setup_tab._browser.toPlainText()
        assert "materials layer" in text
        assert "roadmap" in text

    def test_it_says_why_no_control_is_offered_yet(self, win):
        """The project's own rule, stated where a user would otherwise
        wonder why the tab is empty."""
        text = win._setup_tab._browser.toPlainText()
        assert "silently does nothing" in text
        assert "relative" in text

    def test_it_offers_no_input_at_all(self, win):
        """The guard. `Problem.coating` exists and is inert; a coating field
        here would look like a working option that changes nothing -- exactly
        the mistake M9 named and rejected. A future 'let's just add a field'
        fails here rather than slipping through review.
        """
        from PySide6.QtWidgets import QComboBox, QLineEdit

        assert win._setup_tab.findChildren(QLineEdit) == []
        assert win._setup_tab.findChildren(QComboBox) == []

    def test_it_is_not_in_the_solver_tab_registry(self, win):
        """Setup implements none of the solve/cancel contract, so it must
        never be somewhere `solve(name)` could reach."""
        assert win._setup_tab not in win.tabs.values()
        assert "setup" not in win.tabs

    def test_it_renders_as_typeset_markdown_not_raw_source(self, win):
        html = win._setup_tab._browser.toHtml()
        assert "<h1" in html
        assert "# Setup" not in win._setup_tab._browser.toPlainText()


class TestHelpMenu:
    """Read through `win._help_menu`, never `menuBar().actions()[i].menu()`.

    Under PySide6 6.11 a QMenu reached back down from the menu bar is not
    reliably usable -- it raises "Internal C++ object already deleted" while
    the menu the window built is alive and correct. Holding the QMenuBar does
    not help; holding the QMenu does. `_build_menu` therefore keeps both, and
    these tests use what the window kept.
    """

    def _labels(self, win, *, separators=False):
        return [
            action.text()
            for action in win._help_menu.actions()
            if separators or not action.isSeparator()
        ]

    def test_the_menu_bar_holds_exactly_view_and_help(self, win):
        """Two menus, and only two. View exists solely because the geometry
        dock is closable and a control with no way back is a trap; Help is
        still the one that matters."""
        titles = [action.text() for action in win._menu_bar.actions()]
        assert titles == ["&View", "&Help"]

    def test_lists_every_registered_solver_and_general_page(self, win):
        from gratinglab.gui.docs import general_pages, theory_pages
        from gratinglab.gui.richtext import menu_label

        labels = self._labels(win)
        for page in list(general_pages()) + list(theory_pages()):
            assert any(menu_label(page.title) in label for label in labels), page.title

    def test_an_ampersand_in_a_title_is_escaped_not_eaten(self, win):
        """Qt reads `&` as a mnemonic, so 'Grating Geometry & Conventions'
        would appear as 'Grating Geometry Conventions' with C underlined."""
        assert any("&&" in label for label in self._labels(win))

    def test_the_escaping_test_is_not_vacuous(self):
        """It only means something while a real page title contains `&`."""
        from gratinglab.gui.docs import general_pages

        assert any("&" in page.title for page in general_pages())

    def test_about_reports_the_real_version(self, win, monkeypatch):
        from gratinglab import __version__

        seen = {}
        monkeypatch.setattr(
            "PySide6.QtWidgets.QMessageBox.about",
            lambda parent, title, text: seen.update(title=title, text=text),
        )
        win.show_about()
        assert seen["title"] == "About GratingLab"
        assert __version__ in seen["text"]


class TestTheoryViewer:
    def _open(self, qtbot, win, name="scalar"):
        from gratinglab.gui.docs import theory_pages

        page = next(p for p in theory_pages() if p.name == name)
        viewer = win.show_theory(page)
        qtbot.addWidget(viewer)
        return viewer

    def test_opens_a_page_with_typeset_math(self, qtbot, win):
        text = self._open(qtbot, win)._browser.toPlainText()
        assert "5. Energy is not conserved" in text
        assert "$" not in text
        assert "\\Phi_m" not in text

    def test_the_equations_are_registered_as_images(self, qtbot, win):
        """The `Text.dump()` image count this replaces, expressed against the
        document's own resources."""
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QTextDocument

        document = self._open(qtbot, win)._browser.document()
        loaded = document.resource(
            QTextDocument.ResourceType.ImageResource, QUrl("gratinglab-math:0")
        )
        assert not loaded.isNull()

    def test_a_table_survives_into_the_document(self, qtbot, win):
        """What this milestone bought: scalar.md section 5 is a table, and it
        used to render as literal pipe characters."""
        assert "<table" in self._open(qtbot, win)._browser.toHtml()

    def test_a_non_rigorous_solver_is_bannered(self, qtbot, win):
        from gratinglab.gui.docs import theory_pages

        page = next(p for p in theory_pages() if p.name == "scalar")
        assert not page.rigorous  # the case this test needs
        assert "Approximate method" in self._open(qtbot, win)._browser.toPlainText()

    def test_a_page_that_does_not_exist_yet_explains_itself(self, qtbot, win):
        from gratinglab.gui.docs import TheoryPage

        stub = TheoryPage(
            name="future", title="Future Method", available=False, path=None,
            text="No theory page has been written yet for 'future'.",
            rigorous=True,
        )
        viewer = win.show_theory(stub)  # must not raise
        qtbot.addWidget(viewer)
        text = viewer._browser.toPlainText()
        assert "No theory page has been written yet" in text
        assert "Approximate method" not in text


class TestToolkitBoundary:
    def test_only_gui_qt_may_import_a_toolkit(self):
        """A package boundary rather than a convention, so the drift that
        always happens -- someone needing 'just a QColor' in a pure module --
        fails here instead of quietly spreading.

        `app.py` is the one file at this level that's expected to *reach into*
        the toolkit -- but only deferred, inside `_require_qt()`, which is
        exactly what keeps `import gratinglab.gui.app` safe without the `gui`
        extra. So this checks module-*level* imports only (`.body`, not a full
        recursive walk) -- a recursive walk would flag that sanctioned
        deferred import as if it were the unconditional kind this test exists
        to catch.
        """
        import ast
        from pathlib import Path

        import gratinglab.gui as gui

        root = Path(gui.__file__).parent
        offenders = {}
        for path in sorted(root.glob("*.py")):
            imported: set[str] = set()
            for node in ast.parse(path.read_text()).body:
                if isinstance(node, ast.Import):
                    imported.update(a.name.split(".")[0] for a in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0])
            found = {"PySide6", "PyQt6", "PyQt5", "tkinter"} & imported
            if found:
                offenders[path.name] = found

        assert offenders == {}, "only gui/qt/ may import a toolkit unconditionally"


class TestEntryPoint:
    """`app.py`'s own two jobs: explain a missing toolkit, and stay safe to
    import without one.

    Ported from the Tk suite this replaces -- PySide6 is pip-installable
    (unlike Tk), so the *message* changed, but the failure mode is still real:
    an install of `[dev]` without `[gui]` hits exactly this.
    """

    def test_gives_an_explanation_not_a_traceback(self, monkeypatch):
        import builtins

        from gratinglab.gui import app as app_module

        real_import = builtins.__import__

        def refuse(name, *args, **kwargs):
            if name == "PySide6.QtWidgets" or name.startswith("PySide6"):
                raise ImportError("No module named 'PySide6'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", refuse)

        with pytest.raises(SystemExit) as excinfo:
            app_module._require_qt()

        message = str(excinfo.value)
        assert 'pip install -e ".[gui]"' in message
        assert "Everything except the GUI works without it" in message

    def test_importing_the_module_does_not_require_pyside6(self):
        """Import must be safe where PySide6 is missing; only calling
        `_require_qt`/`main` needs it."""
        import ast
        import inspect

        from gratinglab.gui import app as app_module

        top_level = {
            node.names[0].name
            for node in ast.parse(inspect.getsource(app_module)).body
            if isinstance(node, ast.Import)
        }
        assert "PySide6" not in top_level
