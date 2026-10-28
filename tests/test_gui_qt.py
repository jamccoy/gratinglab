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

from endstation.testing.boundary import assert_qt_free, modules_in

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


class TestTheIntegralTab:
    """The second solver's tab, and the boundary condition it now offers.

    The tab existed before `conductivity` did, so it asserted a perfectly
    conducting boundary in a static label while the solver had grown a second
    one. These pin the selector, and pin that nothing downstream renames the
    quantity it is showing.
    """

    #: A scan short enough to run twice in a test and coarse enough to clear
    #: the metal-side mesh floor for Au. The point here is the wiring, not the
    #: physics -- `test_gui_integral_options.py` owns the guards.
    TABULATED_FORM = dict(
        coating="Au", wavelength_start="4.0", wavelength_stop="4.5",
        wavelength_count="2",
    )
    BOUNDARY_POINTS = "128"

    def test_the_conductivity_selector_offers_the_solvers_modes(self, win):
        from gratinglab.gui.integral_options import CONDUCTIVITY_MODES

        combo = win.tabs["integral"]._fields["conductivity"]
        offered = [combo.itemText(i) for i in range(combo.count())]
        assert tuple(offered) == CONDUCTIVITY_MODES
        assert combo.currentText() == "perfect"

    def test_the_note_follows_the_selector(self, win):
        """The label it replaced was a claim about the result that a control
        beside it could falsify."""
        tab = win.tabs["integral"]
        before = tab._conductivity_note.text()
        assert "Perfectly conducting" in before

        tab._fields["conductivity"].setCurrentText("tabulated")
        after = tab._conductivity_note.text()
        assert after != before
        assert "absolute" in after

    def test_the_selector_reaches_build_options(self, win):
        tab = win.tabs["integral"]
        tab._fields["conductivity"].setCurrentText("tabulated")
        parsed = build(FormState(coating="Au"))
        options = tab.build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths
        )
        assert options["conductivity"] == "tabulated"

    def test_the_corner_mesh_selector_offers_off_and_the_orders(self, win):
        from gratinglab.gui.integral_options import GRADING_CHOICES

        combo = win.tabs["integral"]._fields["corner_grading"]
        assert tuple(combo.itemText(i) for i in range(combo.count())) == GRADING_CHOICES
        assert combo.currentText() == "off"

    def test_the_corner_mesh_note_follows_the_selector(self, win):
        tab = win.tabs["integral"]
        assert "Equal-arc" in tab._grading_note.text()
        tab._fields["corner_grading"].setCurrentText("4")
        assert "order 4" in tab._grading_note.text()

    def test_the_corner_mesh_reaches_build_options(self, win):
        """Regression: the field and the state agreed on the name, so the
        consistency check passed while `build_options` never read it, and the
        selector did nothing."""
        tab = win.tabs["integral"]
        tab._fields["corner_grading"].setCurrentText("4")
        tab._fields["boundary_points"].setText("1200")
        parsed = build(FormState())
        options = tab.build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths
        )
        assert options["corner_grading"] == 4

    def test_a_tabulated_solve_equals_a_direct_one(self, qtbot, win):
        """`TestNotASecondSourceOfTruth`, for the coupled system: the tab must
        not compute, reshape or rescale anything on the way to the plot."""
        from gratinglab.gui.integral_options import (
            IntegralOptionsState,
            build_options as build_integral_options,
        )
        from gratinglab.solvers import integral

        for field, value in self.TABULATED_FORM.items():
            widget = win.geometry._fields[field]
            if hasattr(widget, "setCurrentIndex"):
                widget.setCurrentIndex(widget.findData(value))
            else:
                widget.setText(value)
        tab = win.tabs["integral"]
        tab._fields["boundary_points"].setText(self.BOUNDARY_POINTS)
        tab._fields["conductivity"].setCurrentText("tabulated")
        resolve(qtbot, win, "integral")

        parsed = build(FormState(**self.TABULATED_FORM))
        options = build_integral_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            IntegralOptionsState(
                boundary_points=self.BOUNDARY_POINTS, conductivity="tabulated"
            ),
        )
        expected = integral.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(tab._scan.efficiency, expected.efficiency)
        assert np.array_equal(tab._scan.absorption, expected.absorption)

    def test_the_panel_names_the_material_and_the_conserved_quantity(
        self, qtbot, win
    ):
        """A finite-conductivity run whose Σ sits at 0.7 is not a solver
        losing power -- the rest is absorbed, and the panel has to say so
        rather than leaving a correct result looking like a deficit."""
        for field, value in self.TABULATED_FORM.items():
            widget = win.geometry._fields[field]
            if hasattr(widget, "setCurrentIndex"):
                widget.setCurrentIndex(widget.findData(value))
            else:
                widget.setText(value)
        tab = win.tabs["integral"]
        tab._fields["boundary_points"].setText(self.BOUNDARY_POINTS)
        tab._fields["conductivity"].setCurrentText("tabulated")
        resolve(qtbot, win, "integral")

        shown = tab._provenance.toPlainText()
        assert "boundary points" in shown
        assert "absolute (Au)" in shown
        assert "absorption: A" in shown
        assert "Σ + A" in shown

    def test_without_a_coating_tabulated_is_refused_on_the_form(self, win):
        """And refused before anything reaches the worker: the previous
        result stays on screen and the panel names the field."""
        tab = win.tabs["integral"]
        tab._fields["conductivity"].setCurrentText("tabulated")
        win.solve("integral")
        shown = tab._provenance.toPlainText()
        assert "coating" in shown
        assert "field(s) need attention" in shown


class TestTheRcwaTab:
    """The third solver's tab. The wiring only: `test_gui_rcwa_options.py`
    owns the guards and `test_rcwa.py` the physics."""

    #: Soft X-ray, off-plane, lamellar, two wavelengths, Au: the bundled table,
    #: and small enough to solve twice in a test.
    FORM = dict(
        profile_kind="Lamellar", coating="Au",
        wavelength_start="4.0", wavelength_stop="4.5", wavelength_count="2",
    )
    TRUNCATION = "10"

    def set_geometry(self, win, form):
        for field, value in form.items():
            widget = win.geometry._fields[field]
            if hasattr(widget, "setCurrentIndex"):
                index = widget.findData(value)
                widget.setCurrentIndex(index if index >= 0 else widget.findText(value))
            else:
                widget.setText(value)

    def test_it_is_a_tab_with_its_short_title(self, win):
        from gratinglab.gui.docs import display_title, short_title

        tab = win.tabs["rcwa"]
        index = win._tab_widget.indexOf(tab)
        assert index >= 0
        assert win._tab_widget.tabText(index) == short_title("rcwa")
        assert win._tab_widget.tabToolTip(index) == display_title("rcwa")

    def test_a_solve_equals_a_direct_one(self, qtbot, win):
        """`TestNotASecondSourceOfTruth`, for RCWA: the tab must not compute,
        reshape or rescale anything on the way to the plot."""
        from gratinglab.gui.rcwa_options import RcwaOptionsState
        from gratinglab.gui.rcwa_options import build_options as build_rcwa_options
        from gratinglab.solvers import rcwa

        self.set_geometry(win, self.FORM)
        tab = win.tabs["rcwa"]
        tab._fields["truncation"].setText(self.TRUNCATION)
        resolve(qtbot, win, "rcwa")

        parsed = build(FormState(**self.FORM))
        options = build_rcwa_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            RcwaOptionsState(truncation=self.TRUNCATION),
        )
        expected = rcwa.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(tab._scan.efficiency, expected.efficiency)
        assert np.array_equal(tab._scan.absorption, expected.absorption)
        shown = tab._provenance.toPlainText()
        assert "absolute (Au)" in shown
        assert "exact" in tab._profile_note.text()

    def test_the_profile_note_follows_the_geometry(self, qtbot, win):
        self.set_geometry(win, {**self.FORM, "profile_kind": "Blazed"})
        tab = win.tabs["rcwa"]
        tab._fields["truncation"].setText(self.TRUNCATION)
        tab._fields["slices"].setText("6")
        resolve(qtbot, win, "rcwa")
        assert "6 slabs" in tab._profile_note.text()

    def test_without_a_coating_it_is_refused_on_the_form(self, win):
        tab = win.tabs["rcwa"]
        win.solve("rcwa")
        shown = tab._provenance.toPlainText()
        assert "coating" in shown
        assert "field(s) need attention" in shown


class TestTheLayersField:
    """Layers on the geometry panel: the RCWA tab models them, the others
    refuse them on the field rather than from the worker."""

    FORM = dict(
        profile_kind="Lamellar", coating="Au",
        wavelength_start="4.0", wavelength_stop="4.5", wavelength_count="2",
    )

    def set_geometry(self, win, form):
        for field, value in form.items():
            widget = win.geometry._fields[field]
            if hasattr(widget, "setCurrentIndex"):
                index = widget.findData(value)
                widget.setCurrentIndex(index if index >= 0 else widget.findText(value))
            else:
                widget.setText(value)

    def test_it_reaches_the_problem(self, win):
        self.set_geometry(win, {**self.FORM, "layers": "Au 2"})
        problem = build(win.geometry.read_form()).problem
        assert [(x.material, x.thickness) for x in problem.layers] == [("Au", 2.0)]

    @pytest.mark.parametrize("name", ["scalar", "integral"])
    def test_a_solver_without_layers_refuses_on_the_field(self, win, name):
        self.set_geometry(win, {**self.FORM, "layers": "Au 2"})
        before = win.tabs[name]._scan
        win.solve(name)
        shown = win.tabs[name]._provenance.toPlainText()
        assert "layers" in shown and "cannot model conformal layers" in shown
        assert win.tabs[name]._scan is before  # nothing reached the worker

    def test_the_rcwa_tab_solves_them_as_a_direct_solve_would(self, qtbot, win):
        from gratinglab.gui.rcwa_options import RcwaOptionsState
        from gratinglab.gui.rcwa_options import build_options as build_rcwa_options
        from gratinglab.solvers import rcwa

        form = {**self.FORM, "layers": "Au 2"}
        self.set_geometry(win, form)
        tab = win.tabs["rcwa"]
        tab._fields["truncation"].setText("8")
        tab._fields["slices"].setText("4")
        resolve(qtbot, win, "rcwa")

        parsed = build(FormState(**form))
        options = build_rcwa_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            RcwaOptionsState(truncation="8", slices="4"),
        )
        expected = rcwa.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(tab._scan.efficiency, expected.efficiency)
        assert tab._scan.provenance.notes["layers"] == ["Au 2 nm"]
        assert "staircase" in tab._profile_note.text()


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

        def record(self, token, method, geometry, options, cancel):
            seen["thread"] = threading.current_thread().ident
            return original_run(self, token, method, geometry, options, cancel)

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(worker_module.SolveWorker, "run", record)
            resolve(qtbot, win)

        assert "thread" in seen, "the worker slot never ran"
        assert seen["thread"] != threading.current_thread().ident

    def test_a_two_sided_theorem_is_not_checked_at_machine_precision(self):
        """`R + A = 1` is pinned from both sides, so its discretisation error
        lands above unity as often as below. A converged finite-conductivity
        run measured +2.7e-6 and rendered as "EXCEEDS UNITY" in red -- the
        exact "nothing correct looks wrong" failure `gui/provenance.py` exists
        to prevent, reachable only once the tabulated mode was selectable."""
        from gratinglab.gui.qt.worker import _energy
        from gratinglab.result import EfficiencyScan, Provenance

        def scan(total, absorption):
            return EfficiencyScan(
                wavelengths=np.array([1.0]),
                orders=np.array([0]),
                efficiency=np.array([[total]]),
                propagating=np.array([[True]]),
                provenance=Provenance("integral"),
                absorption=None if absorption is None else np.array([absorption]),
            )

        assert _energy(scan(0.7, 0.3 + 2.7e-6)).passed
        # Non-vacuity in both directions: a real violation still fails, and a
        # one-sided quantity keeps the strict default, since an excess there is
        # a modelling error rather than rounding.
        assert not _energy(scan(0.7, 0.8)).passed
        assert not _energy(scan(1.0 + 2.7e-6, None)).passed

    def test_the_worker_thread_gets_a_full_stack(self, win):
        """A QThread starts with a far smaller stack than the main thread's,
        and LAPACK's dense complex solve overruns it -- `np.linalg.solve` on a
        few-hundred-square matrix took down the whole process, so the integral
        tab returned nothing in the real window under either boundary
        condition while every headless test passed.

        Pinned as a number as well as behaviourally (`TestTheIntegralTab`
        actually solves through the worker) because the behavioural failure is
        a segfault, not an assertion: it takes the suite with it rather than
        reporting, so it is worth catching here first."""
        from gratinglab.gui.qt.main_window import _WORKER_STACK_BYTES

        assert _WORKER_STACK_BYTES >= 8 * 1024 * 1024
        assert win._thread.stackSize() == _WORKER_STACK_BYTES

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

    def test_cancel_says_the_calculation_stopped(self, win):
        """It used to say "still finishing", which was true when cancelling
        could only stop the waiting. The solver now checks at every
        wavelength, so that sentence would be the panel's first lie."""
        win.solve("scalar")
        win.cancel()
        text = win.tabs["scalar"]._provenance.toPlainText()
        assert "the calculation stopped" in text
        assert "still finishing" not in text

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


class TestCancellationReallyStops:
    """M14-D. Before this, Cancel stopped the *waiting* and the CPU kept
    going. These measure the difference rather than asserting it."""

    def _long_solve(self, win, points="4000"):
        """A scan slow enough that Cancel lands mid-flight."""
        win.geometry._fields["wavelength_count"].setText(points)

    def test_cancel_sets_the_event_the_solver_reads(self, win):
        win.solve("scalar")
        assert not win._cancel.is_set()
        win.cancel()
        assert win._cancel.is_set()

    def test_and_the_worker_stops_where_it_was_told(self, qtbot, win):
        """The measured half: it stops *mid-flight*, not merely before it
        started.

        Waiting for real progress first is what makes this a test of
        cancellation rather than of racing the thread to the starting line --
        a fixed sleep let the flag be set before the worker's first callback,
        and the assertion passed on a solve that never began.
        """
        reached = []
        stopped = []
        win._worker.progress.connect(lambda t, d, n: reached.append(d))
        win._worker.cancelled.connect(lambda t, m: stopped.append(t))

        self._long_solve(win)
        win.solve("scalar")
        qtbot.waitUntil(lambda: any(d > 0 for d in reached), timeout=SOLVE_TIMEOUT_MS)
        under_way = max(reached)
        win.cancel()

        qtbot.waitUntil(lambda: bool(stopped), timeout=SOLVE_TIMEOUT_MS)
        assert under_way > 0, "it had not started, so stopping proves nothing"
        # Reports are throttled to whole percent, so the last one can trail the
        # true stopping point by up to 40 wavelengths here. Half the scan is
        # far outside that and still fails loudly if nothing was interrupted.
        assert max(reached) < 2000, (
            f"reached {max(reached)} of 4000 -- it ran on despite being cancelled"
        )

    def test_a_solve_left_alone_runs_to_the_end(self, qtbot, win):
        """Non-vacuity for the test above: without it, a solver that stopped
        early for some entirely different reason would pass."""
        reached = []
        win._worker.progress.connect(lambda t, d, n: reached.append(d))
        resolve(qtbot, win)
        assert max(reached) == 200

    def test_cancelling_is_reported_as_cancelled_not_failed(self, qtbot, win):
        """A stop the user asked for is not an error, and rendering it as one
        would be the category mistake M8 fixed for the no-coating default."""
        events = []
        win._worker.cancelled.connect(lambda t, m: events.append("cancelled"))
        win._worker.failed.connect(lambda t, m, msg: events.append("failed"))

        self._long_solve(win)
        win.solve("scalar")
        qtbot.wait(30)
        win.cancel()

        qtbot.waitUntil(lambda: bool(events), timeout=SOLVE_TIMEOUT_MS)
        assert events == ["cancelled"]

    def test_but_a_real_solver_failure_still_is(self, win):
        """Non-vacuity: the distinction only means something if `failed` can
        still happen."""
        seen = []
        win._worker.failed.connect(lambda t, m, msg: seen.append(msg))
        win._active_name = "scalar"
        win._on_failed(win._token, "scalar", "ValueError: boom")
        assert "solve failed" in win.tabs["scalar"]._provenance.toPlainText()

    def test_each_solve_gets_its_own_event(self, win):
        """A worker-owned event cleared at the top of `run` loses a Cancel
        clicked between the request and the clear. A fresh one per solve
        cannot."""
        win.solve("scalar")
        first = win._cancel
        win.cancel()
        win.solve("scalar")
        assert win._cancel is not first
        assert not win._cancel.is_set()


class TestTheProgressBar:
    def test_it_is_determinate_once_the_solver_reports(self, qtbot, win):
        resolve(qtbot, win)
        bar = win.tabs["scalar"]._progress
        assert bar.maximum() == 200

    def test_and_indeterminate_again_for_the_next_run(self, qtbot, win):
        """The previous run's range is meaningless for a solve that has not
        started, and leaving it would show a full bar for no work."""
        resolve(qtbot, win)
        win.tabs["scalar"].show_solving(200)
        assert win.tabs["scalar"]._progress.maximum() == 0

    def test_reports_are_throttled_to_whole_percent(self, qtbot, win):
        """200 queued cross-thread signals for a bar with ~100 visible states
        is 100 of them nobody could have seen."""
        seen = []
        win._worker.progress.connect(lambda t, d, n: seen.append(d))
        resolve(qtbot, win)
        assert len(seen) <= 101
        assert len(seen) > 1, "and it does report"

    def test_a_stale_report_does_not_move_the_bar(self, win):
        """A cancelled solve keeps reporting for the wavelength it is in.
        Letting that walk the bar over a result the window has moved on from
        would be the same staleness the token exists to catch."""
        win._active_name = "scalar"
        bar = win.tabs["scalar"]._progress
        bar.setRange(0, 100)
        bar.setValue(50)
        win._on_progress(win._token - 1, 7, 100)
        assert bar.value() == 50


class TestTheConvergenceButton:
    """M14-E. The sweep is the first operation here long enough for a progress
    bar to be worth having, and the first a user is likely to want out of."""

    def _sweep(self, qtbot, win):
        with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
            win.tabs["scalar"]._converge_button.click()

    def test_it_runs_the_harness_and_the_panel_says_so(self, qtbot, win):
        """No new rendering: `provenance_lines` already reads the study off
        `notes["convergence"]`, which is what M3-C put there."""
        self._sweep(qtbot, win)
        text = win.tabs["scalar"]._provenance.toPlainText()
        assert "convergence: yes" in text
        assert "quadrature_points=4096 is enough" in text

    def test_a_plain_solve_still_says_it_did_not_sweep(self, qtbot, win):
        """Non-vacuity: pressing Solve must not silently sweep."""
        resolve(qtbot, win)
        assert "not checked" in win.tabs["scalar"]._provenance.toPlainText()

    def test_the_knob_field_ends_at_the_coarsest_defensible_setting(self, qtbot, win):
        """The actionable number, put where the next plain Solve will use it."""
        field = win.tabs["scalar"]._fields["quadrature_points"]
        assert field.text() == "2048"
        self._sweep(qtbot, win)
        assert field.text() == "4096"

    def test_the_knob_field_is_dead_while_the_sweep_chooses_it(self, win):
        """An input that looks live and is ignored is the M9 mistake."""
        field = win.tabs["scalar"]._fields["quadrature_points"]
        assert field.isEnabled()
        win.solve("scalar", sweep=True)
        assert not field.isEnabled()
        assert not win.tabs["scalar"]._converge_button.isEnabled()

    def test_and_comes_back_afterwards(self, qtbot, win):
        self._sweep(qtbot, win)
        assert win.tabs["scalar"]._fields["quadrature_points"].isEnabled()
        assert win.tabs["scalar"]._converge_button.isEnabled()

    def test_it_shares_the_one_at_a_time_rule(self, win):
        """Same token, same event, same window-wide lock as an ordinary solve
        -- it is the same kind of thing, only longer."""
        win.solve("scalar", sweep=True)
        assert win._running
        assert not win.tabs["scalar"]._solve_button.isEnabled()
        token = win._token
        win.solve("scalar")  # ignored while one is in flight
        assert win._token == token

    def test_the_bar_counts_rungs_not_wavelengths(self, qtbot, win):
        """`check_convergence` reports per rung: the numbers a bar can show
        without jumping backwards."""
        self._sweep(qtbot, win)
        assert win.tabs["scalar"]._progress.maximum() == 10

    def test_cancelling_a_sweep_stops_it(self, qtbot, win):
        stopped = []
        rungs = []
        win._worker.cancelled.connect(lambda t, m: stopped.append(t))
        win._worker.progress.connect(lambda t, d, n: rungs.append(d))

        win.solve("scalar", sweep=True)
        qtbot.waitUntil(lambda: bool(rungs), timeout=SOLVE_TIMEOUT_MS)
        win.cancel()
        qtbot.waitUntil(lambda: bool(stopped), timeout=SOLVE_TIMEOUT_MS)
        assert max(rungs) < 9, "it worked through the whole ladder anyway"

    def test_the_result_survives_a_cancelled_sweep(self, qtbot, win):
        previous = win.tabs["scalar"]._scan
        win.solve("scalar", sweep=True)
        win.cancel()
        assert win.tabs["scalar"]._scan is previous


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

    def test_the_view_menu_holds_the_docks_own_toggle_and_the_themes(self, win):
        """The dock's own action, so its checked state and its label cannot
        drift from the dock -- which a hand-rolled checkable action would, the
        moment someone clicked the dock's own close button. The only other
        entry is the Theme submenu, added with themes."""
        actions = [a for a in win._view_menu.actions() if not a.isSeparator()]
        assert actions[0] is win.geometry_dock.toggleViewAction()
        assert [a.menu() is not None for a in actions[1:]] == [True]
        assert actions[1].text() == "&Theme"

    def test_the_toggle_is_reachable_by_keyboard(self, win):
        from PySide6.QtGui import QKeySequence

        assert win._view_menu.actions()[0].shortcut() == QKeySequence("Ctrl+G")

    def test_the_dock_is_outside_the_central_widget(self, win):
        """What lets closing it hand the reclaimed width to the tabs. (The
        central widget still carries the profile plot above the tabs until
        M13-C moves it into the geometry tab.)"""
        assert not win.centralWidget().isAncestorOf(win.geometry)


class TestDockAlignment:
    """The four groups share one label column, so the inputs line up."""

    def test_every_label_has_one_width(self, win):
        from PySide6.QtWidgets import QFormLayout, QLabel

        widths = set()
        for form in win.geometry._forms:
            for row in range(form.rowCount()):
                item = form.itemAt(row, QFormLayout.ItemRole.LabelRole)
                if item is not None and isinstance(item.widget(), QLabel):
                    widths.add(item.widget().width() if item.widget().isVisible()
                               else item.widget().maximumWidth())
        assert len(widths) == 1

    def test_switching_mount_keeps_the_column(self, win):
        """Measured over every text a mount label can take, so relabelling
        the angle fields never widens the column."""
        label = win.geometry._angle_labels["alpha"]
        before = label.maximumWidth()
        for mount in ("Classical", "Conical", "Off-plane"):
            win.geometry._fields["mount"].setCurrentText(mount)
            assert label.maximumWidth() == before
            assert label.fontMetrics().horizontalAdvance(label.text()) <= before


class TestArrowLabels:
    def test_no_two_labels_in_the_cross_section_overlap(self, qtbot, win):
        """The order tips crowd one arc; `_declutter` pushes labels apart."""
        tab = win.geometry_tab
        # Shown and at a real window size: an unshown canvas keeps the
        # figure's default size, which leaves the labels room, and then this
        # test passes with or without the declutter (measured).
        win.resize(1400, 860)
        win.show()
        win._tab_widget.setCurrentWidget(tab)
        qtbot.waitExposed(win)
        qtbot.waitUntil(lambda: tab._diagram is not None, timeout=5000)
        qtbot.wait(200)
        tab._paint(tab._diagram)
        renderer = tab._canvas.get_renderer()
        axes = tab._axes["main"]
        boxes = [
            text.get_window_extent(renderer)
            for text in axes.texts
            if text.get_text().startswith(("m=", "β"))
        ]
        assert len(boxes) >= 3
        for i, a in enumerate(boxes):
            for b in boxes[i + 1:]:
                assert not a.overlaps(b)
        title = axes.title.get_window_extent(renderer)
        assert not any(box.overlaps(title) for box in boxes)


class TestGeometryTab:
    """The tab and its wiring. What the *drawing* contains is decided and
    tested in `tests/test_gui_diagram.py`, without a window."""

    def test_it_sits_between_setup_and_the_solver_tabs(self, win):
        from gratinglab.gui.docs import display_title, short_title

        titles = [win._tab_widget.tabText(i) for i in range(win._tab_widget.count())]
        assert titles[0] == "Setup"
        assert titles[1] == "Grating Geometry"
        # Everything after the boundary is a solver tab -- however many
        # solvers are registered, and in whatever order.
        # Short labels on the tabs, so the bar does not overflow as solvers
        # are added; the full name is the tooltip, and the Help menu's.
        assert titles[2:] == [short_title(name) for name in win.tabs]
        tips = [win._tab_widget.tabToolTip(i) for i in range(2, win._tab_widget.count())]
        assert tips == [display_title(name) for name in win.tabs]

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

        # Back to the real scene, so a failure above is reported against the
        # geometry the rest of this class assumes rather than against a
        # leftover in-plane one.
        panel.show_scene(tab._scene)
        assert panel._focus_button.isEnabled()


class TestLiveGeometry:
    """M13-H. The picture follows the form without a solve.

    Nothing here reaches a worker: `state.build` is microseconds of parsing
    and the redraw is matplotlib. The token is the proof -- it never moves.
    """

    def _edit(self, qtbot, win, key, text):
        """Type into a field the way a user does.

        `setText` fires `textChanged`, not `textEdited`, and the panel
        deliberately listens to the latter -- so this drives the real path
        rather than the one tests happen to find convenient.
        """
        field = win.geometry._fields[key]
        field.setText(text)
        field.textEdited.emit(text)
        qtbot.wait(250)  # past the 120 ms debounce

    def test_editing_a_field_redraws_the_geometry(self, qtbot, win):
        before = win.geometry_tab._parsed.problem.period
        self._edit(qtbot, win, "period", "500")
        assert win.geometry_tab._parsed.problem.period == pytest.approx(500.0)
        assert before != pytest.approx(500.0)

    def test_and_does_not_solve(self, qtbot, win):
        """The distinction the whole step rests on. A redraw that quietly
        started a solve would be the live re-solve ruled out of scope."""
        token = win._token
        self._edit(qtbot, win, "period", "480")
        assert win._token == token
        assert not win._running

    def test_the_3d_scene_follows_too(self, qtbot, win):
        """Non-vacuity in the other direction: `_parsed` moving proves the
        parse ran, not that anything was drawn from it."""
        self._edit(qtbot, win, "gamma", "12")
        assert win.geometry_tab._scene.gamma == pytest.approx(np.radians(12.0))

    def test_a_burst_of_keystrokes_costs_one_redraw(self, qtbot, win):
        """What the debounce is for. Without it, typing "315.15" draws six
        times -- each ~50 ms of matplotlib, while the user is still typing."""
        redraws = []
        original = win.geometry_tab.show_geometry

        def counted(parsed):
            redraws.append(parsed)
            original(parsed)

        win.geometry_tab.show_geometry = counted
        # Free earlier tests' windows now. Freeing one can take hundreds of ms,
        # and since gc_guard runs the collector on this thread, a collection
        # landing mid-burst would outlast the debounce and count as a redraw:
        # a measurement of garbage, not of the debounce.
        import gc

        gc.collect()
        try:
            field = win.geometry._fields["period"]
            for text in ("3", "31", "315", "315.", "315.1", "315.15"):
                field.setText(text)
                field.textEdited.emit(text)
                qtbot.wait(10)
            qtbot.wait(250)
        finally:
            win.geometry_tab.show_geometry = original
        assert len(redraws) == 1

    def test_a_half_typed_field_keeps_the_last_good_drawing(self, qtbot, win):
        """`"-"` fails to parse for exactly one keystroke. Blanking the canvas
        for it would style ordinary typing as a problem."""
        self._edit(qtbot, win, "period", "420")
        good = win.geometry_tab._scene

        self._edit(qtbot, win, "period", "-")
        assert win.geometry_tab._scene is good
        assert win.geometry_tab._parsed.problem.period == pytest.approx(420.0)

    def test_and_is_not_styled_as_an_error(self, qtbot, win):
        """The absent-assertion: `warn` and `bad` are reserved for something
        actually wrong. A form mid-edit is not that. The colours are the
        theme's, which is what the panel renders with."""
        from gratinglab.gui.qt import theming
        from gratinglab.gui.theme import provenance_colors

        colors = provenance_colors(theming.current())
        self._edit(qtbot, win, "period", "-")
        html = win.geometry_tab._captions.toHtml()
        assert colors["warn"] not in html
        assert colors["bad"] not in html
        assert "waiting" in win.geometry_tab._captions.toPlainText()
        # Its own line: `to_html` concatenates runs, so a missing newline
        # would run the note straight into the caption beneath it.
        assert win.geometry_tab._captions.toPlainText().splitlines()[0].endswith(
            "the drawing below is from the last complete geometry"
        )

    def test_and_the_captions_come_back_when_it_parses(self, qtbot, win):
        """Non-vacuity for the pending line: it has to clear."""
        self._edit(qtbot, win, "period", "-")
        assert "waiting" in win.geometry_tab._captions.toPlainText()
        self._edit(qtbot, win, "period", "315.15")
        assert "waiting" not in win.geometry_tab._captions.toPlainText()

    def test_programmatic_setText_schedules_nothing(self, qtbot, win):
        """Why the panel listens to `textEdited` rather than `textChanged`:
        six existing tests set fields programmatically, and every one of them
        would otherwise start a debounce timer."""
        win.geometry._fields["period"].setText("999")
        assert not win._redraw_timer.isActive()
        win.geometry._fields["period"].setText("315.15")


class TestSetupTab:
    """No longer a stub. It waited on the materials layer, which now exists --
    but the control it was waiting for did *not* land here: `coating` and
    `roughness` are `Problem` fields and belong in the geometry dock beside
    the period, because one `Problem` feeds every solver tab. What this tab
    owns is the material *library*, which duplicates no state."""

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

    def test_it_lists_the_installed_materials(self, win):
        """Generated from `materials.available()`, not written out -- a
        hand-maintained list is right on the day it is written."""
        text = win._setup_tab._browser.toPlainText()
        assert "Au" in text
        assert "nm" in text

    def test_and_says_where_they_came_from(self, win):
        text = win._setup_tab._browser.toPlainText()
        assert "Henke" in text
        assert "CXRO" in text

    def test_it_explains_what_naming_a_coating_does(self, win):
        """The distinction the whole milestone turns on, stated where the
        dropdown's one-word entries cannot."""
        text = win._setup_tab._browser.toPlainText()
        assert "absolute" in text
        assert "relative" in text

    def test_and_that_a_table_has_edges(self, win):
        """Out-of-range raises, which a user meets the first time they scan
        past 6.2 nm with gold selected."""
        assert "extrapolat" in win._setup_tab._browser.toPlainText()

    def test_it_still_offers_no_input_at_all(self, win):
        """The guard survives the tab gaining a purpose, and now protects
        something subtler: a coating control *here* would be a second place to
        set a `Problem` field the geometry dock already owns. Two sources of
        truth for one value is worse than the inert field this guard was
        originally written against.
        """
        from PySide6.QtWidgets import QComboBox, QLineEdit

        assert win._setup_tab.findChildren(QLineEdit) == []
        assert win._setup_tab.findChildren(QComboBox) == []

    def test_and_the_dock_is_where_the_control_actually_is(self, win):
        """Non-vacuity for the guard: the control exists, just not here."""
        assert "coating" in win.geometry._fields

    @pytest.fixture
    def library(self, tmp_path, monkeypatch):
        from gratinglab import materials

        directory = tmp_path / "library"
        monkeypatch.setenv(materials.LIBRARY_ENV, str(directory))
        return directory

    @staticmethod
    def _answer(monkeypatch, tab, path, name, replace=False):
        monkeypatch.setattr(tab, "_choose_file", lambda: str(path))
        monkeypatch.setattr(tab, "_choose_name", lambda default: name)
        monkeypatch.setattr(tab, "_confirm_replace", lambda _name: replace)

    @staticmethod
    def _coatings(win):
        combo = win.geometry._fields["coating"]
        return [combo.itemData(i) for i in range(combo.count())]

    def test_an_added_material_is_a_coating_without_a_restart(
        self, win, library, monkeypatch
    ):
        from gratinglab import materials

        tab = win._setup_tab
        self._answer(monkeypatch, tab, materials.data_dir() / "Au.txt", "Pt")
        assert "Pt" not in self._coatings(win)
        mentions = tab._browser.toPlainText().count("Pt")
        tab.add_material()
        assert (library / "Pt.txt").is_file()
        assert "Pt" in self._coatings(win)
        # The page re-rendered with a row for it; "Pt" was already the
        # example name in the prose, so count rather than test presence.
        assert tab._browser.toPlainText().count("Pt") > mentions
        assert "Added Pt" in tab._status.text()

    def test_an_unreadable_table_is_reported_and_not_added(
        self, win, library, tmp_path, monkeypatch
    ):
        broken = tmp_path / "Pt.txt"
        broken.write_text("not a CXRO export\n", encoding="utf-8")
        self._answer(monkeypatch, win._setup_tab, broken, "Pt")
        win._setup_tab.add_material()
        assert "was not added" in win._setup_tab._status.text()
        assert "Pt" not in self._coatings(win)

    def test_a_declined_replace_leaves_the_table_alone(
        self, win, library, monkeypatch
    ):
        from gratinglab import materials

        tab = win._setup_tab
        self._answer(monkeypatch, tab, materials.data_dir() / "Au.txt", "Pt")
        tab.add_material()
        before = (library / "Pt.txt").stat().st_mtime_ns
        tab.add_material()
        assert "was not added" in tab._status.text()
        assert (library / "Pt.txt").stat().st_mtime_ns == before

    def test_cancelling_either_dialog_does_nothing(self, win, library, monkeypatch):
        tab = win._setup_tab
        self._answer(monkeypatch, tab, "", "Pt")
        tab.add_material()
        monkeypatch.setattr(tab, "_choose_file", lambda: "Pt.txt")
        monkeypatch.setattr(tab, "_choose_name", lambda default: "")
        tab.add_material()
        assert not library.exists()
        assert tab._status.text() == ""

    def test_rescan_finds_a_file_added_by_hand(self, win, library):
        from gratinglab import materials

        library.mkdir()
        (library / "Ni.txt").write_bytes((materials.data_dir() / "Au.txt").read_bytes())
        assert "Ni" not in self._coatings(win)
        win._setup_tab.rescan()
        assert "Ni" in self._coatings(win)
        assert "Library read" in win._setup_tab._status.text()

    def test_a_removed_coating_falls_back_to_none_and_says_so(
        self, win, library, qtbot
    ):
        from gratinglab import materials

        library.mkdir()
        table = library / "Ni.txt"
        table.write_bytes((materials.data_dir() / "Au.txt").read_bytes())
        win._setup_tab.rescan()
        combo = win.geometry._fields["coating"]
        combo.setCurrentIndex(combo.findData("Ni"))
        table.unlink()
        with qtbot.waitSignal(win.geometry.changed, timeout=1000):
            win._setup_tab.rescan()
        assert combo.currentData() == ""

    def test_a_kept_coating_stays_selected(self, win, library):
        combo = win.geometry._fields["coating"]
        combo.setCurrentIndex(combo.findData("Au"))
        win._setup_tab.rescan()
        assert combo.currentData() == "Au"

    def test_open_library_folder_creates_it(self, win, library, monkeypatch):
        from gratinglab.gui.qt import setup_tab

        opened = []
        monkeypatch.setattr(
            setup_tab.QDesktopServices, "openUrl", lambda url: opened.append(url) or True
        )
        win._setup_tab.open_library_folder()
        assert library.is_dir()
        assert opened and opened[0].toLocalFile().endswith("library")

    def test_web_links_open_outside_and_nothing_else_does(self, win, monkeypatch):
        from PySide6.QtCore import QUrl

        from gratinglab.gui.qt import setup_tab

        opened = []
        monkeypatch.setattr(
            setup_tab.QDesktopServices, "openUrl", lambda url: opened.append(url) or True
        )
        win._setup_tab._follow(QUrl("https://henke.lbl.gov/optical_constants/"))
        win._setup_tab._follow(QUrl("file:///etc/passwd"))
        assert [u.toString() for u in opened] == ["https://henke.lbl.gov/optical_constants/"]

    def test_provenance_cites_what_is_installed(self, win):
        """Only the bundled CXRO gold is installed in the suite, so the page
        cites Henke and not the refractiveindex.info or PCGrate notes."""
        text = win._setup_tab._browser.toPlainText()
        assert "Atomic Data and Nuclear Data Tables" in text
        assert "Cite that measurement" not in text
        assert "carry no citation" not in text
        assert "roughness" in win.geometry._fields

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

        This used to be an AST scan of module-*level* imports. It is now a
        subprocess import per module, which is strictly stronger: it catches a
        *transitive* toolkit import, where a pure module grows an innocent
        import of something that itself pulls in Qt. A source scan cannot see
        that, and it is the version of this bug that would actually get past
        review.

        The deferred import in `app.py`'s `_require_qt()` was the reason the
        scan had to be shallow -- a recursive walk would flag it as if it were
        the unconditional kind. A subprocess import has no such problem: the
        function never runs at import time, so Qt never appears in `sys.modules`
        and the sanctioned deferral is simply invisible here.

        Modules are discovered rather than listed, so a new pure module is
        covered the day it lands.
        """
        assert_qt_free(
            modules_in("gratinglab.gui", exclude=("gratinglab.gui.qt",)),
            reason="only gratinglab/gui/qt may import a toolkit.",
        )

    def test_the_walk_finds_something(self):
        """Guard against discovery silently emptying and passing vacuously."""
        found = modules_in("gratinglab.gui", exclude=("gratinglab.gui.qt",))
        assert "gratinglab.gui.state" in found


class TestEntryPoint:
    """`app.py`'s own two jobs: explain a missing toolkit, and stay safe to
    import without one.

    Ported from the Tk suite this replaces -- PySide6 is pip-installable
    (unlike Tk), so the *message* changed, but the failure mode is still real:
    an install of `[dev]` without `[gui]` hits exactly this.
    """

    @pytest.mark.parametrize("module", ["app", "launcher"])
    def test_gives_an_explanation_not_a_traceback(self, monkeypatch, module):
        """The check is endstation's `require_extra`, which asks `find_spec`
        rather than importing, so that is what a missing PySide6 is faked
        through."""
        import importlib

        import endstation.optional

        entry = importlib.import_module(f"gratinglab.gui.{module}")
        real_find_spec = endstation.optional.find_spec

        def absent(name, *args, **kwargs):
            if name == "PySide6" or name.startswith("PySide6."):
                return None
            return real_find_spec(name, *args, **kwargs)

        monkeypatch.setattr(endstation.optional, "find_spec", absent)

        with pytest.raises(SystemExit) as excinfo:
            entry._require_qt()

        message = str(excinfo.value)
        assert "pip install -e '.[gui]'" in message
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


class TestOpeningASavedBoundary:
    """A gratinglab .json brings its measured period with it."""

    def _saved(self, tmp_path, period=314.97):
        from types import SimpleNamespace

        import numpy as np

        from gratinglab.io.boundary import write_boundary_json

        t = np.linspace(0.0, 1.0, 41)
        y = np.where(t < 0.7, t / 0.7, (1 - t) / 0.3) * 0.31
        y[0] = y[-1] = 0.0
        profile = SimpleNamespace(x_norm=t, y_norm=y, period_nm=period, n_grooves=6,
                                  n_used=5, n_edge_rejected=1, metrics={})
        return write_boundary_json(tmp_path / "measured.json", profile)

    def test_the_period_and_profile_both_arrive(self, qtbot, win, tmp_path):
        from gratinglab.gui.state import build

        path = self._saved(tmp_path)
        with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
            win.open_boundary(path)
        form = win.geometry.read_form()
        assert form.profile_kind == "From file"
        assert float(form.period) == pytest.approx(314.97)
        parsed = build(form)
        assert parsed.problem.period == pytest.approx(314.97)
        assert parsed.problem.profile.depth == pytest.approx(0.31)

    def test_the_label_says_where_the_period_came_from(self, qtbot, win, tmp_path):
        with qtbot.waitSignal(win.solved, timeout=SOLVE_TIMEOUT_MS):
            win.open_boundary(self._saved(tmp_path))
        assert "from the file" in win.geometry._path_label.text()


class TestGarbageCollectionStaysOnTheGuiThread:
    """See gratinglab.gui.qt.gc_guard: a collection on a worker thread can
    destroy Qt objects there, which segfaulted the macOS CI leg."""

    def test_solves_never_collect_on_the_worker(self, qtbot, win):
        import gc
        import threading

        gui_thread = threading.get_ident()
        seen = []

        def record(phase, info):
            if phase == "start":
                seen.append(threading.get_ident())

        gc.callbacks.append(record)
        try:
            for _ in range(3):
                junk = [[] for _ in range(5000)]
                for a, b in zip(junk, junk[1:]):
                    a.append(b)
                    b.append(a)
                del junk
                resolve(qtbot, win)
            qtbot.wait(1200)
        finally:
            gc.callbacks.remove(record)
        assert seen, "no collection ran at all, so the test proves nothing"
        assert set(seen) == {gui_thread}

    def test_automatic_collection_is_off_while_a_window_exists(self, win):
        import gc

        from gratinglab.gui.qt import gc_guard

        assert gc_guard.is_active()
        assert not gc.isenabled()

    def test_the_last_window_going_restores_it(self, qtbot):
        import gc

        from PySide6.QtWidgets import QApplication

        from gratinglab.gui.qt import gc_guard
        from gratinglab.gui.qt.main_window import MainWindow

        from PySide6.QtCore import QCoreApplication, QEvent

        # An earlier test's window may still be queued for deletion; let it
        # go first, or it releases the guard mid-test and the count races.
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        QCoreApplication.processEvents()
        before = gc_guard._holders
        window = MainWindow()
        with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
            pass
        assert gc_guard._holders == before + 1
        # close() first: it joins the solve thread, and deleting a window
        # whose QThread is still running aborts the process.
        window.close()
        window.deleteLater()
        qtbot.waitUntil(lambda: gc_guard._holders == before, timeout=5000)
        QApplication.processEvents()
        assert gc.isenabled() == (before == 0)


def test_exit_restores_collection_even_with_a_window_still_open(qtbot, win):
    """A closed-but-undeleted window never releases the guard; the atexit hook
    must still hand collection back, or Python's final pass frees Qt garbage
    after Qt is gone (the py3.11 bool_dealloc abort)."""
    import gc

    from gratinglab.gui.qt import gc_guard

    assert gc_guard.is_active() and not gc.isenabled()
    try:
        gc_guard._restore_at_exit()
        assert gc.isenabled()
    finally:
        # Put back what a live window expects, or the rest of the session
        # runs with collection on while windows still hold the guard.
        gc.disable()
