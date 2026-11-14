"""The main window.

The window is thin, as its Tk predecessor was. Every number shown comes from
:func:`gratinglab.compare.sweep`, every validation from
:mod:`gratinglab.gui.state`, and every line of a tab's provenance panel from
:mod:`gratinglab.gui.provenance`. Nothing here computes physics, so a bug in
this file can misplace a widget but cannot produce a wrong answer, and
`tests/test_gui_qt.py` checks this by comparing what a tab plotted against a
direct solver call.

The geometry inputs (:class:`~.geometry_panel.GeometryPanel`) live in a
closable dock, because they apply to every tab but should not take up space
while they are not being edited; Ctrl+G hides them and gives the tabs the full
width. The tabs alone form the central widget.

Each tab holds whatever differs by method. A solver tab owns its options, its
result and its own provenance. The geometry tab owns every picture of the
grating itself, since the shape of a groove and the directions in which light
leaves it do not depend on a solver. The Compare tab, last, owns no solver:
it runs several in turn, each with its own tab's options, and shows them
together (`compare`). Only one solve runs in the whole window
at a time (see `solve`/`_set_running`), because there is not yet a second,
slow solver to justify concurrent workers per tab.
"""

from __future__ import annotations

from functools import partial

from endstation.qt.jobs import JobRunner
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (
    QDockWidget,
    QMainWindow,
    QScrollArea,
    QTabWidget,
    QWidget,
)

from ..failures import describe_failure
from ..state import FieldError, FormErrors, build, validate
from . import recent, theming
from .cmethod_tab import CMethodTab
from .compare_tab import CompareTab
from .geometry_panel import GeometryPanel
from .geometry_tab import GeometryTab
from .integral_tab import IntegralTab
from .rcwa_tab import RcwaTab
from .scalar_tab import ScalarTab
from .setup_tab import SetupTab
from .sizing import fit_width_to_contents
from . import worker

__all__ = ["MainWindow"]

#: Stack for the solve thread, matching the platform's own main-thread default.
#: Not a tuning knob: a rigorous solver factors dense matrices through LAPACK,
#: and LAPACK expects the stack the main thread has. See `__init__`.
_WORKER_STACK_BYTES = 8 * 1024 * 1024

#: How long typing must pause before the geometry redraws. Long enough that a
#: burst of keystrokes costs one redraw rather than six, short enough that the
#: picture feels attached to the field you are editing.
_REDRAW_DEBOUNCE_MS = 120

#: How a registered solver name finds a tab class. A new backend adds its own
#: entry here and nothing else about this dispatch changes (RCWA did, as the
#: third) -- tabs are
#: generated from `available_solvers()`, not hardcoded, so a second solver
#: needs no change to `_build_tabs` itself.
#:
#: Every value here is expected to implement a small, duck-typed contract
#: rather than a formal base class (see `scalar_tab.py`'s module docstring for
#: why): a `name` class attribute matching its registry key; `solve_requested`
#: and `cancel_requested` signals; `build_options(problem, illumination,
#: wavelengths) -> dict` (raising `FormErrors`); `set_running(bool)`;
#: `show_progress_value(done, total)`; `show_solving(wavelength_count)`; `show_result(scan, energy,
#: lambda_over_period)`; `show_cancelled()`; `show_error(message)`;
#: `show_field_errors(errors)`. Optionally a `convergence_requested` signal and
#: `show_convergence(scan, energy, lambda_over_period, report)` -- optional
#: because a backend with no `accuracy_knob` has nothing to sweep.
_TAB_FACTORIES = {"scalar": ScalarTab, "integral": IntegralTab, "rcwa": RcwaTab}
_TAB_FACTORIES["cmethod"] = CMethodTab


def _refuse_layers_unless_handled(name: str, problem) -> None:
    """Put a solver's refusal of layers on the Layers field, not the worker.

    ``Capabilities.check`` would refuse the same problem inside the solve, but
    as an exception from the worker thread naming no field. Asked here, it
    lands beside every other form error, on the field that caused it.
    """
    from ...solvers import get_solver
    from ..docs import short_title

    if problem.layers and not get_solver(name).capabilities.handles_layers:
        raise FormErrors((FieldError(
            "layers",
            f"the {short_title(name)} solver cannot model conformal layers; "
            "the RCWA tab can. Clear Layers to solve the bare grating here.",
        ),))


class MainWindow(QMainWindow):
    """One geometry, a tab per modeling approach."""

    #: Emitted after a result is drawn. Tests wait on this rather than
    #: sleeping, which is what keeps the widget suite fast and deterministic.
    solved = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("gratinglab")
        self.resize(1180, 820)
        # Before any widget exists, so every figure is born in the theme.
        from PySide6.QtWidgets import QApplication

        self._themes = theming.install(QApplication.instance())
        # Solves run on a worker thread; keep the cycle collector off it.
        from endstation.qt.gc_guard import guard

        guard(self)

        self._parsed = None
        self._running = False
        self._active_name: str | None = None
        # One worker thread, shared by every tab, since only one solve is ever
        # in flight in the window. Its stack is the main thread's: a QThread's
        # default is far smaller, and LAPACK's dense complex solve overran it,
        # so `np.linalg.solve` on a few-hundred-square matrix segfaulted the
        # process and the integral tab returned nothing in the real window
        # while every headless test passed. The runner drops a cancelled job's
        # late result, and words a refusal as a refusal rather than a bug.
        self._runner = JobRunner(
            self, stack_size=_WORKER_STACK_BYTES, describe=describe_failure
        )

        # Coalesces a burst of keystrokes into one redraw. Restartable and
        # single-shot, parented to the window so it dies with it -- a bare
        # `QTimer.singleShot` cannot be restarted, so every keystroke would
        # queue its own redraw and typing "315.15" would draw six times.
        self._redraw_timer = QTimer(self)
        self._redraw_timer.setSingleShot(True)
        self._redraw_timer.setInterval(_REDRAW_DEBOUNCE_MS)
        self._redraw_timer.timeout.connect(self._refresh_geometry_tab)

        # Layout first: `_build_menu` puts the dock's own toggleViewAction in
        # the View menu, so the dock has to exist by then.
        self._build_layout()
        self._build_menu()
        recent.restore_layout(self, "efficiency")
        # Geometry first: it needs no solver, so it can be right before any
        # result exists rather than showing an empty panel until one does.
        self._refresh_geometry_tab()
        self._solve_active_tab()

    # -- construction ----------------------------------------------------

    def _build_layout(self) -> None:
        self._build_dock()
        # The tabs *are* the central widget. With the dock closable, that is
        # what lets Ctrl+G hand the whole window width to whichever tab is
        # showing.
        self.setCentralWidget(self._build_tabs())

    def _build_dock(self) -> None:
        """Geometry inputs, in a panel the user is allowed to close.

        A `QDockWidget` is used instead of a splitter pane because this is
        already a `QMainWindow`, so a title bar, a close button, drag-to-float
        and a restorable action are provided without extra code. Without the
        dock, the panel would permanently occupy 300 px of a 1180 px window.
        """
        self.geometry = GeometryPanel()
        self.geometry.solve_requested.connect(self._solve_active_tab)
        self.geometry.changed.connect(self._geometry_edited)

        scroller = QScrollArea()
        scroller.setWidget(self.geometry)
        scroller.setWidgetResizable(True)
        # Derived, with the old literal kept as a floor so a platform whose
        # fonts fit in 280 px renders exactly as before. The literal alone
        # clipped the dock on Windows, whose form needs 388 px; see sizing.py.
        scroller.setMinimumWidth(max(280, fit_width_to_contents(scroller)))

        dock = self.geometry_dock = QDockWidget("Grating geometry", self)
        dock.setObjectName("geometry_dock")  # so a future saveState() works
        dock.setWidget(scroller)
        dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )
        # Left or right only. Docked top or bottom, a tall column of form rows
        # becomes a 1180-px-wide strip of three fields -- allowed by Qt,
        # useless in practice.
        dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
        # resize() does nothing to a dock; this is the supported way.
        self.resizeDocks([dock], [300], Qt.Orientation.Horizontal)

    def _build_tabs(self) -> QWidget:
        from ..docs import display_title, short_title
        from ...solvers import available_solvers

        self.tabs: dict[str, QWidget] = {}
        self._tab_widget = QTabWidget()

        # The two non-solver tabs come first, so the bar has a single boundary
        # between "about the problem" and "about a method". Neither is in
        # `self.tabs`: neither implements the solve/cancel contract, and
        # `_solve_active_tab`'s `getattr(..., "name", None)` skips both.
        self._setup_tab = SetupTab()
        self._tab_widget.addTab(self._setup_tab, "Setup")
        # The dock is built first, so its dropdown can follow the library.
        self._setup_tab.library_changed.connect(self.geometry.refresh_coatings)

        self.geometry_tab = GeometryTab()
        self._tab_widget.addTab(self.geometry_tab, "Grating Geometry")

        for name in available_solvers():
            factory = _TAB_FACTORIES.get(name)
            if factory is None:
                continue  # registered, but no tab written for it yet
            tab = factory()
            tab.solve_requested.connect(lambda _checked=False, n=name: self.solve(n))
            tab.cancel_requested.connect(self.cancel)
            converge = getattr(tab, "convergence_requested", None)
            if converge is not None:
                # Optional: a future tab for a backend with no accuracy_knob
                # has nothing to sweep and should not carry the button.
                converge.connect(lambda _checked=False, n=name: self.solve(n, sweep=True))
            self.tabs[name] = tab
            index = self._tab_widget.addTab(tab, short_title(name))
            self._tab_widget.setTabToolTip(index, display_title(name))

        # Last, after every solver it draws on. Not in `self.tabs`: it has no
        # solver of its own, and everything that loops over `self.tabs`
        # (the solve contract, `_set_running`'s Solve buttons) is about one.
        self.compare_tab = CompareTab(tuple(self.tabs))
        self.compare_tab.compare_requested.connect(self.compare)
        self.compare_tab.cancel_requested.connect(self.cancel)
        index = self._tab_widget.addTab(self.compare_tab, "Compare")
        self._tab_widget.setTabToolTip(
            index, "Solve the geometry with several solvers and show them together"
        )

        # Setup is leftmost so it is easy to find, but a solver tab -- not
        # the explanatory stub -- is what should be showing, and solving,
        # when the window opens: a user wants a result, not a "nothing here
        # yet" page. `_solve_active_tab` (used for both the construction-time
        # solve and Enter-in-a-geometry-field) reads `currentWidget()`, so
        # this is also what makes the window open already solving.
        #
        # Scalar specifically, not the alphabetically first: the opening
        # solve should be the ~70 ms survey, and the default geometry is a
        # 200-wavelength scan that the integral solver takes *minutes* over.
        # A rigorous method is one deliberate click away, behind a progress
        # bar and a working Cancel, which is where minutes belong.
        if self.tabs:
            first = self.tabs.get("scalar", next(iter(self.tabs.values())))
            self._tab_widget.setCurrentWidget(first)

        return self._tab_widget

    def _build_menu(self) -> None:
        """Build the menus: File, View and Help.

        File holds the two actions worth a shortcut: Ctrl+O loads a boundary
        profile and Ctrl+Return solves whichever tab is showing (or runs the
        comparison on the Compare tab), the same as Enter in a geometry field.

        Help is the original menu and remains the main one. This is not a
        general-purpose application with File/Edit concerns, and the one thing
        that merits a menu is the explanation of the mathematics, which the
        window itself never gives. It is built by
        :func:`gratinglab.gui.qt.help_menu.install_help_menu`, which the AFM
        window and the chooser call too, so the pages are reachable from
        wherever the user happens to be.

        View began with a single purpose and a single action, and now holds
        a second: the Theme submenu, which is about how the window looks and
        so belongs nowhere else. The geometry panel can now be closed (otherwise it would occupy 300 px of
        a 1180 px window), and a control that the user can close with no way
        to reopen it would leave the user stuck. Placing the reopening action
        under Help would be worse than adding a View menu, since Help here
        means "explain the math", and diluting that meaning has a real cost.

        `QDockWidget.toggleViewAction()` is used instead of a hand-written
        checkable action, because a hand-written one falls out of
        synchronisation as soon as the user clicks the dock's own close
        button.
        """
        from .help_menu import install_help_menu

        # Both the bar and the menus are held on the instance rather than
        # merely added. Under PySide6 6.11 a QMenu reached back through
        # `menuBar().actions()[i].menu()` is not reliably usable -- touching
        # it raises "Internal C++ object already deleted" even while the menu
        # the window built is perfectly alive. Measured, not guessed: holding
        # the QMenuBar makes no difference, holding the QMenu does.
        #
        # So the window owns what it builds, and nothing (tests included)
        # should navigate back down from the menu bar to get at it.
        self._menu_bar = self.menuBar()

        # File first, so the bar reads File | View | Help.
        self._file_menu = self._menu_bar.addMenu("&File")
        load = self._file_menu.addAction("&Open boundary profile…")
        load.setShortcut(QKeySequence.StandardKey.Open)
        load.triggered.connect(self.geometry.load_profile)
        solve = self._file_menu.addAction("&Solve this tab")
        solve.setShortcut(QKeySequence("Ctrl+Return"))
        solve.triggered.connect(self._solve_active_tab)
        self._file_actions = {"open": load, "solve": solve}

        self._view_menu = self._menu_bar.addMenu("&View")
        toggle = self.geometry_dock.toggleViewAction()
        toggle.setText("&Geometry inputs")
        toggle.setShortcut(QKeySequence("Ctrl+G"))
        self._view_menu.addAction(toggle)
        self._view_menu.addSeparator()

        # One exclusive group: System follows the OS, the rest are fixed.
        # Shared with the AFM window, so both offer it and both stay checked
        # in step whichever one the choice was made in.
        self._theme_menu = theming.ThemeMenu(self)
        self._view_menu.addMenu(self._theme_menu)
        self._theme_actions = self._theme_menu._group

        self._help_menu = install_help_menu(self, self._menu_bar)

    def set_theme(self, choice: str) -> None:
        """Switch theme by key, or ``"system"``; the menu calls this."""
        self._themes.set_choice(choice)

    # -- actions ---------------------------------------------------------

    def open_boundary(self, path) -> None:
        """Load a boundary profile file as the grating, as if chosen by hand.

        The metrology window's "Send to efficiency window" calls this with a
        gratinglab ``.json``, so the measured period arrives with the shape.
        """
        self.geometry.load_boundary_file(path)
        if not self.geometry_dock.isVisible():
            self.geometry_dock.show()

    def _geometry_edited(self) -> None:
        """Schedule a redraw after a field was edited, without solving.

        Restarting the timer implements the debounce, so that only a pause in
        typing reaches `_refresh_geometry_tab`.
        """
        self._redraw_timer.start()

    def _refresh_geometry_tab(self) -> None:
        """Redraw the geometry tab from the form as it stands.

        It is not driven from `_on_solved`. The tab holds no solver output, so
        waiting for a solve before redrawing a picture to which no solver
        contributed would introduce a false dependency. It runs at
        construction, after each solve, and (debounced) on every edit.

        This is not the live re-solve that was ruled out of scope. Nothing
        here reaches a worker or a solver: `state.build` takes microseconds of
        parsing, and the redraw is done by matplotlib, which the tab defers in
        any case while it is hidden.
        """
        errors = validate(self.geometry.read_form())
        if errors:
            # A form mid-edit is not a form with an error. The tab keeps its
            # last good drawing and says, dimly, that it is one edit behind.
            self.geometry_tab.show_pending(errors)
            return
        self.geometry_tab.show_geometry(build(self.geometry.read_form()))

    def _solve_active_tab(self) -> None:
        """Solve whichever tab is showing, in response to Enter in a geometry field.

        The geometry panel does not know solver names; it only knows that
        something changed. A tab with no `name` (Setup, once it exists) is
        silently skipped, because there is nothing to solve there.
        """
        current = self._tab_widget.currentWidget()
        if current is self.compare_tab:
            self.compare()
            return
        name = getattr(current, "name", None)
        if name is not None:
            self.solve(name)

    def solve(self, name: str, *, sweep: bool = False) -> None:
        """Validate here and solve on the worker thread.

        The geometry and the tab's own options are validated in two separate
        steps, but errors from both are caught before anything reaches the
        worker, so a rejected form is never sent to it, whatever kind of error
        it contained.

        ``sweep`` runs the convergence harness instead of a single solve. It
        shares all of the machinery below (the one runner, Cancel, and the
        rule that only one solve runs in the window at a time), because a
        sweep is the same kind of operation, only longer. Only the function
        the runner calls differs.
        """
        if self._running:
            return
        tab = self.tabs[name]

        try:
            geometry = build(self.geometry.read_form())
            _refuse_layers_unless_handled(name, geometry.problem)
            options = tab.build_options(
                geometry.problem, geometry.illumination, geometry.wavelengths
            )
        except FormErrors as exc:
            tab.show_field_errors(exc.errors)
            return

        self._active_name = name
        self._parsed = geometry
        self._set_running(True)
        # The tab's bar appears only if the solve is still running after its
        # delay, and only on this tab; `_set_running(False)` ends every tab's.
        tab.show_solving(len(geometry.wavelengths))
        if sweep:
            self._start(worker.converge, name, geometry, options,
                        on_done=partial(self._on_converged, name))
        else:
            self._start(worker.solve, name, geometry, options,
                        on_done=partial(self._on_solved, name))

    def compare(self) -> None:
        """Solve the geometry with every solver ticked on the Compare tab.

        Validated here, as a single solve is: a geometry error stops the
        whole comparison and lands on the Compare tab's panel. A solver whose
        *own* tab rejects its settings, or that cannot model the layers, is
        different. That is one solver declining, so it is named among the
        refusals with the reason, its tab shows the error on the field, and
        the others still run. A refusal only the solver itself can detect
        (the C-method on a vertical facet) is caught the same way on the
        worker; see `worker.compare`.

        It is one operation to the rest of the window: one job, one Cancel,
        every Solve button disabled until it ends.
        """
        if self._running:
            return
        tab = self.compare_tab
        names = tab.selected()
        if not names:
            return
        try:
            geometry = build(self.geometry.read_form())
        except FormErrors as exc:
            tab.show_field_errors(exc.errors)
            return

        from ..docs import short_title

        options: dict[str, dict] = {}
        refused: dict[str, str] = {}
        for name in names:
            try:
                _refuse_layers_unless_handled(name, geometry.problem)
                options[name] = self.tabs[name].build_options(
                    geometry.problem, geometry.illumination, geometry.wavelengths
                )
            except FormErrors as exc:
                self.tabs[name].show_field_errors(exc.errors)
                refused[name] = "; ".join(
                    f"{e.field}: {e.message}" for e in exc.errors
                ) + f" (on the {short_title(name)} tab)"

        self._active_name = tab.name
        self._parsed = geometry
        self._set_running(True)
        tab.show_solving(tuple(options), len(geometry.wavelengths))
        self._start(worker.compare, geometry, names, options, refused,
                    on_done=self._on_compared)

    def _start(self, fn, *args, on_done) -> None:
        """Hand one piece of work to the runner, with the window's handlers.

        A failure is reported against the tab that asked, which is the
        active one. A solve started straight after a Cancel queues behind the
        cancelled one while it unwinds.
        """
        self._runner.start(
            fn, *args,
            on_done=on_done,
            on_fail=partial(self._on_failed, self._active_name),
            on_progress=self._on_progress,
            on_cancel=self._on_cancelled,
        )

    def _on_compared(self, outcome) -> None:
        self._set_running(False)
        self.geometry_tab.show_geometry(self._parsed)
        self.compare_tab.show_outcome(outcome)
        self.solved.emit()

    def _tab_for(self, name: str):
        """The tab a reply belongs to: a solver's, or the Compare tab."""
        if name == self.compare_tab.name:
            return self.compare_tab
        return self.tabs[name]

    def cancel(self) -> None:
        """Stop the running solve.

        The runner's progress callback raises at the solver's next
        wavelength, so the work itself stops (see `worker.py`). A result that
        was already on its way when Cancel landed is dropped by the runner
        rather than drawn.

        There is only ever one thing that could be running, so it does not
        matter which tab's Cancel button was clicked.
        """
        if not self._running:
            return
        self._runner.cancel()
        self._set_running(False)
        if self._active_name is not None:
            self._tab_for(self._active_name).show_cancelled()

    def _on_solved(self, method: str, result) -> None:
        self._set_running(False)
        # No `profile_panel.draw` here any more: the groove's shape owes
        # nothing to a solver, so redrawing it on solve completion was a false
        # dependency. It lives in the geometry tab now and redraws with it.
        self.geometry_tab.show_geometry(self._parsed)
        self.tabs[method].show_result(
            result.scan, result.energy, self._parsed.lambda_over_period
        )
        self.solved.emit()

    def _on_progress(self, done: int, total: int) -> None:
        """Advance the originating tab's bar, and nobody else's.

        A cancelled solve's last reports, from the wavelength it was in, never
        arrive here: the runner drops them, so they cannot walk the bar over a
        result the window has moved on from.
        """
        if self._active_name is None:
            return
        self._tab_for(self._active_name).show_progress_value(done, total)

    def _on_cancelled(self) -> None:
        """Acknowledge that the solve has stopped.

        There is nothing to show: `cancel()` has already updated the panel and
        re-enabled the buttons, and the previous result is still on screen.
        This handler exists so that a cancelled job reports somewhere other
        than `on_fail`, since a stop that the user requested is not an error.
        """

    def _on_converged(self, method: str, outcome) -> None:
        result, report = outcome
        self._set_running(False)
        self.geometry_tab.show_geometry(self._parsed)
        self.tabs[method].show_convergence(
            result.scan, result.energy, self._parsed.lambda_over_period, report
        )
        self.solved.emit()

    def _on_failed(self, method: str, message: str) -> None:
        self._set_running(False)
        self._tab_for(method).show_error(message)

    def _set_running(self, running: bool) -> None:
        """Enable or disable every tab's Solve button together.

        All of them change state, not only the one that was pressed, because
        only one solve can be in flight in the whole window."""
        self._running = running
        for tab in self.tabs.values():
            tab.set_running(running)
        self.compare_tab.set_running(running)

    def show_about(self) -> None:
        """The About box. The menu calls `help_menu.show_about` directly; this
        stays as the window's own way in, which is what a test drives."""
        from .help_menu import show_about

        show_about(self)

    def show_theory(self, page):
        """Open the read-only viewer for one theory or reference page."""
        from .help_menu import show_theory

        return show_theory(page, self)

    # -- shutdown --------------------------------------------------------

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        """Stop the solve in flight and wait for the worker to finish.

        A solve still running when the window goes would deliver into widgets
        that no longer exist, most visibly in the offscreen CI runs, where
        windows are created and destroyed in quick succession.
        """
        recent.save_layout(self, "efficiency")
        self._redraw_timer.stop()  # a redraw firing into a closing window
        self._runner.cancel()
        self._runner.wait(5000)
        super().closeEvent(event)
