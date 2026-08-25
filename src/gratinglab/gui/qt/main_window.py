"""The window.

Deliberately thin, as its Tk predecessor was. Every number shown comes from
:func:`gratinglab.compare.sweep`; every validation from
:mod:`gratinglab.gui.state`; every line of a tab's provenance panel from
:mod:`gratinglab.gui.provenance`. Nothing here computes physics, so a bug in
this file misplaces a widget rather than producing a wrong answer -- and
`tests/test_gui_qt.py` pins that by comparing what a tab plotted against a
direct solver call.

Geometry inputs (:class:`~.geometry_panel.GeometryPanel`) live in a closable
dock, because they apply to every tab but should not tax the window when you
are not editing them -- Ctrl+G, and the tabs get the full width. The tabs are
the central widget, and nothing else is.

Per tab: whatever genuinely differs by method. A solver tab owns its options,
its result and its own provenance; the geometry tab owns every picture of the
grating itself, since a groove's shape and the directions light leaves in owe
nothing to a solver. Only one solve runs window-wide at a time (see
`solve`/`_set_running`) -- there is no second, slow solver yet to justify
concurrent per-tab workers.
"""

from __future__ import annotations

import threading

from PySide6.QtCore import Qt, QThread, QTimer, Signal
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (
    QDockWidget,
    QMainWindow,
    QMessageBox,
    QScrollArea,
    QTabWidget,
    QWidget,
)

from .. import provenance
from ..state import FormErrors, build, validate
from .geometry_panel import GeometryPanel
from .geometry_tab import GeometryTab
from .integral_tab import IntegralTab
from .scalar_tab import ScalarTab
from .setup_tab import SetupTab
from .worker import SolveWorker

__all__ = ["MainWindow"]

#: Stack for the solve thread, matching the platform's own main-thread default.
#: Not a tuning knob: a rigorous solver factors dense matrices through LAPACK,
#: and LAPACK expects the stack the main thread has. See `_start_worker`.
_WORKER_STACK_BYTES = 8 * 1024 * 1024

#: A solve faster than this never shows a progress bar. The scalar solver takes
#: about 70 ms, and flashing a bar on and off within one frame reads as a
#: glitch rather than as feedback.
_PROGRESS_DELAY_MS = 150

#: How long typing must pause before the geometry redraws. Long enough that a
#: burst of keystrokes costs one redraw rather than six, short enough that the
#: picture feels attached to the field you are editing.
_REDRAW_DEBOUNCE_MS = 120

#: How a registered solver name finds a tab class. A future RCWA backend adds
#: its own entry here and nothing else about this dispatch changes -- tabs are
#: generated from `available_solvers()`, not hardcoded, so a second solver
#: needs no change to `_build_tabs` itself.
#:
#: Every value here is expected to implement a small, duck-typed contract
#: rather than a formal base class (see `scalar_tab.py`'s module docstring for
#: why): a `name` class attribute matching its registry key; `solve_requested`
#: and `cancel_requested` signals; `build_options(problem, illumination,
#: wavelengths) -> dict` (raising `FormErrors`); `set_running(bool)`;
#: `show_progress()`; `show_progress_value(done, total)`;
#: `show_solving(wavelength_count)`; `show_result(scan, energy,
#: lambda_over_period)`; `show_cancelled()`; `show_error(message)`;
#: `show_field_errors(errors)`. Optionally a `convergence_requested` signal and
#: `show_convergence(scan, energy, lambda_over_period, report)` -- optional
#: because a backend with no `accuracy_knob` has nothing to sweep.
_TAB_FACTORIES = {"scalar": ScalarTab, "integral": IntegralTab}


class MainWindow(QMainWindow):
    """One geometry, a tab per modeling approach."""

    #: Emitted after a result is drawn. Tests wait on this rather than
    #: sleeping, which is what keeps the widget suite fast and deterministic.
    solved = Signal()

    #: How work reaches the worker. It must be a signal, not a method call:
    #: `moveToThread` only redirects *signal* delivery, so calling
    #: `worker.run(...)` directly would run the solve on the UI thread and
    #: freeze the window -- the exact thing this design exists to prevent.
    _requested = Signal(int, str, object, dict, object)
    #: The same payload, routed to the worker's sweep slot instead. A second
    #: signal rather than a mode flag: `moveToThread` dispatches per slot, and
    #: a flag would put the branch on the worker thread instead of here.
    _convergence_requested = Signal(int, str, object, dict, object)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("GratingLab")
        self.resize(1180, 820)

        self._parsed = None
        self._token = 0
        self._running = False
        self._active_name: str | None = None
        # Replaced per solve rather than reused, so a Cancel clicked between
        # the request and the worker picking it up cannot be lost to a clear().
        self._cancel = threading.Event()

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
        self._start_worker()
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

        A `QDockWidget` rather than a splitter pane because this is already a
        `QMainWindow`, so a title bar, a close button, drag-to-float and a
        restorable action all come for free -- and the panel is a permanent
        300 px tax on a 1180 px window otherwise.
        """
        self.geometry = GeometryPanel()
        self.geometry.solve_requested.connect(self._solve_active_tab)
        self.geometry.changed.connect(self._geometry_edited)

        scroller = QScrollArea()
        scroller.setWidget(self.geometry)
        scroller.setWidgetResizable(True)
        scroller.setMinimumWidth(280)

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
        from ..docs import display_title
        from ...solvers import available_solvers

        self.tabs: dict[str, QWidget] = {}
        self._tab_widget = QTabWidget()

        # The two non-solver tabs come first, so the bar has a single boundary
        # between "about the problem" and "about a method". Neither is in
        # `self.tabs`: neither implements the solve/cancel contract, and
        # `_solve_active_tab`'s `getattr(..., "name", None)` skips both.
        self._setup_tab = SetupTab()
        self._tab_widget.addTab(self._setup_tab, "Setup")

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
            self._tab_widget.addTab(tab, display_title(name))

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
        """Two menus, and only two.

        Help is the original and still the main one: this is not a
        general-purpose application with File/Edit concerns, and the single
        thing worth a menu for is explaining the math, which the window itself
        never does.

        View exists for exactly one reason and holds exactly one action. The
        geometry panel is closable now -- it costs 300 px of a 1180 px window
        otherwise -- and a control the user can close with no way to reopen it
        is a trap. Putting the way back under *Help* would be worse than a
        View menu: Help here means "explain the math", and diluting that has a
        real cost.

        `QDockWidget.toggleViewAction()` rather than a hand-rolled checkable
        action, because a hand-rolled one desynchronises the moment the user
        clicks the dock's own close button.
        """
        from ..docs import general_pages, theory_pages
        from ..richtext import menu_label

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

        # View first, so the bar reads View | Help.
        self._view_menu = self._menu_bar.addMenu("&View")
        toggle = self.geometry_dock.toggleViewAction()
        toggle.setText("&Geometry inputs")
        toggle.setShortcut(QKeySequence("Ctrl+G"))
        self._view_menu.addAction(toggle)

        menu = self._help_menu = self._menu_bar.addMenu("&Help")
        menu.addAction("About GratingLab", self.show_about)
        menu.addSeparator()

        # Foundational, solver-independent reference first -- the generalized
        # grating equation and the angle conventions -- then one page per
        # method.
        for page in general_pages():
            suffix = "" if page.available else " (not written yet)"
            menu.addAction(
                menu_label(page.title, suffix), lambda p=page: self.show_theory(p)
            )
        menu.addSeparator()
        for page in theory_pages():
            suffix = " Theory" + ("" if page.available else " (not written yet)")
            menu.addAction(
                menu_label(page.title, suffix), lambda p=page: self.show_theory(p)
            )

    def _start_worker(self) -> None:
        self._thread = QThread(self)
        # A QThread starts with a much smaller stack than the main thread's,
        # and LAPACK's dense complex solve overruns it: `np.linalg.solve` on a
        # few-hundred-square matrix segfaults the *process* here while running
        # fine on the main thread and on a plain `threading.Thread`. That made
        # the integral tab unusable in the real window -- both boundary
        # conditions -- while every headless test passed, because nothing but
        # the GUI ran a solve off the main thread. The scalar solver never
        # tripped it: it factors nothing.
        self._thread.setStackSize(_WORKER_STACK_BYTES)
        self._worker = SolveWorker()
        self._worker.moveToThread(self._thread)
        # Cross-thread, so both directions are queued automatically: requests
        # arrive on the worker's thread, results back on this one.
        self._requested.connect(self._worker.run)
        self._worker.finished.connect(self._on_solved)
        self._worker.failed.connect(self._on_failed)
        self._worker.cancelled.connect(self._on_cancelled)
        self._worker.progress.connect(self._on_progress)
        self._worker.converged.connect(self._on_converged)
        self._convergence_requested.connect(self._worker.run_convergence)
        self._thread.start()

    # -- actions ---------------------------------------------------------

    def _geometry_edited(self) -> None:
        """A field was edited. Schedule a redraw; solve nothing.

        Restarting the timer is the debounce: only a pause in typing reaches
        `_refresh_geometry_tab`.
        """
        self._redraw_timer.start()

    def _refresh_geometry_tab(self) -> None:
        """Redraw the geometry tab from the form as it stands.

        Deliberately not driven off `_on_solved`: the tab holds no solver
        output, so waiting for a solve to redraw a picture no solver
        contributed to would be a false dependency. It runs at construction,
        after each solve, and -- debounced -- on every edit.

        This is **not** the live re-solve ruled out of scope. Nothing here
        reaches a worker or a solver: `state.build` is microseconds of parsing
        and the redraw is matplotlib, which the tab defers anyway when it is
        hidden.
        """
        errors = validate(self.geometry.read_form())
        if errors:
            # A form mid-edit is not a form with an error. The tab keeps its
            # last good drawing and says, dimly, that it is one edit behind.
            self.geometry_tab.show_pending(errors)
            return
        self.geometry_tab.show_geometry(build(self.geometry.read_form()))

    def _solve_active_tab(self) -> None:
        """Enter in a geometry field solves whichever tab is showing.

        The geometry panel does not know solver names; it only knows
        something changed. A tab with no `name` (Setup, once it exists) is
        silently not solvable -- there is nothing to solve there.
        """
        name = getattr(self._tab_widget.currentWidget(), "name", None)
        if name is not None:
            self.solve(name)

    def solve(self, name: str, *, sweep: bool = False) -> None:
        """Validate here, solve on the worker thread.

        Geometry and the tab's own options are validated as two separate
        steps, but both are caught before anything reaches the worker: a
        rejected form is never sent to the worker, whichever kind of error it
        was.

        ``sweep`` runs the convergence harness instead of a single solve. It
        shares every bit of the machinery below -- one token, one cancel
        event, one "only one at a time window-wide" rule -- because it is the
        same kind of thing, only longer. Only the destination slot differs.
        """
        if self._running:
            return
        tab = self.tabs[name]

        try:
            geometry = build(self.geometry.read_form())
            options = tab.build_options(
                geometry.problem, geometry.illumination, geometry.wavelengths
            )
        except FormErrors as exc:
            tab.show_field_errors(exc.errors)
            return

        self._active_name = name
        self._parsed = geometry
        self._token += 1
        self._cancel = threading.Event()
        self._set_running(True)
        tab.show_solving(len(geometry.wavelengths))
        # `tab` as the context object, not a bare singleShot: Qt then cancels
        # the timer if the tab is destroyed first. Without it, closing the
        # window inside the 150 ms delay leaves the timer holding a deleted
        # QProgressBar, and an exception escaping a slot under PySide6 can
        # abort the process rather than surfacing. Found by macOS CI, where the
        # timing landed the stale timer inside the *next* test's setup.
        QTimer.singleShot(
            _PROGRESS_DELAY_MS, tab, lambda: self._show_progress_if_still_running(tab)
        )
        signal = self._convergence_requested if sweep else self._requested
        signal.emit(self._token, name, geometry, options, self._cancel)

    def cancel(self) -> None:
        """Stop the running solve.

        Setting the event is what actually stops it: the solver checks it at
        every wavelength and raises out of its own loop (see `worker.py`).

        Bumping the token as well is not redundant. The worker may already
        have finished when Cancel is clicked, in which case a result is
        already in flight and the event will never be read -- the token is
        what makes that one stale on arrival. Two different situations, two
        mechanisms.

        There is only ever one thing that could be running, so it does not
        matter which tab's Cancel button was clicked.
        """
        if not self._running:
            return
        self._cancel.set()
        self._token += 1
        self._set_running(False)
        if self._active_name is not None:
            self.tabs[self._active_name].show_cancelled()

    def _on_solved(self, token: int, method: str, result) -> None:
        if token != self._token:
            return  # cancelled or superseded; the window has moved on
        self._set_running(False)
        # No `profile_panel.draw` here any more: the groove's shape owes
        # nothing to a solver, so redrawing it on solve completion was a false
        # dependency. It lives in the geometry tab now and redraws with it.
        self.geometry_tab.show_geometry(self._parsed)
        self.tabs[method].show_result(
            result.scan, result.energy, self._parsed.lambda_over_period
        )
        self.solved.emit()

    def _on_progress(self, token: int, done: int, total: int) -> None:
        """Advance the originating tab's bar, and nobody else's.

        A stale token means a cancelled solve is still reporting on its way
        out of its own loop -- a few wavelengths at most, but enough to walk a
        bar backwards over a result the window has already moved on from.
        """
        if token != self._token or self._active_name is None:
            return
        self.tabs[self._active_name].show_progress_value(done, total)

    def _on_cancelled(self, token: int, method: str) -> None:
        """The solve really stopped.

        Nothing to show: `cancel()` already updated the panel and re-enabled
        the buttons, and the previous result is still on screen. This exists
        so the worker has somewhere to report to other than `failed` -- a stop
        the user asked for is not an error.
        """

    def _on_converged(self, token: int, method: str, result, report) -> None:
        if token != self._token:
            return
        self._set_running(False)
        self.geometry_tab.show_geometry(self._parsed)
        self.tabs[method].show_convergence(
            result.scan, result.energy, self._parsed.lambda_over_period, report
        )
        self.solved.emit()

    def _on_failed(self, token: int, method: str, message: str) -> None:
        if token != self._token:
            return
        self._set_running(False)
        self.tabs[method].show_error(message)

    def _set_running(self, running: bool) -> None:
        """Every tab's Solve button disables together, not just the one that
        was pressed -- only one solve can be in flight window-wide."""
        self._running = running
        for tab in self.tabs.values():
            tab.set_running(running)

    def _show_progress_if_still_running(self, tab) -> None:
        """Progress is shown only on the tab that started the solve.

        Showing it on an idle tab would claim that tab is computing something
        it isn't -- the same shape of dishonesty a provenance warning must
        never carry either.
        """
        if self._running:
            tab.show_progress()

    def show_about(self) -> None:
        from ... import __version__

        QMessageBox.about(self, "About GratingLab", provenance.about_text(__version__))

    def show_theory(self, page) -> None:
        """Open the read-only viewer for one theory or reference page."""
        from .theory_viewer import TheoryViewer

        viewer = TheoryViewer(page, parent=self)
        viewer.show()
        return viewer

    # -- shutdown --------------------------------------------------------

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        """Join the worker thread.

        Without this Qt reports "QThread: Destroyed while thread is still
        running" and aborts intermittently -- most visibly in the offscreen
        CI runs, where windows are created and destroyed back to back.
        """
        self._redraw_timer.stop()  # a redraw firing into a closing window
        self._thread.quit()
        self._thread.wait(5000)
        super().closeEvent(event)
