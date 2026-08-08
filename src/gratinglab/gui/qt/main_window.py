"""The window.

Deliberately thin, as its Tk predecessor was. Every number shown comes from
:func:`gratinglab.compare.sweep`; every validation from
:mod:`gratinglab.gui.state`; every line of a tab's provenance panel from
:mod:`gratinglab.gui.provenance`. Nothing here computes physics, so a bug in
this file misplaces a widget rather than producing a wrong answer -- and
`tests/test_gui_qt.py` pins that by comparing what a tab plotted against a
direct solver call.

Shared, above the tabs: geometry (:class:`~.geometry_panel.GeometryPanel`) and
the groove-profile plot (:class:`~.profile_plot_panel.ProfilePlotPanel`),
neither of which depends on which solver is about to run -- a
`Problem`/`Illumination` and a groove's own shape mean the same thing
regardless. Per tab: a solver's own options, its result, and its own
provenance, since those genuinely differ by method. Only one solve runs
window-wide at a time (see `solve`/`_set_running`) -- there is no second, slow
solver yet to justify concurrent per-tab workers.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QThread, QTimer, Signal
from PySide6.QtWidgets import (
    QMainWindow,
    QMessageBox,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QWidget,
)

from .. import provenance
from ..state import FormErrors, build
from .geometry_panel import GeometryPanel
from .profile_plot_panel import ProfilePlotPanel
from .scalar_tab import ScalarTab
from .setup_tab import SetupTab
from .worker import SolveWorker

__all__ = ["MainWindow"]

#: A solve faster than this never shows a progress bar. The scalar solver takes
#: about 70 ms, and flashing a bar on and off within one frame reads as a
#: glitch rather than as feedback.
_PROGRESS_DELAY_MS = 150

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
#: `show_progress()`; `show_solving(wavelength_count)`; `show_result(scan,
#: energy, lambda_over_period)`; `show_cancelled()`; `show_error(message)`;
#: `show_field_errors(errors)`.
_TAB_FACTORIES = {"scalar": ScalarTab}


class MainWindow(QMainWindow):
    """One geometry, a tab per modeling approach."""

    #: Emitted after a result is drawn. Tests wait on this rather than
    #: sleeping, which is what keeps the widget suite fast and deterministic.
    solved = Signal()

    #: How work reaches the worker. It must be a signal, not a method call:
    #: `moveToThread` only redirects *signal* delivery, so calling
    #: `worker.run(...)` directly would run the solve on the UI thread and
    #: freeze the window -- the exact thing this design exists to prevent.
    _requested = Signal(int, str, object, dict)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("GratingLab")
        self.resize(1180, 820)

        self._parsed = None
        self._token = 0
        self._running = False
        self._active_name: str | None = None

        self._build_menu()
        self._build_layout()
        self._start_worker()
        self._solve_active_tab()

    # -- construction ----------------------------------------------------

    def _build_layout(self) -> None:
        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(self._build_left_column())
        body.addWidget(self._build_right_column())
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([300, 880])
        self.setCentralWidget(body)

    def _build_left_column(self) -> QWidget:
        self.geometry = GeometryPanel()
        self.geometry.solve_requested.connect(self._solve_active_tab)

        scroller = QScrollArea()
        scroller.setWidget(self.geometry)
        scroller.setWidgetResizable(True)
        scroller.setMinimumWidth(280)
        return scroller

    def _build_right_column(self) -> QWidget:
        """Profile above the tabs -- the profile is shared, so it sits once,
        outside whichever tab is active."""
        self.profile_panel = ProfilePlotPanel()

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self.profile_panel)
        right.addWidget(self._build_tabs())
        right.setStretchFactor(0, 0)
        right.setStretchFactor(1, 1)
        return right

    def _build_tabs(self) -> QWidget:
        from ..docs import display_title
        from ...solvers import available_solvers

        self.tabs: dict[str, QWidget] = {}
        self._tab_widget = QTabWidget()

        # Leftmost, always, and never in `self.tabs`: Setup is not a solver,
        # has no solve/cancel contract, and MainWindow never routes to it.
        self._setup_tab = SetupTab()
        self._tab_widget.addTab(self._setup_tab, "Setup")

        for name in available_solvers():
            factory = _TAB_FACTORIES.get(name)
            if factory is None:
                continue  # registered, but no tab written for it yet
            tab = factory()
            tab.solve_requested.connect(lambda _checked=False, n=name: self.solve(n))
            tab.cancel_requested.connect(self.cancel)
            self.tabs[name] = tab
            self._tab_widget.addTab(tab, display_title(name))

        # Setup is leftmost so it is easy to find, but a solver tab -- not
        # the explanatory stub -- is what should be showing, and solving,
        # when the window opens: a user wants a result, not a "nothing here
        # yet" page. `_solve_active_tab` (used for both the construction-time
        # solve and Enter-in-a-geometry-field) reads `currentWidget()`, so
        # this is also what makes the window open already solving.
        if self.tabs:
            self._tab_widget.setCurrentWidget(next(iter(self.tabs.values())))

        return self._tab_widget

    def _build_menu(self) -> None:
        """Help only, and deliberately so.

        This is not a general-purpose application with File/Edit/View
        concerns. It is one window with one job, and the single thing worth a
        menu for is explaining the math, which the window itself never does.
        """
        from ..docs import general_pages, theory_pages
        from ..richtext import menu_label

        # Both the bar and the menu are held on the instance rather than
        # merely added. Under PySide6 6.11 a QMenu reached back through
        # `menuBar().actions()[i].menu()` is not reliably usable -- touching
        # it raises "Internal C++ object already deleted" even while the menu
        # the window built is perfectly alive. Measured, not guessed: holding
        # the QMenuBar makes no difference, holding the QMenu does.
        #
        # So the window owns what it builds, and nothing (tests included)
        # should navigate back down from the menu bar to get at it.
        self._menu_bar = self.menuBar()
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
        self._worker = SolveWorker()
        self._worker.moveToThread(self._thread)
        # Cross-thread, so both directions are queued automatically: requests
        # arrive on the worker's thread, results back on this one.
        self._requested.connect(self._worker.run)
        self._worker.finished.connect(self._on_solved)
        self._worker.failed.connect(self._on_failed)
        self._thread.start()

    # -- actions ---------------------------------------------------------

    def _solve_active_tab(self) -> None:
        """Enter in a geometry field solves whichever tab is showing.

        The geometry panel does not know solver names; it only knows
        something changed. A tab with no `name` (Setup, once it exists) is
        silently not solvable -- there is nothing to solve there.
        """
        name = getattr(self._tab_widget.currentWidget(), "name", None)
        if name is not None:
            self.solve(name)

    def solve(self, name: str) -> None:
        """Validate here, solve on the worker thread.

        Geometry and the tab's own options are validated as two separate
        steps, but both are caught before anything reaches the worker: a
        rejected form is never sent to the worker, whichever kind of error it
        was.
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
        self._requested.emit(self._token, name, geometry, options)

    def cancel(self) -> None:
        """Stop waiting for the running solve.

        It keeps running -- see `worker.py`. Bumping the token is what makes
        its result stale on arrival. There is only ever one thing that could
        be running, so it does not matter which tab's Cancel button was
        clicked; all of them mean the same thing.
        """
        if not self._running:
            return
        self._token += 1
        self._set_running(False)
        if self._active_name is not None:
            self.tabs[self._active_name].show_cancelled()

    def _on_solved(self, token: int, method: str, result) -> None:
        if token != self._token:
            return  # cancelled or superseded; the window has moved on
        self._set_running(False)
        self.profile_panel.draw(self._parsed)
        self.tabs[method].show_result(
            result.scan, result.energy, self._parsed.lambda_over_period
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
        self._thread.quit()
        self._thread.wait(5000)
        super().closeEvent(event)
