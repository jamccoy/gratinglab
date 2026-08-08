"""The window.

Deliberately thin, as its Tk predecessor was. Every number shown comes from
:func:`gratinglab.compare.sweep`; every validation from
:mod:`gratinglab.gui.state`; every line of the provenance panel from
:mod:`gratinglab.gui.provenance`. Nothing here computes physics, so a bug in
this file misplaces a widget rather than producing a wrong answer -- and
`tests/test_gui_qt.py` pins that by comparing what the window plotted against a
direct solver call.

Regions, unchanged in spirit from the Tk version: geometry inputs on the left,
the groove profile and efficiency stacked on the right, provenance along the
bottom. Geometry (:class:`~.geometry_panel.GeometryPanel`) and the profile plot
(:class:`~.profile_plot_panel.ProfilePlotPanel`) are shared, standalone
widgets now -- neither depends on which solver is about to run, since a
`Problem`/`Illumination` and a groove's own shape mean the same thing
regardless. Scalar's own controls (its options, Solve/Cancel/Export, the
efficiency plot) stay directly on this window for now; they move into their
own tab once a second solver exists to prove what "a solver's own tab" should
actually contain.
"""

from __future__ import annotations

import csv
import dataclasses
from pathlib import Path

from PySide6.QtCore import Qt, QThread, QTimer, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLineEdit,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import provenance
from ..scalar_options import ScalarOptionsState, build_options
from ..state import FormErrors, build
from .geometry_panel import GeometryPanel
from .profile_plot_panel import ProfilePlotPanel
from .worker import SolveWorker

__all__ = ["MainWindow"]

#: A solve faster than this never shows a progress bar. The scalar solver takes
#: about 70 ms, and flashing a bar on and off within one frame reads as a
#: glitch rather than as feedback.
_PROGRESS_DELAY_MS = 150


class MainWindow(QMainWindow):
    """Single-window explorer for one solver over one wavelength scan."""

    #: Emitted after a result is drawn. Tests wait on this rather than
    #: sleeping, which is what keeps the widget suite fast and deterministic.
    solved = Signal()

    #: How work reaches the worker. It must be a signal, not a method call:
    #: `moveToThread` only redirects *signal* delivery, so calling
    #: `worker.run(...)` directly would run the solve on the UI thread and
    #: freeze the window -- the exact thing this design exists to prevent.
    #: (token, method, geometry, options) -- method and options are still
    #: always "scalar" and its options dict until a tab exists to vary them.
    _requested = Signal(int, str, object, dict)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("GratingLab")
        self.resize(1180, 820)

        self._fields: dict[str, QWidget] = {}
        self._scan = None
        self._energy = None
        self._parsed = None
        self._token = 0
        self._running = False

        self._build_menu()
        self._build_layout()
        self._start_worker()
        self.solve()

    # -- construction ----------------------------------------------------

    def _build_layout(self) -> None:
        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(self._build_left_column())
        body.addWidget(self._build_right_column())
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([300, 880])

        outer = QSplitter(Qt.Orientation.Vertical)
        outer.addWidget(body)
        outer.addWidget(self._build_provenance())
        outer.setStretchFactor(0, 1)
        outer.setStretchFactor(1, 0)
        outer.setSizes([620, 200])

        self.setCentralWidget(outer)

    def _build_left_column(self) -> QWidget:
        """Geometry, then scalar's own controls, stacked -- the same visual
        grouping the flat window had before the split, just now assembled
        from a shared `GeometryPanel` plus this window's own remaining
        solver-specific widgets."""
        self.geometry = GeometryPanel()
        self.geometry.solve_requested.connect(self.solve)

        panel = QWidget()
        column = QVBoxLayout(panel)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self.geometry)
        column.addWidget(self._build_solver_controls())

        # The column outgrows a short window once the order panel joins it.
        scroller = QScrollArea()
        scroller.setWidget(panel)
        scroller.setWidgetResizable(True)
        scroller.setMinimumWidth(280)
        return scroller

    def _build_solver_controls(self) -> QWidget:
        """Scalar's own options and actions.

        Everything here moves into its own tab once a second solver exists to
        prove what that tab should actually contain -- see M11-C.
        """
        scalar_defaults = ScalarOptionsState()
        panel = QWidget()
        column = QVBoxLayout(panel)
        column.setContentsMargins(0, 0, 0, 0)

        solver = QGroupBox("Scalar solver")
        solver_form = QFormLayout(solver)
        quadrature = QLineEdit(scalar_defaults.quadrature_points)
        quadrature.returnPressed.connect(self.solve)
        self._fields["quadrature_points"] = quadrature
        solver_form.addRow("Quadrature pts", quadrature)

        self._solve_button = QPushButton("Solve")
        self._solve_button.clicked.connect(self.solve)
        self._cancel_button = QPushButton("Cancel")
        self._cancel_button.clicked.connect(self.cancel)
        self._cancel_button.setEnabled(False)
        self._export_button = QPushButton("Export CSV…")
        self._export_button.clicked.connect(self.export_csv)

        self._progress = QProgressBar()
        self._progress.setRange(0, 0)  # indeterminate: sweep reports nothing
        self._progress.hide()

        buttons = QHBoxLayout()
        buttons.addWidget(self._solve_button)
        buttons.addWidget(self._cancel_button)

        column.addWidget(solver)
        column.addLayout(buttons)
        column.addWidget(self._export_button)
        column.addWidget(self._progress)
        column.addStretch(1)

        self._check_fields_match_scalar_options()
        return panel

    def _check_fields_match_scalar_options(self) -> None:
        """Refuse to open rather than fail later on one code path.

        `_read_scalar_options` builds `ScalarOptionsState` by keyword from
        this dict, so a renamed or forgotten field is a TypeError at solve
        time -- visible only to whoever presses Solve. Checking here makes it
        immediate and says which field. `GeometryPanel` runs the equivalent
        check for its own fields against `FormState`.
        """
        declared = {f.name for f in dataclasses.fields(ScalarOptionsState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "scalar-option fields do not match ScalarOptionsState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    def _build_right_column(self) -> QWidget:
        """Profile above efficiency -- the same stacked arrangement the flat
        window had as two subplots of one figure, now two separate plots so
        the profile can be shared once solver tabs exist."""
        self.profile_panel = ProfilePlotPanel()

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self.profile_panel)
        right.addWidget(self._build_efficiency_plot())
        right.setStretchFactor(0, 0)
        right.setStretchFactor(1, 1)
        return right

    def _build_efficiency_plot(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure

        self._efficiency_figure = Figure(figsize=(8, 3), layout="constrained")
        self._efficiency_axes = self._efficiency_figure.add_subplot(1, 1, 1)
        self._efficiency_canvas = FigureCanvasQTAgg(self._efficiency_figure)

        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self._efficiency_canvas, panel))
        layout.addWidget(self._efficiency_canvas)
        return panel

    def _build_provenance(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 4, 8, 8)
        heading = QLabel("Provenance")
        heading.setStyleSheet("font-weight: bold")
        self._provenance = QTextEdit()
        self._provenance.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._provenance)
        return panel

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

    # -- reactions -------------------------------------------------------

    def _read_scalar_options(self) -> ScalarOptionsState:
        return ScalarOptionsState(
            quadrature_points=_value(self._fields["quadrature_points"])
        )

    # -- actions ---------------------------------------------------------

    def solve(self) -> None:
        """Validate here, solve on the worker thread.

        Geometry and scalar's own options are validated as two separate
        steps, but both are caught before anything reaches the worker --
        preserving "rejected synchronously, never reaches the worker" for
        either kind of error, exactly as it was when there was only one
        dataclass to validate.
        """
        if self._running:
            return

        try:
            geometry = build(self.geometry.read_form())
            options = build_options(
                geometry.problem, geometry.illumination, geometry.wavelengths,
                self._read_scalar_options(),
            )
        except FormErrors as exc:
            self._paint(provenance.error_lines(exc.errors))
            return

        self._parsed = geometry
        self._token += 1
        self._set_running(True)
        self._paint(provenance.solving_lines("scalar", len(geometry.wavelengths)))
        QTimer.singleShot(_PROGRESS_DELAY_MS, self._show_progress_if_still_running)
        self._requested.emit(self._token, "scalar", geometry, options)

    def cancel(self) -> None:
        """Stop waiting for the running solve.

        It keeps running -- see `worker.py`. Bumping the token is what makes
        its result stale on arrival.
        """
        if not self._running:
            return
        self._token += 1
        self._set_running(False)
        if self._scan is not None:
            self._paint(
                provenance.provenance_lines(
                    self._scan,
                    self._energy,
                    self._parsed.lambda_over_period,
                    cancelled=True,
                )
            )

    def _on_solved(self, token: int, method: str, result) -> None:
        # `method` is unused until a second solver tab exists to route
        # between -- for now there is exactly one tab and it is always
        # "scalar". Accepted here so the signal signature is already what
        # M11-C needs.
        if token != self._token:
            return  # cancelled or superseded; the window has moved on
        self._set_running(False)
        self._scan, self._energy = result.scan, result.energy

        self.profile_panel.draw(self._parsed)
        self._draw_efficiency(result.scan)
        self._efficiency_canvas.draw_idle()
        self._paint(
            provenance.provenance_lines(
                result.scan, result.energy, self._parsed.lambda_over_period
            )
        )
        self.solved.emit()

    def _on_failed(self, token: int, method: str, message: str) -> None:
        if token != self._token:
            return
        self._set_running(False)
        self._paint((provenance.Line(f"solve failed: {message}\n", "bad"),))

    def _set_running(self, running: bool) -> None:
        self._running = running
        self._solve_button.setEnabled(not running)
        self._cancel_button.setEnabled(running)
        if not running:
            self._progress.hide()

    def _show_progress_if_still_running(self) -> None:
        if self._running:
            self._progress.show()

    def export_csv(self) -> None:
        if self._scan is None:
            QMessageBox.information(self, "Nothing to export", "Solve first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export efficiencies", "efficiency.csv", "CSV (*.csv)"
        )
        if not path:
            return
        rows = self._scan.to_records()
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        QMessageBox.information(
            self, "Exported", f"{len(rows)} rows to {Path(path).name}"
        )

    def show_about(self) -> None:
        from ... import __version__

        QMessageBox.about(self, "About GratingLab", provenance.about_text(__version__))

    def show_theory(self, page) -> None:
        """Open the read-only viewer for one theory or reference page."""
        from .theory_viewer import TheoryViewer

        viewer = TheoryViewer(page, parent=self)
        viewer.show()
        return viewer

    # -- drawing ---------------------------------------------------------

    def _paint(self, lines) -> None:
        self._provenance.setHtml(provenance.to_html(lines))

    def _draw_efficiency(self, scan) -> None:
        axes = self._efficiency_axes
        axes.clear()
        for index, order in enumerate(scan.orders):
            values = scan.efficiency[:, index]
            if values.max() < 1e-4:
                continue
            axes.plot(scan.wavelengths, values, lw=1.3, label=f"m={int(order):+d}")
        axes.plot(scan.wavelengths, scan.total, "k--", lw=1.0, alpha=0.7, label="Σ")
        axes.axhline(1.0, color="#b00020", lw=0.8, ls=":", alpha=0.7)
        axes.set_xlabel("wavelength (nm)")
        axes.set_ylabel("efficiency")
        axes.set_ylim(bottom=0)
        axes.grid(alpha=0.25)
        handles, _ = axes.get_legend_handles_labels()
        if handles:
            axes.legend(fontsize=7, ncol=max(1, len(handles) // 8), loc="upper right")

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


def _value(widget: QWidget) -> str:
    """The text of an input, whichever kind it is."""
    if isinstance(widget, QComboBox):
        return widget.currentText()
    return widget.text()


