"""The window.

Deliberately thin, as its Tk predecessor was. Every number shown comes from
:func:`gratinglab.compare.sweep`; every validation from
:mod:`gratinglab.gui.state`; every line of the provenance panel from
:mod:`gratinglab.gui.provenance`. Nothing here computes physics, so a bug in
this file misplaces a widget rather than producing a wrong answer -- and
`tests/test_gui_qt.py` pins that by comparing what the window plotted against a
direct solver call.

Three regions, unchanged from the Tk version: inputs on the left, the groove
profile and efficiency stacked on the right, provenance along the bottom. The
splitters are new, and are the one thing `pack()` could not give -- the panel
holds warnings several lines long and used to be a fixed seven lines tall.
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
    QLabel,
    QLineEdit,
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
from ..state import (
    ANGLE_LABELS,
    MOUNTS,
    PROFILE_FIELDS,
    PROFILE_KINDS,
    FormErrors,
    FormState,
    build,
)
from .worker import SolveWorker

__all__ = ["MainWindow"]

#: Fields whose visibility depends on the selected profile kind.
_PROFILE_ROW_KEYS = ("blaze_angle", "antiblaze_angle", "depth_fraction", "duty_cycle")

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
    _requested = Signal(int, object)

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

        self._on_mount_change()
        self._on_profile_change()
        self.solve()

    # -- construction ----------------------------------------------------

    def _build_layout(self) -> None:
        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(self._build_inputs())
        body.addWidget(self._build_plots())
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

    def _build_inputs(self) -> QWidget:
        defaults = FormState()
        panel = QWidget()
        column = QVBoxLayout(panel)

        def entry(form: QFormLayout, label: str, key: str) -> QLineEdit:
            widget = QLineEdit(getattr(defaults, key))
            widget.returnPressed.connect(self.solve)
            self._fields[key] = widget
            form.addRow(label, widget)
            return widget

        def combo(form: QFormLayout, label: str, key: str, values, on_change) -> QComboBox:
            widget = QComboBox()
            widget.addItems(list(values))
            widget.setCurrentText(getattr(defaults, key))
            widget.currentTextChanged.connect(lambda _t: on_change())
            self._fields[key] = widget
            form.addRow(label, widget)
            return widget

        grating = QGroupBox("Grating")
        grating_form = self._grating_form = QFormLayout(grating)
        entry(grating_form, "Period (nm)", "period")
        combo(grating_form, "Profile", "profile_kind", PROFILE_KINDS,
              self._on_profile_change)

        self._profile_rows: dict[str, QWidget] = {}
        for key, label in (
            ("blaze_angle", "Blaze δ (deg)"),
            ("antiblaze_angle", "Anti-blaze (deg)"),
            ("depth_fraction", "Depth / period"),
            ("duty_cycle", "Duty cycle"),
        ):
            self._profile_rows[key] = entry(grating_form, label, key)

        # Not a visible row: the path is set by the file dialog, and shown by
        # the label beneath it. It still has to be a field, because _read_form
        # builds FormState from exactly this dict.
        path_field = QLineEdit("")
        path_field.hide()
        self._fields["profile_path"] = path_field

        self._path_label = QLabel("")
        self._path_label.setWordWrap(True)
        self._path_label.setStyleSheet("color: #555")
        self._load_button = QPushButton("Load profile…")
        self._load_button.clicked.connect(self.load_profile)
        grating_form.addRow(self._load_button)
        grating_form.addRow(self._path_label)

        mount = QGroupBox("Mount")
        mount_form = QFormLayout(mount)
        combo(mount_form, "Geometry", "mount", MOUNTS, self._on_mount_change)
        self._angle_labels: dict[str, QLabel] = {}
        for key in ("alpha", "gamma"):
            widget = QLineEdit(getattr(defaults, key))
            widget.returnPressed.connect(self.solve)
            self._fields[key] = widget
            label = QLabel("")
            mount_form.addRow(label, widget)
            self._angle_labels[key] = label

        scan = QGroupBox("Wavelengths (nm)")
        scan_form = QFormLayout(scan)
        entry(scan_form, "Start", "wavelength_start")
        entry(scan_form, "Stop", "wavelength_stop")
        entry(scan_form, "Points", "wavelength_count")

        solver = QGroupBox("Scalar solver")
        solver_form = QFormLayout(solver)
        entry(solver_form, "Quadrature pts", "quadrature_points")

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

        for group in (grating, mount, scan, solver):
            column.addWidget(group)
        column.addLayout(buttons)
        column.addWidget(self._export_button)
        column.addWidget(self._progress)
        column.addStretch(1)

        self._check_fields_match_formstate()

        # The column outgrows a short window once the order panel joins it.
        scroller = QScrollArea()
        scroller.setWidget(panel)
        scroller.setWidgetResizable(True)
        scroller.setMinimumWidth(280)
        return scroller

    def _check_fields_match_formstate(self) -> None:
        """Refuse to open rather than fail later on one code path.

        `_read_form` builds a `FormState` by keyword from this dict, so a
        renamed or forgotten field is a TypeError at solve time -- visible only
        to whoever presses Solve. Checking here makes it immediate and says
        which field.
        """
        declared = {f.name for f in dataclasses.fields(FormState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "form fields do not match FormState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    def _build_plots(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 6), layout="constrained")
        self._profile_axes = self._figure.add_subplot(2, 1, 1)
        self._efficiency_axes = self._figure.add_subplot(2, 1, 2)
        self._canvas = FigureCanvasQTAgg(self._figure)

        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self._canvas, panel))
        layout.addWidget(self._canvas)
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

    def _on_profile_change(self) -> None:
        """Show only the parameters the selected profile actually uses."""
        needed = PROFILE_FIELDS.get(
            self._fields["profile_kind"].currentText(), frozenset()
        )
        for key, widget in self._profile_rows.items():
            # setRowVisible hides the label with the field. Hiding the field
            # alone would leave a caption for a control that is not there.
            self._grating_form.setRowVisible(widget, key in needed)
        from_file = "profile_path" in needed
        self._load_button.setVisible(from_file)
        self._path_label.setVisible(from_file)

    def _on_mount_change(self) -> None:
        """Relabel the two angle fields; the mount decides what they mean."""
        primary, secondary = ANGLE_LABELS[self._fields["mount"].currentText()]
        self._angle_labels["alpha"].setText(primary)
        self._angle_labels["alpha"].show()
        self._fields["alpha"].show()

        if secondary is None:
            self._angle_labels["gamma"].hide()
            self._fields["gamma"].hide()
        else:
            self._angle_labels["gamma"].setText(secondary)
            self._angle_labels["gamma"].show()
            self._fields["gamma"].show()

    def _read_form(self) -> FormState:
        return FormState(**{k: _value(w) for k, w in self._fields.items()})

    # -- actions ---------------------------------------------------------

    def solve(self) -> None:
        """Validate here, solve on the worker thread."""
        if self._running:
            return

        try:
            parsed = build(self._read_form())
        except FormErrors as exc:
            self._paint(provenance.error_lines(exc.errors))
            return

        self._parsed = parsed
        self._token += 1
        self._set_running(True)
        self._paint(provenance.solving_lines("scalar", len(parsed.wavelengths)))
        QTimer.singleShot(_PROGRESS_DELAY_MS, self._show_progress_if_still_running)
        self._requested.emit(self._token, parsed)

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

    def _on_solved(self, token: int, result) -> None:
        if token != self._token:
            return  # cancelled or superseded; the window has moved on
        self._set_running(False)
        self._scan, self._energy = result.scan, result.energy

        self._draw_profile(self._parsed)
        self._draw_efficiency(result.scan)
        self._canvas.draw_idle()
        self._paint(
            provenance.provenance_lines(
                result.scan, result.energy, self._parsed.lambda_over_period
            )
        )
        self.solved.emit()

    def _on_failed(self, token: int, message: str) -> None:
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

    def load_profile(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load boundary profile", "", "PCGrate boundary (*.ggp);;All files (*)"
        )
        if path:
            self._fields["profile_path"].setText(path)
            self._path_label.setText(Path(path).name)
            self.solve()

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

    def _draw_profile(self, parsed) -> None:
        import numpy as np

        axes = self._profile_axes
        axes.clear()
        # Drawn from the same Profile the solver integrates, so what is shown
        # is literally what is computed.
        t = np.linspace(0.0, 1.0, 600, endpoint=False)
        two_periods = np.concatenate([t, t + 1.0])
        height = parsed.problem.height_nm(two_periods)
        axes.plot(two_periods, height, color="#1f3b73", lw=1.8)
        axes.fill_between(two_periods, 0, height, color="#1f3b73", alpha=0.12)
        axes.set_xlabel("position / period  (two periods shown)")
        axes.set_ylabel("height (nm)")
        axes.set_title(
            f"{type(parsed.problem.profile).__name__} — "
            f"period {parsed.problem.period:g} nm, depth {parsed.problem.depth:.4g} nm",
            fontsize=10,
        )
        axes.grid(alpha=0.25)

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


