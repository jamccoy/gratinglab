"""
The window.

Thin by intention: it reads controls into a `FormState`, asks `gui.state` to turn
that into settings, hands the work to a thread, and draws what comes back. Any
question that can be answered without a window belongs in `gui/state.py`.
"""
from __future__ import annotations

import contextlib
import io
import os

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFrame,
    QGroupBox, QHBoxLayout, QLabel, QMainWindow, QPushButton, QScrollArea,
    QSpinBox, QStatusBar, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

from ...config import PROJECT_ROOT
from ...core.processing import load_afm_data, raw_data, raw_data_multi_group
from ...settings import AnalysisSettings, MAX_FACET_TRIM, VALID_BLAZE_SIDES
from ..state import FormState, build, summarize_result
from .canvas import PlotCanvas
from .worker import AnalysisWorker

__all__ = ["MainWindow"]


class MainWindow(QMainWindow):
    """Interactive blaze-angle analysis of a single AFM scan."""

    _request_analysis = Signal(str, object, int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("AFM Blaze Angle Analysis")
        self.resize(1280, 760)

        self._defaults = AnalysisSettings.from_config()
        self._data = None
        self._disp_um = None
        self._profile_nm = None
        self._scan_size = None
        self._filename = None
        self._result = None
        self._settings = None
        self._token = 0          # bumped to abandon an in-flight result

        self._build_ui()
        self._start_worker()

    # ── Worker ───────────────────────────────────────────────────────────────

    def _start_worker(self):
        self._thread = QThread(self)
        self._worker = AnalysisWorker()
        self._worker.moveToThread(self._thread)
        self._request_analysis.connect(self._worker.run)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._thread.start()

    def closeEvent(self, event):
        self._token += 1
        self._thread.quit()
        self._thread.wait(3000)
        super().closeEvent(event)

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        left = QWidget()
        left.setFixedWidth(272)
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(8)
        left_layout.addWidget(self._build_file_group())
        left_layout.addWidget(self._build_analysis_group())
        left_layout.addWidget(self._build_view_group())
        left_layout.addStretch()

        scroll = QScrollArea()
        scroll.setWidget(left)
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(292)
        scroll.setFrameShape(QFrame.NoFrame)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = PlotCanvas(self)
        right_layout.addWidget(NavigationToolbar(self.canvas, self))
        right_layout.addWidget(self.canvas, stretch=1)
        right_layout.addWidget(self._build_results_panel())

        root.addWidget(scroll)
        divider = QFrame()
        divider.setFrameShape(QFrame.VLine)
        divider.setFrameShadow(QFrame.Sunken)
        root.addWidget(divider)
        root.addWidget(right, stretch=1)

        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self._set_status("Ready. Open an AFM file to begin.")

    def _build_file_group(self):
        group = QGroupBox("File")
        layout = QVBoxLayout(group)
        self.file_label = QLabel("No file loaded")
        self.file_label.setWordWrap(True)
        self.file_label.setStyleSheet("color: gray; font-size: 11px;")
        layout.addWidget(self.file_label)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse_file)
        layout.addWidget(browse)

        layout.addWidget(QLabel("Fallback scan width (µm):"))
        self.scan_size_spin = QDoubleSpinBox()
        self.scan_size_spin.setRange(0.1, 1000.0)
        self.scan_size_spin.setDecimals(3)
        self.scan_size_spin.setSingleStep(0.1)
        self.scan_size_spin.setValue(self._defaults.scan_x_size)
        self.scan_size_spin.setToolTip(
            "Used only when the scan width cannot be read from the file header.")
        layout.addWidget(self.scan_size_spin)

        self.info_label = QLabel("—")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("font-size: 11px;")
        layout.addWidget(self.info_label)
        return group

    def _build_analysis_group(self):
        group = QGroupBox("Analysis parameters")
        layout = QVBoxLayout(group)

        layout.addWidget(QLabel("Estimated period (nm):"))
        self.period_spin = QDoubleSpinBox()
        self.period_spin.setRange(1.0, 100000.0)
        self.period_spin.setDecimals(2)
        self.period_spin.setValue(self._defaults.period_est)
        self.period_spin.setToolTip("Groove spacing. Must match your grating.")
        layout.addWidget(self.period_spin)

        layout.addWidget(QLabel("Facet trim (fraction):"))
        self.trim_spin = QDoubleSpinBox()
        # The hard limit is enforced in gui.state, not here; the spin box range
        # only keeps the common case convenient.
        self.trim_spin.setRange(0.0, MAX_FACET_TRIM)
        self.trim_spin.setSingleStep(0.05)
        self.trim_spin.setDecimals(2)
        self.trim_spin.setValue(self._defaults.facet_trim)
        self.trim_spin.setToolTip(
            "Fraction trimmed from each end of the facet before fitting.\n"
            "Too small risks measuring the rounded groove top and reading low.\n"
            "This parameter moves the answer more than any other.")
        layout.addWidget(self.trim_spin)

        layout.addWidget(QLabel("Blaze side:"))
        self.side_combo = QComboBox()
        self.side_combo.addItems(list(VALID_BLAZE_SIDES))
        self.side_combo.setCurrentText(self._defaults.blaze_side)
        self.side_combo.setToolTip(
            "Which facet to measure, chosen by slope sign rather than position.")
        layout.addWidget(self.side_combo)

        layout.addWidget(QLabel("Edge exclusion (periods):"))
        self.edge_spin = QDoubleSpinBox()
        self.edge_spin.setRange(0.0, 3.0)
        self.edge_spin.setSingleStep(0.1)
        self.edge_spin.setDecimals(2)
        self.edge_spin.setValue(self._defaults.edge_exclusion_periods)
        self.edge_spin.setToolTip(
            "Reject grooves this close to either end of the scan line.\n"
            "Their facet is clipped by the edge, so the fitted angle is\n"
            "meaningless. 0 disables the check.")
        layout.addWidget(self.edge_spin)

        self.row_groups_check = QCheckBox("Use row groups")
        self.row_groups_check.setChecked(self._defaults.use_row_groups)
        self.row_groups_check.setToolTip(
            "Analyse N horizontal bands separately instead of averaging the\n"
            "whole image. Many more measurements, but they re-measure the same\n"
            "grooves and are not independent - see the ICC diagnostic.")
        layout.addWidget(self.row_groups_check)

        n_layout = QHBoxLayout()
        n_layout.addWidget(QLabel("N groups:"))
        self.n_groups_spin = QSpinBox()
        self.n_groups_spin.setRange(2, 200)
        self.n_groups_spin.setValue(self._defaults.n_row_groups)
        n_layout.addWidget(self.n_groups_spin)
        layout.addLayout(n_layout)

        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setEnabled(False)
        font = QFont()
        font.setBold(True)
        self.run_btn.setFont(font)
        self.run_btn.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_btn)
        return group

    def _build_view_group(self):
        group = QGroupBox("View")
        layout = QVBoxLayout(group)
        self.raw_profile_btn = QPushButton("Raw Profile (1D)")
        self.raw_profile_btn.clicked.connect(self.show_raw_profile)
        self.topography_btn = QPushButton("2D Topography")
        self.topography_btn.clicked.connect(self.show_2d)
        self.row_groups_btn = QPushButton("Row Groups")
        self.row_groups_btn.clicked.connect(self.show_row_groups)
        self.detection_btn = QPushButton("Groove Detection")
        self.detection_btn.clicked.connect(self.show_detection)
        self.angles_btn = QPushButton("Blaze Angles")
        self.angles_btn.clicked.connect(self.show_angles)

        self._file_buttons = [self.raw_profile_btn, self.topography_btn,
                              self.row_groups_btn]
        self._result_buttons = [self.detection_btn, self.angles_btn]
        for b in self._file_buttons + self._result_buttons:
            b.setEnabled(False)
            layout.addWidget(b)
        return group

    def _build_results_panel(self):
        group = QGroupBox("Results")
        layout = QVBoxLayout(group)
        self.results_label = QLabel("Load a file and press Run Analysis.")
        self.results_label.setStyleSheet("font-family: monospace; font-size: 11px;")
        self.results_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.results_label)
        group.setMaximumHeight(130)
        return group

    # ── Form ─────────────────────────────────────────────────────────────────

    def form_state(self) -> FormState:
        """Current control values, as plain data"""
        return FormState(
            period_est=self.period_spin.value(),
            facet_trim=self.trim_spin.value(),
            blaze_side=self.side_combo.currentText(),
            edge_exclusion_periods=self.edge_spin.value(),
            use_row_groups=self.row_groups_check.isChecked(),
            n_row_groups=self.n_groups_spin.value(),
            scan_x_size=self.scan_size_spin.value(),
        )

    # ── Actions ──────────────────────────────────────────────────────────────

    def _browse_file(self):
        start_dir = os.path.join(PROJECT_ROOT, 'data')
        if not os.path.isdir(start_dir):
            start_dir = ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open AFM File", start_dir,
            "AFM Data Files (*.txt *.dat *.asc);;All Files (*)")
        if path:
            self.load(path)

    def load(self, path):
        try:
            self._set_status(f"Loading {os.path.basename(path)}…")
            QApplication.processEvents()

            with contextlib.redirect_stdout(io.StringIO()):
                data, scan_size = load_afm_data(
                    path, default_scan_size=self.scan_size_spin.value())
                disp_um, profile_nm = raw_data(data, scan_size)

            self._data, self._scan_size = data, scan_size
            self._disp_um, self._profile_nm = disp_um, profile_nm
            self._filename = path
            self._result = self._settings = None

            self.file_label.setText(os.path.basename(path))
            self.file_label.setStyleSheet("font-size: 11px;")
            self.scan_size_spin.setValue(scan_size)

            rows, cols = data.shape
            self.info_label.setText(
                f"Shape: {rows} × {cols} px\n"
                f"Scan width: {scan_size:.3f} µm\n"
                f"Height range: {profile_nm.min():.1f} – {profile_nm.max():.1f} nm")

            for b in self._file_buttons:
                b.setEnabled(True)
            for b in self._result_buttons:
                b.setEnabled(False)
            self.run_btn.setEnabled(True)
            self.results_label.setText("Press Run Analysis.")

            self.n_groups_spin.setMaximum(max(2, rows // 3))
            self.show_raw_profile()
            self._set_status(f"Loaded: {os.path.basename(path)}")

        except Exception as exc:
            self._set_status(f"Error loading file: {exc}")
            self.info_label.setText(f"Error:\n{exc}")

    def run_analysis(self):
        """Validate the form, then hand the work to the thread."""
        if self._filename is None:
            return

        settings, errors = build(self.form_state())
        if errors:
            self.results_label.setText(
                "\n".join(f"{e.field}: {e.message}" for e in errors))
            self._set_status("Fix the highlighted settings and try again.")
            return

        self._token += 1
        self.run_btn.setEnabled(False)
        self._set_status("Running analysis…")
        self._request_analysis.emit(self._filename, settings, self._token)

    def _on_finished(self, result, settings, token):
        if token != self._token:
            return  # a later run superseded this one
        self.run_btn.setEnabled(True)

        if result is None:
            self._result = self._settings = None
            for b in self._result_buttons:
                b.setEnabled(False)
            self.results_label.setText(
                "No blaze angles could be extracted.\n"
                "Check that the estimated period matches your grating.")
            self._set_status("Analysis produced no measurements.")
            return

        self._result, self._settings = result, settings
        for b in self._result_buttons:
            b.setEnabled(True)
        self.results_label.setText(summarize_result(result))
        self.show_angles()
        self._set_status(f"Analysis complete: {result['n_grooves']} measurements, "
                         f"mean {result['mean_angle']:.2f}°")

    def _on_failed(self, message, token):
        if token != self._token:
            return
        self.run_btn.setEnabled(True)
        self._result = self._settings = None
        self.results_label.setText(f"Analysis failed:\n{message}")
        self._set_status("Analysis failed.")

    # ── Views ────────────────────────────────────────────────────────────────

    def show_raw_profile(self):
        if self._data is not None:
            self.canvas.plot_raw(self._disp_um, self._profile_nm,
                                 self._filename, self._scan_size)

    def show_2d(self):
        if self._data is not None:
            self.canvas.plot_2d(self._data, self._scan_size, self._filename)

    def show_row_groups(self):
        if self._data is None:
            return
        with contextlib.redirect_stdout(io.StringIO()):
            disp_um, profiles, info = raw_data_multi_group(
                self._data, self._scan_size, n_groups=self.n_groups_spin.value())
        self.canvas.plot_row_groups(disp_um, profiles, info,
                                    self._filename, self._scan_size)
        self._set_status(f"{info['n_groups']} groups × {info['rows_per_group']} rows")

    def show_detection(self):
        if self._result is not None:
            self.canvas.plot_groove_detection(self._result, self._filename,
                                              self._settings)

    def show_angles(self):
        if self._result is not None:
            self.canvas.plot_blaze_angles(self._result, self._filename)

    def _set_status(self, msg):
        self.status.showMessage(msg)
