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
    QSpinBox, QStatusBar, QTabWidget, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

from ...config import PROJECT_ROOT
from ...core.processing import load_afm_data, raw_data, raw_data_multi_group
from ...settings import (AnalysisSettings, MAX_FACET_TRIM, VALID_BLAZE_SIDES,
                         VALID_SPM_DIRECTIONS)
from ..state import FormState, build, summarize_result
from .canvas import PlotCanvas
from .boundary_view import BoundaryView
from .import_view import ImportView
from .wiki_view import WikiView
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
        # Tabs rather than a single pane: the analysis and the prose explaining
        # it are both things a user wants to return to. The Analysis tab is the
        # previous central widget, moved wholesale - no control, view or signal
        # changed, which is what lets the existing GUI tests prove the move was
        # clean.
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Import first: a scan enters the application here, and the Analysis tab
        # consumes whatever this produces.
        self.importer = ImportView()
        self.importer.dataChanged.connect(self._on_import_changed)
        self.tabs.addTab(self.importer, "Import")

        central = QWidget()
        root = QHBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        left = QWidget()
        left.setFixedWidth(272)
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(8)
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

        self.tabs.addTab(central, "Analysis")

        self.boundary = BoundaryView()
        self.tabs.addTab(self.boundary, "Boundary")

        self.wiki = WikiView()
        self.tabs.addTab(self.wiki, "Wiki")

        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self._set_status("Ready. Open an AFM file to begin.")

    def _build_analysis_group(self):
        group = QGroupBox("Analysis parameters")
        layout = QVBoxLayout(group)

        layout.addWidget(QLabel("Estimated period (nm):"))
        self.period_spin = QDoubleSpinBox()
        self.period_spin.setRange(1.0, 100000.0)
        self.period_spin.setDecimals(2)
        self.period_spin.setValue(self._defaults.period_est)
        self.period_spin.setToolTip("Groove spacing. Must match your grating.")
        # level_grooves needs the period to find what it levels on, so the
        # Import preview has to know when this changes.
        self.period_spin.valueChanged.connect(
            lambda v: self.importer.set_period_est(v))
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
        overrides = self.importer.settings_overrides()
        return FormState(
            period_est=self.period_spin.value(),
            facet_trim=self.trim_spin.value(),
            blaze_side=self.side_combo.currentText(),
            edge_exclusion_periods=self.edge_spin.value(),
            use_row_groups=self.row_groups_check.isChecked(),
            n_row_groups=self.n_groups_spin.value(),
            scan_x_size=overrides['scan_x_size'],
            spm_direction=overrides['spm_direction'],
            image_flatten_method=overrides['image_flatten_method'],
            flatten_method=overrides['flatten_method'],
            flatten_poly_order=overrides['flatten_poly_order'],
            flatten_feature=overrides['flatten_feature'],
            flatten_exclude_edges=overrides['flatten_exclude_edges'],
        )

    # ── Actions ──────────────────────────────────────────────────────────────

    def _on_import_changed(self):
        """
        Take whatever the Import tab produced.

        The Import tab owns loading and both flattening stages; this window keeps
        one copy of the current scan so the Analysis views and the worker all see
        the same thing.
        """
        self._data = self.importer.data
        self._scan_size = self.importer.scan_size
        self._filename = self.importer.filename
        self._result = self._settings = None

        if self._data is None:
            for button in self._file_buttons + self._result_buttons:
                button.setEnabled(False)
            self.run_btn.setEnabled(False)
            self.canvas.show_placeholder("Load an AFM file in the Import tab")
            self.boundary.set_scan(None, None, None, None)
            self._set_status("No data loaded.")
            return

        self._disp_um, self._profile_nm = raw_data(self._data, self._scan_size)

        for button in self._file_buttons:
            button.setEnabled(True)
        for button in self._result_buttons:
            button.setEnabled(False)
        self.run_btn.setEnabled(True)
        self.results_label.setText("Press Run Analysis.")

        # The boundary export works from the same scan and the same settings,
        # so it sees whatever Import produced too.
        settings, errors = build(self.form_state())
        if settings is not None:
            self.boundary.set_scan(self._data, self._scan_size,
                                   self._filename, settings)

        self.n_groups_spin.setMaximum(max(2, self._data.shape[0] // 3))
        self.show_raw_profile()
        self._set_status(f"Loaded: {os.path.basename(self._filename)}  "
                         f"(image flattening: "
                         f"{self.importer.settings_overrides()['image_flatten_method']})")

    def load(self, path):
        """Load a file. Delegates to the Import tab, which owns loading."""
        self.importer.load(path)

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
