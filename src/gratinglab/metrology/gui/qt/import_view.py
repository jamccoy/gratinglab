"""
The Import tab: load a scan, flatten it, see what that did.

Owns file loading and both flattening stages. The Analysis tab consumes whatever
this produces, so there is one place a scan enters the application.

The preview is the reason this tab exists. Profile flattening moves the mean
blaze angle by about 0.49 degrees across its four methods - comparable to the
differences between samples - and until now it was a config-file setting with no
way to see what it was removing. Here the fitted background is drawn over the
profile it is about to be subtracted from.
"""
from __future__ import annotations

import contextlib
import io
import os

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

import numpy as np
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QFileDialog, QFrame, QGroupBox,
    QHBoxLayout, QLabel, QPushButton, QScrollArea, QSpinBox, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ...config import PROJECT_ROOT
from ...core.image_flatten import (VALID_IMAGE_FLATTEN_METHODS, flatten_image,
                                   row_offset_spread)
from ...core.processing import flatten_profile, load_afm_data, raw_data
from ...settings import (AnalysisSettings, VALID_FLATTEN_FEATURES,
                         VALID_FLATTEN_METHODS, VALID_SPM_DIRECTIONS)

__all__ = ["ImportView"]


class ImportView(QWidget):
    """Load a scan, choose how it is flattened, and see the result."""

    #: A scan was loaded, or a flattening choice changed. The window re-reads
    #: this tab's state rather than the signal carrying it, so there is one
    #: source of truth.
    dataChanged = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self._defaults = AnalysisSettings.from_config()
        self._raw = None          # as loaded, before any flattening
        self._flat = None         # after image flattening; what analysis sees
        self._scan_size = None
        self._filename = None
        self._period_est = self._defaults.period_est

        self._build_ui()

    # ── UI ───────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        left = QWidget()
        left.setFixedWidth(272)
        column = QVBoxLayout(left)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(8)
        column.addWidget(self._build_file_group())
        column.addWidget(self._build_image_group())
        column.addWidget(self._build_profile_group())
        column.addStretch()

        scroll = QScrollArea()
        scroll.setWidget(left)
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(292)
        scroll.setFrameShape(QFrame.NoFrame)

        self.figure = Figure(figsize=(8, 6), tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self._show_placeholder()

        root.addWidget(scroll)
        divider = QFrame()
        divider.setFrameShape(QFrame.VLine)
        divider.setFrameShadow(QFrame.Sunken)
        root.addWidget(divider)
        root.addWidget(self.canvas, stretch=1)

    def _build_file_group(self):
        group = QGroupBox("File")
        layout = QVBoxLayout(group)

        self.file_label = QLabel("No file loaded")
        self.file_label.setWordWrap(True)
        self.file_label.setStyleSheet("color: gray; font-size: 11px;")
        layout.addWidget(self.file_label)

        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
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

        layout.addWidget(QLabel("Scan direction (.spm only):"))
        self.direction_combo = QComboBox()
        self.direction_combo.addItems(list(VALID_SPM_DIRECTIONS))
        self.direction_combo.setCurrentText(self._defaults.spm_direction)
        self.direction_combo.setEnabled(False)
        self.direction_combo.setToolTip(
            "Which pass of the tip to analyse. A Nanoscope file records both.\n"
            "Retrace matches the project's existing Gwyddion exports.\n"
            "Disabled for text exports, which contain a single plane.")
        self.direction_combo.currentTextChanged.connect(self._reload)
        layout.addWidget(self.direction_combo)

        self.info_label = QLabel("—")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("font-size: 11px;")
        layout.addWidget(self.info_label)
        return group

    def _build_image_group(self):
        group = QGroupBox("Image flattening (2-D)")
        layout = QVBoxLayout(group)

        layout.addWidget(QLabel("Method:"))
        self.image_method_combo = QComboBox()
        self.image_method_combo.addItems(list(VALID_IMAGE_FLATTEN_METHODS))
        self.image_method_combo.setCurrentText(self._defaults.image_flatten_method)
        self.image_method_combo.setToolTip(
            "Applied to the whole image, before rows are averaged.\n\n"
            "These correct scan lines relative to one another and make the\n"
            "image readable. They do NOT change the measured blaze angle:\n"
            "the profile flattening below removes any constant or linear\n"
            "term again. Measured at 0.0000° across all samples.")
        self.image_method_combo.currentTextChanged.connect(self._reflatten)
        layout.addWidget(self.image_method_combo)

        self.image_note = QLabel("—")
        self.image_note.setWordWrap(True)
        self.image_note.setStyleSheet("font-size: 11px; color: gray;")
        layout.addWidget(self.image_note)
        return group

    def _build_profile_group(self):
        group = QGroupBox("Profile flattening (1-D)")
        layout = QVBoxLayout(group)

        caption = QLabel("Applied per row group, after averaging.\n"
                         "This one changes the answer.")
        caption.setWordWrap(True)
        caption.setStyleSheet("font-size: 11px; color: gray;")
        layout.addWidget(caption)

        layout.addWidget(QLabel("Method:"))
        self.profile_method_combo = QComboBox()
        self.profile_method_combo.addItems(list(VALID_FLATTEN_METHODS))
        self.profile_method_combo.setCurrentText(self._defaults.flatten_method)
        self.profile_method_combo.setToolTip(
            "Across the four methods the mean blaze angle varies by about\n"
            "0.49°, which is comparable to the differences between samples.\n"
            "The preview draws the background each one would remove.")
        self.profile_method_combo.currentTextChanged.connect(self._redraw)
        layout.addWidget(self.profile_method_combo)

        row = QHBoxLayout()
        row.addWidget(QLabel("Poly order:"))
        self.poly_order_spin = QSpinBox()
        self.poly_order_spin.setRange(1, 5)
        self.poly_order_spin.setValue(self._defaults.flatten_poly_order)
        self.poly_order_spin.valueChanged.connect(self._redraw)
        row.addWidget(self.poly_order_spin)
        layout.addLayout(row)

        layout.addWidget(QLabel("Feature (level_grooves):"))
        self.feature_combo = QComboBox()
        self.feature_combo.addItems(list(VALID_FLATTEN_FEATURES))
        self.feature_combo.setCurrentText(self._defaults.flatten_feature)
        self.feature_combo.setToolTip(
            "Which features to level on: 'peaks' (the lands between grooves),\n"
            "'troughs' (groove bottoms), or 'both'.")
        self.feature_combo.currentTextChanged.connect(self._redraw)
        layout.addWidget(self.feature_combo)

        layout.addWidget(QLabel("Exclude edges (fraction):"))
        self.exclude_spin = QDoubleSpinBox()
        self.exclude_spin.setRange(0.0, 0.45)
        self.exclude_spin.setSingleStep(0.05)
        self.exclude_spin.setDecimals(2)
        self.exclude_spin.setValue(self._defaults.flatten_exclude_edges)
        self.exclude_spin.valueChanged.connect(self._redraw)
        layout.addWidget(self.exclude_spin)
        return group

    # ── State the window reads ───────────────────────────────────────────────

    @property
    def data(self):
        """The flattened image the analysis should use, or None."""
        return self._flat

    @property
    def raw_data_array(self):
        """The image as loaded, before flattening. Kept for comparison."""
        return self._raw

    @property
    def scan_size(self):
        return self._scan_size

    @property
    def filename(self):
        return self._filename

    def settings_overrides(self) -> dict:
        """The settings this tab owns, for the window to merge into a run."""
        return dict(
            scan_x_size=self.scan_size_spin.value(),
            spm_direction=self.direction_combo.currentText(),
            image_flatten_method=self.image_method_combo.currentText(),
            flatten_method=self.profile_method_combo.currentText(),
            flatten_poly_order=self.poly_order_spin.value(),
            flatten_feature=self.feature_combo.currentText(),
            flatten_exclude_edges=self.exclude_spin.value(),
        )

    def set_period_est(self, period_nm: float) -> None:
        """
        The period lives with the analysis controls, but `level_grooves` needs it
        to find the features it levels on. The window pushes it here so the
        preview matches what the analysis will actually do.
        """
        if period_nm and period_nm != self._period_est:
            self._period_est = period_nm
            self._redraw()

    # ── Loading ──────────────────────────────────────────────────────────────

    def _browse(self):
        start = os.path.join(PROJECT_ROOT, 'data')
        if not os.path.isdir(start):
            start = ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open AFM File", start,
            "AFM data (*.txt *.dat *.asc *.spm);;"
            "Nanoscope (*.spm);;Text export (*.txt);;All files (*)")
        if path:
            self.load(path)

    def _reload(self):
        if self._filename:
            self.load(self._filename)

    def load(self, path):
        """Read a file, flatten it, and tell the window."""
        try:
            from ...io.spm import is_nanoscope_file
            nanoscope = is_nanoscope_file(path)

            self.direction_combo.blockSignals(True)
            self.direction_combo.setEnabled(nanoscope)
            self.direction_combo.blockSignals(False)

            settings = self._defaults.with_(**self.settings_overrides())
            with contextlib.redirect_stdout(io.StringIO()):
                data, scan_size = load_afm_data(
                    path, default_scan_size=self.scan_size_spin.value(),
                    settings=settings)

            self._raw = data
            self._scan_size = scan_size
            self._filename = path
            self.scan_size_spin.setValue(scan_size)

            self.file_label.setText(os.path.basename(path))
            self.file_label.setStyleSheet("font-size: 11px;")

            rows, cols = data.shape
            source = (f"Nanoscope: {settings.spm_channel} / "
                      f"{settings.spm_direction}" if nanoscope else "Text export")
            self.info_label.setText(
                f"{source}\nShape: {rows} × {cols} px\n"
                f"Scan width: {scan_size:.3f} µm")

            self._reflatten()

        except Exception as exc:
            self.info_label.setText(f"Error:\n{exc}")
            self._raw = self._flat = None
            self.dataChanged.emit()

    def _reflatten(self):
        """Re-apply image flattening and redraw. Cheap; no thread needed."""
        if self._raw is None:
            return
        method = self.image_method_combo.currentText()
        before = row_offset_spread(self._raw)
        self._flat = flatten_image(self._raw, method)
        after = row_offset_spread(self._flat)
        self.image_note.setText(
            f"Row-offset spread\n  {before*1e9:.2f} → {after*1e9:.2f} nm")
        self._redraw()
        self.dataChanged.emit()

    # ── Drawing ──────────────────────────────────────────────────────────────

    def _show_placeholder(self, message="Load an AFM file to begin"):
        self.figure.clear()
        axes = self.figure.add_subplot(111)
        axes.text(0.5, 0.5, message, ha='center', va='center', fontsize=12,
                  color='gray', transform=axes.transAxes)
        axes.set_axis_off()
        self.canvas.draw()

    def _redraw(self):
        if self._flat is None:
            self._show_placeholder()
            return

        QApplication.processEvents()
        self.figure.clear()
        image_ax = self.figure.add_subplot(2, 1, 1)
        profile_ax = self.figure.add_subplot(2, 1, 2)

        image = image_ax.imshow(
            self._flat * 1e9, aspect='auto', cmap='viridis',
            extent=[0, self._scan_size, 0, self._flat.shape[0]], origin='upper')
        self.figure.colorbar(image, ax=image_ax, label='Height (nm)')
        image_ax.set_xlabel("X (µm)", fontsize=10)
        image_ax.set_ylabel("Scan line", fontsize=10)
        image_ax.set_title(
            f"{os.path.basename(self._filename)} — image flattening: "
            f"{self.image_method_combo.currentText()}", fontsize=10)

        x_um, profile_nm = raw_data(self._flat, self._scan_size)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                flat, background = flatten_profile(
                    x_um, profile_nm,
                    method=self.profile_method_combo.currentText(),
                    poly_order=self.poly_order_spin.value(),
                    exclude_edges=self.exclude_spin.value(),
                    feature=self.feature_combo.currentText(),
                    period_nm=self._period_est)
        except Exception as exc:
            profile_ax.text(0.5, 0.5, f"Profile flattening failed:\n{exc}",
                            ha='center', va='center', color='crimson',
                            transform=profile_ax.transAxes)
            profile_ax.set_axis_off()
            self.canvas.draw()
            return

        profile_ax.plot(x_um, profile_nm, color='0.6', linewidth=0.9,
                        label='Averaged profile')
        profile_ax.plot(x_um, background, color='crimson', linewidth=1.4,
                        linestyle='--', label='Background to remove')
        profile_ax.plot(x_um, flat, color='steelblue', linewidth=1.0,
                        label='After flattening')
        profile_ax.set_xlabel("Position (µm)", fontsize=10)
        profile_ax.set_ylabel("Height (nm)", fontsize=10)
        profile_ax.set_title(
            f"Profile flattening: {self.profile_method_combo.currentText()}"
            f"  (removes {np.ptp(background):.1f} nm of background)", fontsize=10)
        profile_ax.legend(fontsize=8)
        profile_ax.grid(True, alpha=0.3)

        self.canvas.draw()
