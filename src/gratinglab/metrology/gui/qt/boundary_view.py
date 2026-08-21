"""
The Boundary tab: build and export a PCGrate .ggp profile.

Takes the scan the Import tab loaded, averages its grooves into one, and shows
what will be written before writing it.

The preview is the point. A boundary profile feeds a grating-efficiency model,
and the two ways it goes wrong - one groove dragging the average, or a step at
the period boundary where the profile fails to tile - are obvious in a plot and
invisible in the file. The numbers alone will not tell you.

Computation lives in ``afm_analysis.boundary.build_boundary_profile``, which the
CLI mode also uses, so the panel cannot drift from ``ANALYSIS_MODE = 'ggp'``.
"""
from __future__ import annotations

import os

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtWidgets import (
    QCheckBox, QFileDialog, QFrame, QGroupBox, QHBoxLayout, QLabel,
    QMessageBox, QPushButton, QScrollArea, QSpinBox, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ...boundary import build_boundary_profile
from ...config import RESULTS_DIR
from ...io.ggp import write_ggp, write_profile_metrics
from ...settings import AnalysisSettings

__all__ = ["BoundaryView"]


class BoundaryView(QWidget):
    """Average the grooves of the loaded scan into a PCGrate boundary profile."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self._defaults = AnalysisSettings.from_config()
        self._data = None
        self._scan_size = None
        self._filename = None
        self._settings = None
        self._profile = None

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
        column.addWidget(self._build_export_group())
        column.addWidget(self._build_metrics_group())
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

    def _build_export_group(self):
        group = QGroupBox("Boundary profile")
        layout = QVBoxLayout(group)

        caption = QLabel("Grooves averaged into one, normalised to a single "
                         "period for PCGrate.")
        caption.setWordWrap(True)
        caption.setStyleSheet("font-size: 11px; color: gray;")
        layout.addWidget(caption)

        layout.addWidget(QLabel("Points in the exported profile:"))
        self.n_points_spin = QSpinBox()
        self.n_points_spin.setRange(50, 20000)
        self.n_points_spin.setSingleStep(100)
        self.n_points_spin.setValue(self._defaults.ggp_n_points)
        self.n_points_spin.setToolTip(
            "How finely the averaged groove is resampled.\n"
            "This is the number of x y pairs in the .ggp file.")
        self.n_points_spin.valueChanged.connect(self.refresh)
        layout.addWidget(self.n_points_spin)

        self.smoothing_check = QCheckBox("Smooth the profile")
        self.smoothing_check.setChecked(self._defaults.ggp_apply_smoothing)
        self.smoothing_check.setToolTip(
            "Light smoothing that wraps at the period boundary, to remove\n"
            "kinks left by interpolating each groove onto a common axis.")
        self.smoothing_check.stateChanged.connect(self.refresh)
        layout.addWidget(self.smoothing_check)

        row = QHBoxLayout()
        row.addWidget(QLabel("Window:"))
        self.smoothing_spin = QSpinBox()
        self.smoothing_spin.setRange(2, 101)
        self.smoothing_spin.setValue(self._defaults.ggp_smoothing_window)
        self.smoothing_spin.valueChanged.connect(self.refresh)
        row.addWidget(self.smoothing_spin)
        layout.addLayout(row)

        layout.addWidget(QLabel("Minimum half-width (samples):"))
        self.min_half_spin = QSpinBox()
        self.min_half_spin.setRange(1, 200)
        self.min_half_spin.setValue(self._defaults.ggp_min_half_width)
        self.min_half_spin.setToolTip(
            "Grooves whose symmetric extent is at or below this are skipped.\n"
            "A groove near a scan edge would otherwise narrow the window every\n"
            "groove shares, and normalisation would then stretch the exported\n"
            "profile to fill a period.")
        self.min_half_spin.valueChanged.connect(self.refresh)
        layout.addWidget(self.min_half_spin)

        self.export_btn = QPushButton("Export .ggp…")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self.export)
        layout.addWidget(self.export_btn)
        return group

    def _build_metrics_group(self):
        group = QGroupBox("Profile metrics")
        layout = QVBoxLayout(group)
        self.metrics_label = QLabel("Load a scan in the Import tab.")
        self.metrics_label.setWordWrap(True)
        self.metrics_label.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(self.metrics_label)
        return group

    # ── Data in ──────────────────────────────────────────────────────────────

    def set_scan(self, data, scan_size, filename, settings):
        """Take the current scan from the window. Called on every Import change."""
        self._data = data
        self._scan_size = scan_size
        self._filename = filename
        self._settings = settings
        self.refresh()

    def settings_overrides(self) -> dict:
        """The export settings this tab owns."""
        return dict(
            ggp_n_points=self.n_points_spin.value(),
            ggp_apply_smoothing=self.smoothing_check.isChecked(),
            ggp_smoothing_window=self.smoothing_spin.value(),
            ggp_min_half_width=self.min_half_spin.value(),
        )

    @property
    def profile(self):
        """The current BoundaryProfile, or None."""
        return self._profile

    # ── Compute and draw ─────────────────────────────────────────────────────

    def refresh(self):
        """
        Rebuild the profile and redraw.

        No worker thread: averaging a few grooves is well under a second, unlike
        the blaze-angle analysis which has one.
        """
        if self._data is None or self._settings is None:
            self._profile = None
            self.export_btn.setEnabled(False)
            self._show_placeholder()
            return

        settings = self._settings.with_(**self.settings_overrides())
        try:
            self._profile = build_boundary_profile(
                self._data, self._scan_size, settings)
        except Exception as exc:
            self._profile = None
            self.export_btn.setEnabled(False)
            self.metrics_label.setText(f"Could not build a profile:\n{exc}")
            self._show_placeholder("No boundary profile")
            return

        self.export_btn.setEnabled(True)
        self._update_metrics()
        self._draw()

    def _update_metrics(self):
        p = self._profile
        m = p.metrics
        edge = (f"\n{'':<18}({p.n_edge_rejected} rejected at a scan edge)"
                if p.n_edge_rejected else "")
        self.metrics_label.setText(
            f"Period          : {p.period_nm:.2f} nm\n"
            f"Grooves averaged: {p.n_used} of {p.n_grooves}{edge}\n"
            f"Groove depth    : {m['groove_depth']:.4f} of period\n"
            f"Peak-to-valley  : {m['peak_to_valley']:.4f} of period\n"
            f"RMS slope       : {m['rms_slope']:.4f}\n"
            f"Max sidewall    : {m['max_angle_deg']:.2f} deg\n"
            f"Max curvature   : {m['max_curvature']:.4f}")

    def _show_placeholder(self, message="Load an AFM file in the Import tab"):
        self.figure.clear()
        axes = self.figure.add_subplot(111)
        axes.text(0.5, 0.5, message, ha='center', va='center', fontsize=12,
                  color='gray', transform=axes.transAxes)
        axes.set_axis_off()
        self.canvas.draw()

    def _draw(self):
        p = self._profile
        self.figure.clear()
        avg_ax = self.figure.add_subplot(2, 1, 1)
        norm_ax = self.figure.add_subplot(2, 1, 2)

        # Averaged groove with the spread across the grooves that made it. A
        # wide band means the grooves disagree - either real variation, or one
        # groove dragging the average.
        x_nm = p.x_avg_um * 1000
        avg_ax.plot(x_nm, p.y_avg_nm, color='crimson', linewidth=1.8,
                    label=f'Average of {p.n_used} grooves')
        avg_ax.fill_between(x_nm, p.y_avg_nm - p.y_std_nm,
                            p.y_avg_nm + p.y_std_nm,
                            color='crimson', alpha=0.2, label='±1σ')
        avg_ax.set_xlabel("Position (nm)", fontsize=10)
        avg_ax.set_ylabel("Height (nm)", fontsize=10)
        avg_ax.set_title(
            f"Averaged groove — {os.path.basename(self._filename)}", fontsize=10)
        avg_ax.legend(fontsize=8)
        avg_ax.grid(True, alpha=0.3)

        # What actually lands in the file.
        norm_ax.plot(p.x_norm, p.y_norm, color='steelblue', linewidth=1.6)
        norm_ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
        norm_ax.set_xlabel("Normalised position (0 to 1 = one period)", fontsize=10)
        norm_ax.set_ylabel("Height (fraction of period)", fontsize=10)
        norm_ax.set_title(
            f"Exported profile — {len(p.x_norm)} points, "
            f"depth {p.metrics['groove_depth']:.4f} of period", fontsize=10)
        norm_ax.grid(True, alpha=0.3)

        self.canvas.draw()

    # ── Export ───────────────────────────────────────────────────────────────

    def export(self):
        """Write the .ggp, and the metrics file beside it."""
        if self._profile is None:
            return

        stem = os.path.splitext(os.path.basename(self._filename))[0]
        suggested = os.path.join(RESULTS_DIR,
                                 f'averaged_groove_profile_{stem}.ggp')
        os.makedirs(RESULTS_DIR, exist_ok=True)

        path, _ = QFileDialog.getSaveFileName(
            self, "Export PCGrate boundary profile", suggested,
            "PCGrate boundary profile (*.ggp);;All files (*)")
        if not path:
            return

        try:
            write_ggp(path, self._profile.x_norm, self._profile.y_norm)
            metrics_path = os.path.splitext(path)[0] + '_metrics.txt'
            write_profile_metrics(metrics_path, self._profile.metrics)
        except OSError as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
            return

        QMessageBox.information(
            self, "Exported",
            f"Boundary profile written to:\n{path}\n\n"
            f"Metrics written to:\n{metrics_path}")
