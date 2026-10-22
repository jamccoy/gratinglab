"""
The Boundary tab: build and export a PCGrate .ggp profile.

Takes the scan the Import tab loaded, averages its grooves into one, and shows
what will be written before writing it.

The preview is the point. A boundary profile feeds a grating-efficiency model,
and the two ways it goes wrong - one groove dragging the average, or a step at
the period boundary where the profile fails to tile - are obvious in a plot and
invisible in the file. The numbers alone will not tell you.

Computation lives in ``gratinglab.metrology.boundary.build_boundary_profile``, which the
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
from ...core.tip import format_angles
from ...config import RESULTS_DIR
from gratinglab.io.boundary import write_boundary_csv, write_boundary_json
from gratinglab.io.ggp import write_ggp

from ...io.metrics import write_profile_metrics
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
            "Grooves sampled at or below this many points per half period are\n"
            "skipped. A groove the scan edge cuts into is always skipped,\n"
            "whatever this says: it would narrow the window every groove\n"
            "shares, and normalisation would then stretch the profile.")
        self.min_half_spin.valueChanged.connect(self.refresh)
        layout.addWidget(self.min_half_spin)

        self.export_btn = QPushButton("Export…")
        self.export_btn.setEnabled(False)
        self.export_btn.setToolTip(
            "PCGrate .ggp (shape only), gratinglab .json (shape, measured\n"
            "period and provenance, reloads into the efficiency window), or\n"
            ".csv (columns in period fractions and nm, for other software).")
        self.export_btn.clicked.connect(self.export)
        layout.addWidget(self.export_btn)

        # The same profile, handed to a solver instead of to a file. Enabled on
        # the same condition, because both need exactly one thing: a profile.
        self.efficiency_btn = QPushButton("Efficiency…")
        self.efficiency_btn.setEnabled(False)
        self.efficiency_btn.setToolTip(
            "Diffraction efficiency of this groove, with the period this scan\n"
            "measured. Exporting a .ggp and loading it elsewhere loses that\n"
            "number -- the file format cannot carry it.")
        self.efficiency_btn.clicked.connect(self.show_efficiency)
        layout.addWidget(self.efficiency_btn)

        # The full efficiency application rather than the preview above: every
        # solver, mount and convergence check. Handed over as a gratinglab
        # .json, so the measured period arrives with the shape.
        self.send_btn = QPushButton("Send to efficiency window")
        self.send_btn.setEnabled(False)
        self.send_btn.setToolTip(
            "Open this groove in the main gratinglab window, with the period\n"
            "this scan measured already filled in.")
        self.send_btn.clicked.connect(self.send_to_efficiency)
        layout.addWidget(self.send_btn)

        self.swaths_btn = QPushButton("Swath statistics…")
        self.swaths_btn.setEnabled(False)
        self.swaths_btn.setToolTip(
            "Build the groove separately for each horizontal band of rows and "
            "show how much the shape varies down the scan.")
        self.swaths_btn.clicked.connect(self.show_swaths)
        layout.addWidget(self.swaths_btn)
        self._efficiency_window = None
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
            self._set_actions_enabled(False)
            self._show_placeholder()
            return

        settings = self._settings.with_(**self.settings_overrides())
        try:
            self._profile = build_boundary_profile(
                self._data, self._scan_size, settings)
        except Exception as exc:
            self._profile = None
            self._set_actions_enabled(False)
            self.metrics_label.setText(f"Could not build a profile:\n{exc}")
            self._show_placeholder("No boundary profile")
            return

        self._set_actions_enabled(True)
        self._update_metrics()
        self._draw()

    def _update_metrics(self):
        p = self._profile
        m = p.metrics
        edge = (f"\n{'':<18}({p.n_edge_rejected} rejected at a scan edge)"
                if p.n_edge_rejected else "")
        tip = ""
        if m.get('tip_correction', 'none') != 'none':
            # The uncertain remainder is an upper bound, never a measurement,
            # so the fraction recovered belongs beside the depth it qualifies.
            tip = (f"\nTip correction  : {m['tip_correction']}, "
                   f"R {m.get('tip_radius_nm', float('nan')):g} nm, "
                   f"{format_angles(m.get('tip_half_angle_deg', float('nan')), m.get('tip_half_angle_right_deg'))}"
                   f"\n{'':<18}{100 * m.get('tip_certain_fraction', float('nan')):.1f}% "
                   f"of pixels recovered exactly")
        self.metrics_label.setText(
            f"Period          : {p.period_nm:.2f} nm\n"
            f"Grooves averaged: {p.n_used} of {p.n_grooves}{edge}\n"
            f"Groove depth    : {m['groove_depth']:.4f} of period\n"
            f"Peak-to-valley  : {m['peak_to_valley']:.4f} of period\n"
            f"RMS slope       : {m['rms_slope']:.4f}\n"
            f"Max sidewall    : {m['max_angle_deg']:.2f} deg\n"
            f"Max curvature   : {m['max_curvature']:.4f}{tip}")

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
            f"Averaged groove: {os.path.basename(self._filename)}", fontsize=10)
        avg_ax.legend(fontsize=8)
        avg_ax.grid(True, alpha=0.3)

        # What lands in the file.
        norm_ax.plot(p.x_norm, p.y_norm, color='steelblue', linewidth=1.6)
        norm_ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
        norm_ax.set_xlabel("Normalised position (0 to 1 = one period)", fontsize=10)
        norm_ax.set_ylabel("Height (fraction of period)", fontsize=10)
        norm_ax.set_title(
            f"Exported profile: {len(p.x_norm)} points, "
            f"depth {p.metrics['groove_depth']:.4f} of period", fontsize=10)
        norm_ax.grid(True, alpha=0.3)

        self.canvas.draw()

    def _set_actions_enabled(self, enabled: bool) -> None:
        """Both actions need one thing -- a profile -- so they move together."""
        self.export_btn.setEnabled(enabled)
        self.efficiency_btn.setEnabled(enabled)
        self.send_btn.setEnabled(enabled)
        self.swaths_btn.setEnabled(enabled)

    # ── Out to a solver ──────────────────────────────────────────────────────

    def show_efficiency(self):
        """Hand the previewed profile straight to the scalar solver.

        No file in between, because a .ggp round trip keeps
        the shape and drops the period, which then has to be retyped from the
        metrics panel above.
        """
        if self._profile is None:
            return
        # Imported here rather than at module scope so that a Boundary tab that
        # is only ever used to export keeps loading the solver stack off its
        # start-up path.
        from .efficiency_dialog import EfficiencyDialog

        EfficiencyDialog(self._profile, self).exec()

    def swath_dialog(self):
        """A SwathDialog for the current scan and settings, or None."""
        if self._profile is None:
            return None
        from .swath_dialog import SwathDialog

        settings = self._settings.with_(**self.settings_overrides())
        return SwathDialog(self._data, self._scan_size, settings, self._filename, self)

    def show_swaths(self):
        """How much the groove varies from swath to swath down the scan."""
        dialog = self.swath_dialog()
        if dialog is not None:
            dialog.exec()

    # ── Export ───────────────────────────────────────────────────────────────

    #: File-dialog filters, and the suffix each one writes.
    EXPORT_FORMATS = (
        ("PCGrate boundary profile (*.ggp)", ".ggp"),
        ("gratinglab boundary, with period and provenance (*.json)", ".json"),
        ("CSV table, fractions of the period and nm (*.csv)", ".csv"),
    )

    def export(self):
        """Write the profile as .ggp, .json or .csv, by the chosen file type.

        A .ggp or .csv also gets the metrics file beside it; the JSON carries
        the metrics inside itself.
        """
        if self._profile is None:
            return

        stem = os.path.splitext(os.path.basename(self._filename))[0]
        suggested = os.path.join(RESULTS_DIR,
                                 f'averaged_groove_profile_{stem}.ggp')
        os.makedirs(RESULTS_DIR, exist_ok=True)

        filters = ";;".join(f for f, _ in self.EXPORT_FORMATS) + ";;All files (*)"
        path, chosen = QFileDialog.getSaveFileName(
            self, "Export boundary profile", suggested, filters)
        if not path:
            return

        try:
            path, metrics_path = self.write_profile(path, chosen)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
            return

        message = f"Boundary profile written to:\n{path}"
        if metrics_path:
            message += f"\n\nMetrics written to:\n{metrics_path}"
        QMessageBox.information(self, "Exported", message)

    def write_profile(self, path, chosen_filter=""):
        """Write the current profile; return (path written, metrics path or None).

        The format follows the file's suffix, or the dialog's filter when the
        name has no recognised suffix, in which case that suffix is added.
        Separate from :meth:`export` so it can be driven without a dialog.
        """
        suffix = os.path.splitext(path)[1].lower()
        known = {sfx for _, sfx in self.EXPORT_FORMATS}
        if suffix not in known:
            suffix = next((sfx for f, sfx in self.EXPORT_FORMATS if f == chosen_filter), ".ggp")
            path += suffix

        settings = self._settings.with_(**self.settings_overrides()) if self._settings else None
        metrics_path = None
        if suffix == ".json":
            write_boundary_json(path, self._profile, source=self._filename, settings=settings)
        elif suffix == ".csv":
            write_boundary_csv(path, self._profile, source=self._filename)
            metrics_path = os.path.splitext(path)[0] + '_metrics.txt'
            write_profile_metrics(metrics_path, self._profile.metrics)
        else:
            write_ggp(path, t=self._profile.x_norm, y=self._profile.y_norm)
            metrics_path = os.path.splitext(path)[0] + '_metrics.txt'
            write_profile_metrics(metrics_path, self._profile.metrics)
        return path, metrics_path

    def _forget_efficiency_window(self, *_):
        self._efficiency_window = None

    def send_to_efficiency(self):
        """Open this groove in the main efficiency window, period included.

        Written to a gratinglab .json in the temporary directory first, because
        the efficiency window's form takes a profile *file*. That also leaves
        the user a file the period travels with, should they want to keep it.
        The same window is reused on a second send rather than stacking copies.
        """
        if self._profile is None:
            return
        try:
            from gratinglab.gui.qt.main_window import MainWindow
        except ImportError as exc:
            QMessageBox.warning(
                self, "Efficiency window unavailable",
                f"The efficiency window could not be loaded:\n{exc}\n\n"
                'It needs the gui extra: pip install "gratinglab[gui]"')
            return

        import tempfile

        stem = os.path.splitext(os.path.basename(self._filename or "scan"))[0]
        folder = os.path.join(tempfile.gettempdir(), "gratinglab")
        os.makedirs(folder, exist_ok=True)
        path, _ = self.write_profile(os.path.join(folder, f"{stem}_boundary.json"))

        if self._efficiency_window is None:
            self._efficiency_window = MainWindow()
            # Forget it once Qt destroys it, so a later send opens a fresh one
            # rather than touching a deleted C++ object, and no dangling
            # wrapper outlives the window into interpreter shutdown.
            self._efficiency_window.destroyed.connect(self._forget_efficiency_window)
        window = self._efficiency_window
        window.open_boundary(path)
        window.show()
        window.raise_()
        window.activateWindow()
        return window
