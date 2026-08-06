"""
AFM Blaze Angle Analysis - GUI
v0.2: file viewing plus interactive blaze-angle analysis

Run:
    .venv/bin/python afm_gui.py

The analysis parameters exposed here are the ones that actually move the answer.
FACET_TRIM in particular is worth exploring: on the master sample, changing it
from 0.10 to 0.25 shifts the mean blaze angle by over 3 degrees, because a larger
trim keeps the fit away from the rounded groove top. Being able to see that
interactively is the point of this window.
"""

import contextlib
import io
import os
import sys

import numpy as np

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QDoubleSpinBox, QSpinBox, QGroupBox,
    QStatusBar, QFrame, QComboBox, QCheckBox, QScrollArea
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_DIR)

# Imported eagerly and without a fallback. Earlier versions carried private
# copies of the loader in case the package could not be imported; those copies
# then drifted from the real ones. The package is right next to this file, so a
# failure here is a broken install and should say so loudly.
from afm_analysis.core.processing import (
    load_afm_data, raw_data, raw_data_multi_group
)
import afm_analysis.analyzer as analyzer
from afm_analysis.settings import AnalysisSettings, MAX_FACET_TRIM

# Control defaults come from config.py, read once at import.
_DEFAULTS = AnalysisSettings.from_config()


# ── Matplotlib canvas ─────────────────────────────────────────────────────────

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(8, 5), tight_layout=True)
        super().__init__(self.fig)
        self.setParent(parent)
        self._draw_placeholder()

    def _draw_placeholder(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.5, "Load an AFM file to see the raw profile",
                ha='center', va='center', fontsize=12, color='gray',
                transform=ax.transAxes)
        ax.set_axis_off()
        self.draw()

    def plot_raw(self, disp_um, profile_nm, filename, scan_size_um):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.plot(disp_um, profile_nm, color='steelblue', linewidth=1.0)
        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        ax.set_title(f"Raw AFM Profile — {os.path.basename(filename)}\n"
                     f"Scan width: {scan_size_um:.3f} µm", fontsize=11)
        ax.grid(True, alpha=0.3)
        self.draw()

    def plot_row_groups(self, disp_um, profiles_nm, group_info, filename, scan_size_um):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        n_groups = group_info['n_groups']
        cmap = matplotlib.colormaps['coolwarm']

        for i, profile in enumerate(profiles_nm):
            ax.plot(disp_um, profile, color=cmap(i / max(n_groups - 1, 1)),
                    linewidth=0.8, alpha=0.6)

        mean_profile = np.mean(np.array(profiles_nm), axis=0)
        ax.plot(disp_um, mean_profile, color='black', linewidth=1.8,
                label=f'Mean (n={n_groups})', zorder=5)

        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        ax.set_title(f"Row-Group Profiles — {os.path.basename(filename)}\n"
                     f"{n_groups} groups × {group_info['rows_per_group']} rows each  |  "
                     f"Scan width: {scan_size_um:.3f} µm", fontsize=10)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        sm = matplotlib.cm.ScalarMappable(
            cmap=cmap,
            norm=matplotlib.colors.Normalize(vmin=0, vmax=group_info['n_rows']))
        sm.set_array([])
        cb = self.fig.colorbar(sm, ax=ax, pad=0.02)
        cb.set_label("Image row (top → bottom)", fontsize=9)
        self.draw()

    def plot_2d(self, data, scan_size_um, filename):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        im = ax.imshow(data * 1e9, aspect='auto', cmap='viridis',
                       extent=[0, scan_size_um, 0, data.shape[0]], origin='upper')
        self.fig.colorbar(im, ax=ax, label='Height (nm)')
        ax.set_xlabel("X (µm)", fontsize=11)
        ax.set_ylabel("Scan line", fontsize=11)
        ax.set_title(f"2D AFM Topography — {os.path.basename(filename)}", fontsize=11)
        self.draw()

    def plot_groove_detection(self, result, filename, settings):
        """Flattened profile with the grooves the analysis actually used"""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        x, y = result['raw_x'], result['flat_y']
        centers = result['groove_centers']

        ax.plot(x, y, 'k-', linewidth=0.7, label='Flattened profile')
        ax.plot(x[centers], y[centers], 'ro', markersize=9,
                label=f'Detected grooves (N={len(centers)})')

        # Shade the zone excluded near each scan edge, so it is obvious why a
        # groove at the boundary was dropped.
        edge = settings.edge_exclusion_periods
        if edge > 0 and len(centers) > 0:
            period_um = result['period_nm'] / 1000.0
            margin = edge * period_um
            for lo, hi in ((x[0], x[0] + margin), (x[-1] - margin, x[-1])):
                ax.axvspan(lo, hi, color='red', alpha=0.10)
            ax.axvspan(x[0], x[0], color='red', alpha=0.10,
                       label=f'Edge exclusion ({edge}× period)')

        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        ax.set_title(f"Groove Detection — {os.path.basename(filename)}", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        self.draw()

    def plot_blaze_angles(self, result, filename):
        """Angle distribution, fit quality, and angle against facet width"""
        self.fig.clear()
        ax1, ax2, ax3 = (self.fig.add_subplot(1, 3, i) for i in (1, 2, 3))

        angles = np.asarray(result['all_angles'])
        quality = result['quality']
        mean, sem = result['mean_angle'], result.get('sem', 0.0)

        ax1.hist(angles, bins=20, edgecolor='black', alpha=0.75, color='steelblue')
        ax1.axvline(mean, color='r', linestyle='--', linewidth=2,
                    label=f'Mean {mean:.2f}°')
        ax1.set_xlabel("Blaze angle (°)", fontsize=10)
        ax1.set_ylabel("Count", fontsize=10)
        ax1.set_title(f"Distribution (N={len(angles)})\nσ={result['std_angle']:.2f}°",
                      fontsize=10)
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)

        r2 = np.array([q.get('blaze_r2', np.nan) for q in quality])
        ax2.hist(r2, bins=20, edgecolor='black', alpha=0.75, color='seagreen')
        ax2.set_xlabel("Fit R²", fontsize=10)
        ax2.set_ylabel("Count", fontsize=10)
        ax2.set_title(f"Fit quality\nworst R² = {np.nanmin(r2):.4f}", fontsize=10)
        ax2.grid(True, alpha=0.3)

        widths = np.array([q.get('blaze_width_nm', np.nan) for q in quality])
        ax3.plot(widths, angles, 'o', markersize=4, alpha=0.45, color='darkorange')
        ax3.set_xlabel("Blaze facet width (nm)", fontsize=10)
        ax3.set_ylabel("Blaze angle (°)", fontsize=10)
        ax3.set_title("Angle vs facet width\n(narrow facets = unreliable fits)",
                      fontsize=10)
        ax3.grid(True, alpha=0.3)

        self.fig.suptitle(f"{os.path.basename(filename)} — "
                          f"{mean:.2f}° ± {sem:.3f}° (SEM)", fontsize=11)
        self.draw()


# ── Main window ───────────────────────────────────────────────────────────────

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AFM Blaze Angle Analysis")
        self.resize(1280, 760)

        self._data = None
        self._disp_um = None
        self._profile_nm = None
        self._scan_size = None
        self._filename = None
        self._group_profiles = None
        self._group_info = None
        self._result = None
        self._settings = None

        self._build_ui()

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
        left_layout.addWidget(self._build_scan_group())
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
        self.info_label = QLabel("—")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("font-size: 11px;")
        layout.addWidget(self.info_label)
        return group

    def _build_scan_group(self):
        group = QGroupBox("Scan")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Default scan width (µm):"))
        self.scan_size_spin = QDoubleSpinBox()
        self.scan_size_spin.setRange(0.1, 1000.0)
        self.scan_size_spin.setValue(2.0)
        self.scan_size_spin.setSingleStep(0.1)
        self.scan_size_spin.setDecimals(3)
        self.scan_size_spin.setToolTip(
            "Used only when the scan width cannot be read from the file header.")
        layout.addWidget(self.scan_size_spin)
        reload_btn = QPushButton("Reload with these settings")
        reload_btn.clicked.connect(self._reload_file)
        layout.addWidget(reload_btn)
        return group

    def _build_analysis_group(self):
        group = QGroupBox("Analysis parameters")
        layout = QVBoxLayout(group)

        layout.addWidget(QLabel("Estimated period (nm):"))
        self.period_spin = QDoubleSpinBox()
        self.period_spin.setRange(1.0, 100000.0)
        self.period_spin.setDecimals(2)
        self.period_spin.setValue(_DEFAULTS.period_est)
        self.period_spin.setToolTip("Groove spacing. Must match your grating.")
        layout.addWidget(self.period_spin)

        layout.addWidget(QLabel("Facet trim (fraction):"))
        self.trim_spin = QDoubleSpinBox()
        # Capped at 0.28, not 0.5: the blaze facet is trimmed 2.5x harder on the
        # trough side, so trim x 3.5 is removed in total and anything above
        # ~0.286 empties the facet and the analysis returns nothing.
        self.trim_spin.setRange(0.0, MAX_FACET_TRIM)
        self.trim_spin.setSingleStep(0.05)
        self.trim_spin.setDecimals(2)
        self.trim_spin.setValue(_DEFAULTS.facet_trim)
        self.trim_spin.setToolTip(
            "Fraction trimmed from each end of the facet before fitting.\n"
            "Too small risks measuring the rounded groove top and reading low.\n"
            "This parameter moves the answer more than any other.\n"
            "Capped at 0.28 — above that the trim consumes the whole facet.")
        layout.addWidget(self.trim_spin)

        layout.addWidget(QLabel("Blaze side:"))
        self.side_combo = QComboBox()
        self.side_combo.addItems(['negative_slope', 'positive_slope', 'longer'])
        self.side_combo.setCurrentText(_DEFAULTS.blaze_side)
        self.side_combo.setToolTip(
            "Which facet to measure, chosen by slope sign rather than position.")
        layout.addWidget(self.side_combo)

        layout.addWidget(QLabel("Edge exclusion (periods):"))
        self.edge_spin = QDoubleSpinBox()
        self.edge_spin.setRange(0.0, 3.0)
        self.edge_spin.setSingleStep(0.1)
        self.edge_spin.setDecimals(2)
        self.edge_spin.setValue(_DEFAULTS.edge_exclusion_periods)
        self.edge_spin.setToolTip(
            "Reject grooves this close to either end of the scan line.\n"
            "Their facet is clipped by the edge, so the fitted angle is\n"
            "meaningless. 0 disables the check.")
        layout.addWidget(self.edge_spin)

        self.row_groups_check = QCheckBox("Use row groups")
        self.row_groups_check.setChecked(_DEFAULTS.use_row_groups)
        self.row_groups_check.setToolTip(
            "Analyse N horizontal bands separately instead of averaging the\n"
            "whole image into one profile. Many more measurements, but they\n"
            "re-measure the same grooves and are not independent.")
        layout.addWidget(self.row_groups_check)

        n_layout = QHBoxLayout()
        n_layout.addWidget(QLabel("N groups:"))
        self.n_groups_spin = QSpinBox()
        self.n_groups_spin.setRange(2, 200)
        self.n_groups_spin.setValue(_DEFAULTS.n_row_groups)
        n_layout.addWidget(self.n_groups_spin)
        layout.addLayout(n_layout)

        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setEnabled(False)
        font = QFont()
        font.setBold(True)
        self.run_btn.setFont(font)
        self.run_btn.clicked.connect(self._run_analysis)
        layout.addWidget(self.run_btn)

        return group

    def _build_view_group(self):
        group = QGroupBox("View")
        layout = QVBoxLayout(group)

        self.raw_profile_btn = QPushButton("Raw Profile (1D)")
        self.raw_profile_btn.clicked.connect(self._show_raw_profile)
        self.topography_btn = QPushButton("2D Topography")
        self.topography_btn.clicked.connect(self._show_2d)
        self.row_groups_btn = QPushButton("Row Groups")
        self.row_groups_btn.clicked.connect(self._show_row_groups)
        self.detection_btn = QPushButton("Groove Detection")
        self.detection_btn.clicked.connect(self._show_detection)
        self.angles_btn = QPushButton("Blaze Angles")
        self.angles_btn.clicked.connect(self._show_angles)

        self._file_buttons = [self.raw_profile_btn, self.topography_btn,
                              self.row_groups_btn]
        self._result_buttons = [self.detection_btn, self.angles_btn]
        for b in self._file_buttons + self._result_buttons:
            b.setEnabled(False)
            layout.addWidget(b)
        return group

    def _build_results_panel(self):
        self.results_group = QGroupBox("Results")
        layout = QVBoxLayout(self.results_group)
        self.results_label = QLabel("Load a file and press Run Analysis.")
        self.results_label.setStyleSheet("font-family: monospace; font-size: 11px;")
        self.results_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.results_label)
        self.results_group.setMaximumHeight(130)
        return self.results_group

    # ── Actions ──────────────────────────────────────────────────────────────

    def _browse_file(self):
        start_dir = os.path.join(PROJECT_DIR, 'data')
        if not os.path.isdir(start_dir):
            start_dir = ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open AFM File", start_dir,
            "AFM Data Files (*.txt *.dat *.asc);;All Files (*)")
        if path:
            self._load(path)

    def _reload_file(self):
        if self._filename:
            self._load(self._filename)

    def _load(self, path):
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
            self._group_profiles = self._group_info = None
            self._result = None

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

            max_groups = max(2, data.shape[0] // 3)
            self.n_groups_spin.setMaximum(max_groups)
            if self.n_groups_spin.value() > max_groups:
                self.n_groups_spin.setValue(max_groups)

            self._show_raw_profile()
            self._set_status(f"Loaded: {os.path.basename(path)}")

        except Exception as exc:
            self._set_status(f"Error loading file: {exc}")
            self.info_label.setText(f"Error:\n{exc}")

    def _current_settings(self):
        """
        Build an AnalysisSettings from the controls.

        This used to rebind module-level names inside analyzer, which worked only
        because analyzer bound its configuration at import. Settings are now a
        value that gets passed in, so the window mutates nothing.
        """
        return AnalysisSettings.from_config().with_(
            period_est=self.period_spin.value(),
            facet_trim=self.trim_spin.value(),
            blaze_side=self.side_combo.currentText(),
            edge_exclusion_periods=self.edge_spin.value(),
            use_row_groups=self.row_groups_check.isChecked(),
            n_row_groups=self.n_groups_spin.value(),
        )

    def _run_analysis(self):
        if self._filename is None:
            return
        try:
            self._set_status("Running analysis…")
            QApplication.processEvents()
            settings = self._current_settings()

            log = io.StringIO()
            with contextlib.redirect_stdout(log):
                result = analyzer.analyze_single_file(
                    self._filename, show_plots=False, settings=settings)

            if result is None:
                # Clear the previous result too, or the panel and the view
                # buttons keep describing a run that is no longer current.
                self._result = None
                self._settings = None
                self._set_status("Analysis produced no measurements.")
                self.results_label.setText(
                    "No blaze angles could be extracted.\n"
                    "Facet trim above ~0.28 removes the whole facet — try lowering it,\n"
                    "or check that the estimated period matches your grating.")
                for b in self._result_buttons:
                    b.setEnabled(False)
                return

            self._result = result
            self._settings = settings
            for b in self._result_buttons:
                b.setEnabled(True)
            self._update_results_panel(result)
            self._show_angles()
            self._set_status(
                f"Analysis complete: {result['n_grooves']} measurements, "
                f"mean {result['mean_angle']:.2f}°")

        except Exception as exc:
            self._set_status(f"Analysis failed: {exc}")
            self.results_label.setText(f"Analysis failed:\n{exc}")

    def _update_results_panel(self, r):
        angles = np.asarray(r['all_angles'])
        r2 = np.array([q.get('blaze_r2', np.nan) for q in r['quality']])
        mode = f"row groups ×{r.get('n_groups', 1)}" if r.get('n_groups') else "averaged profile"
        self.results_label.setText(
            f"Mean blaze angle : {r['mean_angle']:.3f}°  ± {r.get('sem', 0):.3f}° (SEM)\n"
            f"Spread           : σ = {r['std_angle']:.3f}°   "
            f"range {r['min_angle']:.2f}–{r['max_angle']:.2f}°\n"
            f"Measurements     : N = {r['n_grooves']}   ({mode})\n"
            f"Period           : {r['period_nm']:.2f} ± {r.get('period_std', 0):.2f} nm   "
            f"worst fit R² = {np.nanmin(r2):.4f}")

    # ── Views ────────────────────────────────────────────────────────────────

    def _show_raw_profile(self):
        if self._data is not None:
            self.canvas.plot_raw(self._disp_um, self._profile_nm,
                                 self._filename, self._scan_size)

    def _show_2d(self):
        if self._data is not None:
            self.canvas.plot_2d(self._data, self._scan_size, self._filename)

    def _show_row_groups(self):
        if self._data is None:
            return
        n = self.n_groups_spin.value()
        with contextlib.redirect_stdout(io.StringIO()):
            disp_um, profiles, info = raw_data_multi_group(
                self._data, self._scan_size, n_groups=n)
        self._group_profiles, self._group_info = profiles, info
        self.canvas.plot_row_groups(disp_um, profiles, info,
                                    self._filename, self._scan_size)
        self._set_status(f"{info['n_groups']} groups × {info['rows_per_group']} rows")

    def _show_detection(self):
        if self._result is not None:
            self.canvas.plot_groove_detection(self._result, self._filename,
                                             self._settings)

    def _show_angles(self):
        if self._result is not None:
            self.canvas.plot_blaze_angles(self._result, self._filename)

    def _set_status(self, msg):
        self.status.showMessage(msg)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
