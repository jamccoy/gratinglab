"""
AFM Blaze Angle Analysis - GUI
v0.1: File loader + raw data plot

Requirements:
    pip install PyQt5 matplotlib numpy scipy
"""

import sys
import os
import numpy as np

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QDoubleSpinBox, QSpinBox, QGroupBox,
    QSplitter, QStatusBar, QFrame
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


# ── Attempt to import your project's data loader ──────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_DIR)

try:
    from afm_analysis.core.processing import load_afm_data, raw_data
    LOADER_AVAILABLE = True
except ImportError:
    LOADER_AVAILABLE = False


def _fallback_load(filename, scan_x_size=2.0):
    """
    Minimal fallback loader — mirrors what load_afm_data + raw_data do,
    used when the project package isn't on the path.
    """
    import re
    detected_size = None

    try:
        with open(filename, 'r') as f:
            header_lines = [f.readline() for _ in range(10)]
        for line in header_lines:
            if 'width' in line.lower() or ('scan' in line.lower() and 'size' in line.lower()):
                m = re.search(r'(\d+\.?\d*)\s*(um|µm|micron)', line, re.IGNORECASE)
                if m:
                    detected_size = float(m.group(1))
                    break
                m = re.search(r'(\d+\.?\d*)\s*nm', line, re.IGNORECASE)
                if m:
                    detected_size = float(m.group(1)) / 1000
                    break
    except Exception:
        pass

    actual_size = detected_size if detected_size else scan_x_size
    data = np.genfromtxt(filename, skip_header=4)

    disp_um = actual_size * np.arange(data.shape[1]) / (data.shape[1] - 1)
    profile_nm = 1e9 * (np.mean(data, axis=0) - np.min(np.mean(data, axis=0)))

    return data, actual_size, disp_um, profile_nm


def load_file(filename, scan_x_size=2.0):
    """Unified loader — uses project code if available, otherwise fallback."""
    if LOADER_AVAILABLE:
        data, actual_size = load_afm_data(filename, default_scan_size=scan_x_size)
        disp_um, profile_nm = raw_data(data, actual_size)
        return data, actual_size, disp_um, profile_nm
    else:
        return _fallback_load(filename, scan_x_size)


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
        ax.set_title(
            f"Raw AFM Profile — {os.path.basename(filename)}\n"
            f"Scan width: {scan_size_um:.3f} µm",
            fontsize=11
        )
        ax.grid(True, alpha=0.3)
        self.draw()

    def plot_row_groups(self, disp_um, profiles_nm, group_info, filename, scan_size_um):
        self.fig.clear()
        ax = self.fig.add_subplot(111)

        n_groups = group_info['n_groups']
        cmap = matplotlib.colormaps['coolwarm']

        # Plot each group profile with a color gradient (top→bottom of image)
        for i, profile in enumerate(profiles_nm):
            color = cmap(i / max(n_groups - 1, 1))
            ax.plot(disp_um, profile, color=color, linewidth=0.8, alpha=0.6)

        # Mean profile on top
        mean_profile = np.mean(np.array(profiles_nm), axis=0)
        ax.plot(disp_um, mean_profile, color='black', linewidth=1.8,
                label=f'Mean (n={n_groups})', zorder=5)

        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        ax.set_title(
            f"Row-Group Profiles — {os.path.basename(filename)}\n"
            f"{n_groups} groups × {group_info['rows_per_group']} rows each  |  "
            f"Scan width: {scan_size_um:.3f} µm",
            fontsize=10
        )
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Colorbar to show spatial position in image
        sm = matplotlib.cm.ScalarMappable(
            cmap=cmap,
            norm=matplotlib.colors.Normalize(vmin=0, vmax=group_info['n_rows'])
        )
        sm.set_array([])
        cb = self.fig.colorbar(sm, ax=ax, pad=0.02)
        cb.set_label("Image row (top → bottom)", fontsize=9)

        self.draw()

    def plot_2d(self, data, scan_size_um, filename):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        im = ax.imshow(
            data * 1e9,
            aspect='auto',
            cmap='viridis',
            extent=[0, scan_size_um, 0, data.shape[0]],
            origin='upper'
        )
        self.fig.colorbar(im, ax=ax, label='Height (nm)')
        ax.set_xlabel("X (µm)", fontsize=11)
        ax.set_ylabel("Scan line", fontsize=11)
        ax.set_title(
            f"2D AFM Topography — {os.path.basename(filename)}",
            fontsize=11
        )
        self.draw()


# ── Main window ───────────────────────────────────────────────────────────────

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AFM Blaze Angle Analysis")
        self.resize(1100, 680)

        self._data = None
        self._disp_um = None
        self._profile_nm = None
        self._scan_size = None
        self._filename = None
        self._group_profiles = None   # list of per-group 1D arrays
        self._group_info = None       # metadata dict from raw_data_multi_group

        self._build_ui()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # ── Left panel (controls) ──────────────────────────────────────────
        left = QWidget()
        left.setFixedWidth(240)
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(10)

        # File group
        file_group = QGroupBox("File")
        file_layout = QVBoxLayout(file_group)
        self.file_label = QLabel("No file loaded")
        self.file_label.setWordWrap(True)
        self.file_label.setStyleSheet("color: gray; font-size: 11px;")
        file_layout.addWidget(self.file_label)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._browse_file)
        file_layout.addWidget(browse_btn)
        left_layout.addWidget(file_group)

        # Scan parameters group
        params_group = QGroupBox("Scan Parameters")
        params_layout = QVBoxLayout(params_group)
        params_layout.addWidget(QLabel("Default scan width (µm):"))
        self.scan_size_spin = QDoubleSpinBox()
        self.scan_size_spin.setRange(0.1, 1000.0)
        self.scan_size_spin.setValue(2.0)
        self.scan_size_spin.setSingleStep(0.1)
        self.scan_size_spin.setDecimals(3)
        self.scan_size_spin.setToolTip(
            "Used only when scan width cannot be detected from the file header."
        )
        params_layout.addWidget(self.scan_size_spin)
        reload_btn = QPushButton("Reload with these settings")
        reload_btn.clicked.connect(self._reload_file)
        params_layout.addWidget(reload_btn)
        left_layout.addWidget(params_group)

        # View options group
        view_group = QGroupBox("View")
        view_layout = QVBoxLayout(view_group)
        self.raw_profile_btn = QPushButton("Raw Profile (1D)")
        self.raw_profile_btn.setEnabled(False)
        self.raw_profile_btn.clicked.connect(self._show_raw_profile)
        view_layout.addWidget(self.raw_profile_btn)
        self.topography_btn = QPushButton("2D Topography")
        self.topography_btn.setEnabled(False)
        self.topography_btn.clicked.connect(self._show_2d)
        view_layout.addWidget(self.topography_btn)

        self.row_groups_btn = QPushButton("Row Groups")
        self.row_groups_btn.setEnabled(False)
        self.row_groups_btn.clicked.connect(self._show_row_groups)
        view_layout.addWidget(self.row_groups_btn)

        # N groups spinner
        n_groups_layout = QHBoxLayout()
        n_groups_layout.addWidget(QLabel("N groups:"))
        self.n_groups_spin = QSpinBox()
        self.n_groups_spin.setRange(2, 200)
        self.n_groups_spin.setValue(20)
        self.n_groups_spin.setEnabled(False)
        self.n_groups_spin.setToolTip(
            "Number of horizontal bands to divide the image into.\n"
            "Each band is averaged separately to give one profile."
        )
        self.n_groups_spin.valueChanged.connect(self._on_n_groups_changed)
        n_groups_layout.addWidget(self.n_groups_spin)
        view_layout.addLayout(n_groups_layout)

        left_layout.addWidget(view_group)

        # File info group
        self.info_group = QGroupBox("File Info")
        info_layout = QVBoxLayout(self.info_group)
        self.info_label = QLabel("—")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("font-size: 11px;")
        info_layout.addWidget(self.info_label)
        left_layout.addWidget(self.info_group)

        left_layout.addStretch()

        # ── Right panel (plot) ─────────────────────────────────────────────
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = PlotCanvas(self)
        toolbar = NavigationToolbar(self.canvas, self)
        right_layout.addWidget(toolbar)
        right_layout.addWidget(self.canvas)

        # ── Assemble ───────────────────────────────────────────────────────
        root.addWidget(left)
        divider = QFrame()
        divider.setFrameShape(QFrame.VLine)    # PyQt5: no Shape. prefix
        divider.setFrameShadow(QFrame.Sunken)  # PyQt5: no Shadow. prefix
        root.addWidget(divider)
        root.addWidget(right, stretch=1)

        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self._set_status("Ready. Open an AFM file to begin.")

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _browse_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open AFM File", "",
            "AFM Data Files (*.txt *.dat *.asc);;All Files (*)"
        )
        if path:
            self._load(path)

    def _reload_file(self):
        if self._filename:
            self._load(self._filename)

    def _load(self, path):
        try:
            self._set_status(f"Loading {os.path.basename(path)}…")
            QApplication.processEvents()

            data, scan_size, disp_um, profile_nm = load_file(
                path, self.scan_size_spin.value()
            )

            self._data = data
            self._scan_size = scan_size
            self._disp_um = disp_um
            self._profile_nm = profile_nm
            self._filename = path

            self.file_label.setText(os.path.basename(path))
            self.file_label.setStyleSheet("font-size: 11px;")
            self.scan_size_spin.setValue(scan_size)

            rows, cols = data.shape
            self.info_label.setText(
                f"Shape: {rows} × {cols} px\n"
                f"Scan width: {scan_size:.3f} µm\n"
                f"Height range: {profile_nm.min():.1f} – {profile_nm.max():.1f} nm\n"
                f"Profile points: {len(profile_nm)}"
            )

            self.raw_profile_btn.setEnabled(True)
            self.topography_btn.setEnabled(True)
            self.row_groups_btn.setEnabled(True)
            self.n_groups_spin.setEnabled(True)

            # Cap N groups to max sensible value for this image
            max_groups = data.shape[0] // 3
            self.n_groups_spin.setMaximum(max_groups)
            if self.n_groups_spin.value() > max_groups:
                self.n_groups_spin.setValue(max_groups)

            # Clear cached group data from any previous file
            self._group_profiles = None
            self._group_info = None

            self._show_raw_profile()
            self._set_status(f"Loaded: {os.path.basename(path)}")

        except Exception as e:
            self._set_status(f"Error loading file: {e}")
            self.info_label.setText(f"Error:\n{e}")

    def _show_raw_profile(self):
        if self._data is not None:
            self.canvas.plot_raw(
                self._disp_um, self._profile_nm,
                self._filename, self._scan_size
            )

    def _show_2d(self):
        if self._data is not None:
            self.canvas.plot_2d(self._data, self._scan_size, self._filename)

    def _show_row_groups(self):
        if self._data is None:
            return
        self._compute_row_groups()
        if self._group_profiles is not None:
            self.canvas.plot_row_groups(
                self._disp_um, self._group_profiles,
                self._group_info, self._filename, self._scan_size
            )

    def _on_n_groups_changed(self):
        """Re-compute and re-plot when N groups spinner changes, if that view is active."""
        self._group_profiles = None
        self._group_info = None
        # Only re-plot if row groups view is currently shown
        # (checking canvas title is fragile; just always re-plot if data loaded)
        if self._data is not None and self.row_groups_btn.isEnabled():
            self._show_row_groups()

    def _compute_row_groups(self):
        """Run raw_data_multi_group (or fallback) and cache the result."""
        if self._group_profiles is not None:
            return  # already computed for current settings

        n = self.n_groups_spin.value()

        try:
            from afm_analysis.core.processing import raw_data_multi_group
            disp_um, profiles_nm, group_info = raw_data_multi_group(
                self._data, self._scan_size, n_groups=n
            )
        except ImportError:
            # Inline fallback mirroring raw_data_multi_group logic
            disp_um, profiles_nm, group_info = self._fallback_multi_group(n)

        self._group_profiles = profiles_nm
        self._group_info = group_info

        # Update info label with group stats
        rows, cols = self._data.shape
        gi = group_info
        self.info_label.setText(
            f"Shape: {rows} × {cols} px\n"
            f"Scan width: {self._scan_size:.3f} µm\n"
            f"Groups: {gi['n_groups']}  ×  {gi['rows_per_group']} rows\n"
            f"Height range: {self._profile_nm.min():.1f} – {self._profile_nm.max():.1f} nm"
        )
        self._set_status(
            f"Row groups computed: {gi['n_groups']} groups × {gi['rows_per_group']} rows"
        )

    def _fallback_multi_group(self, n_groups):
        """Inline fallback when project package is not importable."""
        data = self._data
        scan_x_size = self._scan_size
        n_rows, n_cols = data.shape
        disp_um = scan_x_size * np.arange(n_cols) / (n_cols - 1)

        rows_per_group = n_rows // n_groups
        if rows_per_group < 3:
            n_groups = max(1, n_rows // 3)
            rows_per_group = n_rows // n_groups

        profiles_nm, group_ranges = [], []
        for i in range(n_groups):
            start = i * rows_per_group
            end = n_rows if i == n_groups - 1 else (i + 1) * rows_per_group
            avg = np.mean(data[start:end, :], axis=0)
            profile = 1e9 * (avg - np.min(avg))
            profiles_nm.append(profile)
            group_ranges.append((start, end))

        group_info = {
            'n_groups': n_groups,
            'n_rows': n_rows,
            'rows_per_group': rows_per_group,
            'group_ranges': group_ranges
        }
        return disp_um, profiles_nm, group_info

    def _set_status(self, msg):
        self.status.showMessage(msg)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())   # PyQt5 uses exec_() not exec()


if __name__ == "__main__":
    main()
