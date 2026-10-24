"""
The Import tab: load a scan, flatten it, see what that did.

This tab owns file loading and both flattening stages. The Analysis tab
consumes whatever this produces, so there is one place a scan enters the
application.

The tab exists mainly for its preview. Profile flattening moves the mean blaze
angle by about 0.49 degrees across its four methods, which is comparable to
the differences between samples, and until now it was a config-file setting
with no way to see what it was removing. Here the fitted background is drawn
over the profile it is about to be subtracted from.
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
from matplotlib.patches import Rectangle
from matplotlib.widgets import RectangleSelector

from ...config import PROJECT_ROOT
from ...core.image_flatten import (VALID_IMAGE_FLATTEN_METHODS, flatten_image,
                                   row_offset_spread)
from ...core.processing import flatten_profile, load_afm_data, raw_data
from ...core.region_mask import (apply_region_mask, bridge_nan, describe_mask,
                                 pixels_from_drag)
from ...core.scars import VALID_SCAR_REPAIRS, repair_scars_for
from ...settings import (AnalysisSettings, VALID_FLATTEN_FEATURES,
                         VALID_FLATTEN_METHODS, VALID_SPM_DIRECTIONS)

__all__ = ["ImportView"]


#: What a drag on the image does. A rectangle replaces any earlier one; line
#: and column ranges accumulate until the mask is cleared.
DRAG_MODES = ('rectangle', 'rows', 'columns')


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
        self._scars = None        # the ScarReport behind _flat

        # The region mask, in pixels of the image as loaded. Settings rather
        # than drawing state, so a batch run reproduces what was drawn here.
        self._mask_rect = self._defaults.mask_rect
        self._mask_rows = list(self._defaults.mask_rows)
        self._mask_cols = list(self._defaults.mask_cols)
        self._selector = None     # kept alive: matplotlib holds it weakly

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
        column.addWidget(self._build_cleanup_group())
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

    def _build_cleanup_group(self):
        group = QGroupBox("Scars and mask (2-D)")
        layout = QVBoxLayout(group)

        layout.addWidget(QLabel("Scar repair:"))
        self.scar_combo = QComboBox()
        self.scar_combo.addItems(list(VALID_SCAR_REPAIRS))
        self.scar_combo.setCurrentText(self._defaults.scar_repair)
        self.scar_combo.setToolTip(
            "Scan lines that stand off the lines either side of them: tip\n"
            "skips and feedback ringing, Gwyddion's 'scars'. Found after\n"
            "image flattening.\n\n"
            "interpolate: refill them down each column from the lines\n"
            "either side, which is where the grooves run.\n"
            "drop: leave those lines out of every row group instead.")
        self.scar_combo.currentTextChanged.connect(self._reflatten)
        layout.addWidget(self.scar_combo)

        row = QHBoxLayout()
        row.addWidget(QLabel("Threshold (sigma):"))
        self.scar_threshold_spin = QDoubleSpinBox()
        self.scar_threshold_spin.setRange(1.0, 50.0)
        self.scar_threshold_spin.setSingleStep(0.5)
        self.scar_threshold_spin.setDecimals(1)
        self.scar_threshold_spin.setValue(self._defaults.scar_threshold)
        self.scar_threshold_spin.setToolTip(
            "How far a scar must stand off both neighbouring lines, in units\n"
            "of the scan's own line-to-line noise. Lower finds more.")
        self.scar_threshold_spin.valueChanged.connect(self._reflatten)
        row.addWidget(self.scar_threshold_spin)
        layout.addLayout(row)

        self.scar_note = QLabel("Scar repair off")
        self.scar_note.setWordWrap(True)
        self.scar_note.setStyleSheet("font-size: 11px; color: gray;")
        layout.addWidget(self.scar_note)

        layout.addWidget(QLabel("Drag on the image to exclude:"))
        self.drag_mode_combo = QComboBox()
        self.drag_mode_combo.addItems(list(DRAG_MODES))
        self.drag_mode_combo.setToolTip(
            "rectangle: a patch such as a particle (one at a time).\n"
            "rows: whole scan lines; they join no row group.\n"
            "columns: whole columns; no groove is measured across them.")
        layout.addWidget(self.drag_mode_combo)

        self.mask_label = QLabel("No mask")
        self.mask_label.setWordWrap(True)
        self.mask_label.setStyleSheet("font-size: 11px;")
        layout.addWidget(self.mask_label)

        clear = QPushButton("Clear mask")
        clear.clicked.connect(self.clear_mask)
        layout.addWidget(clear)
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
            scar_repair=self.scar_combo.currentText(),
            scar_threshold=self.scar_threshold_spin.value(),
            mask_rect=self._mask_rect,
            mask_rows=tuple(self._mask_rows),
            mask_cols=tuple(self._mask_cols),
        )

    @property
    def scar_report(self):
        """What scar repair found on the current scan, or None."""
        return self._scars

    # ── Region mask ──────────────────────────────────────────────────────────

    def set_mask(self, rect=None, rows=(), cols=()):
        """Replace the whole mask. Pixel indices, half-open, as in settings."""
        self._mask_rect = tuple(int(v) for v in rect) if rect else None
        self._mask_rows = [(int(a), int(b)) for a, b in rows]
        self._mask_cols = [(int(a), int(b)) for a, b in cols]
        self._mask_changed()

    def clear_mask(self):
        self.set_mask()

    def add_drag(self, x0, y0, x1, y1):
        """
        A drag from (x0, y0) to (x1, y1) in image coordinates (micrometres,
        scan lines), added to the mask according to the drag mode.
        """
        if self._flat is None:
            return
        r0, r1, c0, c1 = pixels_from_drag(x0, y0, x1, y1, self._flat.shape,
                                          self._scan_size)
        mode = self.drag_mode_combo.currentText()
        if mode == 'rows':
            self._mask_rows.append((r0, r1))
        elif mode == 'columns':
            self._mask_cols.append((c0, c1))
        else:
            self._mask_rect = (r0, r1, c0, c1)
        self._mask_changed()

    def _on_select(self, press, release):
        if None in (press.xdata, press.ydata, release.xdata, release.ydata):
            return
        self.add_drag(press.xdata, press.ydata, release.xdata, release.ydata)

    def _mask_changed(self):
        text = describe_mask(self._mask_settings())
        self.mask_label.setText(
            f"Excluded: {text}\n(pixels, half-open; MASK_* in config.py)"
            if text else "No mask")
        self._reflatten()

    def _mask_settings(self):
        return self._defaults.with_(mask_rect=self._mask_rect,
                                    mask_rows=tuple(self._mask_rows),
                                    mask_cols=tuple(self._mask_cols))

    def set_period_est(self, period_nm: float) -> None:
        """
        The period lives with the analysis controls, but `level_grooves` needs it
        to find the features it levels on. The window pushes it here so the
        preview matches what the analysis will do.
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
            "AFM data (*.txt *.dat *.asc *.spm *.ibw);;"
            "Nanoscope (*.spm);;Asylum (*.ibw);;Text export (*.txt);;"
            "All files (*)")
        if path:
            self.load(path)

    def _reload(self):
        if self._filename:
            self.load(self._filename)

    def load(self, path):
        """Read a file, flatten it, and tell the window."""
        try:
            from ...io.ibw import is_igor_wave_file
            from ...io.spm import is_nanoscope_file
            nanoscope = is_nanoscope_file(path)
            asylum = not nanoscope and is_igor_wave_file(path)

            self.direction_combo.blockSignals(True)
            self.direction_combo.setEnabled(nanoscope or asylum)
            self.direction_combo.blockSignals(False)

            if path != self._filename:
                # A mask belongs to the scan it was drawn on.
                self._mask_rect, self._mask_rows, self._mask_cols = None, [], []
                self.mask_label.setText("No mask")

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
            if nanoscope:
                source = (f"Nanoscope: {settings.spm_channel} / "
                          f"{settings.spm_direction}")
            elif asylum:
                from ...io.ibw import _choose, list_channels
                chosen = _choose(list(list_channels(path)), None,
                                 settings.spm_direction)
                source = f"Asylum: {chosen.channel} / {chosen.direction}"
            else:
                source = "Text export"
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
        # The same order as the analysis: mask, level, repair scars. The mask
        # comes first because levelling reads every pixel of a line, so what
        # this tab hands on, NaN where excluded, is what the analysis sees.
        masked, _ = apply_region_mask(self._raw, self._mask_settings())
        before = row_offset_spread(masked)
        flat = flatten_image(masked, method)
        after = row_offset_spread(flat)
        self.image_note.setText(
            f"Row-offset spread\n  {before*1e9:.2f} → {after*1e9:.2f} nm")

        settings = self._defaults.with_(
            scar_repair=self.scar_combo.currentText(),
            scar_threshold=self.scar_threshold_spin.value())
        self._flat, self._scars = repair_scars_for(flat, settings)
        self.scar_note.setText(self._scars.summary.capitalize())
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

        # Line 0 at the top and the y axis counting lines, so what the axis
        # reads is what the mask stores (core.region_mask.pixels_from_drag).
        n_rows, n_cols = self._flat.shape
        image = image_ax.imshow(
            self._flat * 1e9, aspect='auto', cmap='viridis',
            extent=[0, self._scan_size, n_rows, 0], origin='upper')
        self.figure.colorbar(image, ax=image_ax, label='Height (nm)')
        image_ax.set_xlabel("X (µm)", fontsize=10)
        image_ax.set_ylabel("Scan line", fontsize=10)
        self._draw_mask(image_ax, n_rows, n_cols)
        self._selector = RectangleSelector(
            image_ax, self._on_select, useblit=False, button=[1],
            minspanx=0, minspany=0, spancoords='data', interactive=False,
            props=dict(facecolor='crimson', edgecolor='crimson', alpha=0.25,
                       fill=True))
        image_ax.set_title(
            f"{os.path.basename(self._filename)}: image flattening: "
            f"{self.image_method_combo.currentText()}", fontsize=10)

        try:
            x_um, profile_nm = raw_data(self._flat, self._scan_size)
            with contextlib.redirect_stdout(io.StringIO()):
                flat, background = flatten_profile(
                    x_um, bridge_nan(profile_nm),
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

        # Masked columns are gaps, not the straight bridge the fit used.
        gaps = ~np.isfinite(profile_nm)
        profile_ax.plot(x_um, profile_nm, color='0.6', linewidth=0.9,
                        label='Averaged profile')
        profile_ax.plot(x_um, background, color='crimson', linewidth=1.4,
                        linestyle='--', label='Background to remove')
        profile_ax.plot(x_um, np.where(gaps, np.nan, flat), color='steelblue',
                        linewidth=1.0, label='After flattening')
        profile_ax.set_xlabel("Position (µm)", fontsize=10)
        profile_ax.set_ylabel("Height (nm)", fontsize=10)
        profile_ax.set_title(
            f"Profile flattening: {self.profile_method_combo.currentText()}"
            f"  (removes {np.ptp(background):.1f} nm of background)", fontsize=10)
        profile_ax.legend(fontsize=8)
        profile_ax.grid(True, alpha=0.3)

        self.canvas.draw()

    def _draw_mask(self, axes, n_rows, n_cols):
        """Hatch the excluded pixels and tick the scar lines."""
        width_um = self._scan_size
        px = width_um / n_cols

        def patch(r0, r1, c0, c1):
            axes.add_patch(Rectangle(
                (c0 * px, r0), (c1 - c0) * px, r1 - r0, fill=False,
                hatch='///', edgecolor='crimson', linewidth=1.0))

        if self._mask_rect:
            patch(*self._mask_rect)
        for r0, r1 in self._mask_rows:
            patch(r0, r1, 0, n_cols)
        for c0, c1 in self._mask_cols:
            patch(0, n_rows, c0, c1)

        if self._scars is not None and self._scars.n_lines:
            for line in self._scars.lines:
                axes.plot([width_um * 0.985, width_um], [line + 0.5] * 2,
                          color='orange', linewidth=2.0)
