"""
How much the averaged groove varies from swath to swath down the scan.

Opened from the Boundary tab. The tab averages every row into one groove; this
splits the rows into swaths, the same bands the blaze-angle analysis calls row
groups, builds the groove for each, and shows how far they disagree.

Computation lives in ``gratinglab.metrology.boundary.swath_statistics``.
"""
from __future__ import annotations

import os

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

import numpy as np
from PySide6.QtWidgets import (
    QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton,
    QSpinBox, QVBoxLayout, QWidget,
)
from matplotlib import colormaps
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from gratinglab.gui.qt import theming
from gratinglab.gui.theme import style_figure

from ...boundary import swath_statistics
from ...config import RESULTS_DIR

__all__ = ["SwathDialog"]


class SwathDialog(QDialog):
    """Per-swath groove shapes, their spread, and the scalars down the scan."""

    def __init__(self, data, scan_size, settings, filename=None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._data = data
        self._scan_size = scan_size
        self._settings = settings
        self._filename = filename
        self._stats = None

        self.setWindowTitle("Groove shape across swaths")
        self.resize(980, 640)
        self._build_ui()
        self.recompute()

    def _build_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)

        left = QWidget()
        left.setFixedWidth(300)
        column = QVBoxLayout(left)
        column.setContentsMargins(0, 0, 0, 0)

        row = QHBoxLayout()
        row.addWidget(QLabel("Swaths:"))
        self.n_spin = QSpinBox()
        self.n_spin.setRange(2, 100)
        self.n_spin.setValue(max(2, int(self._settings.n_row_groups)))
        self.n_spin.setToolTip(
            "Horizontal bands of rows, each averaged into its own groove. "
            "At least 3 rows each, so a short scan gets fewer.")
        self.n_spin.valueChanged.connect(self.recompute)
        row.addWidget(self.n_spin)
        column.addLayout(row)

        self.summary_label = QLabel("")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-family: monospace; font-size: 11px;")
        column.addWidget(self.summary_label)

        note = QLabel(
            "Between swaths: how much the groove changes along its length, down "
            "the image. Within one: how much grooves differ side by side. "
            "Swaths re-measure the same grooves, so these describe variation, "
            "not independent samples.")
        note.setWordWrap(True)
        note.setStyleSheet("font-size: 11px; color: gray;")
        column.addWidget(note)

        self.export_btn = QPushButton("Export per-swath CSV…")
        self.export_btn.clicked.connect(self.export)
        column.addWidget(self.export_btn)
        column.addStretch(1)
        root.addWidget(left)

        self.figure = Figure(figsize=(7, 6), layout="constrained")
        self.canvas = FigureCanvas(self.figure)
        root.addWidget(self.canvas, stretch=1)
        theming.subscribe(self._on_theme)

    @property
    def stats(self):
        """The current SwathStatistics, or None if the scan gave too few swaths."""
        return self._stats

    def recompute(self):
        try:
            self._stats = swath_statistics(
                self._data, self._scan_size, self._settings, n_swaths=self.n_spin.value())
        except ValueError as exc:
            self._stats = None
            self.export_btn.setEnabled(False)
            self.summary_label.setText(f"No swath statistics:\n{exc}")
            self.figure.clear()
            self.canvas.draw()
            return
        self.export_btn.setEnabled(True)
        text = self._stats.summary()
        if self._stats.skipped:
            text += "\n\nSkipped rows " + ", ".join(
                f"{a}-{b - 1}" for (a, b), _ in self._stats.skipped)
        self.summary_label.setText(text)
        self._draw()

    def _on_theme(self, theme):
        """Draw the figure again in the new theme.

        Each draw clears the figure and builds fresh axes, which take their
        colours from matplotlib's settings, so redrawing is the restyle. The
        figure predates the change and needs its own background set first.
        """
        style_figure(self.figure, theme)
        if self._stats is not None:
            self._draw()
        else:
            self.canvas.draw_idle()

    def _draw(self):
        s = self._stats
        self.figure.clear()
        grid = self.figure.add_gridspec(2, 3, height_ratios=[1.6, 1])
        shapes = self.figure.add_subplot(grid[0, :])
        colours = colormaps["viridis"](np.linspace(0.0, 0.9, s.n_swaths))
        for profile, colour, (first, end) in zip(s.profiles, colours, s.row_ranges):
            shapes.plot(s.t, profile, color=colour, lw=0.8, alpha=0.8)
        shapes.fill_between(s.t, s.mean - 2 * s.std, s.mean + 2 * s.std,
                            color="0.6", alpha=0.35, lw=0, label="±2σ between swaths")
        shapes.plot(s.t, s.whole.y_norm, color=theming.current().ink, lw=1.4, ls="--", label="whole scan")
        shapes.set_xlabel("position (fraction of period)")
        shapes.set_ylabel("height (fraction of period)")
        shapes.set_title(f"{s.n_swaths} swaths, coloured top to bottom of the scan")
        shapes.legend(loc="best", fontsize=8)

        centre = [(a + b - 1) / 2 for a, b in s.row_ranges]
        panels = (
            ("groove_depth", "depth (frac. of period)"),
            ("max_angle_deg", "max sidewall (deg)"),
        )
        for i, (key, label) in enumerate(panels):
            ax = self.figure.add_subplot(grid[1, i])
            ax.plot(centre, [p[key] for p in s.per_swath], "o-", ms=3)
            ax.set_xlabel("scan row")
            ax.set_ylabel(label)
        ax = self.figure.add_subplot(grid[1, 2])
        ax.plot(centre, s.shift_nm, "o-", ms=3)
        ax.set_xlabel("scan row")
        ax.set_ylabel("lateral shift (nm)")
        self.canvas.draw()

    def export(self):
        if self._stats is None:
            return
        stem = os.path.splitext(os.path.basename(self._filename or "scan"))[0]
        start = os.path.join(str(RESULTS_DIR), f"{stem}_swaths.csv")
        path, _ = QFileDialog.getSaveFileName(
            self, "Export per-swath statistics", start, "CSV table (*.csv)")
        if not path:
            return
        try:
            self._stats.write_csv(path)
        except OSError as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
