"""The Mirror tab: the coating stack as a flat mirror, before any grooves.

Reads the stack from the geometry dock through ``read_form`` (Coating is the
substrate, then Layers and Roughness) and its own fields for the scan, and
shows the reflectivity, how deep the field reaches into each material, and
whether the outermost layer is thick enough to treat as bulk. Everything it
computes and says is in :mod:`gratinglab.gui.mirror_view`; this is widgets.

Built by :class:`~gratinglab.gui.qt.main_window.MainWindow`, beside the
Grating Geometry tab. It solves in the GUI thread, since a scan takes about a
second; it is not a solver and has no worker.
"""

from __future__ import annotations

from typing import Callable

from endstation.qt.canvas import restyle_toolbar
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import provenance
from ..mirror_view import MirrorResult, MirrorSettings, compute, draw, lines
from ..state import FormErrors, FormState
from ..theme import provenance_colors, style_figure
from . import theming
from .sizing import fit_width_to_contents

__all__ = ["MirrorTab"]


class MirrorTab(QWidget):
    """Reflectivity and penetration of the dock's coating stack."""

    def __init__(self, read_form: Callable[[], FormState], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._read_form = read_form
        self._result: MirrorResult | None = None
        self._lines: tuple = ()
        self._fields: dict[str, QLineEdit] = {}
        self._build()

    # -- construction ----------------------------------------------------

    def _build(self) -> None:
        controls = QWidget()
        column = QVBoxLayout(controls)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self._build_scan_group())
        self._calculate = QPushButton("Calculate")
        self._calculate.setToolTip("Reflect from the dock's Surface stack as a flat mirror.")
        self._calculate.clicked.connect(self.calculate)
        column.addWidget(self._calculate)
        column.addStretch(1)

        scroller = QScrollArea()
        scroller.setWidget(controls)
        scroller.setWidgetResizable(True)
        fit_width_to_contents(scroller)

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self._build_plot())
        right.addWidget(self._build_panel())
        right.setStretchFactor(0, 1)
        right.setSizes([560, 160])

        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(scroller)
        body.addWidget(right)
        body.setStretchFactor(1, 1)
        body.setSizes([260, 700])

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(body)

    def _build_scan_group(self) -> QGroupBox:
        defaults = MirrorSettings()
        group = QGroupBox("Scan")
        form = QFormLayout(group)

        self._scan = QComboBox()
        self._scan.addItem("Graze angle", "angle")
        self._scan.addItem("Wavelength (dock's scan)", "wavelength")
        self._scan.currentIndexChanged.connect(self._on_scan_change)
        form.addRow("Scan", self._scan)

        def entry(label: str, key: str, tip: str) -> None:
            field = QLineEdit(getattr(defaults, key))
            field.setToolTip(tip)
            field.returnPressed.connect(self.calculate)
            self._fields[key] = field
            form.addRow(label, field)

        entry("From (deg)", "graze_start", "Graze angle from the surface; 0 is grazing.")
        entry("To (deg)", "graze_stop", "Graze angle from the surface.")
        entry("Points", "graze_count", "Angles in the scan, at most 1000.")
        entry("At wavelength (nm)", "wavelength", "The one wavelength an angle scan uses.")
        entry("At graze (deg)", "graze", "The one graze angle a wavelength scan uses.")

        self._polarization = QComboBox()
        for label, value in (("Unpolarized", "unpolarized"), ("s", "s"), ("p", "p")):
            self._polarization.addItem(label, value)
        self._polarization.setToolTip("s and p are local to the mirror's plane of incidence.")
        form.addRow("Polarization", self._polarization)

        entry("Bulk within", "tolerance",
              "The outermost layer reads as bulk once the layers beneath it\n"
              "change the reflectivity by less than this, relative.")
        self._log = QCheckBox("Log reflectivity")
        self._log.setChecked(True)
        self._log.toggled.connect(self._redraw)
        form.addRow(self._log)
        self._on_scan_change()
        return group

    def _build_plot(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 5), layout="constrained")
        grid = self._figure.add_gridspec(2, 1, height_ratios=(1, 1))
        self._r_axes = self._figure.add_subplot(grid[0])
        self._depth_axes = self._figure.add_subplot(grid[1], sharex=self._r_axes)
        self._canvas = FigureCanvasQTAgg(self._figure)
        style_figure(self._figure, theming.current())
        theming.subscribe(self._on_theme)

        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        self._toolbar = NavigationToolbar2QT(self._canvas, panel)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas)
        return panel

    def _build_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 4, 8, 8)
        heading = QLabel("Mirror")
        bold = heading.font()
        bold.setBold(True)
        heading.setFont(bold)
        self._panel = QTextEdit()
        self._panel.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._panel)
        self._paint((provenance.Line(
            "The dock's Surface stack, taken flat: its Grating material is the "
            "substrate. Press Calculate to see how far the light goes in, and "
            "whether the top layer can be treated as bulk.\n", "dim",
        ),))
        return panel

    # -- what the window and tests ask ---------------------------------------

    def settings(self) -> MirrorSettings:
        values = {key: field.text() for key, field in self._fields.items()}
        return MirrorSettings(scan=self._scan.currentData(), polarization=self._polarization.currentData(), **values)

    def calculate(self) -> None:
        try:
            self._result = compute(self._read_form(), self.settings())
        except FormErrors as exc:
            self._result = None
            self._paint(provenance.error_lines(exc.errors))
            return
        except (ValueError, KeyError) as exc:
            # A material that does not cover the scan, or one not installed.
            self._result = None
            self._paint((provenance.Line(f"{exc}\n", "bad"),))
            return
        self._redraw()

    # -- internals -------------------------------------------------------

    def _on_scan_change(self, *_args) -> None:
        angle = self._scan.currentData() == "angle"
        for key in ("graze_start", "graze_stop", "graze_count", "wavelength"):
            self._fields[key].setEnabled(angle)
        self._fields["graze"].setEnabled(not angle)

    def _paint(self, lines_) -> None:
        self._lines = tuple(lines_)
        self._panel.setHtml(provenance.to_html(self._lines, provenance_colors(theming.current())))

    def _on_theme(self, theme) -> None:
        restyle_toolbar(self._toolbar, theme)
        style_figure(self._figure, theme)
        self._redraw()
        if self._lines:
            self._paint(self._lines)

    def _redraw(self, *_args) -> None:
        if self._result is not None:
            draw(self._r_axes, self._depth_axes, self._result, theming.current(), self._log.isChecked())
            self._paint(lines(self._result))
        self._canvas.draw_idle()
