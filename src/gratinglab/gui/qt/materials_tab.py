"""The Materials tab: name a material, see its index, send it to the dock.

Three fields build the name (the kind of material, its table name or formula,
and an optional density), so ``factors:SiO2@2.2`` never has to be typed from
memory. The plot and the panel follow as the fields change: the index against
wavelength with the dock's scan shaded, where the numbers came from, and
whether they cover that scan. Two buttons put the material into the dock's
Surface stack, as the grating material or as a new top layer.

The library is here too: a list of its tables, where a click previews one,
above the three actions that change it (:mod:`.library_bar`).

What is resolved, drawn and said is in :mod:`gratinglab.gui.materials_view`;
this is widgets. Built by
:class:`~gratinglab.gui.qt.coatings_window.CoatingsWindow`, which hands it the
dock.
"""

from __future__ import annotations

from endstation.qt.canvas import restyle_toolbar
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import provenance
from ..materials_view import KINDS, MaterialPreview, draw, lines, material_name, preview
from ..theme import provenance_colors, style_figure
from . import theming
from .library_bar import LibraryBar
from .sizing import fit_width_to_contents

__all__ = ["MaterialsTab"]


class MaterialsTab(QWidget):
    """Preview one material and hand it to the geometry dock."""

    def __init__(self, dock, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        #: The window's dock: read for its scan, written by the two buttons.
        self._dock = dock
        self._preview: MaterialPreview | None = None
        self._lines: tuple = ()
        self._build()
        self.refresh_library()
        dock.changed.connect(self.update_preview)  # the scan may have moved

    # -- construction ----------------------------------------------------

    def _build(self) -> None:
        controls = QWidget()
        column = QVBoxLayout(controls)
        column.setContentsMargins(0, 0, 0, 0)

        group = QGroupBox("Material")
        form = QFormLayout(group)
        self._kind = QComboBox()
        for key, label, _field in KINDS:
            self._kind.addItem(label, key)
        self._kind.currentIndexChanged.connect(self._on_kind_change)
        form.addRow("Kind", self._kind)

        self._name_label = QLabel("Name")
        self._name = QComboBox()
        self._name.setEditable(True)
        self._name.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._name.editTextChanged.connect(self.update_preview)
        form.addRow(self._name_label, self._name)

        self._density = QLineEdit()
        self._density.setToolTip(
            "g/cm^3. Optional: a library table is rescaled to it (CXRO exports\n"
            "only), and a formula without one uses its default where there is one."
        )
        self._density.textEdited.connect(self.update_preview)
        form.addRow("Density (g/cm³)", self._density)

        self._result = QLineEdit()
        self._result.setReadOnly(True)
        self._result.setToolTip("The name to use as a coating, here or in Python.")
        form.addRow("Coating name", self._result)
        column.addWidget(group)

        self._use = QPushButton("Use as substrate")
        self._use.setToolTip("Set the dock's Substrate row to this material.")
        self._use.clicked.connect(self._use_as_grating)
        self._layer = QPushButton("Add as top layer")
        self._layer.setToolTip("Add a layer of this material on top of the dock's stack.")
        self._layer.clicked.connect(self._add_as_layer)
        column.addWidget(self._use)
        column.addWidget(self._layer)

        shelf = QGroupBox("Your library")
        shelf_column = QVBoxLayout(shelf)
        self._list = QListWidget()
        self._list.setToolTip("The tables in your library. Click one to preview it.")
        self._list.itemClicked.connect(lambda item: self.set_fields("table", item.text()))
        shelf_column.addWidget(self._list, 1)
        self._elements = QLabel("")
        self._elements.setWordWrap(True)
        self._elements.setToolTip("The elements a scattering-factors formula can use here.")
        shelf_column.addWidget(self._elements)
        #: The library's actions; the window listens to its ``library_changed``.
        self.library = LibraryBar()
        shelf_column.addWidget(self.library)
        column.addWidget(shelf, 1)

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
        body.setSizes([280, 700])

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(body)

    def _build_plot(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 5), layout="constrained")
        self._axes = self._figure.add_subplot()
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
        heading = QLabel("Material")
        bold = heading.font()
        bold.setBold(True)
        heading.setFont(bold)
        self._panel = QTextEdit()
        self._panel.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._panel)
        return panel

    # -- what the window and tests ask ---------------------------------------

    def refresh_library(self) -> None:
        """Re-read the library's names into the Name dropdown."""
        from ...materials import available

        try:
            self._library = available()
        except ValueError:  # two tables for one name; Setup reports it
            self._library = ()
        self._list.clear()
        self._list.addItems(self._library)
        from ...materials import factor_elements

        elements = factor_elements()
        self._elements.setText(
            "Scattering factors: " + (", ".join(elements) if elements else "none installed")
        )
        self._on_kind_change()

    def coating_name(self) -> str:
        return material_name(self._kind.currentData(), self._name.currentText(), self._density.text())

    def set_fields(self, kind: str, name: str, density: str = "") -> None:
        self._kind.setCurrentIndex(self._kind.findData(kind))
        self._name.setEditText(name)
        self._density.setText(density)
        self.update_preview()

    def update_preview(self, *_args) -> None:
        name = self.coating_name()
        self._result.setText(name)
        self._preview = None
        if not name:
            self._paint((provenance.Line(
                "Choose a table from your library, or enter a formula, to see "
                "its index and whether it covers the dock's wavelengths.\n", "dim"),))
        else:
            form = self._dock.read_form()
            try:
                self._preview = preview(name, form.wavelength_start, form.wavelength_stop)
            except (ValueError, KeyError, OSError, ImportError) as refused:
                self._paint((provenance.Line(f"{refused}\n", "bad"),))
        for button in (self._use, self._layer):
            button.setEnabled(self._preview is not None)
        self._redraw()

    # -- internals -------------------------------------------------------

    def _on_kind_change(self, *_args) -> None:
        from ...materials.typical import DENSITIES

        kind = self._kind.currentData()
        self._name_label.setText(next(field for key, _l, field in KINDS if key == kind))
        typed = self._name.currentText()
        self._name.blockSignals(True)
        self._name.clear()
        # A table's name is one of the library's; a formula is free, with
        # the ones that have a default density offered as a start.
        self._name.addItems(self._library if kind == "table" else tuple(DENSITIES))
        self._name.setEditText(typed)
        self._name.blockSignals(False)
        self._density.setPlaceholderText("as tabulated" if kind == "table" else "default")
        self.update_preview()

    def _use_as_grating(self) -> None:
        self._dock.set_grating_material(self.coating_name())

    def _add_as_layer(self) -> None:
        self._dock.add_top_layer(self.coating_name())

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
        if self._preview is not None:
            draw(self._axes, self._preview, theming.current())
            self._paint(lines(self._preview))
        else:
            self._axes.clear()
            style_figure(self._figure, theming.current())
        self._canvas.draw_idle()
