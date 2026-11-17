"""The coating stack as rows, drawn the way the light meets it.

One row per conformal layer, outermost at the top: a colour swatch, a
material (the library's names, or a typed model such as ``typical:Au``), a
thickness in nm, and buttons to move or remove it. "Add layer" puts a new one
on top, which is what depositing one does. Beneath the rows the caller places
the material the grating is cut into, so the whole column reads vacuum,
layers, grating.

The form still carries the stack as text, ``"Au 20, Cr 5"`` with the
innermost layer first (:class:`~gratinglab.gui.state.FormState.layers`), and
this widget speaks that through :meth:`text` and :meth:`setText`, the two
methods the geometry dock asks of every field. The order is reversed here and
nowhere else. Nothing is validated: a row's fields go into the text as typed,
and ``state.build`` reports what is wrong with them.

Built by :class:`~gratinglab.gui.qt.geometry_panel.GeometryPanel`.
"""

from __future__ import annotations

import zlib

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import (
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

__all__ = ["StackEditor", "material_color", "Swatch"]


def material_color(name: str) -> QColor:
    """One colour per material name, the same every run, so the same material
    is the same colour in every row; grey for none."""
    if not name:
        return QColor(150, 150, 150)
    # Twelve hues by three lightnesses, so two names rarely land close: a
    # hue taken straight from the hash put C and Au four degrees apart.
    code = zlib.crc32(name.encode("utf-8"))
    return QColor.fromHsl((code % 12) * 30, 120, (105, 140, 175)[(code // 12) % 3])


class Swatch(QFrame):
    """A small block of a material's colour, at the left of its row."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedWidth(10)
        self.setMinimumHeight(18)
        self.setAutoFillBackground(True)
        self.show_material("")

    def show_material(self, name: str) -> None:
        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, material_color(name.strip()))
        self.setPalette(palette)


class _Row(QWidget):
    """One layer: swatch, material, thickness, up, down, remove."""

    def __init__(self, editor: "StackEditor", material: str, thickness: str) -> None:
        super().__init__(editor)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)

        self.swatch = Swatch()
        self.material = QComboBox()
        self.material.setEditable(True)
        self.material.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.material.setMinimumContentsLength(6)
        self.material.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.material.addItems(editor.materials())
        self.material.setEditText(material)
        self.material.setToolTip(
            "A table from your library, or a typed model: typical:Au,\n"
            "factors:SiO2@2.2, or Au@17.5 for a table at another density."
        )
        self.thickness = QLineEdit(thickness)
        self.thickness.setFixedWidth(46)
        self.thickness.setToolTip("Thickness along the surface normal, nm.")
        unit = QLabel("nm")
        unit.setForegroundRole(QPalette.ColorRole.Mid)

        layout.addWidget(self.swatch)
        layout.addWidget(self.material, 1)
        layout.addWidget(self.thickness)
        layout.addWidget(unit)
        self.buttons = {}
        for key, arrow, tip in (
            ("up", Qt.ArrowType.UpArrow, "Move towards the surface"),
            ("down", Qt.ArrowType.DownArrow, "Move towards the grating"),
        ):
            button = QToolButton()
            button.setArrowType(arrow)
            button.setToolTip(tip)
            button.setAutoRaise(True)
            layout.addWidget(button)
            self.buttons[key] = button
        remove = QToolButton()
        remove.setText("×")
        remove.setToolTip("Remove this layer")
        remove.setAutoRaise(True)
        layout.addWidget(remove)
        self.buttons["remove"] = remove
        self.swatch.show_material(material)
        self.material.editTextChanged.connect(self.swatch.show_material)

    def entry(self) -> str:
        """``"Au 20"``, or just what was typed when a field is empty."""
        return " ".join(
            part for part in (self.material.currentText().strip(), self.thickness.text().strip()) if part
        )


class StackEditor(QWidget):
    """Rows of layers, outermost first, speaking the form's layers text."""

    #: A field changed; the dock redraws on it.
    changed = Signal()
    #: Enter in a thickness; the dock solves on it, as for any field.
    solve_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._materials: tuple[str, ...] = ()
        self._rows: list[_Row] = []

        column = QVBoxLayout(self)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(3)
        # The button is at the top because that is where a new layer goes.
        head = QHBoxLayout()
        self._top = QLabel("vacuum (light arrives)")
        self._top.setForegroundRole(QPalette.ColorRole.Mid)
        self._add = QPushButton("Add layer")
        head.addWidget(self._top, 1)
        head.addWidget(self._add)
        column.addLayout(head)
        self._rows_layout = QVBoxLayout()
        self._rows_layout.setSpacing(3)
        column.addLayout(self._rows_layout)
        self._add.setToolTip(
            "A conformal layer on top of the stack: it follows the surface\n"
            "beneath it as its true offset, so convex corners are rounded and\n"
            "concave ones stay sharp. Only the RCWA tab models layers on a\n"
            "grating; the Mirror tab reflects from the flat stack."
        )
        self._add.clicked.connect(lambda: self.add_layer())

    # -- the form's side -----------------------------------------------------

    def text(self) -> str:
        """The stack as ``FormState.layers`` carries it: innermost first."""
        return ", ".join(entry for row in reversed(self._rows) if (entry := row.entry()))

    def setText(self, text: str) -> None:  # noqa: N802 - the field protocol's name
        """Replace the rows from layers text, without announcing a change."""
        for row in self._rows:
            self._rows_layout.removeWidget(row)
            row.deleteLater()
        self._rows = []
        for entry in (part.strip() for part in text.split(",")):
            if not entry:
                continue
            material, _, thickness = entry.partition(" ")
            self._insert(0, material, thickness.strip())
        self._sync_buttons()

    # -- editing -------------------------------------------------------------

    def materials(self) -> tuple[str, ...]:
        return self._materials

    def set_materials(self, names) -> None:
        """The library's names, for every row's dropdown; typed text is kept."""
        self._materials = tuple(names)
        for row in self._rows:
            typed = row.material.currentText()
            row.material.blockSignals(True)
            row.material.clear()
            row.material.addItems(self._materials)
            row.material.setEditText(typed)
            row.material.blockSignals(False)

    def add_layer(self, material: str = "", thickness: str = "") -> None:
        """A new outermost layer, with the cursor in it."""
        row = self._insert(0, material, thickness)
        self._sync_buttons()
        row.material.setFocus()
        self.changed.emit()

    def layer_count(self) -> int:
        return len(self._rows)

    def _insert(self, index: int, material: str, thickness: str) -> _Row:
        row = _Row(self, material, thickness)
        row.material.editTextChanged.connect(lambda _t: self.changed.emit())
        row.thickness.textEdited.connect(lambda _t: self.changed.emit())
        row.thickness.returnPressed.connect(self.solve_requested)
        row.buttons["up"].clicked.connect(lambda: self._move(row, -1))
        row.buttons["down"].clicked.connect(lambda: self._move(row, +1))
        row.buttons["remove"].clicked.connect(lambda: self._remove(row))
        self._rows.insert(index, row)
        self._rows_layout.insertWidget(index, row)
        return row

    def _move(self, row: _Row, step: int) -> None:
        index = self._rows.index(row)
        target = index + step
        if not 0 <= target < len(self._rows):
            return
        self._rows.insert(target, self._rows.pop(index))
        self._rows_layout.removeWidget(row)
        self._rows_layout.insertWidget(target, row)
        self._sync_buttons()
        self.changed.emit()

    def _remove(self, row: _Row) -> None:
        self._rows.remove(row)
        self._rows_layout.removeWidget(row)
        row.deleteLater()
        self._sync_buttons()
        self.changed.emit()

    def _sync_buttons(self) -> None:
        last = len(self._rows) - 1
        for index, row in enumerate(self._rows):
            row.buttons["up"].setEnabled(index > 0)
            row.buttons["down"].setEnabled(index < last)
