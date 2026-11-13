"""The provenance strip under each solver tab's plot: one summary row, expandable.

Collapsed, it is a single row (a toggle, then :func:`gratinglab.gui.provenance.summary`
in the panel's colours), so the plot above takes the height the full text used
to. Expanded, the full panel shows below that row. The four solver tabs build
one each in ``_build_provenance`` and fill it from ``_paint``.

Whether it is expanded is one choice for all four tabs, kept in the settings
store (:func:`gratinglab.gui.qt.theming.settings`) and applied to every open
panel when it changes. A run that failed opens its panel without changing that
choice, because a failure has to be read.
"""

from __future__ import annotations

import weakref

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QSplitter,
    QTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .. import provenance
from .theming import settings

__all__ = ["ProvenancePanel"]

_KEY = "panels/provenance_expanded"

#: Height the full text takes when expanded, if the splitter gave it less.
_EXPANDED_HEIGHT = 220

_PANELS: "weakref.WeakSet[ProvenancePanel]" = weakref.WeakSet()


class ProvenancePanel(QWidget):
    """A toggle and summary row over a read-only text panel."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 2, 8, 4)
        layout.setSpacing(2)

        row = QHBoxLayout()
        row.setSpacing(6)
        self.toggle = QToolButton()
        self.toggle.setText("Provenance")
        self.toggle.setCheckable(True)
        self.toggle.setAutoRaise(True)
        self.toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        # Bold by font, not stylesheet: a stylesheet freezes the palette at
        # polish time, and the label would keep the old theme's colour.
        bold = self.toggle.font()
        bold.setBold(True)
        self.toggle.setFont(bold)
        self.toggle.toggled.connect(self._on_toggled)
        row.addWidget(self.toggle)

        self.summary = QLabel()
        self.summary.setTextFormat(Qt.TextFormat.RichText)
        self.summary.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        row.addWidget(self.summary, 1)
        layout.addLayout(row)

        self.text = QTextEdit()
        self.text.setReadOnly(True)
        layout.addWidget(self.text)

        self._lines: tuple = ()
        self._colors: dict | None = None
        _PANELS.add(self)
        self.set_expanded(settings().value(_KEY, False, type=bool), remember=False)

    # -- content ---------------------------------------------------------

    def show_lines(self, lines, colors: "dict[str, str] | None" = None) -> None:
        """Fill the panel and its summary row. A failed run opens the panel."""
        self._lines = tuple(lines)
        self._colors = colors
        self.text.setHtml(provenance.to_html(self._lines, colors))
        self.text.verticalScrollBar().setValue(0)
        row = provenance.summary(self._lines)
        self.summary.setText(
            provenance.to_html(row, colors).replace("white-space:pre-wrap", "white-space:pre")
        )
        self.summary.setToolTip("".join(line.text for line in row))
        if any(line.tag == "bad" for line in row) and not self.is_expanded():
            self.set_expanded(True, remember=False)

    # -- expanding -------------------------------------------------------

    def is_expanded(self) -> bool:
        return self.toggle.isChecked()

    def set_expanded(self, expanded: bool, *, remember: bool = True) -> None:
        """Show or hide the full text. ``remember`` makes it the choice for every panel."""
        self.toggle.blockSignals(True)
        self.toggle.setChecked(expanded)
        self.toggle.blockSignals(False)
        self._apply(expanded)
        if remember:
            settings().setValue(_KEY, expanded)
            for other in list(_PANELS):
                if other is not self:
                    other.set_expanded(expanded, remember=False)

    def _on_toggled(self, expanded: bool) -> None:
        self.set_expanded(expanded)

    def _apply(self, expanded: bool) -> None:
        self.toggle.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self.text.setVisible(expanded)
        # Open, the panel's own first row says what the summary would.
        self.summary.setVisible(not expanded)
        if expanded:
            self.text.verticalScrollBar().setValue(0)
            self.setMaximumHeight(16777215)
            splitter = self.parentWidget()
            if isinstance(splitter, QSplitter):
                sizes = splitter.sizes()
                index = splitter.indexOf(self)
                total = sum(sizes)
                if index >= 0 and total and sizes[index] < _EXPANDED_HEIGHT:
                    sizes = [total - _EXPANDED_HEIGHT if i != index else _EXPANDED_HEIGHT
                             for i in range(len(sizes))]
                    splitter.setSizes(sizes)
        else:
            margins = self.layout().contentsMargins()
            self.setMaximumHeight(
                self.toggle.sizeHint().height() + margins.top() + margins.bottom()
            )
