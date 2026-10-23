"""The chooser's widget. Which windows exist is :mod:`gratinglab.gui.launcher`.

Thin, by the rule for this package: it draws a button per mode and imports the
chosen window. The one piece of real behaviour is keeping a reference to each
window it opens -- see :data:`_OPEN_WINDOWS`.
"""

from __future__ import annotations

from functools import partial
from importlib import import_module

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..launcher import MODES

__all__ = ["LauncherWindow", "open_window"]

#: Windows opened from the chooser, kept alive deliberately.
#:
#: A top-level Qt window with no Python reference is garbage collected and
#: vanishes a moment after it appears, which reads as "the button did nothing".
#: The list cannot live on the chooser, because the chooser closes itself as
#: soon as a window is up and would take the only reference with it.
_OPEN_WINDOWS: list[QWidget] = []


def open_window(chosen):
    """Import the window a mode names, show it, and keep it alive."""
    module_name, _, attribute = chosen.window.partition(":")
    window_class = getattr(import_module(module_name), attribute)

    window = window_class()
    _OPEN_WINDOWS.append(window)
    window.show()
    window.raise_()
    window.activateWindow()
    return window


class LauncherWindow(QWidget):
    """Two buttons and a sentence each, so neither half can be missed."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("GratingLab")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(28, 24, 28, 24)
        layout.setSpacing(18)

        heading = QLabel("GratingLab")
        font = heading.font()
        font.setPointSize(font.pointSize() + 8)
        font.setBold(True)
        heading.setFont(font)
        layout.addWidget(heading)

        subtitle = QLabel(
            "Groove metrology and diffraction efficiency. Both halves are this "
            "one package; pick where to start."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        layout.addWidget(subtitle)

        self.buttons = {}
        for chosen in MODES:
            layout.addWidget(self._card(chosen))

        layout.addStretch(1)

        hint = QLabel(
            "Each window is also a command of its own: "
            + ", ".join(f"<code>{m.script}</code>" for m in MODES)
            + "."
        )
        hint.setWordWrap(True)
        hint.setTextFormat(Qt.RichText)
        hint.setStyleSheet("color: palette(mid);")
        layout.addWidget(hint)

        self.setMinimumWidth(460)

    def _card(self, chosen):
        """One mode: its button, and the sentence saying what it is for."""
        frame = QFrame()
        frame.setFrameShape(QFrame.StyledPanel)
        inner = QVBoxLayout(frame)
        inner.setContentsMargins(16, 14, 16, 14)
        inner.setSpacing(6)

        row = QHBoxLayout()
        button = QPushButton(chosen.title)
        button.setDefault(chosen is MODES[0])
        button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button.clicked.connect(partial(self.choose, chosen))
        button.setObjectName(f"mode-{chosen.key}")
        row.addWidget(button)
        inner.addLayout(row)

        blurb = QLabel(chosen.blurb)
        blurb.setWordWrap(True)
        inner.addWidget(blurb)

        self.buttons[chosen.key] = button
        return frame

    def choose(self, chosen):
        """Open a mode's window, then close the chooser.

        In that order: Qt quits when the last window closes, so closing first
        would end the process before the window it was asked for exists.
        """
        window = open_window(chosen)
        self.close()
        return window
