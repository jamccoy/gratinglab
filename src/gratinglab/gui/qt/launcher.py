"""The chooser's widget. Which windows exist is :mod:`gratinglab.gui.launcher`.

Thin, by the rule for this package: it draws a button per mode and imports the
chosen window. The one piece of real behaviour is keeping a reference to each
window it opens -- see :data:`_OPEN_WINDOWS`.
"""

from __future__ import annotations

from functools import partial
from importlib import import_module

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence
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

__all__ = [
    "LauncherWindow",
    "SWITCH_MENU_TITLE",
    "install_switch_menu",
    "open_window",
    "show_chooser",
    "window_for",
]

#: Title of the menu each window gets, through which the other mode and the
#: chooser are reachable. Not "Window": macOS has conventions for a menu of
#: that name and these are modes, not windows of one document.
SWITCH_MENU_TITLE = "&Go"

#: Windows opened from the chooser, kept alive deliberately.
#:
#: A top-level Qt window with no Python reference is garbage collected and
#: vanishes a moment after it appears, which reads as "the button did nothing".
#: The list cannot live on the chooser, because the chooser closes itself as
#: soon as a window is up and would take the only reference with it.
_OPEN_WINDOWS: list[QWidget] = []


#: The chooser, kept so that "Main Menu" shows the same one rather than a
#: second copy. Closing a window does not destroy it while a reference is
#: held, so this one can be shown again after it closed itself.
_CHOOSER: list = []


def window_for(chosen):
    """The window already open for a mode, or None.

    Opening a second window for a mode would be worse than useless: the
    metrology window holds a loaded scan and its settings, so a "switch" that
    quietly produced an empty second copy would look like losing the work.
    """
    for window in _OPEN_WINDOWS:
        if getattr(window, "_gratinglab_mode", None) == chosen.key:
            return window
    return None


def _raise(window):
    """Show a window that may be open behind something, or closed."""
    window.show()
    window.raise_()
    window.activateWindow()
    return window


def install_switch_menu(window, current):
    """Give a window a menu leading to the other mode and to the chooser.

    Added from here rather than by each window, so that neither window needs
    to know the other exists and a third mode is still one entry in MODES.

    The menu and its actions are held on the window. `gui/qt/main_window.py`
    records why: under PySide6 a QMenu reached back through
    ``menuBar().actions()[i].menu()`` is not reliably usable, so whatever
    builds a menu must keep it.
    """
    menu = window.menuBar().addMenu(SWITCH_MENU_TITLE)
    window._gratinglab_switch_menu = menu

    actions = []
    for other in MODES:
        if other.key == current.key:
            continue
        action = menu.addAction(other.title)
        action.triggered.connect(partial(open_window, other))
        actions.append(action)

    menu.addSeparator()
    back = menu.addAction("Main Menu")
    back.setShortcut(QKeySequence("Ctrl+Shift+M"))
    back.triggered.connect(show_chooser)
    actions.append(back)

    window._gratinglab_switch_actions = actions
    return menu


def open_window(chosen):
    """Show the window a mode names, reusing it if it is already open."""
    existing = window_for(chosen)
    if existing is not None:
        return _raise(existing)

    module_name, _, attribute = chosen.window.partition(":")
    window_class = getattr(import_module(module_name), attribute)

    window = window_class()
    window._gratinglab_mode = chosen.key
    install_switch_menu(window, chosen)
    _OPEN_WINDOWS.append(window)
    return _raise(window)


def show_chooser():
    """Bring the chooser back, the same one if it has been up before."""
    if not _CHOOSER:
        _CHOOSER.append(LauncherWindow())
    return _raise(_CHOOSER[0])


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
            "Either window's <b>Go</b> menu switches to the other and comes "
            "back here, so this is not a one-way door. Each is also a command "
            "of its own: "
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
