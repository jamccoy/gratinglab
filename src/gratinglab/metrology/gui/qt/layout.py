"""
The controls-and-plot layout every AFM tab shares.

Each tab was a fixed 272 px column of controls in a 292 px scroll area beside
its plot, so long labels were crammed and nothing could be widened. A splitter
replaces that: the column starts wide enough for its contents, cannot be
dragged narrower than its widest control needs, and gives the plot every pixel
the window gains.

The width is one setting, not one per tab. Dragging the column in any tab
moves it in the others, so switching tabs never makes the plot jump, and the
width is saved as it is dragged, so the window opens at it next time. It is
kept in the same store as the theme (``theming.settings``), which the test
suite redirects to a scratch file.
"""
from __future__ import annotations

import weakref

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFrame, QScrollArea, QSplitter, QWidget

from gratinglab.gui.qt.theming import settings

__all__ = ["controls_beside", "MIN_CONTROLS_WIDTH", "START_CONTROLS_WIDTH",
           "WIDTH_KEY"]

#: Narrower than this and the group boxes' labels start to wrap mid-word.
MIN_CONTROLS_WIDTH = 260
#: Where the column starts, unless its contents ask for more.
START_CONTROLS_WIDTH = 320
#: Where the chosen width is remembered.
WIDTH_KEY = "afm/controls_width"

# Every live splitter, so a drag in one tab can move the others. Weak, so a
# closed window's splitters drop out on their own.
_SPLITTERS: "weakref.WeakSet[QSplitter]" = weakref.WeakSet()


def _saved_width() -> int | None:
    value = settings().value(WIDTH_KEY)
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None     # a hand-edited or corrupt value: fall back to the default


def _follow(source: QSplitter) -> None:
    """Remember ``source``'s column width and give it to every other tab."""
    width = source.sizes()[0]
    store = settings()
    store.setValue(WIDTH_KEY, width)
    store.sync()
    for other in list(_SPLITTERS):
        if other is source:
            continue
        try:
            total = sum(other.sizes())
            other.blockSignals(True)
            other.setSizes([width, max(1, total - width)])
            other.blockSignals(False)
        except RuntimeError:
            # Its window was closed and Qt has deleted it, while the Python
            # wrapper waits for collection.
            _SPLITTERS.discard(other)


def controls_beside(controls: QWidget, plot: QWidget) -> QSplitter:
    """A splitter holding the scrolling ``controls`` column and the ``plot``."""
    scroll = QScrollArea()
    scroll.setWidget(controls)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.NoFrame)
    # Width is the splitter's to decide; only the height should ever scroll.
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    gutter = scroll.verticalScrollBar().sizeHint().width()
    floor = max(MIN_CONTROLS_WIDTH, controls.minimumSizeHint().width() + gutter)
    scroll.setMinimumWidth(floor)

    splitter = QSplitter(Qt.Horizontal)
    splitter.setChildrenCollapsible(False)
    splitter.addWidget(scroll)
    splitter.addWidget(plot)
    splitter.setStretchFactor(0, 0)
    splitter.setStretchFactor(1, 1)
    # A remembered width wins, but never below what this tab's controls need.
    start = max(_saved_width() or START_CONTROLS_WIDTH, floor)
    splitter.setSizes([start, max(1, 4 * start)])

    _SPLITTERS.add(splitter)
    splitter.splitterMoved.connect(lambda *_: _follow(splitter))
    return splitter
