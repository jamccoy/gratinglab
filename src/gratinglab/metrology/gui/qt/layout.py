"""
The controls-and-plot layout every AFM tab shares.

Each tab was a fixed 272 px column of controls in a 292 px scroll area beside
its plot, so long labels were crammed and nothing could be widened. A splitter
replaces that: the column starts wide enough for its contents, cannot be
dragged narrower than its widest control needs, and gives the plot every pixel
the window gains.
"""
from __future__ import annotations

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFrame, QScrollArea, QSplitter, QWidget

__all__ = ["controls_beside", "MIN_CONTROLS_WIDTH", "START_CONTROLS_WIDTH"]

#: Narrower than this and the group boxes' labels start to wrap mid-word.
MIN_CONTROLS_WIDTH = 260
#: Where the column starts, unless its contents ask for more.
START_CONTROLS_WIDTH = 320


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
    start = max(START_CONTROLS_WIDTH, floor)
    splitter.setSizes([start, max(1, 4 * start)])
    return splitter
