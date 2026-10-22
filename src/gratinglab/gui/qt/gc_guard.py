"""Run Python's cycle collector on the GUI thread, never on a worker's.

Both windows do their heavy work on a ``QThread``: the metrology window's
analysis and the main window's solves. Python's automatic collector runs on
whichever thread happens to allocate past its threshold, so it regularly ran
*inside* those workers, and when the garbage it found held Qt objects (a
closed window, a figure from an earlier preview) it destroyed them off the GUI
thread. Qt does not allow that, and on macOS it segfaults: CI run 39835186938
died with ``Fatal Python error: Segmentation fault`` and "Garbage-collecting"
at the top of the worker thread's stack, in the middle of an ordinary
``genfromtxt``. The same thing can happen to a user: close some windows, press
Run, and the application crashes.

The remedy is the usual one for PySide and PyQt applications. While any window
holds the guard, automatic collection is off, and a timer on the GUI thread
runs the collector whenever the generation counts say it is due, so garbage is
still collected promptly and always on the one thread allowed to free a Qt
object. When the last window releases it, automatic collection is restored,
so code that runs with no window open behaves exactly as before.
"""

from __future__ import annotations

import gc

from PySide6.QtCore import QObject, QTimer
from PySide6.QtWidgets import QApplication

__all__ = ["guard", "is_active"]

#: How often the GUI thread checks whether a collection is due. Checking is
#: three integer comparisons; collecting happens only past the thresholds.
INTERVAL_MS = 500

_collector: "_Collector | None" = None
_holders = 0


class _Collector(QObject):
    def __init__(self, parent: QObject) -> None:
        super().__init__(parent)
        self._timer = QTimer(self)
        self._timer.timeout.connect(self.collect_if_due)
        self._timer.start(INTERVAL_MS)

    def collect_if_due(self) -> None:
        """What the automatic collector would have done, on this thread."""
        count0, count1, count2 = gc.get_count()
        limit0, limit1, limit2 = gc.get_threshold()
        if count0 > limit0:
            generation = 2 if count2 > limit2 else 1 if count1 > limit1 else 0
            gc.collect(generation)


def guard(window: QObject) -> None:
    """Keep the collector on the GUI thread for as long as ``window`` exists.

    Call from the window's constructor, on the GUI thread. Each call is
    released by that window's ``destroyed`` signal; the last release restores
    automatic collection.
    """
    global _collector, _holders
    if _holders == 0:
        app = QApplication.instance()
        # Whatever garbage already exists is collected here and now, on this
        # thread, before any worker can start allocating.
        gc.collect()
        gc.disable()
        _collector = _Collector(app if app is not None else window)
    _holders += 1
    window.destroyed.connect(_release)


def _release(*_args) -> None:
    global _collector, _holders
    _holders = max(0, _holders - 1)
    if _holders == 0 and _collector is not None:
        _collector._timer.stop()
        _collector.deleteLater()
        _collector = None
        gc.enable()


def is_active() -> bool:
    """True while some window holds the guard, so automatic collection is off."""
    return _holders > 0
