"""
Running the analysis off the UI thread.

A single scan takes about 0.3 s, so today this buys little a user would notice.
It matters as soon as the window grows a compare mode: eight scans is seconds,
and retrofitting an event flow around a window that has already frozen is worse
than building for it once.

**Cancel abandons a result; it does not interrupt one.** `analyze_single_file`
reaches NumPy-bound loops with no check point, and a Python thread cannot be
killed, so there is nothing to interrupt from out here. Cancelling bumps a
token: the in-flight result still arrives, is recognised as stale, and is
dropped. The window stops waiting; the CPU does not stop working.

The worker returns plain data. It never touches a figure - see the rule in
`afm_analysis.gui`.
"""
from __future__ import annotations

import contextlib
import io

from PySide6.QtCore import QObject, Signal, Slot

__all__ = ["AnalysisWorker"]


class AnalysisWorker(QObject):
    """Runs one analysis and emits the result dictionary."""

    finished = Signal(object, object, int)   # result, settings, token
    failed = Signal(str, int)                # message, token
    log = Signal(str)                        # captured stdout

    def __init__(self, parent=None):
        super().__init__(parent)

    @Slot(str, object, int)
    def run(self, filename, settings, token):
        """
        Analyse one file.

        `token` is echoed back untouched so the window can tell whether the
        result it receives is the one it is still waiting for.
        """
        from ...analyzer import analyze_single_file

        buffer = io.StringIO()
        try:
            with contextlib.redirect_stdout(buffer):
                result = analyze_single_file(filename, show_plots=False,
                                             settings=settings)
        except Exception as exc:  # surfaced in the window, not a traceback
            self.log.emit(buffer.getvalue())
            self.failed.emit(str(exc), token)
            return

        self.log.emit(buffer.getvalue())
        self.finished.emit(result, settings, token)
