"""
Running the analysis off the UI thread.

A single scan takes about 0.3 s, so at present this gives little benefit a
user would notice. It becomes relevant once the window gains a compare mode,
where eight scans take seconds, and retrofitting an event flow around a window
that already freezes is harder than building for it once.

Cancel abandons a result without interrupting the computation.
`analyze_single_file` reaches NumPy-bound loops with no check point, and a
Python thread cannot be killed, so there is nothing to interrupt from out
here. Cancelling increments a token: the in-flight result still arrives, is
recognised as stale, and is dropped. The window stops waiting, while the CPU
continues the computation to completion.

The worker returns plain data and never touches a figure; see the rule in
`gratinglab.metrology.gui`.
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
