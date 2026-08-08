"""Console-script entry point.

The rest of the window lives in :mod:`gratinglab.gui.qt.main_window` and its
siblings. This module stays separate, and its PySide6 import stays deferred
into :func:`_require_qt`, so that ``import gratinglab.gui.app`` -- and by
extension ``import gratinglab`` -- is safe on an install that never asked for
the ``gui`` extra. The rest of :mod:`gratinglab.gui.qt` gets that same
property from never being imported except from in here or from a test that
already guarded for PySide6's absence.

Run with ``gratinglab-gui``.
"""

from __future__ import annotations

import sys

_QT_HELP = """\
gratinglab-gui needs PySide6, which this install does not have.

  pip install -e ".[gui]"        (from a checkout)
  pip install gratinglab[gui]    (once published)

Everything except the GUI works without it.
"""


def _require_qt():
    """Import PySide6 with an explanation instead of a traceback."""
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError as exc:  # pragma: no cover - depends on the install
        raise SystemExit(_QT_HELP) from exc
    return QApplication


def _bring_to_front(window) -> None:
    """Put the window in front on launch.

    A window started from a shell stub opens *behind* whatever is frontmost on
    macOS, because the process has no bundle identity of its own. `raise_()` +
    `activateWindow()` is the ordinary Qt remedy; the `osascript` call beneath
    it is the same belt-and-braces fallback the Tk window used, kept for now
    and worth re-testing once the app bundle carries real identity.
    """
    window.raise_()
    window.activateWindow()
    try:
        import os
        from subprocess import run

        run(
            ["osascript", "-e",
             'tell application "System Events" to set frontmost of the first '
             f'process whose unix id is {os.getpid()} to true'],
            capture_output=True, timeout=3,
        )
    except Exception:  # pragma: no cover - best effort only
        pass


def main() -> int:
    """Console-script entry point."""
    QApplication = _require_qt()
    from .qt.main_window import MainWindow

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    _bring_to_front(window)
    return app.exec()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
