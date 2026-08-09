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

    There used to be an `osascript`/System Events fallback here, inherited
    from the Tk window. M10-H measured what it actually bought, and the answer
    was nothing:

    ============================  =============  ==============
    Launched from                 with the hack  without it
    ============================  =============  ==============
    ``GratingLab.app``            front          front
    ``gratinglab-gui`` in a shell not front      not front
    ============================  =============  ==============

    The bundle is what does the work: it carries a real `CFBundleIdentifier`,
    so macOS treats the process as an application ("GratingLab", not
    "python") and LaunchServices activates it. From a bare shell the process
    has no such identity and stays behind -- and the AppleScript does not
    rescue it, because driving System Events needs Accessibility permission
    the launching terminal has usually never been granted.

    So it was a `subprocess.run` on every launch, wrapped in a bare `except`,
    that silently no-opped in the one case it existed for. Removed rather than
    kept as a talisman. The supported fix for the shell path is to launch the
    bundle, which `tools/make_app.py` builds and the README points at.
    """
    window.raise_()
    window.activateWindow()


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
