"""Console-script entry point.

The rest of the window lives in :mod:`gratinglab.gui.qt.main_window` and its
siblings. This module stays separate, and its PySide6 import stays deferred
into :func:`_require_qt`, so that ``import gratinglab.gui.app`` (and by
extension ``import gratinglab``) is safe on an install that never asked for
the ``gui`` extra. The rest of :mod:`gratinglab.gui.qt` has the same
property because it is imported only from this module or from a test that
has already guarded against PySide6 being absent.

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
    from the Tk window. M10-H measured its effect and found none:

    ============================  =============  ==============
    Launched from                 with the hack  without it
    ============================  =============  ==============
    ``GratingLab.app``            front          front
    ``gratinglab-gui`` in a shell not front      not front
    ============================  =============  ==============

    The bundle brings the window forward. It carries a real
    `CFBundleIdentifier`, so macOS treats the process as an application
    named "GratingLab" instead of "python", and LaunchServices activates it.
    Launched from a bare shell, the process has no such identity and stays
    behind. The AppleScript does not bring it forward either, because driving
    System Events requires Accessibility permission, which the launching
    terminal has usually not been granted.

    The fallback was therefore a `subprocess.run` on every launch, wrapped in
    a bare `except`, which did nothing in the one case it was written for, and
    it has been removed. The supported way to bring the window forward when
    starting from a shell is to launch the bundle, which `tools/apps.py bundle`
    builds and the README describes.
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
