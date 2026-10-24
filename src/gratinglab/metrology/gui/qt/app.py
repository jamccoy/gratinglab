"""
Application entry point for the Qt front-end.

Reached via ``gratinglab-metrology-gui`` or ``python -m gratinglab.metrology.gui.qt.app``.
"""
from __future__ import annotations

import sys

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow

__all__ = ["main"]


def main(argv=None):
    app = QApplication.instance() or QApplication(argv or sys.argv)
    # Fusion plus the saved theme, shared with the efficiency window.
    from gratinglab.gui.qt.theming import install

    install(app)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
