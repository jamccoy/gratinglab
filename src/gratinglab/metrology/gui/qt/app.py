"""
Application entry point for the Qt front-end.

Reached via ``afm-gui`` or ``python -m afm_analysis.gui.qt.app``.
"""
from __future__ import annotations

import sys

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow

__all__ = ["main"]


def main(argv=None):
    app = QApplication.instance() or QApplication(argv or sys.argv)
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
