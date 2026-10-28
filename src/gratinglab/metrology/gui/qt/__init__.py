"""The widget layer, and the only part of the GUI allowed to import a toolkit.

That boundary is drawn as a package boundary so that a test can enforce it:
`tests/test_gui_boundary.py` walks every module under `gui/` outside this
subpackage and fails if PySide6 appears. The test guards against ordinary
drift, such as someone needing "just a QColor" in a module that was pure.

Everything here should be thin. If a question can be answered without a window,
it belongs in `gui/state.py`, which is pure and tested headlessly.
"""
from __future__ import annotations

# matplotlib's Qt shim picks a binding from sys.modules first and QT_API second.
# Without QT_API, a stray PyQt5 in the environment would win and its canvas
# would refuse to parent into a PySide6 window. This project shipped PyQt5
# until recently, so that stray binding is a realistic thing to find installed.
# Importing endstation.qt sets it, with setdefault; the modules here run
# `from . import *` first thing to get this import, so it exports nothing.
import endstation.qt  # noqa: F401

__all__: list[str] = []
