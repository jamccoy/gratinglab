"""The widget layer, and the only part of the GUI allowed to import a toolkit.

That boundary is a package boundary rather than a convention so a test can
enforce it: `tests/test_gui_qt.py` walks every module under `gui/` outside this
subpackage and fails if PySide6 appears. The drift it prevents is the ordinary
kind -- someone needing "just a QColor" in a module that was pure.

Everything here should be thin. If a question can be answered without a window
it belongs in `gui/state.py`, `gui/provenance.py`, `gui/richtext.py` or
`gui/orders.py`, all of which are pure and tested headlessly.
"""

from __future__ import annotations

import os

# matplotlib's Qt shim picks a binding from sys.modules first and QT_API
# second. Without this, a stray PyQt5 in the environment would win and its
# canvas would refuse to parent into a PySide6 window.
os.environ.setdefault("QT_API", "PySide6")
