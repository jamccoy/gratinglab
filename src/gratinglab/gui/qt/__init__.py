"""The widget layer, and the only part of the GUI allowed to import a toolkit.

The boundary is drawn as a package boundary, so that a test can enforce it:
`tests/test_gui_qt.py` walks every module under `gui/` outside this
subpackage and fails if PySide6 appears. It guards against the common case of
a toolkit import, such as a single QColor, being added to a module that was
pure.

Everything here should be thin. If a question can be answered without a
window, its answer belongs in `gui/state.py`, `gui/provenance.py`,
`gui/richtext.py` or `gui/orders.py`, all of which are pure and tested
headlessly.
"""

from __future__ import annotations

# matplotlib's Qt shim picks a binding from sys.modules first and QT_API
# second. Without QT_API, a stray PyQt5 in the environment would win and its
# canvas would refuse to parent into a PySide6 window. Importing endstation.qt
# sets it, with setdefault, before anything here imports matplotlib.
import endstation.qt  # noqa: F401
