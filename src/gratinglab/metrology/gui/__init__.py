"""
A PySide6 front-end for the blaze-angle analysis.

**This package carries no analysis logic.** It computes nothing itself: every
number it shows comes from :func:`afm_analysis.analyzer.analyze_single_file`,
which knows nothing about Qt, so behaviour stays testable without a display.

Two rules make that boundary real rather than aspirational:

1. **Only `gui.qt` may import a toolkit.** Everything else here -- `state.py` in
   particular -- is pure and tested headlessly. `tests/test_gui_boundary.py`
   walks this package and fails if PySide6 appears outside `gui/qt`.

2. **Only the main thread may touch a matplotlib figure.** The worker returns
   the result dictionary; the window draws from it. Drawing from a worker is the
   classic way a Qt + matplotlib application crashes, usually intermittently and
   far from the cause.

Importing `gui.qt` requires PySide6, which is an optional extra::

    pip install -e '.[gui]'
"""
from __future__ import annotations

__all__ = ["main", "require_qt", "QT_MISSING_MESSAGE"]

QT_MISSING_MESSAGE = (
    "the Qt interface needs PySide6, which is an optional extra.\n"
    "Install it with:\n\n"
    "    pip install -e '.[gui]'\n\n"
    "The analysis, the .ggp export and the ICC diagnostic all work without it."
)


def require_qt() -> None:
    """Raise something useful rather than a bare ImportError.

    ``ModuleNotFoundError: No module named 'PySide6'`` tells a user nothing about
    what to do, and this is the single most likely first-run failure.

    Uses find_spec rather than importing: this module sits outside `gui/qt`, and
    the rule that only `gui/qt` touches a toolkit is enforced by
    `tests/test_gui_boundary.py`. Checking availability should not cost a full
    toolkit import either - the point is to fail fast and clearly.
    """
    from importlib.util import find_spec

    if find_spec("PySide6") is None:  # pragma: no cover - depends on environment
        raise ModuleNotFoundError(QT_MISSING_MESSAGE)


def main(argv=None):
    """Launch the window. Entry point for the ``afm-gui`` console script."""
    require_qt()
    from .qt.app import main as _main
    return _main(argv)
