"""
A PySide6 front-end for the blaze-angle analysis.

This package carries no analysis logic and computes nothing itself. Every
number it shows comes from
:func:`gratinglab.metrology.analyzer.analyze_single_file`, which knows nothing
about Qt, so the behaviour stays testable without a display.

Two rules enforce that boundary:

1. **Only `gui.qt` may import a toolkit.** Everything else here, `state.py` in
   particular, is pure and tested headlessly. `tests/test_gui_boundary.py`
   walks this package and fails if PySide6 appears outside `gui/qt`.

2. **Only the main thread may touch a matplotlib figure.** The worker returns
   the result dictionary, and the window draws from it. Drawing from a worker
   is a well-known cause of crashes in Qt + matplotlib applications, which
   usually occur intermittently and far from the cause.

Importing `gui.qt` requires PySide6, which is an optional extra::

    pip install -e '.[gui]'
"""
from __future__ import annotations

from endstation.optional import missing_extra_message, require_extra

__all__ = ["main", "require_qt", "QT_MISSING_MESSAGE"]

_STILL_WORKS = (
    "The analysis, the .ggp export and the ICC diagnostic all work without it."
)

QT_MISSING_MESSAGE = missing_extra_message("PySide6", "gui", still_works=_STILL_WORKS)


def require_qt() -> None:
    """Raise something useful rather than a bare ImportError.

    ``ModuleNotFoundError: No module named 'PySide6'`` tells a user nothing about
    what to do, and this is the most likely first-run failure.

    The check is :func:`endstation.optional.require_extra`, a ``find_spec``
    call rather than an import, because this module sits outside `gui/qt` and
    the rule that only `gui/qt` touches a toolkit is enforced by
    `tests/test_gui_boundary.py`. It raises
    :class:`~endstation.optional.MissingExtra`, a ``ModuleNotFoundError``.
    """
    require_extra("PySide6", "gui", still_works=_STILL_WORKS)


def main(argv=None):
    """Launch the window. Entry point for the ``gratinglab-metrology-gui`` console script."""
    require_qt()
    from .qt.app import main as _main
    return _main(argv)
