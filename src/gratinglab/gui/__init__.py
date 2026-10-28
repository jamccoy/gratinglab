"""The desktop layer: two windows and the chooser that opens them.

Only :mod:`gratinglab.gui.qt` may import a toolkit. This module stays
toolkit-free so the console-script entry points, which live beside it, can
explain a missing PySide6 before anything tries to import it.
"""

from __future__ import annotations

from endstation.optional import MissingExtra, require_extra

__all__ = ["STILL_WORKS", "exit_without_qt"]

#: What a missing toolkit leaves working. Without it, a user who wanted only
#: the solvers reads an absent GUI extra as a broken install.
STILL_WORKS = "Everything except the GUI works without it."


def exit_without_qt() -> None:
    """Exit with an install instruction if PySide6 is not installed.

    The check is :func:`endstation.optional.require_extra`, a ``find_spec``
    call rather than an import. Its :class:`~endstation.optional.MissingExtra`
    becomes ``SystemExit`` here, so a console script prints the message and
    not a traceback naming a package the user never asked for.
    """
    try:
        require_extra("PySide6", "gui", still_works=STILL_WORKS)
    except MissingExtra as exc:
        raise SystemExit(str(exc)) from exc
