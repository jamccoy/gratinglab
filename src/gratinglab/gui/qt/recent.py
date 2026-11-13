"""What the windows remember between sessions: the last folder per kind of file, and window layout.

Both live in the settings store from :func:`gratinglab.gui.qt.theming.settings`,
so ``GRATINGLAB_SETTINGS`` redirects them with the theme. File dialogs call
:func:`start_in` for their starting folder and :func:`remember` with the path
the user chose; the main windows call :func:`restore_layout` when built and
:func:`save_layout` when closed.

A folder is keyed by the kind of file (``"scans"``, ``"profiles"``,
``"tables"``, ``"exports"``), so opening a scan does not move where the next
export goes.
"""

from __future__ import annotations

from pathlib import Path

from .theming import settings

__all__ = ["start_in", "remember", "restore_layout", "save_layout"]


def start_in(kind: str, filename: str = "", fallback: str = "") -> str:
    """The folder last used for ``kind``, joined to ``filename`` if given.

    A remembered folder that no longer exists is ignored, and ``fallback`` is
    used instead (an empty string lets Qt pick its own default).
    """
    folder = settings().value(f"folders/{kind}", "", type=str)
    if not folder or not Path(folder).is_dir():
        folder = fallback
    if not filename:
        return folder
    return str(Path(folder) / filename) if folder else filename


def remember(kind: str, path: "str | Path") -> None:
    """Record the folder holding ``path`` as the one to start in for ``kind``."""
    if path:
        settings().setValue(f"folders/{kind}", str(Path(path).resolve().parent))


def restore_layout(window, key: str) -> bool:
    """Restore a main window's size, position and dock layout. True if any was saved."""
    store = settings()
    geometry = store.value(f"windows/{key}/geometry")
    state = store.value(f"windows/{key}/state")
    if geometry is not None:
        window.restoreGeometry(geometry)
    if state is not None:
        window.restoreState(state)
    return geometry is not None


def save_layout(window, key: str) -> None:
    """Save a main window's size, position and dock layout."""
    store = settings()
    store.setValue(f"windows/{key}/geometry", window.saveGeometry())
    store.setValue(f"windows/{key}/state", window.saveState())
