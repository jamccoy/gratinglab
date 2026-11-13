"""The last folder per kind of file, and window layout: :mod:`endstation.qt.recent`.

Kept as a module so the windows and dialogs import it from here. It writes to
the store :mod:`gratinglab.gui.qt.theming` names, so ``GRATINGLAB_SETTINGS``
redirects it with the theme.
"""

from __future__ import annotations

from endstation.qt.recent import remember, restore_layout, save_layout, start_in

from . import theming as _theming  # noqa: F401 - names the settings store

__all__ = ["start_in", "remember", "restore_layout", "save_layout"]
