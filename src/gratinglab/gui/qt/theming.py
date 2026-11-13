"""gratinglab's themes, applied by :mod:`endstation.qt.theming`.

The manager, the View > Theme menu, the palette and the settings store are
endstation's. This module names the store and passes gratinglab's themes
(:mod:`gratinglab.gui.theme`, which add the diagram and order colours), so
the rest of the package keeps importing from here.

The choice is stored with ``QSettings`` under ``gratinglab``.
``GRATINGLAB_SETTINGS`` names an INI file to use instead, which is what the
test suite sets, the same arrangement as ``GRATINGLAB_MATERIALS``: a run of
the tests must never write to the user's own settings.
"""

from __future__ import annotations

from endstation.qt import settings as _settings_store
from endstation.qt import theming as _theming
from endstation.qt.theming import ThemeManager, ThemeMenu, manager, palette, subscribe, swatch
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication

from ..theme import LIGHT, THEMES, Theme

__all__ = [
    "QIcon",
    "SETTINGS_ENV",
    "ThemeMenu",
    "ThemeManager",
    "current",
    "install",
    "manager",
    "palette",
    "settings",
    "subscribe",
    "swatch",
]

#: Names an INI file to keep settings in instead of the platform's store.
SETTINGS_ENV = "GRATINGLAB_SETTINGS"

# At import, so every entry point (both windows, the chooser, a test) reads
# and writes the same store. The organization and application names are the
# ones the store had before it moved to endstation, so saved choices carry over.
_settings_store.configure("gratinglab", "gratinglab", env=SETTINGS_ENV)

settings = _settings_store.settings
_settings = settings


def install(app: QApplication) -> ThemeManager:
    """The application's manager with gratinglab's themes, created on first call."""
    return _theming.install(app, THEMES)


def current() -> Theme:
    """The theme in effect, or :data:`gratinglab.gui.theme.LIGHT` before any is installed."""
    return _theming.current(LIGHT)
