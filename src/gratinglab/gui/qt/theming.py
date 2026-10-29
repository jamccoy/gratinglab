"""Apply a :mod:`gratinglab.gui.theme` to the running application.

One :class:`ThemeManager` per ``QApplication``. It turns the chosen theme
into a ``QPalette`` (under Fusion, which is what every entry point already
asked for, so the palette governs every widget) and into matplotlib
settings, remembers the choice, and emits ``changed`` so that each figure
and each rendered page can restyle itself. A window that draws something
connects to ``changed`` through :func:`subscribe`.

``"system"`` follows the operating system's light or dark setting, live:
Qt reports a change of colour scheme and the manager re-applies.

The choice is stored with ``QSettings``. ``GRATINGLAB_SETTINGS`` names an
INI file to use instead, which is what the test suite sets, the same
arrangement as ``GRATINGLAB_MATERIALS``: a run of the tests must never
write to the user's own settings.
"""

from __future__ import annotations

import os
from collections.abc import Callable

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtGui import (QAction, QActionGroup, QColor, QGuiApplication,
                           QIcon, QPalette, QPixmap)
from PySide6.QtWidgets import QApplication, QMenu

from ..theme import DARK, LIGHT, SYSTEM, THEMES, Theme, rc_params

__all__ = [
    "QIcon",
    "SETTINGS_ENV",
    "ThemeMenu",
    "ThemeManager",
    "current",
    "install",
    "manager",
    "settings",
    "subscribe",
    "swatch",
]

#: Names an INI file to keep settings in instead of the platform's store.
SETTINGS_ENV = "GRATINGLAB_SETTINGS"

_KEY = "appearance/theme"
_PROPERTY = "_gratinglab_theme_manager"


def settings() -> QSettings:
    """The application's settings store, or the INI file ``GRATINGLAB_SETTINGS`` names.

    Shared by everything the windows remember (the theme, the AFM column
    width), so the test suite's redirect covers all of it at once.
    """
    path = os.environ.get(SETTINGS_ENV)
    if path:
        return QSettings(path, QSettings.Format.IniFormat)
    return QSettings("gratinglab", "gratinglab")


_settings = settings


def palette(theme: Theme) -> QPalette:
    """A ``QPalette`` for ``theme``, with disabled text dimmed to its muted colour."""
    role = QPalette.ColorRole
    group = QPalette.ColorGroup
    result = QPalette()
    entries = {
        role.Window: theme.window,
        role.WindowText: theme.text,
        role.Base: theme.base,
        role.AlternateBase: theme.alt_base,
        role.ToolTipBase: theme.base,
        role.ToolTipText: theme.text,
        role.PlaceholderText: theme.muted,
        role.Text: theme.text,
        role.Button: theme.button,
        role.ButtonText: theme.text,
        role.BrightText: theme.bad,
        role.Highlight: theme.highlight,
        role.HighlightedText: theme.highlight_text,
        role.Link: theme.link,
        role.LinkVisited: theme.link,
        role.Mid: theme.muted,
    }
    for which, color in entries.items():
        result.setColor(which, QColor(color))
    for which in (role.WindowText, role.Text, role.ButtonText):
        result.setColor(group.Disabled, which, QColor(theme.muted))
    # Light and Dark shading for Fusion's bevels, derived rather than listed.
    base = QColor(theme.button)
    result.setColor(role.Light, base.lighter(130))
    result.setColor(role.Midlight, base.lighter(115))
    result.setColor(role.Dark, base.darker(160))
    result.setColor(role.Shadow, base.darker(250))
    return result


class ThemeManager(QObject):
    """Holds the chosen theme and applies it to one application."""

    changed = Signal(object)  # the resolved Theme

    def __init__(self, app: QApplication) -> None:
        super().__init__(app)
        self._app = app
        stored = str(_settings().value(_KEY, SYSTEM))
        self._choice = stored if stored == SYSTEM or stored in THEMES else SYSTEM
        hints = QGuiApplication.styleHints()
        hints.colorSchemeChanged.connect(self._on_scheme_changed)

    @property
    def choice(self) -> str:
        """The stored choice: ``"system"`` or a key of :data:`THEMES`."""
        return self._choice

    @property
    def theme(self) -> Theme:
        """The theme in effect, with ``"system"`` resolved."""
        if self._choice != SYSTEM:
            return THEMES[self._choice]
        scheme = QGuiApplication.styleHints().colorScheme()
        return DARK if scheme == Qt.ColorScheme.Dark else LIGHT

    def set_choice(self, choice: str) -> None:
        """Choose a theme by key, or ``"system"``, and apply it."""
        if choice != SYSTEM and choice not in THEMES:
            raise ValueError(
                f"unknown theme {choice!r}; choose {SYSTEM!r} or one of "
                f"{', '.join(THEMES)}"
            )
        self._choice = choice
        settings = _settings()
        settings.setValue(_KEY, choice)
        settings.sync()
        self.apply()

    def apply(self) -> None:
        """Push the current theme into Qt and matplotlib, then announce it."""
        import matplotlib

        theme = self.theme
        self._app.setPalette(palette(theme))
        matplotlib.rcParams.update(rc_params(theme))
        self.changed.emit(theme)

    def _on_scheme_changed(self, _scheme) -> None:
        if self._choice == SYSTEM:
            self.apply()


def install(app: QApplication) -> ThemeManager:
    """The application's manager, created and applied on first call.

    Sets the Fusion style too, which every entry point already asked for:
    under the native styles a palette is only a suggestion.
    """
    existing = manager(app)
    if existing is not None:
        return existing
    app.setStyle("Fusion")
    created = ThemeManager(app)
    app.setProperty(_PROPERTY, created)
    created.apply()
    return created


def manager(app: QApplication | None = None) -> ThemeManager | None:
    """The manager installed on ``app`` (the running one by default), if any."""
    app = app or QApplication.instance()
    if app is None:
        return None
    found = app.property(_PROPERTY)
    return found if isinstance(found, ThemeManager) else None


def current() -> Theme:
    """The theme in effect, or :data:`LIGHT` before any is installed."""
    found = manager()
    return found.theme if found is not None else LIGHT


def subscribe(callback: Callable[[Theme], None]) -> None:
    """Call ``callback(theme)`` whenever the theme changes.

    A no-op when no manager is installed, so a widget built on its own (in a
    test, or embedded elsewhere) still works and simply keeps its colours.

    Pass a bound method of a ``QObject``, never a lambda: Qt drops a
    connection to a method when its object is destroyed, but keeps calling
    a lambda, into a widget that no longer exists.
    """
    found = manager()
    if found is not None:
        found.changed.connect(callback)


class ThemeMenu(QMenu):
    """The View > Theme submenu, one per window, all kept in step.

    Every window that has a menu bar carries one, built here so that none can
    be left without it: the AFM window had no theme choice at all while the
    efficiency window did. Choosing in one window re-checks the others,
    because each menu follows the manager's ``changed`` signal rather than
    its own clicks.
    """

    def __init__(self, parent=None) -> None:
        super().__init__("&Theme", parent)
        self._group = QActionGroup(self)
        self._group.setExclusive(True)
        choices = [(SYSTEM, "System")] + [(k, t.label) for k, t in THEMES.items()]
        for key, label in choices:
            action = self.addAction(label)
            action.setCheckable(True)
            action.setData(key)
            action.triggered.connect(self._chosen)
            self._group.addAction(action)
            if key == SYSTEM:
                self.addSeparator()
        found = manager()
        if found is not None:
            found.changed.connect(self._sync)
        self._sync()

    def actions_by_key(self) -> dict[str, QAction]:
        return {a.data(): a for a in self._group.actions()}

    def _chosen(self) -> None:
        action = self.sender()
        found = manager()
        if found is not None and isinstance(action, QAction):
            found.set_choice(action.data())

    def _sync(self, _theme=None) -> None:
        found = manager()
        choice = found.choice if found is not None else SYSTEM
        for key, action in self.actions_by_key().items():
            action.setChecked(key == choice)


def swatch(color: str, size: int = 12) -> QIcon:
    """A square of ``color``, for a list entry that names a plotted line."""
    pixmap = QPixmap(size, size)
    pixmap.fill(QColor(color))
    return QIcon(pixmap)
