"""The chooser's widget. Which windows exist is :mod:`gratinglab.gui.launcher`.

Thin, by the rule for this package: it draws a button per mode and imports the
chosen window. Two pieces of real behaviour: keeping a reference to each
window it opens (see :data:`_OPEN_WINDOWS`), and importing both windows off
the GUI thread while the chooser is on screen (see :func:`prewarm`, which
exists because the first open used to freeze the application).

It also carries the View and Help menus the windows carry, built by the
same code. The chooser is the main menu: it stays open while the windows it
opened are up, as Gwiddion's toolbox does, and closing it closes them all.
There is one window per mode, so a second click on a card raises the window
already open rather than opening another.
"""

from __future__ import annotations

import atexit
import sys
import threading
from functools import partial
from importlib import import_module

import shiboken6
from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QMenuBar,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..launcher import MODES
from . import theming

__all__ = [
    "LauncherWindow",
    "SWITCH_MENU_TITLE",
    "install_switch_menu",
    "open_window",
    "prewarm",
    "show_chooser",
    "window_for",
]

#: Title of the menu each window gets, through which the other mode and the
#: chooser are reachable. Not "Window": macOS has conventions for a menu of
#: that name and these are modes, not windows of one document.
SWITCH_MENU_TITLE = "&Go"

#: Windows opened from the chooser, kept alive deliberately.
#:
#: A top-level Qt window with no Python reference is garbage collected and
#: vanishes a moment after it appears, which reads as "the button did nothing".
#: The list is not on the chooser, because a window opened from a Go menu,
#: in a session started by a direct script, has no chooser to hold it.
_OPEN_WINDOWS: list[QWidget] = []


#: The toolbox's width, and its buttons' height. Wide enough for the
#: sentence to wrap to three lines and the longest title to sit in its
#: button with room either side; a toolbox wider than this crowds the
#: windows it sits beside.
_TOOLBOX_WIDTH = 260
_BUTTON_HEIGHT = 34

#: The chooser, kept so that "Main Menu" shows the same one rather than a
#: second copy. Closing a window does not destroy it while a reference is
#: held, so this one can be shown again after it was closed.
_CHOOSER: list = []

#: The background import started by :func:`prewarm`, at most one per process.
_PREWARM: list[threading.Thread] = []

#: How long exit waits for that import to finish. It is seconds of work, and
#: waiting is only to keep the interpreter from finalising mid-import; a
#: thread that is somehow stuck must not hold the application open.
_PREWARM_JOIN_S = 10.0


def window_for(chosen):
    """The window already open for a mode, or None.

    Opening a second window for a mode would be worse than useless: the
    metrology window holds a loaded scan and its settings, so a "switch" that
    quietly produced an empty second copy would look like losing the work.
    """
    # A window Qt has deleted is dropped here, so the next open builds a
    # fresh one instead of touching a dead C++ object.
    _OPEN_WINDOWS[:] = [w for w in _OPEN_WINDOWS if shiboken6.isValid(w)]
    for window in _OPEN_WINDOWS:
        if getattr(window, "_gratinglab_mode", None) == chosen.key:
            return window
    return None


def _raise(window):
    """Show a window that may be open behind something, or closed."""
    window.show()
    window.raise_()
    window.activateWindow()
    return window


def install_switch_menu(window, current):
    """Give a window a menu leading to the other mode and to the chooser.

    Added from here rather than by each window, so that neither window needs
    to know the other exists and a third mode is still one entry in MODES.

    The menu and its actions are held on the window. `gui/qt/main_window.py`
    records why: under PySide6 a QMenu reached back through
    ``menuBar().actions()[i].menu()`` is not reliably usable, so whatever
    builds a menu must keep it.
    """
    menu = window.menuBar().addMenu(SWITCH_MENU_TITLE)
    window._gratinglab_switch_menu = menu

    actions = []
    for other in MODES:
        if other.key == current.key:
            continue
        action = menu.addAction(other.title)
        action.triggered.connect(partial(open_window, other))
        actions.append(action)

    menu.addSeparator()
    back = menu.addAction("Main Menu")
    back.setShortcut(QKeySequence("Ctrl+Shift+M"))
    back.triggered.connect(show_chooser)
    actions.append(back)

    window._gratinglab_switch_actions = actions
    return menu


def open_window(chosen):
    """Show the window a mode names, reusing it if it is already open."""
    existing = window_for(chosen)
    if existing is not None:
        return _raise(existing)

    module_name, _, attribute = chosen.window.partition(":")

    # A wait cursor around the part that cannot leave this thread. `prewarm`
    # should have made the import free, but it has not necessarily finished
    # when the user clicks, and the window itself still has to be built here.
    QApplication.setOverrideCursor(Qt.WaitCursor)
    try:
        window_class = getattr(import_module(module_name), attribute)
        window = window_class()
    finally:
        QApplication.restoreOverrideCursor()

    window._gratinglab_mode = chosen.key
    install_switch_menu(window, chosen)
    _OPEN_WINDOWS.append(window)
    return _raise(window)


def prewarm():
    """Import the window modules off the GUI thread, while the chooser waits.

    Opening a window for the first time used to freeze the application, badly
    enough on Windows to be reported as "almost non-responsive". The work is
    all import: measured on macOS, opening the metrology window blocks the GUI
    thread for 1.7 s with everything cached and 8.8 s on a machine that has
    not built matplotlib's font cache yet -- which is once per machine, and so
    is exactly the first open. The efficiency window costs 6.6 s on that same
    cold cache. Windows adds a multiplier these numbers do not show, because
    Defender scans each of scipy's and matplotlib's binaries the first time it
    is touched, and it paints "Not Responding" over any window that stops
    pumping events for about five seconds.

    None of that work needs the GUI thread. Importing a module binds classes
    and creates no QObject, so it runs here instead, in the seconds the user
    spends reading the chooser. The same measurement with this in place: the
    background import takes 7.6 s with the GUI live throughout, and the click
    that follows costs 0.2 s where it cost 8.8 s.

    A click part-way through is no worse than before: the GUI thread blocks on
    Python's per-module import lock, which is the same wait it used to do
    itself. The thread is joined at exit rather than merely left daemonic, so
    the interpreter cannot begin finalising in the middle of an import.

    Returns the thread, so a test can wait for it, or None when there is
    nothing left to import.
    """
    if _PREWARM:
        return None
    pending = [
        name
        for name in (chosen.window.partition(":")[0] for chosen in MODES)
        if name not in sys.modules
    ]
    if not pending:
        return None

    def warm():
        for module_name in pending:
            try:
                import_module(module_name)
            except Exception:
                # This is an optimisation and nothing more. A module that
                # cannot be imported has to raise where the user asked for
                # it -- from the click -- and not from a thread nobody
                # started, as a traceback on a console they may not have.
                # `continue`, not `return`: one broken window is no reason to
                # leave the other one cold.
                continue

    thread = threading.Thread(target=warm, name="gratinglab-prewarm", daemon=True)
    _PREWARM.append(thread)
    thread.start()
    atexit.register(thread.join, _PREWARM_JOIN_S)
    return thread


def show_chooser():
    """Bring the chooser back, the same one if it has been up before."""
    if not _CHOOSER:
        _CHOOSER.append(LauncherWindow())
    shown = _raise(_CHOOSER[0])
    # After the chooser is up, never before: the point is to use the time the
    # user spends reading it, and an import started first would delay it.
    prewarm()
    return shown


def _either(items):
    """``a``, ``a or b``, ``a, b, or c``: a list read as alternatives."""
    items = list(items)
    if len(items) < 3:
        return " or ".join(items)
    return ", ".join(items[:-1]) + ", or " + items[-1]


class LauncherWindow(QWidget):
    """A compact toolbox: one sentence, and a button per window.

    Kept small because it stays open beside the windows, as Gwiddion's
    toolbox does. What each window is for, and the command that opens it
    directly, are in its button's tooltip; the sentence names every window
    in a line, so none of them, the AFM half least of all, can be missed.

    It carries the same View and Help menus the windows do, since the theme
    is a setting for the whole application and not for whichever window
    happens to be open.

    It is the main menu, and stays open while windows are: see `choose` and
    `closeEvent`.

    A plain ``QWidget`` has no ``menuBar()``, so the bar is built here and put
    at the top of the layout. `install_help_menu` takes one for exactly this
    reason.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("gratinglab")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 8, 12, 12)
        layout.setSpacing(6)
        layout.setMenuBar(self._build_menu())

        opening = _either([m.summary for m in MODES if m.summary])
        # Only the first letter: str.capitalize() would lower "AFM" too.
        subtitle = QLabel(opening[:1].upper() + opening[1:] + ".")
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        subtitle.setToolTip(
            "The Go menu in each window switches to another window or back "
            "to this one, and closing this one closes them all."
        )
        layout.addWidget(subtitle)
        layout.addSpacing(4)

        self.buttons = {}
        for chosen in MODES:
            layout.addWidget(self._button(chosen))

        layout.addStretch(1)
        self.setFixedWidth(_TOOLBOX_WIDTH)

    def _build_menu(self):
        """View > Theme and Help, held on the instance and returned.

        Held for the reason the windows hold theirs: under PySide6 a QMenu
        reached back through the menu bar is not reliably usable, so whatever
        builds one keeps it. Returned because the caller has to hand the bar
        to the layout, which is how a QWidget gets a menu bar at all.
        """
        from .help_menu import install_help_menu

        # Install the manager here as well as in `gui.launcher.main`, for the
        # same reason the efficiency window does: a chooser built any other
        # way -- a test, an embedding, a future entry point -- would otherwise
        # show a Theme menu that silently did nothing. It is idempotent, so
        # the normal path still installs exactly one.
        self._themes = theming.install(QApplication.instance())

        self._menu_bar = QMenuBar(self)
        self._view_menu = self._menu_bar.addMenu("&View")
        self._theme_menu = theming.ThemeMenu(self)
        self._view_menu.addMenu(self._theme_menu)
        self._help_menu = install_help_menu(self, self._menu_bar)
        return self._menu_bar

    def _button(self, chosen):
        """One mode's button, with what it is for and its command as a tooltip."""
        button = QPushButton(chosen.title)
        button.setDefault(chosen is MODES[0])
        button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button.setMinimumHeight(_BUTTON_HEIGHT)
        button.setToolTip(
            f"<p>{chosen.blurb}</p>"
            f"<p>To open it directly, run <code>{chosen.script}</code>.</p>"
        )
        button.clicked.connect(partial(self.choose, chosen))
        button.setObjectName(f"mode-{chosen.key}")
        self.buttons[chosen.key] = button
        return button

    def choose(self, chosen):
        """Open a mode's window, or raise it if it is open, and stay open.

        The main menu stays up so it is where any window is reopened from,
        as Gwiddion's toolbox is. A window covering it is brought back with
        Go > Main Menu.
        """
        return open_window(chosen)

    def closeEvent(self, event):  # noqa: N802 - Qt's spelling
        """Closing the main menu closes every window it opened.

        Each window gets its own close, so it saves its layout and stops a
        job in flight as it would if closed by hand. A window that refuses
        keeps the main menu open too, so the menu is never gone while a
        window is still up. Nothing calls quit: once the last window is
        closed, Qt's quit-on-last-window ends the application, which also
        leaves a test's other widgets alone.
        """
        for window in list(_OPEN_WINDOWS):
            if not shiboken6.isValid(window) or not window.isVisible():
                continue
            if not window.close():
                event.ignore()
                return
        super().closeEvent(event)
