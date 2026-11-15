"""The geometry window: a grating to look at, with no solver attached.

The chooser's "Explore a geometry" card opens it. It holds the same geometry
form and the same Grating Geometry tab the efficiency window does, and
nothing else: no solver tab, no worker thread, no solve. So it is a
reference to turn a grating over in, and a place to settle a geometry before
spending a solve on it.

It owns its form. Editing here changes nothing in the efficiency window
until "Send to efficiency window" copies the form across, which is the one
link between the two and goes one way. The efficiency window keeps its own
Grating Geometry tab, drawn from its own form.

Reads and writes: the form, the tab, and the window layout in the settings
store under ``"geometry"``. Called by: :mod:`gratinglab.gui.qt.launcher`,
through the ``geometry`` entry in :data:`gratinglab.gui.launcher.MODES`.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ..state import build, validate
from . import recent, theming
from .geometry_panel import GeometryPanel
from .geometry_tab import GeometryTab
from .sizing import fit_width_to_contents

__all__ = ["GeometryWindow"]

#: The efficiency window's debounce, for the same reason: one redraw per
#: pause in typing, not one per keystroke.
_REDRAW_DEBOUNCE_MS = 120

#: The key the window layout is saved under, and this window's mode key.
_LAYOUT_KEY = "geometry"


class GeometryWindow(QMainWindow):
    """A geometry form and the drawings of it, with a way to the solvers."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("gratinglab: geometry")
        self.resize(1180, 820)
        # Before any widget exists, so every figure is born in the theme.
        self._themes = theming.install(QApplication.instance())

        self._redraw_timer = QTimer(self)
        self._redraw_timer.setSingleShot(True)
        self._redraw_timer.setInterval(_REDRAW_DEBOUNCE_MS)
        self._redraw_timer.timeout.connect(self.refresh)

        self._build_dock()
        self.geometry_tab = GeometryTab()
        self.setCentralWidget(self.geometry_tab)
        self._build_menu()
        recent.restore_layout(self, _LAYOUT_KEY)
        self.refresh()

    # -- construction ----------------------------------------------------

    def _build_dock(self) -> None:
        """The form, with the send button under it, in a closable dock.

        The efficiency window's dock, for its reasons: Ctrl+G hands the
        whole width to the drawings. Enter in a field redraws at once, where
        in the efficiency window it would solve, since a redraw is all there
        is to do here.
        """
        self.geometry = GeometryPanel()
        self.geometry.solve_requested.connect(self.refresh)
        self.geometry.changed.connect(self._redraw_timer.start)

        self.send_button = QPushButton("Send to efficiency window")
        self.send_button.setToolTip(
            "Copy this geometry into the efficiency window's form, opening "
            "the window if it is not open. Nothing is solved until you ask."
        )
        self.send_button.clicked.connect(self.send_to_efficiency)

        column = QWidget()
        layout = QVBoxLayout(column)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.geometry)
        layout.addWidget(self.send_button)

        scroller = QScrollArea()
        scroller.setWidget(column)
        scroller.setWidgetResizable(True)
        scroller.setMinimumWidth(max(280, fit_width_to_contents(scroller)))

        dock = self.geometry_dock = QDockWidget("Grating geometry", self)
        dock.setObjectName("geometry_dock")
        dock.setWidget(scroller)
        dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )
        dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
        self.resizeDocks([dock], [300], Qt.Orientation.Horizontal)

    def _build_menu(self) -> None:
        """File, View and Help, each menu held on the window.

        Held for the reason `main_window.MainWindow._build_menu` gives: under
        PySide6 a QMenu reached back through the menu bar is not reliably
        usable.
        """
        from .help_menu import install_help_menu

        self._menu_bar = self.menuBar()

        self._file_menu = self._menu_bar.addMenu("&File")
        load = self._file_menu.addAction("&Open boundary profile…")
        load.setShortcut(QKeySequence.StandardKey.Open)
        load.triggered.connect(self.geometry.load_profile)
        send = self._file_menu.addAction("&Send to efficiency window")
        send.setShortcut(QKeySequence("Ctrl+E"))
        send.triggered.connect(self.send_to_efficiency)
        self._file_actions = {"open": load, "send": send}

        self._view_menu = self._menu_bar.addMenu("&View")
        toggle = self.geometry_dock.toggleViewAction()
        toggle.setText("&Geometry inputs")
        toggle.setShortcut(QKeySequence("Ctrl+G"))
        self._view_menu.addAction(toggle)
        self._view_menu.addSeparator()
        self._theme_menu = theming.ThemeMenu(self)
        self._view_menu.addMenu(self._theme_menu)

        self._help_menu = install_help_menu(self, self._menu_bar)

    # -- actions ---------------------------------------------------------

    def set_theme(self, choice: str) -> None:
        """Switch theme by key, or ``"system"``; the menu calls this."""
        self._themes.set_choice(choice)

    def refresh(self) -> None:
        """Redraw from the form, or keep the last drawing while it does not parse."""
        form = self.geometry.read_form()
        errors = validate(form)
        if errors:
            self.geometry_tab.show_pending(errors)
            return
        self.geometry_tab.show_geometry(build(form))

    def send_to_efficiency(self):
        """Copy this geometry into the efficiency window, and return that window.

        A form that does not parse is not sent: it would replace a working
        geometry over there with one no solver can take. The window is the
        one the Go menu opens, so sending never leaves two of them open.
        """
        form = self.geometry.read_form()
        errors = validate(form)
        if errors:
            QMessageBox.warning(
                self,
                "Geometry not sent",
                "Fix the geometry before sending it:\n\n"
                + "\n".join(f"{e.field}: {e.message}" for e in errors),
            )
            return None

        from ..launcher import mode
        from .launcher import open_window

        window = open_window(mode("efficiency"))
        window.open_geometry(form)
        return window

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        """Save the layout, and stop a redraw from firing into a closed window."""
        recent.save_layout(self, _LAYOUT_KEY)
        self._redraw_timer.stop()
        super().closeEvent(event)
