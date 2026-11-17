"""The coatings window: a stack of materials, with no grooves in it.

The chooser's "Design a coating" card opens it. It answers what a coating
does before a grating is cut into it: the Materials tab names a material and
shows its index, and the Mirror tab reflects the stack as a flat mirror and
says how deep the light goes and whether the top layer may be taken as bulk.
No solver, no worker thread.

It owns its stack. The dock is the efficiency window's geometry panel showing
its Surface and Wavelengths groups alone, with the base material called the
Substrate. Editing here changes nothing in the efficiency window until "Send
to efficiency window" copies the stack across; the efficiency window's "Send
coating to coatings window" is the way back, and brings its mount's graze
angles with it for the Mirror tab to mark.

The material library is one for the whole application. Its actions are on the
Materials tab here, and every window follows
:data:`~gratinglab.gui.qt.library_bar.library_events`.

Reads and writes: the dock, the two tabs, and the window layout in the
settings store under ``"coatings"``. Called by:
:mod:`gratinglab.gui.qt.launcher`, through the ``coatings`` entry in
:data:`gratinglab.gui.launcher.MODES`.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QMainWindow,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from . import recent, theming
from .geometry_panel import GeometryPanel
from .library_bar import library_events
from .materials_tab import MaterialsTab
from .mirror_tab import MirrorTab
from .sizing import fit_width_to_contents

__all__ = ["CoatingsWindow"]

#: The key the window layout is saved under, and this window's mode key.
_LAYOUT_KEY = "coatings"


class CoatingsWindow(QMainWindow):
    """A coating stack, its materials, and the stack as a mirror."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("gratinglab: coatings")
        self.resize(1180, 820)
        # Before any widget exists, so every figure is born in the theme.
        self._themes = theming.install(QApplication.instance())

        self._build_dock()
        self.materials_tab = MaterialsTab(self.stack)
        self.mirror_tab = MirrorTab(self.stack)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.materials_tab, "Materials")
        self.tabs.addTab(self.mirror_tab, "Mirror")
        self.tabs.setTabToolTip(0, "Name a material, see its index, and put it in the stack")
        self.tabs.setTabToolTip(1, "The stack as a flat mirror: reflectivity and penetration depth")
        self.setCentralWidget(self.tabs)
        self._build_menu()

        # The library is the application's; follow it wherever it changed.
        library_events.changed.connect(self.stack.refresh_coatings)
        library_events.changed.connect(self.materials_tab.refresh_library)
        recent.restore_layout(self, _LAYOUT_KEY)

    # -- construction ----------------------------------------------------

    def _build_dock(self) -> None:
        """The stack and the wavelengths, with the send button under them."""
        self.stack = GeometryPanel(sections=("surface", "scan"), base_label="Substrate")
        # Enter in a field has nothing to solve here; the tabs follow `changed`.

        self.send_button = QPushButton("Send to efficiency window")
        self.send_button.setToolTip(
            "Copy this stack into the efficiency window's Surface group:\n"
            "the layers, the substrate as its Grating material, and the\n"
            "roughness. Its grating, mount and wavelengths are left alone,\n"
            "and nothing is solved until you ask."
        )
        self.send_button.clicked.connect(self.send_to_efficiency)

        column = QWidget()
        layout = QVBoxLayout(column)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.stack, 0)
        layout.addWidget(self.send_button)
        layout.addStretch(1)

        scroller = QScrollArea()
        scroller.setWidget(column)
        scroller.setWidgetResizable(True)
        scroller.setMinimumWidth(max(280, fit_width_to_contents(scroller)))

        dock = self.stack_dock = QDockWidget("Coating stack", self)
        dock.setObjectName("stack_dock")
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
        """File, View and Help, each menu held on the window, as the other
        windows hold theirs: a QMenu reached back through the menu bar is
        not reliably usable under PySide6."""
        from .help_menu import install_help_menu

        self._menu_bar = self.menuBar()

        self._file_menu = self._menu_bar.addMenu("&File")
        send = self._file_menu.addAction("&Send to efficiency window")
        send.setShortcut(QKeySequence("Ctrl+E"))
        send.triggered.connect(self.send_to_efficiency)
        script = self._file_menu.addAction("Copy mirror as &Python")
        script.setShortcut(QKeySequence("Ctrl+Shift+C"))
        script.setToolTip("Copy a script that builds this stack and reflects it as the Mirror tab does")
        script.triggered.connect(self.copy_as_python)
        self._file_actions = {"send": send, "python": script}

        self._view_menu = self._menu_bar.addMenu("&View")
        toggle = self.stack_dock.toggleViewAction()
        toggle.setText("&Coating stack")
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

    def open_stack(self, form, marks: str = "") -> None:
        """Take a coating from the efficiency window: its Surface group and
        wavelengths out of ``form``, and its mount's angles as ``marks``."""
        self.stack.write_form(form)
        self.mirror_tab.set_marks(marks)
        if not self.stack_dock.isVisible():
            self.stack_dock.show()

    def python_script(self) -> "str | None":
        """The Mirror tab's calculation as a script, computed afresh from the
        dock so it is what is on screen; ``None`` when it does not compute,
        which the tab's own panel then explains."""
        from ..script import mirror_script

        self.mirror_tab.calculate()
        result = self.mirror_tab._result
        return None if result is None else mirror_script(result)

    def copy_as_python(self) -> None:
        """Put `python_script` on the clipboard, or show why there is none."""
        script = self.python_script()
        if script is None:
            self.tabs.setCurrentWidget(self.mirror_tab)
            self.statusBar().showMessage("Nothing copied: the Mirror tab says what to fix", 6000)
            return
        QApplication.clipboard().setText(script)
        self.statusBar().showMessage(
            f"Copied {len(script.splitlines())} lines of Python to the clipboard", 6000
        )

    def send_to_efficiency(self):
        """Copy this stack into the efficiency window, and return that window.

        Only the Surface group goes: a mirror has no grating or mount to
        send, and the efficiency window's wavelengths are a scan of its own.
        The window is the one the Go menu opens, so there is never a second.
        """
        from ..launcher import mode
        from .launcher import open_window

        form = self.stack.read_form()
        window = open_window(mode("efficiency"))
        window.open_surface(form.coating, form.layers, form.roughness)
        return window

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        recent.save_layout(self, _LAYOUT_KEY)
        super().closeEvent(event)
