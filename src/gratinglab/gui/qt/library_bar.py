"""The material library's three actions, and the line that reports on them.

Adding a table, re-reading the directories, and opening the one you own.
They act on the library and never on the `Problem`. Each one that may have
changed what is installed emits ``library_changed``, which the window routes
to everything that lists materials: the dock's dropdowns, the Materials tab's
list, and the Setup page's table.

These lived on the Setup tab until the Materials tab existed; they are here so
that tab can hold them beside the list they change. The dialogs are one
method each so a test can answer them.

Built by :class:`~gratinglab.gui.qt.materials_tab.MaterialsTab`.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QFileDialog,
    QInputDialog,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from . import recent

__all__ = ["LibraryBar", "library_events"]


class _LibraryEvents(QObject):
    """One ``changed`` signal for the whole application.

    The library is the same for every window, and the bar that changes it is
    in the coatings window while the efficiency window's dropdowns list it
    too. A window connects to this rather than to a bar it may not hold.
    """

    changed = Signal()


#: Connect to ``library_events.changed`` to follow the library.
library_events = _LibraryEvents()


class LibraryBar(QWidget):
    """Add material, Rescan library, Open library folder, and a status line."""

    library_changed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._add_button = QPushButton("Add material…")
        self._add_button.setToolTip(
            "Pick a CXRO .txt, refractiveindex.info .yml, PCGrate .ari, plain\n"
            "n,k .nk or .csv table, or one element's .nff scattering factors,\n"
            "and copy it into your library"
        )
        self._add_button.clicked.connect(self.add_material)
        self._rescan_button = QPushButton("Rescan library")
        self._rescan_button.setToolTip(
            "Re-read the library directories, after adding or removing a file "
            "outside this window"
        )
        self._rescan_button.clicked.connect(self.rescan)
        self._folder_button = QPushButton("Open library folder")
        self._folder_button.clicked.connect(self.open_library_folder)
        self.library_changed.connect(library_events.changed)
        self._status = QLabel("")
        self._status.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        for widget in (self._add_button, self._rescan_button, self._folder_button, self._status):
            layout.addWidget(widget)

    # -- library actions --------------------------------------------------

    def rescan(self) -> None:
        """Re-read every library directory into everything that lists it."""
        from ...materials import available

        try:
            count = len(available())
        except ValueError:
            self._status.setText("A library could not be read; the Setup tab says why.")
        else:
            self._status.setText(
                f"Library read: {count} material{'s' if count != 1 else ''}."
            )
        self.library_changed.emit()

    def add_material(self) -> None:
        """Pick a table, name it, validate it, and install it."""
        from ...materials import install

        chosen = self._choose_file()
        if not chosen:
            return
        # One element's scattering factors are filed under its symbol, which
        # the file name already gives, so there is no coating name to ask for.
        factors = Path(chosen).suffix.lower() == ".nff"
        name = Path(chosen).stem.capitalize() if factors else self._choose_name(Path(chosen).stem)
        if not name:
            return
        try:
            try:
                written = install(chosen, name)
            except FileExistsError:
                if not self._confirm_replace(name):
                    raise
                written = install(chosen, name, overwrite=True)
        except (OSError, ValueError, ImportError) as refused:
            # A table that does not read is reported and never installed:
            # the library is checked here rather than inside a later solve.
            self._status.setText(f"{Path(chosen).name} was not added: {refused}")
            return
        if factors:
            self._status.setText(
                f"Added {name}'s scattering factors as {written}. Use them as "
                f"factors:<formula>@<density>, e.g. factors:{name}@<g/cm^3>."
            )
        else:
            self._status.setText(
                f"Added {name} from {Path(chosen).name} as {written}. It is in the "
                "library list and the dock's material dropdowns now."
            )
        self.library_changed.emit()

    def open_library_folder(self) -> None:
        """Open your library directory, creating it first if need be."""
        from ...materials import library_dirs

        directory = library_dirs()[0]
        try:
            directory.mkdir(parents=True, exist_ok=True)
        except OSError as refused:
            self._status.setText(f"Could not create {directory}: {refused}")
            return
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(directory))):
            self._status.setText(f"Your library is {directory}")

    # -- dialogs, one method each so a test can answer them ----------------

    def _choose_file(self) -> str:
        path, _filter = QFileDialog.getOpenFileName(
            self,
            "Add material",
            recent.start_in("tables", fallback=str(Path.home())),
            "Optical constants or scattering factors "
            "(*.txt *.yml *.yaml *.ari *.nk *.csv *.nff);;All files (*)",
        )
        recent.remember("tables", path)
        return path

    def _choose_name(self, default: str) -> str:
        name, accepted = QInputDialog.getText(
            self,
            "Material name",
            "Name as a coating (refractiveindex.info files are named for "
            "their author, so rename those):",
            text=default,
        )
        return name.strip() if accepted else ""

    def _confirm_replace(self, name: str) -> bool:
        answer = QMessageBox.question(
            self,
            "Replace material",
            f"Your library already has {name}. Replace it?",
        )
        return answer == QMessageBox.StandardButton.Yes
