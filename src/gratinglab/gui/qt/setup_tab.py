"""The Setup tab: what optical constants are installed, and where from.

Until the materials layer existed, this tab was a stub that stated as much
instead of offering a coating field that did nothing, because an error or a
no-op must never be presented as a real choice (`docs/theory/scalar.md` §5,
on why the obliquity factor was removed instead of being kept as an "opt-in
variant").

The layer now exists, but the control it was waiting for was not placed
here. `coating` and `roughness` are :class:`~gratinglab.problem.Problem`
fields, so they belong in the geometry dock beside the period: one `Problem`
feeds every solver tab, and a second place to set it would be a second source
of truth for the same value.

This tab instead owns the library: which materials are available, what they
were derived from, and how to add one. That content suits a tab, since it is
reference material and not an input, it duplicates no state, and it is what a
user needs when the coating dropdown offers one entry and they wanted nickel.

Its three buttons act on the library and never on the `Problem`: adding a
table, re-reading the directories, and opening the one you own. Each change
emits `library_changed`, which the window routes to the geometry dock's
Coating dropdown, so a table added here is a coating there without a restart.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMessageBox,
    QPushButton,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from ..docs import TheoryPage
from .theory_viewer import render_into

#: What each reader's ``source`` string starts with, and what to cite for it.
_CITATIONS = {
    "CXRO/Henke": (
        "**CXRO/Henke tables** are the Henke–Gullikson–Davis compilation. Work "
        "using them should cite *X-ray interactions: photoabsorption, "
        "scattering, transmission, and reflection at E = 50–30000 eV, "
        "Z = 1–92*, Atomic Data and Nuclear Data Tables **54** (2), 181–342 "
        "(1993)."
    ),
    "refractiveindex.info": (
        "**refractiveindex.info tables** each come from a measurement of their "
        "own, named in the Source column above and carried into every "
        "result's provenance. Cite that measurement, not the database."
    ),
    "PCGrate": (
        "**PCGrate `.ari` tables** carry no citation of their own: the "
        "provenance names the file and nothing more, so record where it came "
        "from beside it."
    ),
}

__all__ = ["SetupTab", "library_page"]


def library_page() -> str:
    """The tab's markdown, built from the installed tables.

    It is generated instead of written out by hand, so it cannot drift from
    what `materials.available()` returns, whereas a hand-maintained list is
    correct only on the day it is written.
    """
    from ...materials import LIBRARY_ENV, available, library_dirs, lookup

    try:
        names = available()
        unlisted = ""
    except ValueError as refused:  # two tables for one name in one directory
        names, unlisted = (), f"**A library could not be read:** {refused}\n\n"
    kinds: set[str] = set()
    if not names:
        listing = unlisted + (
            "**No tables are installed.** That is what an install starts with, "
            "not a broken one: the package ships none, and any CXRO/Henke "
            "export you download becomes a coating when you add it to a "
            "library, below.\n"
        )
    else:
        rows = ["| Material | Range | Points | Source |", "|---|---|---|---|"]
        for name in names:
            # One row per table, and an unreadable file is a row saying so
            # rather than a page that fails: a library is the user's own
            # directory and can hold anything, including a .yml file on an
            # install without the refractiveindex extra (MissingExtra).
            try:
                table = lookup(name)
            except (OSError, ValueError, ImportError) as unreadable:
                rows.append(f"| `{name}` | unreadable | | {unreadable} |")
                continue
            kinds.update(k for k in _CITATIONS if table.source.startswith(k))
            low, high = table.range_nm
            rows.append(
                f"| `{name}` | {low:.3g} – {high:.3g} nm | "
                f"{len(table.wavelength_nm)} | {table.source} |"
            )
        listing = unlisted + "\n".join(rows) + "\n"

    searched = "\n".join(f"    {d}" for d in library_dirs())
    # Cite what is installed; with nothing installed, every kind a library
    # could hold, since that is what the next table will be.
    citations = "\n\n".join(
        text for kind, text in _CITATIONS.items() if kind in kinds or not kinds
    )

    return f"""\
# Materials

Optical constants for the **Coating** dropdown in the geometry panel. On the
Scalar tab, naming one applies its Fresnel reflectivity and makes the result
*absolute*; leaving it empty gives *relative* efficiency, which is a correct
answer to a different question rather than a missing input. The provenance
panel on each solver tab says which you got.

{listing}
Asking for a wavelength outside a table's range **raises** rather than
extrapolating. An extrapolated optical constant is a plausible number with
nothing behind it, and every efficiency computed from one would inherit that.

## What the coating means to the integral solver

The **Conductivity** selector on the Integral tab decides whether the coating
is read at all, because there the material is the boundary condition rather
than a factor applied afterwards.

- **perfect** (default) is a perfectly conducting boundary. It has no material
  in it, so a named coating is *ignored*, and the provenance panel says so rather
  than letting the setting look consulted. Efficiencies are relative to a
  perfect reflector and sum to 1, which is a theorem, so the summed deviation
  is a direct read on discretisation error.
- **tabulated** reads the coating's complex index at every wavelength and
  solves the coupled system of Goray & Schmidt (2010). Efficiencies are
  absolute, the absorbed fraction is computed independently, and `R + A = 1`
  is the theorem that replaces the sum rule. It needs a coating: without one
  there is nothing rigorous to compute, and the form says so instead of
  quietly defaulting.

Folding a Fresnel reflectivity onto a perfect-conductivity result is **not**
equivalent to either of these, which is why the integral tab offers no
reflectivity model.

Tabulated costs roughly twice as much per wavelength, and demands a finer
mesh: the field oscillates and decays inside the material on the scale
`λ/|n|`, and the **Boundary pts** guard refuses a mesh that cannot resolve it,
naming the number it needs.

## Adding a material

**Add material…** above picks a table, asks the name it should have as a
coating, reads it to prove it is usable, and copies it into your library; the
Coating dropdown has it at once. **Rescan library** re-reads the directories
after you add or remove a file yourself, and **Open library folder** opens the
first of them, creating it if it does not exist yet.

A library is a directory of tables, one material per file and named for it.
Any of these, saved as `Pt` plus its suffix, is the coating `Pt`:

- `Pt.txt`: an `(Energy, Delta, Beta)` export from
  [CXRO](https://henke.lbl.gov/optical_constants/), for X-rays.
- `Pt.yml`: a file from [refractiveindex.info](https://refractiveindex.info),
  ultraviolet to infrared. The database names files for their authors
  (`Au/nk/Johnson.yml`), so rename it to the material. Needs the
  `refractiveindex` extra.
- `Pt.ari`: a PCGrate table.

These directories are searched in this order, so a table of your own replaces
a bundled one of the same name:

{searched}

The first is yours. Set `{LIBRARY_ENV}` to use other directories instead,
separated as your platform separates `PATH`. The list above reads them, so no
code changes. Or load a table from anywhere without adding it:

    from gratinglab import materials
    nickel = materials.from_cxro_file("Ni_CXRO.txt")

## Provenance

{citations}

The package ships no tables: CXRO's terms for redistributing them inside a
BSD-3 package have not been confirmed, and `SOURCES.md` beside the bundled
directory records why. Each table's source names the file it was read from,
so a result's provenance says which library answered.

## Surface roughness

The **Roughness σ** field damps the reflectivity through the Névot–Croce
factor. It does nothing without a coating, because there is no reflectivity to
damp, and the solver does not pretend otherwise. Debye–Waller is offered
alongside it in the **Roughness** selector on the Scalar tab; the two differ by
about 1% near the critical angle, reverse order well below it, and converge far
above it.

## Reflectivity across the groove

The **Reflectivity** selector on the Scalar tab chooses how much of the groove
the reflection calculation looks at, and also does nothing without a coating.

- **local** (default) evaluates the Fresnel amplitude at every quadrature point
  from the local facet tilt and carries it inside the diffraction integral. A
  groove whose reflectivity varies across the cycle is an amplitude grating as
  well as a phase grating, so reflectivity becomes order-dependent.
- **average** takes the groove-cycle mean of the intensity: one factor per
  wavelength, but one that sees the shadowing.
- **facet** is the older treatment: a single reflectivity at the active-facet
  angle. Kept so an earlier run can be reproduced; it breaks reciprocity, so it
  is not the default.

A facet turned away from the beam contributes nothing under the two resolved
models. On a blazed groove that can be a sixth of the period, and it can drive
a weak order to exactly zero, which the provenance panel names, so it is never
mistaken for an order passing off.
"""


class SetupTab(QWidget):
    """Reference material, presented as a page, with three library actions.

    It has no ``solve_requested`` or ``build_options``, because it is not
    part of the solve/cancel contract that every solver tab implements, and
    ``MainWindow`` never routes a solve to it. ``library_changed`` fires
    after any action that may have changed what is installed.
    """

    library_changed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(False)
        self._browser.setOpenLinks(False)
        # Links leave for the system browser; the page itself never navigates.
        self._browser.anchorClicked.connect(self._follow)

        self._add_button = QPushButton("Add material…")
        self._add_button.setToolTip(
            "Pick a CXRO .txt, refractiveindex.info .yml or PCGrate .ari "
            "table, name it, and copy it into your library"
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

        buttons = QHBoxLayout()
        for button in (self._add_button, self._rescan_button, self._folder_button):
            buttons.addWidget(button)
        buttons.addStretch(1)
        self._status = QLabel("")
        self._status.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.addLayout(buttons)
        layout.addWidget(self._status)
        layout.addWidget(self._browser)
        self._render()
        from . import theming

        theming.subscribe(self._on_theme)

    def _on_theme(self, _theme) -> None:
        # A bound method, not a lambda: Qt drops the connection when this
        # widget is destroyed, and would keep calling a lambda into a corpse.
        self._render()

    # -- library actions --------------------------------------------------

    def rescan(self) -> None:
        """Re-read every library directory into the page and the dropdown."""
        from ...materials import available

        self._render()
        try:
            count = len(available())
        except ValueError:
            self._status.setText("A library could not be read; see below.")
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
        name = self._choose_name(Path(chosen).stem)
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
        self._render()
        self._status.setText(
            f"Added {name} from {Path(chosen).name} as {written}. It is in the "
            "Coating dropdown now."
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
            str(Path.home()),
            "Optical constants (*.txt *.yml *.yaml *.ari);;All files (*)",
        )
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

    # -- page -------------------------------------------------------------

    def _follow(self, url: QUrl) -> None:
        if url.scheme() in ("http", "https"):
            QDesktopServices.openUrl(url)

    def _render(self) -> None:
        page = TheoryPage(
            name="setup",
            title="Setup",
            available=True,
            path=None,
            text=library_page(),
            # Not an approximate method -- suppresses the "approximate" banner,
            # which would be nonsense on a page about where data came from.
            rigorous=True,
        )
        render_into(
            self._browser,
            page,
            device_pixel_ratio=self.devicePixelRatioF(),
            text_color=self.palette().windowText().color().name(),
        )
