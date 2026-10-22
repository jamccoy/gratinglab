"""The Setup tab: what optical constants are installed, and where from.

It was a stub until the materials layer existed, saying so plainly rather than
offering a coating field that did nothing -- an error or a no-op must never be
dressed as a real choice (`docs/theory/scalar.md` §5, on why the obliquity
factor was removed rather than kept as an "opt-in variant").

The layer exists now, and the control it was waiting for did **not** land
here. `coating` and `roughness` are :class:`~gratinglab.problem.Problem`
fields, so they belong in the geometry dock beside the period: one `Problem`
feeds every solver tab, and a second place to set it would be a second source
of truth for the same value.

What this tab owns instead is the *library* -- which materials are available,
what they were derived from, and how to add one. That is genuinely tab-shaped
(it is reference material, not an input), it duplicates no state, and it is
the thing a user needs when the coating dropdown offers one entry and
they wanted nickel.
"""

from __future__ import annotations

from PySide6.QtWidgets import QTextBrowser, QVBoxLayout, QWidget

from ..docs import TheoryPage
from .theory_viewer import render_into

__all__ = ["SetupTab", "library_page"]


def library_page() -> str:
    """The tab's markdown, built from the installed tables.

    Generated rather than written out, so it cannot drift from what
    `materials.available()` returns -- the failure mode of a
    hand-maintained list is that it is right on the day it is written.
    """
    from ...materials import LIBRARY_ENV, available, library_dirs, lookup

    try:
        names = available()
        unlisted = ""
    except ValueError as refused:  # two tables for one name in one directory
        names, unlisted = (), f"**A library could not be read:** {refused}\n\n"
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
            low, high = table.range_nm
            rows.append(
                f"| `{name}` | {low:.3g} – {high:.3g} nm | "
                f"{len(table.wavelength_nm)} | {table.source} |"
            )
        listing = unlisted + "\n".join(rows) + "\n"

    searched = "\n".join(f"    {d}" for d in library_dirs())

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

These are the Henke–Gullikson–Davis compilation. Work using them should cite
*X-ray interactions: photoabsorption, scattering, transmission, and reflection
at E = 50–30000 eV, Z = 1–92*, Atomic Data and Nuclear Data Tables **54** (2),
181–342 (1993).

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
    """Reference material, not a form.

    Has no ``solve_requested`` or ``build_options`` -- it is not
    part of the solve/cancel contract every solver tab implements, and
    ``MainWindow`` never routes to it.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(False)
        self._browser.setOpenLinks(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.addWidget(self._browser)

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
