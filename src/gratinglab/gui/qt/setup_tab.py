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
the thing a user actually needs when the coating dropdown offers one entry and
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
    `materials.available()` actually returns -- the failure mode of a
    hand-maintained list is that it is right on the day it is written.
    """
    from ...materials import available, data_dir, lookup

    names = available()
    if not names:
        listing = (
            "**No tables are installed.** That is a supported configuration, "
            "not a broken one: the reader works on any CXRO/Henke export you "
            "download, and nothing in `gratinglab` requires a bundled file.\n"
        )
    else:
        rows = ["| Material | Range | Points | Source |", "|---|---|---|---|"]
        for name in names:
            table = lookup(name)
            low, high = table.range_nm
            rows.append(
                f"| `{name}` | {low:.3g} – {high:.3g} nm | "
                f"{len(table.wavelength_nm)} | {table.source} |"
            )
        listing = "\n".join(rows) + "\n"

    return f"""\
# Materials

Optical constants for the **Coating** dropdown in the geometry panel. Naming
one applies its Fresnel reflectivity and makes the result *absolute*; leaving
it empty gives *relative* efficiency, which is a correct answer to a different
question rather than a missing input. The provenance panel on each solver tab
says which you got.

{listing}
Asking for a wavelength outside a table's range **raises** rather than
extrapolating. An extrapolated optical constant is a plausible number with
nothing behind it, and every efficiency computed from one would inherit that.

## Adding a material

Download an `(Energy, Delta, Beta)` table from
[CXRO](https://henke.lbl.gov/optical_constants/) and drop it in as
`<Element>.txt` here:

    {data_dir()}

The list above reads that directory, so no code changes. Or load one from
anywhere without installing it:

    from gratinglab import materials
    nickel = materials.from_cxro_file("Ni_CXRO.txt")

## Provenance

These are the Henke–Gullikson–Davis compilation. Work using them should cite
*X-ray interactions: photoabsorption, scattering, transmission, and reflection
at E = 50–30000 eV, Z = 1–92*, Atomic Data and Nuclear Data Tables **54** (2),
181–342 (1993).

`SOURCES.md` in the directory above records per-file provenance, and states
plainly that redistribution terms have not been formally confirmed for
bundling inside a BSD-3 package. If they do not hold up, deleting that
directory is the whole fix -- the reader keeps working on user-supplied files.

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
- **average** takes the groove-cycle mean of the intensity — one factor per
  wavelength, but one that sees the shadowing.
- **facet** is the older treatment: a single reflectivity at the active-facet
  angle. Kept so an earlier run can be reproduced; it breaks reciprocity, so it
  is not the default.

A facet turned away from the beam contributes nothing under the two resolved
models. On a blazed groove that can be a sixth of the period, and it can drive
a weak order to exactly zero — which the provenance panel names, so it is never
mistaken for an order passing off.
"""


class SetupTab(QWidget):
    """Reference material, not a form.

    Deliberately has no ``solve_requested`` or ``build_options`` -- it is not
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
