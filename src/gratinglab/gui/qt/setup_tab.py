"""The Setup tab: an honest stub, not a fake control.

`Problem.coating` exists as a field today and is completely inert -- nothing
reads it. The real materials layer (CXRO optical constants, Fresnel
reflectivity, roughness factors) is deferred milestone "Materials layer" in
`docs/roadmap.md`, not started. Putting a coating text field on this tab
before that layer exists would be exactly the mistake the scalar formulation
work already named and rejected once: an error or a no-op must never be
offered as though it were a real choice (see `docs/theory/scalar.md` §5 on
why the obliquity factor was removed rather than kept as an "opt-in
variant"). A control that looks like it does something and does nothing is
worse than no control at all.

So this tab explains itself instead, reusing the exact rendering path
:class:`~.theory_viewer.TheoryViewer` uses for an unwritten theory page --
same typeset markdown, same "not yet" honesty, no new UI invented for it.
"""

from __future__ import annotations

from PySide6.QtWidgets import QTextBrowser, QVBoxLayout, QWidget

from ..docs import TheoryPage
from .theory_viewer import render_into

__all__ = ["SetupTab"]

_TEXT = """\
# Setup

Nothing configurable lives here yet.

This tab will hold project-level and materials configuration once the
**materials layer** lands (CXRO/Henke optical constants, Fresnel
reflectivity, and the Névot–Croce and Debye–Waller roughness factors --
see `docs/roadmap.md`). `Problem.coating` already exists as a field, but
nothing in the solver reads it today.

No coating control is offered in the meantime. A field that looks like it
sets something, and silently does nothing, is worse than no field at all --
the scalar formulation work reached exactly that conclusion once already,
about an option that turned out to be an error rather than a real choice
(see `docs/theory/scalar.md` section 5). Every result today is **relative
efficiency**, with no coating applied, which the provenance panel on each
solver tab already states.
"""


class SetupTab(QWidget):
    """A stub, not a form. Deliberately has no `solve_requested` or
    `build_options` -- it is not part of the solve/cancel contract every
    solver tab implements, and `MainWindow` never routes to it."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(False)
        self._browser.setOpenLinks(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.addWidget(self._browser)

        page = TheoryPage(
            name="setup", title="Setup", available=False, path=None,
            text=_TEXT, rigorous=True,
        )
        render_into(
            self._browser,
            page,
            device_pixel_ratio=self.devicePixelRatioF(),
            text_color=self.palette().windowText().color().name(),
        )
