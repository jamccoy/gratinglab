"""
The Wiki tab: a page list beside a rendered page.

The view is thin. `gratinglab.metrology.wiki` owns the content and the
ordering, and this module only displays what it is given, so the pages stay
testable without a window and a command-line ``--explain`` could print the
same text.

A page is shown in endstation's reading pane (`endstation.qt.reading`), the
same pane and renderer as the Help menu's theory pages, so the two read alike:
a reading font, a column of readable width, and tables with room around them.
"""
from __future__ import annotations

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from endstation.qt.reading import ReadingBrowser, render_into, text_color
from PySide6.QtWidgets import QHBoxLayout, QListWidget, QListWidgetItem, QWidget

from gratinglab.gui.qt import theming

from ...wiki import pages

__all__ = ["WikiView"]

#: Shown when the pages are missing, so an installed copy without its package
#: data says so rather than presenting an empty panel.
NO_PAGES = (
    "# No wiki pages found\n\n"
    "The markdown files that ship with `gratinglab.metrology.wiki` could not "
    "be located. In an installed copy this means the package data "
    "was not included."
)


class WikiView(QWidget):
    """How the analysis works, in prose."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self._pages = pages()

        self._contents = QListWidget()
        self._contents.setFixedWidth(210)
        for wiki_page in self._pages:
            item = QListWidgetItem(wiki_page.title)
            item.setData(0x0100, wiki_page.slug)  # Qt.UserRole
            self._contents.addItem(item)
        self._contents.currentRowChanged.connect(self._show_row)

        self._browser = ReadingBrowser()
        # These pages cross-reference each other by title in prose rather than by
        # link, so there is nothing for a click to resolve. Opening links would
        # only ever blank the view.
        self._browser.setOpenExternalLinks(False)
        self._browser.setOpenLinks(False)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        layout.addWidget(self._contents)
        layout.addWidget(self._browser, stretch=1)

        if self._pages:
            self._contents.setCurrentRow(0)
        else:
            self._render(NO_PAGES)

        # Any equation is an image in the text colour, so a theme change
        # re-renders the page rather than only repainting it.
        theming.subscribe(self._on_theme)

    # ── State, for tests and callers ─────────────────────────────────────────

    @property
    def slugs(self) -> tuple[str, ...]:
        """Slugs currently listed, in display order."""
        return tuple(p.slug for p in self._pages)

    @property
    def current_slug(self) -> str | None:
        row = self._contents.currentRow()
        if 0 <= row < len(self._pages):
            return self._pages[row].slug
        return None

    def rendered_html(self) -> str:
        """The page as currently displayed."""
        return self._browser.toHtml()

    def show_page(self, slug: str) -> bool:
        """Select a page by slug. False if there is no such page."""
        for row, wiki_page in enumerate(self._pages):
            if wiki_page.slug == slug:
                self._contents.setCurrentRow(row)
                return True
        return False

    # ── Internals ────────────────────────────────────────────────────────────

    def _render(self, markdown: str) -> None:
        render_into(
            self._browser,
            markdown,
            device_pixel_ratio=self.devicePixelRatioF(),
            text_color=text_color(self),
        )

    def _show_row(self, row: int) -> None:
        if 0 <= row < len(self._pages):
            self._render(self._pages[row].markdown)
            self._browser.verticalScrollBar().setValue(0)

    def _on_theme(self, _theme) -> None:
        # A bound method, so Qt disconnects it when the view is destroyed.
        # The reader keeps their place: only the colours changed.
        row = self._contents.currentRow()
        if 0 <= row < len(self._pages):
            scroll = self._browser.verticalScrollBar()
            position = scroll.value()
            self._render(self._pages[row].markdown)
            scroll.setValue(position)
