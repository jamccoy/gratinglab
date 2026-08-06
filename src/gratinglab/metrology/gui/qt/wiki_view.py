"""
The Wiki tab: a page list beside a rendered page.

Thin. `afm_analysis.wiki` owns the content and the ordering; this does nothing
but display what it is given, so the pages stay testable without a window and a
command-line ``--explain`` could print the same text.

Markdown is rendered by Qt itself (`QTextBrowser.setMarkdown`), which handles
headings, tables, lists, code spans and emphasis - everything these pages use.
GratingLab carries a LaTeX rasteriser for its theory pages because diffraction
algebra needs typesetting; the formulas here are one-liners and read fine as
text, so that machinery is not worth the ~670 lines here.
"""
from __future__ import annotations

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

from PySide6.QtWidgets import (
    QHBoxLayout, QListWidget, QListWidgetItem, QTextBrowser, QWidget,
)

from ...wiki import pages

__all__ = ["WikiView"]


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

        self._browser = QTextBrowser()
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
            # A checkout missing its package data should say so rather than
            # presenting an empty panel that looks like a rendering bug.
            self._browser.setMarkdown(
                "# No wiki pages found\n\n"
                "The markdown files that ship with `afm_analysis.wiki` could not "
                "be located. In an installed copy this means the package data "
                "was not included."
            )

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

    def _show_row(self, row: int) -> None:
        if 0 <= row < len(self._pages):
            self._browser.setMarkdown(self._pages[row].markdown)
            self._browser.verticalScrollBar().setValue(0)
