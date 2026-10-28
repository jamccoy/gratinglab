"""
The wiki: how this analysis works, in prose.

The module is pure, with no Qt, no plotting and no file writing. The GUI
displays these pages but does not own them: they describe the analysis and not
the interface, so a command-line ``--explain`` could print exactly the same
text.

Pages are markdown files sitting next to this module and shipped as package
data, so they are present in an installed copy and not only in a checkout. The
ordering is explicit instead of alphabetical, because a reader should meet the
overview before the correction it motivates.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from endstation.docs.index import load_pages

__all__ = ["WikiPage", "pages", "page", "PAGE_ORDER"]

WIKI_DIR = Path(__file__).resolve().parent

#: Reading order. A slug missing from disk is skipped rather than raising - a
#: half-written wiki should not stop the window opening - and any page found on
#: disk but absent here is appended, so a new file shows up without needing to
#: be registered in two places.
PAGE_ORDER = (
    "overview",
    "file-formats",
    "flattening",
    "groove-orientation",
    "scars-and-masking",
    "tip-correction",
    "row-groups",
    "icc-correction",
    "swath-statistics",
    "facet-fitting",
    "reading-outputs",
)

_TITLE_RE = re.compile(r"^#\s+(.+)$", re.MULTILINE)


@dataclass(frozen=True)
class WikiPage:
    """One page: its slug, its display title, and its markdown."""

    slug: str
    title: str
    markdown: str

    @property
    def is_empty(self) -> bool:
        """True when the page has a heading but no prose under it."""
        body = _TITLE_RE.sub("", self.markdown, count=1)
        return not body.strip()


def pages() -> tuple[WikiPage, ...]:
    """
    Every wiki page, in reading order.

    The pages are read from disk at call time instead of being cached at
    import, so editing a page and reopening the tab shows the edit. The pages
    are prose under active revision, and a stale cache would be confusing to
    debug.

    Which pages exist, their order and their titles come from
    :func:`endstation.docs.index.load_pages`, which applies the rules
    `PAGE_ORDER` describes.
    """
    return tuple(
        WikiPage(slug=entry.slug, title=entry.title, markdown=entry.read())
        for entry in load_pages(WIKI_DIR, PAGE_ORDER)
    )


def page(slug: str) -> WikiPage | None:
    """One page by slug, or None if there is no such page."""
    for candidate in pages():
        if candidate.slug == slug:
            return candidate
    return None
