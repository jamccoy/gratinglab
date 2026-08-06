"""
The wiki: how this analysis works, in prose.

Pure - no Qt, no plotting, no file writing. The GUI displays these pages, but it
does not own them: they describe the analysis rather than the interface, so a
command-line ``--explain`` could print exactly the same text.

Pages are markdown files sitting next to this module and shipped as package data,
so they are present in an installed copy and not only in a checkout. Ordering is
explicit rather than alphabetical, because a reader should meet the overview
before the correction it motivates.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

__all__ = ["WikiPage", "pages", "page", "PAGE_ORDER"]

WIKI_DIR = Path(__file__).resolve().parent

#: Reading order. A slug missing from disk is skipped rather than raising - a
#: half-written wiki should not stop the window opening - and any page found on
#: disk but absent here is appended, so a new file shows up without needing to
#: be registered in two places.
PAGE_ORDER = (
    "overview",
    "row-groups",
    "icc-correction",
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


def _title_of(markdown: str, slug: str) -> str:
    """The leading `# ` heading, or a readable fallback built from the slug."""
    match = _TITLE_RE.search(markdown)
    if match:
        return match.group(1).strip()
    return slug.replace("-", " ").capitalize()


def _load(path: Path) -> WikiPage:
    markdown = path.read_text(encoding="utf-8")
    slug = path.stem
    return WikiPage(slug=slug, title=_title_of(markdown, slug), markdown=markdown)


def pages() -> tuple[WikiPage, ...]:
    """
    Every wiki page, in reading order.

    Read from disk at call time rather than cached at import, so editing a page
    and reopening the tab shows the edit - the pages are prose under active
    revision, and a stale cache would be a confusing thing to debug.
    """
    found = {path.stem: path for path in sorted(WIKI_DIR.glob("*.md"))}

    ordered = [found.pop(slug) for slug in PAGE_ORDER if slug in found]
    ordered.extend(found[slug] for slug in sorted(found))

    return tuple(_load(path) for path in ordered)


def page(slug: str) -> WikiPage | None:
    """One page by slug, or None if there is no such page."""
    for candidate in pages():
        if candidate.slug == slug:
            return candidate
    return None
