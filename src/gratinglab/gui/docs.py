"""Locate each solver's theory page.

This module is pure and independent of any toolkit, like
:mod:`gratinglab.gui.state`, because the logic that matches pages to solvers
is worth testing on its own.

``docs/theory/`` lives at the repository root, outside ``src/gratinglab/``,
because so much of the tree names the pages by that path. ``setup.py`` copies
the pages the menu shows into ``gratinglab/_docs/`` when a wheel is built, in
the same layout, so an install finds them there and a checkout finds
``docs/`` by walking up. When neither is found, :func:`theory_pages` returns
an unavailable entry with an explanation instead of raising.

A link on a page is resolved by :func:`resolve_link`: to another page the
menu shows, opened in the viewer, or else to the file on GitHub at the
installed version's tag, so the page and the code it names are one release.

Every registered solver receives an entry, whether or not its page is
available, so the Help menu always reflects which solvers exist. A solver
without a page still appears in the menu and states that its page is missing.

Alongside the per-solver pages, :func:`general_pages` returns entries that are
not tied to any solver. At present this includes the normative
``conventions.md``, which every solver's page cross-references for the shared
geometry (the generalized grating equation and the angle conventions) instead
of duplicating it. These entries are returned by a separate function so that
the contract of :func:`theory_pages`, exactly one entry per registered solver,
continues to hold; a test asserts it.
"""

from __future__ import annotations

import posixpath
from dataclasses import dataclass
from pathlib import Path

from .. import __version__
from ..solvers.base import available_solvers, get_solver

__all__ = [
    "TheoryPage",
    "find_theory_root",
    "theory_pages",
    "general_pages",
    "display_title",
    "resolve_link",
    "REPOSITORY",
]

#: Where the Help pages are in an installed package; setup.py puts them there.
_BUNDLED_THEORY = Path(__file__).resolve().parents[1] / "_docs" / "theory"

#: The repository a link to anything but a Help page opens in.
REPOSITORY = "https://github.com/jamccoy/gratinglab"

#: Human-readable titles for known solvers. Falls back to the bare registry
#: name for anything not listed here, so a new solver still gets a menu entry
#: before anyone remembers to add a nicer title.
_TITLES = {
    "cmethod": "C-method (Chandezon)",
    "scalar": "Scalar (Kirchhoff)",
    "integral": "Integral (boundary integral equation)",
    "rcwa": "RCWA (Fourier modal method)",
}

#: Non-solver reference pages: (key, title, filename relative to docs/).
_GENERAL_PAGES = (
    ("conventions", "Grating Geometry & Conventions", "conventions.md"),
    # Not a solver page -- nothing here approximates a wave equation -- but an
    # efficiency computed from a measured profile inherits every assumption
    # made while producing that profile, and those are invisible in the curve.
    # It belongs beside the conventions for the same reason they do.
    ("metrology", "Groove Metrology", "theory/metrology.md"),
)

#: How many parent directories to search before giving up. gui/docs.py sits at
#: src/gratinglab/gui/, four levels below a repo root that holds docs/theory/,
#: so this comfortably covers the dev-checkout case with room to spare.
_MAX_SEARCH_DEPTH = 6


@dataclass(frozen=True, slots=True)
class TheoryPage:
    """One Help-menu entry.

    Attributes
    ----------
    name
        Solver registry key, e.g. ``"scalar"``.
    title
        Display title for the menu.
    available
        Whether an actual page was found.
    path
        The theory file, if found.
    text
        The page content if available, otherwise a short explanation of why
        it is not. It is always text the viewer can show, and never ``None``.
    rigorous
        The solver's own ``Capabilities.rigorous``. Used to decide whether the
        viewer prefixes an approximate-method banner.
    """

    name: str
    title: str
    available: bool
    path: Path | None
    text: str
    rigorous: bool


def display_title(name: str) -> str:
    """A solver's human-readable name, e.g. ``"scalar"`` -> ``"Scalar
    (Kirchhoff)"``.

    It is public so that the Help menu and a solver's own tab label read the
    same name. A solver's display name is therefore decided in one place, and
    the two cannot drift apart.
    """
    return _TITLES.get(name, name)


def short_title(name: str) -> str:
    """A solver's tab label: its display title without the bracketed gloss.

    ``"Integral (boundary integral equation)"`` becomes ``"Integral"``. The
    full title was the tab label until 2026-10-27, and with four tabs plus the
    geometry dock it already overflowed the bar, so every new solver made it
    worse. The full title stays in the Help menu and the tab's tooltip.
    """
    return display_title(name).split(" (", 1)[0]


def find_theory_root(start: Path | None = None) -> Path | None:
    """The directory holding the theory pages.

    With no ``start``, the copy bundled into an installed package, if there
    is one. Otherwise, or failing that, search upward from ``start`` (or
    this file) for a ``docs/theory`` directory, which is where a checkout
    keeps them. Returns ``None`` instead of raising if neither is found
    within :data:`_MAX_SEARCH_DEPTH` levels.
    """
    if start is None and _BUNDLED_THEORY.is_dir():
        return _BUNDLED_THEORY
    here = (start or Path(__file__)).resolve()
    for candidate in [here, *here.parents][: _MAX_SEARCH_DEPTH + 1]:
        theory = candidate / "docs" / "theory"
        if theory.is_dir():
            return theory
    return None


def _load_page(
    name: str, title: str, path: Path | None, *, rigorous: bool, not_written: str
) -> TheoryPage:
    """Build one entry, handling the three cases every page can be in.

    Shared by :func:`theory_pages` and :func:`general_pages` so the "docs/ not
    bundled" and "file missing" fallback text stays consistent between them.
    """
    if path is None:
        return TheoryPage(
            name=name, title=title, available=False, path=None,
            text=(
                "Documentation is not bundled with this installation "
                "(docs/ ships in a development checkout, not in a built "
                f"package). See the project repository for {not_written}."
            ),
            rigorous=rigorous,
        )
    if path.is_file():
        return TheoryPage(
            name=name, title=title, available=True, path=path,
            text=path.read_text(encoding="utf-8"), rigorous=rigorous,
        )
    return TheoryPage(
        name=name, title=title, available=False, path=path,
        text=f"No page has been written yet for '{name}' ({not_written}).",
        rigorous=rigorous,
    )


def theory_pages() -> tuple[TheoryPage, ...]:
    """One entry per registered solver, in registration order.

    There is exactly one entry per name in ``available_solvers()``. A solver
    with no page yet still appears, marked unavailable, so that it is not
    silently omitted from the menu.
    """
    root = find_theory_root()
    return tuple(
        _load_page(
            name,
            display_title(name),
            None if root is None else root / f"{name}.md",
            rigorous=get_solver(name).capabilities.rigorous,
            not_written=f"docs/theory/{name}.md",
        )
        for name in available_solvers()
    )


def general_pages() -> tuple[TheoryPage, ...]:
    """Reference pages not tied to any one solver.

    ``rigorous=True`` is set for all of them. They are references on geometry
    and conventions and describe no solver's approximation, so the viewer's
    approximate-method banner never applies to them.
    """
    root = find_theory_root()
    return tuple(
        _load_page(
            name,
            title,
            None if root is None else root.parent / filename,
            rigorous=True,
            not_written=f"docs/{filename}",
        )
        for name, title, filename in _GENERAL_PAGES
    )


def resolve_link(
    page: TheoryPage, href: str, pages: tuple[TheoryPage, ...] | None = None
) -> TheoryPage | str | None:
    """What following ``href`` on ``page`` should do.

    Returns another of ``pages`` (by default every page the menu shows) to
    open in the viewer, or a URL for the browser: a web link as written, and
    any other file, a source file, a test or a doc the menu does not show, as
    its page on GitHub at this version's tag. ``None`` when there is nothing
    to follow: an anchor on the same page, which the viewer does not mark, or
    a link from a page with no file behind it.
    """
    if "://" in href or href.startswith("mailto:"):
        return href
    target, _, fragment = href.partition("#")
    if not target or page.path is None:
        return None

    resolved = (page.path.parent / target).resolve()
    for candidate in pages if pages is not None else (*general_pages(), *theory_pages()):
        if candidate.path is not None and candidate.path.resolve() == resolved:
            return candidate

    # The pages sit in docs/ and docs/theory/, in a checkout and bundled alike,
    # so the link's place in the repository follows from the page's.
    folder = "docs/theory" if page.path.parent.name == "theory" else "docs"
    in_repository = posixpath.normpath(f"{folder}/{target}")
    if in_repository.startswith(".."):
        return None
    anchor = f"#{fragment}" if fragment else ""
    return f"{REPOSITORY}/blob/v{__version__}/{in_repository}{anchor}"
