"""Locating each solver's theory page.

Pure and tkinter-free, like :mod:`gratinglab.gui.state`, because this is the one
piece of "which page goes with which solver" logic worth testing on its own.

``docs/theory/`` lives at the **repo root**, outside ``src/gratinglab/``, so it
ships in a development checkout but not in a built wheel. That is a real
limitation of a pre-alpha, locally-run tool -- not something to paper over by
duplicating markdown into package data before there is a second solver to
justify it. :func:`theory_pages` degrades to an explanatory, unavailable entry
rather than raising when the directory cannot be found.

Every *registered* solver gets an entry, available or not, so the Help menu
always reflects what solvers exist. A solver with no page yet still shows up
and says so, rather than silently vanishing from the menu until someone
remembers to write it.

Alongside the per-solver pages, :func:`general_pages` returns entries that are
not tied to any solver -- currently just the normative ``conventions.md``,
which every solver's page cross-references for shared geometry (the
generalized grating equation, the angle conventions) rather than duplicating
it. Kept as a separate function rather than folded into :func:`theory_pages`
so that function's existing contract -- one entry per registered solver, no
more, no less -- stays exactly true; a test asserts it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..solvers.base import available_solvers, get_solver

__all__ = [
    "TheoryPage",
    "find_theory_root",
    "theory_pages",
    "general_pages",
    "display_title",
]

#: Human-readable titles for known solvers. Falls back to the bare registry
#: name for anything not listed here, so a new solver still gets a menu entry
#: before anyone remembers to add a nicer title.
_TITLES = {
    "scalar": "Scalar (Kirchhoff)",
}

#: Non-solver reference pages: (key, title, filename relative to docs/).
_GENERAL_PAGES = (
    ("conventions", "Grating Geometry & Conventions", "conventions.md"),
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
        not -- always something the viewer can show, never ``None``.
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

    Public so both the Help menu and a solver's own tab label read the same
    name -- there is exactly one place that decides what a solver is called,
    not two that could quietly drift apart.
    """
    return _TITLES.get(name, name)


def find_theory_root(start: Path | None = None) -> Path | None:
    """Search upward from ``start`` for a ``docs/theory`` directory.

    Returns ``None`` rather than raising if none is found within
    :data:`_MAX_SEARCH_DEPTH` levels -- e.g. running from an installed wheel,
    where ``docs/`` was never packaged.
    """
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
            text=path.read_text(), rigorous=rigorous,
        )
    return TheoryPage(
        name=name, title=title, available=False, path=path,
        text=f"No page has been written yet for '{name}' ({not_written}).",
        rigorous=rigorous,
    )


def theory_pages() -> tuple[TheoryPage, ...]:
    """One entry per registered solver, in registration order.

    Exactly one entry per name in ``available_solvers()`` -- a solver with no
    page yet still appears, marked unavailable, rather than being silently
    omitted until someone remembers to write it.
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

    ``rigorous=True`` throughout: these are geometry/convention references,
    not a solver's approximation, so the viewer's approximate-method banner
    never applies to them.
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
