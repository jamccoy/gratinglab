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
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..solvers.base import available_solvers, get_solver

__all__ = ["TheoryPage", "find_theory_root", "theory_pages"]

#: Human-readable titles for known solvers. Falls back to the bare registry
#: name for anything not listed here, so a new solver still gets a menu entry
#: before anyone remembers to add a nicer title.
_TITLES = {
    "scalar": "Scalar (Kirchhoff)",
}

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


def theory_pages() -> tuple[TheoryPage, ...]:
    """One entry per registered solver, in registration order."""
    root = find_theory_root()
    pages = []

    for name in available_solvers():
        title = _TITLES.get(name, name)
        rigorous = get_solver(name).capabilities.rigorous

        if root is None:
            pages.append(
                TheoryPage(
                    name=name,
                    title=title,
                    available=False,
                    path=None,
                    text=(
                        "Documentation is not bundled with this installation "
                        "(docs/ ships in a development checkout, not in a "
                        "built package). See the project repository for "
                        f"docs/theory/{name}.md."
                    ),
                    rigorous=rigorous,
                )
            )
            continue

        path = root / f"{name}.md"
        if path.is_file():
            pages.append(
                TheoryPage(
                    name=name, title=title, available=True, path=path,
                    text=path.read_text(), rigorous=rigorous,
                )
            )
        else:
            pages.append(
                TheoryPage(
                    name=name,
                    title=title,
                    available=False,
                    path=path,
                    text=(
                        f"No theory page has been written yet for '{name}'. "
                        "See docs/roadmap.md for what is planned."
                    ),
                    rigorous=rigorous,
                )
            )

    return tuple(pages)
