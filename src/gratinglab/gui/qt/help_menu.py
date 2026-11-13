"""The Help menu, built in one place and installed in every window.

The menu was written into the efficiency window, which is where the theory
pages came from, so the AFM window -- merged in later, from a project that
had no such pages -- offered no route to them at all. That is the wrong way
round for two of them: ``conventions.md`` fixes the angle conventions a
measured profile is reported in, and ``theory/metrology.md`` is about the AFM
half itself. A user measuring a scan was the one user who could not read
either.

So the menu lives here and each window installs it, the same arrangement as
the Go menu in :mod:`gratinglab.gui.qt.launcher`: neither window knows the
other has it, and a page added to :mod:`gratinglab.gui.docs` appears in both
at once.

``menu_bar`` exists because the chooser is a plain ``QWidget`` and has no
``menuBar()`` of its own; it hands in the one it put in its layout. The
windows pass nothing and get theirs.

The menu itself is built by :func:`endstation.qt.help_menu.install_help_menu`;
what this module supplies is the About text and the pages.
"""

from __future__ import annotations

from endstation.qt.help_menu import HELP_MENU_TITLE
from endstation.qt.help_menu import install_help_menu as _install
from PySide6.QtWidgets import QMessageBox

__all__ = ["HELP_MENU_TITLE", "install_help_menu", "show_about", "show_theory"]

def _about_text() -> str:
    from ... import __version__, __version_date__
    from .. import provenance

    return provenance.about_text(__version__, __version_date__)


def show_about(parent=None) -> None:
    """The About box. Its text is :func:`gratinglab.gui.provenance.about_text`,
    which is a plain string so it can be asserted without a dialog."""
    QMessageBox.about(parent, "About gratinglab", _about_text())


def show_theory(page, parent=None):
    """Open the read-only viewer for one theory or reference page.

    ``parent`` is what keeps the viewer alive: a top-level Qt window with no
    reference on either side is collected and vanishes a moment after it
    appears.
    """
    from .theory_viewer import TheoryViewer

    viewer = TheoryViewer(page, parent=parent)
    viewer.show()
    return viewer


def install_help_menu(window, menu_bar=None):
    """Add the Help menu to ``window`` and return it.

    The caller holds the returned QMenu. Under PySide6 6.11 a QMenu reached
    back through ``menuBar().actions()[i].menu()`` is not reliably usable --
    it raises "Internal C++ object already deleted" while the menu the window
    built is alive and correct -- so whatever builds a menu must keep it. It
    is also stored on the window, so a caller that forgets still cannot lose
    it.
    """
    from ..docs import general_pages, theory_pages
    from ..richtext import menu_label

    def entry(page, suffix):
        return (menu_label(page.title, suffix), lambda: show_theory(page, window))

    def missing(page):
        return "" if page.available else " (not written yet)"

    # Foundational, solver-independent reference first -- the generalized
    # grating equation and the angle conventions -- then one page per method.
    menu = _install(
        window,
        about_title="About gratinglab",
        about_text=_about_text,
        pages=[
            [entry(page, missing(page)) for page in general_pages()],
            [entry(page, " Theory" + missing(page)) for page in theory_pages()],
        ],
        menu_bar=menu_bar,
    )
    window._gratinglab_help_menu = menu
    return menu
