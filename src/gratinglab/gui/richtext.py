"""What gratinglab adds to endstation's rendering of a Help page.

The markdown and the equations are rendered by :mod:`endstation.docs.markdown`,
which began as this module. Two things stay here because they are about
gratinglab's pages rather than markdown: the banner above the page of a
solver that approximates, and the escaping a page title needs to become a
menu entry.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from endstation.docs.markdown import RichPage, render_markdown

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from .docs import TheoryPage

__all__ = ["BANNER_COLOR", "banner", "menu_label", "to_html"]

#: Colour of the approximate-method banner. The same as endstation's colour
#: for an equation that failed to render, and kept as its own constant: the
#: two facts are unrelated, and a test that cannot tell them apart could one
#: day pass for the wrong reason.
BANNER_COLOR = "#a5370d"


def banner(page: "TheoryPage") -> str:
    """The HTML placed above the page, empty for a rigorous solver.

    A solver whose ``Capabilities.rigorous`` is False says so at the top, so a
    reader learns about the approximation before section 7 rather than in it.
    """
    if page.rigorous:
        return ""
    return (
        f'<p style="color:{BANNER_COLOR}"><b>'
        "⚠ Approximate method: see Limits below."
        "</b></p>"
    )


def to_html(page: "TheoryPage", **kwargs) -> RichPage:
    """A page rendered by :func:`endstation.docs.markdown.render_markdown`,
    banner included. ``kwargs`` go to ``render_markdown``."""
    rendered = render_markdown(page.text, **kwargs)
    notice = banner(page)
    if not notice:
        return rendered
    return RichPage(
        html=notice + "\n" + rendered.html,
        images=rendered.images,
        unrendered=rendered.unrendered,
    )


def menu_label(title: str, suffix: str = "") -> str:
    """A menu entry's text, with Qt's mnemonic character neutralised.

    Qt reads ``&`` in an action's text as "underline the next character", so
    the page titled "Grating Geometry & Conventions" would appear as "Grating
    Geometry Conventions" with a stray underline. Doubling it escapes it.
    """
    return (title + suffix).replace("&", "&&")
