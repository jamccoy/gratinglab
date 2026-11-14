"""The read-only window for one theory or reference page.

The pane and the rendering are endstation's (:mod:`endstation.qt.reading`):
a ``QTextBrowser`` laid out for reading, with the page's equations sized to
its font. This module adds what is gratinglab's: the approximate-method banner
(:func:`gratinglab.gui.richtext.banner`), links between Help pages opening in
windows of their own, and a re-render when the theme changes.
"""

from __future__ import annotations

from endstation.qt.reading import ReadingBrowser, text_color
from endstation.qt.reading import render_into as _render_into
from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QTextBrowser, QVBoxLayout, QWidget

from .. import richtext
from ..docs import resolve_link

__all__ = ["ReadingBrowser", "TheoryViewer", "render_into"]


def render_into(
    browser: QTextBrowser, page, *, device_pixel_ratio: float, text_color: str
) -> None:
    """Render a page into a browser, banner and images included.

    Shared by :class:`TheoryViewer` and :class:`~.setup_tab.SetupTab`.
    """
    _render_into(
        browser,
        page.text,
        device_pixel_ratio=device_pixel_ratio,
        text_color=text_color,
        before=richtext.banner(page),
    )


class TheoryViewer(QWidget):
    """One page, typeset."""

    def __init__(self, page, parent: QWidget | None = None) -> None:
        # A window in its own right rather than a modal dialog: reading the
        # theory while adjusting the grating is the point.
        super().__init__(parent)
        self.setWindowFlag(Qt.WindowType.Window)
        self.setWindowTitle(page.title)
        self.resize(860, 720)

        self._browser = ReadingBrowser()
        self._browser.setOpenExternalLinks(False)
        # QTextBrowser would load a relative link as a file in this window and
        # show raw markdown, so links are followed here instead.
        self._browser.setOpenLinks(False)
        self._browser.anchorClicked.connect(self._follow)
        #: The last page a link opened, held so its window is not collected.
        self.opened: TheoryViewer | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._browser)

        self.set_page(page)
        from . import theming

        # Equations are images rasterised in the text colour, so a theme
        # change re-renders the page rather than only repainting it.
        theming.subscribe(self._on_theme)

    def _follow(self, url: QUrl) -> None:
        """Open a linked Help page in a window of its own, beside this one,
        the way the Help menu opens them; anything else goes to the browser."""
        target = resolve_link(self._page, url.toString())
        if isinstance(target, str):
            QDesktopServices.openUrl(QUrl(target))
        elif target is not None:
            self.opened = TheoryViewer(target, parent=self.parentWidget())
            self.opened.show()

    def _on_theme(self, _theme) -> None:
        # A bound method so Qt disconnects it with the window (see SetupTab).
        self.set_page(self._page)

    def set_page(self, page) -> None:
        """Render a page into the browser, images and all."""
        self._page = page
        render_into(
            self._browser,
            page,
            device_pixel_ratio=self.devicePixelRatioF(),
            text_color=text_color(self),
        )
