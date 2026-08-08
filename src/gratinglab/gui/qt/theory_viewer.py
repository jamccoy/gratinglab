"""The read-only viewer for one theory or reference page.

Thin by design: :mod:`gratinglab.gui.richtext` turns the page into HTML and a
list of rasterised equations, and this class does nothing but hand both to a
`QTextBrowser`. Everything that could be wrong in an interesting way -- table
structure, escaping, math placement, the approximate-method banner -- is
decided and tested there, without a window.

Images are registered with `QTextDocument.addResource` rather than embedded as
`data:` URIs, because `loadResource` does not handle `data:` by default. That
also disposes of the Tk viewer's ugliest detail: a list monkey-patched onto the
window purely to keep `PhotoImage` objects alive against garbage collection
while Tk was still displaying them.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QImage, QTextDocument
from PySide6.QtWidgets import QTextBrowser, QVBoxLayout, QWidget

from .. import richtext

__all__ = ["TheoryViewer", "render_into"]


def render_into(
    browser: QTextBrowser, page, *, device_pixel_ratio: float, text_color: str
) -> None:
    """Render a page into an existing `QTextBrowser`, images and all.

    Shared by :class:`TheoryViewer` and :class:`~.setup_tab.SetupTab` --
    both are "one markdown-ish page, typeset, in a browser," and the second
    use is exactly the signal that this was worth factoring out rather than
    copied.
    """
    rendered = richtext.to_html(
        page, device_pixel_ratio=device_pixel_ratio, color=text_color
    )
    document = browser.document()
    for image in rendered.images:
        document.addResource(
            QTextDocument.ResourceType.ImageResource,
            QUrl(image.key),
            QImage.fromData(image.png, "PNG"),
        )
    browser.setHtml(rendered.html)


class TheoryViewer(QWidget):
    """One page, typeset."""

    def __init__(self, page, parent: QWidget | None = None) -> None:
        # A window in its own right rather than a modal dialog: reading the
        # theory while adjusting the grating is the point.
        super().__init__(parent)
        self.setWindowFlag(Qt.WindowType.Window)
        self.setWindowTitle(page.title)
        self.resize(860, 720)

        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(False)
        # Relative links between doc pages point at files on disk that this
        # viewer has no route to open. Following them would empty the window.
        self._browser.setOpenLinks(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._browser)

        self.set_page(page)

    def set_page(self, page) -> None:
        """Render a page into the browser, images and all."""
        render_into(
            self._browser,
            page,
            device_pixel_ratio=self.devicePixelRatioF(),
            text_color=self._text_color(),
        )

    def _text_color(self) -> str:
        """Rasterise math in the palette's own text colour.

        Equations are images, so they do not follow a stylesheet the way the
        surrounding prose does. Black math on a dark theme would be unreadable
        for no reason other than that it was baked in.
        """
        return self.palette().windowText().color().name()
