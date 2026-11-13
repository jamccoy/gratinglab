"""The read-only viewer for one theory or reference page.

The viewer is thin by design. :mod:`gratinglab.gui.richtext` turns the page
into HTML and a list of rasterised equations, and this class only passes both
to a `QTextBrowser`. Everything that could be wrong in a non-trivial way
(table structure, escaping, math placement, the approximate-method banner) is
decided and tested in that module, without a window.

Images are registered with `QTextDocument.addResource` instead of being
embedded as `data:` URIs, because `loadResource` does not handle `data:` by
default. This also removes the least elegant detail of the Tk viewer: a list
monkey-patched onto the window solely to keep `PhotoImage` objects alive
against garbage collection while Tk was still displaying them.

:class:`ReadingBrowser` is the page itself: a slightly larger font than the
widgets around it, and a text column kept to a readable width and centred,
so a wide window adds margin rather than longer lines. The math is sized
from that font, so an equation in a sentence matches the words beside it.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import (
    QDesktopServices,
    QFontMetricsF,
    QImage,
    QTextBlockFormat,
    QTextCursor,
    QTextDocument,
    QTextImageFormat,
)
from PySide6.QtWidgets import QTextBrowser, QVBoxLayout, QWidget

from .. import richtext
from ..docs import resolve_link

__all__ = ["ReadingBrowser", "TheoryViewer", "render_into"]

#: Body text relative to the application font. The pages are long prose, and
#: the 9 pt that suits a button is small for reading paragraphs of it.
READING_SCALE = 1.15

#: Widest the text column grows, in ems of the reading font: about 100
#: characters of prose a line. A display equation wider than that widens the
#: column to fit rather than scrolling sideways.
READING_WIDTH_EMS = 50

#: Space either side of the column, in pixels, when the window is narrower
#: than the column.
READING_MARGIN = 24

#: Space around tables, which Qt otherwise sets flush against the paragraph
#: above them.
_STYLE = "table { margin-top: 6px; margin-bottom: 10px; } th { text-align: left; }"

#: Extra space between lines, as a fraction of the font's pixel size. Added
#: as a fixed distance rather than CSS ``line-height: 135%``, which scales
#: each line's own height and so opens a wide gap around every line holding
#: an equation.
LEADING = 0.3


class ReadingBrowser(QTextBrowser):
    """A `QTextBrowser` laid out for reading a page rather than a log."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        font = self.font()
        font.setPointSizeF(font.pointSizeF() * READING_SCALE)
        self.setFont(font)
        self.document().setDefaultFont(font)
        self.document().setDefaultStyleSheet(_STYLE)
        self._widest_image = 0.0

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        super().resizeEvent(event)
        self._centre_column()

    def pixel_size(self) -> float:
        """The reading font's size in pixels.

        From the point size and the screen's DPI rather than `QFontInfo`,
        which reports the matched font and so reads near zero where no system
        fonts are installed, as on the offscreen CI runs.
        """
        return self.font().pointSizeF() * self.logicalDpiY() / 72.0

    def _centre_column(self) -> None:
        column = max(self.pixel_size() * READING_WIDTH_EMS, self._widest_image)
        spare = self.viewport().width() - column
        side = max(READING_MARGIN, spare / 2)
        root = self.document().rootFrame()
        frame = root.frameFormat()
        if frame.leftMargin() == side and frame.rightMargin() == side:
            return
        frame.setLeftMargin(side)
        frame.setRightMargin(side)
        frame.setTopMargin(READING_MARGIN / 2)
        frame.setBottomMargin(READING_MARGIN)
        root.setFrameFormat(frame)

    def setHtml(self, text: str) -> None:  # noqa: N802 - Qt's spelling
        # setHtml replaces the root frame's format along with the content.
        super().setHtml(text)
        leading = QTextBlockFormat()
        leading.setLineHeight(
            round(self.pixel_size() * LEADING),
            QTextBlockFormat.LineHeightTypes.LineDistanceHeight.value,
        )
        cursor = QTextCursor(self.document())
        cursor.select(QTextCursor.SelectionType.Document)
        cursor.mergeBlockFormat(leading)
        self._widest_image = _widest_image(self.document())
        self._centre_column()


def _widest_image(document: QTextDocument) -> float:
    """The widest image in a document, in pixels, tables included."""
    widest = 0.0
    block = document.begin()
    while block.isValid():
        fragments = block.begin()
        while not fragments.atEnd():
            fmt = fragments.fragment().charFormat()
            if fmt.isImageFormat():
                widest = max(widest, QTextImageFormat(fmt).width())
            fragments += 1
        block = block.next()
    return widest


def render_into(
    browser: QTextBrowser, page, *, device_pixel_ratio: float, text_color: str
) -> None:
    """Render a page into an existing `QTextBrowser`, images and all.

    It is shared by :class:`TheoryViewer` and :class:`~.setup_tab.SetupTab`,
    both of which show "one markdown-ish page, typeset, in a browser," and the
    second use is what justified factoring this function out instead of
    copying it.
    """
    font = browser.document().defaultFont()
    pixels = font.pointSizeF() * browser.logicalDpiY() / 72.0
    rendered = richtext.to_html(
        page,
        device_pixel_ratio=device_pixel_ratio,
        color=text_color,
        inline_dpi=richtext.math_dpi(pixels),
        display_dpi=richtext.math_dpi(pixels, display=True),
        **richtext.inline_metrics(pixels, QFontMetricsF(font).xHeight()),
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
            text_color=self._text_color(),
        )

    def _text_color(self) -> str:
        """Rasterise math in the palette's own text colour.

        Equations are images, so they do not follow a stylesheet the way the
        surrounding prose does. Black math on a dark theme would be unreadable
        only because its colour had been fixed at render time.
        """
        return self.palette().windowText().color().name()
