r"""Theory pages as rich text: markdown structure *and* typeset math.

M7 got the equations rendering and deliberately stopped there -- tables and
bullet lists stayed as literal pipes and dashes, because real column layout in a
`tkinter.Text` widget -- the toolkit at the time -- is a great deal of
engineering for content that is readable enough as-is. That trade has since gone bad: `docs/theory/scalar.md`
carries the three-way tension table and the shallow-limit evidence table, and
`docs/conventions.md` carries the angle reconciliation and the errata table.
The tables *are* the content now, in the one place a user goes to understand
what the solver does.

Qt's rich-text engine parses a small HTML subset, which makes this cheap. So
this module converts a page to HTML, and the widget layer only displays it.

Three things worth knowing before editing:

**Math is split out before anything is escaped.** Four math spans in the real
docs contain `<` or `>` (`$t < w$`, `$\zeta < \theta_c$`, `$\lambda > 32 ...$`,
`$k > 0$`). Escaping first would hand mathtext the string `t &lt; w`, which it
renders literally. So :func:`~gratinglab.gui.mathtext.split_segments` runs on
the raw source and only the *text* spans are escaped -- the LaTeX never is.

**Images are referenced by key, not embedded.** `QTextDocument.loadResource`
does not handle `data:` URLs, so each image gets a `gratinglab-math:N` src and
the viewer registers the bytes with `addResource`. That also keeps this module
free of base64 padding it would only have to undo.

**A render failure still degrades rather than crashing**, exactly as the Tk
viewer did: the raw `$...$` source appears in a distinct colour, so a future
edit using an unsupported macro is diagnosable at a glance instead of taking
down the page.

Deliberately not supported, because the two real pages do not use them and
guessing is worse than a visible gap: nested lists, ordered lists inside
blockquotes, inline HTML, reference-style links, and setext (`===`) headings.
A standalone `---` is always a horizontal rule here, never a setext underline.
"""

from __future__ import annotations

import html
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import TYPE_CHECKING, Callable, Sequence

from . import mathtext

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from .docs import TheoryPage

__all__ = [
    "MathImage",
    "RichPage",
    "render_markdown",
    "to_html",
    "png_size",
    "math_dpi",
    "menu_label",
    "INLINE_DPI",
    "DISPLAY_DPI",
]

#: Base rasterisation DPI. Display equations render larger than inline ones,
#: roughly as a printed page distinguishes a headline equation from a symbol in
#: running text. Both are scaled by the device pixel ratio at render time.
INLINE_DPI = 130
DISPLAY_DPI = 170

#: Colour for a math span mathtext could not parse. Same signal the Tk viewer's
#: `unrendered` tag gave.
UNRENDERED_COLOR = "#a5370d"

#: Colour for the approximate-method banner. Deliberately its own constant even
#: though it currently matches: "this solver approximates" and "this equation
#: failed to render" are unrelated facts, and a test that cannot tell them
#: apart is a test that will one day pass for the wrong reason.
BANNER_COLOR = "#a5370d"

#: `# Title` through `###### h6`. The Tk viewer matched `#{2,6}` and so left a
#: literal `#` on the first line of every page, both of which open with a
#: level-1 heading.
_HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")
_RULE_RE = re.compile(r"^(?:-{3,}|_{3,}|\*{3,})$")
_BULLET_RE = re.compile(r"^\s*[-*+]\s+(.*)$")
_TABLE_DELIMITER_RE = re.compile(r"^\|[\s:|-]+\|$")
_FENCE_RE = re.compile(r"^```")

_BOLD_RE = re.compile(r"\*\*(.+?)\*\*")
_ITALIC_RE = re.compile(r"(?<![\w*])\*([^*\s][^*]*?)\*(?![\w*])")
_CODE_RE = re.compile(r"`([^`]+)`")
_LINK_RE = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
_PLACEHOLDER_RE = re.compile(r"\x00(\d+)\x00")


@dataclass(frozen=True, slots=True)
class MathImage:
    """One rasterised equation, ready to register with a `QTextDocument`.

    ``width`` and ``height`` are *logical* pixels -- the PNG's own dimensions
    divided by the device pixel ratio it was rendered at -- so that a 2x image
    is downsampled into a correctly sized box on a Retina panel rather than
    rendering at double size.
    """

    key: str
    png: bytes
    width: int
    height: int


@dataclass(frozen=True, slots=True)
class RichPage:
    """A page's HTML, plus the images its `<img src>` keys refer to."""

    html: str
    images: tuple[MathImage, ...] = ()


@dataclass
class _Collector:
    """Accumulates images while the inline pass walks a page."""

    inline_dpi: int
    display_dpi: int
    device_pixel_ratio: float
    color: str
    render: Callable[..., bytes | None]
    images: list[MathImage] = field(default_factory=list)

    def add(self, latex: str, *, display: bool) -> str:
        """Render one math span; return the HTML that stands for it."""
        base = self.display_dpi if display else self.inline_dpi
        dpi = max(1, round(base * self.device_pixel_ratio))
        png = self.render(latex, dpi=dpi, color=self.color)

        if png is None:
            delimiter = "$$" if display else "$"
            source = html.escape(f"{delimiter}{latex}{delimiter}")
            return f'<span style="color:{UNRENDERED_COLOR}">{source}</span>'

        pixels = png_size(png)
        if pixels is None:  # pragma: no cover - only a truncated PNG gets here
            return ""
        width, height = pixels
        key = f"gratinglab-math:{len(self.images)}"
        self.images.append(
            MathImage(
                key=key,
                png=png,
                width=max(1, round(width / self.device_pixel_ratio)),
                height=max(1, round(height / self.device_pixel_ratio)),
            )
        )
        image = f'<img src="{key}" width="{self.images[-1].width}" height="{self.images[-1].height}">'
        return f"<br>{image}<br>" if display else image


def render_markdown(
    markdown: str,
    *,
    inline_dpi: int = INLINE_DPI,
    display_dpi: int = DISPLAY_DPI,
    device_pixel_ratio: float = 1.0,
    color: str = "black",
    render: Callable[..., bytes | None] | None = None,
) -> RichPage:
    """Convert markdown with `$...$` math to Qt-displayable HTML.

    ``render`` is injectable so the whole pipeline can be tested without
    rasterising anything, including the degradation path where it returns
    ``None``. It defaults to a cached
    :func:`~gratinglab.gui.mathtext.render_math_png`; a page carries roughly
    forty equations and the cache is what keeps reopening the viewer instant.
    """
    collector = _Collector(
        inline_dpi=inline_dpi,
        display_dpi=display_dpi,
        device_pixel_ratio=device_pixel_ratio,
        color=color,
        render=render if render is not None else _cached_render,
    )
    parts = [_block_html(block, collector) for block in _blocks(markdown)]
    return RichPage(html="\n".join(p for p in parts if p), images=tuple(collector.images))


def to_html(
    page: "TheoryPage",
    **kwargs,
) -> RichPage:
    """Render a theory page, banner included.

    A solver whose `Capabilities.rigorous` is False is announced at the top,
    matching the project's practice of surfacing an approximation rather than
    letting a reader discover it in section 7.
    """
    rendered = render_markdown(page.text, **kwargs)
    if page.rigorous:
        return rendered
    banner = (
        f'<p style="color:{BANNER_COLOR}"><b>'
        "⚠ Approximate method — see Limits below."
        "</b></p>"
    )
    return RichPage(html=banner + "\n" + rendered.html, images=rendered.images)


# -- block structure -----------------------------------------------------


@dataclass(frozen=True, slots=True)
class _Block:
    kind: str
    lines: tuple[str, ...] = ()
    rows: tuple[tuple[str, ...], ...] = ()
    level: int = 0


def _blocks(markdown: str) -> tuple[_Block, ...]:
    """Split a page into block-level structures.

    One pass, no lookbehind: every block kind is identified by its first line,
    which is true of the markdown these docs actually use.
    """
    lines = markdown.split("\n")
    blocks: list[_Block] = []
    i, n = 0, len(lines)

    while i < n:
        line = lines[i]

        if not line.strip():
            i += 1
            continue

        if _FENCE_RE.match(line):
            i += 1
            body: list[str] = []
            while i < n and not _FENCE_RE.match(lines[i]):
                body.append(lines[i])
                i += 1
            i += 1  # closing fence, or end of page
            blocks.append(_Block("code", tuple(body)))
            continue

        heading = _HEADING_RE.match(line)
        if heading:
            blocks.append(
                _Block("heading", (heading.group(2),), level=len(heading.group(1)))
            )
            i += 1
            continue

        if _RULE_RE.match(line.strip()):
            blocks.append(_Block("rule"))
            i += 1
            continue

        if line.lstrip().startswith("|"):
            rows: list[tuple[str, ...]] = []
            while i < n and lines[i].lstrip().startswith("|"):
                if not _TABLE_DELIMITER_RE.match(lines[i].strip()):
                    rows.append(_table_cells(lines[i]))
                i += 1
            blocks.append(_Block("table", rows=tuple(rows)))
            continue

        if line.lstrip().startswith(">"):
            quoted: list[str] = []
            while i < n and lines[i].strip():
                stripped = lines[i].lstrip()
                quoted.append(
                    stripped[1:].strip() if stripped.startswith(">") else stripped
                )
                i += 1
            blocks.append(_Block("quote", (" ".join(quoted),)))
            continue

        bullet = _BULLET_RE.match(line)
        if bullet:
            items: list[str] = []
            while i < n and lines[i].strip():
                match = _BULLET_RE.match(lines[i])
                if match:
                    items.append(match.group(1))
                elif items:
                    # An indented continuation of the item above. Markdown
                    # soft-wraps, and these docs wrap at about 80 columns.
                    items[-1] += " " + lines[i].strip()
                else:  # pragma: no cover - unreachable: the loop starts on a bullet
                    break
                i += 1
            blocks.append(_Block("bullets", tuple(items)))
            continue

        paragraph: list[str] = []
        while i < n and lines[i].strip() and not _starts_a_block(lines[i]):
            paragraph.append(lines[i].strip())
            i += 1
        blocks.append(_Block("paragraph", (" ".join(paragraph),)))

    return tuple(blocks)


def _starts_a_block(line: str) -> bool:
    """Whether this line begins something other than more paragraph text."""
    stripped = line.lstrip()
    return bool(
        _HEADING_RE.match(line)
        or _RULE_RE.match(line.strip())
        or _FENCE_RE.match(line)
        or stripped.startswith("|")
        or stripped.startswith(">")
        or _BULLET_RE.match(line)
    )


def _table_cells(line: str) -> tuple[str, ...]:
    """Cells of one table row, outer pipes discarded."""
    return tuple(cell.strip() for cell in line.strip().strip("|").split("|"))


def _block_html(block: _Block, collector: _Collector) -> str:
    if block.kind == "heading":
        level = min(block.level, 6)
        return f"<h{level}>{_inline(block.lines[0], collector)}</h{level}>"

    if block.kind == "paragraph":
        return f"<p>{_inline(block.lines[0], collector)}</p>"

    if block.kind == "bullets":
        items = "".join(f"<li>{_inline(item, collector)}</li>" for item in block.lines)
        return f"<ul>{items}</ul>"

    if block.kind == "quote":
        return f"<blockquote><p>{_inline(block.lines[0], collector)}</p></blockquote>"

    if block.kind == "rule":
        return "<hr>"

    if block.kind == "code":
        # Escaped and left literal: code is the one place markdown and math
        # syntax must not be interpreted.
        body = html.escape("\n".join(block.lines))
        return f'<pre style="white-space:pre">{body}</pre>'

    if block.kind == "table":
        if not block.rows:
            return ""
        header, *body = block.rows
        cells = "".join(f"<th>{_inline(c, collector)}</th>" for c in header)
        out = [f"<tr>{cells}</tr>"]
        for row in body:
            cells = "".join(f"<td>{_inline(c, collector)}</td>" for c in row)
            out.append(f"<tr>{cells}</tr>")
        return (
            '<table border="1" cellpadding="4" cellspacing="0">'
            + "".join(out)
            + "</table>"
        )

    raise AssertionError(f"unhandled block kind {block.kind!r}")  # pragma: no cover


# -- inline structure ----------------------------------------------------


def _inline(raw: str, collector: _Collector) -> str:
    """Render one run of markdown text, math spans included.

    Math is separated **before** escaping: several math spans in the real docs
    contain `<` or `>`, and handing mathtext `&lt;` would render it literally.
    """
    out: list[str] = []
    for segment in mathtext.split_segments(raw):
        if segment.kind == "math":
            out.append(collector.add(segment.content, display=segment.display))
        else:
            out.append(_markup(html.escape(segment.content)))
    return "".join(out)


def _markup(escaped: str) -> str:
    """Bold, italic, code and links, applied to already-escaped text.

    Order matters: code spans are extracted first so that markup characters
    inside them stay literal, and links before emphasis so a link's text can be
    emphasised without the URL being scanned for asterisks.
    """
    placeholders: list[str] = []

    def stash(html_fragment: str) -> str:
        placeholders.append(html_fragment)
        return f"\x00{len(placeholders) - 1}\x00"

    text = _CODE_RE.sub(lambda m: stash(f"<code>{m.group(1)}</code>"), escaped)
    text = _LINK_RE.sub(
        lambda m: stash(f'<a href="{m.group(2)}">{_emphasis(m.group(1))}</a>'), text
    )
    text = _emphasis(text)

    # Repeatedly, because placeholders nest: scalar.md's first line is
    # [`scalar.py`](../../src/gratinglab/solvers/scalar.py), so the link's
    # placeholder contains the code span's. `re.sub` does not rescan what it
    # substituted, and a single pass left the inner one as a raw null byte.
    # Terminates because a placeholder can only ever contain lower indices.
    while _PLACEHOLDER_RE.search(text):
        text = _PLACEHOLDER_RE.sub(lambda m: placeholders[int(m.group(1))], text)
    return text


def _emphasis(text: str) -> str:
    text = _BOLD_RE.sub(r"<b>\1</b>", text)
    return _ITALIC_RE.sub(r"<i>\1</i>", text)


# -- helpers -------------------------------------------------------------


def png_size(data: bytes) -> tuple[int, int] | None:
    """Pixel dimensions from a PNG's IHDR chunk, without a decoder.

    Layout is fixed by the format: an 8-byte signature, a 4-byte chunk length,
    the literal ``IHDR``, then width and height as big-endian uint32. Reading
    them directly keeps the size calculation pure, so the HiDPI arithmetic is
    testable without Qt or a display.
    """
    if len(data) < 24 or data[12:16] != b"IHDR":
        return None
    return (
        int.from_bytes(data[16:20], "big"),
        int.from_bytes(data[20:24], "big"),
    )


def math_dpi(font_pixel_size: float, *, display: bool = False) -> int:
    """DPI at which mathtext output matches text of the given size.

    matplotlib renders mathtext at ``rcParams["font.size"]`` (10 pt by
    default), so 72 dpi would give ten-pixel-tall glyphs. Scaling by the
    viewer's actual font size keeps equations the same visual weight as the
    prose around them at any zoom or platform default, which the fixed 130/170
    constants only manage on the machine they were chosen on.
    """
    base = 72.0 * font_pixel_size / 10.0
    return max(1, round(base * (1.3 if display else 1.0)))


def menu_label(title: str, suffix: str = "") -> str:
    """A menu entry's text, with Qt's mnemonic character neutralised.

    Qt reads ``&`` in an action's text as "underline the next character", so
    the page titled "Grating Geometry & Conventions" would appear as "Grating
    Geometry Conventions" with a stray underline. Doubling escapes it.
    """
    return (title + suffix).replace("&", "&&")


@lru_cache(maxsize=512)
def _cached_render(latex: str, *, dpi: int, color: str) -> bytes | None:
    """Rasterise, remembering. Reopening the viewer re-renders every equation
    on the page otherwise -- roughly forty of them, and more at 2x."""
    return mathtext.render_math_png(latex, dpi=dpi, color=color)
