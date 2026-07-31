r"""Rendering LaTeX math from a theory page, without a markdown dependency.

The theory viewer previously showed `$$\Phi_m(x) = ...$$` as literal source --
backslashes and all -- because it was a plain :class:`tkinter.Text` insertion.
This module is the fix's pure half: splitting a markdown string into text and
math spans, and rasterizing math spans to PNG bytes via
:func:`matplotlib.mathtext.math_to_image`, which is already a dependency (it is
what draws every axis label in the plots) and needs no live display to run --
the whole pipeline down to PNG bytes is headlessly testable. Only the final
step of turning those bytes into a `tkinter.PhotoImage` and inserting it into a
`Text` widget lives in :mod:`gratinglab.gui.app`, where it needs a real Tk root.

Two things worth knowing before touching this file:

**`\boxed{...}` is not supported by mathtext** (confirmed by testing --
it raises `ValueError` wrapping a `ParseFatalException`). `scalar.md`'s two
headline equations use it. :func:`strip_boxed` removes the wrapper by counting
braces rather than by regex, because the boxed content itself contains nested
braces (`\boxed{\;G_m = \int_0^1 e^{i\Phi_m(t)}\, ...\;}`) and a naive
`\boxed\{(.*)\}` pattern would either stop at the first inner `}` or, greedy,
swallow past the intended end. The visual rectangle is not reproduced --
mathtext has no boxed macro, and rebuilding one would mean splitting a single
display equation into multiple sub-images for a cosmetic rule. The equation
still renders correctly and reads as its own centered line, which is the part
that matters.

**A render failure degrades, it does not crash the viewer.** `render_math_png`
returns `None` rather than propagating the parser's exception, so one
unsupported macro in a future edit shows as visible raw source for that one
span -- diagnosable at a glance -- rather than taking down the whole page.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

__all__ = ["Segment", "split_segments", "strip_boxed", "render_math_png"]


@dataclass(frozen=True, slots=True)
class Segment:
    """One piece of a theory page: either plain text or a math span.

    Attributes
    ----------
    kind
        ``"text"`` or ``"math"``.
    content
        For ``"text"``, the literal markdown text. For ``"math"``, the LaTeX
        source with its `$` or `$$` delimiters already stripped.
    display
        Only meaningful for ``"math"``: ``True`` for `$$...$$` (its own
        centered line), ``False`` for inline `$...$`.
    """

    kind: Literal["text", "math"]
    content: str
    display: bool = False


def split_segments(markdown: str) -> tuple[Segment, ...]:
    r"""Split markdown into alternating text and math segments.

    Recognises `$$...$$` (display) and `$...$` (inline). A backslash-escaped
    `\$` is never a delimiter and is emitted as a literal `$` in the
    surrounding text segment. An odd, unterminated `$` (or `$$`) at the end of
    the string is treated as literal text rather than raising -- a markdown
    source with a typo should degrade to showing the stray dollar sign, not
    crash the viewer that is trying to fix rendering problems.
    """
    segments: list[Segment] = []
    buffer: list[str] = []
    i, n = 0, len(markdown)

    def flush_text() -> None:
        if buffer:
            segments.append(Segment(kind="text", content="".join(buffer)))
            buffer.clear()

    while i < n:
        char = markdown[i]

        if char == "\\" and i + 1 < n and markdown[i + 1] == "$":
            buffer.append("$")
            i += 2
            continue

        if char == "$":
            display = markdown[i : i + 2] == "$$"
            delimiter = "$$" if display else "$"
            start = i + len(delimiter)
            end = _find_unescaped(markdown, delimiter, start)
            if end == -1:
                # No closing delimiter: treat the rest as literal text rather
                # than silently dropping content or raising.
                buffer.append(markdown[i:])
                i = n
                continue
            flush_text()
            segments.append(
                Segment(kind="math", content=markdown[start:end], display=display)
            )
            i = end + len(delimiter)
            continue

        buffer.append(char)
        i += 1

    flush_text()
    return tuple(segments)


def _find_unescaped(text: str, delimiter: str, start: int) -> int:
    """Index of the next ``delimiter`` in ``text`` from ``start``, skipping
    escaped occurrences. -1 if none."""
    i = start
    while True:
        i = text.find(delimiter, i)
        if i == -1:
            return -1
        if i > 0 and text[i - 1] == "\\":
            i += len(delimiter)
            continue
        return i


def strip_boxed(latex: str) -> str:
    r"""Remove every ``\boxed{...}`` wrapper, keeping its inner content.

    Brace-counting, not regex: the content inside ``\boxed{}`` routinely
    contains its own nested braces (subscripts, ``\left\{`` groups, fractions),
    so a pattern like ``\\boxed\{(.*)\}`` cannot reliably find the matching
    close brace. Handles more than one ``\boxed{}`` in the same string, as in
    the side-by-side pair in `scalar.md` section 3.

    A ``\boxed`` with no matching close brace is left as-is -- the string is
    passed through unchanged from that point, rather than raising on malformed
    input this function did not create.
    """
    marker = "\\boxed{"
    result: list[str] = []
    i, n = 0, len(latex)

    while i < n:
        start = latex.find(marker, i)
        if start == -1:
            result.append(latex[i:])
            break

        result.append(latex[i:start])
        depth = 1
        j = start + len(marker)
        while j < n and depth > 0:
            if latex[j] == "{":
                depth += 1
            elif latex[j] == "}":
                depth -= 1
            j += 1

        if depth != 0:
            # No matching close brace: not well-formed \boxed{}; leave the
            # remainder untouched rather than guessing.
            result.append(latex[start:])
            i = n
            break

        result.append(latex[start + len(marker) : j - 1])
        i = j

    return "".join(result)


def render_math_png(
    latex: str, *, dpi: int = 140, color: str = "black"
) -> bytes | None:
    r"""Rasterize a LaTeX string (no `$` delimiters) to PNG bytes.

    Returns ``None`` on any parser failure -- an unsupported macro, malformed
    LaTeX -- instead of raising, so the caller can fall back to showing the raw
    source for that one span. ``\boxed{}`` wrappers are stripped first (see
    module docstring); everything else the project's docs use
    (`\operatorname`, `\text`, `\qquad`, Greek letters, `\left[...\right]`,
    sub/superscripts) is supported by mathtext directly, confirmed by testing.

    A display equation split across source lines for readability (e.g. a long
    ``$$...$$`` wrapped onto a second, indented line) carries a literal
    newline and leading whitespace into its LaTeX source once markdown
    delimiters are stripped. Ordinary TeX ignores such whitespace, but
    mathtext's parser does not tolerate an embedded ``\n`` -- confirmed by
    testing -- so any run of whitespace is collapsed to a single space before
    parsing, matching what plain LaTeX would do with the same source.
    """
    import io
    import re

    import matplotlib

    matplotlib.use("Agg")  # no display needed; safe to call repeatedly
    from matplotlib.mathtext import math_to_image

    normalized = re.sub(r"\s+", " ", strip_boxed(latex)).strip()

    buffer = io.BytesIO()
    try:
        math_to_image(f"${normalized}$", buffer, dpi=dpi, format="png", color=color)
    except ValueError:
        return None
    return buffer.getvalue()
