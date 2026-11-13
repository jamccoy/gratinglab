r"""Render LaTeX math from a theory page, without a markdown dependency.

A theory page's `$$\Phi_m(x) = ...$$` should not reach a reader as literal
source with its backslashes. This module is the pure half of the solution. It
splits a markdown string into text and math spans and rasterizes the math
spans to PNG bytes with :func:`matplotlib.mathtext.math_to_image`. That
function is already a dependency (it draws every axis label in the plots) and
needs no live display, so the whole pipeline down to PNG bytes can be tested
headlessly. Only the final step, turning those bytes into an embeddable image,
lives in :mod:`gratinglab.gui.richtext` and
:mod:`gratinglab.gui.qt.theory_viewer`, where a real Qt document is needed.

Two points should be understood before this file is changed.

First, `\boxed{...}` is not supported by mathtext; testing confirmed that
it raises `ValueError` wrapping a `ParseFatalException`. The two headline
equations of `scalar.md` use it. :func:`strip_boxed` removes the wrapper by
counting braces instead of by regex, because the boxed content itself
contains nested braces (`\boxed{\;G_m = \int_0^1 e^{i\Phi_m(t)}\, ...\;}`),
and a naive `\boxed\{(.*)\}` pattern would either stop at the first inner `}`
or, if greedy, continue past the intended end. The visual rectangle is not
reproduced. Mathtext has no boxed macro, and reconstructing one would require
splitting a single display equation into several sub-images for a cosmetic
rule. The equation still renders correctly on its own centered line.

Second, a render failure does not crash the viewer. `render_math_png`
returns `None` instead of propagating the parser's exception, so an
unsupported macro introduced in a future edit appears as raw source for that
one span, where it is easy to diagnose, and the rest of the page still
renders.
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
    the string is treated as literal text instead of raising. A markdown
    source with a typo should therefore show the stray dollar sign, and should
    not crash a viewer whose purpose is to improve how pages render.
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

    The wrapper is found by counting braces. A regex is not used because the
    content inside ``\boxed{}`` routinely contains its own nested braces
    (subscripts, ``\left\{`` groups, fractions), so a pattern like
    ``\\boxed\{(.*)\}`` cannot reliably find the matching close brace. More
    than one ``\boxed{}`` in the same string is handled, as in the
    side-by-side pair in `scalar.md` section 3.

    A ``\boxed`` with no matching close brace is left as it is: the string is
    passed through unchanged from that point, so the function does not raise
    on malformed input that it did not create.
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
    latex: str,
    *,
    dpi: int = 140,
    color: str = "black",
    axis_pt: float | None = None,
    min_height_pt: float = 0.0,
) -> bytes | None:
    r"""Rasterize a LaTeX string (no `$` delimiters) to PNG bytes.

    Returns ``None`` on any parser failure, such as an unsupported macro or
    malformed LaTeX, instead of raising, so the caller can fall back to
    showing the raw source for that one span. ``\boxed{}`` wrappers are
    stripped first (see the module docstring). Testing confirmed that
    everything else the project's docs use (`\operatorname`, `\text`,
    `\qquad`, Greek letters, `\left[...\right]`, sub/superscripts) is
    supported by mathtext directly.

    A display equation split across source lines for readability (e.g. a long
    ``$$...$$`` wrapped onto a second, indented line) carries a literal
    newline and leading whitespace into its LaTeX source once the markdown
    delimiters are stripped. Ordinary TeX ignores such whitespace, but testing
    confirmed that mathtext's parser does not tolerate an embedded ``\n``. Any
    run of whitespace is therefore collapsed to a single space before parsing,
    which matches what plain LaTeX would do with the same source.

    ``axis_pt`` is for an equation set inside a line of text. Qt places an
    image aligned ``vertical-align: middle`` with its centre a quarter of the
    font's x-height above the baseline, so given that height in points of the
    10 pt render, space is added above or below until the image's centre lands
    there exactly when the equation's own baseline sits on the text's. Without
    it the image's bottom edge sits on the baseline, and every symbol with a
    descender rides above the line.

    That rule holds only for an image at least about three x-heights tall;
    Qt sets a shorter one higher, by an amount that depends on its height
    (measured with PySide6 6.11 on 2026-11-12, by scanning solid images from
    8 to 32 px against a 14 px font). ``min_height_pt`` pads the image evenly
    above and below up to that height, which keeps the centre where it was.
    """
    import io
    import re

    import matplotlib

    matplotlib.use("Agg")  # no display needed; safe to call repeatedly
    from matplotlib.mathtext import math_to_image

    normalized = re.sub(r"\s+", " ", strip_boxed(latex)).strip()

    buffer = io.BytesIO()
    try:
        if axis_pt is None:
            math_to_image(f"${normalized}$", buffer, dpi=dpi, format="png", color=color)
            return buffer.getvalue()
        return _inline_png(normalized, dpi, color, axis_pt, min_height_pt)
    except ValueError:
        return None


def _inline_png(
    latex: str, dpi: int, color: str, axis_pt: float, min_height_pt: float
) -> bytes:
    """An inline equation as a transparent PNG, padded about its baseline.

    The baseline is measured rather than taken from mathtext's box depth,
    which is the layout box and not where the glyphs land: rasterised at
    101 dpi the two differ by a pixel in either direction depending on the
    expression. So the expression is rendered a second time with an upright
    ``x`` after it, which sits on the baseline, and the bottom of that ``x``
    is read off. The two renders share one vertical frame whenever the ``x``
    is no taller and no deeper than the expression, which is any expression
    holding a letter or digit; otherwise the box depth is used.
    """
    import io

    import numpy as np
    from matplotlib.colors import to_rgb
    from matplotlib.image import imsave
    from matplotlib.mathtext import MathTextParser

    parser = MathTextParser("agg")
    ink = np.asarray(parser.parse(f"${latex}$", dpi=dpi).image)
    probe = parser.parse(f"${latex}" + r"\mathrm{x}$", dpi=dpi)
    probe_ink = np.asarray(probe.image)
    rows = ink.shape[0]
    if probe_ink.shape[0] == rows:
        # The rightmost couple of points of the x: its right leg ends on
        # the baseline, and nothing of the expression reaches that far.
        columns = np.flatnonzero(probe_ink.max(axis=0))
        leg = probe_ink[:, columns[-1] - max(1, round(2 * dpi / 72)) + 1 : columns[-1] + 1]
        baseline = int(np.flatnonzero(leg.max(axis=1) > 100)[-1]) + 1
    else:
        rendered = parser.parse(f"${latex}$", dpi=dpi)
        baseline = int(rendered.height - rendered.depth) - 1

    # Centre on the axis: (above - below) / 2 must equal the axis height.
    scale = dpi / 72.0
    above, below = float(baseline), float(rows - baseline)
    excess = (above - below) - 2 * axis_pt * scale
    if excess > 0:
        below += excess
    else:
        above -= excess
    short = min_height_pt * scale - (above + below)
    if short > 0:
        above += short / 2
        below += short / 2
    top = round(above) - baseline
    bottom = round(below) - (rows - baseline)

    alpha = np.pad(ink, ((top, bottom), (0, 0))).astype(float) / 255.0
    rgba = np.empty(alpha.shape + (4,))
    rgba[..., :3] = to_rgb(color)
    rgba[..., 3] = alpha
    buffer = io.BytesIO()
    imsave(buffer, rgba, format="png")
    return buffer.getvalue()
