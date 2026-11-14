"""Every equation in gratinglab's pages renders.

The renderer is endstation's (``endstation.docs.mathtext``) and is tested
there. What is tested here is the pages: the spans are parsed out of the real
files, so an edit that adds a macro mathtext cannot draw fails here rather
than showing up as raw source in the app.
"""

import pytest
from endstation.docs.mathtext import render_math_png, split_segments

from gratinglab.gui.docs import find_theory_root, general_pages, theory_pages

THEORY_ROOT = find_theory_root()
requires_theory_docs = pytest.mark.skipif(
    THEORY_ROOT is None, reason="docs/theory not found from this checkout"
)


def _failures(text):
    spans = [s.content for s in split_segments(text) if s.kind == "math"]
    return spans, [s for s in spans if render_math_png(s) is None]


@requires_theory_docs
def test_every_math_span_in_scalar_md_renders():
    spans, failures = _failures((THEORY_ROOT / "scalar.md").read_text(encoding="utf-8"))
    assert len(spans) > 10, "sanity check: did the file load correctly?"
    assert not failures, f"{len(failures)} span(s) failed to render: {failures}"


@requires_theory_docs
@pytest.mark.parametrize(
    "page",
    [p for p in (*theory_pages(), *general_pages()) if p.available],
    ids=lambda p: p.name,
)
def test_every_math_span_in_every_help_page_renders(page):
    _, failures = _failures(page.text)
    assert not failures, f"{page.name}: {len(failures)} span(s) failed: {failures}"
