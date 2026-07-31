r"""LaTeX splitting and rendering.

Headless -- matplotlib's Agg backend needs no display, so this whole pipeline
down to PNG bytes is testable in CI. Only turning those bytes into a
tkinter.PhotoImage needs a live root, which is tested separately in
test_gui_app.py behind the display gate.
"""

import pytest

from gratinglab.gui.docs import find_theory_root
from gratinglab.gui.mathtext import (
    Segment,
    render_math_png,
    split_segments,
    strip_boxed,
)


class TestSplitSegments:
    def test_plain_text_is_a_single_segment(self):
        segments = split_segments("no math here at all")
        assert segments == (Segment(kind="text", content="no math here at all"),)

    def test_inline_math(self):
        segments = split_segments("the value $\\alpha$ is an angle")
        assert segments == (
            Segment(kind="text", content="the value "),
            Segment(kind="math", content="\\alpha", display=False),
            Segment(kind="text", content=" is an angle"),
        )

    def test_display_math(self):
        segments = split_segments("before\n\n$$x = 1$$\n\nafter")
        math = [s for s in segments if s.kind == "math"]
        assert len(math) == 1
        assert math[0].display is True
        assert math[0].content == "x = 1"

    def test_adjacent_math_spans_preserve_order_and_text_between(self):
        segments = split_segments("$a$ and $b$")
        assert [s.content for s in segments] == ["a", " and ", "b"]

    def test_escaped_dollar_is_literal_not_a_delimiter(self):
        segments = split_segments("costs \\$5 total")
        assert segments == (Segment(kind="text", content="costs $5 total"),)

    def test_escaped_dollar_inside_surrounding_text_of_real_math(self):
        segments = split_segments("\\$5 then $x$ then \\$10")
        assert [s.content for s in segments] == ["$5 then ", "x", " then $10"]

    def test_unterminated_dollar_is_literal_text_not_an_exception(self):
        segments = split_segments("oops $unterminated")
        assert segments == (Segment(kind="text", content="oops $unterminated"),)

    def test_unterminated_display_math_is_literal_text(self):
        segments = split_segments("oops $$unterminated")
        assert segments == (Segment(kind="text", content="oops $$unterminated"),)

    def test_empty_string(self):
        assert split_segments("") == ()

    def test_math_at_the_very_start_or_end(self):
        segments = split_segments("$x$ trailing")
        assert segments[0] == Segment(kind="math", content="x", display=False)
        segments = split_segments("leading $x$")
        assert segments[-1] == Segment(kind="math", content="x", display=False)

    def test_display_and_inline_do_not_confuse_each_other(self):
        """A display block right after inline math must not be misparsed as
        an empty inline span plus stray dollars."""
        segments = split_segments("inline $a$ then display $$b$$ done")
        kinds_and_display = [(s.kind, s.display) for s in segments if s.kind == "math"]
        assert kinds_and_display == [("math", False), ("math", True)]


class TestStripBoxed:
    def test_no_boxed_is_unchanged(self):
        assert strip_boxed("x = 1") == "x = 1"

    def test_simple_boxed(self):
        assert strip_boxed(r"\boxed{x = 1}") == "x = 1"

    def test_boxed_with_nested_braces(self):
        latex = r"\boxed{G_m = \int_0^1 e^{i\Phi_m(t)}\, dt}"
        assert strip_boxed(latex) == r"G_m = \int_0^1 e^{i\Phi_m(t)}\, dt"

    def test_two_boxed_expressions_side_by_side(self):
        """The exact shape used in scalar.md section 3."""
        latex = r"\boxed{\;G_m = \int_0^1 e^{i\Phi_m(t)}\,dt\;} \qquad \boxed{\;\mathscr{E}_m = |G_m|^2\;}"
        stripped = strip_boxed(latex)
        assert "\\boxed" not in stripped
        assert stripped == r"\;G_m = \int_0^1 e^{i\Phi_m(t)}\,dt\; \qquad \;\mathscr{E}_m = |G_m|^2\;"

    def test_unmatched_boxed_is_left_untouched(self):
        assert strip_boxed(r"\boxed{unterminated") == r"\boxed{unterminated"

    def test_text_outside_boxed_is_preserved(self):
        assert strip_boxed(r"before \boxed{x} after") == "before x after"


class TestRenderMathPng:
    def test_a_simple_expression_renders(self):
        png = render_math_png(r"\alpha + \beta")
        assert png is not None
        assert png[:8] == b"\x89PNG\r\n\x1a\n"

    def test_boxed_expressions_render_via_stripping(self):
        """The two headline equations in scalar.md must not silently fail."""
        png = render_math_png(r"\boxed{G_m = \int_0^1 e^{i\Phi_m(t)}\,dt}")
        assert png is not None

    def test_operatorname_text_and_qquad_all_render(self):
        """Confirmed by direct testing to be the macros scalar.md actually uses."""
        for latex in [
            r"\operatorname{sinc}(x)",
            r"\text{hello}",
            r"a \qquad b",
            r"\Phi_m(x) = k\,g(x)\,\sin\gamma\left[\cos\alpha + \cos\beta_m\right]",
        ]:
            assert render_math_png(latex) is not None, latex

    def test_an_unsupported_macro_returns_none_not_an_exception(self):
        assert render_math_png(r"\definitelynotarealmacro{x}") is None

    def test_render_is_deterministic(self):
        a = render_math_png(r"\alpha")
        b = render_math_png(r"\alpha")
        assert a == b


THEORY_ROOT = find_theory_root()
requires_theory_docs = pytest.mark.skipif(
    THEORY_ROOT is None, reason="docs/theory not found from this checkout"
)


@requires_theory_docs
class TestAgainstTheRealDocs:
    """Parsed out of the actual files, not hand-copied -- so an edit that
    introduces an unsupported macro fails here, not silently in the app."""

    def test_every_math_span_in_scalar_md_renders(self):
        text = (THEORY_ROOT / "scalar.md").read_text()
        math_segments = [s for s in split_segments(text) if s.kind == "math"]
        assert len(math_segments) > 10, "sanity check: did the file load correctly?"

        failures = [s.content for s in math_segments if render_math_png(s.content) is None]
        assert not failures, f"{len(failures)} span(s) failed to render: {failures}"

    def test_every_math_span_in_conventions_md_renders(self):
        conventions = THEORY_ROOT.parent / "conventions.md"
        if not conventions.is_file():
            pytest.skip("conventions.md not found")
        text = conventions.read_text()
        math_segments = [s for s in split_segments(text) if s.kind == "math"]

        failures = [s.content for s in math_segments if render_math_png(s.content) is None]
        assert not failures, f"{len(failures)} span(s) failed to render: {failures}"
