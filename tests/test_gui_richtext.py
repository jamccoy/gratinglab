r"""Theory pages as rich text.

Headless -- no Qt, no display. The pipeline down to an HTML string plus PNG
bytes is entirely pure, which is the point: whether the conventions table
renders *as a table* is checkable without opening a window.

Tests that matter run against the **real** doc files rather than hand-written
markdown, so an edit to `scalar.md` that breaks rendering fails here rather
than being discovered by eye.
"""

import re

import pytest

from gratinglab.gui import richtext
from gratinglab.gui.docs import TheoryPage, general_pages, theory_pages
from gratinglab.gui.mathtext import split_segments
from gratinglab.gui.richtext import (
    DISPLAY_DPI,
    INLINE_DPI,
    RichPage,
    math_dpi,
    menu_label,
    png_size,
    render_markdown,
    to_html,
)


def fake_render(latex, *, dpi, color):
    """A 4x3 PNG, so rendering can be exercised without matplotlib."""
    import struct
    import zlib

    def chunk(kind, payload):
        body = kind + payload
        return (
            struct.pack(">I", len(payload))
            + body
            + struct.pack(">I", zlib.crc32(body))
        )

    ihdr = struct.pack(">IIBBBBB", 4, 3, 8, 0, 0, 0, 0)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IEND", b"")


def markup(source, **kwargs):
    """Render markdown without rasterising anything."""
    kwargs.setdefault("render", fake_render)
    return render_markdown(source, **kwargs)


@pytest.fixture(scope="module")
def scalar_page():
    return next(p for p in theory_pages() if p.name == "scalar")


@pytest.fixture(scope="module")
def conventions_page():
    return next(p for p in general_pages() if p.name == "conventions")


@pytest.fixture(scope="module")
def scalar_html(scalar_page):
    return to_html(scalar_page, render=fake_render).html


@pytest.fixture(scope="module")
def conventions_html(conventions_page):
    return to_html(conventions_page, render=fake_render).html


class TestTheRealPages:
    """What this module was written for: the docs are mostly tables now."""

    def test_the_tension_table_is_a_table(self, scalar_html):
        assert "<table" in scalar_html
        assert scalar_html.count("<tr>") > 10

    def test_table_headers_are_headers(self, scalar_html):
        assert "<th>" in scalar_html

    def test_math_renders_inside_table_cells(self, scalar_html):
        """The specific thing that rules out rendering markdown wholesale and
        patching images in afterwards -- the tension table's own column
        headings are equations."""
        cells_with_math = [
            cell
            for cell in re.findall(r"<t[hd]>.*?</t[hd]>", scalar_html, re.S)
            if "<img " in cell
        ]
        assert len(cells_with_math) > 3

    def test_bullet_lists_are_lists(self, conventions_html):
        assert "<ul>" in conventions_html and "<li>" in conventions_html

    def test_the_blockquote_is_a_blockquote(self, scalar_html):
        """scalar.md section 5's warning that referencing the blaze direction
        is circular. It reads as an aside and must look like one."""
        assert "<blockquote>" in scalar_html

    def test_horizontal_rules_survive(self, conventions_html):
        assert "<hr>" in conventions_html

    def test_the_level_one_title_is_a_heading(self, conventions_html):
        """The Tk viewer matched `#{2,6}` and so showed a literal `#` on the
        first line of every page. Both real pages open with `# Title`."""
        assert "<h1>Conventions</h1>" in conventions_html

    def test_no_raw_dollar_delimiters_reach_the_output(
        self, scalar_html, conventions_html
    ):
        assert "$" not in scalar_html
        assert "$" not in conventions_html

    def test_no_raw_latex_reaches_the_output(self, scalar_html):
        assert "\\Phi_m" not in scalar_html
        assert "\\alpha" not in scalar_html

    def test_every_equation_becomes_an_image(self, scalar_page):
        rendered = to_html(scalar_page, render=fake_render)
        expected = sum(
            1 for s in split_segments(scalar_page.text) if s.kind == "math"
        )
        assert len(rendered.images) == expected
        assert expected > 50

    def test_image_keys_are_unique_and_referenced(self, scalar_page):
        rendered = to_html(scalar_page, render=fake_render)
        keys = [image.key for image in rendered.images]
        assert len(set(keys)) == len(keys)
        for key in keys:
            assert f'src="{key}"' in rendered.html

    def test_every_equation_actually_rasterises(self, scalar_page):
        """The real renderer, not the fake one. An unsupported macro added in
        a future edit degrades visibly, and this is what notices.

        `render_markdown`, not `to_html`, so the approximate-method banner is
        not in the string being searched.
        """
        rendered = render_markdown(scalar_page.text)
        assert richtext.UNRENDERED_COLOR not in rendered.html

    def test_a_code_span_inside_a_link_survives(self, scalar_html):
        """scalar.md's first line is [`scalar.py`](../../src/...). The
        placeholder that protects the code span from emphasis nests inside the
        link's, and `re.sub` does not rescan its own substitutions -- one pass
        left a raw null byte where the filename should be.
        """
        assert "<code>solvers/scalar.py</code></a>" in scalar_html
        assert "\x00" not in scalar_html

    def test_no_placeholder_bytes_escape_into_any_page(
        self, scalar_html, conventions_html
    ):
        assert "\x00" not in scalar_html
        assert "\x00" not in conventions_html


class TestEscaping:
    def test_an_ampersand_in_prose_is_escaped(self, conventions_html):
        """conventions.md cites 'Heilmann ... & McEntaffer'."""
        assert "&amp;" in conventions_html

    def test_math_is_split_out_before_anything_is_escaped(self):
        r"""Four math spans in the real docs contain `<` or `>`. Escaping
        first would hand mathtext `t &lt; w`, which it renders literally.

        Asserted on what reaches the renderer, since with images the
        corruption would otherwise be invisible in the HTML.
        """
        seen = []

        def spy(latex, *, dpi, color):
            seen.append(latex)
            return fake_render(latex, dpi=dpi, color=color)

        markup("Valid for $t < w$ and $k > 0$.", render=spy)
        assert seen == ["t < w", "k > 0"]

    def test_a_bracket_in_prose_cannot_become_a_tag(self):
        html = markup("Compare <script>alert(1)</script> here.").html
        assert "<script>" not in html
        assert "&lt;script&gt;" in html

    def test_escaping_does_not_break_bold(self):
        assert "<b>bold</b>" in markup("**bold** & plain").html


class TestInlineMarkup:
    def test_bold(self):
        assert "<b>x</b>" in markup("**x**").html

    def test_italic(self):
        assert "<i>x</i>" in markup("an *x* here").html

    def test_code_stays_literal(self):
        """Markup characters inside a code span must not be interpreted."""
        html = markup("use `a**b**c` now").html
        assert "<code>a**b**c</code>" in html

    def test_links_become_anchors(self):
        html = markup("see [the docs](conventions.md) for more").html
        assert '<a href="conventions.md">the docs</a>' in html

    def test_a_link_url_is_not_scanned_for_emphasis(self):
        html = markup("[x](a_b_c*d*e.md)").html
        assert "a_b_c*d*e.md" in html

    def test_bold_inside_a_bullet(self):
        assert "<li><b>A deficit is ordinary.</b>" in markup(
            "- **A deficit is ordinary.** Power goes elsewhere."
        ).html

    def test_a_lone_asterisk_is_not_emphasis(self):
        assert "<i>" not in markup("2 * 3 = 6").html


class TestBlockParsing:
    def test_soft_wrapped_paragraphs_are_joined(self):
        """The docs hard-wrap at about 80 columns; rendering each source line
        as its own paragraph would double-space the whole page."""
        html = markup("one two\nthree four\n\nnext").html
        assert "<p>one two three four</p>" in html
        assert html.count("<p>") == 2

    def test_a_bullet_continuation_line_joins_its_item(self):
        html = markup("- first part\n  second part\n- other").html
        assert "<li>first part second part</li>" in html
        assert html.count("<li>") == 2

    def test_a_blockquote_joins_its_lines(self):
        html = markup("> one\n> two\n\nafter").html
        assert "<blockquote><p>one two</p></blockquote>" in html

    def test_a_table_delimiter_row_is_not_a_data_row(self):
        html = markup("| a | b |\n|---|---|\n| 1 | 2 |").html
        assert html.count("<tr>") == 2
        assert "---" not in html

    def test_every_heading_level(self):
        html = markup("\n\n".join(f"{'#' * n} h{n}" for n in range(1, 7))).html
        for n in range(1, 7):
            assert f"<h{n}>h{n}</h{n}>" in html

    def test_a_rule_is_a_rule_not_a_heading_underline(self):
        """These docs use `---` as a horizontal rule throughout, never as a
        setext underline."""
        assert markup("text\n\n---\n\nmore").html.count("<hr>") == 1

    def test_a_fenced_code_block_stays_literal(self):
        """No page uses one today, but rcwa.md will. Interpreting markdown
        inside code would mangle it."""
        html = markup("```\na = **b**\n```").html
        assert "<pre" in html
        assert "**b**" in html
        assert "<b>" not in html

    def test_an_unclosed_fence_does_not_swallow_the_parser(self):
        assert "<pre" in markup("```\nx = 1").html

    def test_an_empty_page_is_empty_not_an_error(self):
        assert markup("").html == ""

    def test_blank_lines_produce_no_empty_paragraphs(self):
        assert "<p></p>" not in markup("a\n\n\n\nb").html


class TestDegradation:
    def test_an_unparseable_span_shows_its_source_in_colour(self):
        """Same signal the Tk viewer's `unrendered` tag gave: a rendering gap
        is visible and diagnosable, never silently swallowed."""
        html = render_markdown("before $\\notamacro{x}$ after", render=lambda *a, **k: None).html
        assert richtext.UNRENDERED_COLOR in html
        assert "\\notamacro{x}" in html

    def test_a_failure_does_not_lose_the_rest_of_the_page(self):
        page = render_markdown(
            "## Heading\n\ntext $bad$ more\n\n| a |\n|---|\n| b |",
            render=lambda *a, **k: None,
        )
        assert "<h2>Heading</h2>" in page.html
        assert "<table" in page.html

    def test_a_failed_span_registers_no_image(self):
        page = render_markdown("$x$", render=lambda *a, **k: None)
        assert page.images == ()

    def test_the_failed_source_is_escaped(self):
        """The fallback shows raw LaTeX, which is exactly the text most likely
        to contain a bracket."""
        html = render_markdown("$a < b$", render=lambda *a, **k: None).html
        assert "&lt;" in html


class TestHiDpi:
    def test_display_equations_render_larger_than_inline_ones(self):
        seen = {}

        def spy(latex, *, dpi, color):
            seen[latex] = dpi
            return fake_render(latex, dpi=dpi, color=color)

        markup("inline $a$ and\n\n$$b$$\n", render=spy)
        assert seen["a"] == INLINE_DPI
        assert seen["b"] == DISPLAY_DPI
        assert seen["b"] > seen["a"]

    def test_a_retina_ratio_rasterises_at_double_resolution(self):
        seen = []

        def spy(latex, *, dpi, color):
            seen.append(dpi)
            return fake_render(latex, dpi=dpi, color=color)

        markup("$a$", device_pixel_ratio=2.0, render=spy)
        assert seen == [2 * INLINE_DPI]

    def test_but_is_displayed_at_logical_size(self):
        """Otherwise a 2x image would render at double size on a Retina panel
        instead of crisply at the right one."""
        one = markup("$a$", device_pixel_ratio=1.0).images[0]
        two = markup("$a$", device_pixel_ratio=2.0).images[0]
        # fake_render always returns a 4x3 PNG, so the logical box halves.
        assert (one.width, one.height) == (4, 3)
        assert (two.width, two.height) == (2, 2)  # 3/2 rounds to 2

    def test_a_logical_size_never_rounds_to_zero(self):
        assert markup("$a$", device_pixel_ratio=8.0).images[0].width >= 1

    def test_the_color_is_passed_through(self):
        """So a dark-mode viewer can render light-on-dark math."""
        seen = []
        markup(
            "$a$",
            color="#eeeeee",
            render=lambda latex, *, dpi, color: seen.append(color)
            or fake_render(latex, dpi=dpi, color=color),
        )
        assert seen == ["#eeeeee"]


class TestPngSize:
    def test_reads_the_ihdr_dimensions(self):
        assert png_size(fake_render("x", dpi=1, color="black")) == (4, 3)

    def test_matches_a_real_matplotlib_png(self):
        from gratinglab.gui.mathtext import render_math_png

        png = render_math_png("x", dpi=100)
        assert png is not None
        width, height = png_size(png)
        assert width > 1 and height > 1

    def test_truncated_data_returns_none_rather_than_raising(self):
        assert png_size(b"\x89PNG") is None

    def test_data_that_is_not_a_png_returns_none(self):
        assert png_size(b"x" * 64) is None


class TestMathDpi:
    def test_ten_point_text_wants_seventy_two_dpi(self):
        """mathtext renders at rcParams['font.size'], 10 pt by default."""
        assert math_dpi(10.0) == 72

    def test_larger_text_wants_proportionally_more(self):
        assert math_dpi(20.0) == 2 * math_dpi(10.0)

    def test_display_equations_get_a_boost(self):
        assert math_dpi(10.0, display=True) > math_dpi(10.0)

    def test_never_returns_zero(self):
        assert math_dpi(0.0) >= 1


class TestMenuLabel:
    def test_an_ampersand_is_escaped_for_qt(self):
        """Qt reads `&` as 'underline the next character', so the real page
        title would appear as 'Grating Geometry Conventions'."""
        assert menu_label("Grating Geometry & Conventions") == (
            "Grating Geometry && Conventions"
        )

    def test_the_real_page_title_needs_it(self):
        titles = [p.title for p in general_pages()]
        assert any("&" in title for title in titles), "escaping test went vacuous"

    def test_a_suffix_is_appended_before_escaping(self):
        assert menu_label("A & B", " Theory") == "A && B Theory"

    def test_a_title_without_an_ampersand_is_untouched(self):
        assert menu_label("Scalar (Kirchhoff)", " Theory") == (
            "Scalar (Kirchhoff) Theory"
        )


class TestBanner:
    def test_a_non_rigorous_solver_is_announced(self, scalar_page):
        assert not scalar_page.rigorous  # the case this test needs
        assert "Approximate method" in to_html(scalar_page, render=fake_render).html

    def test_a_rigorous_one_is_not(self):
        page = TheoryPage(
            name="integral", title="Integral", available=True, path=None,
            text="## Theory\n\nExact.", rigorous=True,
        )
        assert "Approximate method" not in to_html(page, render=fake_render).html

    def test_a_page_that_does_not_exist_yet_still_renders(self):
        page = TheoryPage(
            name="future", title="Future Method", available=False, path=None,
            text="No theory page has been written yet for 'future'.",
            rigorous=True,
        )
        assert "No theory page has been written yet" in to_html(page).html


class TestPurity:
    def test_imports_no_toolkit(self):
        import ast
        import inspect

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(richtext))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6"} & imported

    def test_returns_a_richpage_not_a_widget(self):
        assert isinstance(markup("x"), RichPage)


class TestCaching:
    def test_the_same_equation_is_rasterised_once(self):
        """A page carries dozens of equations and the viewer reopens often;
        without this, every open re-renders the lot."""
        from gratinglab.gui.richtext import _cached_render

        _cached_render.cache_clear()
        render_markdown("$x$ and $x$ and $y$")
        info = _cached_render.cache_info()
        assert info.misses == 2
        assert info.hits == 1
