"""gratinglab's Help pages as rich text.

The markdown renderer is endstation's and is tested there. These run it on the
real pages, so an edit to ``scalar.md`` that breaks rendering fails here
rather than being found by eye, and they test the two things gratinglab adds:
the approximate-method banner and menu labels.

Headless: no Qt, no display.
"""

import ast
import inspect
import re
import struct
import zlib

import pytest
from endstation.docs.markdown import render_markdown
from endstation.docs.mathtext import split_segments

from gratinglab.gui import richtext
from gratinglab.gui.docs import TheoryPage, general_pages, theory_pages
from gratinglab.gui.richtext import menu_label, to_html


def fake_render(latex, *, dpi, color, **_):
    """A 4x3 PNG, so a page renders without matplotlib."""
    def chunk(kind, payload):
        body = kind + payload
        return struct.pack(">I", len(payload)) + body + struct.pack(">I", zlib.crc32(body))

    ihdr = struct.pack(">IIBBBBB", 4, 3, 8, 0, 0, 0, 0)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IEND", b"")


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
    """The pages are mostly tables, with equations in the cells."""

    def test_the_tension_table_is_a_table(self, scalar_html):
        assert "<table" in scalar_html and scalar_html.count("<tr>") > 10
        assert "<th>" in scalar_html

    def test_math_renders_inside_table_cells(self, scalar_html):
        """The tension table's column headings are equations."""
        cells = re.findall(r"<t[hd]>.*?</t[hd]>", scalar_html, re.S)
        assert len([c for c in cells if "<img " in c]) > 3

    def test_lists_rules_quotes_and_the_title(self, scalar_html, conventions_html):
        assert "<ul>" in conventions_html and "<hr>" in conventions_html
        assert "<blockquote>" in scalar_html
        assert "<h1>Conventions</h1>" in conventions_html

    def test_no_raw_delimiters_latex_or_placeholders_reach_the_output(
        self, scalar_html, conventions_html
    ):
        for html in (scalar_html, conventions_html):
            assert "$" not in html and "\x00" not in html
        assert "\\Phi_m" not in scalar_html and "\\alpha" not in scalar_html

    def test_an_ampersand_in_prose_is_escaped(self, conventions_html):
        """conventions.md cites 'Heilmann ... & McEntaffer'."""
        assert "&amp;" in conventions_html

    def test_a_code_span_inside_a_link_survives(self, scalar_html):
        """scalar.md's first line is [`scalar.py`](../../src/...)."""
        assert "<code>solvers/scalar.py</code></a>" in scalar_html

    def test_every_equation_becomes_a_referenced_image(self, scalar_page):
        rendered = to_html(scalar_page, render=fake_render)
        expected = sum(1 for s in split_segments(scalar_page.text) if s.kind == "math")
        assert len(rendered.images) == expected > 50
        for image in rendered.images:
            assert f'src="{image.key}"' in rendered.html

    def test_every_equation_actually_rasterises(self, scalar_page):
        """The real renderer: an unsupported macro added later degrades
        visibly, and this is what notices."""
        assert render_markdown(scalar_page.text).unrendered == ()

    def test_the_metrology_steps_are_a_numbered_list(self):
        page = next(p for p in general_pages() if p.name == "metrology")
        html = to_html(page, render=fake_render).html
        assert "<ol><li><b>Endpoint flattening.</b>" in html
        assert "<p>1. " not in html


class TestBanner:
    def test_a_non_rigorous_solver_is_announced(self, scalar_page):
        assert not scalar_page.rigorous  # the case this test needs
        assert "Approximate method" in to_html(scalar_page, render=fake_render).html

    def test_a_rigorous_one_is_not(self):
        page = TheoryPage(name="integral", title="Integral", available=True, path=None,
                          text="## Theory\n\nExact.", rigorous=True)
        assert richtext.banner(page) == ""
        assert "Approximate method" not in to_html(page, render=fake_render).html

    def test_a_page_that_does_not_exist_yet_still_renders(self):
        page = TheoryPage(name="future", title="Future Method", available=False,
                          path=None, rigorous=True,
                          text="No theory page has been written yet for 'future'.")
        assert "No theory page has been written yet" in to_html(page).html


class TestMenuLabel:
    def test_an_ampersand_is_escaped_for_qt(self):
        assert menu_label("Grating Geometry & Conventions") == "Grating Geometry && Conventions"

    def test_the_real_page_title_needs_it(self):
        assert any("&" in p.title for p in general_pages()), "escaping test went vacuous"

    def test_a_suffix_is_appended_before_escaping(self):
        assert menu_label("A & B", " Theory") == "A && B Theory"

    def test_a_title_without_an_ampersand_is_untouched(self):
        assert menu_label("Scalar (Kirchhoff)", " Theory") == "Scalar (Kirchhoff) Theory"


def test_richtext_imports_no_toolkit():
    imported = set()
    for node in ast.walk(ast.parse(inspect.getsource(richtext))):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert not {"tkinter", "PySide6", "PyQt6"} & imported
