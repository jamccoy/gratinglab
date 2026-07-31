"""Theory-page lookup.

Headless, no tkinter -- this module's whole purpose is to be the testable part
of "which page goes with which solver."
"""

from pathlib import Path

import pytest

from gratinglab.gui.docs import find_theory_root, theory_pages
from gratinglab.solvers.base import Capabilities, available_solvers, register


class TestFindTheoryRoot:
    def test_finds_the_real_repo_root(self):
        root = find_theory_root()
        assert root is not None
        assert root.name == "theory"
        assert (root / "scalar.md").is_file()

    def test_returns_none_when_nothing_is_found(self, tmp_path):
        """Simulates running from somewhere with no docs/theory above it --
        e.g. an installed wheel, where docs/ was never packaged."""
        isolated = tmp_path / "a" / "b" / "c"
        isolated.mkdir(parents=True)
        assert find_theory_root(isolated) is None

    def test_finds_a_root_several_levels_up(self, tmp_path):
        (tmp_path / "docs" / "theory").mkdir(parents=True)
        deep = tmp_path / "src" / "pkg" / "gui"
        deep.mkdir(parents=True)
        assert find_theory_root(deep) == tmp_path / "docs" / "theory"


class TestTheoryPages:
    def test_one_entry_per_registered_solver(self):
        pages = theory_pages()
        assert {p.name for p in pages} == set(available_solvers())

    def test_scalar_page_is_available_and_has_real_content(self):
        page = next(p for p in theory_pages() if p.name == "scalar")
        assert page.available
        assert page.path is not None and page.path.name == "scalar.md"
        # A heading known to exist in the real file, so this fails if the file
        # is ever emptied or renamed rather than passing on an empty string.
        assert "## 5. Energy is not conserved" in page.text

    def test_scalar_title_is_friendly_not_the_bare_key(self):
        page = next(p for p in theory_pages() if p.name == "scalar")
        assert page.title == "Scalar (Kirchhoff)"

    def test_scalar_reports_its_rigorous_flag(self):
        """Scalar theory is approximate; the viewer decides on the banner from this."""
        page = next(p for p in theory_pages() if p.name == "scalar")
        assert page.rigorous is False

    def test_a_solver_with_no_page_yet_is_still_listed(self):
        """The whole point: report missing docs, do not hide the solver."""

        class Stub:
            capabilities = Capabilities(name="stub_no_docs", rigorous=True)

            @staticmethod
            def solve(problem, illumination, wavelengths, **options):
                raise NotImplementedError

        register(Stub())
        try:
            page = next(p for p in theory_pages() if p.name == "stub_no_docs")
            assert not page.available
            assert "No theory page has been written yet" in page.text
            assert "stub_no_docs" in page.text
            assert page.rigorous is True
        finally:
            from gratinglab.solvers import base

            del base._REGISTRY["stub_no_docs"]

    def test_missing_theory_root_gives_every_solver_an_explanation(self, monkeypatch):
        """When docs/ is not bundled (e.g. an installed wheel), every entry
        must still carry a usable message rather than raising."""
        monkeypatch.setattr(
            "gratinglab.gui.docs.find_theory_root", lambda *a, **k: None
        )
        for page in theory_pages():
            assert not page.available
            assert page.path is None
            assert "not bundled" in page.text

    def test_text_is_never_empty(self):
        """The viewer always has something to show."""
        for page in theory_pages():
            assert page.text.strip()
