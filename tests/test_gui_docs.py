"""Theory-page lookup.

Headless, no tkinter -- this module's whole purpose is to be the testable part
of "which page goes with which solver."
"""

from pathlib import Path

import pytest

from gratinglab.gui.docs import find_theory_root, general_pages, theory_pages
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
            assert "No page has been written yet" in page.text
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


class TestGeneralPages:
    @staticmethod
    def _page(name):
        """By name rather than by index -- the tuple has grown once already."""
        return next(p for p in general_pages() if p.name == name)

    def test_conventions_page_is_available_with_real_content(self):
        page = self._page("conventions")
        assert page.title == "Grating Geometry & Conventions"
        assert page.available
        assert page.path is not None and page.path.name == "conventions.md"

    def test_conventions_page_contains_the_grating_equation(self):
        """The whole point of adding this page."""
        page = self._page("conventions")
        assert "generalized (conical) grating equation" in page.text
        assert "sin\\gamma" in page.text or "sin\\,\\gamma" in page.text

    def test_metrology_page_is_available_with_real_content(self):
        page = self._page("metrology")
        assert page.title == "Groove Metrology"
        assert page.available
        assert page.path is not None and page.path.name == "metrology.md"

    def test_metrology_page_states_that_roughness_is_not_measured(self):
        """The assumption most likely to be silently assumed away."""
        page = self._page("metrology")
        assert "RMS surface roughness" in page.text
        assert "not yet computed" in page.text

    def test_general_pages_are_never_banner_flagged(self):
        """These are geometry references, not a solver's approximation --
        the approximate-method banner must never apply to them."""
        for page in general_pages():
            assert page.rigorous is True

    def test_does_not_collide_with_solver_names(self):
        """general_pages() and theory_pages() must be safely concatenable in
        the Help menu without a name clash."""
        general_names = {p.name for p in general_pages()}
        solver_names = {p.name for p in theory_pages()}
        assert general_names.isdisjoint(solver_names)

    def test_missing_theory_root_gives_a_usable_explanation(self, monkeypatch):
        monkeypatch.setattr(
            "gratinglab.gui.docs.find_theory_root", lambda *a, **k: None
        )
        page = general_pages()[0]
        assert not page.available
        assert "not bundled" in page.text
