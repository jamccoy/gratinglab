"""
The wiki content, checked without a window.

These pages are the project's explanation of its own statistics, so the failure
worth guarding against is not a crash - it is a page that quietly went missing,
emptied, or stopped saying what the code does.
"""
import os
import re
import sys

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from afm_analysis.wiki import PAGE_ORDER, WIKI_DIR, page, pages  # noqa: E402


def test_every_expected_page_exists():
    found = {p.slug for p in pages()}
    missing = set(PAGE_ORDER) - found
    assert not missing, f"missing wiki pages: {sorted(missing)}"


def test_pages_come_back_in_reading_order():
    """
    Ordering is explicit, not alphabetical.

    A reader should meet the overview before the correction it motivates;
    alphabetical order would open on facet-fitting.
    """
    slugs = [p.slug for p in pages()]
    expected = [s for s in PAGE_ORDER if s in slugs]
    assert slugs[:len(expected)] == expected


def test_no_page_is_empty():
    for wiki_page in pages():
        assert not wiki_page.is_empty, f"{wiki_page.slug} has a title but no body"
        assert len(wiki_page.markdown) > 500, (
            f"{wiki_page.slug} is suspiciously short ({len(wiki_page.markdown)} chars)")


def test_every_page_has_a_title_from_its_heading():
    for wiki_page in pages():
        assert wiki_page.title
        assert not wiki_page.title.startswith('#')
        first_heading = re.search(r'^#\s+(.+)$', wiki_page.markdown, re.MULTILINE)
        assert first_heading, f"{wiki_page.slug} has no level-1 heading"
        assert wiki_page.title == first_heading.group(1).strip()


def test_slugs_are_unique():
    slugs = [p.slug for p in pages()]
    assert len(slugs) == len(set(slugs))


def test_lookup_by_slug():
    assert page('icc-correction') is not None
    assert page('icc-correction').slug == 'icc-correction'
    assert page('no-such-page') is None


def test_pages_are_read_fresh_rather_than_cached():
    """
    Editing a page and reopening the tab should show the edit.

    A cache populated at import would make that a confusing thing to debug, so
    pages() reads from disk each call.
    """
    first = pages()
    second = pages()
    assert first == second          # same content
    assert first is not second      # but not the same object


class TestContentStaysTrue:
    """
    A handful of specific claims, pinned.

    These pages state numbers the code produces. If the code changes and the
    prose does not, the wiki becomes the most authoritative-looking wrong
    document in the project - worse than no wiki. Each assertion here names a
    fact that a real change would invalidate.
    """

    def _markdown(self, slug):
        found = page(slug)
        assert found is not None, f"{slug} missing"
        return found.markdown

    def test_icc_page_states_the_measured_range(self):
        text = self._markdown('icc-correction')
        assert '0.097' in text and '0.429' in text, "measured ICC range not stated"
        assert '0.244' in text, "median ICC not stated"

    def test_icc_page_gives_the_design_effect_formula(self):
        text = self._markdown('icc-correction')
        assert 'N_eff = N / (1 + (m − 1) · ICC)' in text

    def test_icc_page_records_that_no_conclusion_changed(self):
        """The honest headline. If this ever stops being true it must be edited."""
        assert 'No conclusion' in self._markdown('icc-correction')

    def test_facet_page_states_the_trim_limit(self):
        text = self._markdown('facet-fitting')
        assert '0.286' in text, "the hard trim limit is not stated"
        assert '3.5' in text, "the trim multiplier is not explained"

    def test_edge_page_uses_the_real_failure_as_the_example(self):
        text = self._markdown('row-groups')
        assert '2.40' in text, "the worked example angles are missing"
        assert '25.58' in text, "the post-fix minimum is not stated"

    def test_outputs_page_covers_both_sem_columns(self):
        text = self._markdown('reading-outputs')
        assert 'SEM_deg' in text and 'SEM_corrected_deg' in text
        assert 'Row_group' in text, "the column that makes ICC computable is unlisted"


def test_markdown_files_are_package_data():
    """
    The pages must live inside the package, not at the repo root.

    Content outside src/ ships in a checkout but not in a wheel, which would
    leave the Wiki tab empty for anyone who pip installed. pyproject.toml
    declares afm_analysis.wiki = ["*.md"] to match.
    """
    assert WIKI_DIR.name == 'wiki'
    assert WIKI_DIR.parent.name == 'afm_analysis'
    assert list(WIKI_DIR.glob('*.md')), "no markdown found beside the module"


if __name__ == '__main__':
    failures = 0
    targets = list(globals().items())
    for name, obj in sorted(targets):
        if name.startswith('test_'):
            try:
                obj()
                print(f"PASS  {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL  {name}: {exc}")
        elif name.startswith('Test'):
            instance = obj()
            for method in sorted(dir(instance)):
                if method.startswith('test_'):
                    try:
                        getattr(instance, method)()
                        print(f"PASS  {name}.{method}")
                    except AssertionError as exc:
                        failures += 1
                        print(f"FAIL  {name}.{method}: {exc}")
    print(f"\n{'all tests passed' if not failures else f'{failures} failure(s)'}")
    sys.exit(1 if failures else 0)
