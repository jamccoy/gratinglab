"""
The widget layer.

The single most important property here is that the window reports what the
analysis computes. Everything else it does is presentation; if the number on
screen can drift from the number `analyze_single_file` returns, the window is
worse than useless because it looks authoritative.

Skipped entirely when PySide6 is absent, since the GUI is an optional extra.
"""
import contextlib
import io
import os

import pytest

pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"')

from afm_analysis.analyzer import analyze_single_file          # noqa: E402
from afm_analysis.config import PROJECT_ROOT                   # noqa: E402
from afm_analysis.gui.qt.main_window import MainWindow         # noqa: E402

SAMPLE = os.path.join(PROJECT_ROOT, 'data', 'ALD_master_1p5um_flatten.txt')


@pytest.fixture
def window(qtbot):
    win = MainWindow()
    qtbot.addWidget(win)
    yield win
    win.close()


def _run_and_wait(qtbot, window, timeout=60000):
    """
    Trigger an analysis and wait until the *window* has taken the result.

    Deliberately not `waitSignal(worker.finished)`. That returns when the signal
    is emitted, which across a thread boundary is not the same moment the
    window's slot has run - the assertion that follows then reads state the
    window has not written yet. It passed in isolation and failed in the full
    suite, which is the signature of exactly that race.
    """
    previous = window._result
    window.run_analysis()
    qtbot.waitUntil(lambda: window._result is not previous, timeout=timeout)


def test_qt_api_is_pinned_to_pyside6():
    """
    matplotlib's Qt shim picks a binding from QT_API when none is imported yet.
    This project shipped PyQt5 until recently, so a stray install is realistic
    and would otherwise produce a canvas that cannot parent into the window.
    """
    import afm_analysis.gui.qt  # noqa: F401
    assert os.environ.get("QT_API") == "PySide6"


class TestLoading:
    def test_loading_enables_analysis(self, qtbot, window):
        window.load(SAMPLE)
        assert window.run_btn.isEnabled()
        assert os.path.basename(SAMPLE) in window.file_label.text()

    def test_result_views_stay_disabled_until_there_is_a_result(self, qtbot, window):
        window.load(SAMPLE)
        assert not any(b.isEnabled() for b in window._result_buttons)


class TestTheNumbersOnScreen:
    """The property that makes the window trustworthy."""

    def test_reported_values_equal_a_direct_analysis(self, qtbot, window):
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)

        assert window._result is not None, "worker returned no result"
        with contextlib.redirect_stdout(io.StringIO()):
            direct = analyze_single_file(SAMPLE, show_plots=False,
                                         settings=window._settings)

        assert window._result['mean_angle'] == direct['mean_angle']
        assert window._result['n_grooves'] == direct['n_grooves']
        assert window._result['period_nm'] == direct['period_nm']

    def test_the_panel_shows_the_computed_mean(self, qtbot, window):
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)
        assert f"{window._result['mean_angle']:.3f}" in window.results_label.text()

    def test_changing_a_control_changes_the_result(self, qtbot, window):
        """Controls must actually reach the analysis, not just look adjustable"""
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)
        first = window._result['mean_angle']

        window.trim_spin.setValue(0.20)
        _run_and_wait(qtbot, window)
        second = window._result['mean_angle']

        assert first != second, (
            "facet trim 0.10 -> 0.20 should move the mean; if it does not, the "
            "control is not reaching the analysis")


class TestViews:
    def test_every_view_draws_without_error(self, qtbot, window):
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)
        for show in (window.show_raw_profile, window.show_2d,
                     window.show_row_groups, window.show_detection,
                     window.show_angles):
            show()

    def test_result_views_become_available(self, qtbot, window):
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)
        assert all(b.isEnabled() for b in window._result_buttons)


class TestWikiTab:
    """The Wiki tab, and that adding it disturbed nothing."""

    def test_both_tabs_exist(self, qtbot, window):
        assert [window.tabs.tabText(i) for i in range(window.tabs.count())] == \
            ["Analysis", "Wiki"]

    def test_every_page_renders(self, qtbot, window):
        for slug in window.wiki.slugs:
            assert window.wiki.show_page(slug), f"{slug} not selectable"
            assert window.wiki.current_slug == slug
            html = window.wiki.rendered_html()
            assert len(html) > 500, f"{slug} rendered almost nothing"

    def test_tables_survive_markdown_rendering(self, qtbot, window):
        """The ICC page's measured-values table is the reason tables matter"""
        window.wiki.show_page('icc-correction')
        assert '<table' in window.wiki.rendered_html()

    def test_switching_tabs_leaves_a_result_intact(self, qtbot, window):
        """The tab move must not disturb analysis state"""
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)
        before = window._result['mean_angle']

        window.tabs.setCurrentIndex(1)
        window.wiki.show_page('facet-fitting')
        window.tabs.setCurrentIndex(0)

        assert window._result['mean_angle'] == before
        assert all(b.isEnabled() for b in window._result_buttons)


class TestInvalidInput:
    def test_invalid_settings_report_without_running(self, qtbot, window):
        window.load(SAMPLE)
        window.trim_spin.setRange(0.0, 1.0)      # bypass the convenience cap
        window.trim_spin.setValue(0.9)
        window.run_analysis()                     # must not reach the worker
        assert window._result is None
        assert 'facet_trim' in window.results_label.text()
