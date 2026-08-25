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

from gratinglab.metrology.analyzer import analyze_single_file          # noqa: E402
from gratinglab.metrology.gui.qt.main_window import MainWindow         # noqa: E402
from scans import SYNTHETIC, real_scan                                 # noqa: E402

SAMPLE = str(SYNTHETIC)


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
    import gratinglab.metrology.gui.qt  # noqa: F401
    assert os.environ.get("QT_API") == "PySide6"


class TestLoading:
    def test_loading_enables_analysis(self, qtbot, window):
        window.load(SAMPLE)
        assert window.run_btn.isEnabled()
        assert os.path.basename(SAMPLE) in window.importer.file_label.text()

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


class TestImportTab:
    """Import owns loading; Analysis consumes what it produces."""

    def test_loading_in_import_enables_analysis(self, qtbot, window):
        window.load(SAMPLE)
        assert window.run_btn.isEnabled()
        assert window._data is not None
        assert window.importer.filename == SAMPLE

    def test_image_flattening_choice_reaches_the_form(self, qtbot, window):
        window.load(SAMPLE)
        window.importer.image_method_combo.setCurrentText('none')
        assert window.form_state().image_flatten_method == 'none'

    def test_profile_flattening_choice_changes_the_result(self, qtbot, window):
        """
        The knob that matters. Roughly 0.5 degrees between methods, so a control
        that leaves the answer alone is not wired up.

        Needs a real scan. The methods differ in how they estimate a background
        to subtract, and on the synthetic fixture -- whose background is one
        exact plane -- they agree to 0.04 deg, which cannot distinguish "wired
        up" from "ignored". That agreement is a fact about ideal data, not
        evidence about the control.
        """
        window.load(str(real_scan('ALD_master_1p5um_flatten.txt')))
        _run_and_wait(qtbot, window)
        before = window._result['mean_angle']

        window.importer.profile_method_combo.setCurrentText('groove_peaks')
        _run_and_wait(qtbot, window)
        after = window._result['mean_angle']

        assert abs(after - before) > 0.05, (
            f"profile flattening did not reach the analysis: {before} -> {after}")

    def test_affine_image_flattening_leaves_the_answer_alone(self, qtbot, window):
        """
        Free by construction - the reason align_rows could be the default.

        Equal to within floating point, not bit-for-bit: subtracting row medians
        perturbs the last bits, and the averaging that follows is not
        associative. Measured at 5e-15 degrees, against a sigma near 2.
        """
        window.load(SAMPLE)
        _run_and_wait(qtbot, window)
        with_align = window._result['mean_angle']

        window.importer.image_method_combo.setCurrentText('none')
        _run_and_wait(qtbot, window)
        assert abs(window._result['mean_angle'] - with_align) < 1e-9


class TestBoundaryTab:
    """PCGrate export, driven from the window."""

    def test_tab_order(self, qtbot, window):
        assert [window.tabs.tabText(i) for i in range(window.tabs.count())] == \
            ["Import", "Analysis", "Boundary", "Wiki"]

    def test_loading_populates_the_boundary_panel(self, qtbot, window):
        assert window.boundary.profile is None
        window.load(SAMPLE)
        assert window.boundary.profile is not None
        assert window.boundary.export_btn.isEnabled()

    def test_metrics_are_shown(self, qtbot, window):
        window.load(SAMPLE)
        text = window.boundary.metrics_label.text()
        for field in ("Period", "Grooves averaged", "Groove depth",
                      "Max sidewall"):
            assert field in text, f"{field} missing from the metrics panel"

    def test_point_count_control_reaches_the_profile(self, qtbot, window):
        window.load(SAMPLE)
        window.boundary.n_points_spin.setValue(500)
        assert len(window.boundary.profile.x_norm) == 500

    def test_export_writes_a_valid_ggp(self, qtbot, window, tmp_path, monkeypatch):
        """The header must be the uncommented two-line form PCGrate accepts"""
        from PySide6.QtWidgets import QFileDialog
        window.load(SAMPLE)
        target = str(tmp_path / "out.ggp")
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            staticmethod(lambda *a, **k: (target, "")))
        monkeypatch.setattr("gratinglab.metrology.gui.qt.boundary_view.QMessageBox",
                            type("Stub", (), {"information": staticmethod(lambda *a: None),
                                              "warning": staticmethod(lambda *a: None)}))
        window.boundary.export()

        with open(target) as handle:
            lines = handle.read().splitlines()
        assert lines[0] == "3 0 - Polygonal type"
        assert lines[1] == "Period: 1 PSC: 1"
        assert os.path.exists(str(tmp_path / "out_metrics.txt"))

    def test_efficiency_button_tracks_the_export_button(self, qtbot, window):
        """Both need one thing -- a profile -- so they must not drift apart."""
        assert not window.boundary.export_btn.isEnabled()
        assert not window.boundary.efficiency_btn.isEnabled()
        window.load(SAMPLE)
        assert window.boundary.export_btn.isEnabled()
        assert window.boundary.efficiency_btn.isEnabled()

    def test_panel_matches_a_direct_build(self, qtbot, window):
        """The window must show what the pipeline computes"""
        import numpy as np
        from gratinglab.metrology.boundary import build_boundary_profile
        window.load(SAMPLE)
        direct = build_boundary_profile(
            window._data, window._scan_size,
            window._defaults.with_(**window.boundary.settings_overrides()))
        assert np.array_equal(direct.y_norm, window.boundary.profile.y_norm)


class TestWikiTab:
    """The Wiki tab, and that adding it disturbed nothing."""

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

        window.tabs.setCurrentIndex(2)          # Wiki
        window.wiki.show_page('facet-fitting')
        window.tabs.setCurrentIndex(1)          # back to Analysis

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


class TestEfficiencyDialog:
    """The handoff that makes the two halves one package."""

    def _dialog(self, qtbot, window):
        from gratinglab.metrology.gui.qt.efficiency_dialog import EfficiencyDialog
        window.load(SAMPLE)
        dialog = EfficiencyDialog(window.boundary.profile)
        qtbot.addWidget(dialog)
        return dialog

    def test_it_uses_the_measured_period_without_being_told(self, qtbot, window):
        """The number a .ggp cannot carry, and would otherwise be retyped."""
        dialog = self._dialog(qtbot, window)
        assert dialog._to_problem().period == window.boundary.profile.period_nm

    def test_computing_produces_a_finite_scan(self, qtbot, window):
        import numpy as np
        dialog = self._dialog(qtbot, window)
        dialog.count_spin.setValue(12)
        dialog.compute()
        assert dialog._scan is not None
        assert np.all(np.isfinite(dialog._scan.efficiency))

    def test_the_facet_angle_reaches_the_problem_only_when_enabled(self, qtbot, window):
        dialog = self._dialog(qtbot, window)
        assert dialog._to_problem().profile.blaze_angle is None
        dialog.blaze_check.setChecked(True)
        dialog.blaze_spin.setValue(29.5)
        assert dialog._to_problem().profile.blaze_angle == pytest.approx(29.5)

    def test_it_says_the_result_is_not_defensible(self, qtbot, window):
        """No convergence check runs here, and the panel must not imply one."""
        dialog = self._dialog(qtbot, window)
        dialog.count_spin.setValue(8)
        dialog.compute()
        assert "not convergence-checked" in dialog.status.text()

    def test_it_reports_the_resolving_power_of_the_measured_grooves(
            self, qtbot, window):
        """R = |m|*N from the groove count the scan measured, no retyping."""
        dialog = self._dialog(qtbot, window)
        dialog.count_spin.setValue(8)
        dialog.compute()
        n = window.boundary.profile.n_used
        assert f"Resolving power from the {n} averaged grooves" in dialog.status.text()
        assert f"m = +1: R = {n:,.0f}" in dialog.status.text()

    def test_a_backwards_wavelength_range_is_refused_not_raised(self, qtbot, window):
        dialog = self._dialog(qtbot, window)
        dialog.start_spin.setValue(5.0)
        dialog.stop_spin.setValue(1.0)
        dialog.compute()
        assert dialog._scan is None
        assert "must increase" in dialog.status.text()

    def test_optical_constants_that_do_not_cover_the_range_are_reported(
            self, qtbot, window):
        """Au is tabulated 0.62-6.2 nm; asking outside it must not traceback."""
        dialog = self._dialog(qtbot, window)
        dialog.start_spin.setValue(400.0)
        dialog.stop_spin.setValue(700.0)
        dialog.count_spin.setValue(8)
        dialog.compute()
        assert dialog._scan is None
        assert dialog.status.text()
