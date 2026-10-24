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

    def test_an_asylum_ibw_loads_and_names_its_channel(self, qtbot, window, tmp_path):
        """Asylum files open in the window, with the direction selector live
        and the channel the reader chose shown, as a .spm's is."""
        from test_ibw_reader import _sawtooth, _write_ibw

        image = _sawtooth(period_nm=400.0, blaze_deg=25.0, antiblaze_deg=70.0,
                          scan_um=4.0, cols=512, rows=32)
        path = _write_ibw(tmp_path / "scan.ibw",
                          {"ZSensorRetrace": image, "ZSensorTrace": image},
                          dx_m=4e-6 / 512)
        window.load(str(path))
        assert window.run_btn.isEnabled()
        assert "Asylum: ZSensor / Retrace" in window.importer.info_label.text()
        assert window.importer.direction_combo.isEnabled()

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

        with open(target, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        assert lines[0] == "3 0 - Polygonal type"
        assert lines[1] == "Period: 1 PSC: 1"
        assert os.path.exists(str(tmp_path / "out_metrics.txt"))

    def test_efficiency_button_tracks_the_export_button(self, qtbot, window):
        """All three need one thing -- a profile -- so they must not drift apart."""
        assert not window.boundary.export_btn.isEnabled()
        assert not window.boundary.efficiency_btn.isEnabled()
        assert not window.boundary.send_btn.isEnabled()
        window.load(SAMPLE)
        assert window.boundary.export_btn.isEnabled()
        assert window.boundary.efficiency_btn.isEnabled()
        assert window.boundary.send_btn.isEnabled()

    def test_a_json_export_keeps_the_measured_period(self, qtbot, window, tmp_path):
        from gratinglab.io.boundary import read_boundary_json
        window.load(SAMPLE)
        path, metrics = window.boundary.write_profile(str(tmp_path / "out.json"))
        saved = read_boundary_json(path)
        assert saved.period_nm == window.boundary.profile.period_nm
        assert saved.provenance["source"] == window.boundary._filename
        assert metrics is None, "the JSON carries its metrics inside"

    def test_a_name_without_a_suffix_takes_the_chosen_format(self, qtbot, window, tmp_path):
        window.load(SAMPLE)
        json_filter = window.boundary.EXPORT_FORMATS[1][0]
        path, _ = window.boundary.write_profile(str(tmp_path / "groove"), json_filter)
        assert path.endswith(".json") and os.path.exists(path)

    def test_a_csv_export_writes_its_metrics_beside_it(self, qtbot, window, tmp_path):
        window.load(SAMPLE)
        path, metrics = window.boundary.write_profile(str(tmp_path / "out.csv"))
        with open(path, encoding="utf-8") as handle:
            assert handle.readline().strip() == "t,y,x_nm,y_nm"
        assert metrics and os.path.exists(metrics)

    def test_send_to_efficiency_opens_the_groove_with_its_period(self, qtbot, window):
        from gratinglab.gui.qt.main_window import MainWindow
        window.load(SAMPLE)
        efficiency = MainWindow()
        qtbot.addWidget(efficiency)
        with qtbot.waitSignal(efficiency.solved, timeout=60_000):
            pass  # the solve its constructor starts
        window.boundary._efficiency_window = efficiency
        efficiency.destroyed.connect(window.boundary._forget_efficiency_window)
        with qtbot.waitSignal(efficiency.solved, timeout=60_000):
            opened = window.boundary.send_to_efficiency()
        assert opened is efficiency
        form = efficiency.geometry.read_form()
        assert form.profile_kind == "From file"
        assert float(form.period) == pytest.approx(window.boundary.profile.period_nm, rel=1e-5)

        with qtbot.waitSignal(efficiency.solved, timeout=60_000):
            again = window.boundary.send_to_efficiency()
        assert again is efficiency, "a second send reuses the window"

    def test_a_destroyed_efficiency_window_is_forgotten(self, qtbot, window):
        """Closing and deleting it must not leave a dead reference behind."""
        from PySide6.QtCore import QCoreApplication, QEvent
        window.load(SAMPLE)
        opened = window.boundary.send_to_efficiency()
        with qtbot.waitSignal(opened.solved, timeout=60_000):
            pass
        opened.close()
        opened.deleteLater()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        qtbot.waitUntil(lambda: window.boundary._efficiency_window is None, timeout=5000)

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


class TestSwathDialog:
    """Per-swath groove statistics, opened from the Boundary tab."""

    def test_button_follows_the_profile(self, qtbot, window):
        assert not window.boundary.swaths_btn.isEnabled()
        window.load(SAMPLE)
        assert window.boundary.swaths_btn.isEnabled()

    def test_dialog_shows_statistics_for_the_loaded_scan(self, qtbot, window):
        window.load(SAMPLE)
        dialog = window.boundary.swath_dialog()
        qtbot.addWidget(dialog)
        assert dialog.stats is not None and dialog.stats.n_swaths >= 2
        assert "nm RMS" in dialog.summary_label.text()
        assert dialog.export_btn.isEnabled()

    def test_swath_count_reaches_the_statistics(self, qtbot, window):
        window.load(SAMPLE)
        dialog = window.boundary.swath_dialog()
        qtbot.addWidget(dialog)
        dialog.n_spin.setValue(4)
        assert dialog.stats.n_swaths == 4

    def test_boundary_settings_reach_the_swaths(self, qtbot, window):
        window.load(SAMPLE)
        window.boundary.n_points_spin.setValue(500)
        dialog = window.boundary.swath_dialog()
        qtbot.addWidget(dialog)
        assert dialog.stats.profiles.shape[1] == 500

    def test_export_writes_one_row_per_swath(self, qtbot, window, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QFileDialog
        window.load(SAMPLE)
        dialog = window.boundary.swath_dialog()
        qtbot.addWidget(dialog)
        dialog.n_spin.setValue(5)
        target = tmp_path / "swaths.csv"
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            lambda *a, **k: (str(target), "CSV table (*.csv)"))
        dialog.export()
        rows = [line for line in target.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")]
        assert len(rows) == 1 + 5


class TestTipCorrectionControls:
    """Tip deconvolution, from the window rather than config.py."""

    def test_the_tip_fields_follow_the_choice(self, qtbot, window):
        assert window.tip_combo.currentText() == 'none'
        assert not window.tip_radius_spin.isEnabled()
        window.tip_combo.setCurrentText('erosion')
        assert window.tip_radius_spin.isEnabled() and window.tip_angle_spin.isEnabled()

    def test_the_form_carries_them_into_the_settings(self, qtbot, window):
        from gratinglab.metrology.gui.state import build
        window.tip_combo.setCurrentText('erosion')
        window.tip_radius_spin.setValue(20.0)
        window.tip_angle_spin.setValue(18.0)
        settings, errors = build(window.form_state())
        assert errors == ()
        assert (settings.tip_correction, settings.tip_radius_nm,
                settings.tip_half_angle_deg) == ('erosion', 20.0, 18.0)

    def test_choosing_erosion_corrects_the_boundary_profile_now(self, qtbot, window):
        """The dilated fixture is the sawtooth seen through a 20 nm, 18 deg tip.

        Choosing erosion must reach the Boundary profile without a re-load and
        report how much of the image it recovered. It does NOT assert the depth
        recovers: the trough the tip never reached is exactly the part erosion
        cannot give back (core/tip.py), which an earlier draft of this test got
        wrong. The pointwise recovery is tests/metrology/test_tip.py's job.
        """
        import numpy as np
        from scans import SYNTHETIC_DILATED
        window.load(str(SYNTHETIC_DILATED))
        before = window.boundary.profile.y_norm.copy()
        window.tip_radius_spin.setValue(20.0)
        window.tip_angle_spin.setValue(18.0)
        window.tip_combo.setCurrentText('erosion')
        after = window.boundary.profile
        assert after.metrics['tip_correction'] == 'erosion'
        assert 0.5 < after.metrics['tip_certain_fraction'] < 1.0
        assert not np.array_equal(after.y_norm, before)
        assert "recovered exactly" in window.boundary.metrics_label.text()

    def test_any_analysis_setting_now_reaches_the_boundary_tab(self, qtbot, window):
        """It used to keep the settings current at load time."""
        window.load(SAMPLE)
        window.edge_spin.setValue(0.0)
        assert window.boundary._settings.edge_exclusion_periods == 0.0

    def test_the_results_say_the_scan_was_corrected(self, qtbot, window):
        window.load(SAMPLE)
        window.tip_combo.setCurrentText('erosion')
        _run_and_wait(qtbot, window)
        assert "Tip correction   : erosion" in window.results_label.text()


class TestTipShapePicker:
    """Named probes and asymmetric flanks, from the window."""

    def _pick(self, window, key):
        combo = window.tip_shape_combo
        combo.setCurrentIndex(combo.findData(key))

    def test_starts_as_a_symmetric_cone_with_the_right_flank_tied(self, qtbot, window):
        window.tip_combo.setCurrentText('erosion')
        assert window.tip_shape_combo.currentText() == window.TIP_SYMMETRIC
        assert not window.tip_angle_right_spin.isEnabled()
        window.tip_angle_spin.setValue(21.0)
        assert window.tip_angle_right_spin.value() == 21.0
        from gratinglab.metrology.gui.state import build
        settings, _ = build(window.form_state())
        assert settings.tip_half_angle_right_deg is None

    def test_a_preset_fills_the_numbers_and_reaches_the_settings(self, qtbot, window):
        from gratinglab.metrology.gui.state import build
        window.tip_combo.setCurrentText('erosion')
        self._pick(window, 'scanasyst-air-along-front-minus')
        assert (window.tip_radius_spin.value(), window.tip_angle_spin.value(),
                window.tip_angle_right_spin.value()) == (2.0, 15.0, 25.0)
        assert window.tip_angle_right_spin.isEnabled()
        settings, errors = build(window.form_state())
        assert errors == ()
        assert (settings.tip_radius_nm, settings.tip_half_angle_deg,
                settings.tip_half_angle_right_deg) == (2.0, 15.0, 25.0)

    def test_the_across_preset_is_symmetric_in_the_settings(self, qtbot, window):
        from gratinglab.metrology.gui.state import build
        window.tip_combo.setCurrentText('erosion')
        self._pick(window, 'scanasyst-air-across')
        settings, _ = build(window.form_state())
        assert (settings.tip_half_angle_deg, settings.tip_half_angle_right_deg) == (17.5, None)

    def test_editing_a_preset_makes_it_custom(self, qtbot, window):
        window.tip_combo.setCurrentText('erosion')
        self._pick(window, 'scanasyst-air-along-front-plus')
        window.tip_radius_spin.setValue(12.0)    # a worn tip, Bruker's maximum
        assert window.tip_shape_combo.currentText() == window.TIP_CUSTOM
        assert window.tip_angle_spin.value() == 25.0

    def test_an_asymmetric_tip_reaches_the_boundary_tab(self, qtbot, window):
        window.load(SAMPLE)
        window.tip_combo.setCurrentText('erosion')
        self._pick(window, 'scanasyst-air-along-front-minus')
        m = window.boundary.profile.metrics
        assert (m['tip_half_angle_deg'], m['tip_half_angle_right_deg']) == (15.0, 25.0)
        assert "15 deg (-x) / 25 deg (+x)" in window.boundary.metrics_label.text()


class TestGarbageCollectionStaysOnTheGuiThread:
    """CI run 39835186938 segfaulted on macOS with "Garbage-collecting" at the
    top of the analysis worker's stack: the cycle collector ran on the worker
    thread and destroyed Qt objects there. gratinglab.gui.qt.gc_guard keeps it
    on the GUI thread while a window exists."""

    def test_no_collection_runs_on_the_worker_during_an_analysis(self, qtbot, window):
        import gc
        import threading

        gui_thread = threading.get_ident()
        seen = []

        def record(phase, info):
            if phase == "start":
                seen.append(threading.get_ident())

        window.load(SAMPLE)
        gc.callbacks.append(record)
        try:
            # Plenty of cyclic garbage, so the collector has reason to run
            # somewhere while the worker allocates.
            for _ in range(3):
                junk = [[] for _ in range(5000)]
                for a, b in zip(junk, junk[1:]):
                    a.append(b)
                    b.append(a)
                del junk
                _run_and_wait(qtbot, window)
            qtbot.wait(1200)  # let the GUI-thread collector take its turns
        finally:
            gc.callbacks.remove(record)
        assert seen, "no collection ran at all, so the test proves nothing"
        assert set(seen) == {gui_thread}, (
            f"{sum(t != gui_thread for t in seen)} of {len(seen)} collections "
            "ran off the GUI thread")
