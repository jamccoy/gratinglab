"""
Profile flattening, the stage that does move the answer.

`level_grooves` fits its background through detected features, which makes the
size of that feature set load-bearing. Two ways it used to go wrong silently:

- Edge exclusion ran *after* detection and overwrote the feature arrays, so a
  trace whose features all sat near the ends was fitted on whatever remained -
  sometimes a single point - while the code printed "using all".
- The guard was `len(feature_x) < 2`, but `poly_order` defaults to 2 and a
  degree-k polynomial needs k+1 points. Two features and a quadratic gave a
  rank-deficient background behind a numpy RankWarning.

Both produce a background that looks like a fit and is not one, and because
flattening sets the plane every blaze angle on the scan is quoted against, a
bad background is a bias on every groove rather than noise on one. These tests
pin the handling.

Each fallback raises a `UserWarning`. The suite runs with
`filterwarnings = error`, so a test that provokes one has to say so, and any
other test that starts flattening by fallback fails rather than passing
quietly.

The profiles here are built rather than measured: a few smooth bumps at known
positions, so the feature set is exactly what the test intends and the
expected background can be written down.
"""
from pathlib import Path

import numpy as np
import pytest

from gratinglab.metrology.analyzer import _aggregate_flattening
from gratinglab.metrology.core.processing import (
    FlatteningReport,
    flatten_profile,
)

#: x in micrometres, 1 nm per sample - `flatten_profile` reads dx in nm as
#: `(x[1] - x[0]) * 1000`, and the feature separation gate is in periods.
DX_UM = 0.001
PERIOD_NM = 100.0


def _profile(centres_um, heights=None, n=1001, width_um=0.01):
    """A baseline of zero with a smooth bump at each centre."""
    x = np.arange(n) * DX_UM
    heights = [1.0] * len(centres_um) if heights is None else heights
    y = np.zeros_like(x)
    for centre, height in zip(centres_um, heights):
        y = y + height * np.exp(-(((x - centre) / width_um) ** 2))
    return x, y


def _background_through(x, feature_x, feature_y, order):
    """The polynomial a correct implementation must have fitted."""
    return np.polyval(np.polyfit(feature_x, feature_y, order), x)


class TestEdgeExclusionDoesNotStarveTheFit:
    """
    Excluding the edges is an improvement to the fit, not a licence to break it.

    When exclusion would leave too few features for the requested order, the
    full set is used and the warning says so. The old code said so and then
    fitted the remnant anyway.
    """

    CENTRES = (0.05, 0.5, 0.95)

    def test_all_features_are_used_when_exclusion_would_leave_too_few(self):
        x, y = _profile(self.CENTRES)

        with pytest.warns(UserWarning, match='using all of them'):
            flat, background = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                exclude_edges=0.2, feature='peaks', period_nm=PERIOD_NM)

        # exclude_edges=0.2 keeps only 0.2 <= x <= 0.8, i.e. the middle bump.
        # A degree-2 background needs three points, so all three are used.
        expected = _background_through(
            x, np.array(self.CENTRES), np.ones(3), 2)

        assert np.allclose(background, expected, atol=1e-9), (
            "the background was not fitted through every detected feature")
        assert np.allclose(flat, y - expected, atol=1e-9)

    def test_the_warning_says_what_actually_happened(self):
        x, y = _profile(self.CENTRES)

        with pytest.warns(UserWarning) as caught:
            flatten_profile(x, y, method='level_grooves', poly_order=2,
                            exclude_edges=0.2, feature='peaks',
                            period_nm=PERIOD_NM)

        message = str(caught[0].message)
        assert 'edge exclusion leaves 1 feature(s)' in message, message
        assert 'using all of them' in message, message

    def test_exclusion_still_bites_when_enough_features_survive(self):
        """
        The fix must not quietly disable the setting it guards.

        Eleven evenly spaced bumps, the outer ones raised: excluding the edges
        has to change the background, or edge exclusion has stopped working.
        """
        centres = [0.05 + 0.09 * i for i in range(11)]
        heights = [2.0] + [1.0] * 9 + [2.0]
        x, y = _profile(centres, heights)

        # Neither call may warn: eleven features survive either way, so a
        # warning here would mean the fit was starved after all.
        _, excluded = flatten_profile(
            x, y, method='level_grooves', poly_order=2,
            exclude_edges=0.2, feature='peaks', period_nm=PERIOD_NM)
        _, everything = flatten_profile(
            x, y, method='level_grooves', poly_order=2,
            exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)

        assert not np.allclose(excluded, everything, atol=1e-6), (
            "excluding the edges made no difference, so it was not applied")


class TestOrderIsReducedRatherThanRankDeficient:
    """
    Two features cannot determine a quadratic.

    numpy will return something regardless, behind a RankWarning that nothing
    in this pipeline was reading. Reducing the order gives the best background
    the data actually supports.
    """

    def test_two_features_give_the_line_through_them(self):
        x, y = _profile((0.25, 0.75))

        with pytest.warns(UserWarning, match='degree-1 background'):
            _, background = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)

        expected = _background_through(x, np.array((0.25, 0.75)), np.ones(2), 1)
        assert np.allclose(background, expected, atol=1e-9)

        # A straight line, not a curve that happened to pass through two points.
        assert np.allclose(np.diff(background, 2), 0.0, atol=1e-12)

    def test_the_reduction_is_reported(self):
        x, y = _profile((0.25, 0.75))

        with pytest.warns(UserWarning) as caught:
            flatten_profile(x, y, method='level_grooves', poly_order=2,
                            exclude_edges=0.0, feature='peaks',
                            period_nm=PERIOD_NM)

        message = str(caught[0].message)
        assert 'support only a degree-1 background' in message, message


class TestTheUnaffectedPathsStayPut:
    """
    The fix touches one branch. These pin the behaviour around it.
    """

    def test_ample_features_fit_the_requested_order(self):
        centres = [0.05 + 0.09 * i for i in range(11)]
        x, y = _profile(centres)

        _, background = flatten_profile(
            x, y, method='level_grooves', poly_order=2,
            exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)

        expected = _background_through(
            x, np.array(centres), np.ones(len(centres)), 2)
        assert np.allclose(background, expected, atol=1e-9)

    def test_too_few_features_still_falls_back_to_linear(self):
        """
        One bump cannot be levelled on its features, and the outer guard
        catches that before the branch this change touches.
        """
        x, y = _profile((0.5,))

        with pytest.warns(UserWarning, match='using linear'):
            _, background = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)
        assert np.allclose(background, np.polyval(np.polyfit(x, y, 1), x))

    def test_level_grooves_still_demands_a_period(self):
        x, y = _profile((0.25, 0.75))
        with pytest.raises(ValueError, match='period_nm'):
            flatten_profile(x, y, method='level_grooves', period_nm=None)


class TestTheReport:
    """
    `return_report=True` says what was actually done.

    The warnings announce a fallback as it happens; the report is the same
    information in a form the caller can store beside the result, which is what
    makes it survive a batch run.
    """

    def test_a_clean_fit_reports_no_fallback(self):
        centres = [0.05 + 0.09 * i for i in range(11)]
        x, y = _profile(centres)

        _, _, report = flatten_profile(
            x, y, method='level_grooves', poly_order=2, exclude_edges=0.0,
            feature='peaks', period_nm=PERIOD_NM, return_report=True)

        assert report.method_used == report.method_requested == 'level_grooves'
        assert report.poly_order_used == report.poly_order_requested == 2
        assert report.n_features == 11
        assert not report.edge_exclusion_skipped
        assert not report.fell_back
        assert report.dof == 11 - 3

    def test_the_residual_measures_the_lands_not_the_grooves(self):
        """
        Equal-height lands lie on a background exactly; uneven ones do not.

        This is the quantity worth having: how far from a smooth low-order
        surface the features the reference was fitted through actually lie.
        """
        centres = [0.05 + 0.09 * i for i in range(11)]
        x_flat, y_flat_lands = _profile(centres)
        _, _, flat_report = flatten_profile(
            x_flat, y_flat_lands, method='level_grooves', poly_order=2,
            feature='peaks', period_nm=PERIOD_NM, return_report=True)

        # Same lands, one of them 5 nm high.
        heights = [1.0] * 11
        heights[5] = 1.05
        x_bumpy, y_bumpy = _profile(centres, heights)
        _, _, bumpy_report = flatten_profile(
            x_bumpy, y_bumpy, method='level_grooves', poly_order=2,
            feature='peaks', period_nm=PERIOD_NM, return_report=True)

        assert flat_report.residual_rms_nm < 1e-6, "coplanar lands should fit"
        assert bumpy_report.residual_rms_nm > flat_report.residual_rms_nm

    def test_whole_trace_methods_report_no_residual(self):
        """
        `linear` and `polynomial` fit the trace, not features.

        Their residual would be the profile about its own background, dominated
        by the groove depth, and would say nothing about how flat the reference
        is. None is the honest answer.
        """
        x, y = _profile([0.05 + 0.09 * i for i in range(11)])

        for method in ('linear', 'polynomial'):
            _, _, report = flatten_profile(
                x, y, method=method, poly_order=2, return_report=True)
            assert report.residual_rms_nm is None, method
            assert report.n_features == 0, method
            assert not report.fell_back, method

    def test_an_exact_fit_has_no_degrees_of_freedom(self):
        """
        Three points and a quadratic give a residual of zero by arithmetic.

        Reporting that as flatness would be wrong, so `dof` is zero and
        `summary` says the number is meaningless.
        """
        x, y = _profile((0.05, 0.5, 0.95))

        _, _, report = flatten_profile(
            x, y, method='level_grooves', poly_order=2, exclude_edges=0.0,
            feature='peaks', period_nm=PERIOD_NM, return_report=True)

        assert report.n_features == 3
        assert report.dof == 0
        assert report.residual_rms_nm == pytest.approx(0.0, abs=1e-9)
        assert 'meaningless' in report.summary()

    def test_a_skipped_exclusion_is_a_fallback(self):
        x, y = _profile((0.05, 0.5, 0.95))

        with pytest.warns(UserWarning, match='using all of them'):
            _, _, report = flatten_profile(
                x, y, method='level_grooves', poly_order=2, exclude_edges=0.2,
                feature='peaks', period_nm=PERIOD_NM, return_report=True)

        assert report.edge_exclusion_skipped
        assert report.fell_back
        assert 'edge exclusion skipped' in report.summary()

    def test_a_method_fallback_is_named_in_the_summary(self):
        x, y = _profile((0.5,))

        with pytest.warns(UserWarning, match='using linear'):
            _, _, report = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                feature='peaks', period_nm=PERIOD_NM, return_report=True)

        assert report.method_used == 'linear'
        assert report.fell_back
        assert 'asked for level_grooves' in report.summary()

    def test_two_values_still_unpack(self):
        """
        The default is unchanged, which is why the GUI needed no edit.
        """
        x, y = _profile([0.05 + 0.09 * i for i in range(11)])
        result = flatten_profile(x, y, method='level_grooves',
                                 feature='peaks', period_nm=PERIOD_NM)
        assert len(result) == 2


class TestAggregationAcrossRowGroups:
    """One sample's flattening, summarised from its row groups."""

    def _report(self, **kw):
        defaults = dict(method_requested='level_grooves',
                        method_used='level_grooves',
                        poly_order_requested=2, poly_order_used=2,
                        n_features=11, edge_exclusion_skipped=False,
                        residual_rms_nm=0.5)
        return FlatteningReport(**{**defaults, **kw})

    def test_the_worst_group_sets_the_residual(self):
        """Trust is set by the weakest reference, not the average one."""
        got = _aggregate_flattening([
            self._report(residual_rms_nm=0.2),
            self._report(residual_rms_nm=1.7),
            self._report(residual_rms_nm=0.4),
        ])
        assert got['flatten_residual_rms_nm'] == pytest.approx(1.7)
        assert got['flatten_residual_groups'] == 3
        assert got['flatten_fallback_groups'] == 0

    def test_exact_fits_are_left_out_of_the_worst_case(self):
        """
        A group with no spare features reports zero, which is not flatness.

        Including it would drag the reported worst case towards a number
        nothing measured.
        """
        got = _aggregate_flattening([
            self._report(n_features=3, residual_rms_nm=0.0),
            self._report(residual_rms_nm=0.9),
        ])
        assert got['flatten_residual_rms_nm'] == pytest.approx(0.9)
        assert got['flatten_residual_groups'] == 1

    def test_all_exact_fits_report_no_residual_at_all(self):
        got = _aggregate_flattening([self._report(n_features=3,
                                                  residual_rms_nm=0.0)])
        assert got['flatten_residual_rms_nm'] is None
        assert got['flatten_residual_groups'] == 0

    def test_fallbacks_are_counted(self):
        got = _aggregate_flattening([
            self._report(),
            self._report(method_used='linear', poly_order_used=1),
            self._report(edge_exclusion_skipped=True),
        ])
        assert got['flatten_fallback_groups'] == 2
        assert got['flatten_method'] == 'level_grooves'

    def test_nothing_to_aggregate_adds_no_keys(self):
        """
        An older result dict carries no reports, and must not gain empty
        columns that would read as measurements.
        """
        assert _aggregate_flattening([]) == {}
        assert _aggregate_flattening([None, None]) == {}


#: Enough of a result dict for the writers to run. Only the flattening keys
#: matter here; the rest are plausible values the summary needs present.
def _result(**flattening):
    quality = [{'groove_depth_nm': 100.0, 'blaze_width_nm': 250.0,
                'steep_width_nm': 60.0, 'blaze_r2': 0.999,
                'local_period_used_nm': 315.0}]
    result = {
        'filename': 'synthetic.txt', 'n_grooves': 1, 'n_groups': 2,
        'mean_angle': 30.0, 'std_angle': 0.1,
        'min_angle': 29.9, 'max_angle': 30.1,
        'mean_slope': 0.577, 'mean_steep': -1.7, 'period_nm': 315.0,
        'period_std': 0.4, 'groove_periods': [315.0], 'all_angles': [30.0],
        'groove_row_groups': [0], 'quality': quality,
        'local_angle_std': 0.05, 'local_angle_range': 0.2,
        'all_local_angles': [[30.0]], 'total_std': 0.1, 'sem': 0.05,
        'measurement_variance': 0.001, 'physical_variance': 0.002,
        'mean_measurement_uncertainty': 0.01, 'within_image_std': 0.03,
        'icc': 0.24, 'design_effect': 1.2, 'n_effective': 0.8,
        'sem_corrected': 0.07,
    }
    result.update(flattening)
    return result


class TestItReachesTheOutput:
    """
    The report is only worth having if it survives the run.

    A warning in a console that has scrolled is what this replaces, so these
    tests run the writers and read the files back.
    """

    def _write(self, tmp_path, result):
        from gratinglab.metrology.io.file_io import save_results_to_file

        txt, csv_path, _ = save_results_to_file(
            [result], labels=['S1'], output_dir=str(tmp_path))
        return Path(txt).read_text(encoding='utf-8'), Path(csv_path).read_text(
            encoding='utf-8')

    def test_the_csv_carries_the_columns_and_the_values(self, tmp_path):
        summary, csv_text = self._write(tmp_path, _result(
            flatten_method='level_grooves', flatten_fallback_groups=1,
            flatten_residual_rms_nm=0.42, flatten_residual_groups=2))

        header, row = csv_text.splitlines()[:2]
        assert header.endswith("Flatten_method,Flatten_fallback_groups,"
                               "Flatten_residual_nm,Flatten_residual_groups")
        assert row.split(',')[-4:] == ['level_grooves', '1', '0.4200', '2']

    def test_the_columns_are_appended_not_inserted(self, tmp_path):
        """
        Existing columns keep their positions, so a baseline recorded before
        this change still diffs against one recorded after.
        """
        _, csv_text = self._write(tmp_path, _result(
            flatten_method='level_grooves', flatten_fallback_groups=0,
            flatten_residual_rms_nm=0.1, flatten_residual_groups=2))

        columns = csv_text.splitlines()[0].split(',')
        assert columns[:4] == ['Sample', 'File', 'N_scans', 'N_grooves']
        assert columns.index('SEM_corrected_deg') == columns.index(
            'Flatten_method') - 1

    def test_the_summary_says_whether_every_group_got_the_method(self, tmp_path):
        clean, _ = self._write(tmp_path, _result(
            flatten_method='level_grooves', flatten_fallback_groups=0,
            flatten_residual_rms_nm=0.2, flatten_residual_groups=3))
        assert 'every row group as asked' in clean

        fell, _ = self._write(tmp_path, _result(
            flatten_method='level_grooves', flatten_fallback_groups=2,
            flatten_residual_rms_nm=0.2, flatten_residual_groups=3))
        assert '2 row group(s) fell back' in fell

    def test_the_summary_reports_the_residual_with_its_group_count(self, tmp_path):
        summary, _ = self._write(tmp_path, _result(
            flatten_method='level_grooves', flatten_fallback_groups=0,
            flatten_residual_rms_nm=0.42, flatten_residual_groups=2))

        assert 'Lands about the background: 0.420 nm RMS, worst of 2 group(s)' \
            in summary

    def test_a_result_without_the_fields_writes_blanks_not_zeros(self, tmp_path):
        """
        An older result dict, or the boundary path, carries no report. A zero
        in these columns would read as "no fallbacks, perfectly flat", which is
        a measurement nobody made.
        """
        summary, csv_text = self._write(tmp_path, _result())

        assert csv_text.splitlines()[1].endswith(',,,')
        assert 'Flattening:' not in summary
        assert 'Lands about the background' not in summary
