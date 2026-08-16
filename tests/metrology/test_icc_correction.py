"""
The correlation correction applied to real analysis output.

`test_icc.py` covers the ICC arithmetic itself on synthetic data. This file
checks the wiring: that the analyzer's corrected SEM follows from the ICC it
reports, that the correction vanishes when it should, and that it can only ever
make a result less significant.
"""
import contextlib
import io
import os
import sys

import numpy as np

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

import matplotlib
matplotlib.use('Agg')

from afm_analysis.analyzer import analyze_single_file          # noqa: E402
from afm_analysis.config import PROJECT_ROOT                   # noqa: E402
from afm_analysis.settings import AnalysisSettings             # noqa: E402
from afm_analysis.stats.icc import (                           # noqa: E402
    effective_sample_size, sem_inflation_factor)
from afm_analysis.stats.analysis import _effective_n           # noqa: E402

SAMPLE = os.path.join(PROJECT_ROOT, 'data', 'ALD_master_1p5um_flatten.txt')

_cache = {}


def _result(path=SAMPLE, **overrides):
    key = (path, tuple(sorted(overrides.items())))
    if key not in _cache:
        settings = AnalysisSettings.from_config().with_(**overrides)
        with contextlib.redirect_stdout(io.StringIO()):
            _cache[key] = analyze_single_file(path, show_plots=False,
                                              settings=settings)
    return _cache[key]


# ── The correction follows from the ICC ──────────────────────────────────────

def test_result_carries_the_correction_fields():
    r = _result()
    for key in ('icc', 'design_effect', 'n_effective', 'sem_corrected'):
        assert key in r, f"{key} missing from the result"
        assert r[key] is not None


def test_effective_size_matches_the_design_effect_formula():
    r = _result()
    expected = effective_sample_size(r['n_grooves'],
                                     r['n_grooves'] / r['n_groups'],
                                     r['icc'])
    assert np.isclose(r['n_effective'], expected)


def test_corrected_sem_is_the_plain_sem_times_the_inflation_factor():
    """
    The two code paths must agree.

    stats/icc.py derives an inflation factor from the ICC; the analyzer divides
    total_std by sqrt(n_eff). If those ever disagree the reported SEM and the ICC
    report would tell different stories about the same data.
    """
    r = _result()
    factor = sem_inflation_factor(r['n_grooves'] / r['n_groups'], r['icc'])
    assert np.isclose(r['sem_corrected'], r['sem'] * factor, rtol=1e-9)


def test_correction_never_shrinks_the_error_bar():
    """ICC is non-negative here, so the corrected SEM cannot be the smaller one"""
    for path in ('ALD_master_1p5um_flatten.txt', '20250820_215C_00001.txt',
                 '500C_N2_flatten.txt'):
        r = _result(os.path.join(PROJECT_ROOT, 'data', path))
        assert r['sem_corrected'] >= r['sem'], path


def test_effective_size_never_exceeds_the_measurement_count():
    r = _result()
    assert r['n_effective'] <= r['n_grooves']
    assert r['n_effective'] >= r['n_groups'] * 0.5, (
        "effective size collapsed below the group count - suspect a bad ICC")


# ── The correction vanishes when it is not needed ────────────────────────────

def test_zero_icc_leaves_the_sem_untouched():
    """
    A correction that does not disappear at ICC = 0 is wrong.

    Checked on the formula rather than by manufacturing uncorrelated AFM data:
    n_eff = N / (1 + (m-1)*0) = N, so sem_corrected = total_std/sqrt(N) = sem.
    """
    n, m = 100, 5
    assert effective_sample_size(n, m, 0.0) == n
    assert sem_inflation_factor(m, 0.0) == 1.0


def test_traditional_mode_falls_back_to_the_raw_count():
    """
    Without row groups there is no clustering to correct for.

    _effective_n must return the real count rather than nan, or every downstream
    standard error becomes nan.
    """
    r = _result(use_row_groups=False)
    n_eff = _effective_n(r)
    assert np.isfinite(n_eff)
    assert n_eff == float(r['n_grooves'])


def test_effective_n_handles_results_predating_the_field():
    """Old result dicts must not break the comparison code"""
    legacy = {'n_grooves': 42, 'std_angle': 1.0, 'mean_angle': 30.0}
    assert _effective_n(legacy) == 42.0
    assert _effective_n({**legacy, 'n_effective': float('nan')}) == 42.0
    assert _effective_n({**legacy, 'n_effective': 0}) == 42.0


# ── Effect on inference ──────────────────────────────────────────────────────

def test_no_p_value_gets_smaller():
    """
    The correction removes information; it can never add confidence.

    Any comparison whose p-value falls after correcting is a sign the effective
    sizes went the wrong way.
    """
    from scipy import stats as scipy_stats
    from afm_analysis.stats.analysis import _calculate_welch_df

    a = _result(os.path.join(PROJECT_ROOT, 'data', 'ALD_master_1p5um_flatten.txt'))
    b = _result(os.path.join(PROJECT_ROOT, 'data', '500C_N2_flatten.txt'))

    def p_value(use_effective):
        s1, s2 = a['total_std'], b['total_std']
        n1 = a['n_effective'] if use_effective else a['n_grooves']
        n2 = b['n_effective'] if use_effective else b['n_grooves']
        se = np.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2)
        t = (b['mean_angle'] - a['mean_angle']) / se
        df = _calculate_welch_df(s1, n1, s2, n2)
        return 2 * (1 - scipy_stats.t.cdf(abs(t), df))

    assert p_value(True) >= p_value(False)


def test_means_are_untouched_by_the_correction():
    """The whole point: uncertainties change, the measurement does not"""
    r = _result()
    assert np.isclose(r['mean_angle'], 33.2337, atol=1e-4), (
        "the master sample's mean moved; the correction must only affect "
        "uncertainty")


if __name__ == '__main__':
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith('test_'):
            continue
        try:
            fn()
            print(f"PASS  {name}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL  {name}: {exc}")
    print(f"\n{'all tests passed' if not failures else f'{failures} failure(s)'}")
    sys.exit(1 if failures else 0)
