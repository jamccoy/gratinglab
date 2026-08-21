"""
Tests for the ICC calculation and for row-group label integrity.

The ICC decides whether the project's reported uncertainties need reworking, so
the arithmetic is checked against synthetic data with known structure before any
real number is trusted.

Run directly:
    .venv/bin/python tests/test_icc.py
"""
import contextlib
import io
import os
import sys

import numpy as np

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import matplotlib
matplotlib.use('Agg')

from gratinglab.metrology.stats.icc import (
    compute_icc, effective_sample_size, sem_inflation_factor, interpret_icc)


# ── The maths, on data whose structure is known by construction ──────────────

def test_independent_data_gives_icc_near_zero():
    """Draws from one distribution, arbitrarily grouped, are uncorrelated"""
    rng = np.random.default_rng(42)
    values = rng.normal(30.0, 1.0, 400)
    labels = np.repeat(np.arange(20), 20)

    icc = compute_icc(values, labels)['icc']
    assert icc < 0.1, f"independent data should give ICC near 0, got {icc:.4f}"


def test_strong_group_offsets_give_high_icc():
    """A large per-group offset puts nearly all variance between groups"""
    rng = np.random.default_rng(0)
    offsets = rng.normal(0, 5.0, 20)          # between-group spread
    values, labels = [], []
    for g, off in enumerate(offsets):
        values.extend(30.0 + off + rng.normal(0, 0.1, 20))  # tiny within-group
        labels.extend([g] * 20)

    icc = compute_icc(values, labels)['icc']
    assert icc > 0.9, f"strongly grouped data should give ICC near 1, got {icc:.4f}"


def test_icc_matches_hand_computed_value():
    """Exact check against the variance ratio computed by hand"""
    values = np.array([1.0, 2.0, 11.0, 12.0])
    labels = np.array([0, 0, 1, 1])
    # group means 1.5 and 11.5 -> between_var = var([1.5, 11.5], ddof=1) = 50
    # each group var (ddof=1) = 0.5 -> within_var = 0.5
    # icc = 50 / 50.5
    got = compute_icc(values, labels)
    assert np.isclose(got['between_var'], 50.0), got['between_var']
    assert np.isclose(got['within_var'], 0.5), got['within_var']
    assert np.isclose(got['icc'], 50.0 / 50.5), got['icc']


def test_single_group_is_undefined_not_zero():
    """One group cannot yield a between-group variance; must not report 0"""
    got = compute_icc([1.0, 2.0, 3.0], [0, 0, 0])
    assert np.isnan(got['icc']), \
        "a single group must give nan, not 0 - 0 would read as 'independent'"


def test_length_mismatch_raises():
    try:
        compute_icc([1.0, 2.0, 3.0], [0, 0])
    except ValueError:
        return
    raise AssertionError("mismatched lengths should raise ValueError")


def test_design_effect_relationships():
    """n_eff falls to the group count as ICC rises; inflation is its sqrt"""
    n, m = 400, 20
    assert np.isclose(effective_sample_size(n, m, 0.0), n)
    assert np.isclose(sem_inflation_factor(m, 0.0), 1.0)

    # At ICC = 1 every group contributes one independent observation
    assert np.isclose(effective_sample_size(n, m, 1.0), n / m)

    # Monotonic in ICC
    assert (effective_sample_size(n, m, 0.05) >
            effective_sample_size(n, m, 0.30))
    assert sem_inflation_factor(m, 0.30) > sem_inflation_factor(m, 0.05)


def test_interpretation_thresholds():
    assert 'negligible' in interpret_icc(0.05)
    assert 'moderate' in interpret_icc(0.15)
    assert 'substantial' in interpret_icc(0.40)


# ── Label integrity on real data ─────────────────────────────────────────────

def _analyze(path, trim):
    """
    Run an analysis at a given facet trim.

    Settings are passed, not patched. This used to rebind analyzer.FACET_TRIM and
    restore it in a finally block, which worked only because the module bound its
    configuration at import - the same coupling that made a CLI impossible.
    """
    from gratinglab.metrology.analyzer import analyze_single_file
    from gratinglab.metrology.settings import AnalysisSettings

    settings = AnalysisSettings.from_config().with_(facet_trim=trim)
    with contextlib.redirect_stdout(io.StringIO()):
        return analyze_single_file(path, show_plots=False, settings=settings)


def test_row_group_labels_align_with_angles():
    """Every measurement carries a label, at the default settings"""
    from gratinglab.metrology.config import PROJECT_ROOT
    path = os.path.join(PROJECT_ROOT, 'data', '20250820_280C_00004.txt')
    r = _analyze(path, 0.10)
    assert len(r['groove_row_groups']) == len(r['all_angles'])


def test_labels_stay_aligned_when_fits_fail():
    """
    The regression test for the bug this work replaced.

    Group membership used to be reconstructed by slicing the angle list using
    each group's detected-centre count. Angles are only appended when a fit
    succeeds, so any failure shifted every later group. At FACET_TRIM = 0.28 this
    file detects 102 centres but only ~70 fits succeed - exactly the case the old
    code got wrong.
    """
    from gratinglab.metrology.config import PROJECT_ROOT
    path = os.path.join(PROJECT_ROOT, 'data', '20250820_280C_00004.txt')
    r = _analyze(path, 0.28)
    assert r is not None, "expected some measurements to survive at trim 0.28"
    n_angles, n_labels = len(r['all_angles']), len(r['groove_row_groups'])
    assert n_angles == n_labels, f"{n_angles} angles but {n_labels} labels"
    assert n_angles < 102, ("expected some fits to fail at trim 0.28, so this "
                            "test exercises the misalignment case")


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
