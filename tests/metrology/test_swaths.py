"""
Groove shape statistics across swaths.

The committed fixture has the same groove in every row, so it can only show
that identical swaths read as identical. The scans built here put a known
change down the image, a blaze angle that varies from top to bottom or a
grating that drifts sideways, so each statistic has a known answer to recover.
"""
import contextlib
import io
import os
import sys

import numpy as np
import pytest

from scans import SYNTHETIC

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

from gratinglab.metrology.boundary import SwathStatistics, swath_statistics    # noqa: E402
from gratinglab.metrology.core.image_flatten import flatten_image               # noqa: E402
from gratinglab.metrology.core.processing import load_afm_data                  # noqa: E402
from gratinglab.metrology.settings import AnalysisSettings                      # noqa: E402

PERIOD_NM = 315.0
SCAN_X_UM = 2.0
N_ROWS, N_COLS = 60, 384


def _sawtooth(t, blaze_deg, antiblaze_deg=70.0):
    cot_b = 1.0 / np.tan(np.radians(blaze_deg))
    cot_a = 1.0 / np.tan(np.radians(antiblaze_deg))
    apex = cot_b / (cot_b + cot_a)
    depth = 1.0 / (cot_b + cot_a)
    t = np.mod(-t, 1.0)
    return np.where(t <= apex, t / apex * depth, (1.0 - t) / (1.0 - apex) * depth)


def _scan(blaze_top=30.0, blaze_bottom=30.0, drift_nm=0.0, noise_nm=0.15, seed=1):
    """Heights in metres; blaze and lateral offset vary linearly down the rows."""
    rng = np.random.default_rng(seed)
    x_nm = np.linspace(0.0, SCAN_X_UM * 1000.0, N_COLS, endpoint=False)
    rows = []
    for i, frac in enumerate(np.linspace(0.0, 1.0, N_ROWS)):
        blaze = blaze_top + (blaze_bottom - blaze_top) * frac
        shift = drift_nm * frac
        rows.append(_sawtooth((x_nm - shift) / PERIOD_NM, blaze) * PERIOD_NM)
    rows = np.array(rows) + rng.normal(0.0, noise_nm, size=(N_ROWS, N_COLS))
    return rows * 1e-9


SETTINGS = AnalysisSettings.from_config()


def test_identical_swaths_read_as_identical():
    with contextlib.redirect_stdout(io.StringIO()):
        data, size = load_afm_data(str(SYNTHETIC), default_scan_size=2.0)
    stats = swath_statistics(flatten_image(data, SETTINGS.image_flatten_method), size, SETTINGS)
    assert isinstance(stats, SwathStatistics)
    assert stats.n_swaths >= 2 and not stats.skipped
    # Every row is the same groove, so swaths differ only by noise.
    assert stats.between_rms_nm < 0.5
    assert stats.metric_summary()['max_angle_deg']['std'] < 0.5
    assert np.ptp(stats.shift_nm) < 1.0


def test_a_blaze_gradient_down_the_image_is_found():
    stats = swath_statistics(_scan(blaze_top=26.0, blaze_bottom=34.0), SCAN_X_UM, SETTINGS,
                             n_swaths=6)
    depth = [s['groove_depth'] for s in stats.per_swath]
    # A steeper blaze is a deeper groove, so depth rises monotonically down the scan.
    assert all(b > a for a, b in zip(depth, depth[1:]))
    uniform = swath_statistics(_scan(), SCAN_X_UM, SETTINGS, n_swaths=6)
    assert stats.between_rms_nm > 5 * uniform.between_rms_nm
    # The change is along the grooves, so it dominates the groove-to-groove spread.
    assert stats.between_rms_nm > stats.within_rms_nm


def test_a_lateral_drift_is_reported_as_shift_not_shape():
    drift = 60.0
    stats = swath_statistics(_scan(drift_nm=drift), SCAN_X_UM, SETTINGS, n_swaths=6)
    shifts = stats.shift_nm
    assert abs(abs(shifts[-1] - shifts[0]) - drift * (1 - 1 / 6)) < 15.0
    assert np.all(np.diff(shifts) * np.sign(shifts[-1] - shifts[0]) >= -1.0)
    # Aligned, the drifting grooves are the same shape.
    assert stats.between_rms_nm < 1.0
    # But averaging every row blurs them, so the whole-scan groove is rounder.
    still = swath_statistics(_scan(), SCAN_X_UM, SETTINGS, n_swaths=6)
    assert still.whole_vs_swaths_rms_nm < 1.0
    assert stats.whole_vs_swaths_rms_nm > 5 * max(still.whole_vs_swaths_rms_nm, 0.2)


def test_default_swath_count_follows_the_row_group_setting():
    stats = swath_statistics(_scan(), SCAN_X_UM, SETTINGS.with_(n_row_groups=4))
    assert stats.n_swaths == 4
    assert stats.row_ranges[0][0] == 0 and stats.row_ranges[-1][1] == N_ROWS


def test_aligned_profiles_stay_closed():
    stats = swath_statistics(_scan(drift_nm=40.0), SCAN_X_UM, SETTINGS, n_swaths=5)
    assert stats.profiles.shape == (5, len(stats.t))
    np.testing.assert_array_equal(stats.profiles[:, 0], stats.profiles[:, -1])


def test_too_few_usable_swaths_is_an_error():
    with pytest.raises(ValueError, match="nothing to compare"):
        swath_statistics(_scan(), SCAN_X_UM, SETTINGS, n_swaths=1)


def test_summary_and_csv(tmp_path):
    stats = swath_statistics(_scan(blaze_top=28.0, blaze_bottom=32.0), SCAN_X_UM, SETTINGS,
                             n_swaths=5)
    text = stats.summary()
    assert "5 swaths" in text and "nm RMS" in text
    path = stats.write_csv(tmp_path / "swaths.csv")
    raw = path.read_bytes().decode("utf-8")
    assert raw.splitlines()[0].startswith("first_row,end_row,shift_nm")
    table = np.genfromtxt(path, delimiter=",", names=True, comments="#")
    assert len(table) == 5
    np.testing.assert_allclose(table['groove_depth'],
                               [s['groove_depth'] for s in stats.per_swath])
