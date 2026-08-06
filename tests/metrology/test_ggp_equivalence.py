"""
Regression test: the ported boundary-profile code must still reproduce the
original standalone script bit-for-bit.

fixtures/rev3_reference.ggp was produced by afm_scan_avg_profile_rev3.py (the
pre-port standalone script, since retired) run on
data/TASTE_ALS_A205_Ti_Pt_flatten.txt with its original settings. That script
applied no edge-exclusion rule, so the equivalence only holds with
edge_exclusion=0 - which is exactly the point: it separates "the port moved the
code" from "the edge rule changed the result".

Run directly:
    .venv/bin/python tests/test_ggp_equivalence.py
"""
import os
import sys

import numpy as np

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from afm_analysis.config import PROJECT_ROOT
from afm_analysis.core.processing import load_afm_data, raw_data, find_groove_positions
from afm_analysis.boundary import (flatten_endpoints, average_grooves,
                                   normalize_profile, profile_metrics)

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'fixtures', 'rev3_reference.ggp')
SOURCE = os.path.join(PROJECT_ROOT, 'data', 'TASTE_ALS_A205_Ti_Pt_flatten.txt')


def build_profile(edge_exclusion):
    """Reproduce the original script's pipeline at a given edge-exclusion setting"""
    data, scan_x_size = load_afm_data(SOURCE, default_scan_size=2.0)
    raw_x, raw_y = raw_data(data, scan_x_size)

    flat_y = flatten_endpoints(raw_x, raw_y)
    flat_y = flat_y - np.polyval(np.polyfit(raw_x, flat_y, 1), raw_x)

    scan_width_nm = scan_x_size * 1000
    period_nm = scan_width_nm / max(2, int(scan_width_nm / 315.0))

    centers = find_groove_positions(raw_x, flat_y, period_nm,
                                    prominence_factor=0.01, distance_factor=0.3,
                                    edge_exclusion=edge_exclusion)
    if len(centers) > 1:
        period_nm = np.mean(np.diff(centers) * (raw_x[1] - raw_x[0]) * 1000)

    x_avg, y_avg, _, n_used = average_grooves(raw_x, flat_y, centers, period_nm,
                                              margin=0.0, n_points=2000)
    x_norm, y_norm, _ = normalize_profile(x_avg, y_avg, period_nm)
    return x_norm, y_norm, period_nm, n_used


def _as_written(values, fmt='%f'):
    """
    Round to the precision actually written to disk.

    .ggp files are written with fmt='%f', i.e. 6 decimal places, so the stored
    reference is rounded. Comparing full-precision floats against it fails by
    ~5e-7 - half a unit in the last written place - which is rounding, not a
    behaviour change. The file is what PCGrate consumes, so the file is what
    this compares.
    """
    return np.array([float(fmt % v) for v in values])


def test_matches_original_script():
    """With no edge rule, the port must reproduce the original script exactly"""
    reference = np.loadtxt(FIXTURE)  # numpy '#' header is skipped automatically
    x_norm, y_norm, _, _ = build_profile(edge_exclusion=0.0)

    assert reference.shape[0] == len(x_norm), (
        f"point count changed: {reference.shape[0]} -> {len(x_norm)}")
    assert np.array_equal(reference[:, 0], _as_written(x_norm)), \
        "x values differ from the original"
    assert np.array_equal(reference[:, 1], _as_written(y_norm)), \
        "y values differ from the original"


def test_edge_rule_reduces_x_stretch():
    """
    The edge rule should shrink, not introduce, horizontal distortion.

    normalize_profile stretches the averaged window to span exactly one period.
    A groove near the scan edge narrows the window shared by all grooves, so the
    profile gets stretched to compensate. Excluding it should bring the window
    back in line with the true period.
    """
    def stretch(edge_exclusion):
        data, scan_x_size = load_afm_data(SOURCE, default_scan_size=2.0)
        raw_x, raw_y = raw_data(data, scan_x_size)
        flat_y = flatten_endpoints(raw_x, raw_y)
        flat_y = flat_y - np.polyval(np.polyfit(raw_x, flat_y, 1), raw_x)
        dx_nm = (raw_x[1] - raw_x[0]) * 1000
        period_nm = (scan_x_size * 1000) / max(2, int(scan_x_size * 1000 / 315.0))
        centers = find_groove_positions(raw_x, flat_y, period_nm,
                                        prominence_factor=0.01, distance_factor=0.3,
                                        edge_exclusion=edge_exclusion)
        period_nm = np.mean(np.diff(centers) * dx_nm)
        half_width = int(round((period_nm / 2) / dx_nm))
        common = min(min(half_width, c, len(raw_x) - 1 - c) for c in centers)
        window_nm = 2 * common * dx_nm
        return abs(period_nm - window_nm) / period_nm

    assert stretch(0.6) < stretch(0.0), "edge rule should reduce x-axis stretch"
    assert stretch(0.6) < 0.01, "edge-excluded window should be within 1% of a period"


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
