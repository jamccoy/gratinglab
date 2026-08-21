"""
Image flattening: the 2-D stage, before rows are averaged.

The property worth protecting here is a surprising one: the affine methods must
give *identical* blaze angles. That is what makes defaulting to `align_rows`
free, and if the pipeline order ever changes such that it stops holding, these
tests are the warning.
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

from gratinglab.metrology.config import PROJECT_ROOT                      # noqa: E402
from gratinglab.metrology.core.image_flatten import (                     # noqa: E402
    VALID_IMAGE_FLATTEN_METHODS, flatten_image, row_offset_spread)

# The scan with the worst row-offset spread in the dataset, 2.74 nm - if any
# file would show a difference, it is this one.
WORST = os.path.join(PROJECT_ROOT, 'data', '20250820_280C_00004.txt')

_cache = {}


def _image():
    if 'image' not in _cache:
        from gratinglab.metrology.core.processing import load_afm_data
        with contextlib.redirect_stdout(io.StringIO()):
            _cache['image'] = load_afm_data(WORST, default_scan_size=2.0)
    return _cache['image']


# ── The methods themselves ───────────────────────────────────────────────────

def test_every_method_preserves_shape_and_stays_finite():
    data, _ = _image()
    for method in VALID_IMAGE_FLATTEN_METHODS:
        out = flatten_image(data, method)
        assert out.shape == data.shape, method
        assert np.all(np.isfinite(out)), f"{method} produced non-finite values"


def test_none_is_the_identity():
    data, _ = _image()
    assert np.array_equal(flatten_image(data, 'none'), data)


def test_align_rows_removes_the_row_offsets():
    """That is the whole job: every scan line on a common level."""
    data, _ = _image()
    before = row_offset_spread(data)
    after = row_offset_spread(flatten_image(data, 'align_rows'))
    assert before > 1e-12, "this fixture was supposed to have row offsets"
    assert after < 1e-15, f"row-offset spread still {after*1e9:.4f} nm"


def test_plane_removes_a_tilt():
    """A deliberately tilted image should come back level"""
    data, _ = _image()
    rows, cols = data.shape
    y, x = np.mgrid[0:rows, 0:cols]
    tilted = data + 1e-9 * (0.05 * x + 0.03 * y)

    def residual_tilt(z):
        basis = np.column_stack([x.ravel(), y.ravel(), np.ones(z.size)])
        coefficients, *_ = np.linalg.lstsq(basis, z.ravel(), rcond=None)
        return abs(coefficients[0]) + abs(coefficients[1])

    assert residual_tilt(flatten_image(tilted, 'plane')) < \
        residual_tilt(tilted) * 1e-6


def test_the_input_is_not_modified():
    """The Import tab keeps the raw image to show a before/after comparison"""
    data, _ = _image()
    original = data.copy()
    for method in VALID_IMAGE_FLATTEN_METHODS:
        flatten_image(data, method)
    assert np.array_equal(data, original)


def test_unknown_method_lists_the_registry():
    data, _ = _image()
    try:
        flatten_image(data, 'gwyddion_magic')
    except ValueError as exc:
        assert 'align_rows' in str(exc), "the error should list what is available"
        return
    raise AssertionError("expected a ValueError for an unknown method")


def test_a_one_dimensional_array_is_refused():
    try:
        flatten_image(np.arange(10.0), 'plane')
    except ValueError as exc:
        assert '2-D' in str(exc) or '2D' in str(exc)
        return
    raise AssertionError("expected a ValueError for a 1-D input")


# ── The property that makes the default free ─────────────────────────────────

def test_affine_methods_give_identical_blaze_angles():
    """
    none, plane and align_rows must agree to within floating point.

    Every profile-flattening method fits at least a first-order polynomial to the
    averaged profile, so any constant or linear term removed from the image is
    removed again downstream. This is why `align_rows` could become the default
    without moving a single stored number.

    Not bit-for-bit: subtracting row medians perturbs the last bits and the
    averaging that follows is not associative, so the results differ by around
    5e-15 degrees against a sigma near 2. The tolerance below is far tighter than
    any real effect and far looser than exact equality, which does fail.

    If this ever fails, the pipeline order has changed - most likely profile
    flattening was disabled or moved - and the defaults need revisiting rather
    than the test relaxing.
    """
    from gratinglab.metrology.analyzer import analyze_single_file
    from gratinglab.metrology.settings import AnalysisSettings

    base = AnalysisSettings.from_config()
    angles = {}
    for method in ('none', 'plane', 'align_rows'):
        with contextlib.redirect_stdout(io.StringIO()):
            result = analyze_single_file(
                WORST, show_plots=False,
                settings=base.with_(image_flatten_method=method))
        angles[method] = result['mean_angle']

    spread = max(angles.values()) - min(angles.values())
    assert spread < 1e-9, (
        f"affine image flattening changed the answer by {spread:.2e} deg: "
        f"{angles}. Profile flattening should have removed these terms again.")


def test_profile_flattening_does_change_the_answer():
    """
    The contrast that makes the point.

    Image flattening is free; profile flattening is not. About 0.49 degrees
    across the four methods on this scan - comparable to the differences between
    samples this software exists to detect.
    """
    from gratinglab.metrology.analyzer import analyze_single_file
    from gratinglab.metrology.settings import AnalysisSettings

    base = AnalysisSettings.from_config()
    angles = {}
    for method in ('linear', 'polynomial', 'groove_peaks', 'level_grooves'):
        with contextlib.redirect_stdout(io.StringIO()):
            result = analyze_single_file(
                WORST, show_plots=False, settings=base.with_(flatten_method=method))
        angles[method] = result['mean_angle']

    spread = max(angles.values()) - min(angles.values())
    assert spread > 0.1, (
        f"profile flattening methods differed by only {spread:.4f} deg; "
        f"measured at ~0.49. If this collapsed, the methods may no longer be "
        f"doing different things: {angles}")


def test_default_is_align_rows_and_costs_nothing():
    """The default in config.py, and the reason it was safe to change"""
    from gratinglab.metrology.settings import AnalysisSettings
    assert AnalysisSettings.from_config().image_flatten_method == 'align_rows'


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
