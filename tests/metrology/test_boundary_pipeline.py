"""
Building a boundary profile, and the panel that previews one.

`test_ggp_equivalence.py` already pins the exported numbers against the
standalone script this replaced. What is checked here is the seam introduced so
the GUI could preview it: that computing and writing stayed equivalent, and that
every control actually reaches the output.
"""
import contextlib
import io
import os
import sys
import tempfile

import numpy as np

from scans import SYNTHETIC, SYNTHETIC_ANTIBLAZE_DEG, SYNTHETIC_BLAZE_DEG, SYNTHETIC_PERIOD_NM
import pytest

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import matplotlib
matplotlib.use('Agg')

from gratinglab.metrology.boundary import build_boundary_profile                # noqa: E402
from gratinglab.metrology.core.image_flatten import flatten_image               # noqa: E402
from gratinglab.metrology.core.processing import load_afm_data                  # noqa: E402
from gratinglab.io.ggp import GGP_HEADER, write_ggp                             # noqa: E402
from gratinglab.metrology.settings import AnalysisSettings                      # noqa: E402

SOURCE = str(SYNTHETIC)

_cache = {}


def _scan():
    if 'scan' not in _cache:
        with contextlib.redirect_stdout(io.StringIO()):
            _cache['scan'] = load_afm_data(SOURCE, default_scan_size=2.0)
    return _cache['scan']


def _profile(**overrides):
    data, scan_size = _scan()
    settings = AnalysisSettings.from_config().with_(**overrides)
    if settings.image_flatten_method != 'none':
        data = flatten_image(data, settings.image_flatten_method)
    return build_boundary_profile(data, scan_size, settings)


# ── The profile itself ───────────────────────────────────────────────────────

def test_builds_a_normalised_profile():
    p = _profile()
    assert len(p.x_norm) == len(p.y_norm)
    assert np.isclose(p.x_norm[0], 0.0) and np.isclose(p.x_norm[-1], 1.0), \
        "x must span exactly one period, 0 to 1"
    assert np.all(np.isfinite(p.y_norm))


def test_endpoints_are_exactly_zero_so_periods_tile():
    """
    A step at the period boundary is the defect that makes a PCGrate result
    wrong while the file still looks fine.
    """
    p = _profile()
    assert p.y_norm[0] == 0.0
    assert p.y_norm[-1] == 0.0


def test_reports_how_many_grooves_were_used():
    p = _profile()
    assert p.n_used > 0
    assert p.n_used <= p.n_grooves
    assert 'grooves averaged' in p.summary


def test_sigma_band_accompanies_the_average():
    """The panel draws this; an empty or wrong-length band would be silent"""
    p = _profile()
    assert p.y_std_nm.shape == p.y_avg_nm.shape
    assert np.all(p.y_std_nm >= 0)


def test_no_grooves_raises_rather_than_returning_something_empty():
    """
    Detecting nothing should say so, not return an empty profile.

    Triggered with an impossible prominence rather than a wrong period: a period
    far from reality still finds grooves, because detection keys on prominence
    and spacing rather than on the estimate being right.
    """
    try:
        _profile(prominence_factor=2.0)
    except ValueError as exc:
        assert 'no grooves detected' in str(exc).lower()
        assert 'prominence' in str(exc).lower(), \
            "the message should point at what to change"
        return
    raise AssertionError("expected a ValueError when no grooves are detected")


# ── The controls the panel exposes ───────────────────────────────────────────

def test_point_count_reaches_the_output():
    assert len(_profile(ggp_n_points=500).x_norm) == 500
    assert len(_profile(ggp_n_points=2000).x_norm) == 2000


def test_smoothing_changes_the_curve():
    """A control that changes nothing is not wired up"""
    smoothed = _profile(ggp_apply_smoothing=True).y_norm
    raw = _profile(ggp_apply_smoothing=False).y_norm
    assert not np.array_equal(smoothed, raw)


def test_minimum_half_width_excludes_grooves_once_it_exceeds_their_extent():
    """
    On this scan the exclusion is all-or-nothing, which is worth knowing.

    Every surviving groove sits far enough from both edges to get the full
    half-width, so they share one extent and cross the threshold together. On
    the synthetic scan that extent is about 30 samples -- a 315 nm period at
    5.21 nm per sample, 384 columns over 2 um -- so 29 keeps all five and 30
    keeps none. The control is a guard against clipped grooves, not a gradual
    filter.

    The numbers are sampling-dependent and would change if the fixture were
    regenerated at a different width; the all-or-nothing *behaviour* is the
    claim, and the constants only pin where the step falls.
    """
    assert _profile(ggp_min_half_width=10).n_used == 5
    assert _profile(ggp_min_half_width=29).n_used == 5

    try:
        _profile(ggp_min_half_width=31)
    except ValueError as exc:
        assert 'extracted' in str(exc).lower()
        return
    raise AssertionError(
        "a minimum half-width above every groove's extent should leave nothing "
        "to average, and say so")


# ── Consistency with the command-line mode ───────────────────────────────────

def test_image_flattening_is_a_no_op_at_written_precision():
    """
    The boundary path applies image flattening for consistency with the analysis
    path. It must not change the exported file - and does not, for the same
    structural reason as everywhere else: the endpoint flattening and linear
    detrend downstream remove constant and linear terms again.

    Measured at 8e-16 in full precision. The file is written at 6 decimals, so
    that is what this compares - the file is what PCGrate consumes.
    """
    def written(method):
        y = _profile(image_flatten_method=method).y_norm
        return np.array([float('%f' % v) for v in y])

    baseline = written('none')
    for method in ('align_rows', 'plane'):
        assert np.array_equal(baseline, written(method)), \
            f"image flattening '{method}' changed the exported profile"


def test_the_panel_and_the_cli_compute_the_same_profile():
    """
    ANALYSIS_MODE = 'ggp' and the Boundary tab must not drift apart. They share
    build_boundary_profile precisely so they cannot, and this is the check that
    the workflow really does use it.
    """
    from gratinglab.metrology.workflows import run_boundary_profile_export

    # Named explicitly rather than leaning on config.GGP_SOURCE_FILE, which
    # points at a real scan this machine may not have. The claim is that the two
    # paths agree on one input, not which input config happens to name.
    with contextlib.redirect_stdout(io.StringIO()):
        via_cli = run_boundary_profile_export(filename=SOURCE)
    via_panel = _profile()

    assert np.array_equal(via_cli['x'], via_panel.x_norm)
    assert np.array_equal(via_cli['y'], via_panel.y_norm)
    assert via_cli['metrics']['groove_depth'] == \
        via_panel.metrics['groove_depth']


# ── The written file ─────────────────────────────────────────────────────────

def test_written_file_has_the_uncommented_two_line_header():
    """
    PCGrate rejects a commented header. Three different header variants were
    found among the hand-made files before this was automated.
    """
    p = _profile(ggp_n_points=100)
    with tempfile.TemporaryDirectory() as scratch:
        path = os.path.join(scratch, 'out.ggp')
        write_ggp(path, t=p.x_norm, y=p.y_norm)
        with open(path) as handle:
            lines = handle.read().splitlines()

    assert lines[0] == GGP_HEADER.splitlines()[0] == "3 0 - Polygonal type"
    assert lines[1] == "Period: 1 PSC: 1"
    assert not lines[0].startswith('#')
    assert len(lines) == 102, "two header lines plus one per point"


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


# ── The handoff to the solvers ───────────────────────────────────────────────
#
# The point of these living in one package. A .ggp round trip preserves the
# shape and loses the scale; to_problem() keeps the period this scan measured.

def test_to_problem_carries_the_measured_period():
    """The number a .ggp cannot hold, and that corpus.toml exists to restore."""
    p = _profile()
    problem = p.to_problem()
    assert problem.period == p.period_nm
    assert problem.period > 0


def test_to_problem_scales_the_depth_by_that_period():
    """y_norm is a fraction of the period, so the period is inside the physics."""
    p = _profile()
    problem = p.to_problem()
    assert problem.depth == pytest.approx(
        p.metrics['peak_to_valley'] * p.period_nm, rel=1e-9)


def test_to_problem_produces_the_same_shape_as_the_ggp_export():
    """Whichever route a profile takes, it must be the same groove."""
    p = _profile()
    profile = p.to_problem().profile
    assert np.allclose(profile.t, p.x_norm)
    assert np.allclose(profile.y, p.y_norm)


def test_to_problem_carries_a_fitted_blaze_angle_when_given_one():
    p = _profile()
    assert p.to_problem().profile.blaze_angle is None
    assert p.to_problem(blaze_angle=29.5).profile.blaze_angle == pytest.approx(29.5)


def test_to_problem_reports_the_grooves_actually_averaged():
    p = _profile()
    assert p.to_problem().n_grooves == p.n_used


def test_to_problem_does_not_invent_a_roughness():
    """y_std_nm is groove-to-groove form spread, not surface microroughness.

    Wiring one into the other would produce a confident wrong efficiency
    instead of an obviously missing one. 0.0 here means "not supplied".
    """
    p = _profile()
    assert p.to_problem().roughness == 0.0
    assert np.any(p.y_std_nm > 0), "the sample does have groove-to-groove spread"


def test_the_measured_groove_solves():
    """End to end: an AFM scan reaches an efficiency without touching disk."""
    from gratinglab.illumination import Illumination
    from gratinglab.solvers import get_solver

    problem = _profile().to_problem(coating='Au', blaze_angle=29.5)
    scan = get_solver('scalar').solve(
        problem,
        Illumination.offplane(graze=1.5, azimuth=25.0),
        np.linspace(1.0, 5.0, 5),
        quadrature_points=1024)
    assert np.all(np.isfinite(scan.efficiency))
    assert scan.efficiency.max() > 0


# ── The synthetic fixture as a control ───────────────────────────────────────
#
# These are the tests a real scan cannot support, because on a real scan nobody
# knows the true answer. The fixture is an ideal sawtooth of stated period and
# blaze angle, so the pipeline can be asked to recover them.

def test_the_measured_period_matches_the_one_that_was_built_in():
    p = _profile()
    assert p.period_nm == pytest.approx(SYNTHETIC_PERIOD_NM, rel=0.01), \
        "period is measured from groove spacing; it must recover what was drawn"


def test_the_depth_matches_the_blaze_geometry_that_was_built_in():
    """Depth and facet angle are locked together on a sharp groove."""
    from gratinglab.profiles import Blazed

    ideal = Blazed(blaze_angle=SYNTHETIC_BLAZE_DEG,
                   antiblaze_angle=SYNTHETIC_ANTIBLAZE_DEG)
    p = _profile()
    # A few percent low: the apex falls between samples at 384 columns, so
    # discretisation clips it. Real scans lose far more, and for a reason that
    # is not discretisation -- see findings.md, "The measured groove is rounded".
    assert p.metrics['peak_to_valley'] == pytest.approx(ideal.depth, rel=0.03)


def test_depth_and_facet_fit_agree_on_a_sharp_groove():
    """The control the TASTE finding lacked.

    On a sharp sawtooth the blaze angle implied by the depth and the blaze
    angle from a facet fit are the same number by construction. They disagree
    by 7.6 deg on the one real scan in evidence, and the whole weight of that
    finding rests on this case behaving differently. Without this test, a
    reader has only the assertion that a sharp groove would agree.

    Note what the tolerance is *not* expressed in. The fitted angle here has a
    scatter of ~0.1 deg, so even this 0.3 deg discretisation gap is several
    sigma -- the same "3.6 sigma" the real scan shows. Sigma counts the noise,
    not the disagreement. The absolute gap is what separates the two cases:
    0.3 deg here against 7.6 deg there.
    """
    import contextlib
    import io

    from scipy.optimize import brentq

    from gratinglab.metrology.analyzer import analyze_single_file
    from gratinglab.profiles import Blazed

    p = _profile()
    with contextlib.redirect_stdout(io.StringIO()):
        fit = analyze_single_file(SOURCE, settings=AnalysisSettings.from_config())

    implied = brentq(
        lambda a: Blazed(blaze_angle=a,
                         antiblaze_angle=SYNTHETIC_ANTIBLAZE_DEG).depth
        - p.metrics['peak_to_valley'],
        0.5, 89.0)

    assert fit['mean_angle'] == pytest.approx(SYNTHETIC_BLAZE_DEG, abs=0.5)
    assert abs(fit['mean_angle'] - implied) < 1.0, (
        f"a sharp groove must agree with itself: fitted {fit['mean_angle']:.2f} "
        f"deg against {implied:.2f} deg implied by the depth")
