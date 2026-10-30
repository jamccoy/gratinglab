"""
Groove orientation: a sample rotated on the stage, measured and undone.

The synthetic scans here are the fixture's grating (315 nm, 30 deg blaze, 70 deg
anti-blaze, long facet descending) drawn at a known rotation, so every test
asserts recovery of a number it was given. Two effects are pinned, both in the
direction the module docstring derives:

* uncorrected, ``tan(blaze)`` reads ``cos(theta)`` short and the period
  ``1 / cos(theta)`` long;
* corrected, both come back to the truth.

The period half is asserted at 3 deg only. At 1 deg the effect is 0.048 nm,
and the analysis measures the period from whole-pixel groove centres, which on
the finest grid used here resolves about 0.05 nm: a test there would be
asserting noise. The blaze angle resolves 1 deg comfortably.
"""
import contextlib
import io
import math
import os
import sys

import numpy as np
import pytest

from scans import real_scan

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import matplotlib  # noqa: E402
matplotlib.use('Agg')

from gratinglab.metrology.analyzer import analyze_single_file  # noqa: E402
from gratinglab.metrology.boundary import build_boundary_profile  # noqa: E402
from gratinglab.metrology.core.orientation import (  # noqa: E402
    OrientationReport, align_grooves, apply_groove_orientation,
    measure_groove_orientation)
from gratinglab.metrology.core.processing import load_scan_y_size  # noqa: E402
from gratinglab.metrology.gui.state import (  # noqa: E402
    FormState, build, summarize_result)
from gratinglab.metrology.io.metrics import write_profile_metrics  # noqa: E402
from gratinglab.metrology.settings import AnalysisSettings  # noqa: E402

PERIOD_NM = 315.0
BLAZE_DEG = 30.0
ANTIBLAZE_DEG = 70.0


def _sawtooth(u_nm):
    """The fixture's groove, heights in nm, the long facet descending."""
    cot_b = 1.0 / math.tan(math.radians(BLAZE_DEG))
    cot_a = 1.0 / math.tan(math.radians(ANTIBLAZE_DEG))
    apex = cot_b / (cot_b + cot_a)
    depth = 1.0 / (cot_b + cot_a)
    t = np.mod(-u_nm / PERIOD_NM, 1.0)
    return np.where(t < apex, t * depth / apex,
                    depth * (1.0 - t) / (1.0 - apex)) * PERIOD_NM


def _image(theta_deg, *, rows=96, cols=768, width_um=8.0, height_um=2.0,
           noise_nm=0.15, seed=1):
    """A rotated grating in metres. Positive theta moves grooves to +x down the rows."""
    x = np.linspace(0.0, width_um * 1000.0, cols)
    y = np.linspace(0.0, height_um * 1000.0, rows)
    xx, yy = np.meshgrid(x, y)
    theta = math.radians(theta_deg)
    z = _sawtooth(xx * math.cos(theta) - yy * math.sin(theta))
    z = z + np.random.default_rng(seed).normal(0.0, noise_nm, z.shape)
    return z * 1e-9


def _drifting_image(theta_deg, *, rows, cols=3072, width_um=8.0, height_um=2.0,
                    blaze_from=29.0, blaze_to=31.0, noise_nm=0.15, seed=3):
    """A rotated grating whose blaze drifts down the rows, 29 deg to 31 deg.

    The drift is the point. On a grating that is the same all the way down,
    which band a given line lands in makes almost no difference, so a scan cut
    into bands two different ways measures the same either way and a test on
    `_image` cannot see the band placement at all. A real scan varies down the
    image -- the TASTE scan's groups spread over 2 deg -- and then the
    placement shows. 2 deg over 97 lines is that scan's spread, quietened.
    """
    x = np.linspace(0.0, width_um * 1000.0, cols)
    y = np.linspace(0.0, height_um * 1000.0, rows)
    theta = math.radians(theta_deg)
    z = np.empty((rows, cols))
    for row in range(rows):
        blaze = blaze_from + (blaze_to - blaze_from) * row / max(rows - 1, 1)
        cot_b = 1.0 / math.tan(math.radians(blaze))
        cot_a = 1.0 / math.tan(math.radians(ANTIBLAZE_DEG))
        apex = cot_b / (cot_b + cot_a)
        depth = 1.0 / (cot_b + cot_a)
        t = np.mod(-(x * math.cos(theta) - y[row] * math.sin(theta)) / PERIOD_NM, 1.0)
        z[row] = np.where(t < apex, t * depth / apex,
                          depth * (1.0 - t) / (1.0 - apex)) * PERIOD_NM
    z = z + np.random.default_rng(seed).normal(0.0, noise_nm, z.shape)
    return z * 1e-9


def _write_text(path, data, width_um, height_um=None):
    """As Gwyddion writes it; ``height_um=None`` leaves the Height line out."""
    head = ["# Channel: ZSensor", f"# Width: {width_um:.3f} µm"]
    head.append(f"# Height: {height_um:.3f} µm" if height_um is not None
                else "# Note: no height recorded")
    head.append("# Value units: m")
    body = ["\t".join(f"{v:.6e}" for v in row) for row in data]
    path.write_text("\n".join(head + body) + "\n", encoding="utf-8")
    return path


def _analyse(path, mode, n_row_groups=16):
    settings = AnalysisSettings(groove_orientation=mode, n_row_groups=n_row_groups)
    with contextlib.redirect_stdout(io.StringIO()):
        return analyze_single_file(str(path), show_plots=False, settings=settings)


def _apparent_blaze(theta_deg):
    return math.degrees(math.atan(math.tan(math.radians(BLAZE_DEG))
                                  * math.cos(math.radians(theta_deg))))


# ── Measuring ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("theta", [0.0, 1.0, 3.0, -2.0])
def test_the_angle_is_recovered_with_its_sign(theta):
    measured = measure_groove_orientation(_image(theta), 8.0, 2.0, PERIOD_NM)
    assert measured.angle_deg == pytest.approx(theta, abs=0.01)
    assert 0 < measured.stderr_deg < 0.01
    assert abs(measured.angle_deg - theta) < 5 * measured.stderr_deg + 1e-3


def test_the_slow_axis_pixel_is_used_not_assumed_square():
    """
    The same 3 deg grating on a scan half as tall, with the same row count.

    Each row is then half as far from the next, so the groove moves half as
    many pixels per row. A measurement that took the pixel to be square
    would report about half the angle.
    """
    short = _image(3.0, height_um=1.0)
    assert measure_groove_orientation(short, 8.0, 1.0, PERIOD_NM).angle_deg == \
        pytest.approx(3.0, abs=0.01)
    wrong = measure_groove_orientation(short, 8.0, 2.0, PERIOD_NM).angle_deg
    assert wrong == pytest.approx(math.degrees(math.atan(
        math.tan(math.radians(3.0)) / 2)), abs=0.01)


def test_a_scan_with_too_few_rows_is_reported_not_guessed():
    reason = measure_groove_orientation(_image(1.0, rows=6), 8.0, 2.0, PERIOD_NM)
    assert isinstance(reason, str) and "bands" in reason


# ── Correcting ───────────────────────────────────────────────────────────────

def test_alignment_crops_only_what_the_shift_slides_off():
    rows, cols = 96, 768
    data = _image(3.0, rows=rows, cols=cols)
    aligned, width, cropped = align_grooves(data, 8.0, 2.0, 3.0)
    pixel_x = 8000.0 / (cols - 1)
    pixel_y = 2000.0 / (rows - 1)
    travel = math.tan(math.radians(3.0)) * (rows - 1) * pixel_y / pixel_x
    assert cropped == pytest.approx(travel, abs=2)
    assert aligned.shape == (rows, cols - cropped)
    assert width == pytest.approx(
        pixel_x * math.cos(math.radians(3.0)) * (cols - cropped - 1) / 1000.0)


def test_alignment_removes_the_drift():
    """After correcting, the grooves no longer move from row to row."""
    data = _image(3.0)
    aligned, width, _ = align_grooves(data, 8.0, 2.0, 3.0)
    again = measure_groove_orientation(aligned, width, 2.0, PERIOD_NM)
    assert again.angle_deg == pytest.approx(0.0, abs=0.01)


def test_an_absurd_angle_is_refused():
    with pytest.raises(ValueError, match="columns"):
        align_grooves(_image(0.0), 8.0, 2.0, 70.0)


# ── Together with the region mask and scar repair ───────────────────────────

def test_a_masked_rotated_scan_is_still_measured_and_corrected(tmp_path):
    """
    A column mask leaves a gap in every band. The bands are searched across it
    and the grooves beside it dropped, so the angle is still found, and the
    masked analysis still recovers the true facet.
    """
    path = _write_text(tmp_path / "rot3.txt", _image(3.0, cols=3072), 8.0, 2.0)
    settings = AnalysisSettings(groove_orientation='always', n_row_groups=16,
                                mask_cols=((1200, 1500),))
    with contextlib.redirect_stdout(io.StringIO()):
        result = analyze_single_file(str(path), show_plots=False,
                                     settings=settings)
    report = result['groove_orientation']
    assert report['applied']
    assert report['angle_deg'] == pytest.approx(3.0, abs=0.02)
    assert result['mean_angle'] == pytest.approx(BLAZE_DEG, abs=0.02)
    assert result['mask_rejected_grooves'] > 0


def test_scars_are_repaired_before_the_rows_are_shifted(tmp_path):
    """Repair first: the correction would otherwise carry each scar sideways."""
    data = _image(3.0, cols=3072)
    for line in (20, 21, 60):
        data[line, 400:1400] += 60e-9
    path = _write_text(tmp_path / "rot3_scars.txt", data, 8.0, 2.0)
    settings = AnalysisSettings(groove_orientation='always', n_row_groups=16,
                                scar_repair='interpolate')
    with contextlib.redirect_stdout(io.StringIO()):
        result = analyze_single_file(str(path), show_plots=False,
                                     settings=settings)
    assert result['scar_lines'] == 3
    assert result['groove_orientation']['angle_deg'] == pytest.approx(3.0, abs=0.02)
    assert result['mean_angle'] == pytest.approx(BLAZE_DEG, abs=0.02)


# ── The analysis, end to end ─────────────────────────────────────────────────

@pytest.fixture(scope="module")
def rotated_scans(tmp_path_factory):
    """Each angle written once, at a fine grid so the period resolves."""
    folder = tmp_path_factory.mktemp("rotated")
    scans = {}
    for theta in (0.0, 1.0, 3.0):
        data = _image(theta, cols=3072)
        scans[theta] = _write_text(folder / f"rot_{theta:g}.txt", data, 8.0, 2.0)
    return scans


@pytest.mark.parametrize("theta", [0.0, 1.0, 3.0])
def test_corrected_blaze_and_period_recover_the_truth(rotated_scans, theta):
    result = _analyse(rotated_scans[theta], 'always')
    report = result['groove_orientation']
    assert report['applied']
    assert report['angle_deg'] == pytest.approx(theta, abs=0.01)
    assert result['mean_angle'] == pytest.approx(BLAZE_DEG, abs=0.01)
    assert result['period_nm'] == pytest.approx(PERIOD_NM, abs=0.1)


@pytest.mark.parametrize("theta", [1.0, 3.0])
def test_uncorrected_blaze_reads_shallow_by_cos_theta(rotated_scans, theta):
    result = _analyse(rotated_scans[theta], 'off')
    assert result['groove_orientation']['applied'] is False
    assert result['mean_angle'] == pytest.approx(_apparent_blaze(theta), abs=0.002)
    # The error is real, and the correction is what removes it.
    corrected = _analyse(rotated_scans[theta], 'always')['mean_angle']
    assert corrected - result['mean_angle'] == pytest.approx(
        BLAZE_DEG - _apparent_blaze(theta), abs=0.004)


def test_uncorrected_period_reads_long_by_one_over_cos_theta(rotated_scans):
    """0.43 nm at 3 deg, against about 0.05 nm of period resolution."""
    expected = PERIOD_NM / math.cos(math.radians(3.0))
    uncorrected = _analyse(rotated_scans[3.0], 'off')['period_nm']
    corrected = _analyse(rotated_scans[3.0], 'always')['period_nm']
    assert uncorrected == pytest.approx(expected, abs=0.05)
    assert uncorrected - corrected == pytest.approx(expected - PERIOD_NM, abs=0.1)


def test_auto_corrects_a_real_rotation_and_leaves_none_alone(rotated_scans):
    assert _analyse(rotated_scans[3.0], 'auto')['groove_orientation']['applied']
    level = _analyse(rotated_scans[0.0], 'auto')['groove_orientation']
    assert not level['applied']
    assert "threshold" in level['reason']


def test_detect_reports_the_angle_and_changes_nothing(rotated_scans):
    detected = _analyse(rotated_scans[3.0], 'detect')
    untouched = _analyse(rotated_scans[3.0], 'off')
    assert detected['groove_orientation']['angle_deg'] == pytest.approx(3.0, abs=0.01)
    assert not detected['groove_orientation']['applied']
    assert detected['mean_angle'] == untouched['mean_angle']
    assert detected['period_nm'] == untouched['period_nm']


def test_a_scan_read_the_other_way_up_measures_the_same(tmp_path):
    """
    Same grating, slow axis reversed: the correction must give the same answer.

    Not a hypothetical. A Nanoscope `.spm` and the Gwyddion text export of the
    same scan come out of the two readers row for row in opposite order, so one
    grating reaches this pipeline both ways up, and the only thing that should
    flip is the sign of the measured rotation.

    97 lines in 16 bands leaves one line over, which is the case that used to
    hurt: the last band took it, so reversing the rows moved every boundary.
    Here that reads 29.990 deg one way up and 30.009 the other; on the TASTE
    scan, whose bands spread four times as far, it was worth 0.067 deg, against
    the 0.005 deg its two readings now agree to.
    """
    data = _drifting_image(1.5, rows=97)
    upright = _analyse(_write_text(tmp_path / "up.txt", data, 8.0, 2.0), 'always')
    flipped = _analyse(_write_text(tmp_path / "down.txt", data[::-1], 8.0, 2.0),
                       'always')

    assert flipped['groove_orientation']['angle_deg'] == pytest.approx(
        -upright['groove_orientation']['angle_deg'], abs=0.005), \
        "the rotation is the one thing that should change sign"
    assert flipped['mean_angle'] == pytest.approx(upright['mean_angle'], abs=0.001)
    assert flipped['period_nm'] == pytest.approx(upright['period_nm'], abs=0.01)
    assert flipped['n_grooves'] == upright['n_grooves']
    # The centred bands are also the more accurate pair: the drift averages to
    # 30 deg, and the old placement read 29.990 / 30.009 either side of it.
    assert upright['mean_angle'] == pytest.approx(30.0, abs=0.005)


def test_the_traditional_path_carries_the_report_too(rotated_scans):
    settings = AnalysisSettings(groove_orientation='always', use_row_groups=False)
    with contextlib.redirect_stdout(io.StringIO()):
        result = analyze_single_file(str(rotated_scans[3.0]), show_plots=False,
                                     settings=settings)
    assert result['groove_orientation']['applied']
    assert result['mean_angle'] == pytest.approx(BLAZE_DEG, abs=0.01)


# ── Refusing rather than guessing ────────────────────────────────────────────

def test_a_file_without_a_slow_axis_size_is_not_measured(tmp_path):
    path = _write_text(tmp_path / "no_height.txt", _image(3.0), 8.0, None)
    assert load_scan_y_size(path) is None
    result = _analyse(path, 'auto')
    report = result['groove_orientation']
    assert report['angle_deg'] is None and not report['applied']
    assert "slow-axis" in report['reason']


def test_always_refuses_when_it_cannot_measure(tmp_path):
    data = _image(3.0)
    settings = AnalysisSettings(groove_orientation='always')
    with pytest.raises(ValueError, match="slow-axis"):
        apply_groove_orientation(data, 8.0, None, settings)


def test_off_does_nothing_and_says_so():
    data = _image(3.0)
    out, width, report = apply_groove_orientation(
        data, 8.0, 2.0, AnalysisSettings(groove_orientation='off'))
    assert out is data and width == 8.0
    assert not report.measured and report.reason == "switched off"


# ── Where the number lands ───────────────────────────────────────────────────

def test_the_whole_image_average_is_no_longer_smeared():
    """
    The boundary profile averages every row, which is where a rotation hurts
    most: at 3 deg the grooves drift 105 nm, a third of a period, from the
    first row to the last, and the averaged groove loses depth to it.

    Measured on the fine grid (2.6 nm pixels, 120 per period). The row shift
    interpolates linearly, which blunts this ideal groove's knife-edge apex
    and trough by about one pixel over one period: 3% of the depth at 33
    pixels per period. That is the correction's own cost, stated in the module
    docstring, and it is a tenth of the 32% the smearing takes.
    """
    settings = AnalysisSettings(groove_orientation='always')
    level = build_boundary_profile(_image(0.0, cols=3072), 8.0, settings)
    rotated = _image(3.0, cols=3072)
    data, width, report = apply_groove_orientation(rotated, 8.0, 2.0, settings)
    corrected = build_boundary_profile(data, width, settings, orientation=report)
    smeared = build_boundary_profile(rotated, 8.0, settings)

    true_depth = level.metrics['groove_depth']
    assert corrected.metrics['groove_depth'] == pytest.approx(true_depth, rel=0.01)
    assert smeared.metrics['groove_depth'] < 0.9 * true_depth
    assert corrected.metrics['groove_rotation_deg'] == pytest.approx(3.0, abs=0.01)
    assert corrected.metrics['groove_rotation_corrected'] is True


def test_the_metrics_sidecar_states_the_rotation(tmp_path):
    settings = AnalysisSettings(groove_orientation='always')
    data, width, report = apply_groove_orientation(_image(3.0), 8.0, 2.0, settings)
    profile = build_boundary_profile(data, width, settings, orientation=report)
    text = write_profile_metrics(tmp_path / "m.txt", profile.metrics).read_text(
        encoding="utf-8")
    assert "Groove rotation: +3.0" in text and "(corrected)" in text


def test_the_results_panel_states_the_rotation(rotated_scans):
    text = summarize_result(_analyse(rotated_scans[3.0], 'always'))
    assert "Groove rotation  : +3.0" in text and "corrected" in text


def test_the_window_passes_the_mode_through():
    settings, errors = build(FormState(groove_orientation='detect'))
    assert not errors and settings.groove_orientation == 'detect'
    settings, errors = build(FormState(groove_orientation='sideways'))
    assert settings is None
    assert any(e.field == 'groove_orientation' for e in errors)


def test_the_report_summary_reads_both_ways():
    report = OrientationReport(mode='detect', applied=False, reason="detect only",
                               angle_deg=1.25, stderr_deg=0.01, n_tracks=4,
                               n_points=40)
    assert "+1.250" in report.summary() and "not corrected" in report.summary()
    missing = OrientationReport(mode='auto', applied=False, reason="no rows")
    assert missing.summary() == "groove orientation not measured: no rows"


# ── Slow-axis size, per format ───────────────────────────────────────────────

def test_text_height_line_is_read(tmp_path):
    path = _write_text(tmp_path / "h.txt", _image(0.0, rows=12, cols=40), 8.0, 0.333)
    assert load_scan_y_size(path) == pytest.approx(0.333)


def test_asylum_slow_size_comes_from_the_y_scaling(tmp_path):
    """Metres per line times the lines recorded, so an early stop is honest."""
    from test_ibw_reader import _ramp, _write_ibw

    path = _write_ibw(tmp_path / "scan.ibw", {"ZSensorRetrace": _ramp(rows=8, cols=16)},
                      dx_m=5e-9)
    assert load_scan_y_size(path) == pytest.approx(5e-9 * 8 * 1e6)

    unscaled = _write_ibw(tmp_path / "bare.ibw", {"ZSensorRetrace": _ramp(rows=8, cols=16)},
                          dx_m=1.0, note="ScanSize: 5e-06\rScanRate: 1\r")
    assert load_scan_y_size(unscaled) is None, "a size the note does not give was guessed"


def test_nanoscope_slow_size_is_scaled_to_the_lines_recorded(tmp_path):
    """
    The second Scan Size value is the size the scan was set to cover in
    ``Lines`` lines; a scan stopped early holds fewer, so the recorded height
    is that size times rows over lines.
    """
    from gratinglab.metrology.io.spm import scan_y_size_um

    bs = chr(92)
    header = "\r\n".join([
        f"{bs}*File list",
        f"{bs}*Ciao scan list",
        f"{bs}Lines: 512",
        f"{bs}Scan Size: 2000 nm",
        f"{bs}*Ciao image list",
        f"{bs}Samps/line: 512",
        f"{bs}Number of lines: 401",
        f"{bs}Line Direction: Retrace",
        f"{bs}Scan Size: 2.000 2.000 ~m",
        f'{bs}@2:Image Data: S [Height] "Height Sensor"',
        f"{bs}*File list end",
    ])
    path = tmp_path / "fake.spm"
    path.write_bytes(header.encode("latin-1"))
    assert scan_y_size_um(path) == pytest.approx(2.0 * 401 / 512)

    single = path.read_bytes().replace(b"Scan Size: 2.000 2.000", b"Scan Size: 2.000")
    path.write_bytes(single)
    assert scan_y_size_um(path) is None


def test_real_nanoscope_scan_has_a_slow_axis_and_an_angle():
    """On the group's own scan: a size from the header, and a measured angle."""
    path = real_scan('TASTE_ALS_A205_Ti_Pt.0_00003.spm')
    from gratinglab.metrology.core.processing import load_afm_data
    from gratinglab.metrology.core.image_flatten import flatten_image

    size_y = load_scan_y_size(path)
    assert size_y is not None and 0 < size_y <= 2.0
    with contextlib.redirect_stdout(io.StringIO()):
        data, size_x = load_afm_data(path)
    measured = measure_groove_orientation(flatten_image(data, 'align_rows'),
                                          size_x, size_y, PERIOD_NM)
    assert not isinstance(measured, str), measured
    assert abs(measured.angle_deg) < 10 and measured.stderr_deg < 1.0
