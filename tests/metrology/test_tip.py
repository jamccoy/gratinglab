"""
Tip convolution: the erosion recovers what it can and admits what it cannot.

The validation lever is the dilated fixture: the same synthetic grating as
``synthetic_blazed_scan.txt`` with the *identical* noise field, imaged through
a known worn tip (R = 20 nm, 18 deg) by ``core.tip.dilate``. Because the true
surface is committed next to the image of it, recovery can be asserted
pointwise -- which no real scan allows -- and the failure mode can be pinned
just as hard: the trough wedge the sphere cannot enter stays lost, the depth
stays short, and the certainty map is what says so.

Dilation and erosion are one module but the tests do not lean on the
round-trip alone: image >= surface, surface <= reconstruction <= image, exact
redilation, and the spike-images-the-inverted-tip closed form each pin one
operator against geometry rather than against the other operator.
"""
import contextlib
import io
import os
import sys

import numpy as np
import pytest

from scans import (SYNTHETIC, SYNTHETIC_DILATED, SYNTHETIC_ANTIBLAZE_DEG,
                   SYNTHETIC_BLAZE_DEG, SYNTHETIC_TIP_HALF_ANGLE_DEG,
                   SYNTHETIC_TIP_RADIUS_NM)

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import matplotlib
matplotlib.use('Agg')

from gratinglab.metrology.boundary import build_boundary_profile               # noqa: E402
from gratinglab.metrology.core.image_flatten import flatten_image              # noqa: E402
from gratinglab.metrology.core.processing import load_afm_data                 # noqa: E402
from gratinglab.metrology.core.tip import (apply_tip_correction, dilate,       # noqa: E402
                                           erode, tip_cross_section)
from gratinglab.metrology.settings import AnalysisSettings                     # noqa: E402

SCAN_X_UM = 2.0
TIP = dict(radius_nm=SYNTHETIC_TIP_RADIUS_NM,
           half_angle_deg=SYNTHETIC_TIP_HALF_ANGLE_DEG)

_cache = {}


def _fixture(path):
    if path not in _cache:
        with contextlib.redirect_stdout(io.StringIO()):
            _cache[path] = load_afm_data(str(path), default_scan_size=SCAN_X_UM)
    return _cache[path]


def _triangle():
    """A clean, noiseless triangle wave in metres: 100 nm deep, 315 nm period."""
    x = np.linspace(0.0, 2000.0, 384, endpoint=False)
    y = np.abs((x % 315.0) - 157.5) / 157.5 * 100.0
    return np.tile(y, (4, 1)) * 1e-9


# ── The tip itself ───────────────────────────────────────────────────────────

class TestTheTipShape:
    def test_apex_at_zero_symmetric_and_monotone(self):
        t = tip_cross_section(5.0, 20.0, 18.0, 150.0)
        k = len(t) // 2
        assert t[k] == 0.0
        assert np.array_equal(t, t[::-1])
        assert (np.diff(t[k:]) > 0).all()

    def test_sphere_meets_cone_with_a_continuous_slope(self):
        # A slope jump at the tangency would spike the second difference far
        # above the sphere's own curvature floor, 1/(R sin^3(theta)) * dx^2 =
        # 1.7e-4 nm at this sampling. The tangency construction is what this
        # pins: a kink there would imprint a phantom edge on every image.
        dx = 0.01
        t = tip_cross_section(dx, 20.0, 18.0, 30.0)
        curvature_floor = dx**2 / (20.0 * np.sin(np.radians(18.0))**3)
        assert np.abs(np.diff(t, 2)).max() < 1.5 * curvature_floor

    def test_it_reaches_the_stated_height(self):
        t = tip_cross_section(5.0, 20.0, 18.0, 150.0)
        assert t.max() >= 150.0

    def test_bad_parameters_are_refused(self):
        with pytest.raises(ValueError, match="radius"):
            tip_cross_section(5.0, 0.0, 18.0, 100.0)
        with pytest.raises(ValueError, match="half angle"):
            tip_cross_section(5.0, 20.0, 90.0, 100.0)

    def test_a_sub_pixel_tip_still_has_a_window(self):
        t = tip_cross_section(5.0, 0.5, 18.0, 0.2)
        assert len(t) == 3


# ── The operators, against geometry ──────────────────────────────────────────

class TestTheOperators:
    def test_a_tip_can_hide_a_trough_but_never_dig_one(self):
        s = _triangle()
        img = dilate(s, SCAN_X_UM, **TIP)
        assert (img >= s - 1e-15).all()
        assert img.max() == pytest.approx(s.max(), abs=1e-15)  # peaks are touched

    def test_the_reconstruction_sits_between_truth_and_image(self):
        s = _triangle()
        img = dilate(s, SCAN_X_UM, **TIP)
        rec = erode(img, SCAN_X_UM, **TIP).data
        assert (rec >= s - 1e-15).all()
        assert (rec <= img + 1e-15).all()

    def test_certain_pixels_are_recovered_to_machine_precision(self):
        s = _triangle()
        c = erode(dilate(s, SCAN_X_UM, **TIP), SCAN_X_UM, **TIP)
        assert 0.5 < c.certain_fraction < 1.0
        assert np.abs((c.data - s)[c.certain]).max() < 1e-20  # metres

    def test_redilating_the_reconstruction_reproduces_the_image_exactly(self):
        # The reconstruction is the least upper bound *consistent with the
        # image*: push it back through the tip and the image must come back.
        s = _triangle()
        img = dilate(s, SCAN_X_UM, **TIP)
        rec = erode(img, SCAN_X_UM, **TIP).data
        assert np.array_equal(dilate(rec, SCAN_X_UM, **TIP), img)

    def test_a_spike_images_the_inverted_tip(self):
        # The classic closed form, and the one test of `dilate` that involves
        # no morphology at all: a delta-function surface is narrower than any
        # tip, so its image *is* the tip, upside down.
        spike = np.zeros((2, 384))
        spike[:, 200] = 120.0e-9
        img = dilate(spike, SCAN_X_UM, **TIP)
        dx_nm = SCAN_X_UM * 1000.0 / 383
        t = tip_cross_section(dx_nm, TIP['radius_nm'], TIP['half_angle_deg'], 120.0)
        k = len(t) // 2
        expected = np.maximum(120.0 - t, 0.0) * 1e-9
        assert np.allclose(img[0, 200 - k:200 + k + 1], expected, atol=1e-18)

    def test_uncertainty_concentrates_at_the_trough_corners(self):
        s = _triangle()
        c = erode(dilate(s, SCAN_X_UM, **TIP), SCAN_X_UM, **TIP)
        troughs = np.flatnonzero(s[0] == s[0].min())
        mids = (troughs[:-1] + np.diff(troughs) // 4)  # on the open facet
        assert not c.certain[0, troughs].any(), \
            "the sphere cannot enter the trough wedge; those pixels are bounds"
        assert c.certain[0, mids].all(), \
            "mid-facet pixels are apex contacts and must be certain"


# ── The fixtures: recovery with the truth committed next door ────────────────

class TestTheDilatedFixture:
    def test_the_dilation_is_the_only_difference(self):
        # Identical noise by construction, so the pointwise gap *is* the tip.
        sharp, _ = _fixture(SYNTHETIC)
        dil, _ = _fixture(SYNTHETIC_DILATED)
        gap_nm = (dil - sharp) * 1e9
        assert gap_nm.min() > -1e-3, "dilation never digs"
        assert gap_nm.max() > 20.0, "and here it visibly buried the troughs"

    def test_erosion_recovers_the_certain_pixels_to_the_noise_floor(self):
        sharp, _ = _fixture(SYNTHETIC)
        dil, scan = _fixture(SYNTHETIC_DILATED)
        c = erode(dil, scan, **TIP)
        err_nm = np.abs((c.data - sharp) * 1e9)
        # Calibrated on the fixture: uncorrected mean error 9.1 nm; corrected,
        # over certain pixels, 0.18 nm mean / 1.0 nm max -- the residual is
        # the 0.15 nm noise field interacting with the min-filter, not the tip.
        assert np.abs((dil - sharp) * 1e9).mean() > 5.0
        assert err_nm[c.certain].mean() < 0.5
        assert err_nm[c.certain].max() < 2.0
        assert 0.75 < c.certain_fraction < 0.95


class TestThePipeline:
    def _profile(self, path, **overrides):
        data, scan = _fixture(path)
        settings = AnalysisSettings.from_config().with_(**overrides)
        data = flatten_image(data, settings.image_flatten_method)
        return build_boundary_profile(data, scan, settings)

    CORRECTED = dict(tip_correction='erosion',
                     tip_radius_nm=SYNTHETIC_TIP_RADIUS_NM,
                     tip_half_angle_deg=SYNTHETIC_TIP_HALF_ANGLE_DEG)

    def test_the_metrics_record_the_correction_and_its_reach(self):
        p = self._profile(SYNTHETIC_DILATED, **self.CORRECTED)
        m = p.metrics
        assert m['tip_correction'] == 'erosion'
        assert m['tip_radius_nm'] == SYNTHETIC_TIP_RADIUS_NM
        assert m['tip_half_angle_deg'] == SYNTHETIC_TIP_HALF_ANGLE_DEG
        assert 0.75 < m['tip_certain_fraction'] < 0.95

    def test_uncorrected_metrics_say_so(self):
        p = self._profile(SYNTHETIC_DILATED)
        assert p.metrics['tip_correction'] == 'none'
        assert 'tip_certain_fraction' not in p.metrics

    def test_the_depth_deficit_is_not_resurrected_and_the_warning_exists(self):
        """Erosion is a bound, not a resurrection -- pinned, not hoped.

        The sphere never enters the trough wedge (geometric standoff ~9 nm at
        this radius in this 80 deg wedge), so the corrected profile is still
        shallow against the sharp truth. If this test ever fails in the
        deeper-than-ideal direction, the erosion has started inventing
        surface, which is the one thing it must never do. The pipeline's
        honesty about the residual lives in `tip_certain_fraction`, asserted
        above -- that number, not the corrected depth, is the warning the
        rounded-groove finding asked for.
        """
        from gratinglab.profiles import Blazed
        ideal = Blazed(blaze_angle=SYNTHETIC_BLAZE_DEG,
                       antiblaze_angle=SYNTHETIC_ANTIBLAZE_DEG)
        sharp = self._profile(SYNTHETIC).metrics['peak_to_valley']
        buried = self._profile(SYNTHETIC_DILATED).metrics['peak_to_valley']
        corrected = self._profile(SYNTHETIC_DILATED,
                                  **self.CORRECTED).metrics['peak_to_valley']
        assert buried < sharp - 0.02, "the fixture visibly lost depth"
        assert corrected < ideal.depth - 0.02, "and erosion does not restore it"
        assert corrected == pytest.approx(buried, abs=0.01)

    def test_the_blaze_angle_survives_the_corrected_path(self):
        # Facet trimming already keeps the fit off the rounded corners, so the
        # angle is insensitive to this tip either way; what this pins is that
        # the corrected path changes nothing it should not.
        from gratinglab.metrology.analyzer import analyze_single_file
        settings = AnalysisSettings.from_config().with_(**self.CORRECTED)
        with contextlib.redirect_stdout(io.StringIO()):
            fit = analyze_single_file(str(SYNTHETIC_DILATED), show_plots=False,
                                      settings=settings)
        assert fit['mean_angle'] == pytest.approx(SYNTHETIC_BLAZE_DEG, abs=0.5)

    def test_the_sidecar_states_the_correction(self, tmp_path):
        from gratinglab.metrology.io.metrics import write_profile_metrics
        p = self._profile(SYNTHETIC_DILATED, **self.CORRECTED)
        path = write_profile_metrics(tmp_path / "metrics.txt", p.metrics)
        text = path.read_text(encoding="utf-8")
        assert "Tip correction: erosion" in text
        assert "upper bounds" in text

    def test_the_sidecar_states_the_absence_too(self, tmp_path):
        from gratinglab.metrology.io.metrics import write_profile_metrics
        p = self._profile(SYNTHETIC_DILATED)
        text = write_profile_metrics(tmp_path / "m.txt", p.metrics).read_text(encoding="utf-8")
        assert "Tip correction: none" in text


# ── The real scan: the finding, pinned ───────────────────────────────────────

class TestTheRealScan:
    def test_nominal_tip_erosion_does_not_move_the_taste_depth(self):
        """`findings.md`, "A nominal tip does not explain the rounded groove".

        The measured TASTE surface is already reachable by the nominal probe
        (R ~ 1-2 nm, 18 deg), so erosion recovers no depth: the rounding is
        not a sharp tip's convolution artefact. If this ever fails in the
        deeper direction, either the scan on this machine changed or the
        erosion has started inventing surface; both need investigating, not a
        tolerance bump.
        """
        from scans import real_scan
        path = real_scan("TASTE_ALS_A205_Ti_Pt_flatten.txt")

        with contextlib.redirect_stdout(io.StringIO()):
            data, scan = load_afm_data(str(path), default_scan_size=2.0)
        base = AnalysisSettings.from_config()
        flat = flatten_image(data, base.image_flatten_method)
        with contextlib.redirect_stdout(io.StringIO()):
            plain = build_boundary_profile(flat, scan, base)
            eroded = build_boundary_profile(
                flat, scan, base.with_(tip_correction='erosion',
                                       tip_radius_nm=2.0,
                                       tip_half_angle_deg=18.0))
        assert abs(eroded.metrics['peak_to_valley']
                   - plain.metrics['peak_to_valley']) < 1e-3
        assert eroded.metrics['tip_certain_fraction'] > 0.9


# ── Settings and dispatch ────────────────────────────────────────────────────

class TestTheSettings:
    def test_none_is_a_no_op_that_returns_the_same_array(self):
        data = _triangle()
        settings = AnalysisSettings.from_config()
        out, correction = apply_tip_correction(data, SCAN_X_UM, settings)
        assert out is data
        assert correction is None

    def test_erosion_dispatches_and_reports(self):
        settings = AnalysisSettings.from_config().with_(
            tip_correction='erosion', tip_radius_nm=20.0,
            tip_half_angle_deg=18.0)
        out, correction = apply_tip_correction(_triangle(), SCAN_X_UM, settings)
        assert correction is not None
        assert "erosion" in correction.summary
        assert out is correction.data

    def test_validation_catches_the_nonsense(self):
        base = AnalysisSettings.from_config()
        fields = [f for f, _ in base.with_(tip_correction='sandpaper').validate()]
        assert 'tip_correction' in fields
        bad = base.with_(tip_correction='erosion', tip_radius_nm=-1.0,
                         tip_half_angle_deg=95.0)
        fields = [f for f, _ in bad.validate()]
        assert 'tip_radius_nm' in fields and 'tip_half_angle_deg' in fields

    def test_tip_values_are_not_policed_when_the_correction_is_off(self):
        base = AnalysisSettings.from_config().with_(tip_radius_nm=-1.0)
        assert 'tip_radius_nm' not in [f for f, _ in base.validate()]


# ── An asymmetric tip: a pyramid seen along its cantilever ───────────────────

ASYM = dict(radius_nm=2.0, half_angle_deg=15.0, half_angle_right_deg=25.0)


def _wall(rise_deg=70.0, fall_deg=30.0, rows=3, cols=384):
    """A sawtooth whose steep wall RISES towards +x, in metres.

    Going right, the surface climbs a ``rise_deg`` wall and then descends a
    ``fall_deg`` facet. The wall faces -x, so it is the tip's +x flank that has
    to follow it.
    """
    period = 315.0
    cot_r = 1.0 / np.tan(np.radians(rise_deg))
    cot_f = 1.0 / np.tan(np.radians(fall_deg))
    depth = period / (cot_r + cot_f)
    t = np.linspace(0.0, 2000.0, cols, endpoint=False) % period
    run = depth * cot_r
    y = np.where(t <= run, t / cot_r, depth - (t - run) / cot_f)
    return np.tile(y, (rows, 1)) * 1e-9


class TestAnAsymmetricTip:
    def test_each_side_has_its_own_flank(self):
        dx = 0.5
        t = tip_cross_section(dx, 2.0, 15.0, 200.0, half_angle_right_deg=25.0)
        k = len(t) // 2
        # Far from the apex each side is a straight line of slope cot(angle)
        # per unit lateral distance; -x is the 15 deg side, +x the 25 deg one.
        left = (t[k - 80] - t[k - 81]) / dx
        right = (t[k + 81] - t[k + 80]) / dx
        assert left == pytest.approx(-1.0 / np.tan(np.radians(15.0)), rel=1e-9)
        assert right == pytest.approx(1.0 / np.tan(np.radians(25.0)), rel=1e-9)
        assert t[k] == 0.0

    def test_swapping_the_angles_mirrors_the_tip(self):
        a = tip_cross_section(5.0, 2.0, 15.0, 150.0, half_angle_right_deg=25.0)
        b = tip_cross_section(5.0, 2.0, 25.0, 150.0, half_angle_right_deg=15.0)
        assert np.array_equal(a, b[::-1])

    def test_omitting_the_right_angle_is_the_old_symmetric_cone(self):
        a = tip_cross_section(5.0, 20.0, 18.0, 150.0)
        b = tip_cross_section(5.0, 20.0, 18.0, 150.0, half_angle_right_deg=18.0)
        assert np.array_equal(a, b)

    def test_a_spike_images_the_tip_reflected(self):
        # The image of a delta is the tip upside down AND left-right reversed:
        # with the apex at 200 + j the spike sits at offset -j under the tip.
        # Symmetric tests cannot see the reversal; this one exists for it.
        spike = np.zeros((2, 384))
        spike[:, 200] = 120.0e-9
        img = dilate(spike, SCAN_X_UM, **ASYM)
        dx_nm = SCAN_X_UM * 1000.0 / 383
        t = tip_cross_section(dx_nm, 2.0, 15.0, 120.0, half_angle_right_deg=25.0)
        k = len(t) // 2
        expected = np.maximum(120.0 - t[::-1], 0.0) * 1e-9
        assert np.allclose(img[0, 200 - k:200 + k + 1], expected, atol=1e-18)

    def test_the_reconstruction_still_sits_between_truth_and_image(self):
        # Necessary and not sufficient: eroding by the UNREFLECTED tip passes
        # this too, as a looser bound. The three tests after it are the ones
        # that catch that mistake, checked by putting it back.
        surface = _wall()
        img = dilate(surface, SCAN_X_UM, **ASYM)
        rec = erode(img, SCAN_X_UM, **ASYM).data
        assert (rec >= surface - 1e-15).all()
        assert (rec <= img + 1e-15).all()

    def test_redilating_the_reconstruction_reproduces_the_image(self):
        img = dilate(_wall(), SCAN_X_UM, **ASYM)
        rec = erode(img, SCAN_X_UM, **ASYM).data
        assert np.allclose(dilate(rec, SCAN_X_UM, **ASYM), img, atol=1e-15, rtol=0)

    def test_certain_pixels_are_the_surface(self):
        surface = _wall()
        img = dilate(surface, SCAN_X_UM, **ASYM)
        tip = erode(img, SCAN_X_UM, **ASYM)
        assert tip.certain.any()
        assert np.allclose(tip.data[tip.certain], surface[tip.certain],
                           atol=1e-15, rtol=0)

    def test_which_face_meets_the_steep_wall_decides_whether_it_is_recovered(self):
        # The ScanAsyst-Air case. A 70 deg wall rising to +x is followed by the
        # tip's +x flank: the 15 deg front face (a 75 deg flank) clears it and
        # the 25 deg back face (65 deg) cannot, so the same tip, rotated, either
        # recovers the anti-blaze facet or loses it.
        surface = _wall(rise_deg=70.0)
        slope = np.diff(surface[0])
        on_wall = np.zeros(surface.shape[1], dtype=bool)
        on_wall[1:-1] = (slope[:-1] > 0) & (slope[1:] > 0)
        idx = np.nonzero(on_wall)[0]
        walls = np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)

        def certain_per_wall(left, right):
            tip = dict(radius_nm=2.0, half_angle_deg=left, half_angle_right_deg=right)
            img = dilate(surface, SCAN_X_UM, **tip)
            certain = erode(img, SCAN_X_UM, **tip).certain[0]
            return [int(certain[w].sum()) for w in walls]

        # 75 deg flank: nearly every wall pixel, measured 8 or 9 of 10.
        assert all(n >= 0.7 * len(w)
                   for n, w in zip(certain_per_wall(25.0, 15.0), walls))
        # 65 deg flank: the tip rides the top corner and nothing below it, so
        # at most that one pixel per wall. Measured 0 or 1, by sampling phase.
        assert all(n <= 1 for n in certain_per_wall(15.0, 25.0))

    def test_the_summary_names_both_flanks(self):
        tip = erode(dilate(_wall(), SCAN_X_UM, **ASYM), SCAN_X_UM, **ASYM)
        assert "15 deg (-x) / 25 deg (+x)" in tip.summary
        sym = erode(_triangle(), SCAN_X_UM, radius_nm=2.0, half_angle_deg=18.0)
        assert "half angle = 18 deg" in sym.summary


class TestThePresets:
    def test_scanasyst_air_matches_bruker(self):
        from gratinglab.metrology.core.tip import TIP_PRESETS
        across = TIP_PRESETS['scanasyst-air-across']
        assert (across.radius_nm, across.half_angle_deg,
                across.half_angle_right_deg) == (2.0, 17.5, 17.5)
        minus = TIP_PRESETS['scanasyst-air-along-front-minus']
        plus = TIP_PRESETS['scanasyst-air-along-front-plus']
        assert (minus.half_angle_deg, minus.half_angle_right_deg) == (15.0, 25.0)
        assert (plus.half_angle_deg, plus.half_angle_right_deg) == (25.0, 15.0)

    def test_every_preset_runs_through_the_pipeline(self):
        from gratinglab.metrology.core.tip import TIP_PRESETS
        base = AnalysisSettings.from_config().with_(tip_correction='erosion')
        for preset in TIP_PRESETS.values():
            settings = base.with_(tip_radius_nm=preset.radius_nm,
                                  tip_half_angle_deg=preset.half_angle_deg,
                                  tip_half_angle_right_deg=preset.half_angle_right_deg)
            assert not settings.validate()
            _, correction = apply_tip_correction(_triangle(), SCAN_X_UM, settings)
            assert correction.half_angle_right_deg == preset.half_angle_right_deg

    def test_a_bad_right_angle_is_refused(self):
        bad = AnalysisSettings.from_config().with_(
            tip_correction='erosion', tip_half_angle_right_deg=0.0)
        assert 'tip_half_angle_right_deg' in [f for f, _ in bad.validate()]

    def test_the_records_carry_the_right_flank(self, tmp_path):
        from gratinglab.metrology.io.metrics import write_profile_metrics
        data = flatten_image(_fixture(SYNTHETIC_DILATED)[0], 'align_rows')
        settings = AnalysisSettings.from_config().with_(
            tip_correction='erosion', tip_radius_nm=2.0,
            tip_half_angle_deg=15.0, tip_half_angle_right_deg=25.0)
        profile = build_boundary_profile(data, SCAN_X_UM, settings)
        assert profile.metrics['tip_half_angle_right_deg'] == 25.0
        path = tmp_path / "m.txt"
        write_profile_metrics(path, profile.metrics)
        assert "15 deg (-x) / 25 deg (+x)" in path.read_text(encoding="utf-8")
