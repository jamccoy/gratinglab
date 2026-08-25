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
        text = path.read_text()
        assert "Tip correction: erosion" in text
        assert "upper bounds" in text

    def test_the_sidecar_states_the_absence_too(self, tmp_path):
        from gratinglab.metrology.io.metrics import write_profile_metrics
        p = self._profile(SYNTHETIC_DILATED)
        text = write_profile_metrics(tmp_path / "m.txt", p.metrics).read_text()
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
