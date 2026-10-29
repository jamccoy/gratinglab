"""
Tracing a groove off an SEM cross-section, and fitting a trapezoid to it.

Every image here is rendered from a :class:`Trapezoidal` whose shape is known
exactly (``xsec_synth.render``), so each test can state the right answer in
advance, and the suite needs none of the group's micrographs, which are not
in the repository and must not be.

The tolerances are set from a measured sweep rather than by hope. On a clean
render at the OGRE scale (2.233 nm pixels, a 163.3 nm period) the edge sits
within 0.1 nm of the truth per wall and the depth within 0.11 nm. What moves
the answer is physics the method cannot undo, and those tests pin the size of
it: a bright rim widens each wall by about a nanometre, and the beam rounds
every corner so the fitted radius reads high, by roughly
``R_fit^2 = R^2 + (2.5 sigma)^2``.
"""
import numpy as np
import pytest

from gratinglab.metrology.boundary import CrossSection, trace_cross_section
from gratinglab.metrology.io.sem import read_sem
from gratinglab.profiles import FromProfileData, Trapezoidal
from test_sem_reader import BASE_FIELDS, with_databar, write_tiff, zeiss_tag
from xsec_synth import render

PERIOD = 163.3
PIXEL = 2.233

#: The etched OGRE groove of docs/findings.md: 68 nm deep, a trench 0.70 of
#: the period at its mouth, a 72 degree sidewall, and a fillet at each corner.
OGRE = Trapezoidal(
    depth_fraction=68.0 / PERIOD, duty_cycle=0.30, sidewall_angle=72.0,
    corner_radius=6.0 / PERIOD,
)


def trace(truth=OGRE, **kwargs):
    render_kw = {k: kwargs.pop(k) for k in list(kwargs) if k in RENDER_KEYS}
    return trace_cross_section(render(truth, **render_kw), pixel_size_nm=PIXEL, **kwargs)


RENDER_KEYS = {
    "noise", "rim", "rim_nm", "beam_px", "slope", "height_scale", "seed",
    "width_px", "height_px", "surface_row", "phase_nm",
}


@pytest.fixture(scope="module")
def clean():
    return trace()


@pytest.fixture(scope="module")
def noisy():
    return trace(noise=6.0, seed=1)


# ── The shape comes back ─────────────────────────────────────────────────────


class TestRecoversAKnownTrapezoid:
    def test_depth(self, clean):
        assert clean.fit.depth_nm == pytest.approx(68.0, abs=0.2)

    def test_sidewall_angle(self, clean):
        assert clean.fit.sidewall_angle_deg == pytest.approx(72.0, abs=0.2)

    def test_crest_and_base_widths(self, clean):
        assert clean.fit.ridge_top_width == pytest.approx(0.30, abs=0.002)
        assert clean.fit.ridge_base_width == pytest.approx(
            OGRE.ridge_base_width, abs=0.002)

    def test_trench_widths_are_what_a_micrograph_is_read_for(self, clean):
        # findings.md reads 0.70 -> 0.42 off xsec_050; this truth is 0.70 -> 0.43.
        assert clean.fit.trench_width_mouth == pytest.approx(0.70, abs=0.002)
        assert clean.fit.trench_width_floor == pytest.approx(
            1 - OGRE.ridge_base_width, abs=0.002)

    def test_corner_radius_reads_high_by_the_beam_and_no_more(self, clean):
        # A 0.7 px beam alone rounds a sharp corner to about 4 nm, and 6 nm
        # in quadrature with that is 7.2. The fit must not come out sharper
        # than the truth: the beam can only add rounding.
        assert 6.0 <= clean.fit.corner_radius_nm <= 8.5

    def test_period_from_the_ridge_lattice(self, clean):
        assert clean.period_source == "measured"
        assert clean.period_nm == pytest.approx(PERIOD, rel=5e-4)

    def test_every_whole_groove_is_used(self, clean):
        # 640 px of 2.233 nm hold 8.75 periods; 8 ridges have a whole period
        # about them.
        assert clean.fit.n_grooves == 8
        assert len(clean.groove_centres_nm) == 8
        assert np.all(np.diff(clean.groove_centres_nm) == pytest.approx(PERIOD, rel=2e-3))

    def test_a_vertical_wall_is_found_vertical(self):
        resist = Trapezoidal(depth_fraction=223.5 / 169.5, duty_cycle=0.48,
                             sidewall_angle=90.0)
        xs = trace_cross_section(
            render(resist, period_nm=169.5, pixel_nm=2.484, height_px=180,
                   surface_row=40, noise=4.0),
            pixel_size_nm=2.484,
        )
        assert xs.fit.sidewall_angle_deg > 89.5
        assert xs.fit.depth_nm == pytest.approx(223.5, abs=0.5)
        assert xs.fit.ridge_top_width == pytest.approx(0.48, abs=0.003)

    def test_a_shallow_wall(self):
        shallow = Trapezoidal(depth_fraction=40.0 / PERIOD, duty_cycle=0.35,
                              sidewall_angle=45.0, corner_radius=10.0 / PERIOD)
        xs = trace(shallow, noise=4.0)
        assert xs.fit.sidewall_angle_deg == pytest.approx(45.0, abs=0.5)
        assert xs.fit.depth_nm == pytest.approx(40.0, abs=0.3)


class TestFitQuality:
    def test_the_trapezoid_describes_its_own_render(self, clean):
        assert clean.fit.rms_residual_nm < 0.6
        assert clean.fit.mean_profile_rms_nm < 0.2

    def test_a_shape_the_trapezoid_cannot_describe_shows_in_the_averaged_groove(
            self, clean):
        # A sinusoid is not a trapezoid. The averaged groove keeps the shape the
        # fit could not, and its distance from the fitted trapezoid says so.
        # Not by much: a trapezoid with large fillets is within 0.7 nm of a
        # 60 nm sinusoid everywhere, so the test is against the misfit of a
        # true trapezoid traced the same way, which is ten times smaller.
        class Sine(Trapezoidal):
            def height(self, t):
                t = np.asarray(t, dtype=np.float64)
                return 0.5 * self.depth * (1 - np.cos(2 * np.pi * t))

        sine = Sine(depth_fraction=60.0 / PERIOD, duty_cycle=0.3, sidewall_angle=60.0)
        xs = trace(sine)
        assert xs.fit.mean_profile_max_nm > 0.4
        assert xs.fit.mean_profile_rms_nm > 3 * clean.fit.mean_profile_rms_nm
        # ...while the averaged groove itself is still the sinusoid.
        u = np.linspace(0.05, 0.95, 50)
        want = 0.5 * 60.0 / PERIOD * (1 - np.cos(2 * np.pi * u))
        got = np.interp(u, xs.u, xs.y)
        assert np.max(np.abs(got - want)) * PERIOD < 0.6


# ── Uncertainties ────────────────────────────────────────────────────────────


class TestUncertainties:
    def test_errors_are_groove_to_groove(self, noisy):
        assert noisy.fit.uncertainty == "groove-to-groove"
        for key, values in noisy.fit.per_groove.items():
            assert len(values) == noisy.fit.n_grooves, key

    def test_errors_match_the_scatter_between_noise_realisations(self, noisy):
        # The reported standard error should describe how far another image of
        # the same grating would move the answer. Four more realisations.
        depths, angles = [], []
        for seed in range(2, 6):
            xs = trace(noise=6.0, seed=seed)
            depths.append(xs.fit.depth_nm)
            angles.append(xs.fit.sidewall_angle_deg)
        scatter_depth = np.std(depths, ddof=1)
        scatter_angle = np.std(angles, ddof=1)
        assert noisy.fit.depth_nm_err == pytest.approx(scatter_depth, rel=2.5)
        assert noisy.fit.sidewall_angle_deg_err == pytest.approx(scatter_angle, rel=2.5)
        assert noisy.fit.depth_nm_err > 0
        assert noisy.fit.sidewall_angle_deg_err > 0

    def test_noisy_answer_is_within_its_errors_and_the_known_edge_bias(self, noisy):
        f = noisy.fit
        assert abs(f.depth_nm - 68.0) < 4 * f.depth_nm_err + 0.15
        assert abs(f.sidewall_angle_deg - 72.0) < 4 * f.sidewall_angle_deg_err + 0.1
        assert abs(f.ridge_top_width - 0.30) * PERIOD < 4 * f.ridge_top_width_err * PERIOD + 0.15

    def test_one_groove_falls_back_to_the_covariance_and_says_so(self):
        xs = trace(noise=4.0, width_px=110, period_nm=PERIOD)
        assert xs.fit.n_grooves == 1
        assert xs.fit.uncertainty == "covariance"
        assert 0 < xs.fit.depth_nm_err < 1.0
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.5)


# ── What an SEM image does to an edge ────────────────────────────────────────


class TestEdgeContrast:
    def test_a_bright_rim_widens_each_wall_by_about_a_nanometre(self):
        # The rim here is brighter above the bulk than the bulk is above the
        # vacuum, and 4 nm deep: a strong InLens edge. The first crossing keeps
        # it from tracing the rim's inner flank, but the beam spreads the rim
        # into the vacuum and the edge moves out. Pinned, not hidden.
        xs = trace(noise=6.0, rim=80.0, seed=3)
        per_wall = 0.5 * (xs.fit.ridge_top_width - 0.30) * PERIOD
        assert 0.3 < per_wall < 1.5
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.6)
        assert xs.fit.sidewall_angle_deg == pytest.approx(72.0, abs=0.3)

    def test_the_bulk_darker_than_the_rim_is_still_material(self):
        # Material barely brighter than vacuum, with a rim far brighter than
        # either: the bulk falls below the level, and only the first crossing
        # and the fill keep the trench walls from being traced twice.
        img = render(OGRE, vacuum=40.0, material=70.0, rim=120.0, noise=3.0)
        xs = trace_cross_section(img, pixel_size_nm=PIXEL)
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=1.0)
        assert xs.fit.sidewall_angle_deg == pytest.approx(72.0, abs=1.0)

    def test_an_isolated_charging_speck_is_not_the_surface(self):
        img = render(OGRE, noise=3.0).copy()
        img[10, 100:102] = 250  # a bright speck well above the crest
        xs = trace_cross_section(img, pixel_size_nm=PIXEL)
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)
        assert xs.height_nm.max() - np.percentile(xs.height_nm, 90) < 2.0

    def test_material_darker_than_vacuum(self):
        img = 255 - render(OGRE, noise=4.0).astype(np.int64)
        xs = trace_cross_section(img, pixel_size_nm=PIXEL)
        assert xs.material_level < xs.vacuum_level
        assert xs.vacuum_level > xs.level > xs.material_level
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)

    def test_material_above_the_edge(self):
        img = render(OGRE, noise=4.0)[::-1]
        xs = trace_cross_section(img, pixel_size_nm=PIXEL, material_below=False)
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)
        assert xs.fit.sidewall_angle_deg == pytest.approx(72.0, abs=0.3)

    def test_a_tilted_cleave_edge_is_levelled(self):
        xs = trace(noise=4.0, slope=0.03)
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)
        assert xs.fit.sidewall_angle_deg == pytest.approx(72.0, abs=0.4)

    def test_too_little_contrast_is_refused(self):
        img = render(OGRE, vacuum=60.0, material=62.0, noise=4.0)
        with pytest.raises(ValueError, match="vacuum noise"):
            trace_cross_section(img, pixel_size_nm=PIXEL)

    def test_an_explicit_level_must_lie_between_the_two(self):
        with pytest.raises(ValueError, match="between"):
            trace_cross_section(render(OGRE), pixel_size_nm=PIXEL, level=200.0)

    def test_an_explicit_level(self):
        xs = trace_cross_section(render(OGRE), pixel_size_nm=PIXEL, level=75.0)
        assert xs.level == 75.0
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)


# ── Foreshortening is an input, never a guess ────────────────────────────────


class TestForeshortening:
    COS30 = float(np.cos(np.radians(30.0)))

    def test_a_foreshortened_image_is_refused_until_the_axis_is_named(self):
        img = render(OGRE, height_scale=self.COS30, noise=4.0)
        with pytest.raises(ValueError, match="which image axis"):
            trace_cross_section(img, pixel_size_nm=PIXEL, foreshortening=self.COS30)

    def test_vertical_foreshortening_is_undone_on_the_heights(self):
        img = render(OGRE, height_scale=self.COS30, noise=4.0)
        xs = trace_cross_section(img, pixel_size_nm=PIXEL, foreshortening=self.COS30,
                                 foreshortened="vertical")
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)
        assert xs.fit.sidewall_angle_deg == pytest.approx(72.0, abs=0.3)
        assert xs.pixel_y_nm == pytest.approx(PIXEL / self.COS30)
        assert xs.pixel_x_nm == PIXEL

    def test_naming_the_wrong_axis_gives_a_visibly_wrong_period(self):
        img = render(OGRE, height_scale=self.COS30, noise=4.0)
        xs = trace_cross_section(img, pixel_size_nm=PIXEL, foreshortening=self.COS30,
                                 foreshortened="horizontal")
        assert xs.period_nm == pytest.approx(PERIOD / self.COS30, rel=2e-3)

    def test_neither_is_an_explicit_choice(self):
        img = render(OGRE, height_scale=self.COS30, noise=4.0)
        xs = trace_cross_section(img, pixel_size_nm=PIXEL, foreshortening=self.COS30,
                                 foreshortened="neither")
        assert xs.fit.depth_nm == pytest.approx(68.0 * self.COS30, abs=0.3)

    def test_the_micrographs_own_tilt_is_read(self, tmp_path):
        fields = dict(BASE_FIELDS)
        fields["AP_STAGE_AT_T"] = "Stage at T =  30.0 °"
        img = render(OGRE, height_scale=self.COS30, noise=4.0, width_px=640)
        path = write_tiff(tmp_path / "tilted.tif", img, zeiss_tag(fields))
        sem = read_sem(path)
        assert sem.foreshortening == pytest.approx(self.COS30)
        with pytest.raises(ValueError, match="which image axis"):
            trace_cross_section(sem)
        xs = trace_cross_section(sem, foreshortened="vertical")
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)


# ── From a SmartSEM file, and into the package's conventions ─────────────────


@pytest.fixture(scope="module")
def from_file(tmp_path_factory):
    img = render(OGRE, noise=5.0, width_px=1024, height_px=200, surface_row=60)
    # SmartSEM paints its data bar over the foot of the frame.
    img = with_databar(img, rows=30)
    path = write_tiff(tmp_path_factory.mktemp("sem") / "xsec.tif", img,
                      zeiss_tag(BASE_FIELDS))
    sem = read_sem(path)
    assert sem.databar_rows == 30
    return trace_cross_section(sem)


class TestEndToEnd:
    def test_scale_comes_from_the_instrument(self, from_file):
        assert from_file.pixel_x_nm == pytest.approx(PIXEL, rel=1e-3)
        assert from_file.period_nm == pytest.approx(PERIOD, rel=1e-3)

    def test_shape_comes_from_the_file(self, from_file):
        assert from_file.fit.n_grooves >= 12
        assert from_file.fit.depth_nm == pytest.approx(68.0, abs=0.3)
        assert from_file.fit.sidewall_angle_deg == pytest.approx(72.0, abs=0.3)

    def test_a_region_crops_before_tracing(self):
        img = render(OGRE, noise=4.0, height_px=220, surface_row=110)
        img[:40] = 250  # something bright at the top of the frame, not vacuum
        with pytest.raises(ValueError):
            trace_cross_section(img, pixel_size_nm=PIXEL)
        xs = trace_cross_section(img, pixel_size_nm=PIXEL, region=(60, 220, 0, 640))
        assert xs.fit.depth_nm == pytest.approx(68.0, abs=0.3)

    def test_too_narrow_a_region_is_refused(self):
        with pytest.raises(ValueError):
            trace_cross_section(render(OGRE, width_px=40), pixel_size_nm=PIXEL)

    def test_an_array_needs_a_pixel_size(self):
        with pytest.raises(ValueError, match="pixel_size_nm"):
            trace_cross_section(render(OGRE))


def wall_widths(profile):
    """How much of the period the rising and falling walls take, 20% to 80%."""
    t = np.linspace(0.0, 1.0, 4000, endpoint=False)
    h = profile.height(t)
    band = (h > 0.2 * profile.depth) & (h < 0.8 * profile.depth)
    slope = np.gradient(h, t)
    step = t[1] - t[0]
    return (float(np.sum(band & (slope > 0)) * step),
            float(np.sum(band & (slope < 0)) * step))


class TestToProfile:
    def test_is_a_from_profile_data_normalised_to_the_period(self, clean):
        profile = clean.to_profile()
        assert isinstance(profile, FromProfileData)
        assert min(profile.t) >= 0.0 and max(profile.t) <= 1.0
        assert min(profile.y) == 0.0
        assert profile.depth * clean.period_nm == pytest.approx(68.0, abs=0.4)
        assert clean.depth_nm == pytest.approx(profile.depth * clean.period_nm)

    def test_matches_the_true_groove_pointwise(self, clean):
        # Away from the corners, where the beam rounds the truth.
        profile = clean.to_profile()
        t = np.linspace(0.0, 1.0, 400, endpoint=False)
        got = profile.height(t) * PERIOD
        want = OGRE.height(1.0 - t) * PERIOD  # t runs against the image columns
        corners = np.zeros_like(t, dtype=bool)
        for x, _ in OGRE._corners():
            corners |= np.abs(np.mod(1.0 - t - x + 0.5, 1.0) - 0.5) < 0.04
        assert np.max(np.abs(got - want)[~corners]) < 1.0

    def test_t_runs_against_the_image_columns_by_default(self):
        # An asymmetric groove: the right wall shallower than the left. The
        # Trapezoidal cannot hold that, but the averaged groove keeps it.
        class Lopsided(Trapezoidal):
            def height(self, t):
                t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
                up = np.clip((t - 0.15) / 0.05, 0, 1)
                down = np.clip((0.75 - t) / 0.25, 0, 1)
                return self.depth * np.minimum(up, down)

        shape = Lopsided(depth_fraction=50.0 / PERIOD, duty_cycle=0.3, sidewall_angle=60.0)
        along = trace(shape)
        against = trace_cross_section(render(shape), pixel_size_nm=PIXEL,
                                      image_x_along_dispersion=False)
        # The steep wall rises in u. With the columns taken to run against
        # +d-hat, t is u and the steep wall rises in t; with them running
        # along it (the default), t = 1 - u and the steep wall falls.
        assert wall_widths(against.to_profile()) == pytest.approx((0.03, 0.15), abs=0.015)
        assert wall_widths(along.to_profile()) == pytest.approx((0.15, 0.03), abs=0.015)

    def test_a_vertical_wall_stays_single_valued(self):
        resist = Trapezoidal(depth_fraction=1.0, duty_cycle=0.48, sidewall_angle=90.0)
        xs = trace_cross_section(
            render(resist, height_px=150, surface_row=30, noise=4.0),
            pixel_size_nm=PIXEL,
        )
        assert xs.to_profile().is_single_valued()

    def test_result_is_frozen(self, clean):
        assert isinstance(clean, CrossSection)
        with pytest.raises(AttributeError):
            clean.period_nm = 1.0
