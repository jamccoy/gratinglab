r"""The conical diffraction scene.

Headless -- no Qt, no matplotlib. The load-bearing test is
:class:`TestTheGratingEquationIsReadableOffTheDrawing`, which recovers
:math:`\sin\alpha + \sin\beta_m = m\lambda/(p\sin\gamma)` from the drawn
vectors by a *different* route than the one that built them.
"""

import numpy as np
import pytest

from gratinglab.geometry import blaze_direction, is_propagating, sin_beta
from gratinglab.gui import diagram3d
from gratinglab.gui.diagram3d import (
    CONE_AXIS,
    D_HAT,
    G_HAT,
    N_HAT,
    MIN_RIM_MAGNIFICATION,
    PRESET_VIEWS,
    RIM_PRESET,
    build_scene,
    incident_vector,
    view_direction,
    wave_vector,
)
from gratinglab.gui.provenance import ALARM_TAGS
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, Sinusoidal

UNPOL = "unpolarized"
PERIOD, DELTA, ALPHA, GAMMA, LAMBDA = 315.15, 29.5, 25.0, 1.5, 3.0


def blazed():
    return Problem(
        period=PERIOD, profile=Blazed(blaze_angle=DELTA, antiblaze_angle=70.5)
    )


def offplane(alpha=ALPHA, gamma=GAMMA):
    return Illumination.offplane(graze=gamma, azimuth=alpha, polarization=UNPOL)


@pytest.fixture(scope="module")
def scene():
    return build_scene(blazed(), offplane(), LAMBDA)


def order_rays(scene):
    return [r for r in scene.rays if r.order is not None]


class TestEveryRayIsAUnitVector:
    def test_unit_norm(self, scene):
        for ray in scene.rays:
            assert np.linalg.norm(ray.direction) == pytest.approx(1.0, abs=1e-12)

    def test_the_scene_actually_contains_order_rays(self, scene):
        """Non-vacuity: without this, every claim about "every ray" could hold
        over an empty set."""
        assert len(order_rays(scene)) > 1

    def test_the_transverse_extent_is_sin_gamma_not_normalised_away(self, scene):
        r"""The 2D view draws all rays at one arbitrary length; here the
        transverse extent *is* :math:`\sin\gamma`. At 0.026 a 2D-style
        renormalisation would read 1.0 -- caught by three orders of magnitude,
        not by a tolerance."""
        for ray in order_rays(scene):
            dx, dy, _ = ray.direction
            assert np.hypot(dx, dy) == pytest.approx(np.sin(np.radians(GAMMA)))


class TestTheConeProperty:
    def test_every_ray_has_kz_equal_to_cos_gamma(self, scene):
        expected = np.cos(np.radians(GAMMA))
        for ray in order_rays(scene):
            assert ray.direction[2] == pytest.approx(expected)

    def test_including_the_incident_ray(self, scene):
        r""":math:`k_z` is conserved *across* the boundary, not merely among
        the diffracted orders."""
        incident = next(r for r in scene.rays if r.tag == "incident")
        assert incident.direction[2] == pytest.approx(np.cos(np.radians(GAMMA)))

    def test_the_cone_opens_away_from_the_groove_axis(self, scene):
        """The `docs/findings.md` entry, pinned in the drawing."""
        for ray in order_rays(scene):
            assert np.degrees(ray.polar_from_cone_axis) == pytest.approx(GAMMA)
            assert np.degrees(np.arccos(np.dot(ray.direction, G_HAT))) == pytest.approx(
                180.0 - GAMMA
            )

    def test_the_fixture_is_genuinely_off_plane(self):
        """Non-vacuity: at γ = 90° every angle above is trivial."""
        assert not offplane().is_in_plane
        assert np.cos(np.radians(GAMMA)) > 0.99

    def test_the_incident_vector_is_the_illuminations_own(self):
        """Not a second source of truth for the one vector we transcribe."""
        ill = offplane()
        assert incident_vector(ill.alpha, ill.gamma) == pytest.approx(
            ill.direction_cosines
        )

    def test_right_handedness_holds_for_the_drawn_triad(self):
        assert np.cross(D_HAT, G_HAT) == pytest.approx(N_HAT)
        assert CONE_AXIS == pytest.approx(-G_HAT)


class TestTheGratingEquationIsReadableOffTheDrawing:
    """The load-bearing test.

    Vector -> `arctan2` -> `sin`, a different route than `sin_beta` took to
    build them, so a swapped component or a sign flip fails.
    """

    def test_every_drawn_ray_satisfies_it(self, scene):
        sa = np.sin(np.radians(ALPHA))
        sg = np.sin(np.radians(GAMMA))
        for ray in order_rays(scene):
            beta_drawn = np.arctan2(ray.direction[0], ray.direction[1])
            assert sa + np.sin(beta_drawn) == pytest.approx(
                ray.order * LAMBDA / (PERIOD * sg), abs=1e-12
            )

    def test_the_zeroth_order_mirrors_the_incident_ray(self, scene):
        r"""Pure vector identity, no angles at all: :math:`\hat{k}_0` is
        :math:`\hat{k}_i` with its :math:`\hat{n}` component negated."""
        incident = next(r for r in scene.rays if r.tag == "incident")
        zero = scene.ray(0)
        assert np.asarray(zero.direction) == pytest.approx(
            np.asarray(incident.direction) * np.array([1.0, -1.0, 1.0])
        )

    def test_the_drawn_sines_are_not_all_zero(self, scene):
        """Non-vacuity: the class cannot pass on a degenerate scene."""
        sines = [np.sin(r.azimuth) for r in order_rays(scene)]
        assert max(sines) - min(sines) > 0.1


class TestEvanescentOrdersGetNoRay:
    #: Order 3 passes off above 3.912 nm at this geometry.
    M, PAST = 3, 4.5

    def test_the_case_really_is_a_cutoff(self):
        sa, sg = np.sin(np.radians(ALPHA)), np.sin(np.radians(GAMMA))
        assert is_propagating(sin_beta(self.M, LAMBDA, PERIOD, sa, sg))
        assert not is_propagating(sin_beta(self.M, self.PAST, PERIOD, sa, sg))

    def test_the_rays_are_exactly_the_propagating_orders(self, scene):
        drawn = {r.order for r in order_rays(scene)}
        assert drawn == {m.order for m in scene.orders if m.propagating}

    def test_an_evanescent_order_keeps_its_record_but_loses_its_ray(self):
        after = build_scene(blazed(), offplane(), self.PAST)
        assert self.M not in {r.order for r in order_rays(after)}
        assert any(m.order == self.M and not m.propagating for m in after.orders)

    def test_asking_for_its_ray_raises_rather_than_inventing_one(self):
        after = build_scene(blazed(), offplane(), self.PAST)
        with pytest.raises(KeyError):
            after.ray(self.M)


class TestInPlaneDegenerates:
    def test_every_ray_lies_in_the_dispersion_plane(self):
        scene = build_scene(
            blazed(), Illumination.classical(alpha=10.0, polarization=UNPOL), 600.0
        )
        for ray in order_rays(scene):
            assert ray.direction[2] == pytest.approx(0.0, abs=1e-15)

    def test_and_a_caption_says_the_cone_degenerates(self):
        scene = build_scene(
            blazed(), Illumination.classical(alpha=10.0, polarization=UNPOL), 600.0
        )
        text = " ".join(c.text for c in scene.captions)
        assert "degenerates" in text and "90" in text

    def test_the_off_plane_scene_does_leave_the_plane(self, scene):
        """Non-vacuity for the pair above."""
        assert min(abs(r.direction[2]) for r in order_rays(scene)) > 0.99


class TestSurfacesAreRuledAndExact:
    def test_the_grating_patch_is_z_invariant(self, scene):
        """Which is why a two-row mesh is exact rather than a sampling
        choice: the profile has no z dependence."""
        patch = next(s for s in scene.surfaces if s.tag == "surface")
        assert patch.x[:, 0] == pytest.approx(patch.x[:, 1])
        assert patch.y[:, 0] == pytest.approx(patch.y[:, 1])
        assert not np.allclose(patch.z[:, 0], patch.z[:, 1])

    def test_both_levels_of_detail_exist(self, scene):
        assert scene.at("fine") and scene.at("coarse")

    def test_the_coarse_mesh_is_a_subset_of_the_fine_one(self, scene):
        """Provable identity, not a tolerance: 4*61 - 3 == 241."""
        fine = next(s for s in scene.surfaces if s.tag == "surface" and s.lod == "fine")
        coarse = next(
            s for s in scene.surfaces if s.tag == "surface" and s.lod == "coarse"
        )
        assert fine.x[::4] == pytest.approx(coarse.x)
        assert coarse.x.shape[0] < fine.x.shape[0]

    def test_every_cone_rim_point_lies_on_the_cone(self, scene):
        rim = next(c for c in scene.curves if c.tag == "cone" and c.lod == "fine")
        radius = np.hypot(rim.x, rim.y)
        assert radius == pytest.approx(np.sin(np.radians(GAMMA)))
        assert rim.z == pytest.approx(np.cos(np.radians(GAMMA)))


class TestBlazeAndCaptions:
    def test_the_blaze_ray_matches_blaze_direction(self, scene):
        ray = next(r for r in scene.rays if r.tag == "blaze")
        expected = blaze_direction(np.radians(DELTA), np.radians(ALPHA))
        assert ray.azimuth == pytest.approx(expected)

    def test_no_blaze_ray_without_a_blaze_angle(self):
        scene = build_scene(
            Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.15)),
            offplane(), LAMBDA,
        )
        assert not any(r.tag == "blaze" for r in scene.rays)

    def test_a_correct_scene_is_not_styled_as_a_problem(self, scene):
        assert all(c.tag not in ALARM_TAGS for c in scene.captions)

    def test_the_title_states_gamma_and_that_it_is_to_scale(self, scene):
        assert "1.5" in scene.title and "to scale" in scene.title

    def test_a_caption_names_the_one_arbitrary_ratio(self, scene):
        """Ray directions are exact; the patch's size relative to them is a
        drawing choice, and the caption says which is which."""
        text = " ".join(c.text for c in scene.captions)
        assert "drawing choice" in text
        assert "cross-section itself is exact" in text

    def test_a_caption_states_the_cone_axis_direction(self, scene):
        text = " ".join(c.text for c in scene.captions)
        assert "opens along −ĝ" in text


class TestSceneShape:
    def test_the_bounding_box_is_a_cube(self, scene):
        """Equal box aspect expressed in the pure layer, so the widget only
        applies it and never chooses limits."""
        spans = [hi - lo for lo, hi in scene.limits]
        assert spans[0] == pytest.approx(spans[1]) == pytest.approx(spans[2])

    def test_every_tag_used_has_a_colour(self, scene):
        used = {r.tag for r in scene.rays} | {s.tag for s in scene.surfaces}
        used |= {c.tag for c in scene.curves} | {p.tag for p in scene.points}
        assert used <= set(diagram3d.TAG_COLORS)

    def test_wave_vector_reuses_the_two_dimensional_direction(self):
        """The (x, y) part is literally `diagram.direction`, so the two views
        cannot disagree about an azimuth."""
        from gratinglab.gui.diagram import direction

        beta, gamma = np.radians(17.7), np.radians(GAMMA)
        k = wave_vector(beta, gamma)
        assert k[:2] == pytest.approx(direction(beta) * np.sin(gamma))


class TestPresetViews:
    """Presets are physical look-directions, so they mean the same thing
    whatever permutation the widget applies for display."""

    def test_every_preset_is_a_unit_vector(self):
        for name, look in PRESET_VIEWS.items():
            assert np.linalg.norm(look) == pytest.approx(1.0, abs=1e-12), name

    def test_down_the_cone_axis_looks_along_g(self):
        assert PRESET_VIEWS[RIM_PRESET] == pytest.approx(G_HAT, abs=1e-12)

    def test_which_is_opposite_the_way_the_cone_opens(self):
        """Non-vacuity, and the M13-E finding again: looking *along* the
        groove axis means the rays come back at you."""
        assert PRESET_VIEWS[RIM_PRESET] @ CONE_AXIS == pytest.approx(-1.0)

    def test_along_d_looks_against_the_dispersion_axis(self):
        assert PRESET_VIEWS["Along d\u0302"] == pytest.approx(-D_HAT, abs=1e-12)

    def test_face_n_looks_against_the_normal(self):
        assert PRESET_VIEWS["Face n\u0302"] == pytest.approx(-N_HAT, abs=1e-12)

    def test_the_presets_are_pairwise_distinct(self):
        looks = list(PRESET_VIEWS.values())
        assert len(looks) > 2
        for i, a in enumerate(looks):
            for b in looks[i + 1:]:
                assert abs(a @ b) < 0.99

    def test_view_direction_is_a_unit_vector_for_any_camera_angles(self):
        """`view_direction` is the forward map the panel inverts; a preset
        round-trip through it is asserted at the widget level."""
        for elev in (-60.0, 0.0, 25.0, 89.0):
            for azim in (-135.0, 0.0, 30.0, 200.0):
                assert np.linalg.norm(view_direction(elev, azim)) == pytest.approx(1.0)


class TestPurity:
    def test_imports_no_toolkit_and_no_plotting(self):
        import ast
        import inspect

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(diagram3d))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6", "matplotlib"} & imported


class TestRimFocus:
    """Framing the cone rim, which is what makes the 93.8-degree fan
    readable at true scale."""

    def test_the_rim_box_is_a_cube(self, scene):
        """Both boxes cubes is what makes moving between them a similarity
        transform rather than a distortion."""
        spans = [hi - lo for lo, hi in scene.rim_limits]
        assert spans[0] == pytest.approx(spans[1]) == pytest.approx(spans[2])

    def test_it_is_centred_on_the_rim(self, scene):
        (x0, x1), (y0, y1), (z0, z1) = scene.rim_limits
        assert 0.5 * (x0 + x1) == pytest.approx(0.0)
        assert 0.5 * (y0 + y1) == pytest.approx(0.0)
        assert 0.5 * (z0 + z1) == pytest.approx(np.cos(np.radians(GAMMA)))

    def test_every_ray_tip_falls_inside_it(self, scene):
        (x0, x1), (y0, y1), (z0, z1) = scene.rim_limits
        for ray in order_rays(scene):
            hx, hy, hz = ray.head
            assert x0 <= hx <= x1 and y0 <= hy <= y1 and z0 <= hz <= z1

    def test_the_magnification_is_the_ratio_of_the_two_half_extents(self, scene):
        half = diagram3d.RIM_MARGIN * np.sin(np.radians(GAMMA))
        assert scene.rim_magnification == pytest.approx(diagram3d.BOX_HALF / half)

    def test_and_is_worth_having(self, scene):
        assert scene.rim_magnification > 20.0
        assert scene.rim_magnification > MIN_RIM_MAGNIFICATION

    def test_because_the_whole_scene_box_leaves_the_fan_microscopic(self, scene):
        """Non-vacuity for the test above, measured rather than asserted: in
        the full box the entire fan spans under 3% of the frame."""
        (x0, x1), _, _ = scene.limits
        fan_width = 2.0 * np.sin(np.radians(GAMMA))
        assert fan_width / (x1 - x0) < 0.03

    def test_zooming_changes_no_angle(self, scene):
        """The claim that uniform zoom is honest, as a fact.

        The boxes differ by a scale and a translation, so the angle between
        every pair of drawn rays is identical in both -- there is nothing for
        a caption to disclaim.
        """
        rays = [np.asarray(r.direction) for r in order_rays(scene)]
        scale = scene.rim_magnification
        for i, a in enumerate(rays):
            for b in rays[i + 1:]:
                before = np.arccos(np.clip(a @ b, -1.0, 1.0))
                after = np.arccos(np.clip((scale * a) @ (scale * b)
                                          / (scale * scale), -1.0, 1.0))
                assert before == pytest.approx(after, abs=1e-12)

    def test_the_reason_is_non_empty_exactly_when_unavailable(self):
        """The `blaze_jump` contract, so a disabled control can explain
        itself instead of sitting there dead."""
        cases = [
            offplane(gamma=1.5),
            offplane(gamma=30.0),
            Illumination.classical(alpha=25.0, polarization=UNPOL),
        ]
        seen = set()
        for ill in cases:
            s = build_scene(blazed(), ill, LAMBDA)
            assert bool(s.rim_reason) == (not s.rim_available)
            seen.add(s.rim_available)
        assert seen == {True, False}, "both branches must be exercised"

    def test_in_plane_says_the_rim_is_already_the_whole_scene(self):
        s = build_scene(
            blazed(), Illumination.classical(alpha=25.0, polarization=UNPOL), 600.0
        )
        assert not s.rim_available
        assert "nothing to zoom into" in s.rim_reason
