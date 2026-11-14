"""The painted 3D view's camera. Headless: no Qt, no matplotlib.

The load-bearing claim is that every framing the view offers (orbit, roll,
zoom about a point, pan, the rim slider) keeps angles true, because the view
calls itself "to scale" and the matplotlib one refused to zoom to keep that
promise.
"""

from __future__ import annotations

import ast
import inspect
from dataclasses import replace

import numpy as np
import pytest

from gratinglab.gui import camera3d, diagram3d
from gratinglab.gui.camera3d import Camera, basis, fit, orbit_for, project, zoom_about

GAMMA = np.radians(1.5)
RECT = (0.0, 0.0, 900.0, 640.0)
CLOUD = np.array([[0, 0, 0], [0.4, 0, 0], [0, 0.4, 0], [0, 0, 1.0],
                  [0.25, -0.02, 0.25], [-0.25, -0.02, -0.25]])


def screen_angle(frame, a, b, c):
    """The angle at b between a and c, measured on screen."""
    (pa, pb, pc), _ = project(frame, [a, b, c])
    u, v = pa - pb, pc - pb
    return float(np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v))))


class TestBasis:
    @pytest.mark.parametrize("yaw,pitch,roll", [(0, 0, 0), (26.6, 20.4, 0), (200, -60, 33)])
    def test_is_orthonormal_and_right_handed(self, yaw, pitch, roll):
        r, u, f = basis(yaw, pitch, roll)
        m = np.array([r, u, f])
        assert m @ m.T == pytest.approx(np.eye(3), abs=1e-12)
        # right x up points TOWARD the viewer, i.e. against fwd.
        assert np.cross(r, u) == pytest.approx(-f, abs=1e-12)

    def test_keeps_the_grating_normal_up(self):
        """With no roll the normal n-hat never points down the screen."""
        for yaw in range(0, 360, 30):
            for pitch in (-80, -20, 0, 20, 80):
                _r, u, _f = basis(yaw, pitch)
                assert u[1] > 0

    @pytest.mark.parametrize("name", list(diagram3d.PRESET_VIEWS))
    def test_every_preset_round_trips(self, name):
        want = diagram3d.PRESET_VIEWS[name]
        yaw, pitch = orbit_for(want)
        _r, _u, fwd = basis(yaw, pitch)
        # The face-on preset stops PITCH_LIMIT short of the pole.
        tol = 1 - np.cos(np.radians(90 - camera3d.PITCH_LIMIT)) + 1e-9
        assert fwd @ want == pytest.approx(1.0, abs=tol)

    def test_the_app_and_the_package_agree_on_oblique(self):
        """The app's default orbit is (26.6, 20.4); the package's Oblique is a
        direction. They are the same view."""
        yaw, pitch = orbit_for(diagram3d.PRESET_VIEWS["Oblique"])
        assert (yaw, pitch) == pytest.approx((26.6, 20.4), abs=0.1)


class TestAnglesStayTrue:
    """Orthographic, and every reframing a similarity transform."""

    def _angles(self, camera):
        frame = fit(camera, CLOUD, GAMMA, RECT)
        return [screen_angle(frame, CLOUD[i], CLOUD[0], CLOUD[j])
                for i, j in ((1, 2), (1, 3), (2, 3), (4, 5))]

    def test_zoom_pan_and_rim_change_no_angle(self):
        base = Camera(yaw=40, pitch=25)
        frame = fit(base, CLOUD, GAMMA, RECT)
        moved = zoom_about(base, frame, 7.0, (123.0, 456.0)).panned(-50, 80)
        moved = replace(moved, rim=0.6)
        assert self._angles(moved) == pytest.approx(self._angles(base), abs=1e-9)

    def test_roll_changes_no_angle(self):
        base = Camera(yaw=40, pitch=25)
        assert self._angles(base.rolled(73)) == pytest.approx(self._angles(base), abs=1e-9)

    def test_but_an_orbit_does(self):
        """Non-vacuity: a projection does change angles out of the screen
        plane, so the tests above would be empty if this one passed too."""
        base = Camera(yaw=40, pitch=25)
        assert self._angles(base.orbited(30, 0)) != pytest.approx(self._angles(base), abs=1e-3)

    def test_down_the_cone_every_azimuth_is_true(self):
        """Looking along the cone axis, each order's screen angle from n-hat
        is its true azimuth: the readable fan the preset exists for."""
        yaw, pitch = orbit_for(diagram3d.PRESET_VIEWS["Down the cone axis"])
        frame = fit(Camera(yaw=yaw, pitch=pitch, rim=1.0), CLOUD, GAMMA, RECT)
        apex = np.array([0.0, 0.0, 0.0])
        normal_tip = np.array([0.0, np.sin(GAMMA), np.cos(GAMMA)])
        for beta in np.radians([-51.84, -25.0, -3.38, 17.74, 41.94]):
            k = diagram3d.wave_vector(beta, GAMMA)
            assert screen_angle(frame, normal_tip, apex, k) == pytest.approx(
                abs(np.degrees(beta)), abs=1e-6)


class TestFitAndZoom:
    def test_the_fit_fills_the_frame(self):
        frame = fit(Camera(), CLOUD, GAMMA, RECT, pad=30)
        xy, _ = project(frame, CLOUD)
        span = xy.max(axis=0) - xy.min(axis=0)
        assert max(span[0] / (900 - 60), span[1] / (640 - 60)) == pytest.approx(1.0)

    def test_zoom_keeps_the_point_under_the_cursor(self):
        cam = Camera()
        frame = fit(cam, CLOUD, GAMMA, RECT)
        at = (612.0, 211.0)
        zoomed = zoom_about(cam, frame, 3.0, at)
        after = fit(zoomed, CLOUD, GAMMA, RECT)
        # The world point that was under `at`, built from the frame's basis.
        sx = (at[0] - frame.centre[0]) / frame.scale - frame.offset[0]
        sy = (frame.centre[1] - at[1]) / frame.scale - frame.offset[1]
        under = sx * frame.right + sy * frame.up
        (before,), _ = project(frame, [under])
        (q,), _ = project(after, [under])
        assert before == pytest.approx(at, abs=1e-9)
        assert q == pytest.approx(at, abs=1e-9)

    def test_zoom_is_clamped_and_the_pan_follows_the_clamp(self):
        cam = Camera(zoom=camera3d.MAX_ZOOM)
        frame = fit(cam, CLOUD, GAMMA, RECT)
        assert zoom_about(cam, frame, 5.0, (10.0, 10.0)) == cam

    def test_the_rim_is_not_offered_in_plane(self):
        frame = fit(Camera(rim=1.0), CLOUD, np.radians(90.0), RECT, rim_offered=False)
        assert not frame.rim_usable
        assert frame.magnification == pytest.approx(1.0)

    def test_recentre_keeps_the_view_and_the_rim(self):
        cam = Camera(yaw=12, pitch=34, roll=56, zoom=3, pan=(7, 8), rim=0.5)
        assert cam.recentred() == Camera(yaw=12, pitch=34, rim=0.5)


def test_reads_no_physics_and_no_toolkit():
    """The camera knows gamma only to place the rim. It imports neither the
    solvers nor Qt, so it can be tested and reasoned about alone."""
    tree = ast.parse(inspect.getsource(camera3d))
    names = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    names |= {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    assert names <= {"__future__", "dataclasses", "numpy", "numpy.typing"}
