r"""The camera for the painted 3D view, without a GUI toolkit.

:mod:`gratinglab.gui.diagram3d` decides what is in the scene; this module
decides only where it lands on the screen; :mod:`gratinglab.gui.qt.scene3d_panel`
paints it. Nothing here reads physics: the one geometric input is the cone's
half-angle, which places the rim the view can zoom to.

The projection is **orthographic**, and every change of framing is a uniform
scale plus a translation (a similarity transform), so every angle in the
picture is the true angle in the plane of the screen at any zoom. That is
the property the matplotlib view protected by refusing to zoom at all; with
its own camera the view can zoom and keep it.

Frames and names
================

The world is the physical :math:`(\hat{d}, \hat{n}, -\hat{g})` frame of
``conventions.md`` §3: x across the grooves, y the grating normal, z along the
cone axis. No permutation is needed, because this camera keeps y up by
choosing its own basis, unlike `Axes3D`, which always draws its z upward.

- ``yaw``, ``pitch``: the orbit, in degrees. The eye sits at
  :math:`(\sin y \cos p,\ \sin p,\ \cos y \cos p)` and looks at the origin.
- ``roll``: the viewer's twist about the line of sight, in degrees. It moves
  nothing physical and is what lets the fan of orders lie level on screen.
- ``zoom``, ``pan``: the user's free transform on top of the fitted frame.
- ``rim``: 0 frames the whole scene, 1 frames the cone rim, and the scale
  between them goes as ``magnification ** rim`` so the slider feels even.

Called by the panel on every paint; the tests drive it directly.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "PITCH_LIMIT",
    "MIN_ZOOM",
    "MAX_ZOOM",
    "RIM_FILL",
    "Camera",
    "Frame",
    "basis",
    "orbit_for",
    "fit",
    "project",
    "zoom_about",
]

#: The orbit stops short of the poles, where "up" stops being defined. 89.5°
#: is close enough that the face-on preset reads as face-on.
PITCH_LIMIT = 89.5
#: The user's zoom, as a multiple of the fitted frame.
MIN_ZOOM, MAX_ZOOM = 0.3, 60.0
#: The rim circle's diameter as a fraction of the frame's short side, at full
#: rim focus. The app uses the same 0.62, so the two products frame alike.
RIM_FILL = 0.62


@dataclass(frozen=True, slots=True)
class Camera:
    """Where the viewer is. Immutable: every gesture returns a new one."""

    yaw: float = 26.6
    pitch: float = 20.4
    roll: float = 0.0
    zoom: float = 1.0
    pan: tuple[float, float] = (0.0, 0.0)
    rim: float = 0.0

    def orbited(self, dyaw: float, dpitch: float) -> "Camera":
        pitch = float(np.clip(self.pitch + dpitch, -PITCH_LIMIT, PITCH_LIMIT))
        return replace(self, yaw=(self.yaw + dyaw) % 360.0, pitch=pitch)

    def rolled(self, droll: float) -> "Camera":
        return replace(self, roll=(self.roll + droll) % 360.0)

    def panned(self, dx: float, dy: float) -> "Camera":
        return replace(self, pan=(self.pan[0] + dx, self.pan[1] + dy))

    def recentred(self) -> "Camera":
        """Clear the free transform. The orbit and the rim slider stay, as in
        the app: Recentre undoes a wander, not a choice of view."""
        return replace(self, roll=0.0, zoom=1.0, pan=(0.0, 0.0))


@dataclass(frozen=True, slots=True)
class Frame:
    """One frame's mapping from world to screen pixels."""

    right: NDArray[np.float64]
    up: NDArray[np.float64]
    fwd: NDArray[np.float64]
    #: Screen point that world (offset) lands on, in pixels.
    centre: tuple[float, float]
    #: World-plane offset applied before scaling: the fit and the rim focus.
    offset: tuple[float, float]
    #: Pixels per world unit.
    scale: float
    #: Pixels per world unit for the fitted frame alone.
    fit_scale: float
    #: Whether zooming to the rim magnifies enough to be worth offering.
    rim_usable: bool

    @property
    def magnification(self) -> float:
        """How much closer than the fitted frame this is."""
        return self.scale / self.fit_scale


def basis(yaw: float, pitch: float, roll: float = 0.0):
    """``(right, up, fwd)`` for an orbit, as unit world vectors."""
    y, p = np.radians(yaw), np.radians(np.clip(pitch, -PITCH_LIMIT, PITCH_LIMIT))
    eye = np.array([np.sin(y) * np.cos(p), np.sin(p), np.cos(y) * np.cos(p)])
    fwd = -eye
    right = np.cross(fwd, [0.0, 1.0, 0.0])
    right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    if roll:
        c, s = np.cos(np.radians(roll)), np.sin(np.radians(roll))
        right, up = right * c + up * s, up * c - right * s
    return right, up, fwd


def orbit_for(direction) -> tuple[float, float]:
    """``(yaw, pitch)`` that looks along a world ``direction``.

    The inverse of :func:`basis`'s ``fwd``, so a preset that names a physical
    direction is the direction the camera ends up looking along."""
    eye = -np.asarray(direction, dtype=float)
    eye /= np.linalg.norm(eye)
    pitch = float(np.degrees(np.arcsin(np.clip(eye[1], -1.0, 1.0))))
    yaw = float(np.degrees(np.arctan2(eye[0], eye[2])))
    return yaw % 360.0, float(np.clip(pitch, -PITCH_LIMIT, PITCH_LIMIT))


def fit(
    camera: Camera,
    points,
    gamma: float,
    rect: tuple[float, float, float, float],
    *,
    pad: float = 30.0,
    rim_offered: bool = True,
) -> Frame:
    """Fit ``points`` into the pixel ``rect`` ``(x0, y0, x1, y1)``.

    Fitted to the scene's own projected extent rather than to a fixed box:
    the bundle is a long diagonal, not a ball, so a box leaves most of the
    frame empty, which is what made the matplotlib view look small. ``pad``
    is reserved in pixels for labels, which do not scale with the scene.
    """
    right, up, fwd = basis(camera.yaw, camera.pitch, camera.roll)
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    sx, sy = pts @ right, pts @ up
    half_x = max(1e-9, (sx.max() - sx.min()) / 2.0)
    half_y = max(1e-9, (sy.max() - sy.min()) / 2.0)
    x0, y0, x1, y1 = rect
    vis_w, vis_h = max(x1 - x0, 80.0), max(y1 - y0, 80.0)
    fit_scale = min((vis_w - 2 * pad) / (2 * half_x), (vis_h - 2 * pad) / (2 * half_y))
    fit_scale = max(fit_scale, 1e-6)
    fit_off = (-(sx.max() + sx.min()) / 2.0, -(sy.max() + sy.min()) / 2.0)

    # The rim: a circle of radius sin(gamma) about (0, 0, cos(gamma)).
    sg = float(np.sin(gamma))
    rim_scale = RIM_FILL * min(vis_w, vis_h) / (2 * sg) if sg > 0 else np.inf
    mag = rim_scale / fit_scale
    usable = bool(rim_offered and np.isfinite(mag) and mag > 1.0)
    t = float(np.clip(camera.rim, 0.0, 1.0)) if usable else 0.0
    ease = 1.0 - (1.0 - t) ** 2  # the centre settles ahead of the zoom
    rc = np.array([0.0, 0.0, float(np.cos(gamma))])
    rim_off = (-float(rc @ right), -float(rc @ up))
    offset = (
        fit_off[0] + (rim_off[0] - fit_off[0]) * ease,
        fit_off[1] + (rim_off[1] - fit_off[1]) * ease,
    )
    scale = fit_scale * (mag if usable else 1.0) ** t * camera.zoom
    centre = ((x0 + x1) / 2.0 + camera.pan[0], (y0 + y1) / 2.0 + camera.pan[1])
    return Frame(right, up, fwd, centre, offset, scale, fit_scale, usable)


def project(frame: Frame, points) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """World points ``(n, 3)`` to screen ``(n, 2)`` pixels and depth ``(n,)``.

    Depth grows away from the viewer, so painting in decreasing depth puts
    near things over far ones."""
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    sx = (pts @ frame.right + frame.offset[0]) * frame.scale + frame.centre[0]
    sy = frame.centre[1] - (pts @ frame.up + frame.offset[1]) * frame.scale
    return np.column_stack([sx, sy]), pts @ frame.fwd


def zoom_about(camera: Camera, frame: Frame, factor: float, at: tuple[float, float]) -> Camera:
    """Zoom by ``factor`` keeping the scene point under ``at`` where it is.

    Clamp first, then pan by the factor actually applied: using the requested
    one makes the scene slide out from under the cursor at the limits."""
    zoom = float(np.clip(camera.zoom * factor, MIN_ZOOM, MAX_ZOOM))
    k = zoom / camera.zoom
    # screen = c0 + pan + v*scale, held fixed for the v under `at`.
    c0 = (frame.centre[0] - camera.pan[0], frame.centre[1] - camera.pan[1])
    pan = tuple(
        a - c - (a - c - p) * k for a, c, p in zip(at, c0, camera.pan)
    )
    return replace(camera, zoom=zoom, pan=pan)
