r"""The conical diffraction geometry in three dimensions, decided without a
toolkit.

The 2D view in :mod:`gratinglab.gui.diagram` is the projection down the groove
axis. It draws every azimuth truly and hides one thing completely: the
out-of-plane component :math:`k_z = k\cos\gamma`, which at
:math:`\gamma = 1.5°` is 99.97% of :math:`|k|`. This module is the view that
shows it.

Why a narrow cone still makes a good picture
============================================

At the reference geometry (315.15 nm, :math:`\alpha = 25°`,
:math:`\gamma = 1.5°`, :math:`\lambda = 3` nm) the propagating orders are:

===  ==========  ============================================
 m   :math:`\beta_m`   :math:`\hat{k}_m`
===  ==========  ============================================
 -1   −51.84°    (−0.02058, +0.01617, **0.99966**)
  0   −25.00°    (−0.01106, +0.02372, **0.99966**)
 +1    −3.38°    (−0.00154, +0.02613, **0.99966**)
 +2   +17.74°    (+0.00798, +0.02493, **0.99966**)
 +3   +41.94°    (+0.01749, +0.01947, **0.99966**)
===  ==========  ============================================

**Polar spread: exactly 0.0°. Azimuth spread: 93.8°.** The cone is a needle
*and* the fan around it is nearly a right angle. Both at once, and that
combination is precisely what a dispersion-plane drawing cannot convey.

So the default oblique view is not a degenerate picture. The grooves run along
:math:`\hat{z}` and the rays run at 1.5° to them: a beam *skimming along the
groove direction* at grazing incidence, which is exactly what distinguishes an
off-plane mount from a classical one. Looking down the cone axis instead, and
zooming uniformly to the rim, makes the 93.8° fan fully readable -- at true
scale, since uniform zoom is a similarity transform and preserves every angle.

**Nothing here exaggerates γ.** Legibility comes from viewpoint and
magnification only.

Two scales in one box
=====================

The grating patch is in nm (~630 nm across); the rays are unit vectors. There
is no physical ratio between them, so :attr:`Scene.scale_nm` records the one
chosen and a caption states it. The groove *cross-section* is exact and every
ray *direction* is exact; only the patch's size relative to the rays is a
drawing choice -- the same disclaimer the 2D view already carries for ray
length, narrowed to the one thing still arbitrary.

Why a separate module
=====================

:class:`~gratinglab.gui.diagram.Diagram` is 2D by construction: ``limits`` is
typed as two ranges and the widget unpacks it as such, and
:attr:`Arrow.azimuth` would silently keep returning a projection angle if a
``z`` were defaulted onto it. Adding a third dimension there would make the
widget branch on shape, which is what the pure/impure split exists to prevent.
This module reuses :func:`~gratinglab.gui.diagram.direction`,
:func:`~gratinglab.gui.diagram.order_marks`,
:func:`~gratinglab.gui.diagram.strike_point` and
:func:`~gratinglab.gui.diagram.facet_normal` rather than re-deriving anything.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

import numpy as np
from numpy.typing import NDArray

from ..geometry import blaze_direction
from .diagram import (
    TAG_COLORS,
    OrderMark,
    direction,
    facet_normal,
    label_orders,
    order_marks,
    strike_point,
    x_nm,
)
from .provenance import Line

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..illumination import Illumination
    from ..problem import Problem

__all__ = [
    "G_HAT",
    "D_HAT",
    "N_HAT",
    "CONE_AXIS",
    "TAG_COLORS",
    "PRESET_VIEWS",
    "RIM_PRESET",
    "rim_box",
    "Ray3D",
    "Surface3D",
    "Curve3D",
    "Point3D",
    "Scene",
    "wave_vector",
    "incident_vector",
    "view_direction",
    "build_scene",
    "captions3d",
]

#: `conventions.md` §3, forced by :math:`\hat{d} \times \hat{g} = \hat{n}`.
D_HAT = np.array([1.0, 0.0, 0.0])
N_HAT = np.array([0.0, 1.0, 0.0])
G_HAT = np.array([0.0, 0.0, -1.0])

#: The cone opens along :math:`-\hat{g}`, **not** along :math:`\hat{g}`. Every
#: wave vector carries :math:`k_z = +k\cos\gamma`, so at :math:`\gamma = 1.5°`
#: every ray lies 1.5° from :math:`+\hat{z}` and 178.5° from :math:`\hat{g}`.
#: Opening it the other way points every ray 180° out. See
#: ``docs/conventions.md`` §3 and ``docs/findings.md``.
CONE_AXIS = -G_HAT

#: Grating patch width, in ray-length units. The rays are unit vectors and the
#: patch is in nm; this is the ratio between those two scales, and it is a
#: drawing choice that a caption declares.
#:
#: Small on purpose. At gamma = 1.5 deg the rays lie within 1.5 deg of the
#: groove axis, so a patch as long as the rays swallows them; keeping it well
#: under a ray length lets the bundle emerge and be seen for what it is.
PATCH_SPAN: float = 0.5

#: Half-extent of the scene's bounding cube.
BOX_HALF: float = 1.15

#: How far out the rim-focus box reaches, in units of the rim radius
#: :math:`\sin\gamma`. 1.6 puts the rim circle across ~62% of the frame.
RIM_MARGIN: float = 1.6

#: Below this, focusing on the rim is not worth offering -- and at
#: :math:`\gamma = 90°` it would be a zoom *out*.
MIN_RIM_MAGNIFICATION: float = 1.5

#: The direction the camera looks **along**, in the physical
#: :math:`(\hat{d}, \hat{n}, -\hat{g})` frame of ``conventions.md`` §3.
#:
#: Physical, not matplotlib ``(elev, azim)``: the widget permutes the scene to
#: put :math:`\hat{n}` upward for display, so an ``(elev, azim)`` pair here
#: would no longer mean what its name says. The widget converts.
PRESET_VIEWS: dict[str, NDArray[np.float64]] = {
    # Off-axis: the grazing bundle skimming along the grooves.
    "Oblique": np.array([-0.42, -0.35, -0.84]),
    # From +ẑ looking back, so every order projects onto a circle at its true
    # azimuth. Paired with rim focus, this is where the fan is readable.
    "Down the cone axis": G_HAT.copy(),
    # The (n̂, ĝ) plane, where γ itself is the visible angle.
    "Along d̂": -D_HAT.copy(),
    # The grating face-on, grooves running across.
    "Face n̂": -N_HAT.copy(),
}
PRESET_VIEWS = {
    name: v / np.linalg.norm(v) for name, v in PRESET_VIEWS.items()
}

#: The preset whose whole point is the rim, so choosing it switches focus too.
RIM_PRESET = "Down the cone axis"

#: Mesh resolutions. Both surfaces are *ruled* -- the grating patch is the
#: cross-section extruded along ẑ, and every generator of a cone is straight --
#: so the second dimension is 2 and that is exact, not a sampling choice.
_PATCH_NX = {"fine": 241, "coarse": 61}
_CONE_NPHI = {"fine": 289, "coarse": 73}


@dataclass(frozen=True, slots=True)
class Ray3D:
    """A wave vector, drawn from ``origin`` along a **unit** ``direction``.

    Unit, always. The 2D view draws every ray at one arbitrary length because
    only azimuths carry meaning there; here the transverse extent *is* the
    physics -- it is :math:`\\sin\\gamma`, and at 1.5° that is 0.026. A
    renormalised ray would erase exactly the quantity this view exists to
    show.
    """

    origin: tuple[float, float, float]
    direction: tuple[float, float, float]
    length: float
    tag: str
    label: str | None = None
    order: int | None = None
    dashed: bool = False

    @property
    def head(self) -> tuple[float, float, float]:
        o, d, s = np.asarray(self.origin), np.asarray(self.direction), self.length
        return tuple(float(v) for v in o + s * d)

    @property
    def azimuth(self) -> float:
        """Azimuth from :math:`\\hat{n}` toward :math:`\\hat{d}`, radians --
        the convention every angle in §3 uses. Read back off the drawn vector,
        which is what lets a test recover the grating equation from the
        picture rather than recompute it the same way."""
        return float(np.arctan2(self.direction[0], self.direction[1]))

    @property
    def polar_from_cone_axis(self) -> float:
        """Angle to :math:`-\\hat{g}`. Equals :math:`\\gamma` for every wave
        vector, which is the cone property stated as one number."""
        return float(np.arccos(np.dot(self.direction, CONE_AXIS)))


@dataclass(frozen=True, slots=True)
class Surface3D:
    """A ruled surface, as ``(n, 2)`` meshgrid arrays for ``plot_surface``."""

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    z: NDArray[np.float64]
    tag: str
    lod: str = "fine"
    alpha: float = 1.0
    label: str | None = None


@dataclass(frozen=True, slots=True)
class Curve3D:
    """An open or closed 3D path, already sampled."""

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    z: NDArray[np.float64]
    tag: str
    lod: str = "fine"
    dashed: bool = False
    label: str | None = None


@dataclass(frozen=True, slots=True)
class Point3D:
    """A point with a *named* glyph; the widget owns the name-to-marker map."""

    x: float
    y: float
    z: float
    glyph: str
    tag: str
    label: str | None = None


@dataclass(frozen=True, slots=True)
class Scene:
    """One geometry, at one wavelength, as things to draw in three dimensions."""

    rays: tuple[Ray3D, ...]
    surfaces: tuple[Surface3D, ...]
    curves: tuple[Curve3D, ...]
    points: tuple[Point3D, ...]
    orders: tuple[OrderMark, ...]
    captions: tuple[Line, ...]
    title: str
    limits: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    #: The same scene, framed on the cone rim. **Also a cube**, so moving
    #: between the two is a uniform scale plus a translation -- a similarity
    #: transform, which preserves every angle exactly. Nothing is distorted;
    #: only which part of the geometry fills the frame changes.
    rim_limits: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    #: How much closer that is, for the label on the canvas.
    rim_magnification: float
    #: Non-empty **exactly** when the rim focus is not worth offering, so a
    #: caller puts *this* in a disabled control's tooltip rather than
    #: inventing one. Same contract as ``blaze_jump``.
    rim_reason: str
    scale_nm: float
    gamma: float
    wavelength: float
    blaze_order: int | None

    @property
    def rim_available(self) -> bool:
        return not self.rim_reason

    def ray(self, order: int) -> Ray3D:
        """The ray for one order. Raises if it has none -- an evanescent order
        genuinely has no direction, and a default would invent one."""
        for entry in self.rays:
            if entry.order == order:
                return entry
        raise KeyError(f"order {order} has no ray in this scene")

    def at(self, lod: str) -> tuple[object, ...]:
        """Every level-of-detail-bearing primitive at one level."""
        return tuple(
            item
            for group in (self.surfaces, self.curves)
            for item in group
            if item.lod == lod
        )


# -- wave vectors ---------------------------------------------------------


def wave_vector(azimuth: float, gamma: float) -> NDArray[np.float64]:
    r"""A unit wave vector at the given azimuth, on the cone of half-angle
    :math:`\gamma`.

    ``conventions.md`` §3's :math:`\hat{k}_m`. The :math:`(x, y)` part is
    literally :func:`~gratinglab.gui.diagram.direction` scaled by
    :math:`\sin\gamma`, so the 2D and 3D views cannot disagree about a
    direction; only :math:`+\cos\gamma\,\hat{z}` is new here, and that is
    exactly the component the projection drops.
    """
    return np.array([*(direction(azimuth) * np.sin(gamma)), np.cos(gamma)])


def incident_vector(alpha: float, gamma: float) -> NDArray[np.float64]:
    r""":math:`\hat{k}_i`, which travels *toward* the grating (negative
    :math:`\hat{y}`). Equals ``Illumination.direction_cosines``."""
    return np.array([*(-direction(alpha) * np.sin(gamma)), np.cos(gamma)])


def rim_box(
    gamma: float,
) -> tuple[
    tuple[tuple[float, float], tuple[float, float], tuple[float, float]], float, str
]:
    r"""``(limits, magnification, reason)`` for a view framed on the cone rim.

    The rim is a circle of radius :math:`\sin\gamma` at height
    :math:`\cos\gamma`, so the box is a cube centred there. At
    :math:`\gamma = 1.5°` that is a 27x uniform magnification -- which is what
    turns the 93.8° azimuth fan from a point cluster into something readable,
    without touching a single angle.

    ``reason`` is non-empty **exactly** when the magnification is not worth
    offering. At :math:`\gamma = 90°` the rim already *is* the whole scene and
    focusing on it would zoom *out*; a control that did nothing would be the
    mistake this project has named before.
    """
    half = RIM_MARGIN * float(np.sin(gamma))
    magnification = BOX_HALF / half if half > 0.0 else float("inf")
    centre = float(np.cos(gamma))

    limits = (
        (-half, half),
        (-half, half),
        (centre - half, centre + half),
    )

    if magnification < MIN_RIM_MAGNIFICATION:
        degrees = np.degrees(gamma)
        if np.isclose(degrees, 90.0):
            reason = (
                "γ = 90°: the cone rim is the whole scene, so there is nothing "
                "to zoom into."
            )
        else:
            reason = (
                f"γ = {degrees:g}° already separates the orders; focusing on "
                "the rim would change little."
            )
        return limits, magnification, reason

    return limits, magnification, ""


def view_direction(elev_deg: float, azim_deg: float) -> NDArray[np.float64]:
    """The unit vector a matplotlib 3D view looks *along*, from the camera
    toward the origin's far side. Pure, so presets are testable without a
    figure."""
    elev, azim = np.radians(elev_deg), np.radians(azim_deg)
    return np.array(
        [
            -np.cos(elev) * np.cos(azim),
            -np.cos(elev) * np.sin(azim),
            -np.sin(elev),
        ]
    )


# -- captions -------------------------------------------------------------


def scene_title(gamma: float) -> str:
    """Always states γ, and always states that it is to scale."""
    return f"Conical diffraction — γ = {np.degrees(gamma):g}°, to scale"


def captions3d(
    problem: "Problem",
    illumination: "Illumination",
    wavelength: float,
    marks: Sequence[OrderMark],
    scale_nm: float,
) -> tuple[Line, ...]:
    """What the drawing is, and what in it is a drawing choice."""
    live = sum(1 for m in marks if m.propagating)
    dark = len(marks) - live
    gamma = illumination.gamma

    lines: list[Line] = [
        Line(
            f"λ = {wavelength:.3g} nm. {live} of {len(marks)} orders propagate "
            f"and carry a ray; {dark} are evanescent, have no direction, and "
            "appear on the sin β axis beside this view.\n"
        ),
        Line(
            "Every ray is a unit wave vector. Their directions, and the angles "
            "between them, are exact.\n",
            "dim",
        ),
        Line(
            f"The grating patch is in nm and the rays are unit vectors, so "
            f"their relative size is a drawing choice ({scale_nm:.1f} nm per "
            "scene unit). The groove cross-section itself is exact.\n",
            "dim",
        ),
    ]

    if illumination.is_in_plane:
        lines.append(
            Line(
                "γ = 90°: k_z = 0, the cone degenerates to the dispersion "
                "plane, and every ray lies in it. This view and the "
                "cross-section beside it show the same geometry.\n",
                "dim",
            )
        )
    else:
        lines.append(
            Line(
                f"The cone opens along −ĝ, not ĝ. Every wave vector carries "
                f"k_z = +k cos γ = {np.cos(gamma):.5f} k, so each ray lies "
                f"{np.degrees(gamma):g}° from +ẑ and "
                f"{180 - np.degrees(gamma):g}° from ĝ.\n",
                "dim",
            )
        )
        lines.append(
            Line(
                "The orders differ in azimuth, never in polar angle — they all "
                "sit on this one cone. That is why a cone this narrow can still "
                "carry a wide fan.\n",
                "dim",
            )
        )

    lines.append(
        Line(
            "The groove profile has no z dependence, so the surface drawn is "
            "the cross-section extruded exactly along ẑ. That is not an "
            "approximation.\n",
            "dim",
        )
    )
    return tuple(lines)


# -- the scene ------------------------------------------------------------


def _grating_patch(problem: "Problem", scale_nm: float, periods: int, lod: str):
    """The cross-section extruded along ẑ. Ruled, so two rows is exact."""
    t = np.linspace(-periods, 0.0, _PATCH_NX[lod])
    strike = np.asarray(strike_point(problem))
    px = (x_nm(t, problem.period) - strike[0]) / scale_nm
    py = (np.asarray(problem.height_nm(t)) - strike[1]) / scale_nm
    half = 0.5 * PATCH_SPAN

    return Surface3D(
        x=np.column_stack([px, px]),
        y=np.column_stack([py, py]),
        z=np.column_stack([np.full_like(px, -half), np.full_like(px, half)]),
        tag="surface",
        lod=lod,
        alpha=0.55,
        label="grating",
    )


def _cone_mantle(gamma: float, lod: str) -> Surface3D:
    """The cone of half-angle γ about −ĝ. Ruled along its generators."""
    phi = np.linspace(0.0, 2.0 * np.pi, _CONE_NPHI[lod])
    rim = np.array([np.sin(gamma) * np.cos(phi), np.sin(gamma) * np.sin(phi)])
    s = np.array([0.0, 1.0])

    return Surface3D(
        x=np.outer(rim[0], s),
        y=np.outer(rim[1], s),
        z=np.outer(np.full_like(phi, np.cos(gamma)), s),
        tag="cone",
        lod=lod,
        alpha=0.10,
        label="cone of diffraction",
    )


def _cone_rim(gamma: float, lod: str) -> Curve3D:
    phi = np.linspace(0.0, 2.0 * np.pi, _CONE_NPHI[lod])
    return Curve3D(
        x=np.sin(gamma) * np.cos(phi),
        y=np.sin(gamma) * np.sin(phi),
        z=np.full_like(phi, np.cos(gamma)),
        tag="cone",
        lod=lod,
        dashed=True,
    )


def build_scene(
    problem: "Problem",
    illumination: "Illumination",
    wavelength: float,
    *,
    periods: int = 2,
    pad_orders: int = 2,
) -> Scene:
    """Everything one 3D drawing consists of."""
    gamma = illumination.gamma
    marks = order_marks(problem, illumination, wavelength, pad=pad_orders)
    labelled = label_orders(marks)
    scale_nm = periods * problem.period / PATCH_SPAN
    origin = (0.0, 0.0, 0.0)
    rim_limits, rim_magnification, rim_reason = rim_box(gamma)

    rays: list[Ray3D] = []
    surfaces: list[Surface3D] = []
    curves: list[Curve3D] = []
    points: list[Point3D] = []

    for lod in ("fine", "coarse"):
        surfaces.append(_grating_patch(problem, scale_nm, periods, lod))
        surfaces.append(_cone_mantle(gamma, lod))
        curves.append(_cone_rim(gamma, lod))

    # -- axes ------------------------------------------------------------
    rays.append(Ray3D(origin, tuple(G_HAT), 0.85, "axis", label="ĝ (grooves)"))
    rays.append(Ray3D(origin, tuple(CONE_AXIS), 1.0, "cone",
                      label="cone axis −ĝ", dashed=True))
    rays.append(Ray3D(origin, tuple(D_HAT), 0.45, "axis", label="d̂"))
    rays.append(Ray3D(origin, tuple(N_HAT), 0.45, "normal", label="n̂"))

    nx, ny = facet_normal(problem)
    rays.append(
        Ray3D(origin, (float(nx), float(ny), 0.0), 0.5, "normal",
              label="facet normal", dashed=True)
    )

    # -- the waves -------------------------------------------------------
    incident = incident_vector(illumination.alpha, gamma)
    rays.append(
        Ray3D(
            tuple(-incident), tuple(incident), 1.0, "incident", label="k̂ᵢ"
        )
    )

    for mark in marks:
        if mark.beta is None:
            continue  # no direction; it stays on the ladder beside this view
        k = wave_vector(mark.beta, gamma)
        rays.append(
            Ray3D(
                origin, tuple(k), 1.0,
                "zero" if mark.order == 0 else "order",
                # 3D text cannot dodge its neighbours, and every order tip is
                # within 1.5 deg of the same point. Only m=0 is labelled here;
                # the ladder beside this view names them all.
                label=f"m={mark.order:+d}" if mark.order == 0 else None,
                order=mark.order,
            )
        )
        points.append(Point3D(*k, glyph="out-of-page", tag="cone"))
        # Its shadow on the dispersion plane -- the 2D view, literally.
        curves.append(
            Curve3D(
                np.array([0.0, k[0]]), np.array([0.0, k[1]]), np.zeros(2),
                tag="axis", dashed=True,
            )
        )

    blaze_angle = getattr(problem.profile, "blaze_angle", None)
    if blaze_angle is not None:
        beta_b = blaze_direction(np.radians(blaze_angle), illumination.alpha)
        rays.append(
            Ray3D(origin, tuple(wave_vector(beta_b, gamma)), 1.0, "blaze",
                  label=f"β_b = {np.degrees(beta_b):.1f}°")
        )

    points.append(Point3D(0.0, 0.0, 0.0, glyph="strike", tag="incident"))

    return Scene(
        rays=tuple(rays),
        surfaces=tuple(surfaces),
        curves=tuple(curves),
        points=tuple(points),
        orders=marks,
        captions=captions3d(problem, illumination, wavelength, marks, scale_nm),
        title=scene_title(gamma),
        limits=(
            (-BOX_HALF, BOX_HALF),
            (-BOX_HALF, BOX_HALF),
            (-BOX_HALF, BOX_HALF),
        ),
        rim_limits=rim_limits,
        rim_magnification=rim_magnification,
        rim_reason=rim_reason,
        scale_nm=scale_nm,
        gamma=float(gamma),
        wavelength=float(wavelength),
        blaze_order=None,
    )
