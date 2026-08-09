r"""What a grating geometry looks like, decided without a toolkit.

Returns drawable primitives -- polylines, arrows, markers, tagged captions --
in physical coordinates (nm, in the :math:`(\hat{d}, \hat{n})` frame of
``docs/conventions.md`` §3). The widget positions them and colours them by tag;
it computes no angle, no direction and no cutoff, so a bug there misplaces a
line rather than drawing wrong physics. Same split, and the same reason, as
:mod:`gratinglab.gui.provenance`: which line is red is a correctness question.

This module **arranges**; it never re-derives. Every angle comes from
:mod:`gratinglab.geometry`, every height from :meth:`Problem.height_nm`.

Three panels, because the physics needs three
=============================================

**main** -- the view *down the groove axis* :math:`\hat{g}`, i.e. the
projection onto the :math:`\hat{n}`-:math:`\hat{d}` plane. In that projection
:math:`\mathbf{k}_i` and every :math:`\mathbf{k}_m` have the same projected
length :math:`k\sin\gamma`, and their azimuths are exactly :math:`\alpha` and
:math:`\beta_m`. Two things follow, and they are why this is a good view rather
than a compromise:

1. Every propagating order lands on one circle centred at the strike point.
   That circle *is* the diffraction cone, seen end-on -- so "the mount is
   conical" becomes a picture rather than a caption.
2. :math:`\alpha`, :math:`\beta_m`, :math:`\delta` and :math:`\beta_b` are all
   azimuths in this same plane, drawn true at any :math:`\gamma`.

**ladder** -- every order on a :math:`\sin\beta` axis. The main panel can only
draw an order that has a direction; an evanescent one has none, and
``conventions.md`` §4 forbids dropping it ("retained in the record, never
silently dropped, never NaN"). So an evanescent order sits here, outside the
:math:`[-1, 1]` propagating window, with no ray attached. Omitting it from the
picture would be exactly the silent drop that section exists to prevent.

There used to be a third panel here, drawing :math:`\gamma` as a true-to-scale
sliver. It existed because the main panel is a projection: :math:`k_z =
k\cos\gamma` is perpendicular to it and invisible, and at :math:`\gamma = 1.5°`
that is 99.97% of :math:`|k|`, so a reader who sees rays at 25° from the normal
could conclude the mount is not grazing. :mod:`gratinglab.gui.diagram3d` now
draws the actual cone, with :math:`\gamma` as a real angle in a real scene
rather than a stand-in beside one -- so the sliver was retired (M13-I). Two
drawings of one angle are two answers to one question. The captions still say
what this projection cannot show, and now point at the view that can.

What is deliberately *not* drawn
================================

:math:`\zeta`, the facet graze angle, is **text only**.
:func:`~gratinglab.geometry.facet_graze` is the three-dimensional angle between
:math:`\hat{k}_i` and the facet *plane*, :math:`\sin\zeta = \sin\gamma
\cos(\delta-\alpha)`. It is **not** the angle between the projected ray and the
facet line in this drawing -- 1.495° against 85.5° at the reference geometry.
An arc labelled :math:`\zeta` would be the worst lie available here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..geometry import (
    blaze_direction,
    blaze_wavelength,
    facet_graze,
    is_propagating,
    order_range,
    sin_beta,
)
from .provenance import Line

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..illumination import Illumination
    from ..problem import Problem

__all__ = [
    "T_AXIS_SIGN",
    "TAG_COLORS",
    "MAX_ORDER_LABELS",
    "Polyline",
    "Arrow",
    "Marker",
    "OrderMark",
    "Diagram",
    "build",
    "x_nm",
    "direction",
    "order_span",
    "order_marks",
    "strike_point",
    "facet_normal",
    "label_orders",
    "captions",
    "nearest_index",
    "blaze_targets",
    "blaze_jump",
]

#: Which way the profile parameter runs along the dispersion axis.
#:
#: ``Profile.height(t)`` is a shape in ``t``; :mod:`gratinglab.profiles` never
#: says which physical direction ``+t`` points. The wave conventions do fix it:
#: ``blaze_direction`` gives :math:`\beta_b = 2\delta - \alpha`, which requires
#: the active facet's outward normal at azimuth :math:`+\delta`, while
#: ``Blazed.height`` rises with ``t``, putting it at :math:`-\delta` in ``t``.
#:
#: Nothing before this module drew a profile and a ray in the same frame, so
#: nothing before it had to know. See ``docs/conventions.md`` §3 and
#: ``tests/test_geometry.py::TestFacetHandedness``.
T_AXIS_SIGN: float = -1.0

#: Tag -> colour, exactly as :data:`gratinglab.gui.provenance.TAG_COLORS`
#: works: this module decides what a thing *is*, the widget looks up what that
#: colour is.
TAG_COLORS: dict[str, str] = {
    "surface": "#1f3b73",  # the blue ProfilePlotPanel already uses
    "facet": "#2d5aa8",  # the active facet, picked out of the surface
    "incident": "#8fa3c8",  # icon.py's incident-ray grey-blue
    "zero": "#666666",
    "order": "#4dc4ff",
    "blaze": "#ffc233",
    "normal": "#0a7d28",
    "cone": "#999999",
    "evanescent": "#b00020",
    "axis": "#888888",
}

#: At most this many order arrows carry a label before thinning starts. An
#: X-ray geometry at the short end of a scan propagates sixteen orders; sixteen
#: labels is a smear. The **label** is what gets dropped -- never the arrow,
#: which would be a silent omission -- and a caption reports the true total.
MAX_ORDER_LABELS = 7

#: Ray length, in periods. Arbitrary: only directions carry meaning, and a
#: caption says so.
_RAY_PERIODS = 1.15

#: Points used to sample the surface and the cone arc.
_SURFACE_POINTS = 400
_CONE_POINTS = 181


@dataclass(frozen=True, slots=True)
class Polyline:
    """An open path, already sampled.

    Arcs and angle wedges are polylines too: a ``matplotlib.patches.Arc`` takes
    degrees counter-clockwise from ``+x``, which is a *different* angle
    convention from this project's azimuth-from-:math:`\\hat{n}`, and that
    conversion does not belong in a widget.
    """

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    tag: str
    panel: str = "main"
    fill_to: float | None = None
    dashed: bool = False
    label: str | None = None


@dataclass(frozen=True, slots=True)
class Arrow:
    """A directed segment, tail to head, in absolute coordinates.

    Endpoints rather than angle-and-length on purpose: the widget then converts
    nothing, and a sign error has nowhere to hide in the untested layer.
    """

    x0: float
    y0: float
    x1: float
    y1: float
    tag: str
    panel: str = "main"
    label: str | None = None
    order: int | None = None
    dashed: bool = False

    @property
    def azimuth(self) -> float:
        """Direction as an azimuth from :math:`\\hat{n}` toward :math:`\\hat{d}`,
        radians. Derived for tests; the widget never needs it."""
        return float(np.arctan2(self.x1 - self.x0, self.y1 - self.y0))


@dataclass(frozen=True, slots=True)
class Marker:
    """A point with a *named* glyph.

    ``glyph`` is one of ``"out-of-page"``, ``"strike"``. Named, not styled: the
    widget owns the mapping from name to marker.
    """

    x: float
    y: float
    glyph: str
    tag: str
    panel: str = "main"
    label: str | None = None


@dataclass(frozen=True, slots=True)
class OrderMark:
    """One diffraction order, on the :math:`\\sin\\beta` ladder.

    Present whether or not it propagates. ``sin_beta`` is unclipped, straight
    from :func:`~gratinglab.geometry.sin_beta`, so an evanescent order reads
    :math:`|\\sin\\beta| > 1` and lands visibly outside the window rather than
    vanishing.

    ``beta`` is ``None`` **exactly** when the order is evanescent -- not
    ``nan``, not ``0.0``. An evanescent order has no direction, and a number
    here would invite one.
    """

    order: int
    sin_beta: float
    beta: float | None
    propagating: bool
    is_blaze_order: bool


@dataclass(frozen=True, slots=True)
class Diagram:
    """One geometry, at one wavelength, as things to draw."""

    paths: tuple[Polyline, ...]
    arrows: tuple[Arrow, ...]
    markers: tuple[Marker, ...]
    orders: tuple[OrderMark, ...]
    captions: tuple[Line, ...]
    limits: dict[str, tuple[tuple[float, float], tuple[float, float]]]
    wavelength: float
    blaze_order: int | None

    def mark(self, order: int) -> OrderMark:
        """The :class:`OrderMark` for one order. Raises if it is not present --
        an order that should be in the record and isn't is a bug, not a
        default."""
        for entry in self.orders:
            if entry.order == order:
                return entry
        raise KeyError(f"order {order} is not in this diagram")

    def on(self, panel: str) -> tuple[object, ...]:
        """Every primitive belonging to one panel, for the widget's dispatch."""
        return tuple(
            item
            for group in (self.paths, self.arrows, self.markers)
            for item in group
            if item.panel == panel
        )


# -- coordinates ---------------------------------------------------------


def x_nm(t: ArrayLike, period: float) -> NDArray[np.float64]:
    """Normalised profile position to physical position along :math:`\\hat{d}`,
    in nm. See :data:`T_AXIS_SIGN`."""
    return T_AXIS_SIGN * np.asarray(t, dtype=np.float64) * period


def direction(azimuth: float) -> NDArray[np.float64]:
    """Unit vector at an azimuth measured from :math:`\\hat{n}` toward
    :math:`\\hat{d}` -- the convention every angle in §3 uses."""
    return np.array([np.sin(azimuth), np.cos(azimuth)])


# -- orders --------------------------------------------------------------


def order_span(
    wavelength: float,
    period: float,
    sin_alpha: float,
    sin_gamma: float = 1.0,
    *,
    pad: int = 2,
) -> NDArray[np.int64]:
    """Every propagating order, plus ``pad`` evanescent ones on each side.

    :func:`~gratinglab.geometry.order_range` returns exactly the propagating
    set, which is the right answer for a solver and the wrong one for a
    picture: an order that has just passed off leaves no trace, and the passing
    off is the thing worth watching.
    """
    live = order_range(wavelength, period, sin_alpha, sin_gamma)
    return np.arange(int(live.min()) - pad, int(live.max()) + pad + 1, dtype=np.int64)


def order_marks(
    problem: "Problem",
    illumination: "Illumination",
    wavelength: float,
    *,
    pad: int = 2,
) -> tuple[OrderMark, ...]:
    """Describe every order in the span, propagating or not."""
    orders = order_span(
        wavelength,
        problem.period,
        illumination.sin_alpha,
        illumination.sin_gamma,
        pad=pad,
    )
    values = sin_beta(
        orders,
        wavelength,
        problem.period,
        illumination.sin_alpha,
        illumination.sin_gamma,
    )
    live = is_propagating(values)
    blaze = _blaze_order(problem, illumination, wavelength)

    return tuple(
        OrderMark(
            order=int(m),
            sin_beta=float(s),
            beta=float(np.arcsin(s)) if bool(p) else None,
            propagating=bool(p),
            is_blaze_order=(blaze is not None and int(m) == blaze),
        )
        for m, s, p in zip(orders, values, live)
    )


def label_orders(
    marks: Sequence[OrderMark], *, limit: int = MAX_ORDER_LABELS
) -> frozenset[int]:
    """Which orders get a text label.

    Order 0 and the blaze order always, then the rest by descending
    :math:`|\\sin\\beta|` proximity to the centre, so labels thin from the
    crowded edges inward. Only labels are dropped; every propagating order
    keeps its arrow.
    """
    live = [m for m in marks if m.propagating]
    always = {m.order for m in live if m.order == 0 or m.is_blaze_order}
    rest = sorted(
        (m for m in live if m.order not in always), key=lambda m: abs(m.sin_beta)
    )
    room = max(0, limit - len(always))
    return frozenset(always | {m.order for m in rest[:room]})


# -- the grating itself ---------------------------------------------------


def _active_facet_t(problem: "Problem") -> tuple[float, float]:
    """The ``t`` span of the facet the beam is drawn striking, within the
    period that occupies the middle of the window.

    For a sawtooth that is the active (blaze) facet, ``[0, apex]``. For every
    other profile there is no distinguished facet, so the anchor is the crest
    -- the ray *directions* are the physics, and the origin is a visual anchor
    that a caption already declares arbitrary.
    """
    profile = problem.profile
    apex = getattr(profile, "apex", None)
    if apex is None:
        t = np.linspace(0.0, 1.0, 512, endpoint=False)
        apex = float(t[int(np.argmax(profile.height(t)))])
    return -1.0 + 0.0, -1.0 + float(apex)


def strike_point(problem: "Problem") -> tuple[float, float]:
    """Where the incident ray is drawn meeting the grating, in nm."""
    start, end = _active_facet_t(problem)
    t = 0.5 * (start + end)
    return float(x_nm(t, problem.period)), float(problem.height_nm(t))


def facet_normal(problem: "Problem") -> tuple[float, float]:
    """Outward unit normal of the struck facet, in the physical frame.

    Taken from :meth:`Profile.boundary`, which already returns outward normals
    (``ny > 0``), with the ``x`` component flipped by :data:`T_AXIS_SIGN` --
    the normal is a direction in ``t``-space and this frame runs the other way.
    """
    start, end = _active_facet_t(problem)
    t_mid = np.mod(0.5 * (start + end), 1.0)

    curve = problem.profile.boundary(512)
    index = int(np.argmin(np.abs(curve.t - t_mid)))
    nx = T_AXIS_SIGN * float(curve.nx[index])
    ny = float(curve.ny[index])
    norm = float(np.hypot(nx, ny))
    return nx / norm, ny / norm


def _blaze_order(
    problem: "Problem", illumination: "Illumination", wavelength: float
) -> int | None:
    """The order whose blaze wavelength is nearest the one being drawn, if the
    profile has a blaze angle at all."""
    targets = blaze_targets(problem, illumination)
    if not targets:
        return None
    order, _ = min(targets, key=lambda pair: abs(pair[1] - wavelength))
    return order


# -- wavelength selection -------------------------------------------------


def nearest_index(wavelengths: NDArray[np.float64], target: float) -> int:
    """Index of the scan point closest to ``target``."""
    return int(np.argmin(np.abs(np.asarray(wavelengths, dtype=np.float64) - target)))


def blaze_targets(
    problem: "Problem",
    illumination: "Illumination",
    orders: Iterable[int] = range(1, 9),
) -> tuple[tuple[int, float], ...]:
    """``(order, blaze wavelength)`` pairs, or empty for a profile with no
    blaze angle."""
    blaze_angle = getattr(problem.profile, "blaze_angle", None)
    if blaze_angle is None:
        return ()
    delta = np.radians(blaze_angle)
    return tuple(
        (
            int(m),
            float(
                blaze_wavelength(
                    m, problem.period, delta, illumination.alpha, illumination.gamma
                )
            ),
        )
        for m in orders
    )


def blaze_jump(
    wavelengths: NDArray[np.float64],
    problem: "Problem",
    illumination: "Illumination",
) -> tuple[int | None, str]:
    """``(index to jump to, reason)`` for a "go to the blaze wavelength"
    control.

    The reason is non-empty **exactly** when the index is ``None``, so a caller
    puts *this* text in a disabled button's tooltip instead of inventing one --
    or leaving a dead control that says nothing. A disabled control that
    explains itself is honest; an enabled one that silently does nothing is
    not.
    """
    targets = blaze_targets(problem, illumination)
    if not targets:
        kind = type(problem.profile).__name__
        return None, f"{kind} has no blaze angle, so no blaze wavelength."

    lo = float(np.min(wavelengths))
    hi = float(np.max(wavelengths))
    inside = [(m, lam) for m, lam in targets if lo <= lam <= hi]
    if inside:
        order, lam = min(inside, key=lambda pair: pair[0])
        return nearest_index(wavelengths, lam), ""

    nearby = ", ".join(f"{lam:.2f} nm (m={m})" for m, lam in targets[:3])
    return None, (
        f"No blaze wavelength falls inside the {lo:.2f}-{hi:.2f} nm scan: {nearby}."
    )


# -- captions -------------------------------------------------------------


def captions(
    problem: "Problem",
    illumination: "Illumination",
    wavelength: float,
    marks: Sequence[OrderMark],
) -> tuple[Line, ...]:
    """What the drawing cannot show, said in words.

    Every one of these exists because the picture is true but partial, and a
    partial truth left unlabelled is how a reader talks themselves into the
    wrong geometry.
    """
    live = sum(1 for m in marks if m.propagating)
    dark = len(marks) - live
    lines: list[Line] = [
        Line(
            f"λ = {wavelength:.3g} nm. {live} of {len(marks)} orders shown "
            f"propagate; {dark} are evanescent and appear on the sin β axis "
            "with no ray, having no direction.\n"
        ),
        Line(
            "View along the groove axis ĝ. Angles in this plane — α, β_m, δ, "
            "β_b — are true, and the groove cross-section is to scale.\n",
            "dim",
        ),
    ]

    if illumination.is_in_plane:
        lines.append(
            Line(
                "γ = 90°: this is the in-plane case. k_z = 0, so this diagram "
                "is the whole geometry, not a projection.\n",
                "dim",
            )
        )
    else:
        cos_gamma = float(np.cos(illumination.gamma))
        lines.append(
            Line(
                f"γ = {illumination.gamma_deg:g}°: every ray also carries "
                f"k_z = k cos γ = {cos_gamma:.5f} k out of the page, "
                "identically for all orders. It is perpendicular to this "
                "projection and so cannot appear in it — the 3D view beside "
                "it draws γ as a real angle.\n",
                "dim",
            )
        )

    blaze_angle = getattr(problem.profile, "blaze_angle", None)
    if blaze_angle is None:
        lines.append(
            Line(
                f"{type(problem.profile).__name__} has no blaze angle, so no "
                "blaze direction is drawn.\n",
                "dim",
            )
        )
    else:
        zeta = facet_graze(
            illumination.gamma, np.radians(blaze_angle), illumination.alpha
        )
        lines.append(
            Line(
                f"ζ = {np.degrees(zeta):.3g}° is the graze angle onto the facet "
                "in three dimensions. It is not an angle in this plane and is "
                "not drawn.\n",
                "dim",
            )
        )

    labelled = label_orders(marks)
    if labelled != {m.order for m in marks if m.propagating}:
        lines.append(
            Line(
                f"{len(labelled)} of {live} propagating orders are labelled; "
                "every one is still drawn.\n",
                "dim",
            )
        )

    lines.append(
        Line("Ray lengths are arbitrary; only their directions carry meaning.\n", "dim")
    )
    return tuple(lines)


# -- the whole drawing ----------------------------------------------------


def build(
    problem: "Problem",
    illumination: "Illumination",
    wavelength: float,
    *,
    periods: int = 2,
    pad_orders: int = 2,
) -> Diagram:
    """Everything one drawing consists of. The only entry point a widget needs."""
    period = problem.period
    marks = order_marks(problem, illumination, wavelength, pad=pad_orders)
    labelled = label_orders(marks)
    ray = _RAY_PERIODS * period

    strike = np.array(strike_point(problem))
    paths: list[Polyline] = []
    arrows: list[Arrow] = []
    markers: list[Marker] = []

    # -- the grating surface, over `periods` periods --------------------
    t = np.linspace(-periods, 0.0, _SURFACE_POINTS)
    surface_x = x_nm(t, period)
    surface_y = np.asarray(problem.height_nm(t), dtype=np.float64)
    floor = float(surface_y.min()) - 0.35 * period
    paths.append(
        Polyline(surface_x, surface_y, "surface", fill_to=floor, label="grating")
    )

    # The struck facet, picked out of the surface so the eye finds it.
    start, end = _active_facet_t(problem)
    facet_t = np.linspace(start, end, 64)
    paths.append(
        Polyline(x_nm(facet_t, period), np.asarray(problem.height_nm(facet_t)), "facet")
    )

    # -- incident ray ----------------------------------------------------
    # k_i projects to (-sin a, -cos a): it *travels* that way, so the tail sits
    # back along the opposite direction and the head is the strike point.
    incoming = direction(illumination.alpha)
    tail = strike + ray * incoming
    arrows.append(
        Arrow(
            float(tail[0]), float(tail[1]), float(strike[0]), float(strike[1]),
            "incident", label="incident",
        )
    )
    markers.append(Marker(float(strike[0]), float(strike[1]), "strike", "incident"))

    # -- the diffraction cone, seen end-on -------------------------------
    # Every propagating order has the same projected length, so they all land
    # on this circle. That is what makes the mount's conical nature visible.
    cone_az = np.linspace(-np.pi / 2, np.pi / 2, _CONE_POINTS)
    paths.append(
        Polyline(
            strike[0] + ray * np.sin(cone_az),
            strike[1] + ray * np.cos(cone_az),
            "cone",
            dashed=True,
            label="cone of diffraction",
        )
    )

    # -- one arrow per propagating order ---------------------------------
    for mark in marks:
        if mark.beta is None:
            continue  # no direction; it lives on the ladder instead
        head = strike + ray * direction(mark.beta)
        arrows.append(
            Arrow(
                float(strike[0]), float(strike[1]), float(head[0]), float(head[1]),
                "zero" if mark.order == 0 else "order",
                label=f"m={mark.order:+d}" if mark.order in labelled else None,
                order=mark.order,
            )
        )
        markers.append(Marker(float(head[0]), float(head[1]), "out-of-page", "cone"))

    # -- facet normal and blaze direction --------------------------------
    normal = np.array(facet_normal(problem))
    normal_head = strike + 0.55 * ray * normal
    arrows.append(
        Arrow(
            float(strike[0]), float(strike[1]),
            float(normal_head[0]), float(normal_head[1]),
            "normal", label="facet normal", dashed=True,
        )
    )

    blaze_angle = getattr(problem.profile, "blaze_angle", None)
    if blaze_angle is not None:
        beta_b = blaze_direction(np.radians(blaze_angle), illumination.alpha)
        head = strike + ray * direction(beta_b)
        arrows.append(
            Arrow(
                float(strike[0]), float(strike[1]), float(head[0]), float(head[1]),
                "blaze", label=f"β_b = {np.degrees(beta_b):.1f}°",
            )
        )

    # -- the sin(beta) ladder --------------------------------------------
    # The propagating window, so an evanescent order is visibly outside it
    # rather than merely absent.
    paths.append(
        Polyline(
            np.array([-1.0, 1.0]), np.array([0.0, 0.0]), "axis", panel="ladder",
            label="propagating window",
        )
    )
    for mark in marks:
        tag = "zero" if mark.order == 0 else ("order" if mark.propagating else "evanescent")
        markers.append(
            Marker(mark.sin_beta, 0.0, "strike", tag, panel="ladder",
                   label=f"{mark.order:+d}")
        )

    span = periods * period
    pad = 0.12 * span
    limits = {
        "main": (
            (float(surface_x.min()) - pad, float(surface_x.max()) + pad),
            (floor, float(strike[1]) + ray + pad),
        ),
        "ladder": ((-1.6, 1.6), (-1.0, 1.0)),
    }

    return Diagram(
        paths=tuple(paths),
        arrows=tuple(arrows),
        markers=tuple(markers),
        orders=marks,
        captions=captions(problem, illumination, wavelength, marks),
        limits=limits,
        wavelength=float(wavelength),
        blaze_order=_blaze_order(problem, illumination, wavelength),
    )
