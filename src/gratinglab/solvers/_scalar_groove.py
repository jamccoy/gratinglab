"""The groove-resolved reflection behind :class:`~.scalar.ScalarSolver`.

How much of the groove faces the light, at what local angle, and what the
coating reflects there. :func:`_build_reflection` works that out once per solve,
before the wavelength loop, and :mod:`.scalar` calls the efficiency functions
here from inside it, one per reflectivity model. The validity checks the
solver's provenance reports, and the setup it needs before the loop (the
graze, the coating, the span of orders), are here too.

Nothing here is public: :mod:`.scalar` is the module to read first, and
``docs/theory/scalar.md`` section 9 is the physics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from ..geometry import (
    facet_graze,
    horizon_weights,
    order_range,
    sin_facet_graze,
)
from ..illumination import Illumination
from ..problem import Problem
from ..materials.fresnel import RoughnessModel, amplitude, reflectivity

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..materials import OpticalConstants

#: How reflectivity is resolved across the groove cycle. See
#: :class:`ScalarSolver.solve` and ``docs/theory/scalar.md`` section 9.
ReflectivityModel = Literal["local", "average", "facet"]

#: Which shadows the visibility masks see. ``"facet-normal"`` is the local
#: orientation test alone; ``"horizon"`` adds the shadows one part of the
#: groove casts on another. See :class:`ScalarSolver.solve` and
#: ``docs/theory/scalar.md`` section 9.
Visibility = Literal["facet-normal", "horizon"]


# ============ THE GROOVE'S GEOMETRY, BUILT ONCE PER SOLVE ============


@dataclass
class _Reflection:
    r"""Everything the reflectivity treatment knows about one solve.

    It is built once, before the wavelength loop, because the geometry it
    holds (which parts of the groove face the light, and at what local angle)
    depends on :math:`\alpha` and the profile but not on wavelength or order.
    Only the optical constants vary down the scan.

    It also carries the provenance material, so that the reporting and the
    arithmetic cannot drift apart: ``shadowed_fraction`` and ``suppressed``
    record what the solve did and are not a second estimate of it.
    """

    model: ReflectivityModel
    #: Representative graze for the validity guards, radians.
    graze: float
    basis: str
    visibility: "Visibility" = "facet-normal"
    #: Local graze across the visible groove, radians. Empty when unused.
    local_graze: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    #: sin(zeta) toward the incident direction, per quadrature point, unclipped.
    sin_graze_in: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    tilt: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    #: Incident-side quadrature weights: binary orientation under
    #: "facet-normal", sub-cell shadow boundaries from
    #: :func:`~gratinglab.geometry.horizon_weights` under "horizon". Zero
    #: exactly where the beam cannot reach.
    weight_in: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    #: Grid geometry for the exit-side horizon scans. Empty when unused.
    height_nm: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    period: float = 0.0
    shadowed_fraction: float = 0.0
    #: Orders driven to zero purely by the exit-direction visibility mask.
    suppressed: set[int] = field(default_factory=set)
    brewster: bool = False

    @property
    def visible_in(self) -> NDArray[np.bool_]:
        return self.weight_in > 0.0

    @property
    def masked(self) -> bool:
        """Did visibility masks participate in the solve at all?"""
        return self.sin_graze_in.size > 0

    @property
    def guard_graze(self) -> float:
        """The angle the validity guards should be measured against.

        This is the worst local graze when the model resolves the groove,
        because total external reflection fails first wherever the surface is
        steepest toward the beam, and the single facet angle when the model
        knows only one. A guard should describe the model that ran.
        """
        if self.model == "facet" or self.local_graze.size == 0:
            return self.graze
        return float(self.local_graze.max())


def _facet_tilt(profile, t: NDArray[np.float64]) -> NDArray[np.float64]:
    r"""Local facet tilt :math:`\delta(t)`, radians, from the groove slope.

    :math:`\tan\delta = +\,dy/dt` in normalised units. The sign is the
    essential content of this function. ``docs/findings.md`` records that the
    profile parameter runs against the periodicity direction, which makes the
    opposite sign look equally plausible, but the opposite sign is wrong. The
    deciding check is that an ideal sawtooth must return :math:`\delta` equal
    to its own blaze angle everywhere, so that ``sin_facet_graze`` reproduces
    :func:`~gratinglab.geometry.facet_graze` exactly. At
    :math:`\gamma = 1.25°,\ \alpha = 19.99°` on a 29.5° blaze the correct sign
    gives 1.2328° and the incorrect one gives 0.8119°, and both values are
    plausible angles.

    The profile's analytic slope is used where it exists. ``Blazed`` with a
    vertical anti-blaze facet and ``Lamellar`` both refuse to supply one,
    correctly, since the C-method cannot represent a discontinuity, so this
    function falls back to a periodic central difference on ``height``. That
    is legitimate here although it would not be for the C-method: a vertical
    wall has zero horizontal extent, so it carries zero weight in an integral
    over ``dt``, and the finite-difference smear touches only the two samples
    on either side of the discontinuity and shrinks with the grid.
    """
    from ..profiles import ProfileRepresentationError

    try:
        slope = np.asarray(profile.slope(t), dtype=np.float64)
    except ProfileRepresentationError:
        # Requires the uniform periodic grid `solve` builds; it is the only
        # caller, and a non-uniform grid would silently mis-scale the result.
        step = 1.0 / len(t)
        heights = np.asarray(profile.height(t), dtype=np.float64)
        slope = (np.roll(heights, -1) - np.roll(heights, 1)) / (2.0 * step)
    return np.arctan(slope)


def _build_reflection(
    problem: Problem,
    illumination: Illumination,
    coating: "OpticalConstants | None",
    model: "ReflectivityModel",
    visibility: "Visibility",
    t: NDArray[np.float64],
    height: NDArray[np.float64],
) -> _Reflection:
    """Resolve the groove geometry the reflectivity treatment will use."""
    if model not in ("local", "average", "facet"):
        raise ValueError(
            f"reflectivity_model must be 'local', 'average' or 'facet', got "
            f"{model!r}"
        )
    if visibility not in ("facet-normal", "horizon"):
        raise ValueError(
            f"visibility must be 'facet-normal' or 'horizon', got {visibility!r}"
        )
    if visibility == "horizon" and coating is not None and model == "facet":
        # Accepting this and ignoring the horizon would be the silent wrong
        # answer this project is arranged against.
        raise ValueError(
            "visibility='horizon' cannot act under reflectivity_model='facet': "
            "that model applies one reflectivity to the whole groove and "
            "carries no per-point masks for a horizon to narrow. Use 'local' "
            "or 'average', or keep the default visibility."
        )

    graze, basis = _reflecting_graze(problem, illumination)
    reflection = _Reflection(
        model=model, graze=graze, basis=basis, visibility=visibility
    )
    resolved = coating is not None and model != "facet"
    if not resolved and visibility != "horizon":
        # Nothing to resolve. Without a coating there is no reflectivity at all,
        # and `facet` is the model that declines to look.
        return reflection

    tilt = _facet_tilt(problem.profile, t)
    sin_in = sin_facet_graze(illumination.gamma, tilt, illumination.alpha)

    reflection.tilt = tilt
    reflection.sin_graze_in = sin_in
    if visibility == "horizon":
        # The horizon alone. A back-facing point has a falling ray-adapted
        # height, which puts it under the running horizon, so the orientation
        # test is a subset of this one -- and the horizon's boundaries are
        # geometric crossings `horizon_weights` places at sub-cell precision,
        # where the orientation function just jumps at profile corners.
        weight_in = horizon_weights(
            height, problem.period, illumination.alpha
        )
        reflection.height_nm = height
        reflection.period = problem.period
    else:
        # Binary, because the orientation function jumps at profile
        # corners, where a linear crossing estimate is biased rather than
        # sharp, and bit-for-bit reproducibility is this mode's purpose.
        weight_in = (sin_in > 0.0).astype(np.float64)
    reflection.weight_in = weight_in

    visible = reflection.visible_in
    reflection.shadowed_fraction = float(1.0 - weight_in.mean())
    reflection.local_graze = np.arcsin(np.clip(sin_in[visible], -1.0, 1.0))
    if visible.any():
        reflection.basis = (
            f"groove cycle resolved, "
            f"{np.degrees(reflection.local_graze.min()):.4g}"
            f"-{np.degrees(reflection.local_graze.max()):.4g} deg over the "
            f"{100 * (1 - reflection.shadowed_fraction):.3g}% of the period "
            "that faces the beam"
        )
    else:  # pragma: no cover - needs a groove the beam cannot see at all
        reflection.basis = "groove cycle resolved, fully shadowed"
    return reflection


# ============ EFFICIENCY UNDER EACH REFLECTION MODEL ============


def _groove_average_reflectivity(
    reflection: _Reflection,
    coating: "OpticalConstants",
    wavelength: float,
    roughness_nm: float,
    roughness_model: "RoughnessModel",
) -> float:
    r""":math:`\langle R(\zeta(t))\rangle` over the whole period.

    The average is taken over the whole period, including the unlit part. A
    shadowed facet reflects nothing, and averaging only over the visible part
    would quietly remove the shadowing that this model is intended to capture.
    Under ``"horizon"`` the shadow
    boundary enters at sub-cell precision through the mask weights
    (:func:`~gratinglab.geometry.horizon_weights`).
    """
    weights = reflection.weight_in
    graze = np.arcsin(np.clip(reflection.sin_graze_in, -1.0, 1.0))
    local = reflectivity(
        coating.n(wavelength),
        graze,
        polarization="unpolarized",
        roughness_nm=roughness_nm,
        wavelength_nm=wavelength,
        model=roughness_model,
    )
    return float(np.mean(np.where(weights > 0.0, local, 0.0) * weights))


def _local_reflected_efficiency(
    integrand: NDArray[np.complex128],
    reflection: _Reflection,
    coating: "OpticalConstants",
    wavelength: float,
    betas: NDArray[np.float64],
    orders: NDArray[np.int64],
    illumination: Illumination,
    roughness_nm: float,
    roughness_model: "RoughnessModel",
) -> NDArray[np.float64]:
    r"""Reflection resolved across the groove, carried inside the integral.

    .. math::
        \mathscr{E}_m = \tfrac{1}{2}\left(
            \left|\int r_s(t)\,e^{i\Phi_m(t)}e^{-2\pi i m t}dt\right|^2 +
            \left|\int r_p(t)\,\cdots\right|^2\right)

    A geometric mean is used for the following reason. The obvious local
    coefficient, :math:`r(\zeta_i(t))`, depends on :math:`\alpha` alone, which
    destroys the :math:`\alpha\leftrightarrow\beta_m` symmetry and breaks
    Lorentz reciprocity, by 82% as measured on the reference geometry.
    Weighting a facet by both the angle at which it receives and the angle at
    which it emits,

    .. math:: r^{\text{gm}}(t) = \sqrt{r(\zeta_i(t))}\,\sqrt{r(\zeta_{d,m}(t))}

    is symmetric by construction, holds reciprocity to :math:`6\times10^{-15}`,
    and reduces to plain :math:`r(\zeta)` wherever the facet is at specular,
    so a perfectly blazed groove is unaffected by the change.

    The square roots are taken separately, instead of as
    :math:`\sqrt{r_i r_d}`, so that the only discontinuities left in the weight
    are the ones :math:`r` itself has.

    The benefit of this was measured, and it is smaller than it first
    appears. Where the argument of the product wraps uniformly across the
    groove, which is the grazing case because :math:`\arg r` is then nearly
    the same at every :math:`t`, the two forms differ by a global sign, and a
    global sign cancels out of :math:`|\int\cdots|^2`. On Au at 1-6 nm they
    agree to :math:`10^{-15}`.

    The two forms differ only where the wrap is non-uniform, near normal
    incidence on a deep groove, where :math:`\arg r_p` sweeps through the
    Brewster jump and different parts of the same groove land on opposite sides
    of the cut. There the difference reaches 12%, on orders carrying
    :math:`10^{-9}`, and in a regime that `_brewster_crossing` already reports
    as beyond what this symmetrisation can carry.

    Separate roots are therefore an inexpensive defensive choice and are not
    essential. This is stated explicitly because the first draft of this
    docstring claimed that the product wrapped in the grazing regime and that
    the separate roots were therefore essential, which the measurement does
    not support.
    """
    n_c = coating.n(wavelength)
    common = dict(
        roughness_nm=roughness_nm, wavelength_nm=wavelength, model=roughness_model
    )

    graze_in = np.arcsin(np.clip(reflection.sin_graze_in, -1.0, 1.0))
    r_s_in, r_p_in = amplitude(n_c, graze_in, **common)

    sin_out = sin_facet_graze(
        illumination.gamma, reflection.tilt[None, :], betas[:, None]
    )
    # Minimum, not product: the incident and exit shadows share boundaries
    # (every direction's shadow starts at the same apex corner), and where
    # two same-side boundaries land in one cell the lit overlap is the
    # smaller of the two sub-cell weights -- multiplying them would count
    # the shared boundary twice and put an O(1/n) error on the blaze order.
    coverage = np.minimum(
        reflection.weight_in[None, :], _exit_weights(reflection, betas, sin_out)
    )
    visible = coverage > 0.0
    graze_out = np.arcsin(np.clip(sin_out, -1.0, 1.0))
    r_s_out, r_p_out = amplitude(n_c, graze_out, **common)

    zero = np.zeros((), dtype=np.complex128)
    w_s = np.where(visible, np.sqrt(r_s_in)[None, :] * np.sqrt(r_s_out), zero) * coverage
    w_p = np.where(visible, np.sqrt(r_p_in)[None, :] * np.sqrt(r_p_out), zero) * coverage

    # An order no facet can radiate into is zero for a geometric reason, not
    # because it passed off. Record it so the provenance can say which.
    reflection.suppressed.update(
        int(m) for m, lit in zip(orders, visible.any(axis=1)) if not lit
    )
    reflection.brewster |= _brewster_crossing(r_p_in[reflection.visible_in])

    return 0.5 * (
        np.abs(np.mean(w_s * integrand, axis=1)) ** 2
        + np.abs(np.mean(w_p * integrand, axis=1)) ** 2
    )


def _exit_weights(
    reflection: _Reflection,
    betas: NDArray[np.float64],
    sin_out: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Exit-side mask weights, one row per order: orientation, and under
    ``"horizon"`` the cast shadows toward each diffracted direction.

    These are the same scans the incident side ran, at each
    :math:`\\beta_m`. Occlusion along a straight ray reads the same from
    either end, and this is what keeps the horizon model reciprocal.
    """
    if reflection.visibility == "horizon":
        # The horizon alone, for the same two reasons as the incident side:
        # it subsumes the orientation test, and its boundaries are geometric
        # crossings rather than jumps at profile corners.
        return np.stack(
            [
                horizon_weights(
                    reflection.height_nm, reflection.period, float(b)
                )
                for b in betas
            ]
        )
    return (sin_out > 0.0).astype(np.float64)


def _shadow_masked_efficiency(
    integrand: NDArray[np.complex128],
    reflection: _Reflection,
    betas: NDArray[np.float64],
    orders: NDArray[np.int64],
    illumination: Illumination,
) -> NDArray[np.float64]:
    r"""The visibility masks at unit amplitude: shadowing without a material.

    This is the uncoated ``visibility="horizon"`` branch. It uses the same
    geometry as the ``"local"`` model (a point contributes only if it faces,
    and is not cast-shadowed from, both the incident and the diffracted
    direction), but with :math:`r \equiv 1` in place of a Fresnel amplitude.
    The result is a perfect reflector with shadows, which is the
    configuration needed for a perfect-conductor cross-check against the
    integral solver on grooves deep enough to shadow themselves.
    """
    sin_out = sin_facet_graze(
        illumination.gamma, reflection.tilt[None, :], betas[:, None]
    )
    # Minimum, not product -- see `_local_reflected_efficiency`.
    coverage = np.minimum(
        reflection.weight_in[None, :], _exit_weights(reflection, betas, sin_out)
    )
    visible = coverage > 0.0

    reflection.suppressed.update(
        int(m) for m, lit in zip(orders, visible.any(axis=1)) if not lit
    )

    return np.abs(np.mean(np.where(visible, integrand, 0.0) * coverage, axis=1)) ** 2


# ============ CHECKS THE PROVENANCE REPORTS ============


def _sp_split(
    coating: "OpticalConstants",
    graze: float,
    wavelengths: NDArray[np.float64],
) -> float:
    r"""Worst fractional s/p intensity split, :math:`|R_s - R_p| / R_s`.

    The split is evaluated only at the guarded graze and at the extremes of
    the scan's wavelengths. It grows monotonically with graze in this regime,
    so the angle the validity guards already watch is the right single point,
    and two wavelengths bound a table that varies smoothly between them.
    Roughness is
    left out because the damping factors apply to s and p alike and
    would cancel from the ratio.
    """
    split = 0.0
    for wavelength in (float(wavelengths.min()), float(wavelengths.max())):
        r_s, r_p = amplitude(coating.n(wavelength), graze)
        big, small = np.abs(r_s) ** 2, np.abs(r_p) ** 2
        if big > 0.0:
            split = max(split, float(np.abs(big - small) / big))
    return split


def _brewster_crossing(r_p: NDArray[np.complex128]) -> bool:
    r"""Does :math:`r_p` pass through zero somewhere on the lit groove?

    At Brewster's angle the p amplitude changes sign, so its phase jumps by
    :math:`\pi` and the half-phase the geometric mean takes jumps by
    :math:`\pi/2`. The jump is real physics, but the geometric-mean
    symmetrisation is not built to carry it, so the result is reported as
    suspect instead of being returned without comment.

    The condition is reachable only with steep grooves near normal incidence.
    At the grazing angles this project is aimed at, :math:`r_p` stays well
    away from zero.
    """
    if r_p.size == 0:
        return False
    magnitude = np.abs(r_p)
    return bool(magnitude.min() < 0.05 * magnitude.max())


# ============ SETUP: GRAZE, COATING AND ORDERS ============


def _reflecting_graze(problem: Problem, illumination: Illumination):
    r"""``(zeta, description)``: the angle at which the surface reflects.

    A blazed profile has one flat active facet tilted by ``blaze_angle``, and
    :math:`\sin\zeta = \sin\gamma\,\cos(\delta - \alpha)` is exact for it.

    No other profile has a single facet angle, so the reflection is evaluated
    on the mean surface, which is `facet_graze` with a zero tilt and requires
    no separate formula:

    >>> facet_graze(gamma, 0.0, alpha)          # doctest: +SKIP
    ...                                          # == arcsin(|k_i . n|)

    The accuracy of this depends on the profile, and the description records
    which case applies:

    - Lamellar tops and bottoms are parallel to the mean surface, so the
      result is exact for them. It omits the vertical walls, which at grazing
      incidence are nearly edge-on.
    - Trapezoidal has the same flat tops and bottoms, but its walls are tilted
      and so are neither edge-on nor parallel to the mean surface. The result
      is exact on the lands and ignores the walls, and says so.
    - Sinusoidal and measured profiles have a local slope that varies across
      the groove, so the mean is an approximation and is recorded as one.

    A measured profile is the one case that can avoid this. If it arrived from
    metrology with a fitted ``blaze_angle``, it takes the exact facet branch
    above, marked ``fitted`` instead of ``declared``. Without one it falls
    back to the mean surface, which for a sawtooth at grazing incidence is
    wrong by the whole blaze angle, so the fallback is a real loss of
    accuracy.

    The alternative, refusing any profile without a facet angle, would be
    worse than approximating. A sinusoid at grazing incidence does reflect,
    and a solver that declined to estimate how much would be less useful
    without being more accurate.
    """
    kind = type(problem.profile).__name__

    blaze_angle = getattr(problem.profile, "blaze_angle", None)
    if blaze_angle is not None:
        # A measured profile may carry a fitted facet angle. The geometry is
        # used identically either way -- but the reader is entitled to know
        # whether the tilt was specified or estimated, because only one of
        # those has an uncertainty attached to it.
        origin = "fitted" if kind == "FromProfileData" else "declared"
        return (
            facet_graze(illumination.gamma, np.radians(blaze_angle), illumination.alpha),
            f"active facet, {blaze_angle:g} deg tilt ({origin})",
        )

    caveat = {
        "Lamellar": "exact for its flat tops and bottoms, omits the vertical walls",
        "Trapezoidal": "exact for its flat tops and bottoms, omits the sloped "
        "sidewalls",
    }.get(
        kind,
        "the local slope varies across the groove, so this is an approximation",
    )
    return (
        facet_graze(illumination.gamma, 0.0, illumination.alpha),
        f"mean surface ({kind} has no single facet angle; {caveat})",
    )


def _resolve_coating(problem: Problem) -> "OpticalConstants | None":
    """Turn ``Problem.coating`` into optical constants, or refuse.

    ``coating`` is a material name, so that a benchmark case remains a small
    serialisable file instead of embedding a 500-row table. Resolving it here,
    and not in ``Problem.__init__``, keeps a saved case loadable on a machine
    whose vendored tables have been trimmed, since the name has to resolve
    only when someone solves.

    An unknown name raises an error. The string was formerly read by
    nothing, so ``coating="unobtanium"`` silently relabelled a result as
    absolute. A name that cannot be resolved poses a question this solver
    cannot answer, and the solver does not substitute a default for it.
    """
    if problem.coating is None:
        return None
    from ..materials import lookup

    return lookup(problem.coating)


def _spanning_orders(
    wavelengths: NDArray[np.float64],
    period: float,
    sin_alpha: float,
    sin_gamma: float,
) -> NDArray[np.int64]:
    """Every order propagating at any wavelength in the scan, as one range.

    The shortest wavelength admits the most orders, so it bounds the set.
    """
    widest = order_range(float(wavelengths.min()), period, sin_alpha, sin_gamma)
    return widest
