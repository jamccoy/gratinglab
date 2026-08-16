r"""Scalar (Kirchhoff / thin-element) diffraction efficiency.

Implemented as the **general Fourier integral**, not as the closed-form sawtooth
result:

.. math::
    G_m = \int_0^1 e^{i\Phi_m(t)}\, e^{-2\pi i m t}\,dt, \qquad
    \Phi_m(t) = k\,g(t)\,\sin\gamma\,\left[\cos\alpha + \cos\beta_m\right]

where ``t`` is position normalised to the period and ``g(t)`` is the groove
height in nm. Efficiency is :math:`\mathscr{E}_m = |G_m|^2`.

Two consequences of doing it this way are worth stating plainly.

**It accepts any profile**, including a measured AFM boundary, which the closed
forms cannot. That is what lets a scalar calculation and a PCGrate run be driven
from the identical geometry.

**The closed forms become tests, not code paths.** The Appendix D results --
sawtooth ``sinc²(φ/2 − mπ)``, square wave, sinusoid -- are exact answers we
already trust, so they check the numerics rather than duplicating them.

**The phase depends on the order** through :math:`\cos\beta_m` (ISSI eq. 15;
thesis Appendix-D.tex:651), so this is one integral per order rather than a
single FFT. That choice is deliberate and carries a known cost, because no
scalar formulation can satisfy all three of the following at once:

==================  ==================================  =========================
Property            order-dependent (this)              beta-free alternative
==================  ==================================  =========================
Energy, sum <= 1    violated, up to 1.61 measured       exact (Parseval)
Reciprocity         exact, 1e-17                        violated, 0.44 measured
Blaze direction     exact at every alpha                exact at Littrow only
==================  ==================================  =========================

The conflict is structural: reciprocity *requires* the phase to be symmetric
under :math:`\alpha \leftrightarrow \beta_m`, which forces both angles in;
Parseval *requires* a single fixed function, which forbids :math:`\beta_m`.
Physical optics is reciprocal but not energy-conserving; the pure
transmittance-function picture is the reverse.

**Reciprocity is the invariant this solver keeps.** The energy deviation is
therefore expected, is reported rather than hidden, and is never renormalised
away. It is not an implementation defect: in the shallow-groove limit the sum
is 1.00000 exactly, degrading smoothly as groove depth grows (0.0001 -> 1.00000,
0.10 -> 0.92902). The deviation scales with **phase excursion across the
groove** -- depth relative to wavelength, hence working diffraction order --
not with lambda/period or the number of propagating orders.

``docs/theory/scalar.md`` §5 carries the full derivation and evidence;
``docs/conventions.md`` §9 is the reference for every formula here.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..geometry import (
    beta,
    cos_beta,
    facet_graze,
    flux_obliquity,
    is_propagating,
    order_range,
    sin_beta,
    sin_facet_graze,
)
from ..illumination import Illumination
from ..problem import Problem
from ..result import EfficiencyScan, Provenance
from ..materials.fresnel import RoughnessModel, amplitude, reflectivity
from .base import Capabilities, Progress, register

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..materials import OpticalConstants

__all__ = ["ReflectivityModel", "ScalarSolver", "interference_factor", "scalar"]

#: How reflectivity is resolved across the groove cycle. See
#: :class:`ScalarSolver.solve` and ``docs/theory/scalar.md`` section 9.
ReflectivityModel = Literal["local", "average", "facet"]

#: Above this ratio of wavelength to period, scalar theory is losing validity.
#: Soft X-ray work runs ~0.005; visible gratings run ~0.4.
_LAMBDA_OVER_PERIOD_WARN = 0.1


def interference_factor(s: ArrayLike, n_grooves: int) -> NDArray[np.float64]:
    r"""The finite-:math:`N` grating interference function, ISSI eq. (8).

    .. math::
        \left[\frac{\sin(Ns)}{N\sin(s)}\right]^2, \qquad
        s = \left[\sin\alpha + \sin\beta\right]\sin\gamma\,\frac{p\pi}{\lambda}

    **This is deliberately not applied to order efficiencies.** At an exact
    diffraction order ``s = mπ`` and the factor equals 1, so multiplying by it
    would change nothing. Its role is the *angular line shape* around each
    order -- the link to resolving power -- which is a different question from
    how much power lands in order ``m``. Exposed here for that purpose.
    """
    s = np.asarray(s, dtype=np.float64)
    numerator = np.sin(n_grooves * s)
    denominator = n_grooves * np.sin(s)
    # s -> m*pi is removable: the ratio tends to +/-1, so the square tends to 1.
    near_order = np.isclose(np.sin(s), 0.0, atol=1e-12)
    ratio = np.where(near_order, 1.0, numerator / np.where(near_order, 1.0, denominator))
    return ratio**2


class ScalarSolver:
    """Kirchhoff scalar diffraction from a surface-relief grating."""

    capabilities = Capabilities(
        name="scalar",
        conical=True,
        # Scalar theory neglects polarization, but "the TE efficiency" is still
        # a question it can answer -- the answer just happens to equal the TM
        # one. That is an approximation to record, not an unsupported
        # configuration to refuse, so every polarization is accepted and the
        # degeneracy is reported on the provenance instead.
        polarizations=("TE", "TM", "unpolarized"),
        accuracy_knob="quadrature_points",
        rigorous=False,
        handles_undercut=False,
        reports_progress=True,
    )

    def solve(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: ArrayLike,
        *,
        quadrature_points: int = 2048,
        roughness_model: "RoughnessModel" = "nevot-croce",
        reflectivity_model: "ReflectivityModel" = "local",
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        r"""Efficiency over a wavelength scan.

        Efficiency is :math:`\mathscr{E}_m = O_m |G_m|^2` -- the norm-squared
        Fourier coefficient times the **symmetric** flux obliquity
        :math:`4\cos\alpha\cos\beta_m/(\cos\alpha+\cos\beta_m)^2`, and **no
        renormalisation**. Thesis Appendix D carries the *asymmetric* obliquity
        :math:`\cos\beta_m/\cos\alpha`, which breaks reciprocity, and a
        renormalisation to :math:`\sum_m \mathscr{E}_m = 1`, which destroys the
        model's own error signal; both of those are removed. See
        :func:`~gratinglab.geometry.flux_obliquity` for why a flux factor
        nonetheless belongs, and ``tests/test_perturbation.py`` for the
        independent theory that fixes which one.

        Parameters
        ----------
        quadrature_points
            Samples per period for the Fourier integral. The integrand is
            periodic, so the rectangle rule converges spectrally and the
            default is generous. Must exceed twice the highest order to satisfy
            Nyquist; this is checked.
        roughness_model
            How ``Problem.roughness`` damps the reflectivity: ``"nevot-croce"``
            (default), ``"debye-waller"``, or ``"none"``. Both are legitimate
            and they differ, so the choice is exposed rather than buried --
            Debye-Waller is the more pessimistic near the critical angle
            because it knows nothing about the transmitted wave. Ignored
            without a coating, since there is no reflectivity to damp.
        reflectivity_model
            How reflectivity is resolved **across the groove cycle**. Ignored
            without a coating.

            ``"local"`` (default) evaluates the complex Fresnel amplitude at
            every quadrature point, from the local facet tilt, and carries it
            *inside* the diffraction integral. A groove whose reflectivity
            varies across the cycle is an amplitude grating as well as a phase
            grating, so this makes reflectivity **order-dependent** -- see the
            note below, which corrects the module docstring's earlier claim
            that no such mechanism exists.

            ``"average"`` takes the groove-cycle mean of the *intensity*,
            :math:`\langle R(\zeta(t))\rangle`, and applies it as one factor
            per wavelength. It sees the shadowing and the varying local angle
            but not the interference between them, so it stays
            order-independent -- and, like ``"facet"``, it **breaks
            reciprocity**, because its :math:`\zeta(t)` is built from
            :math:`\alpha` alone. Resolving the groove is not what repairs
            reciprocity; symmetrising in the exit direction is, and only
            ``"local"`` does that. Offered as the diagnostic that separates
            shadowing from interference, not as a physical model.

            ``"facet"`` is the M15 treatment: one :math:`R` per wavelength at
            the active-facet angle, applied to every order alike. Kept so a
            prior result can be reproduced exactly, and so the size of the
            change is measurable rather than asserted.

            Both resolved models zero the reflection where a facet is
            back-facing, in either direction. That is honest geometry and it
            can drive a weak order to exactly zero; when it does, the order is
            named in a provenance warning rather than left to look like a
            passing-off.
        progress
            Called ``(0, n)`` before the first wavelength and ``(k, n)`` after
            each one, per the contract in
            :class:`~gratinglab.solvers.base.Solver`. Exceptions from it are
            **not** caught, which is what lets a caller raise
            :class:`~gratinglab.solvers.base.SolveCancelled` and actually stop
            this loop. A scan cancelled that way returns nothing at all: a
            partially filled efficiency array with the remaining rows still
            zero is indistinguishable from a scan whose orders all passed off,
            and handing one back would be the silent-wrong-number failure this
            project is arranged against.

        Notes
        -----
        Summed efficiency will not equal 1, and can exceed it. That is expected
        -- see the module docstring for the three-way tension that makes it
        unavoidable, and why reciprocity is the invariant kept instead. The sum
        is **reported, never rescaled**: renormalising by it (as Appendix D
        does) would erase the model's own signal about how far it has strayed
        from its regime.
        """
        self.capabilities.check(problem, illumination)
        # Before any work. An unknown material must fail loudly rather than
        # produce a scan that quietly means something else -- see
        # `_resolve_coating`.
        coating = _resolve_coating(problem)

        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
        if wavelengths.ndim != 1:
            raise ValueError("wavelengths must be one-dimensional")
        if (wavelengths <= 0).any():
            raise ValueError("wavelengths must be positive")
        if quadrature_points < 16:
            raise ValueError(
                f"quadrature_points must be at least 16, got {quadrature_points}"
            )

        started = time.perf_counter()

        sin_alpha = illumination.sin_alpha
        cos_alpha = illumination.cos_alpha
        sin_gamma = illumination.sin_gamma

        # One column per order, spanning every order that propagates anywhere
        # in the scan, so column j means order j at every wavelength.
        all_orders = _spanning_orders(
            wavelengths, problem.period, sin_alpha, sin_gamma
        )
        max_order = int(np.abs(all_orders).max())
        if quadrature_points <= 2 * max_order:
            raise ValueError(
                f"quadrature_points={quadrature_points} cannot resolve order "
                f"{max_order}; Nyquist needs more than {2 * max_order}"
            )

        t = np.linspace(0.0, 1.0, quadrature_points, endpoint=False)
        height = np.asarray(problem.height_nm(t), dtype=np.float64)
        # exp(-2*pi*i*m*t) for every order, reused at every wavelength.
        kernel = np.exp(-2j * np.pi * all_orders[:, None] * t[None, :])

        reflection = _build_reflection(
            problem, illumination, coating, reflectivity_model, t
        )
        resolved = coating is not None and reflection.model != "facet"

        efficiency = np.zeros((len(wavelengths), len(all_orders)))
        propagating = np.zeros_like(efficiency, dtype=bool)
        # One scale factor per wavelength for the two order-independent models.
        scale = np.ones(len(wavelengths))

        for row, wavelength in enumerate(wavelengths):
            if progress is not None:
                # At the *top* of the row it is about to compute, so `row` is
                # the count already finished. That makes the first call
                # (0, n) -- before any work, as the contract requires -- and
                # puts a cancellation point ahead of every wavelength rather
                # than behind it. The `continue` below cannot skip it.
                progress(row, len(wavelengths))

            sines = sin_beta(
                all_orders, wavelength, problem.period, sin_alpha, sin_gamma
            )
            live = is_propagating(sines)
            if not live.any():
                continue

            cosines = cos_beta(sines[live])
            # Phi_m(t) = k g(t) sin(gamma) [cos(alpha) + cos(beta_m)]
            #
            # Symmetric under alpha <-> beta_m, which is exactly what makes the
            # result reciprocal -- and exactly what stops the G_m from being
            # Fourier coefficients of any single function. See the module
            # docstring; that trade is deliberate.
            phase = (2.0 * np.pi / wavelength) * height * sin_gamma
            phi = phase[None, :] * (cos_alpha + cosines)[:, None]

            integrand = np.exp(1j * phi) * kernel[live]
            # Norm squared of the coefficient, times the flux projection.
            # Still no renormalisation by the sum -- that remains an error.
            #
            # The obliquity factor is *symmetric*, 4 c_a c_b / (c_a + c_b)^2,
            # not the c_b / c_a of thesis Appendix D. That distinction is the
            # whole point: the asymmetric form breaks reciprocity, the
            # symmetric one does not, and only the symmetric one reproduces
            # first-order perturbation theory in the shallow limit. See
            # `geometry.flux_obliquity` and `tests/test_perturbation.py`.
            obliquity = flux_obliquity(cos_alpha, cosines)

            if resolved and reflection.model == "local":
                values = obliquity * _local_reflected_efficiency(
                    integrand,
                    reflection,
                    coating,
                    wavelength,
                    beta(sines[live]),
                    all_orders[live],
                    illumination,
                    problem.roughness,
                    roughness_model,
                )
            else:
                values = obliquity * np.abs(np.mean(integrand, axis=1)) ** 2
                if resolved:  # "average"
                    scale[row] = _groove_average_reflectivity(
                        reflection,
                        coating,
                        wavelength,
                        problem.roughness,
                        roughness_model,
                    )

            efficiency[row, live] = values
            propagating[row, live] = True

        if progress is not None:
            progress(len(wavelengths), len(wavelengths))

        if coating is not None and reflection.model == "facet":
            # The M15 treatment, kept verbatim: one factor per wavelength,
            # applied to every order alike. It is an approximation whose size is
            # now measurable rather than assumed -- `reflectivity_model="local"`
            # is the same run with the groove resolved.
            scale = reflectivity(
                coating.n(wavelengths),
                reflection.graze,
                # Unpolarized deliberately, whatever the illumination says.
                # These s and p are facet-local and `Illumination.polarization`
                # is groove-referenced (conventions.md 7) -- different frames,
                # and in a conical mount the mapping is not the identity.
                # Resolving polarization here would be false precision on a
                # model that already reports TE and TM as identical, which the
                # solver warns about a few lines below.
                polarization="unpolarized",
                roughness_nm=problem.roughness,
                wavelength_nm=wavelengths,
                model=roughness_model,
            )
        efficiency = efficiency * np.asarray(scale)[:, None]

        return EfficiencyScan(
            wavelengths=wavelengths,
            orders=all_orders,
            efficiency=efficiency,
            propagating=propagating,
            provenance=self._provenance(
                problem,
                illumination,
                wavelengths,
                quadrature_points,
                efficiency,
                propagating,
                time.perf_counter() - started,
                coating,
                reflection,
            ),
        )

    def _provenance(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: NDArray[np.float64],
        quadrature_points: int,
        efficiency: NDArray[np.float64],
        propagating: NDArray[np.bool_],
        elapsed: float,
        coating: "OpticalConstants | None" = None,
        reflection: "_Reflection | None" = None,
    ) -> Provenance:
        """Record the run, including every validity guard it tripped."""
        from .. import __version__

        if reflection is None:  # pragma: no cover - defensive
            reflection = _Reflection(model="facet", graze=0.0, basis="none")
        graze = reflection.guard_graze
        warnings: list[str] = []

        ratio = float(wavelengths.max() / problem.period)
        if ratio > _LAMBDA_OVER_PERIOD_WARN:
            warnings.append(
                f"lambda/period reaches {ratio:.3g}; scalar theory assumes "
                f"lambda << period and loses accuracy above "
                f"{_LAMBDA_OVER_PERIOD_WARN}"
            )

        if problem.roughness > 0:
            # Fraunhofer smoothness criterion, ISSI section 4. Uses the same
            # zeta the reflectivity does, rather than recomputing it for blazed
            # profiles only -- an optically rough surface is rough whatever
            # shape its grooves are, and this used to skip every profile
            # without a `blaze_angle`.
            threshold = 32.0 * np.sin(graze) * problem.roughness
            if wavelengths.min() < threshold:
                warnings.append(
                    f"Fraunhofer smoothness criterion violated below "
                    f"{threshold:.4g} nm (32 sin(zeta) sigma); the surface is "
                    "not optically smooth there"
                )

        if coating is not None:
            # Total external reflection, at last. `docs/theory/scalar.md`
            # section 7 has carried this row as "needs a materials layer to
            # evaluate" since it was written; the layer is here.
            #
            # Only with a coating: theta_c comes from the decrement, so without
            # optical constants there is nothing to compare against. That is
            # not a validity concern to warn about, it is a check that does not
            # apply -- the M8 distinction.
            critical = coating.critical_angle(wavelengths)
            past = graze > critical
            if past.any():
                lam = wavelengths[past]
                warnings.append(
                    f"facet graze {np.degrees(graze):.4g} deg exceeds the "
                    f"critical angle for {coating.name} over "
                    f"{lam.min():.4g}-{lam.max():.4g} nm "
                    f"(theta_c falls to {np.degrees(critical[past].min()):.4g} "
                    "deg there); reflectivity collapses above it, so these "
                    "wavelengths are outside the regime a grazing-incidence "
                    "design operates in -- see docs/theory/scalar.md section 7"
                )

            if reflection.suppressed:
                # An order at exactly zero must never be mistaken for one that
                # passed off. Naming it is the whole difference between a
                # modelled result and a silent one.
                named = ", ".join(str(m) for m in sorted(reflection.suppressed))
                warnings.append(
                    f"order(s) {named} are zero at one or more wavelengths in "
                    f"this scan because no part of the groove faces both the "
                    f"incident and the diffracted direction; "
                    f"{100 * reflection.shadowed_fraction:.3g}% of the period "
                    "is shadowed from the beam to begin with. This is geometric "
                    "shadowing by the local facet normal, not passing-off -- "
                    "`propagating` still marks these orders live. Scalar theory "
                    "has no diffraction into a shadowed direction, so zero is "
                    "the model's answer rather than the physical one"
                )

            if reflection.brewster:
                warnings.append(
                    "the p-polarized amplitude passes through zero somewhere on "
                    "the lit groove (Brewster), where its phase jumps by pi; the "
                    "groove-resolved reflectivity model symmetrises with a "
                    "geometric mean and does not carry that jump cleanly. Use "
                    "reflectivity_model='average' or 'facet' here, or read this "
                    "result as indicative -- see docs/theory/scalar.md section 9"
                )

        # No coating is the normal default mode, not a validity concern -- it
        # is reported via notes["normalization"] below, not as a warning.
        # There is nothing wrong with a run that has not been given a coating;
        # warnings are reserved for cases where the model is being pushed
        # outside conditions it can actually answer for.

        if illumination.polarization != "unpolarized":
            warnings.append(
                f"result is labelled {illumination.polarization} but scalar "
                "theory neglects polarization; TE and TM are identical here"
            )

        # Report, never rescale. How far an approximate theory strays from a
        # conservation law is information, not something to hide -- and
        # renormalising by the sum (thesis Appendix D) would destroy exactly
        # that information.
        totals = np.where(propagating, efficiency, 0.0).sum(axis=1)
        deviation = float(np.abs(totals - 1.0).max())
        if deviation > 0.01:
            worst = float(totals[np.argmax(np.abs(totals - 1.0))])
            direction = "above" if worst > 1.0 else "below"
            warnings.append(
                f"summed efficiency reaches {worst:.4f}, {100 * deviation:.1f}% "
                f"{direction} unity. Expected: the thin-element approximation "
                "for a reflection grating does not conserve energy, because the "
                "phase depends on the exit direction and the coefficients are "
                "therefore not a Parseval pair. Kept that way because the "
                "energy-conserving alternative violates Lorentz reciprocity. "
                "The deviation grows with phase excursion across the groove "
                "(depth/wavelength, hence working order), and vanishes in the "
                "shallow-groove limit -- see docs/theory/scalar.md section 5"
                + (
                    ". Note this run is absolute, so part of the deficit is "
                    "ordinary absorption in the coating rather than the "
                    "approximation straying"
                    if coating is not None and worst < 1.0
                    else ""
                )
            )

        return Provenance(
            method="scalar",
            version=__version__,
            truncation=quadrature_points,
            # A closed-form quadrature that has not been swept is still not a
            # demonstrated convergence. The harness sets this.
            converged=None,
            wall_time_s=elapsed,
            warnings=tuple(warnings),
            notes={
                # Keyed on whether a reflectivity was actually applied, never
                # on whether a coating was *named*. `Problem.coating` used to
                # be a free-form string nothing read, so setting it to anything
                # at all relabelled a result as "absolute" while leaving every
                # number unchanged. M15-D is what makes this ever say
                # "absolute"; until then a named coating is recorded and not
                # yet used, which is what this reports.
                "normalization": "absolute" if coating is not None else "relative",
                **(
                    {
                        "coating": f"{coating.name} ({coating.source})",
                        "reflectivity_model": reflection.model,
                        "reflectivity_graze": reflection.basis,
                        "shadowed_fraction": reflection.shadowed_fraction,
                    }
                    if coating is not None
                    else {}
                ),
                "alpha_deg": illumination.alpha_deg,
                "gamma_deg": illumination.gamma_deg,
            },
        )


@dataclass
class _Reflection:
    r"""Everything the reflectivity treatment knows about one solve.

    Built once, before the wavelength loop, because the geometry it holds --
    which parts of the groove face the light, and at what local angle -- depends
    on :math:`\alpha` and the profile but not on wavelength or order. Only the
    optical constants vary down the scan.

    Carries the provenance material too, so the reporting and the arithmetic
    cannot drift apart: ``shadowed_fraction`` and ``suppressed`` are what the
    solve actually did, not a second estimate of it.
    """

    model: ReflectivityModel
    #: Representative graze for the validity guards, radians.
    graze: float
    basis: str
    #: Local graze across the visible groove, radians. Empty when unused.
    local_graze: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    #: sin(zeta) toward the incident direction, per quadrature point, unclipped.
    sin_graze_in: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    tilt: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    shadowed_fraction: float = 0.0
    #: Orders driven to zero purely by the exit-direction visibility mask.
    suppressed: set[int] = field(default_factory=set)
    brewster: bool = False

    @property
    def visible_in(self) -> NDArray[np.bool_]:
        return self.sin_graze_in > 0.0

    @property
    def guard_graze(self) -> float:
        """The angle the validity guards should be measured against.

        The **worst** local graze when the model resolves the groove, because
        total external reflection fails first wherever the surface is steepest
        toward the beam; the single facet angle when the model only knows one.
        A guard should describe the model that ran.
        """
        if self.model == "facet" or self.local_graze.size == 0:
            return self.graze
        return float(self.local_graze.max())


def _facet_tilt(profile, t: NDArray[np.float64]) -> NDArray[np.float64]:
    r"""Local facet tilt :math:`\delta(t)`, radians, from the groove slope.

    :math:`\tan\delta = +\,dy/dt` in normalised units. **The sign is the whole
    content of this function.** ``docs/findings.md`` records that the profile
    parameter runs against the periodicity direction, which makes the opposite
    sign look equally plausible; it is wrong. The check that settles it is that
    an ideal sawtooth must come back with :math:`\delta` equal to its own blaze
    angle everywhere, so that ``sin_facet_graze`` reproduces
    :func:`~gratinglab.geometry.facet_graze` exactly. At
    :math:`\gamma = 1.25°,\ \alpha = 19.99°` on a 29.5° blaze the right sign
    gives 1.2328° and the wrong one gives 0.8119°, and both look like angles.

    Uses the profile's analytic slope where it has one. ``Blazed`` with a
    vertical anti-blaze facet and ``Lamellar`` both *refuse* to supply one --
    correctly, since the C-method cannot represent a discontinuity -- so this
    falls back to a periodic central difference on ``height``. That is
    legitimate here in a way it would not be there: a vertical wall has zero
    horizontal extent, so it carries zero weight in an integral over ``dt``, and
    the finite-difference smear touches the two samples either side of the
    discontinuity and shrinks with the grid.
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
    t: NDArray[np.float64],
) -> _Reflection:
    """Resolve the groove geometry the reflectivity treatment will use."""
    if model not in ("local", "average", "facet"):
        raise ValueError(
            f"reflectivity_model must be 'local', 'average' or 'facet', got "
            f"{model!r}"
        )

    graze, basis = _reflecting_graze(problem, illumination)
    reflection = _Reflection(model=model, graze=graze, basis=basis)
    if coating is None or model == "facet":
        # Nothing to resolve. Without a coating there is no reflectivity at all,
        # and `facet` is the model that deliberately declines to look.
        return reflection

    tilt = _facet_tilt(problem.profile, t)
    sin_in = sin_facet_graze(illumination.gamma, tilt, illumination.alpha)
    visible = sin_in > 0.0

    reflection.tilt = tilt
    reflection.sin_graze_in = sin_in
    reflection.shadowed_fraction = float(1.0 - visible.mean())
    reflection.local_graze = np.arcsin(np.clip(sin_in[visible], -1.0, 1.0))
    reflection.basis = (
        f"groove cycle resolved, {np.degrees(reflection.local_graze.min()):.4g}"
        f"-{np.degrees(reflection.local_graze.max()):.4g} deg over the "
        f"{100 * (1 - reflection.shadowed_fraction):.3g}% of the period that "
        "faces the beam"
    )
    return reflection


def _groove_average_reflectivity(
    reflection: _Reflection,
    coating: "OpticalConstants",
    wavelength: float,
    roughness_nm: float,
    roughness_model: "RoughnessModel",
) -> float:
    r""":math:`\langle R(\zeta(t))\rangle` over the whole period.

    Over the **whole** period, not just the lit part: a shadowed facet reflects
    nothing, and averaging only over what is visible would quietly delete the
    shadowing this model exists to see.
    """
    visible = reflection.visible_in
    graze = np.arcsin(np.clip(reflection.sin_graze_in, -1.0, 1.0))
    local = reflectivity(
        coating.n(wavelength),
        graze,
        polarization="unpolarized",
        roughness_nm=roughness_nm,
        wavelength_nm=wavelength,
        model=roughness_model,
    )
    return float(np.mean(np.where(visible, local, 0.0)))


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

    **Why a geometric mean.** The obvious local coefficient,
    :math:`r(\zeta_i(t))`, depends on :math:`\alpha` alone, which destroys the
    :math:`\alpha\leftrightarrow\beta_m` symmetry and breaks Lorentz
    reciprocity -- by 82% on the reference geometry, measured. Weighting a facet
    by both the angle it receives at and the angle it emits at,

    .. math:: r^{\text{gm}}(t) = \sqrt{r(\zeta_i(t))}\,\sqrt{r(\zeta_{d,m}(t))}

    is symmetric by construction, holds reciprocity to :math:`6\times10^{-15}`,
    and collapses to plain :math:`r(\zeta)` wherever the facet is at specular --
    so a perfectly blazed groove is unaffected by the change.

    **The square roots are taken separately**, rather than as
    :math:`\sqrt{r_i r_d}`, so that the only discontinuities left in the weight
    are the ones :math:`r` itself has.

    How much that is worth was measured rather than assumed, and it is less than
    it first appears. Where the product's argument wraps *uniformly* across the
    groove -- which is the grazing case, since :math:`\arg r` is then nearly the
    same at every :math:`t` -- the two forms differ by a **global** sign, and a
    global sign cancels out of :math:`|\int\cdots|^2`. On Au at 1-6 nm they
    agree to :math:`10^{-15}`.

    They part company only where the wrap is *non-uniform*: near normal
    incidence on a deep groove, where :math:`\arg r_p` sweeps through the
    Brewster jump and different parts of the same groove land on opposite sides
    of the cut. There the difference reaches 12% -- on orders carrying
    :math:`10^{-9}`, and in a regime `_brewster_crossing` already reports as
    outside what this symmetrisation can carry.

    So this is the cheap defensive choice, not a load-bearing one. Stated that
    way because the first draft of this docstring claimed the product wrapped in
    the grazing regime and that the separate roots were therefore essential;
    the measurement says otherwise.
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
    visible = reflection.visible_in[None, :] & (sin_out > 0.0)
    graze_out = np.arcsin(np.clip(sin_out, -1.0, 1.0))
    r_s_out, r_p_out = amplitude(n_c, graze_out, **common)

    zero = np.zeros((), dtype=np.complex128)
    w_s = np.where(visible, np.sqrt(r_s_in)[None, :] * np.sqrt(r_s_out), zero)
    w_p = np.where(visible, np.sqrt(r_p_in)[None, :] * np.sqrt(r_p_out), zero)

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


def _brewster_crossing(r_p: NDArray[np.complex128]) -> bool:
    r"""Does :math:`r_p` pass through zero somewhere on the lit groove?

    At Brewster's angle the p amplitude changes sign, so its phase jumps by
    :math:`\pi` and the half-phase the geometric mean takes jumps by
    :math:`\pi/2`. The jump is real physics, but the geometric-mean
    symmetrisation is not built to carry it, so the result is reported as
    suspect rather than quietly returned.

    Only reachable with steep grooves near normal incidence -- at the grazing
    angles this project is aimed at, :math:`r_p` stays well away from zero.
    """
    if r_p.size == 0:
        return False
    magnitude = np.abs(r_p)
    return bool(magnitude.min() < 0.05 * magnitude.max())


def _reflecting_graze(problem: Problem, illumination: Illumination):
    r"""``(zeta, description)`` -- the angle the surface reflects at.

    A blazed profile has one flat active facet tilted by ``blaze_angle``, and
    :math:`\sin\zeta = \sin\gamma\,\cos(\delta - \alpha)` is exact for it.

    Nothing else has a single facet angle, so the reflection is evaluated on
    the **mean surface** -- which is `facet_graze` with a zero tilt, not a
    separate formula:

    >>> facet_graze(gamma, 0.0, alpha)          # doctest: +SKIP
    ...                                          # == arcsin(|k_i . n|)

    How good that is depends on the profile, and the description says which:

    - **Lamellar** tops and bottoms genuinely *are* parallel to the mean
      surface, so this is exact for them. What it omits is the vertical walls,
      which at grazing incidence are nearly edge-on.
    - **Sinusoidal** and measured profiles have a local slope that varies
      across the groove, so the mean is an approximation and is recorded as
      one.

    The alternative -- refusing any profile without a facet angle -- would be
    worse than approximating: a sinusoid at grazing incidence does reflect,
    and a solver that declines to say how much has not become more honest.
    """
    blaze_angle = getattr(problem.profile, "blaze_angle", None)
    if blaze_angle is not None:
        return (
            facet_graze(illumination.gamma, np.radians(blaze_angle), illumination.alpha),
            f"active facet, {blaze_angle:g} deg tilt",
        )

    kind = type(problem.profile).__name__
    exact = kind == "Lamellar"
    return (
        facet_graze(illumination.gamma, 0.0, illumination.alpha),
        f"mean surface ({kind} has no single facet angle; "
        + (
            "exact for its flat tops and bottoms, omits the vertical walls)"
            if exact
            else "the local slope varies across the groove, so this is an "
            "approximation)"
        ),
    )


def _resolve_coating(problem: Problem) -> "OpticalConstants | None":
    """Turn ``Problem.coating`` into optical constants, or refuse.

    ``coating`` is a material *name*, so that a benchmark case stays a small
    serialisable file rather than embedding a 500-row table. Resolving it here
    rather than in ``Problem.__init__`` keeps a saved case loadable on a
    machine whose vendored tables have been trimmed -- it only has to resolve
    when someone actually solves.

    An unknown name **raises**. It is the whole point: the string used to be
    read by nothing, so ``coating="unobtanium"`` silently relabelled a result
    as absolute. A name that cannot be resolved is a question this solver
    cannot answer, not a default to fall back on.
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


#: The registered singleton.
scalar = register(ScalarSolver())
