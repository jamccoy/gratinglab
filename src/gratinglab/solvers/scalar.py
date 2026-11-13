r"""Scalar (Kirchhoff / thin-element) diffraction efficiency.

The efficiency is computed from the general Fourier integral, and the
closed-form sawtooth result is not used:

.. math::
    G_m = \int_0^1 e^{i\Phi_m(t)}\, e^{-2\pi i m t}\,dt, \qquad
    \Phi_m(t) = k\,g(t)\,\sin\gamma\,\left[\cos\alpha + \cos\beta_m\right]

where ``t`` is position normalised to the period and ``g(t)`` is the groove
height in nm. Efficiency is :math:`\mathscr{E}_m = |G_m|^2`.

This approach has two consequences.

First, it accepts any profile, including a measured AFM boundary, which the
closed forms cannot. This allows a scalar calculation and a PCGrate run to be
driven from identical geometry.

Second, the closed forms serve as tests instead of code paths. The Appendix D
results (sawtooth ``sinc²(φ/2 − mπ)``, square wave, sinusoid) are exact
answers that are already trusted, so they check the numerics instead of
duplicating them.

The phase depends on the order through :math:`\cos\beta_m` (ISSI eq. 15;
thesis Appendix-D.tex:651), so the method requires one integral per order
instead of a single FFT. This choice is deliberate and carries a known cost,
because no scalar formulation can satisfy all three of the following at
once:

==================  ==================================  =========================
Property            order-dependent (this)              beta-free alternative
==================  ==================================  =========================
Energy, sum <= 1    violated, up to 1.61 measured       exact (Parseval)
Reciprocity         exact, 1e-17                        violated, 0.44 measured
Blaze direction     exact at every alpha                exact at Littrow only
==================  ==================================  =========================

The conflict is structural. Reciprocity requires the phase to be symmetric
under :math:`\alpha \leftrightarrow \beta_m`, which brings both angles in,
whereas Parseval requires a single fixed function, which excludes
:math:`\beta_m`. Physical optics is reciprocal but not energy-conserving, and
the pure transmittance-function picture is the reverse.

This solver keeps reciprocity as its invariant. The energy deviation is
therefore expected; it is reported openly and is never renormalised away. It
is not an implementation defect: in the shallow-groove limit the sum is
1.00000 exactly, and it degrades smoothly as groove depth grows
(0.0001 -> 1.00000, 0.10 -> 0.92902). The deviation scales with the phase
excursion across the groove, that is, with depth relative to wavelength and
hence with the working diffraction order. It does not scale with
lambda/period or with the number of propagating orders.

``docs/theory/scalar.md`` §5 gives the full derivation and evidence, and
``docs/conventions.md`` §9 is the reference for every formula here.
"""

from __future__ import annotations

import time
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..geometry import (
    beta,
    cos_beta,
    flux_obliquity,
    interference_factor,
    is_propagating,
    sin_beta,
)
from ..illumination import Illumination
from ..problem import Problem
from ..result import EfficiencyScan, Provenance
from ..materials.fresnel import RoughnessModel, reflectivity
from .base import Capabilities, Progress, register
from ._scalar_groove import (
    ReflectivityModel,
    Visibility,
    _build_reflection,
    _groove_average_reflectivity,
    _local_reflected_efficiency,
    _Reflection,
    _resolve_coating,
    _shadow_masked_efficiency,
    _spanning_orders,
    _sp_split,
)

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..materials import OpticalConstants

__all__ = [
    "ReflectivityModel",
    "ScalarSolver",
    "Visibility",
    "interference_factor",
    "scalar",
]


#: Above this s/p reflectivity split, the unowned groove-TE/TM to facet-s/p
#: polarization mapping stops being free. Below it, any basis rotation moves
#: the unpolarized mean by at most half the split -- under 1% -- which is why
#: the mapping is documented rather than implemented; see
#: ``docs/theory/scalar.md`` section 9 and ``docs/findings.md``.
_SP_SPLIT_WARN = 0.02

#: Above this ratio of the **reduced** wavelength lambda/sin(gamma) to the
#: period, scalar theory is losing validity. In-plane this is plain
#: lambda/period; in a conical mount the reduction is what matters, because
#: the transverse problem the grooves pose lives at the reduced
#: wavelength (M&P eq. 4.65 -- the same decoupling the integral solver is
#: built on). Judged on lambda/period alone, the flagship off-plane geometry
#: reads 0.007 while its effective ratio is 0.32.
_REDUCED_RATIO_WARN = 0.1


# interference_factor lived here until the resolving-power milestone moved it
# to gratinglab.geometry, where its consumers outside the solver stack can
# reach it without importing a solver. Re-exported above (and kept in __all__)
# so `gratinglab.solvers.scalar.interference_factor` keeps resolving.


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
        visibility: "Visibility" = "facet-normal",
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        r"""Efficiency over a wavelength scan.

        Efficiency is :math:`\mathscr{E}_m = O_m |G_m|^2`, the norm-squared
        Fourier coefficient multiplied by the symmetric flux obliquity
        :math:`4\cos\alpha\cos\beta_m/(\cos\alpha+\cos\beta_m)^2`, with no
        renormalisation. Thesis Appendix D carries the asymmetric obliquity
        :math:`\cos\beta_m/\cos\alpha`, which breaks reciprocity, and a
        renormalisation to :math:`\sum_m \mathscr{E}_m = 1`, which destroys the
        model's own error signal; both have been removed. See
        :func:`~gratinglab.geometry.flux_obliquity` for the reason a flux
        factor is nonetheless included, and ``tests/test_perturbation.py`` for
        the independent theory that determines which factor.

        Parameters
        ----------
        quadrature_points
            Samples per period for the Fourier integral. The integrand is
            periodic, so the rectangle rule converges spectrally and the
            default is generous. The value must exceed twice the highest order
            to satisfy Nyquist, and this is checked.
        roughness_model
            How ``Problem.roughness`` damps the reflectivity: ``"nevot-croce"``
            (default), ``"debye-waller"``, or ``"none"``. Both models are
            legitimate and they differ, so the choice is exposed to the
            caller. Debye-Waller is the more pessimistic near the critical
            angle because it takes no account of the transmitted wave. The
            setting is ignored without a coating, since there is then no
            reflectivity to damp.
        reflectivity_model
            How reflectivity is resolved across the groove cycle. Ignored
            without a coating.

            ``"local"`` (default) evaluates the complex Fresnel amplitude at
            every quadrature point, from the local facet tilt, and carries it
            inside the diffraction integral. A groove whose reflectivity
            varies across the cycle is an amplitude grating as well as a phase
            grating, so this makes reflectivity order-dependent. The note
            below corrects the earlier claim in the module docstring that no
            such mechanism exists.

            ``"average"`` takes the groove-cycle mean of the intensity,
            :math:`\langle R(\zeta(t))\rangle`, and applies it as one factor
            per wavelength. It accounts for the shadowing and the varying
            local angle but not for the interference between them, so it
            remains order-independent. Like ``"facet"``, it breaks
            reciprocity, because its :math:`\zeta(t)` is built from
            :math:`\alpha` alone. Reciprocity is restored by symmetrising in
            the exit direction, which only ``"local"`` does; resolving the
            groove alone does not restore it. This model is offered as a
            diagnostic that separates shadowing from interference, and it is
            not intended as a physical model.

            ``"facet"`` is the M15 treatment: one :math:`R` per wavelength at
            the active-facet angle, applied to every order alike. It is kept
            so that a prior result can be reproduced exactly and the size of
            the change can be measured.

            Both resolved models zero the reflection where a facet is
            back-facing, in either direction. This follows from the geometry,
            and it can drive a weak order to exactly zero. When it does, the
            order is named in a provenance warning so that it is not mistaken
            for a passing-off.
        visibility
            Which shadows the masks above can see.

            ``"facet-normal"`` (default) is the local orientation test alone:
            a point is shadowed iff its own facet is turned away from the
            direction in question. This test misses cast shadows, in which the
            groove apex blocks surface beyond the trough that otherwise faces
            the ray. On the reference sawtooth these amount to a fraction of a
            percent of the period on the incident side but 10-50% per order on
            the exit side, moving the blaze order by ~-30%
            (``docs/findings.md``).

            ``"horizon"`` adds cast shadows, via
            :func:`~gratinglab.geometry.horizon_visible` in the transverse
            plane, for the incident direction and each diffracted one.
            Occlusion along a straight ray is the same from either end, so
            this keeps the :math:`\alpha \leftrightarrow \beta_m` symmetry
            and reciprocity is preserved. The treatment is ray optics, and
            near passing-off it overestimates blocking because diffraction
            bends around the apex (the vanishing flux obliquity already
            suppresses those orders). It remains opt-in for now so that
            existing results reproduce bit-for-bit, and, like the
            geometric-mean weight, it awaits validation against a rigorous
            finite-conductivity method.

            With ``"horizon"`` and no coating, the masks are applied at unit
            amplitude inside the integral. This gives geometric shadowing on
            an otherwise perfect reflector, the configuration directly
            comparable to the integral solver on deep grooves. The default
            uncoated run remains the unmodified pure phase integral. The
            ``"average"`` model sees the horizon only on the incident side
            (its exit side is order-blind by construction), and ``"facet"``,
            which has no masks at all, refuses the combination instead of
            ignoring it.
        progress
            Called ``(0, n)`` before the first wavelength and ``(k, n)`` after
            each one, per the contract in
            :class:`~gratinglab.solvers.base.Solver`. Exceptions raised by it
            are not caught, which allows a caller to raise
            :class:`~gratinglab.solvers.base.SolveCancelled` and stop
            this loop. A scan cancelled in this way returns nothing. A
            partially filled efficiency array, with the remaining rows still
            zero, would be indistinguishable from a scan whose orders all
            passed off, and returning one would produce the silently wrong
            number that this project is designed to prevent.

        Notes
        -----
        Summed efficiency will not equal 1, and can exceed it. This is
        expected; the module docstring describes the three-way tension that
        makes it unavoidable and explains why reciprocity is the invariant
        kept. The sum is reported and never rescaled, because renormalising by
        it (as Appendix D does) would erase the model's own signal about how
        far it has strayed from its regime.
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

        # ------------ geometry shared by every wavelength ------------
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
            problem, illumination, coating, reflectivity_model, visibility, t, height
        )
        resolved = coating is not None and reflection.model != "facet"

        efficiency = np.zeros((len(wavelengths), len(all_orders)))
        propagating = np.zeros_like(efficiency, dtype=bool)
        # One scale factor per wavelength for the two order-independent models.
        scale = np.ones(len(wavelengths))

        # ------------ one wavelength at a time ------------
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
            elif coating is None and reflection.visibility == "horizon":
                # Geometric shadowing on a perfect reflector: the masks at
                # unit amplitude, inside the integral. This is the uncoated
                # run made comparable to the integral solver on grooves deep
                # enough to shadow themselves.
                values = obliquity * _shadow_masked_efficiency(
                    integrand,
                    reflection,
                    beta(sines[live]),
                    all_orders[live],
                    illumination,
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

        # ------------ order-independent reflectivity, then the result ------------
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
                # Unpolarized, whatever the illumination says.
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

        # The *reduced* ratio, not lambda/period. A conical problem decouples
        # into an in-plane one at lambda/sin(gamma), so "wavelength small
        # against the structure" has to be judged after the reduction --
        # otherwise the guard passes exactly the extreme off-plane mounts it
        # exists for. tests/test_cross_method.py holds lambda/period fixed and
        # closes the cone to show the scalar-vs-integral discrepancy growing;
        # that measurement is what moved sin(gamma) into this line.
        ratio = float(
            wavelengths.max() / (problem.period * illumination.sin_gamma)
        )
        if ratio > _REDUCED_RATIO_WARN:
            warnings.append(
                f"reduced ratio lambda/(period sin gamma) reaches {ratio:.3g}; "
                "scalar theory needs the wavelength small against the "
                "structure, and in a conical mount that comparison happens at "
                "the reduced wavelength lambda/sin(gamma), so orderwise "
                f"accuracy degrades above {_REDUCED_RATIO_WARN} even where "
                "lambda/period looks safe. Order positions and the blaze "
                "envelope remain useful intuition here; per-order numbers "
                "should be cross-checked against a rigorous method -- see "
                "docs/theory/scalar.md section 7"
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
            #
            # Only where n is near 1, too. That is the regime in which a
            # critical angle exists and reflectivity collapses above it; from
            # the UV to the IR there is usually no critical angle at all, a
            # glass or a metal reflects without grazing, and a near-normal
            # design is ordinary rather than outside anything. Warning there
            # fired on every visible design, so it is the same M8 distinction
            # a second time: a check that does not apply.
            critical = coating.critical_angle(wavelengths)
            regime = coating.near_unity(wavelengths)
            past = regime & (graze > critical)
            if past.any():
                lam = wavelengths[past]
                no_edge = past & (critical == 0.0)
                where = (
                    f"theta_c falls to {np.degrees(critical[past].min()):.4g} "
                    "deg there"
                    if not no_edge.any()
                    else "theta_c is 0 at "
                    f"{no_edge.sum()} of them, where the real part of n is not "
                    "below 1 and there is no total external reflection at all"
                )
                warnings.append(
                    f"facet graze {np.degrees(graze):.4g} deg exceeds the "
                    f"critical angle for {coating.name} over "
                    f"{lam.min():.4g}-{lam.max():.4g} nm ({where}); n is near 1 "
                    "there, so reflectivity collapses above the critical angle "
                    "and these wavelengths are outside the regime a "
                    "grazing-incidence design operates in -- see "
                    "docs/theory/scalar.md section 7"
                )

            # The waived polarization mapping, priced. Every reflectivity
            # model here averages facet-local s and p 50/50 because the
            # groove-referenced TE/TM -> facet s/p rotation is
            # unowned (see materials/fresnel.py). That waiver costs at most
            # half the s/p split, which at the grazes this project runs is
            # under 0.5% -- but it is a measured bound, not a law, so the one
            # angle the guards already track is checked and the waiver is
            # reported exactly when it stops being free.
            split = _sp_split(coating, graze, wavelengths)
            if split > _SP_SPLIT_WARN:
                warnings.append(
                    f"s and p reflectivities differ by {100 * split:.3g}% at "
                    f"the guarded graze ({np.degrees(graze):.4g} deg); the "
                    "50/50 s/p average used here stands in for a "
                    "groove-TE/TM to facet-s/p mapping nobody carries, which "
                    "is free only while the two are degenerate. Read "
                    "polarization-sensitive conclusions from a rigorous "
                    "method -- see docs/theory/scalar.md section 9"
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
        # outside conditions it can answer for.

        if reflection.suppressed:
            # An order at exactly zero must never be mistaken for one that
            # passed off. Naming it is the whole difference between a
            # modelled result and a silent one. Outside the coating block
            # because the uncoated horizon mode suppresses orders too.
            named = ", ".join(str(m) for m in sorted(reflection.suppressed))
            mechanism = (
                "the local facet normal"
                if reflection.visibility == "facet-normal"
                else "the local facet normal and the cast-shadow horizon"
            )
            warnings.append(
                f"order(s) {named} are zero at one or more wavelengths in "
                f"this scan because no part of the groove faces both the "
                f"incident and the diffracted direction; "
                f"{100 * reflection.shadowed_fraction:.3g}% of the period "
                f"is shadowed from the beam to begin with. This is geometric "
                f"shadowing by {mechanism}, not passing-off -- "
                "`propagating` still marks these orders live. Scalar theory "
                "has no diffraction into a shadowed direction, so zero is "
                "the model's answer rather than the physical one"
            )

        if illumination.polarization != "unpolarized":
            warnings.append(
                f"result is labelled {illumination.polarization} but scalar "
                "theory neglects polarization; TE and TM are identical here"
            )

        if coating is not None:
            # A model (materials.typical) names the bands it cannot answer for.
            warnings.extend(coating.caveats(wavelengths))

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
                # Keyed on whether a reflectivity was applied, never
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
                    }
                    if coating is not None
                    else {}
                ),
                **(
                    # Whenever visibility masks participated --
                    # groove-resolved reflectivity, or the uncoated horizon
                    # mode. A "facet" run has no masks and gets neither key.
                    {
                        "visibility": reflection.visibility,
                        "shadowed_fraction": reflection.shadowed_fraction,
                    }
                    if reflection.masked
                    else {}
                ),
                "alpha_deg": illumination.alpha_deg,
                "gamma_deg": illumination.gamma_deg,
            },
        )


#: The registered singleton.
scalar = register(ScalarSolver())
