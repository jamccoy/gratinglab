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

Note that :math:`\Phi_m` depends on the order through :math:`\cos\beta_m`
(ISSI eq. 15; thesis Appendix-D.tex:651). This is one integral per order rather
than a single FFT, and it means the ``G_m`` are *not* Fourier coefficients of one
fixed function -- so Parseval does not give an exact sum rule. See
:func:`solve` for what that implies for energy balance.

``docs/conventions.md`` §9 is the reference for every formula here.
"""

from __future__ import annotations

import time

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..geometry import cos_beta, is_propagating, order_range, sin_beta
from ..illumination import Illumination
from ..problem import Problem
from ..result import EfficiencyScan, Provenance
from .base import Capabilities, register

__all__ = ["ScalarSolver", "interference_factor", "scalar"]

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
    )

    def solve(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: ArrayLike,
        *,
        quadrature_points: int = 2048,
        obliquity: bool = False,
    ) -> EfficiencyScan:
        r"""Efficiency over a wavelength scan.

        Parameters
        ----------
        quadrature_points
            Samples per period for the Fourier integral. The integrand is
            periodic, so the rectangle rule converges spectrally and the
            default is generous. Must exceed twice the highest order to satisfy
            Nyquist; this is checked.
        obliquity
            Include the :math:`\cos\beta_m/\cos\alpha` factor from the thesis
            (Appendix-D.tex:418). Default off, matching ISSI eq. (15) --
            see ``docs/conventions.md`` §5. The two redistribute efficiency
            between orders differently, so this is exposed to make the
            difference measurable rather than buried.

        Notes
        -----
        Summed efficiency over *propagating* orders is close to but not exactly
        1. Two reasons, both physical rather than numerical: power also goes
        into evanescent orders, and :math:`\Phi_m` is order-dependent so the
        coefficients are not a single Parseval pair. The deficit grows as
        scalar theory loses validity, which makes it a useful diagnostic in its
        own right -- it is reported, never rescaled away.
        """
        self.capabilities.check(problem, illumination)

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

        efficiency = np.zeros((len(wavelengths), len(all_orders)))
        propagating = np.zeros_like(efficiency, dtype=bool)

        for row, wavelength in enumerate(wavelengths):
            sines = sin_beta(
                all_orders, wavelength, problem.period, sin_alpha, sin_gamma
            )
            live = is_propagating(sines)
            if not live.any():
                continue

            cosines = cos_beta(sines[live])
            # Phi_m(t) = k g(t) sin(gamma) [cos(alpha) + cos(beta_m)]
            phase = (2.0 * np.pi / wavelength) * height * sin_gamma
            phi = phase[None, :] * (cos_alpha + cosines)[:, None]

            coefficients = np.mean(np.exp(1j * phi) * kernel[live], axis=1)
            values = np.abs(coefficients) ** 2

            if obliquity:
                values = values * cosines / cos_alpha

            efficiency[row, live] = values
            propagating[row, live] = True

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
                obliquity,
                time.perf_counter() - started,
            ),
        )

    def _provenance(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: NDArray[np.float64],
        quadrature_points: int,
        obliquity: bool,
        elapsed: float,
    ) -> Provenance:
        """Record the run, including every validity guard it tripped."""
        from .. import __version__

        warnings: list[str] = []

        ratio = float(wavelengths.max() / problem.period)
        if ratio > _LAMBDA_OVER_PERIOD_WARN:
            warnings.append(
                f"lambda/period reaches {ratio:.3g}; scalar theory assumes "
                f"lambda << period and loses accuracy above "
                f"{_LAMBDA_OVER_PERIOD_WARN}"
            )

        if problem.roughness > 0 and hasattr(problem.profile, "blaze_angle"):
            # Fraunhofer smoothness criterion, ISSI section 4.
            from ..geometry import facet_graze

            zeta = facet_graze(
                illumination.gamma,
                np.radians(problem.profile.blaze_angle),
                illumination.alpha,
            )
            threshold = 32.0 * np.sin(zeta) * problem.roughness
            if wavelengths.min() < threshold:
                warnings.append(
                    f"Fraunhofer smoothness criterion violated below "
                    f"{threshold:.4g} nm (32 sin(zeta) sigma); the surface is "
                    "not optically smooth there"
                )

        if problem.coating is None:
            warnings.append(
                "no coating specified, so the total-external-reflection guard "
                "(zeta < sqrt(2 delta_n)) could not be evaluated; efficiencies "
                "are relative, not absolute"
            )

        if illumination.polarization != "unpolarized":
            warnings.append(
                f"result is labelled {illumination.polarization} but scalar "
                "theory neglects polarization; TE and TM are identical here"
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
                "obliquity": obliquity,
                "normalization": "relative" if problem.coating is None else "absolute",
                "alpha_deg": illumination.alpha_deg,
                "gamma_deg": illumination.gamma_deg,
            },
        )


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
