r"""Fresnel reflectivity at a grazing interface, and what roughness does to it.

This is the factor that turns the scalar solver's relative groove efficiency
into an absolute one: light has to reflect off the facet before the groove can
diffract it.

Angles are measured from the surface
====================================

:math:`\theta = 0` is grazing and :math:`\theta = \pi/2` is normal incidence.
This is the same sense as :func:`gratinglab.geometry.facet_graze`, whose
:math:`\zeta` is what is passed in here, so the two compose without a
conversion. Measuring from the normal instead would introduce a silent
90-degree error that still returns numbers in :math:`[0, 1]`.

Polarization is ``s`` and ``p``, distinct from ``TE`` and ``TM``
================================================================

These are different conventions, and this module does not accept the
project's. ``conventions.md`` §7 defines TE and TM with respect to the
grooves, and notes that the opposite choice is common in the RCWA literature
and a frequent source of cross-code disagreement in conical mounts.

Fresnel reflection off a facet is a local flat-surface problem, so its ``s``
and ``p`` are defined with respect to the local plane of incidence on that
facet. In a conical mount those planes are not the groove plane, so the
mapping from TE/TM to s/p is not the identity. Accepting ``"TE"`` here would
implicitly assert that it is.

Whoever connects this module to a solver is responsible for that mapping and
has to state it. In the soft X-ray grazing regime the question has little
practical consequence: s and p agree to well under a percent below a few
degrees, which `tests/test_fresnel.py` measures directly. A companion test
shows that they separate appreciably at 45°, so the agreement is shown to be a
property of the regime and not of the formula.

Roughness is offered as two models
==================================

Névot–Croce uses the transmitted normal wavevector and is the appropriate model
near and below the critical angle; Debye–Waller uses only the incident one and
is the common approximation. Both are offered, with a stated default, because
they differ. Providing only one and calling it "the roughness factor" would
hide a modelling choice inside a number.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "RoughnessModel",
    "SurfacePolarization",
    "amplitude",
    "critical_angle",
    "debye_waller",
    "nevot_croce",
    "reflectivity",
]

#: Local to the facet, **not** the project's groove-referenced TE/TM. See the
#: module docstring.
SurfacePolarization = Literal["s", "p", "unpolarized"]

RoughnessModel = Literal["nevot-croce", "debye-waller", "none"]


def critical_angle(decrement: ArrayLike) -> NDArray[np.float64]:
    r"""Total-external-reflection graze angle, :math:`\theta_c=\arccos(1-\delta)`.

    Radians, from the surface. Below it a lossless material reflects
    everything, which `tests/test_fresnel.py` checks against :func:`amplitude`
    on both sides of the edge. The familiar :math:`\sqrt{2\delta}` is the
    leading term of this, short by a relative :math:`\delta/12`: about 1e-4 in
    the soft X-ray, and enough at a large decrement to put a graze the
    amplitudes still reflect totally on the wrong side of it.

    :class:`~.optical.OpticalConstants` has a method of the same name that
    reads a table; this is the bare relation, for a decrement you already have.
    It computes only the lossless edge, and does not ask whether the material
    is near enough to n = 1 for an edge to be meaningful, which the method
    does.
    """
    decrement = np.asarray(decrement, dtype=np.float64)
    if (decrement < 0).any():
        raise ValueError(
            "decrement must be non-negative; a negative one puts the real part "
            "of n above 1, which deletes total external reflection entirely"
        )
    if (decrement > 1).any():
        raise ValueError(
            "decrement must not exceed 1; beyond it the real part of n is "
            "negative and there is no lossless edge to compute"
        )
    return np.arccos(1.0 - decrement)


def _normal_components(
    n: ArrayLike, graze: ArrayLike
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    r"""Normal wavevector components either side of the interface, over ``k``.

    :math:`\tilde{k}_{z,1} = \sin\theta` in vacuum and
    :math:`\tilde{k}_{z,2} = \sqrt{n^2 - \cos^2\theta}` in the material, both
    divided by :math:`k = 2\pi/\lambda` so they are dimensionless and the
    caller multiplies back in when it needs a length.

    The branch matters: the transmitted wave must decay into the material,
    which under the ``exp(-i\omega t)`` convention of ``conventions.md`` §2
    means a non-negative imaginary part. NumPy's principal square root delivers
    that here, and `tests/test_fresnel.py` tests it explicitly, since the wrong
    branch gives a wave growing with depth and a reflectivity above one.
    """
    n = np.asarray(n, dtype=np.complex128)
    graze = np.asarray(graze, dtype=np.float64)
    return np.sin(graze).astype(np.complex128), np.sqrt(n**2 - np.cos(graze) ** 2)


def amplitude(
    n: ArrayLike,
    graze: ArrayLike,
    *,
    roughness_nm: float = 0.0,
    wavelength_nm: ArrayLike | None = None,
    model: RoughnessModel = "nevot-croce",
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    r"""Complex amplitude reflection coefficients :math:`(r_s, r_p)`.

    :func:`reflectivity` is this, norm-squared. The amplitudes are exposed
    separately because a caller that puts reflection inside a diffraction
    integral needs the phase: two facets at different local angles reflect with
    different phases as well as different magnitudes, and squaring first
    discards that information. The scalar solver's groove-cycle-resolved
    reflectivity model is the caller in question.

    Roughness enters at amplitude level, as the square root of the intensity
    factor returned by :func:`nevot_croce` or :func:`debye_waller`. That is the
    level at which Névot–Croce was originally formulated, before it was squared
    for :func:`reflectivity`. Applying the intensity factor to an amplitude
    would damp twice.

    Parameters are as for :func:`reflectivity`, without ``polarization``. This
    function returns both coefficients and leaves the caller to combine them,
    because a caller working at amplitude level must decide where to combine
    them, and doing it here would make that decision on its behalf.
    """
    k1, k2 = _normal_components(n, graze)
    n_c = np.asarray(n, dtype=np.complex128)

    # Amplitude coefficients. The `s` form is the one with no n^2 weighting;
    # `p` weights the incident-side term by n^2 because the boundary condition
    # is on H rather than E.
    r_s = (k1 - k2) / (k1 + k2)
    r_p = (n_c**2 * k1 - k2) / (n_c**2 * k1 + k2)

    if roughness_nm == 0.0 or model == "none":
        return r_s, r_p

    damping = np.sqrt(_roughness_factor(k1, k2, roughness_nm, wavelength_nm, model))
    return r_s * damping, r_p * damping


def reflectivity(
    n: ArrayLike,
    graze: ArrayLike,
    *,
    polarization: SurfacePolarization = "unpolarized",
    roughness_nm: float = 0.0,
    wavelength_nm: ArrayLike | None = None,
    model: RoughnessModel = "nevot-croce",
) -> NDArray[np.float64]:
    r"""Specular reflectivity of a flat interface from vacuum into ``n``.

    Parameters
    ----------
    n
        Complex index, :math:`1 - \text{decrement} + i\,\text{absorption}`.
        From :meth:`~.optical.OpticalConstants.n`.
    graze
        Angle from the surface, radians. :math:`\zeta` from
        :func:`~gratinglab.geometry.facet_graze`.
    polarization
        ``"s"``, ``"p"``, or ``"unpolarized"``, local to the facet. These are
        distinct from the project's TE/TM; see the module docstring.
    roughness_nm, wavelength_nm, model
        RMS roughness and the model applied to it. ``wavelength_nm`` is
        required when ``roughness_nm > 0``, because both factors compare
        :math:`\sigma` against a wavelength and neither is meaningful without
        one.

    Returns
    -------
    Real reflectivity in :math:`[0, 1]`.
    """
    if polarization not in ("s", "p", "unpolarized"):
        raise ValueError(
            f"polarization must be 's', 'p' or 'unpolarized', got "
            f"{polarization!r}. Note this is the facet-local convention, not "
            "the groove-referenced TE/TM of conventions.md section 7."
        )

    r_s, r_p = amplitude(
        n,
        graze,
        roughness_nm=roughness_nm,
        wavelength_nm=wavelength_nm,
        model=model,
    )

    if polarization == "s":
        result = np.abs(r_s) ** 2
    elif polarization == "p":
        result = np.abs(r_p) ** 2
    else:
        # conventions.md 7: computed, never assumed.
        result = 0.5 * (np.abs(r_s) ** 2 + np.abs(r_p) ** 2)
    return np.asarray(result, dtype=np.float64)


def _roughness_factor(
    k1: NDArray[np.complex128],
    k2: NDArray[np.complex128],
    roughness_nm: float,
    wavelength_nm: ArrayLike | None,
    model: RoughnessModel,
) -> NDArray[np.float64]:
    if roughness_nm < 0:
        raise ValueError(f"roughness_nm must be non-negative, got {roughness_nm}")
    if wavelength_nm is None:
        raise ValueError(
            "wavelength_nm is required when roughness_nm > 0: both roughness "
            "models compare sigma against a wavelength, and neither means "
            "anything without one"
        )

    k = 2.0 * np.pi / np.asarray(wavelength_nm, dtype=np.float64)
    if model == "nevot-croce":
        return nevot_croce(k * k1, k * k2, roughness_nm)
    if model == "debye-waller":
        return debye_waller(k * k1, roughness_nm)
    raise ValueError(
        f"model must be 'nevot-croce', 'debye-waller' or 'none', got {model!r}"
    )


def nevot_croce(
    k_iz: ArrayLike, k_tz: ArrayLike, roughness_nm: float
) -> NDArray[np.float64]:
    r"""Névot–Croce factor, :math:`\left|e^{-2 k_{iz} k_{tz} \sigma^2}\right|^2`.

    Uses the transmitted normal component as well as the incident one, which
    is what makes it behave correctly through the critical angle: below
    :math:`\theta_c` the transmitted wave is evanescent, :math:`k_{tz}` is
    nearly imaginary, and the product barely damps the reflection. Debye–Waller
    has no such term and over-damps there.

    This returns an intensity factor, which is the reason for the square.
    Névot–Croce is usually quoted on the amplitude,
    :math:`r = r_F e^{-2k_{iz}k_{tz}\sigma^2}`, so the reflectivity carries its
    modulus squared. Both functions in this
    module return intensity factors, because :func:`reflectivity` multiplies
    them into :math:`|r|^2`; :func:`amplitude` takes the square root back out.

    The test that fixes this is the :math:`n \to 1` limit: the transmitted wave
    becomes the incident one, :math:`k_{tz} \to k_{iz}`, and Névot–Croce must
    reduce to Debye–Waller. Without the square it does not, and the discrepancy
    is large, 0.29 against 0.085 at 30° graze. This function returned the
    unsquared value until M16-B, which made every roughened reflectivity in the
    project too high (0.73 where 0.54 was right, at 15° graze with σ = 0.5 nm at
    2 nm).

    The modulus is taken because :math:`k_{tz}` is complex, so the exponential
    is too, while a reflectivity is real. Wavevectors in nm⁻¹, matching
    ``roughness_nm``.
    """
    k_iz = np.asarray(k_iz, dtype=np.complex128)
    k_tz = np.asarray(k_tz, dtype=np.complex128)
    damping = np.abs(np.exp(-2.0 * k_iz * k_tz * roughness_nm**2))
    return (damping**2).astype(np.float64)


def debye_waller(k_iz: ArrayLike, roughness_nm: float) -> NDArray[np.float64]:
    r"""Debye–Waller factor, :math:`e^{-(2 k_{iz}\sigma)^2}`.

    This is the familiar approximation, and the more pessimistic of the two
    near the critical angle. It depends only on the incident wave and not on
    the material, so it damps by the same amount whether or not the wave
    penetrates. It is offered because much of the literature quotes it.
    """
    k_iz = np.asarray(k_iz, dtype=np.complex128)
    return np.exp(-((2.0 * np.abs(k_iz) * roughness_nm) ** 2)).astype(np.float64)
