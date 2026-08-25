r"""The finite-conductivity transverse solve (Goray & Schmidt 2010).

One interface between vacuum and a semi-infinite material of complex
refractive index ``n``, in the conical mount. The two coupled scalars are
``E_z`` and ``B_z = (mu+/eps+)^{1/2} H_z`` (G&S section 2.A); each satisfies
a transverse Helmholtz equation with its own side's wavenumber

.. math::
    k_t^{+} = k \sin\gamma, \qquad
    k_t^{-} = k \sqrt{n^2 - \cos^2\gamma},

so the vacuum side is exactly the reduced-wavelength problem the perfectly
conducting solve already meshes, and the metal side is the same Green's
function code with complex ``k`` (the branch pin in ``_greens``).

The system is derived here independently in this project's conventions
(outward normal up, unscaled operators of ``_nystrom``, ``exp(-i omega t)``)
rather than transcribed -- G&S's normal points into the metal and their
potentials carry a factor 2, and sign-sensitive equations do not survive PDF
extraction (``docs/references.md``). The construction:

- The substrate fields are single-layer potentials, ``u^- = S^- w`` and
  ``v^- = S^- tau``, so their traces are ``V^- w`` and their normal
  derivatives ``(L^- - I/2) w`` (the jump relation ``_nystrom`` tests pin).
- On the vacuum side, Green's representation of the radiating scattered
  field plus the regularity identity of the incident field give the exact
  boundary relation ``(I/2 + K^+) psi - V^+ (d psi/dn) = psi_i`` for the
  *total* field ``psi`` -- whose perfectly conducting limits are literally
  the two milestone-1 equations (``d psi/dn = 0`` gives the TM equation,
  ``psi = 0`` the TE one), which is what pins its signs.
- The jump conditions (G&S eq. 6, mapped) eliminate the vacuum traces:
  continuity of ``E_z, B_z``, and

  .. math::
     \partial_n E_z^+ = c_E\, \partial_n E_z^- + s\, \partial_t B_z, \qquad
     \partial_n B_z^+ = c_B\, \partial_n B_z^- - s\, \partial_t E_z,

  with ``c_E = n^2 k_t^{+2}/k_t^{-2}``, ``c_B = k_t^{+2}/k_t^{-2}``, and
  ``s = \cos\gamma\,(1 - k_t^{+2}/k_t^{-2})`` -- the tangential derivatives
  of the *continuous* traces, hence ``D_t V^-`` on the densities. That
  product is assembled as one operator (``tangential_layer_matrix``), with
  the derivative taken analytically on the kernel: differentiating an
  assembled ``V^-`` numerically instead leaves an ``O(1)`` error that
  refinement does not remove, which destroys the conical solve on any
  boundary carrying appreciable high-frequency content.

The result is a 2x2 block system in ``(w, tau)`` whose diagonal blocks are
``(I/2 + K^+)V^- - c\,V^+(L^- - I/2)`` and whose cross blocks are
``\mp s\,V^+ D_t V^-``. In-plane (``cos gamma = 0``) the cross blocks
vanish and the two equations decouple; both incident polarizations share
the same blocks, differing only in the right-hand side. In the perfectly
conducting limit ``|n| -> inf`` the system degenerates to the milestone-1
equations (``V^-, c_B -> 0``), which the tests measure at ``n = 100i``.

Efficiencies are **absolute**: per reflected order
``(gamma_m/gamma_0)(|E_m|^2 + |B_m|^2)`` for a unit-power incident state
``|p_z|^2 + |q_z|^2 = 1``. Absorption is computed from the solved densities
by an independent boundary integral (G&S eq. 26, mapped and re-derived --
it reduces exactly to ``1 - |r|^2`` on a flat interface, both
polarizations), so ``R + A = 1`` is a genuine two-sided check, the
successor of the perfectly conducting ``sum = 1`` theorem. For a lossless
substrate the same integral equals the transmitted power and ``R + T = 1``
is a theorem again.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...geometry import order_range
from ._boundary import PhysicalBoundary
from ._greens import wavenumbers
from ._nystrom import (
    OperatorAssembler,
    adjoint_layer_matrix,
    double_layer_matrix,
    single_layer_matrix,
    tangential_layer_matrix,
)

__all__ = ["FiniteSolution", "solve_transverse_finite", "solve_finite_states"]


@dataclass(frozen=True, slots=True)
class FiniteSolution:
    """One incident polarization state at one wavelength, absolute units.

    ``e_amplitudes`` and ``b_amplitudes`` are the Rayleigh coefficients of
    ``E_z`` and ``B_z`` for the reflected propagating orders; both carry
    power because conical incidence converts polarization. ``absorption``
    is the boundary-integral value: the power crossing into the substrate,
    which for a lossy material is absorbed there and for a lossless one is
    the transmitted power (then it agrees with ``transmitted_total``).
    """

    orders: NDArray[np.int64]
    e_amplitudes: NDArray[np.complex128]
    b_amplitudes: NDArray[np.complex128]
    efficiencies: NDArray[np.float64]
    absorption: float
    transmitted_orders: NDArray[np.int64]
    transmitted_efficiencies: NDArray[np.float64]

    @property
    def total(self) -> float:
        """Sum of reflected efficiencies."""
        return float(self.efficiencies.sum())

    @property
    def transmitted_total(self) -> float:
        return float(self.transmitted_efficiencies.sum())


def _reflected_amplitudes(
    boundary: PhysicalBoundary,
    single: NDArray[np.complex128],
    double: NDArray[np.complex128],
    *,
    kp: float,
    alpha0: float,
    period: float,
    orders: NDArray[np.int64],
) -> NDArray[np.complex128]:
    """Rayleigh coefficients of ``S^+ single - D^+ double`` (M&P 4.33/4.41).

    With ``single`` the total field's normal derivative and ``double`` its
    trace, these are the reflected amplitudes -- the perfectly conducting
    extractions in ``_core`` are the special cases ``double = 0`` (TE) and
    ``single = 0`` (TM).
    """
    x, y = boundary.x, boundary.y
    spacing = boundary.spacing
    alpha_m, gamma_m = wavenumbers(kp, alpha0, period, orders)
    outgoing = np.exp(-1j * alpha_m[:, None] * x - 1j * gamma_m[:, None] * y)
    direction = (
        (alpha_m[:, None] / gamma_m[:, None]) * boundary.nx[None, :]
        + boundary.ny[None, :]
    )
    return spacing * (
        (outgoing @ single) / (2j * period * gamma_m)
        + ((outgoing * direction) @ double) / (2.0 * period)
    )


def _transmitted_amplitudes(
    boundary: PhysicalBoundary,
    density: NDArray[np.complex128],
    *,
    km: complex,
    alpha0: float,
    period: float,
    orders: NDArray[np.int64],
) -> NDArray[np.complex128]:
    """Rayleigh coefficients of the single-layer field ``S^- density`` below
    the boundary (downgoing spectral form of the metal-side kernel)."""
    x, y = boundary.x, boundary.y
    spacing = boundary.spacing
    alpha_m, gamma_m = wavenumbers(km, alpha0, period, orders)
    downgoing = np.exp(-1j * alpha_m[:, None] * x + 1j * gamma_m[:, None] * y)
    return spacing * (downgoing @ density) / (2j * period * gamma_m)


def solve_transverse_finite(
    boundary: PhysicalBoundary,
    *,
    wavelength: float,
    period: float,
    sin_alpha: float,
    index: complex,
    cos_gamma: float = 0.0,
    incident: tuple[complex, complex] = (1.0, 0.0),
    terms: int,
) -> FiniteSolution:
    """Solve one incident state ``incident = (p_z, q_z)`` at one wavelength.

    ``wavelength`` is the *vacuum* wavelength in nm (both sides' transverse
    wavenumbers derive from it); ``cos_gamma = 0`` is the in-plane mount,
    where the system decouples and ``(1, 0)``/``(0, 1)`` are the classical
    TE/TM problems. Efficiencies are normalized to the incident power
    ``|p_z|^2 + |q_z|^2``.
    """
    return solve_finite_states(
        boundary,
        wavelength=wavelength,
        period=period,
        sin_alpha=sin_alpha,
        index=index,
        cos_gamma=cos_gamma,
        incidents=[incident],
        terms=terms,
    )[0]


def solve_finite_states(
    boundary: PhysicalBoundary,
    *,
    wavelength: float,
    period: float,
    sin_alpha: float,
    index: complex,
    cos_gamma: float = 0.0,
    incidents: "list[tuple[complex, complex]]",
    terms: int,
    assembler: OperatorAssembler | None = None,
) -> "list[FiniteSolution]":
    """Solve several incident states against one assembly.

    Every state shares the same system matrices -- the expensive part (four
    kernel builds, the operator products, one factorization) is paid once
    and the states differ only in the right-hand side, which is how the
    unpolarized average costs two triangular solves rather than two
    factorizations. ``assembler`` shares the kernel machinery across
    wavelengths (the solver passes one per scan); built fresh when omitted.
    """
    index = complex(index)
    if index.imag < 0.0:
        raise ValueError(
            f"index must have a non-negative imaginary part, got {index}: "
            "a gain medium breaks the uniqueness of the transmission problem "
            "(Goray & Schmidt 2010, section 2.B)"
        )
    k = 2.0 * np.pi / wavelength
    sin_gamma = float(np.sqrt(1.0 - cos_gamma**2))
    kp = k * sin_gamma
    km = complex(k * np.sqrt(index**2 - cos_gamma**2))
    if km.imag < 0.0:  # pragma: no cover - principal sqrt of Im >= 0 argument
        km = -km

    cos_alpha = float(np.sqrt(1.0 - sin_alpha**2))
    alpha0 = -kp * sin_alpha
    gamma0 = kp * cos_alpha

    ktp2 = kp**2
    ktm2 = km**2
    c_e = index**2 * ktp2 / ktm2
    c_b = ktp2 / ktm2
    s = cos_gamma * (1.0 - ktp2 / ktm2)

    if assembler is None:
        assembler = OperatorAssembler(boundary, period=period, terms=terms)

    points = len(boundary.x)
    half = 0.5 * np.eye(points)
    v_plus = single_layer_matrix(
        boundary, k=kp, alpha0=alpha0, period=period, terms=terms,
        assembler=assembler,
    )
    v_minus = single_layer_matrix(
        boundary, k=km, alpha0=alpha0, period=period, terms=terms,
        assembler=assembler,
    )
    k_plus = double_layer_matrix(
        boundary, k=kp, alpha0=alpha0, period=period, terms=terms,
        assembler=assembler,
    )
    l_minus = adjoint_layer_matrix(
        boundary, k=km, alpha0=alpha0, period=period, terms=terms,
        assembler=assembler,
    )

    trace_block = (half + k_plus) @ v_minus
    derivative_block = v_plus @ (l_minus - half)

    incident_trace = np.exp(1j * (alpha0 * boundary.x - gamma0 * boundary.y))
    states = [(complex(p), complex(q)) for p, q in incidents]

    conical = cos_gamma != 0.0
    if conical:
        dt_v_minus = tangential_layer_matrix(
            boundary, k=km, alpha0=alpha0, period=period, terms=terms,
            assembler=assembler,
        )
        cross = v_plus @ dt_v_minus
        system = np.block(
            [
                [trace_block - c_e * derivative_block, -s * cross],
                [s * cross, trace_block - c_b * derivative_block],
            ]
        )
        rhs = np.empty((2 * points, len(states)), dtype=np.complex128)
        for column, (p_z, q_z) in enumerate(states):
            rhs[:points, column] = p_z * incident_trace
            rhs[points:, column] = q_z * incident_trace
        stacked = np.linalg.solve(system, rhs)
        densities = [
            (stacked[:points, column], stacked[points:, column])
            for column in range(len(states))
        ]
    else:
        rhs_e = np.stack([p_z * incident_trace for p_z, _ in states], axis=-1)
        rhs_b = np.stack([q_z * incident_trace for _, q_z in states], axis=-1)
        w_all = np.linalg.solve(trace_block - c_e * derivative_block, rhs_e)
        tau_all = np.linalg.solve(trace_block - c_b * derivative_block, rhs_b)
        densities = [
            (w_all[:, column], tau_all[:, column])
            for column in range(len(states))
        ]

    orders = order_range(2.0 * np.pi / kp, period, sin_alpha)
    _, gamma_m = wavenumbers(kp, alpha0, period, orders)
    flux = np.abs(gamma_m.real) / gamma0

    # Transmitted orders exist only for a lossless substrate that admits
    # propagating waves; a complex index sends everything evanescent.
    transmitting = index.imag == 0.0 and ktm2.real > 0.0
    if transmitting:
        km_real = float(km.real)
        transmitted_orders = order_range(
            2.0 * np.pi / km_real, period, kp * sin_alpha / km_real
        )
        _, gamma_minus = wavenumbers(km, alpha0, period, transmitted_orders)
        t_flux = np.abs(gamma_minus.real) / gamma0
    else:
        transmitted_orders = np.array([], dtype=np.int64)

    spacing = boundary.spacing
    solutions = []
    for (p_z, q_z), (w, tau) in zip(states, densities):
        power = abs(p_z) ** 2 + abs(q_z) ** 2

        # Total-field traces and normal derivatives on the vacuum side, via
        # the jump conditions -- the densities of the amplitude extraction.
        e_trace = v_minus @ w
        b_trace = v_minus @ tau
        e_normal_minus = (l_minus - half) @ w
        b_normal_minus = (l_minus - half) @ tau
        if conical:
            dt_e = dt_v_minus @ w
            dt_b = dt_v_minus @ tau
            e_normal_plus = c_e * e_normal_minus + s * dt_b
            b_normal_plus = c_b * b_normal_minus - s * dt_e
        else:
            e_normal_plus = c_e * e_normal_minus
            b_normal_plus = c_b * b_normal_minus

        e_amps = _reflected_amplitudes(
            boundary, e_normal_plus, e_trace,
            kp=kp, alpha0=alpha0, period=period, orders=orders,
        )
        b_amps = _reflected_amplitudes(
            boundary, b_normal_plus, b_trace,
            kp=kp, alpha0=alpha0, period=period, orders=orders,
        )
        efficiencies = (
            flux * (np.abs(e_amps) ** 2 + np.abs(b_amps) ** 2) / power
        )

        if transmitting:
            te_amps = _transmitted_amplitudes(
                boundary, w, km=km, alpha0=alpha0, period=period,
                orders=transmitted_orders,
            )
            tb_amps = _transmitted_amplitudes(
                boundary, tau, km=km, alpha0=alpha0, period=period,
                orders=transmitted_orders,
            )
            transmitted = (
                t_flux
                * (
                    c_e.real * np.abs(te_amps) ** 2
                    + c_b.real * np.abs(tb_amps) ** 2
                )
                / power
            )
        else:
            transmitted = np.array([], dtype=np.float64)

        # Absorption from the solved densities (G&S eq. 26, mapped): an
        # independent post-processing of the same solve, so R + A = 1
        # measures the discretisation rather than restating the order sum.
        # On a flat interface it reduces exactly to 1 - |r|^2.
        integral_e = spacing * np.sum(e_normal_minus * np.conj(e_trace))
        integral_b = spacing * np.sum(b_normal_minus * np.conj(b_trace))
        bracket = -(index**2) * integral_e - integral_b
        if conical:
            integral_t = spacing * np.sum(e_trace * np.conj(dt_b))
            bracket = bracket + 2.0 * cos_gamma * integral_t.real
        absorption = float(
            (ktp2 / (gamma0 * period * power)) * (bracket / ktm2).imag
        )

        solutions.append(
            FiniteSolution(
                orders=orders,
                e_amplitudes=e_amps,
                b_amplitudes=b_amps,
                efficiencies=efficiencies,
                absorption=absorption,
                transmitted_orders=transmitted_orders,
                transmitted_efficiencies=transmitted,
            )
        )
    return solutions
