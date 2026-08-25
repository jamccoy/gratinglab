r"""The in-plane perfectly-conducting solve -- the permanent transverse core.

The conical (off-plane) problem for a perfect conductor decouples exactly
into this in-plane problem at the reduced wavelength
``lambda / sin(gamma)`` (M&P eq. 4.65; thesis Chapter-2.tex:931), so this
module never sees ``gamma``: the caller substitutes the reduced wavelength
and maps the orders straight through ``geometry``.

TE (E along the grooves) is a Dirichlet problem solved with a single-layer
ansatz: ``\int G sigma ds' = -psi_i`` on the boundary (M&P eq. 4.31), then

.. math::
    r_m = \frac{\Delta s}{2 i d \gamma_m}
          \sum_j \sigma_j e^{-i \alpha_m x_j - i \gamma_m y_j}

(eq. 4.33). TM (H along the grooves) is a Neumann problem for the total
boundary field: ``psi/2 + \int (dG/dn') psi ds' = psi_i`` (M&P eq. 4.39
carries the equivalent statement in its own sign convention), then

.. math::
    r_m = \frac{\Delta s}{2 d} \sum_j \psi_j
          \left[\frac{\alpha_m}{\gamma_m} n_{j,x} + n_{j,y}\right]
          e^{-i \alpha_m x_j - i \gamma_m y_j}

(eq. 4.41). Efficiencies are ``(gamma_m / gamma_0) |r_m|^2`` and satisfy
``sum = 1`` exactly in the continuum -- the energy-balance theorem this
solver is tested against (eqs. 4.34, 4.42; thesis eq:prop_order_unity).

Amplitudes are returned alongside efficiencies: the planned
finite-conductivity milestone (Goray & Schmidt 2010) assembles a coupled
system around this core and needs phases, not just moduli.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...geometry import order_range
from ._boundary import PhysicalBoundary
from ._greens import wavenumbers
from ._nystrom import OperatorAssembler, dirichlet_matrix, neumann_matrix

__all__ = ["TransverseSolution", "solve_transverse"]


@dataclass(frozen=True, slots=True)
class TransverseSolution:
    """One polarization at one (reduced) wavelength."""

    orders: NDArray[np.int64]
    amplitudes: NDArray[np.complex128]
    efficiencies: NDArray[np.float64]

    @property
    def total(self) -> float:
        return float(self.efficiencies.sum())


def solve_transverse(
    boundary: PhysicalBoundary,
    *,
    wavelength: float,
    period: float,
    sin_alpha: float,
    polarization: str,
    terms: int,
    assembler: OperatorAssembler | None = None,
) -> TransverseSolution:
    """Solve the in-plane perfectly-conducting problem for one polarization.

    ``wavelength`` is the *transverse* (reduced) wavelength in nm;
    ``polarization`` is ``"TE"`` or ``"TM"`` (the unpolarized average is two
    calls, taken by the solver above). ``terms`` is the spectral truncation
    of the kernels. ``assembler`` shares the kernel machinery across calls
    (the solver passes one per scan); built fresh when omitted.
    """
    k = 2.0 * np.pi / wavelength
    cos_alpha = float(np.sqrt(1.0 - sin_alpha**2))
    alpha0 = -k * sin_alpha
    gamma0 = k * cos_alpha

    orders = order_range(wavelength, period, sin_alpha)
    alpha_m, gamma_m = wavenumbers(k, alpha0, period, orders)

    x, y = boundary.x, boundary.y
    spacing = boundary.spacing
    incident = np.exp(1j * (alpha0 * x - gamma0 * y))
    outgoing = np.exp(-1j * alpha_m[:, None] * x - 1j * gamma_m[:, None] * y)

    if assembler is None:
        assembler = OperatorAssembler(boundary, period=period, terms=terms)

    if polarization == "TE":
        system = dirichlet_matrix(
            boundary, k=k, alpha0=alpha0, period=period, terms=terms,
            assembler=assembler,
        )
        density = np.linalg.solve(system, -incident)
        amplitudes = (
            spacing * (outgoing @ density) / (2j * period * gamma_m)
        )
    elif polarization == "TM":
        system = neumann_matrix(
            boundary, k=k, alpha0=alpha0, period=period, terms=terms,
            assembler=assembler,
        )
        field = np.linalg.solve(system, incident)
        direction = (
            (alpha_m[:, None] / gamma_m[:, None]) * boundary.nx[None, :]
            + boundary.ny[None, :]
        )
        amplitudes = spacing * ((outgoing * direction) @ field) / (2.0 * period)
    else:  # pragma: no cover - the solver above never passes anything else
        raise ValueError(f"polarization must be 'TE' or 'TM', got {polarization!r}")

    # order_range keeps only propagating orders, so gamma_m is real here and
    # the flux ratio is a plain quotient. abs() per conventions.md section 6.
    flux = np.abs(gamma_m.real) / gamma0
    return TransverseSolution(
        orders=orders,
        amplitudes=amplitudes,
        efficiencies=flux * np.abs(amplitudes) ** 2,
    )
