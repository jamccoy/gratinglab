r"""The C-method at one wavelength: modes on each side, then the boundary.

Notation follows ``docs/theory/cmethod.md``. Every length is in nm. Grooves
run along ``z`` with ``k_z = k cos(gamma)`` conserved, so in a homogeneous
medium of permittivity ``eps`` both ``E_z`` and ``H_z`` (``H`` scaled by the
vacuum impedance) satisfy the transverse Helmholtz equation with
``kappa^2 = k^2 eps - k_z^2``. In the translation coordinates
``(x, u = y - a(x))`` a solution ``exp(i rho u) sum_n F_n exp(i alpha_n x)``
obeys the quadratic eigenproblem

    rho^2 (I + A A) F - rho (alpha A + A alpha) F - (kappa^2 - alpha^2) F = 0

with ``A`` the Toeplitz matrix of ``a'``. ``A A`` is the matrix product, not
the Toeplitz matrix of ``a'^2``: the derivation keeps ``E``, ``dE/dy`` and
``dE/dx`` at fixed ``y`` as the unknowns, all three continuous across a facet
corner along ``u = const``, so every product is of a discontinuous ``a'``
with a continuous factor and Laurent's rule is the correct one throughout
(Li's factorization rules; ``docs/theory/cmethod.md`` section 3).

The quantity matched at the boundary is the flux
``Phi = sqrt(1 + a'^2) dF/dn = -a' dF/dx|_y + dF/dy``, which in the same
factorization is ``Phi = -i A alpha F + i rho (I + A A) F``.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as _field

import numpy as np
from numpy.typing import NDArray

from ._series import ProfileSeries

__all__ = ["Harmonics", "Modes", "solve_perfect", "solve_finite", "Solution"]

#: Relative size of ``Im rho`` below which an eigenvalue of the lossless
#: vacuum problem is taken to be real (a propagating order, replaced by its
#: exact plane wave). Measured: the real eigenvalues come back with
#: ``|Im rho| / kappa`` near 1e-13; the smallest genuinely evanescent one is
#: ``sqrt(alpha_m^2 - kappa^2)``, which is small only at a Rayleigh anomaly.
_REAL_TOLERANCE = 1e-8


class ModeCountError(ValueError):
    """The eigenvalues could not be split into up- and downgoing sets."""


@dataclass(frozen=True)
class Harmonics:
    """The retained Fourier harmonics ``center - N .. center + N``.

    ``slope`` is the Toeplitz matrix of ``a'`` over them, built once per
    scan, since it depends on the profile and not on the wavelength.
    """

    series: ProfileSeries
    center: int
    truncation: int
    slope: NDArray[np.complex128]

    @classmethod
    def build(cls, series: ProfileSeries, center: int, truncation: int) -> "Harmonics":
        n = 2 * truncation + 1
        c = series.slope(np.arange(-2 * truncation, 2 * truncation + 1))
        i = np.arange(n)
        slope = c[(i[:, None] - i[None, :]) + 2 * truncation]
        return cls(series, center, truncation, slope)

    @property
    def orders(self) -> NDArray[np.int64]:
        return np.arange(-self.truncation, self.truncation + 1) + self.center

    @property
    def size(self) -> int:
        return 2 * self.truncation + 1


@dataclass(frozen=True)
class Modes:
    """Columns ``F`` and ``Phi`` of the modes kept on one side, ``(n, q)``."""

    field: NDArray[np.complex128]
    flux: NDArray[np.complex128]
    rho: NDArray[np.complex128]


def _eigen(h: Harmonics, alpha, kappa2):
    a = h.slope
    n = h.size
    metric = np.eye(n) + a @ a
    cross = alpha[:, None] * a + a * alpha[None, :]
    lower = np.linalg.solve(metric, np.diag(kappa2 - alpha**2 + 0j))
    right = np.linalg.solve(metric, cross)
    matrix = np.block([[np.zeros((n, n)), np.eye(n)], [lower, right]])
    rho, vectors = np.linalg.eig(matrix)
    field = vectors[:n]
    field = field / np.abs(field).max(axis=0, keepdims=True)
    return rho, field, metric


def _flux(h: Harmonics, alpha, metric, rho, field):
    return -1j * (h.slope @ (alpha[:, None] * field)) + 1j * rho[None, :] * (
        metric @ field
    )


def lossy_modes(h: Harmonics, alpha, kappa2: complex) -> Modes:
    """The ``2N + 1`` modes decaying downward (``Im rho < 0``) in a lossy medium."""
    rho, field, metric = _eigen(h, alpha, kappa2)
    down = rho.imag < 0.0
    if int(down.sum()) != h.size:
        raise ModeCountError(
            f"{int(down.sum())} of {2 * h.size} eigenvalues decay downward in "
            f"the material, where exactly {h.size} must; the material's loss is "
            "too small to separate them at this truncation"
        )
    return Modes(field[:, down], _flux(h, alpha, metric, rho[down], field[:, down]), rho[down])


def lossless_evanescent(
    h: Harmonics, alpha, kappa2: float, n_propagating: int, sign: int
) -> Modes:
    """The evanescent modes of a lossless medium, decaying away from the boundary.

    ``kappa2`` is the real ``kappa^2``, negative for a lossless medium with
    negative permittivity, in which nothing propagates. ``sign`` is ``+1``
    above the boundary (``Im rho > 0``) and ``-1`` below.
    The real eigenvalues are the propagating orders, which the caller
    represents by their exact plane waves instead; they must number exactly
    twice the propagating orders, or the split is not trusted.
    """
    rho, field, metric = _eigen(h, alpha, kappa2)
    tolerance = _REAL_TOLERANCE * max(
        np.sqrt(abs(kappa2)), float(np.abs(alpha).max())
    )
    real = np.abs(rho.imag) <= tolerance
    if int(real.sum()) != 2 * n_propagating:
        raise ModeCountError(
            f"the eigenproblem has {int(real.sum())} real eigenvalues where the "
            f"grating equation gives {2 * n_propagating} ({n_propagating} "
            "propagating orders, each up and down). The truncated problem has "
            "a spurious real mode, which happens on deep or steep profiles; "
            "raise the truncation, and if it persists the profile is beyond "
            "what this solver can do reliably"
        )
    away = sign * rho.imag > tolerance
    return Modes(
        field[:, away], _flux(h, alpha, metric, rho[away], field[:, away]), rho[away]
    )


def plane_waves(h: Harmonics, alpha, kappa: float, columns, betas) -> Modes:
    """Exact plane waves ``exp(i (alpha_m x + beta_m y))``, one per column.

    ``columns`` indexes the retained harmonics, ``betas`` the signed
    ``k_y``. With ``L_j`` the Fourier coefficients of ``exp(i beta a)``, the
    wave's field is ``L_{n-m}`` and its flux
    ``i L_{n-m} (kappa^2 - alpha_m alpha_n) / beta_m``, from
    ``a' exp(i beta a) = (i beta)^-1 d/dx exp(i beta a)``.
    """
    n = h.size
    index = np.arange(n)
    field = np.empty((n, len(columns)), dtype=np.complex128)
    flux = np.empty_like(field)
    lags = np.arange(-2 * h.truncation, 2 * h.truncation + 1)
    for q, (m, beta) in enumerate(zip(columns, betas)):
        coefficients = h.series.plane_wave(float(beta), lags)
        values = coefficients[index - m + 2 * h.truncation]
        field[:, q] = values
        flux[:, q] = 1j * values * (kappa * kappa - alpha[m] * alpha) / beta
    return Modes(field, flux, np.asarray(betas, dtype=np.complex128))


@dataclass(frozen=True)
class Solution:
    """Reflected amplitudes per propagating order, as ``(E_z, H_z)``."""

    orders: NDArray[np.int64]
    e_amplitudes: NDArray[np.complex128]
    h_amplitudes: NDArray[np.complex128]
    efficiencies: NDArray[np.float64]
    absorption: float = 0.0
    condition: float = 1.0
    transmitted_orders: NDArray[np.int64] = _field(
        default_factory=lambda: np.zeros(0, dtype=np.int64)
    )
    transmitted_efficiencies: NDArray[np.float64] = _field(
        default_factory=lambda: np.zeros(0)
    )


@dataclass(frozen=True)
class _Vacuum:
    alpha: NDArray[np.float64]
    kappa: float
    live: NDArray[np.bool_]
    betas: NDArray[np.float64]
    incident: Modes
    outgoing: Modes
    evanescent: Modes

    @property
    def columns(self) -> Modes:
        return Modes(
            np.hstack([self.outgoing.field, self.evanescent.field]),
            np.hstack([self.outgoing.flux, self.evanescent.flux]),
            np.concatenate([self.outgoing.rho, self.evanescent.rho]),
        )


def _vacuum(h: Harmonics, k: float, sin_alpha: float, sin_gamma: float, period: float):
    kappa = k * sin_gamma
    alpha = -kappa * sin_alpha + 2.0 * np.pi * h.orders / period
    beta2 = kappa * kappa - alpha * alpha
    live = beta2 >= 0.0
    betas = np.sqrt(np.where(live, beta2, 0.0))
    zero = int(np.flatnonzero(h.orders == 0)[0])
    columns = np.flatnonzero(live)
    incident = plane_waves(h, alpha, kappa, [zero], [-betas[zero]])
    outgoing = plane_waves(h, alpha, kappa, columns, betas[columns])
    evanescent = lossless_evanescent(h, alpha, kappa * kappa, len(columns), +1)
    return _Vacuum(alpha, kappa, live, betas, incident, outgoing, evanescent)


def _efficiencies(vac: _Vacuum, e, h_):
    """``(beta_m / beta_0)(|E_m|^2 + |H_m|^2)`` per propagating order."""
    beta0 = -vac.incident.rho.real[0]
    weights = vac.betas[vac.live] / beta0
    return weights * (np.abs(e) ** 2 + np.abs(h_) ** 2)


def solve_perfect(
    h: Harmonics,
    *,
    wavelength: float,
    period: float,
    sin_alpha: float,
    sin_gamma: float,
    polarization: str,
) -> Solution:
    """Perfect conductor, any mount: Dirichlet (TE) or Neumann (TM).

    The vacuum side's transverse wavenumber is ``k sin(gamma)``, which is the
    exact conical reduction to the in-plane problem at ``lambda / sin(gamma)``.
    """
    k = 2.0 * np.pi / wavelength
    vac = _vacuum(h, k, sin_alpha, sin_gamma, period)
    cols = vac.columns
    if polarization == "TE":
        matrix, rhs = cols.field, -vac.incident.field[:, 0]
    else:
        matrix, rhs = cols.flux, -vac.incident.flux[:, 0]
    lu = _Factor(matrix)
    x = lu.solve(rhs)
    n_live = int(vac.live.sum())
    amplitudes = x[:n_live]
    zeros = np.zeros_like(amplitudes)
    e, hz = (amplitudes, zeros) if polarization == "TE" else (zeros, amplitudes)
    return Solution(
        orders=h.orders[vac.live],
        e_amplitudes=e,
        h_amplitudes=hz,
        efficiencies=_efficiencies(vac, e, hz),
        condition=lu.condition,
    )


def solve_finite(
    h: Harmonics,
    *,
    wavelength: float,
    period: float,
    sin_alpha: float,
    sin_gamma: float,
    cos_gamma: float,
    index: complex,
    incidents,
) -> list[Solution]:
    """One interface into an absorbing material, coupled ``(E_z, H_z)``.

    ``incidents`` is a list of ``(E_z, H_z)`` amplitude pairs; one
    factorization serves them all. Each solution carries the absorbed
    fraction, computed from the flux through the boundary rather than as
    ``1 - R``, so ``R + A = 1`` is a check and not an identity.
    """
    k = 2.0 * np.pi / wavelength
    vac = _vacuum(h, k, sin_alpha, sin_gamma, period)
    eps = complex(index) ** 2
    kz = k * cos_gamma
    kappa1_2 = vac.kappa**2
    kappa2_2 = k * k * eps - kz * kz
    if eps.imag == 0.0:
        # Lossless: the orders that propagate in the material are its real
        # eigenvalues, replaced by exact downgoing plane waves as on the
        # vacuum side. With eps < cos^2(gamma) (a negative permittivity,
        # say) kappa^2 is negative and nothing propagates.
        k2 = float(kappa2_2.real)
        beta2_sq = k2 - vac.alpha**2
        through = np.flatnonzero(beta2_sq > 0.0)
        betas2 = np.sqrt(beta2_sq[through])
        waves = plane_waves(h, vac.alpha, float(np.sqrt(max(k2, 0.0))), through, -betas2)
        rest = lossless_evanescent(h, vac.alpha, k2, len(through), -1)
        metal = Modes(
            np.hstack([waves.field, rest.field]),
            np.hstack([waves.flux, rest.flux]),
            np.concatenate([waves.rho, rest.rho]),
        )
    else:
        through, betas2 = np.zeros(0, int), np.zeros(0)
        metal = lossy_modes(h, vac.alpha, kappa2_2)

    n = h.size
    up = vac.columns
    tz = (kz / k) * 1j * vac.alpha[:, None]  # (k_z / k) d/dx, on rows
    ratio = kappa1_2 / kappa2_2
    zero = np.zeros((n, n), dtype=np.complex128)
    # Unknowns: [x_E, x_H, c_E, c_H]; rows: E_z, H_z, E_t, H_t, the last two
    # scaled by kappa_1^2 / k so every block is O(1) in the wavenumbers.
    system = np.block(
        [
            [up.field, zero, -metal.field, zero],
            [zero, up.field, zero, -metal.field],
            [tz * up.field, up.flux, -ratio * tz * metal.field, -ratio * metal.flux],
            [-up.flux, tz * up.field, ratio * eps * metal.flux, -ratio * tz * metal.field],
        ]
    )
    lu = _Factor(system)

    f_inc = vac.incident.field[:, 0]
    p_inc = vac.incident.flux[:, 0]
    beta0 = -vac.incident.rho.real[0]
    n_live = int(vac.live.sum())
    solutions = []
    for a, b in incidents:
        rhs = -np.concatenate(
            [
                a * f_inc,
                b * f_inc,
                (kz / k) * 1j * vac.alpha * a * f_inc + b * p_inc,
                -a * p_inc + (kz / k) * 1j * vac.alpha * b * f_inc,
            ]
        )
        x = lu.solve(rhs)
        x_e, x_h, c_e, c_h = np.split(x, 4)
        e, hz = x_e[:n_live], x_h[:n_live]
        # Flux into the material, from the material side's fields on u = 0.
        ez = metal.field @ c_e
        hz_b = metal.field @ c_h
        pe = metal.flux @ c_e
        ph = metal.flux @ c_h
        et = (1j / kappa2_2) * (k * ph + kz * 1j * vac.alpha * ez)
        ht = (1j / kappa2_2) * (kz * 1j * vac.alpha * hz_b - k * eps * pe)
        into = -np.real(np.vdot(ht, ez) - np.vdot(hz_b, et))
        incident_flux = k * beta0 * (abs(a) ** 2 + abs(b) ** 2) / kappa1_2
        n_through = len(through)
        t_e, t_h = c_e[:n_through], c_h[:n_through]
        transmitted = (
            k * betas2 * (eps.real * np.abs(t_e) ** 2 + np.abs(t_h) ** 2)
            / kappa2_2.real
            / incident_flux
            if n_through
            else np.zeros(0)
        )
        solutions.append(
            Solution(
                orders=h.orders[vac.live],
                e_amplitudes=e,
                h_amplitudes=hz,
                efficiencies=_efficiencies(vac, e, hz) / (abs(a) ** 2 + abs(b) ** 2),
                absorption=float(into / incident_flux),
                condition=lu.condition,
                transmitted_orders=h.orders[through],
                transmitted_efficiencies=transmitted,
            )
        )
    return solutions


class _Factor:
    """A square system, row-equilibrated, factorised once for several right sides.

    Rows are scaled to unit largest entry first, so the boundary conditions'
    mixed units (a field row is dimensionless, a flux row a wavenumber) do
    not masquerade as ill-conditioning. ``condition`` is LAPACK's 1-norm
    estimate of the equilibrated matrix: the measure of how many digits the
    solve loses, which on a deep profile is what fails first (module note).
    """

    def __init__(self, matrix):
        import scipy.linalg
        from scipy.linalg import lapack

        self._scale = 1.0 / np.abs(matrix).max(axis=1)
        scaled = matrix * self._scale[:, None]
        self._lu = scipy.linalg.lu_factor(scaled, check_finite=False)
        rcond, _ = lapack.zgecon(
            self._lu[0].astype(np.complex128), float(np.abs(scaled).sum(axis=0).max())
        )
        self.condition = float(np.inf if rcond == 0.0 else 1.0 / rcond)

    def solve(self, rhs):
        import scipy.linalg

        return scipy.linalg.lu_solve(self._lu, rhs * self._scale, check_finite=False)
