r"""The Fourier modal method on one :class:`Stack`, at one wavelength.

Everything is normalised by ``k0 = 2 pi / lambda``. The frame is the RCWA one
of ``docs/theory/rcwa.md`` section 2: ``x`` periodic, ``y`` along the grooves,
``z`` down into the structure. Fields are ``[Ex, Ey]`` and ``[Ux, Uy]`` with
``U = eta0 H``, each a ``2N`` vector over the harmonics ``m = -M..M``.

Propagation is Rumpf's symmetric scattering matrix in a zero-thickness gap
medium (PIER B 35, 241, 2011), section 4 of the theory page. The gap is
chosen with ``kz = 1`` in every harmonic, so ``W0 = I`` and ``V0 = Q0``.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
from numpy.typing import NDArray

from .._conditioning import ConditionLog
from ._stack import Slab, Stack, fourier_matrix

__all__ = ["StackSolution", "solve_stack", "incident_field"]

#: An external region's ``kz`` closer to zero than this is moved off it, so
#: ``V = Q kz^-1`` exists exactly at a Rayleigh anomaly. The backend flags any
#: row near one; this only keeps the measure-zero point finite.
_KZ_FLOOR = 1e-12


@dataclass(frozen=True, slots=True)
class _S:
    """One scattering matrix, four ``2N``-square blocks."""

    s11: NDArray[np.complex128]
    s12: NDArray[np.complex128]
    s21: NDArray[np.complex128]
    s22: NDArray[np.complex128]


@dataclass(frozen=True, slots=True)
class StackSolution:
    """Per-harmonic power fractions for each incident state.

    Arrays are ``(n_incidents, 2M+1)`` over ``harmonics`` unless noted.
    ``reflected`` is formed in the lossless superstrate from ``kz |E|^2``, and
    so is ``transmitted`` when the substrate is lossless, so a propagating
    order cannot come out negative by rounding. In an absorbing substrate
    ``transmitted`` is the Poynting flux of the transverse fields just below
    the interface, as is ``plane_flux`` everywhere. ``plane_flux``
    is ``(n_incidents, n_planes)``: the downward power through each interface
    plane, top first, normalised by the incident power. The first plane is
    just below the superstrate, so ``reflected.sum() + plane_flux[:, 0] = 1``
    for every passive stack, and in a lossless one every plane carries the
    same flux.

    ``condition`` is the worst 1-norm condition estimate per group of dense
    systems the solve factored (:mod:`gratinglab.solvers._conditioning`):
    ``"modes"`` for the slabs' mode matrices ``W`` and ``V``, ``"smatrix"``
    for every S-matrix denominator (each slab's ``A`` and ``D``, the external
    connections' ``A``, the star products' ``I - S S`` and the plane-flux
    solves).
    """

    harmonics: NDArray[np.int64]
    kx: NDArray[np.float64]
    reflected: NDArray[np.float64]
    transmitted: NDArray[np.float64]
    plane_flux: NDArray[np.float64]
    condition: dict[str, float] = dataclasses.field(default_factory=dict)

    @property
    def absorbed(self) -> NDArray[np.float64]:
        """Power entering the structure, per incident state."""
        return self.plane_flux[:, 0]


def _kz(epsilon: complex, kx: NDArray, ky: float) -> NDArray[np.complex128]:
    """``sqrt(eps - kx^2 - ky^2)`` on the branch that decays (or propagates) in +z."""
    kz = np.sqrt(complex(epsilon) - kx.astype(np.complex128) ** 2 - ky**2)
    flip = (kz.imag < 0.0) | ((kz.imag == 0.0) & (kz.real < 0.0))
    kz = np.where(flip, -kz, kz)
    return np.where(np.abs(kz) < _KZ_FLOOR, _KZ_FLOOR + 0j, kz)


def _q_homogeneous(epsilon, kx, ky) -> NDArray[np.complex128]:
    """``Q`` of a homogeneous medium of (diagonal) permittivity ``epsilon``."""
    kx = np.diag(kx.astype(np.complex128))
    eye = np.eye(len(kx))
    ky_ = ky * eye
    epsilon = epsilon * eye if np.isscalar(epsilon) else np.diag(epsilon)
    return np.block([[-kx @ ky_, kx @ kx - epsilon], [epsilon - ky_ @ ky_, ky_ @ kx]])


def _homogeneous_modes(index: complex, kx, ky):
    """``(W, V, kappa)`` of a homogeneous medium: plane waves, no eigensolver."""
    epsilon = complex(index) ** 2
    kz = _kz(epsilon, kx, ky)
    kappa = np.concatenate([kz, kz])
    q = _q_homogeneous(epsilon, kx, ky)
    return np.eye(len(kappa), dtype=np.complex128), q / kappa[None, :], kappa


def _slab_modes(slab: Slab, harmonics: int, kx, ky, factorization):
    """``(W, V, kappa)`` of a patterned slab, from the eigenproblem of ``PQ``.

    Li's rules (``rcwa.md`` section 3): ``Exx`` by the inverse rule, ``Eyy``
    and ``Ezz`` by Laurent's. ``factorization="laurent"`` uses Laurent's rule
    for all three; it exists so a test can show the difference, and the
    backend never passes it.
    """
    if slab.homogeneous:
        return _homogeneous_modes(slab.background, kx, ky)
    laurent = fourier_matrix(slab, harmonics, inverse=False)
    exx = (
        fourier_matrix(slab, harmonics, inverse=True)
        if factorization == "li"
        else laurent
    )
    ezz_inv = np.linalg.inv(laurent)
    kx_ = np.diag(kx.astype(np.complex128))
    eye = np.eye(len(kx))
    ky_ = ky * eye
    p = np.block(
        [
            [kx_ @ ezz_inv @ ky_, eye - kx_ @ ezz_inv @ kx_],
            [ky_ @ ezz_inv @ ky_ - eye, -ky_ @ ezz_inv @ kx_],
        ]
    )
    q = np.block([[-kx_ @ ky_, kx_ @ kx_ - laurent], [exx - ky_ @ ky_, ky_ @ kx_]])
    squared, w = np.linalg.eig(p @ q)
    kappa = np.sqrt(squared.astype(np.complex128))
    # The forward branch: decaying downward, or propagating downward on the
    # real axis. The tolerance keeps rounding in Im from flipping a lossless
    # propagating mode.
    scale = np.abs(kappa)
    flip = (kappa.imag < -1e-12 * scale) | (
        (np.abs(kappa.imag) <= 1e-12 * scale) & (kappa.real < 0.0)
    )
    kappa = np.where(flip, -kappa, kappa)
    kappa = np.where(np.abs(kappa) < _KZ_FLOOR, _KZ_FLOOR + 0j, kappa)
    return w, q @ w / kappa[None, :], kappa


def _layer_s(w, v, kappa, thickness_k0: float, v0, log: ConditionLog) -> _S:
    """Rumpf's symmetric S-matrix of a slab, in gap-mode amplitudes."""
    w_inv = log.inverse("modes", w)
    v_inv_v0 = log.solve("modes", v, v0)
    a = w_inv + v_inv_v0
    b = w_inv - v_inv_v0
    x = np.diag(np.exp(1j * kappa * thickness_k0))
    a_inv = log.inverse("smatrix", a)
    xba = x @ b @ a_inv
    d = a - xba @ x @ b
    both = log.solve(
        "smatrix", d, np.hstack([xba @ x @ a - b, x @ (a - b @ a_inv @ b)])
    )
    s11, s12 = both[:, : len(d)], both[:, len(d) :]
    return _S(s11, s12, s12, s11)


def _external_s(v, v0_inv, log: ConditionLog, *, top: bool) -> _S:
    """The connection of a semi-infinite homogeneous region (``W = I``)."""
    eye = np.eye(len(v))
    a = eye + v0_inv @ v
    b = eye - v0_inv @ v
    a_inv = log.inverse("smatrix", a)
    across = 0.5 * (a - b @ a_inv @ b)
    if top:
        return _S(-a_inv @ b, 2.0 * a_inv, across, b @ a_inv)
    return _S(b @ a_inv, across, 2.0 * a_inv, -a_inv @ b)


def _star(sa: _S, sb: _S, log: ConditionLog) -> _S:
    """The Redheffer star product ``sa (x) sb``, ``sa`` above ``sb``."""
    eye = np.eye(len(sa.s11))
    left = log.inverse("smatrix", eye - sb.s11 @ sa.s22)
    right = log.inverse("smatrix", eye - sa.s22 @ sb.s11)
    return _S(
        sa.s11 + sa.s12 @ left @ sb.s11 @ sa.s21,
        sa.s12 @ left @ sb.s12,
        sb.s21 @ right @ sa.s21,
        sb.s22 + sb.s21 @ right @ sa.s22 @ sb.s12,
    )


def _flux(e, u, size: int) -> NDArray[np.float64]:
    """Downward Poynting flux per harmonic, ``Re(Ex Uy* - Ey Ux*)``.

    ``e`` and ``u`` are ``(2N, k)``; returns ``(k, N)``.
    """
    ex, ey = e[:size], e[size:]
    ux, uy = u[:size], u[size:]
    return np.real(ex * np.conj(uy) - ey * np.conj(ux)).T


def incident_field(
    polarization: Literal["TE", "TM"], kx0: float, ky: float, kz0: float
) -> NDArray[np.complex128]:
    """Unit electric field of an incident plane wave, polarised by the grooves.

    TE has no groove component of ``U`` and TM none of ``E``
    (``docs/conventions.md`` section 7; grooves along ``y`` here). With ``k``
    the wave vector, TM's ``E`` is along ``k x y`` and TE's along
    ``(k x y) x k``.
    """
    k = np.array([kx0, ky, kz0], dtype=np.float64)
    tm = np.cross(k, [0.0, 1.0, 0.0])
    field = np.cross(tm, k) if polarization == "TE" else tm
    return (field / np.linalg.norm(field)).astype(np.complex128)


def solve_stack(
    stack: Stack,
    *,
    wavelength: float,
    period: float,
    kx0: float,
    ky: float,
    incidents: Sequence[NDArray[np.complex128]],
    truncation: int,
    factorization: Literal["li", "laurent"] = "li",
) -> StackSolution:
    """Reflected and transmitted power per harmonic for each incident field.

    ``kx0`` and ``ky`` are the incident wave vector's components normalised by
    ``k0`` (so they already carry the superstrate's index). ``incidents`` are
    full 3-vectors from :func:`incident_field`; only their transverse parts
    enter, in harmonic ``m = 0``. The superstrate must be lossless, since the
    incidence angles and the reflected efficiencies are defined in it.
    """
    if complex(stack.superstrate).imag != 0.0:
        raise ValueError("the superstrate must be lossless")
    if factorization not in ("li", "laurent"):
        raise ValueError(f"unknown factorization {factorization!r}")

    size = 2 * truncation + 1
    harmonics = np.arange(-truncation, truncation + 1)
    kx = kx0 + harmonics * wavelength / period
    k0 = 2.0 * np.pi / wavelength

    # The gap: kz = 1 in every harmonic, so W0 = I and V0 = Q0.
    v0 = _q_homogeneous(1.0 + kx**2 + ky**2, kx, ky)
    v0_inv = np.linalg.inv(v0)

    _, v_ref, kappa_ref = _homogeneous_modes(stack.superstrate, kx, ky)
    _, v_trn, kappa_trn = _homogeneous_modes(stack.substrate, kx, ky)

    log = ConditionLog()
    blocks = [_external_s(v_ref, v0_inv, log, top=True)]
    for slab in stack.slabs:
        w, v, kappa = _slab_modes(slab, truncation, kx, ky, factorization)
        blocks.append(_layer_s(w, v, kappa, slab.thickness * k0, v0, log))
    blocks.append(_external_s(v_trn, v0_inv, log, top=False))

    # Partial products from each end, so the field at every interior plane
    # comes from the S-matrices above and below it.
    above = [blocks[0]]
    for block in blocks[1:]:
        above.append(_star(above[-1], block, log))
    below = [blocks[-1]]
    for block in reversed(blocks[:-1]):
        below.append(_star(block, below[-1], log))
    below.reverse()
    total = above[-1]

    source = np.zeros((2 * size, len(incidents)), dtype=np.complex128)
    kz0 = float(np.real(kappa_ref[truncation]))
    incident_power = np.empty(len(incidents))
    for column, field in enumerate(incidents):
        source[truncation, column] = field[0]
        source[size + truncation, column] = field[1]
        incident_power[column] = kz0 * float(np.sum(np.abs(field) ** 2))

    # Reflected: in the lossless superstrate, flux is kz |E|^2 with the full
    # vector, so a propagating order is non-negative by construction.
    r = total.s11 @ source
    rx, ry = r[:size], r[size:]
    kz_ref = kappa_ref[:size]
    rz = (kx[:, None] * rx + ky * ry) / kz_ref[:, None]  # upward wave: k.E = 0
    reflected = (
        np.real(kz_ref)[:, None] * (np.abs(rx) ** 2 + np.abs(ry) ** 2 + np.abs(rz) ** 2)
    ).T / incident_power[:, None]

    t = total.s21 @ source
    if complex(stack.substrate).imag == 0.0:
        tx, ty = t[:size], t[size:]
        kz_trn = kappa_trn[:size]
        tz = -(kx[:, None] * tx + ky * ty) / kz_trn[:, None]  # downward: k.E = 0
        transmitted = (
            np.real(kz_trn)[:, None]
            * (np.abs(tx) ** 2 + np.abs(ty) ** 2 + np.abs(tz) ** 2)
        ).T / incident_power[:, None]
    else:
        transmitted = _flux(t, v_trn @ t, size) / incident_power[:, None]

    planes = []
    eye = np.eye(2 * size)
    for upper, lower in zip(above[:-1], below[1:]):
        forward = log.solve("smatrix", eye - upper.s22 @ lower.s11, upper.s21 @ source)
        backward = lower.s11 @ forward
        flux = _flux(forward + backward, v0 @ (forward - backward), size)
        planes.append(flux.sum(axis=1) / incident_power)

    return StackSolution(
        harmonics=harmonics,
        kx=kx,
        reflected=reflected,
        transmitted=transmitted,
        plane_flux=np.array(planes).T,
        condition=dict(log.worst),
    )
