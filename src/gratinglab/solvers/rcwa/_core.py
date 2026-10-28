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

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
from numpy.typing import NDArray

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
    ``transmitted`` and ``plane_flux`` from the Poynting flux of the
    transverse fields, which holds in an absorbing medium too. ``plane_flux``
    is ``(n_incidents, n_planes)``: the downward power through each interface
    plane, top first, normalised by the incident power. The first plane is
    just below the superstrate, so ``reflected.sum() + plane_flux[:, 0] = 1``
    for every passive stack, and in a lossless one every plane carries the
    same flux. A solve asked for ``planes="top"`` has that first column only.
    """

    harmonics: NDArray[np.int64]
    kx: NDArray[np.float64]
    reflected: NDArray[np.float64]
    transmitted: NDArray[np.float64]
    plane_flux: NDArray[np.float64]

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
    kappa = _forward_branch(squared)
    return w, q @ w / kappa[None, :], kappa


def _forward_branch(squared) -> NDArray[np.complex128]:
    """``kappa`` from the eigenvalues ``kappa^2``, on the forward branch.

    Decaying downward, or propagating downward on the real axis. The
    tolerance keeps rounding in Im from flipping a lossless propagating mode.
    """
    kappa = np.sqrt(squared.astype(np.complex128))
    scale = np.abs(kappa)
    flip = (kappa.imag < -1e-12 * scale) | (
        (np.abs(kappa.imag) <= 1e-12 * scale) & (kappa.real < 0.0)
    )
    kappa = np.where(flip, -kappa, kappa)
    return np.where(np.abs(kappa) < _KZ_FLOOR, _KZ_FLOOR + 0j, kappa)


# In-plane (``ky = 0``) the ``2N`` problem is two ``N`` problems. ``PQ`` is
# block diagonal, so the modes are pure ``Ex`` (TM) or pure ``Ey`` (TE); ``Q``
# and every ``V`` are anti-diagonal, pairing ``Ex`` with ``Uy`` and ``Ey`` with
# ``Ux``; and so every product the S-matrix algebra forms (``V^-1 V0``,
# ``V0^-1 V``) is block diagonal. Each polarization then runs the same S-matrix
# code on its own ``N``-square blocks: an eigenproblem an eighth the cost and
# matrix products a quarter, and a polarization nothing excites is skipped.
# ``V`` below is the block that maps a system's ``E`` component to its ``U``.


@dataclass(frozen=True, slots=True)
class _System:
    """Which rows of ``[Ex, Ey]`` and ``[Ux, Uy]`` one decoupled system owns."""

    name: Literal["coupled", "TM", "TE"]
    e_rows: NDArray[np.int64]
    u_rows: NDArray[np.int64]


def _systems(size: int, decoupled: bool) -> tuple[_System, ...]:
    if not decoupled:
        rows = np.arange(2 * size)
        return (_System("coupled", rows, rows),)
    x, y = np.arange(size), np.arange(size, 2 * size)
    return (_System("TM", x, y), _System("TE", y, x))


def _homogeneous_v(system: _System, index: complex, kx, ky):
    """``(V, kappa)`` of a homogeneous medium for one system (``W = I``)."""
    if system.name == "coupled":
        _, v, kappa = _homogeneous_modes(index, kx, ky)
        return v, kappa
    epsilon = complex(index) ** 2
    kz = _kz(epsilon, kx, 0.0)
    if system.name == "TM":  # Uy = eps Ex / kz
        return np.diag(epsilon / kz), kz
    return np.diag((kx.astype(np.complex128) ** 2 - epsilon) / kz), kz  # Ux from Ey


def _gap_v(system: _System, kx, ky):
    """``V0`` of the gap (``kz = 1`` in every harmonic) for one system."""
    epsilon = 1.0 + kx**2 + ky**2
    if system.name == "coupled":
        return _q_homogeneous(epsilon, kx, ky)
    if system.name == "TM":
        return np.diag(epsilon.astype(np.complex128))
    return -np.eye(len(kx), dtype=np.complex128)


def _slab_modes_decoupled(system: _System, slab: Slab, harmonics, kx, factorization):
    """``(W, V, kappa)`` of one in-plane polarization in a slab, ``N``-square."""
    if slab.homogeneous:
        v, kappa = _homogeneous_v(system, slab.background, kx, 0.0)
        return np.eye(len(kx), dtype=np.complex128), v, kappa
    laurent = fourier_matrix(slab, harmonics, inverse=False)
    kx_ = np.diag(kx.astype(np.complex128))
    if system.name == "TE":  # Ey: (eps - kx^2) Ey = kappa^2 Ey
        q = kx_ @ kx_ - laurent
        squared, w = np.linalg.eig(-q)
    else:  # Ex: (1 - kx eps_zz^-1 kx) exx Ex = kappa^2 Ex
        exx = (
            fourier_matrix(slab, harmonics, inverse=True)
            if factorization == "li"
            else laurent
        )
        p = np.eye(len(kx)) - kx_ @ np.linalg.solve(laurent, kx_)
        q = exx
        squared, w = np.linalg.eig(p @ q)
    kappa = _forward_branch(squared)
    return w, q @ w / kappa[None, :], kappa


def _layer_s(w, v, kappa, thickness_k0: float, v0) -> _S:
    """Rumpf's symmetric S-matrix of a slab, in gap-mode amplitudes."""
    w_inv = np.linalg.inv(w)
    v_inv_v0 = np.linalg.solve(v, v0)
    a = w_inv + v_inv_v0
    b = w_inv - v_inv_v0
    x = np.diag(np.exp(1j * kappa * thickness_k0))
    a_inv = np.linalg.inv(a)
    xba = x @ b @ a_inv
    d = a - xba @ x @ b
    s11 = np.linalg.solve(d, xba @ x @ a - b)
    s12 = np.linalg.solve(d, x @ (a - b @ a_inv @ b))
    return _S(s11, s12, s12, s11)


def _external_s(v, v0_inv, *, top: bool) -> _S:
    """The connection of a semi-infinite homogeneous region (``W = I``)."""
    eye = np.eye(len(v))
    a = eye + v0_inv @ v
    b = eye - v0_inv @ v
    a_inv = np.linalg.inv(a)
    across = 0.5 * (a - b @ a_inv @ b)
    if top:
        return _S(-a_inv @ b, 2.0 * a_inv, across, b @ a_inv)
    return _S(b @ a_inv, across, 2.0 * a_inv, -a_inv @ b)


def _star(sa: _S, sb: _S) -> _S:
    """The Redheffer star product ``sa (x) sb``, ``sa`` above ``sb``."""
    eye = np.eye(len(sa.s11))
    left = np.linalg.inv(eye - sb.s11 @ sa.s22)
    right = np.linalg.inv(eye - sa.s22 @ sb.s11)
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
    planes: Literal["all", "top"] = "all",
    decouple: bool | None = None,
) -> StackSolution:
    """Reflected and transmitted power per harmonic for each incident field.

    ``kx0`` and ``ky`` are the incident wave vector's components normalised by
    ``k0`` (so they already carry the superstrate's index). ``incidents`` are
    full 3-vectors from :func:`incident_field`; only their transverse parts
    enter, in harmonic ``m = 0``. The superstrate must be lossless, since the
    incidence angles and the reflected efficiencies are defined in it.

    ``planes="top"`` forms ``plane_flux`` at the first plane only, which is
    all ``absorbed`` needs, and skips the products from the top and the solve
    at every interior plane. ``decouple`` splits an in-plane mount into its TE
    and TM systems; ``None`` does so whenever ``ky == 0``, and ``False``
    exists so a test can hold the split against the coupled ``2N`` solve.
    """
    if complex(stack.superstrate).imag != 0.0:
        raise ValueError("the superstrate must be lossless")
    if factorization not in ("li", "laurent"):
        raise ValueError(f"unknown factorization {factorization!r}")
    if planes not in ("all", "top"):
        raise ValueError(f"planes must be 'all' or 'top', got {planes!r}")

    size = 2 * truncation + 1
    harmonics = np.arange(-truncation, truncation + 1)
    kx = kx0 + harmonics * wavelength / period
    k0 = 2.0 * np.pi / wavelength

    # The gap: kz = 1 in every harmonic, so W0 = I and V0 = Q0.
    source = np.zeros((2 * size, len(incidents)), dtype=np.complex128)
    kz_ref = _kz(complex(stack.superstrate) ** 2, kx, ky)
    kz0 = float(np.real(kz_ref[truncation]))
    incident_power = np.empty(len(incidents))
    for column, field in enumerate(incidents):
        source[truncation, column] = field[0]
        source[size + truncation, column] = field[1]
        incident_power[column] = kz0 * float(np.sum(np.abs(field) ** 2))

    if decouple is None:
        decouple = ky == 0.0
    elif decouple and ky != 0.0:
        raise ValueError("decouple=True needs an in-plane mount (ky = 0)")

    # The full fields, filled system by system: each owns disjoint rows.
    r = np.zeros_like(source)
    t = np.zeros_like(source)
    u_trn = np.zeros_like(source)
    n_planes = len(stack.slabs) + 1 if planes == "all" else 1
    plane_e = np.zeros((n_planes, 2 * size, len(incidents)), dtype=np.complex128)
    plane_u = np.zeros_like(plane_e)

    for system in _systems(size, decouple):
        sub_source = source[system.e_rows]
        if not sub_source.any():
            continue  # in-plane and unexcited: its every amplitude is zero
        v0 = _gap_v(system, kx, ky)
        v0_inv = np.linalg.inv(v0)
        v_ref, _ = _homogeneous_v(system, stack.superstrate, kx, ky)
        v_trn, _ = _homogeneous_v(system, stack.substrate, kx, ky)

        blocks = [_external_s(v_ref, v0_inv, top=True)]
        for slab in stack.slabs:
            if system.name == "coupled":
                w, v, kappa = _slab_modes(slab, truncation, kx, ky, factorization)
            else:
                w, v, kappa = _slab_modes_decoupled(
                    system, slab, truncation, kx, factorization
                )
            blocks.append(_layer_s(w, v, kappa, slab.thickness * k0, v0))
        blocks.append(_external_s(v_trn, v0_inv, top=False))

        # Partial products from the bottom, so the field at a plane comes
        # from the S-matrices above and below it; from the top too only when
        # every plane is wanted.
        below = [blocks[-1]]
        for block in reversed(blocks[:-1]):
            below.append(_star(block, below[-1]))
        below.reverse()
        if planes == "all":
            above = [blocks[0]]
            for block in blocks[1:-1]:
                above.append(_star(above[-1], block))
        else:
            above = [blocks[0]]
        total = below[0]

        r[system.e_rows] = total.s11 @ sub_source
        sub_t = total.s21 @ sub_source
        t[system.e_rows] = sub_t
        u_trn[system.u_rows] = v_trn @ sub_t

        eye = np.eye(len(v0))
        for plane, (upper, lower) in enumerate(zip(above, below[1:])):
            forward = np.linalg.solve(
                eye - upper.s22 @ lower.s11, upper.s21 @ sub_source
            )
            backward = lower.s11 @ forward
            plane_e[plane][system.e_rows] = forward + backward
            plane_u[plane][system.u_rows] = v0 @ (forward - backward)

    # Reflected: in the lossless superstrate, flux is kz |E|^2 with the full
    # vector, so a propagating order is non-negative by construction.
    rx, ry = r[:size], r[size:]
    rz = (kx[:, None] * rx + ky * ry) / kz_ref[:, None]  # upward wave: k.E = 0
    reflected = (
        np.real(kz_ref)[:, None] * (np.abs(rx) ** 2 + np.abs(ry) ** 2 + np.abs(rz) ** 2)
    ).T / incident_power[:, None]

    transmitted = _flux(t, u_trn, size) / incident_power[:, None]

    fluxes = [
        _flux(e, u, size).sum(axis=1) / incident_power
        for e, u in zip(plane_e, plane_u)
    ]

    return StackSolution(
        harmonics=harmonics,
        kx=kx,
        reflected=reflected,
        transmitted=transmitted,
        plane_flux=np.array(fluxes).T,
    )
