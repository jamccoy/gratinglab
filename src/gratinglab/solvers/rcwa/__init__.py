r"""RCWA: the Fourier modal method, written in-house.

The grating region is a stack of slabs, each piecewise constant along the
period. In each slab the field is expanded in ``2M + 1`` Floquet harmonics,
and Maxwell's equations become a matrix eigenproblem whose eigenvectors are
the slab's modes; the modes are joined across the slabs by a scattering
matrix. ``truncation`` is ``M``, the accuracy knob.

Four requirements from ``docs/roadmap.md`` section 1 shape it, each argued in
``docs/theory/rcwa.md``:

* **Li's factorization rules** (section 3): the Fourier matrix for ``Ex``,
  the component that jumps at a ridge wall, is formed by the inverse rule,
  and ``Ey`` and ``Ez`` by Laurent's. Without that, metallic TM converges
  badly, and a comparison plot would show the method's failure rather than
  any physics.
* **A scattering matrix** (section 4): Rumpf's symmetric S-matrix in a
  zero-thickness gap, never a plain transmittance matrix.
* **Vectorial and conical** (section 2): the full ``2N`` eigenproblem in
  ``(Ex, Ey)``, so off-plane incidence and TE-TM conversion are exact.
* **Layer stacks** (section 6): the core solves a :class:`Stack` with
  homogeneous layers above and below the grating region. ``Problem`` has no
  field for a layer yet, so ``solve`` builds vacuum / grating / substrate.

A ``Problem`` means here what it means to the integral solver's tabulated
mode: ridges of ``coating`` (or ``substrate``) in vacuum, on a semi-infinite
substrate of the same material. Efficiencies are absolute, and the scan's
``absorption`` is the power that enters the structure.

A lamellar groove is one exact slab. Any other single-valued profile is a
staircase of ``slices`` slabs, which by default follow ``truncation`` so the
one knob refines both; an undercut profile is refused.
"""

from __future__ import annotations

import time

import numpy as np
from numpy.typing import ArrayLike

from ...geometry import cos_beta, is_propagating, order_range, sin_beta
from ...illumination import Illumination
from ...problem import Problem
from ...result import EfficiencyScan, Provenance
from .._execution import capped_blas, check_blas_threads, check_workers, run_rows
from ..base import Capabilities, Progress, UnsupportedConfiguration, register
from ._core import incident_field, solve_stack
from ._stack import groove_geometry

__all__ = ["RcwaSolver", "rcwa"]

#: ``min |cos beta_m|`` below which an order is passing off: the external
#: region's ``V = Q kz^-1`` is nearly singular there (a Rayleigh anomaly).
#: The integral solver's threshold, so the two flag the same rows.
_ANOMALY_COS = 1e-3

#: Fewest slabs a staircase gets by default, however low the truncation.
_MIN_SLICES = 8


def _resolve_index(problem: Problem):
    """The material the ridges and the substrate are made of.

    ``coating`` first, ``substrate`` otherwise, as the integral solver's
    tabulated mode resolves it, so the same ``Problem`` is the same grating
    to both. Naming neither is refused: there is no permittivity to expand.
    """
    from ...materials import lookup

    name = problem.coating or problem.substrate
    if name is None:
        raise UnsupportedConfiguration(
            "rcwa needs a material: name one on Problem.coating (or "
            "substrate). The method expands the permittivity in a Fourier "
            "series, and without a material there is none. A perfect "
            "conductor is the integral solver's conductivity='perfect'."
        )
    return name, lookup(name)


class RcwaSolver:
    """Fourier modal method with Li's rules and an S-matrix, in any mount."""

    capabilities = Capabilities(
        name="rcwa",
        conical=True,
        polarizations=("TE", "TM", "unpolarized"),
        accuracy_knob="truncation",
        rigorous=True,
        handles_undercut=False,
        reports_progress=True,
    )

    def solve(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: ArrayLike,
        *,
        truncation: int = 20,
        slices: int | None = None,
        workers: int = 1,
        blas_threads: int | None = 1,
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        """Efficiency over a wavelength scan at fixed geometry.

        Parameters
        ----------
        truncation
            ``M``: harmonics ``-M..M`` are kept. This is the accuracy knob the
            convergence harness sweeps. It must reach every order the scan
            reports. Cost per wavelength is an eigenproblem of size
            ``2(2M+1)`` per patterned slab.
        slices
            Slabs in the staircase of a non-lamellar profile. ``None`` (the
            default) gives ``max(8, 2 * truncation)``, so sweeping
            ``truncation`` refines the staircase too; pass a number to pin
            it. Not consulted for a :class:`~gratinglab.profiles.Lamellar`
            groove, which is one exact slab.
        workers
            Threads solving wavelengths concurrently. ``1`` (the default) is
            the serial loop. Each wavelength is independent and LAPACK
            releases the GIL, so up to the core count this divides a scan's
            time nearly by ``workers``; the result is identical to the
            serial one. See :mod:`gratinglab.solvers._execution`.
        blas_threads
            The BLAS thread cap for the length of the call, restored on the
            way out. ``1`` (the default) is the measured fastest for these
            matrix sizes, by up to two orders of magnitude on a busy
            machine; ``None`` leaves the process's setting alone.
        progress
            Per-wavelength callback under the contract in
            :class:`~gratinglab.solvers.base.Solver`.
        """
        self.capabilities.check(problem, illumination)

        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
        if wavelengths.ndim != 1:
            raise ValueError("wavelengths must be one-dimensional")
        if (wavelengths <= 0).any():
            raise ValueError("wavelengths must be positive")
        if (
            isinstance(truncation, bool)
            or not isinstance(truncation, (int, np.integer))
            or truncation < 1
        ):
            raise ValueError(f"truncation must be a positive integer, got {truncation!r}")
        truncation = int(truncation)
        if slices is not None and (
            isinstance(slices, bool)
            or not isinstance(slices, (int, np.integer))
            or slices < 1
        ):
            raise ValueError(f"slices must be a positive integer or None, got {slices!r}")
        slice_count = (
            int(slices) if slices is not None else max(_MIN_SLICES, 2 * truncation)
        )
        workers = check_workers(workers)
        threads = check_blas_threads(blas_threads)

        started = time.perf_counter()
        sin_alpha = illumination.sin_alpha
        sin_gamma = illumination.sin_gamma
        cos_gamma = float(np.sqrt(1.0 - sin_gamma**2))
        cos_alpha = float(np.sqrt(1.0 - sin_alpha**2))

        all_orders = order_range(
            float(wavelengths.min()), problem.period, sin_alpha, sin_gamma
        )
        widest = int(all_orders[np.argmax(np.abs(all_orders))])
        if truncation < abs(widest):
            raise ValueError(
                f"truncation={truncation} keeps harmonics -{truncation}..{truncation}, "
                f"but order {widest} propagates at {wavelengths.min():.6g} nm; "
                f"need at least {abs(widest)}"
            )

        material, constants = _resolve_index(problem)
        indices = np.asarray(constants.n(wavelengths), dtype=np.complex128)
        if (indices.imag < 0.0).any():
            raise ValueError(
                f"material {material!r} tabulates a negative absorption (gain); "
                "a passive grating cannot be built from it"
            )

        groove = groove_geometry(problem.profile, problem.period, slice_count)

        polarizations = (
            ("TE", "TM")
            if illumination.polarization == "unpolarized"
            else (illumination.polarization,)
        )
        kx0 = -sin_alpha * sin_gamma
        kz0 = cos_alpha * sin_gamma
        incidents = [incident_field(p, kx0, cos_gamma, kz0) for p in polarizations]

        efficiency = np.zeros((len(wavelengths), len(all_orders)))
        absorption = np.zeros(len(wavelengths))
        propagating = np.zeros_like(efficiency, dtype=bool)
        anomalous: list[float] = []
        columns = all_orders + truncation  # harmonic m sits at index m + M

        def solve_row(row: int, _state: None):
            # Only the top plane: `absorbed` is all the scan records, and the
            # interior planes cost a product chain and a solve each.
            solution = solve_stack(
                groove.stack(complex(indices[row])),
                wavelength=float(wavelengths[row]),
                period=problem.period,
                kx0=kx0,
                ky=cos_gamma,
                incidents=incidents,
                truncation=truncation,
                planes="top",
            )
            return (
                solution.reflected[:, columns].mean(axis=0),
                float(solution.absorbed.mean()),
            )

        # (0, n) before any work, and a cancellation point ahead of every
        # wavelength (see scalar.py); run_rows keeps that with workers too.
        with capped_blas(threads):
            rows = run_rows(
                len(wavelengths),
                solve_row,
                state=lambda: None,
                workers=workers,
                progress=progress,
            )

        for row, (wavelength, (reflected, absorbed)) in enumerate(
            zip(wavelengths, rows)
        ):
            sines = sin_beta(
                all_orders, wavelength, problem.period, sin_alpha, sin_gamma
            )
            live = is_propagating(sines)
            if live.any() and np.abs(cos_beta(sines[live])).min() < _ANOMALY_COS:
                anomalous.append(float(wavelength))
            efficiency[row, live] = reflected[live]
            propagating[row, live] = True
            absorption[row] = absorbed

        warnings: list[str] = []
        if anomalous:
            span = f"{min(anomalous):.6g}-{max(anomalous):.6g} nm"
            warnings.append(
                f"an order passes off within |cos beta| < {_ANOMALY_COS:g} at "
                f"{len(anomalous)} wavelength(s) in {span} (Rayleigh anomaly); "
                "the external region's kz is nearly zero there and those rows "
                "are numerically delicate"
            )
        if (
            problem.coating is not None
            and problem.substrate is not None
            and problem.substrate != problem.coating
        ):
            warnings.append(
                f"substrate {problem.substrate!r} is not consulted: the grating "
                f"is modelled as {problem.coating!r} throughout, since Problem "
                "carries no coating thickness (docs/theory/rcwa.md section 6)"
            )
        if (
            not groove.exact
            and "TM" in polarizations
            and bool((indices.imag > 0.0).any())
        ):
            warnings.append(
                f"the profile is a staircase of {len(groove.slabs)} slabs; for an "
                "absorbing material in TM each step's corners carry field "
                "singularities the true profile lacks, so convergence in slices "
                "is slow and can be non-monotone (Popov, Neviere, Gralak & Tayeb "
                "2002; "
                "docs/theory/rcwa.md section 5)"
            )

        from ... import __version__

        notes: dict = {
            "harmonics": 2 * truncation + 1,
            "slices": len(groove.slabs),
            "profile": "exact (lamellar)" if groove.exact else "staircase",
            "factorization": "Li: inverse rule for Ex, Laurent for Ey and Ez",
            "propagation": "S-matrix (Rumpf 2011)",
            "systems": (
                "TE and TM decoupled (in-plane)"
                if cos_gamma == 0.0
                else "coupled 2N (conical)"
            ),
            "normalization": "absolute",
            "material": material,
            "energy_balance": (
                "R + A = 1 holds to rounding at every truncation (flux "
                "continuity is built into the S-matrix); it checks the "
                "implementation, not convergence"
            ),
        }

        return EfficiencyScan(
            wavelengths=wavelengths,
            orders=all_orders,
            efficiency=efficiency,
            propagating=propagating,
            provenance=Provenance(
                method="rcwa",
                version=__version__,
                truncation=truncation,
                converged=None,  # the convergence harness sets this
                wall_time_s=time.perf_counter() - started,
                warnings=tuple(warnings),
                notes=notes,
            ),
            absorption=absorption,
        )


#: The registered singleton.
rcwa = register(RcwaSolver())
