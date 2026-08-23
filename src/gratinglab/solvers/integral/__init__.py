r"""The integral method: rigorous efficiency from a boundary integral equation.

Milestone 1 scope: a **perfectly conducting** grating in any mount. For a
perfect conductor the conical (off-plane) problem decouples exactly into an
in-plane problem at the reduced wavelength ``lambda / sin(gamma)`` (M&P
eq. 4.65; thesis Chapter-2.tex:931): TE stays a Dirichlet problem, TM a
Neumann problem, and cross-polarization conversion is identically zero. The
decoupling is what lets the off-plane X-ray case -- the reason this project
exists -- run in pure NumPy: the discretisation cost follows the *reduced*
wavelength, which at ``gamma ~ 1.5 deg`` is tens of nm against a
hundreds-of-nm period.

Efficiencies are **relative to a perfect reflector**: the sum over orders is
exactly 1 in the continuum (the energy-balance theorem, M&P eqs. 4.34/4.42;
thesis eq:prop_order_unity via Green's theorem). A named coating is ignored,
with a provenance warning -- this matches the perfect-conductivity PCGrate
reference tables the solver is validated against. Folding a Fresnel
reflectivity on top would repeat the error catalogued in
``docs/conventions.md`` section 10 item 5; absolute efficiency for real
materials is the finite-conductivity milestone (Goray & Schmidt 2010),
reachable through the ``conductivity`` option once implemented.

Known limitation, measured and documented rather than hidden: profile
corners make the TM solve converge at first order (the Meixner edge
condition meets a plain equal-arc mesh), so TM energy balance on cornered
profiles carries a percent-level deficit at moderate ``boundary_points``.
TE is untouched. The energy-balance deviation is reported in provenance;
a graded corner mesh is the known follow-up.
"""

from __future__ import annotations

import time
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from ...geometry import cos_beta, is_propagating, order_range, sin_beta
from ...illumination import Illumination
from ...problem import Problem
from ...result import EfficiencyScan, Provenance
from ..base import Capabilities, Progress, UnsupportedConfiguration, register
from ._boundary import physical_boundary
from ._core import solve_transverse

__all__ = ["IntegralSolver", "integral"]

#: Fewer nodes than this per transverse wavelength along the boundary cannot
#: resolve the field's oscillation; the solve is refused rather than let a
#: coarse rung return a plausible wrong number.
_POINTS_PER_WAVELENGTH = 6.0

#: ``min |cos beta_m|`` below which an order is passing off and the kernel's
#: ``1/gamma_m`` makes the solve numerically delicate (a Rayleigh anomaly).
_ANOMALY_COS = 1e-3


class IntegralSolver:
    """Boundary-integral (Nystrom) solver for perfectly conducting gratings."""

    capabilities = Capabilities(
        name="integral",
        conical=True,
        polarizations=("TE", "TM", "unpolarized"),
        accuracy_knob="boundary_points",
        rigorous=True,
        handles_undercut=True,
        reports_progress=True,
    )

    def solve(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: ArrayLike,
        *,
        boundary_points: int = 400,
        spectral_terms: int | None = None,
        conductivity: Literal["perfect"] = "perfect",
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        """Efficiency over a wavelength scan at fixed geometry.

        Parameters
        ----------
        boundary_points
            Number of equal-arc-length nodes on one groove period -- the
            accuracy knob the convergence harness sweeps. Cost is roughly
            cubic in it (kernel build plus dense solve).
        spectral_terms
            Truncation of the kernel's spectral sums. Defaults to
            ``max(boundary_points // 2, ...)`` so the single knob converges
            both discretisations together (M&P Table 4.1 converges at
            ``P ~ 2.2 M``); override only for truncation studies.
        conductivity
            ``"perfect"`` is the only implemented boundary condition. Any
            other value raises :class:`UnsupportedConfiguration`; the
            finite-conductivity formulation of Goray & Schmidt (2010) is the
            planned successor.
        progress
            Per-wavelength callback under the contract in
            :class:`~gratinglab.solvers.base.Solver`.
        """
        self.capabilities.check(problem, illumination)
        if conductivity != "perfect":
            raise UnsupportedConfiguration(
                f"conductivity={conductivity!r} is not implemented; only the "
                "perfectly conducting boundary condition is. The rigorous "
                "finite-conductivity treatment (Goray & Schmidt 2010) is a "
                "planned milestone."
            )

        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
        if wavelengths.ndim != 1:
            raise ValueError("wavelengths must be one-dimensional")
        if (wavelengths <= 0).any():
            raise ValueError("wavelengths must be positive")
        if boundary_points < 32:
            raise ValueError(
                f"boundary_points must be at least 32, got {boundary_points}"
            )

        started = time.perf_counter()
        sin_alpha = illumination.sin_alpha
        sin_gamma = illumination.sin_gamma

        boundary = physical_boundary(problem, boundary_points)

        # The field oscillates along the boundary at the *reduced* wavelength;
        # a mesh that cannot resolve it would not fail loudly on its own.
        reduced_min = float(wavelengths.min()) / sin_gamma
        needed = int(
            np.ceil(_POINTS_PER_WAVELENGTH * boundary.arc_length / reduced_min)
        )
        if boundary_points < needed:
            raise ValueError(
                f"boundary_points={boundary_points} gives fewer than "
                f"{_POINTS_PER_WAVELENGTH:g} nodes per transverse wavelength "
                f"({reduced_min:.4g} nm) along the {boundary.arc_length:.4g} nm "
                f"boundary; need at least {needed}"
            )

        all_orders = order_range(
            float(wavelengths.min()), problem.period, sin_alpha, sin_gamma
        )
        max_order = int(np.abs(all_orders).max())
        terms = (
            spectral_terms
            if spectral_terms is not None
            else max(boundary_points // 2, max_order + 8)
        )
        if terms <= max_order:
            raise ValueError(
                f"spectral_terms={terms} cannot represent order {max_order}; "
                "the kernel's spectral sum must extend past every propagating "
                "order"
            )

        polarizations = (
            ("TE", "TM")
            if illumination.polarization == "unpolarized"
            else (illumination.polarization,)
        )

        efficiency = np.zeros((len(wavelengths), len(all_orders)))
        propagating = np.zeros_like(efficiency, dtype=bool)
        warnings: list[str] = []
        anomalous: list[float] = []

        for row, wavelength in enumerate(wavelengths):
            if progress is not None:
                # Top of the row: (0, n) before any work, and a cancellation
                # point ahead of every wavelength (see scalar.py).
                progress(row, len(wavelengths))

            sines = sin_beta(
                all_orders, wavelength, problem.period, sin_alpha, sin_gamma
            )
            live = is_propagating(sines)
            if not live.any():
                continue
            if np.abs(cos_beta(sines[live])).min() < _ANOMALY_COS:
                anomalous.append(float(wavelength))

            reduced = float(wavelength) / sin_gamma
            values = np.zeros(int(live.sum()))
            for polarization in polarizations:
                solution = solve_transverse(
                    boundary,
                    wavelength=reduced,
                    period=problem.period,
                    sin_alpha=sin_alpha,
                    polarization=polarization,
                    terms=terms,
                )
                columns = np.searchsorted(all_orders[live], solution.orders)
                values[columns] += solution.efficiencies / len(polarizations)

            efficiency[row, live] = values
            propagating[row, live] = True

        if progress is not None:
            progress(len(wavelengths), len(wavelengths))

        if anomalous:
            span = f"{min(anomalous):.6g}-{max(anomalous):.6g} nm"
            warnings.append(
                f"an order passes off within |cos beta| < {_ANOMALY_COS:g} at "
                f"{len(anomalous)} wavelength(s) in {span} (Rayleigh anomaly); "
                "the kernel's 1/gamma_m is nearly singular there and those "
                "rows are numerically delicate"
            )
        if problem.coating is not None:
            warnings.append(
                f"coating {problem.coating!r} is not consulted: the perfectly "
                "conducting boundary condition has no material in it, and the "
                "sum over orders is relative to a perfect reflector "
                "(docs/conventions.md section 10 item 5)"
            )

        # Report, never rescale: the theorem says 1, the discretisation says
        # how close it got. Corners push TM off at first order.
        totals = np.where(propagating, efficiency, 0.0).sum(axis=1)
        computed = totals[propagating.any(axis=1)]
        deviation = float(np.abs(computed - 1.0).max()) if computed.size else 0.0
        if deviation > 0.005:
            warnings.append(
                f"summed efficiency strays from unity by up to "
                f"{100 * deviation:.2f}% (theorem: exactly 1). Profile corners "
                "limit the TM solve to first-order convergence on an equal-arc "
                "mesh; raise boundary_points, and read the deviation as the "
                "discretisation error estimate"
            )

        from ... import __version__

        return EfficiencyScan(
            wavelengths=wavelengths,
            orders=all_orders,
            efficiency=efficiency,
            propagating=propagating,
            provenance=Provenance(
                method="integral",
                version=__version__,
                truncation=boundary_points,
                converged=None,  # the convergence harness sets this
                wall_time_s=time.perf_counter() - started,
                warnings=tuple(warnings),
                notes={
                    "boundary_condition": "perfectly-conducting",
                    "spectral_terms": terms,
                    "energy_balance_deviation": deviation,
                    "normalization": "relative to a perfect reflector",
                },
            ),
        )


#: The registered singleton.
integral = register(IntegralSolver())
