r"""The integral method: rigorous efficiency from a boundary integral equation.

Two boundary conditions, one solver. With ``conductivity="perfect"`` (the
default) the grating is a **perfect conductor** in any mount: the conical
(off-plane) problem decouples exactly into an in-plane problem at the
reduced wavelength ``lambda / sin(gamma)`` (M&P eq. 4.65; thesis
Chapter-2.tex:931), TE stays a Dirichlet problem, TM a Neumann problem,
cross-polarization conversion is identically zero, and efficiencies are
**relative to a perfect reflector** -- the sum over orders is exactly 1 in
the continuum (M&P eqs. 4.34/4.42; thesis eq:prop_order_unity). This
matches the perfect-conductivity PCGrate reference tables the solver is
validated against; a named coating is ignored, with a provenance warning.

With ``conductivity="tabulated"`` the grating is one interface into a
semi-infinite material whose complex index is read per wavelength from the
named ``coating`` (falling back to ``substrate``) through
``materials.lookup`` -- the coupled conical system of Goray & Schmidt
(JOSA A 27, 585, 2010), assembled in ``_finite``. Efficiencies are
**absolute**: TE and TM couple through the cone angle, and each solve also
returns the absorbed fraction from an independent boundary integral, so
the recorded scan satisfies ``R + A = 1`` as a theorem (the successor of
the perfect-conductor sum rule, and what
:func:`gratinglab.checks.check_energy_balance` then checks). Folding a
Fresnel reflectivity onto the perfect-conductivity result is **not**
equivalent and remains the error catalogued in ``docs/conventions.md``
section 10 item 5.

Known limitation, measured and documented rather than hidden: profile
corners make the TM solve (and with finite conductivity, both
polarizations) converge at first order -- the Meixner edge condition meets
a plain equal-arc mesh -- so energy balance on cornered profiles carries a
percent-level deficit at moderate ``boundary_points``. The deviation is
reported in provenance; a graded corner mesh is the known follow-up.
"""

from __future__ import annotations

import time
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import ArrayLike

from ...geometry import cos_beta, is_propagating, order_range, sin_beta
from ...illumination import Illumination
from ...problem import Problem
from ...result import EfficiencyScan, Provenance
from ..base import Capabilities, Progress, UnsupportedConfiguration, register
from ._boundary import physical_boundary
from ._core import solve_transverse
from ._finite import solve_finite_states

if TYPE_CHECKING:  # pragma: no cover - import cycle guard for annotations
    from ...materials import OpticalConstants

__all__ = ["IntegralSolver", "integral"]

#: Fewer nodes than this per transverse wavelength along the boundary cannot
#: resolve the field's oscillation; the solve is refused rather than let a
#: coarse rung return a plausible wrong number. With finite conductivity the
#: same floor applies to the metal-side transverse wavelength, which is what
#: makes visible-light metals (|n| of a few) cost more nodes than X-rays.
_POINTS_PER_WAVELENGTH = 6.0

#: ``min |cos beta_m|`` below which an order is passing off and the kernel's
#: ``1/gamma_m`` makes the solve numerically delicate (a Rayleigh anomaly).
_ANOMALY_COS = 1e-3


def _resolve_index(problem: Problem) -> "tuple[str, OpticalConstants]":
    """Resolve the material the tabulated boundary condition reads.

    ``coating`` first, else ``substrate`` -- the model is a single interface
    into a semi-infinite material, so whichever is named is that material.
    Naming neither is a refusal, not a default: there is no rigorous answer
    without an index. An unknown name raises ``UnknownMaterial`` listing
    what is available (same policy as ``scalar._resolve_coating``).
    """
    from ...materials import lookup

    name = problem.coating or problem.substrate
    if name is None:
        raise UnsupportedConfiguration(
            'conductivity="tabulated" needs a material: name one on '
            "Problem.coating (or substrate). The finite-conductivity "
            "boundary condition reads its complex index per wavelength; "
            "without a material there is nothing rigorous to compute."
        )
    return name, lookup(name)


class IntegralSolver:
    """Boundary-integral (Nystrom) solver for gratings in any mount.

    Perfectly conducting by default; ``conductivity="tabulated"`` switches
    to the finite-conductivity coupled system of Goray & Schmidt (2010)
    with the material read from the problem's coating.
    """

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
        conductivity: Literal["perfect", "tabulated"] = "perfect",
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        """Efficiency over a wavelength scan at fixed geometry.

        Parameters
        ----------
        boundary_points
            Number of equal-arc-length nodes on one groove period -- the
            accuracy knob the convergence harness sweeps. Cost is roughly
            cubic in it (kernel build plus dense solve); the tabulated
            boundary condition costs ~5-8x the perfect one per wavelength
            (four kernels, operator products, a doubled system off-plane).
        spectral_terms
            Truncation of the kernel's spectral sums. Defaults to
            ``max(boundary_points // 2, ...)`` so the single knob converges
            both discretisations together (M&P Table 4.1 converges at
            ``P ~ 2.2 M``), extended past the metal-side propagating cone
            for the tabulated condition; override only for truncation
            studies.
        conductivity
            ``"perfect"`` -- perfectly conducting boundary, efficiencies
            relative to a perfect reflector. ``"tabulated"`` -- the
            finite-conductivity system of Goray & Schmidt (2010) with the
            complex index read from ``problem.coating`` (or ``substrate``)
            per wavelength; efficiencies absolute, absorption recorded on
            the scan. Anything else raises
            :class:`UnsupportedConfiguration`.
        progress
            Per-wavelength callback under the contract in
            :class:`~gratinglab.solvers.base.Solver`.
        """
        self.capabilities.check(problem, illumination)
        if conductivity not in ("perfect", "tabulated"):
            raise UnsupportedConfiguration(
                f"conductivity={conductivity!r} is not implemented: "
                '"perfect" is the perfectly conducting boundary and '
                '"tabulated" the finite-conductivity system of Goray & '
                "Schmidt (2010). Intermediate models (e.g. a Leontovich "
                "impedance condition) are not offered -- they are invalid "
                "in the soft-X-ray regime this solver targets."
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
        # Exact zero in-plane (sin_gamma == 1), so the coupled system's
        # cross blocks vanish identically rather than to rounding.
        cos_gamma = float(np.sqrt(1.0 - sin_gamma**2))

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

        material: str | None = None
        indices = None
        metal_cone = 0
        if conductivity == "tabulated":
            material, constants = _resolve_index(problem)
            indices = np.asarray(constants.n(wavelengths), dtype=np.complex128)
            if (indices.imag < 0.0).any():
                raise ValueError(
                    f"material {material!r} tabulates a negative absorption "
                    "(gain); the transmission problem loses uniqueness there "
                    "(Goray & Schmidt 2010, section 2.B)"
                )
            # The same resolution floor on the metal side: the transverse
            # wavelength inside the material is lambda / |sqrt(n^2 -
            # cos^2 gamma)|, short for visible-light metals.
            factors = np.abs(np.sqrt(indices**2 - cos_gamma**2))
            metal_reduced_min = float((wavelengths / factors).min())
            needed_metal = int(
                np.ceil(
                    _POINTS_PER_WAVELENGTH
                    * boundary.arc_length
                    / metal_reduced_min
                )
            )
            if boundary_points < needed_metal:
                raise ValueError(
                    f"boundary_points={boundary_points} gives fewer than "
                    f"{_POINTS_PER_WAVELENGTH:g} nodes per *metal-side* "
                    f"transverse wavelength ({metal_reduced_min:.4g} nm for "
                    f"{material}); need at least {needed_metal}. The field "
                    "inside a conductive material oscillates and decays on "
                    "the scale lambda/|n|, and the mesh must resolve it"
                )
            # The kernel's spectral sum must also clear the metal-side
            # propagating cone |k_t^-| d / 2 pi.
            metal_cone = int(
                np.ceil(
                    (factors / wavelengths).max() * problem.period
                )
            )

        all_orders = order_range(
            float(wavelengths.min()), problem.period, sin_alpha, sin_gamma
        )
        max_order = int(np.abs(all_orders).max())
        terms = (
            spectral_terms
            if spectral_terms is not None
            else max(boundary_points // 2, max_order + 8, metal_cone + 8)
        )
        if terms <= max_order or terms <= metal_cone:
            raise ValueError(
                f"spectral_terms={terms} cannot represent order "
                f"{max(max_order, metal_cone)}; the kernel's spectral sum "
                "must extend past every propagating order on both sides of "
                "the boundary"
            )

        polarizations = (
            ("TE", "TM")
            if illumination.polarization == "unpolarized"
            else (illumination.polarization,)
        )
        incidents = {"TE": (1.0 + 0j, 0.0 + 0j), "TM": (0.0 + 0j, 1.0 + 0j)}

        efficiency = np.zeros((len(wavelengths), len(all_orders)))
        absorption = np.zeros(len(wavelengths))
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

            values = np.zeros(int(live.sum()))
            if conductivity == "perfect":
                reduced = float(wavelength) / sin_gamma
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
            else:
                solutions = solve_finite_states(
                    boundary,
                    wavelength=float(wavelength),
                    period=problem.period,
                    sin_alpha=sin_alpha,
                    index=complex(indices[row]),
                    cos_gamma=cos_gamma,
                    incidents=[incidents[p] for p in polarizations],
                    terms=terms,
                )
                for solution in solutions:
                    columns = np.searchsorted(all_orders[live], solution.orders)
                    values[columns] += solution.efficiencies / len(polarizations)
                    absorption[row] += solution.absorption / len(polarizations)

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
        if conductivity == "perfect" and problem.coating is not None:
            warnings.append(
                f"coating {problem.coating!r} is not consulted: the perfectly "
                "conducting boundary condition has no material in it, and the "
                "sum over orders is relative to a perfect reflector "
                "(docs/conventions.md section 10 item 5)"
            )

        # Report, never rescale: the theorem says 1 -- for the sum alone
        # (perfect) or the sum plus absorption (tabulated) -- and the
        # discretisation says how close it got. Corners push the solve off
        # at first order.
        totals = np.where(propagating, efficiency, 0.0).sum(axis=1)
        if conductivity == "tabulated":
            totals = totals + absorption
        computed = totals[propagating.any(axis=1)]
        deviation = float(np.abs(computed - 1.0).max()) if computed.size else 0.0
        if deviation > 0.005:
            warnings.append(
                f"summed efficiency{' plus absorption' if conductivity == 'tabulated' else ''} "
                f"strays from unity by up to {100 * deviation:.2f}% "
                "(theorem: exactly 1). Profile corners limit convergence to "
                "first order on an equal-arc mesh; raise boundary_points, and "
                "read the deviation as the discretisation error estimate"
            )

        from ... import __version__

        notes: dict = {
            "spectral_terms": terms,
            "energy_balance_deviation": deviation,
        }
        if conductivity == "perfect":
            notes["boundary_condition"] = "perfectly-conducting"
            notes["normalization"] = "relative to a perfect reflector"
        else:
            notes["boundary_condition"] = (
                "finite-conductivity (Goray & Schmidt 2010)"
            )
            notes["normalization"] = "absolute"
            notes["material"] = material

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
                notes=notes,
            ),
            absorption=absorption if conductivity == "tabulated" else None,
        )


#: The registered singleton.
integral = register(IntegralSolver())
