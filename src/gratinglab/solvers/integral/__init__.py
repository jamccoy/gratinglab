r"""The integral method: rigorous efficiency from a boundary integral equation.

The solver supports two boundary conditions. With ``conductivity="perfect"``
(the default) the grating is a perfect conductor in any mount. The conical
(off-plane) problem then decouples exactly into an in-plane problem at the
reduced wavelength ``lambda / sin(gamma)`` (M&P eq. 4.65; thesis
Chapter-2.tex:931): TE remains a Dirichlet problem and TM a Neumann problem,
cross-polarization conversion is identically zero, and efficiencies are
relative to a perfect reflector, so that the sum over orders is exactly 1 in
the continuum (M&P eqs. 4.34/4.42; thesis eq:prop_order_unity). This
matches the convention of the perfect-conductivity reference tables against
which the solver is checked. A named coating is ignored, and a provenance
warning records this.

With ``conductivity="tabulated"`` the grating is a single interface into a
semi-infinite material whose complex index is read per wavelength from the
named ``coating`` (falling back to ``substrate``) through
``materials.lookup``. The problem is then the coupled conical system of
Goray & Schmidt (JOSA A 27, 585, 2010), assembled in ``_finite``.
Efficiencies are absolute. TE and TM couple through the cone angle, and each
solve also returns the absorbed fraction from an independent boundary
integral, so the recorded scan satisfies ``R + A = 1`` as a theorem. This
relation succeeds the perfect-conductor sum rule and is what
:func:`gratinglab.checks.check_energy_balance` checks. Folding a Fresnel
reflectivity onto the perfect-conductivity result is not equivalent, and
doing so remains the error catalogued in ``docs/conventions.md`` section 10
item 5.

Profile corners make the TM solve (and, with finite conductivity, both
polarizations) converge at first order on the default equal-arc mesh,
because the Meixner edge singularity of the densities is sampled by evenly
spaced nodes. Energy balance on cornered profiles then carries a
percent-level deficit at moderate ``boundary_points``. ``corner_grading=p``
clusters the nodes at the corners with Kress's sigmoidal transformation of
order ``p`` (M&P section 4.6.6) and restores fast convergence there; see
``docs/theory/integral.md`` section 5 for the measured ladders. It is opt-in
until it has been used on real problems, and on a profile with no corner it
changes nothing. It serves both boundary conditions. The finite-conductivity
system needs two corrections on it that the perfectly conducting one does
not: a singularity subtraction of ``K+`` at the corners
(``_nystrom.corner_subtraction``), without which TM on a visible-light metal
stalled, and the tangential derivative taken on the assembled trace
(``_finite``), without which the conical solve converged erratically.
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
from .._conditioning import ESTIMATOR, ConditionLog
from ..base import Capabilities, Progress, UnsupportedConfiguration, register
from ._boundary import physical_boundary
from ._core import solve_transverse
from ._finite import solve_finite_states
from ._nystrom import OperatorAssembler

if TYPE_CHECKING:  # pragma: no cover - import cycle guard for annotations
    from ...materials import OpticalConstants

__all__ = ["IntegralSolver", "integral"]

#: Fewer nodes than this per transverse wavelength along the boundary cannot
#: resolve the field's oscillation; the solve is refused rather than let a
#: coarse rung return a plausible wrong number. With finite conductivity the
#: same floor applies to the metal-side transverse wavelength, which is what
#: makes visible-light metals (|n| of a few) cost more nodes than X-rays.
_POINTS_PER_WAVELENGTH = 6.0

#: Deviation of the conserved quantity from unity -- the sum alone under the
#: perfect boundary, the sum plus absorption under the tabulated one -- above
#: which the discretisation is worth naming in provenance. Reported, never
#: rescaled: the theorem says 1, and how close the mesh got is the error
#: estimate.
_ENERGY_DEVIATION_WARN = 0.005

#: Highest ``corner_grading`` accepted. At order 8 the two nodes nearest a
#: corner already coincide in double precision at P = 512 on a 1000 nm
#: sawtooth, and the solve returns NaN (measured); 6 still separates them,
#: and 4 converges as well as 6 on every profile measured.
_MAX_GRADING = 6

#: The closest pair of graded nodes must be this many ulps of the largest
#: coordinate apart, so their separation, which the kernels take as a
#: difference of coordinates, carries at least two significant digits.
_MIN_SEPARATION_ULPS = 100.0

def _coincident_nodes(boundary, period: float) -> "tuple[float, float] | None":
    """``(closest, floor)`` in nm when two graded nodes cannot be told apart.

    ``None`` when the closest pair is separable, which an equal-arc mesh
    always is. Shared with the window's form (``gui.integral_options``) so
    the two refuse the same meshes.
    """
    if not boundary.graded:
        return None
    steps = np.hypot(
        np.diff(np.append(boundary.x, boundary.x[0] + period)),
        np.diff(np.append(boundary.y, boundary.y[0])),
    )
    scale = float(max(np.abs(boundary.x).max(), np.abs(boundary.y).max()) + period)
    floor = _MIN_SEPARATION_ULPS * np.finfo(np.float64).eps * scale
    closest = float(steps.min())
    return (closest, floor) if closest <= floor else None


#: ``min |cos beta_m|`` below which an order is passing off and the kernel's
#: ``1/gamma_m`` makes the solve numerically delicate (a Rayleigh anomaly).
_ANOMALY_COS = 1e-3


#: Condition estimate above which the solve warns. The second-kind TM and
#: coupled systems stay small, but TE's first-kind single layer grows with
#: the mesh, and a corner-graded mesh crowds nodes together: measured
#: 2026-10-27 (``docs/theory/integral.md`` section 9), ordinary problems
#: reach ``4.5e12`` (blazed, TE, graded, 3200 points, growing as ``P^4``),
#: and an order passing off pushes it as ``|lambda - lambda*|^-1/2`` to
#: ``9e15`` exactly at grazing. The answer was still right there, so the
#: warning says "unusual", not "wrong".
_WARN_KAPPA = 1e14


def _resolve_index(problem: Problem) -> "tuple[str, OpticalConstants]":
    """Resolve the material the tabulated boundary condition reads.

    ``coating`` is used first, and ``substrate`` otherwise. The model is a
    single interface into a semi-infinite material, so whichever is named is
    that material. Naming neither is refused, since there is no rigorous
    answer without an index. An unknown name raises ``UnknownMaterial``,
    listing what is available (the same policy as
    ``scalar._resolve_coating``).
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
        corner_grading: int = 0,
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        """Efficiency over a wavelength scan at fixed geometry.

        Parameters
        ----------
        boundary_points
            Number of equal-arc-length nodes on one groove period. This is
            the accuracy knob the convergence harness sweeps. Kernels are
            assembled as separable BLAS products (``_kernels``): a few
            ``(P, 2M+1)`` GEMMs per wavelength plus the dense solve, after
            a one-time ``O(P^2 M)`` geometry pass amortised over the scan.
            The tabulated boundary condition costs ~2x the perfect one per
            wavelength (more operators, sharing kernel evaluations per
            medium, and a doubled system off-plane).
        spectral_terms
            Truncation of the kernel's spectral sums. Defaults to
            ``max(boundary_points // 2, ...)`` so the single knob converges
            both discretisations together (M&P Table 4.1 converges at
            ``P ~ 2.2 M``), extended past the metal-side propagating cone
            for the tabulated condition; override only for truncation
            studies.
        conductivity
            ``"perfect"`` selects a perfectly conducting boundary, with
            efficiencies relative to a perfect reflector. ``"tabulated"``
            selects the finite-conductivity system of Goray & Schmidt (2010),
            with the complex index read from ``problem.coating`` (or
            ``substrate``) per wavelength; efficiencies are then absolute and
            absorption is recorded on the scan. Any other value raises
            :class:`UnsupportedConfiguration`.
        corner_grading
            ``0`` (the default) places the nodes at equal arc length.
            ``p`` from 2 to 6 grades them towards the profile's corners with
            Kress's sigmoidal transformation of order ``p``, which is what a
            cornered profile needs for TM to converge faster than first
            order. ``4`` is the recommended order. A profile with no corner
            (a vertex turning by more than
            :data:`gratinglab.profiles.CORNER_TURN_DEGREES`) keeps the
            equal-arc mesh whatever is passed. Both boundary conditions.
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
        if (
            isinstance(corner_grading, bool)
            or not isinstance(corner_grading, (int, np.integer))
            or not (corner_grading == 0 or 2 <= corner_grading <= _MAX_GRADING)
        ):
            raise ValueError(
                f"corner_grading must be 0 (equal arc length) or an integer "
                f"from 2 to {_MAX_GRADING}, got {corner_grading!r}"
            )
        corner_grading = int(corner_grading)

        started = time.perf_counter()
        sin_alpha = illumination.sin_alpha
        sin_gamma = illumination.sin_gamma
        # Exact zero in-plane (sin_gamma == 1), so the coupled system's
        # cross blocks vanish identically rather than to rounding.
        cos_gamma = float(np.sqrt(1.0 - sin_gamma**2))

        boundary = physical_boundary(
            problem, boundary_points, grading=corner_grading
        )
        coincident = _coincident_nodes(boundary, problem.period)
        if coincident is not None:
            closest, floor = coincident
            raise ValueError(
                f"corner_grading={corner_grading} at boundary_points="
                f"{boundary_points} puts two nodes {closest:.3g} nm "
                f"apart, below the {floor:.3g} nm that double precision "
                "can resolve on this groove; lower corner_grading (4 is "
                "recommended) or boundary_points"
            )

        # The field oscillates along the boundary at the *reduced* wavelength;
        # a mesh that cannot resolve it would not fail loudly on its own.
        reduced_min = float(wavelengths.min()) / sin_gamma
        # The coarsest spacing must resolve it, which on a graded mesh is
        # about twice the mean (resolution_length is arc_length otherwise).
        needed = int(
            np.ceil(
                _POINTS_PER_WAVELENGTH * boundary.resolution_length / reduced_min
            )
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
                    * boundary.resolution_length
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
        conditioning = ConditionLog()
        ill_conditioned: list[float] = []

        # One kernel assembly for the whole scan: the geometry-only half
        # (including the cached Kummer tails) is what every wavelength and
        # both media share; the per-wavelength memo is cleared each row.
        assembler = OperatorAssembler(
            boundary, period=problem.period, terms=terms
        )

        for row, wavelength in enumerate(wavelengths):
            if progress is not None:
                # Top of the row: (0, n) before any work, and a cancellation
                # point ahead of every wavelength (see scalar.py).
                progress(row, len(wavelengths))
            assembler.clear()

            sines = sin_beta(
                all_orders, wavelength, problem.period, sin_alpha, sin_gamma
            )
            live = is_propagating(sines)
            if not live.any():
                continue
            if np.abs(cos_beta(sines[live])).min() < _ANOMALY_COS:
                anomalous.append(float(wavelength))

            values = np.zeros(int(live.sum()))
            row_log = ConditionLog()
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
                        assembler=assembler,
                    )
                    row_log.record("nystrom", solution.condition)
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
                    assembler=assembler,
                )
                for solution in solutions:
                    row_log.record("nystrom", solution.condition)
                    columns = np.searchsorted(all_orders[live], solution.orders)
                    values[columns] += solution.efficiencies / len(polarizations)
                    absorption[row] += solution.absorption / len(polarizations)

            conditioning.merge(row_log)
            if row_log.worst.get("nystrom", 0.0) > _WARN_KAPPA:
                ill_conditioned.append(float(wavelength))
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
        conditioning_warning = conditioning.warning(
            "integral", ill_conditioned, _WARN_KAPPA
        )
        if conditioning_warning is not None:
            warnings.append(conditioning_warning)
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
        if deviation > _ENERGY_DEVIATION_WARN:
            warnings.append(
                f"summed efficiency{' plus absorption' if conductivity == 'tabulated' else ''} "
                f"strays from unity by up to {100 * deviation:.2f}% "
                "(theorem: exactly 1). "
                + (
                    "Profile corners limit convergence to first order on an "
                    "equal-arc mesh; pass corner_grading=4 or raise "
                    "boundary_points, and read the deviation as the "
                    "discretisation error estimate"
                    if not boundary.graded
                    else "Raise boundary_points, and read the deviation as "
                    "the discretisation error estimate"
                )
            )

        from ... import __version__

        notes: dict = {
            "spectral_terms": terms,
            "energy_balance_deviation": deviation,
            "condition": {
                "estimator": ESTIMATOR,
                "warn_above": _WARN_KAPPA,
                **conditioning.worst,
            },
            "mesh": (
                f"corner-graded, order {corner_grading}"
                if boundary.graded
                else "equal arc length"
            ),
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
