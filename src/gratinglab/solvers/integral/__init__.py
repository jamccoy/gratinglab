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
solve also returns the power crossing into the material from an independent
boundary integral. Into a lossless material that power leaves as
transmitted orders, recorded as the scan's ``transmission``, and the
absorption keeps only the remainder, which is then the discretisation error;
so the recorded scan satisfies ``R + T + A = 1`` as a theorem. This
relation succeeds the perfect-conductor sum rule and is what
:func:`gratinglab.checks.check_energy_balance` checks. Folding a Fresnel
reflectivity onto the perfect-conductivity result is not equivalent, and
doing so remains the error catalogued in ``docs/conventions.md`` section 10
item 5.

Profile corners make the TM solve (and, with finite conductivity, both
polarizations) converge at first order on an equal-arc mesh,
because the Meixner edge singularity of the densities is sampled by evenly
spaced nodes. Energy balance on cornered profiles then carries a
percent-level deficit at moderate ``boundary_points``. ``corner_grading=p``
clusters the nodes at the corners with Kress's sigmoidal transformation of
order ``p`` (M&P section 4.6.6) and restores fast convergence there; see
``docs/theory/integral.md`` section 5 for the measured ladders. It is the
default, at order 4, since 0.2.7: it was opt-in until then, and the default
equal-arc mesh returned a summed efficiency plus absorption of 2.95 on the
GUI's own default blaze in soft X-rays (3.5 nm, Au, 400 points), with order +2
15% below the C-method, where order 4 closes to unity and matches it to
2e-4. On a profile with no corner it changes nothing, bit for bit, so
``corner_grading=0`` is needed only to reproduce an equal-arc result. It
serves both boundary conditions. The finite-conductivity
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
from ...profiles import (
    CORNER_TURN_DEGREES,
    FromProfileData,
    polyline_corner_count,
    polyline_rounded_turn,
)
from ...result import EfficiencyScan, Provenance, Transmission
from .._conditioning import ESTIMATOR, ConditionLog
from .._execution import capped_blas, check_blas_threads, check_workers, run_rows
from ..base import Capabilities, Progress, UnsupportedConfiguration, register
from ._boundary import physical_boundary, smooths
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

#: ``profile_smoothing="auto"`` low-passes a measured polygon with a Gaussian
#: whose standard deviation is the shortest transverse wavelength in the
#: scan over this. Measured 2026-10-29 on an AFM groove (2018 TASTE scan,
#: 2000 points, 400 nm period, 15.5 nm at |n| = 0.888, where that wavelength
#: is 30.5 nm): 2 nm costs the C-method's own R 3.5e-6 against the polygon,
#: and the integral solver's R + A - 1 falls 1.7e-3, 1.1e-4, 6.4e-6 at 256,
#: 512, 1024 points. 0.5 nm still stood at 2.4e-3 at 1024, because what the
#: mesh must resolve is the finest curvature the smoothing leaves.
_SMOOTHING_PER_WAVELENGTH = 16.0

#: A measured groove whose sub-threshold vertices turn by more than this many
#: degrees within this many nm of arc gets a provenance warning: a corner
#: rounded into slight bends, which smoothing cannot help and the graded mesh
#: cannot see. Measured 2026-10-29 on four grooves: the three that converge
#: carry at most 16.7 deg in any 5 nm, and the one that does not (the
#: repository's synthetic tip-dilated blazed scan through metrology.boundary,
#: whose valley turns 104 deg over +-5 nm) carries 74. A heuristic separating
#: those cases, not a convergence theorem.
_ROUNDED_WINDOW_NM = 5.0
_ROUNDED_TURN_DEGREES = 45.0

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


def _profile_note(problem: Problem, boundary, conductivity, requested) -> dict:
    """What the integral solver meshed for a measured groove, and why.

    Recorded on every :class:`FromProfileData` solve, so a result can be
    compared against PCGrate's linear interpolation knowing which geometry
    produced it.
    """
    profile = problem.profile
    vertices = len(profile.t)
    if boundary.smoothed is not None:
        rule = (
            f"auto: shortest transverse wavelength / {_SMOOTHING_PER_WAVELENGTH:g}"
            if requested == "auto"
            else "given"
        )
        note = {
            "representation": "smoothed",
            "gaussian_sigma_nm": boundary.smoothed,
            "rule": rule,
            "vertices": vertices,
        }
        corners = polyline_corner_count(profile.t, profile.y)
        if corners:
            note["representation"] = (
                "smoothed between its corners, which are kept and graded"
            )
            note["corners"] = corners
        return note
    if conductivity != "tabulated":
        reason = "the perfectly conducting condition meshes the polygon as given"
    elif requested is None:
        reason = "profile_smoothing=None"
    elif not smooths(problem):
        reason = (
            "nothing to smooth: every vertex is straight or turns by more "
            f"than {CORNER_TURN_DEGREES:g} deg, a corner the mesh grades"
        )
    else:  # pragma: no cover - every route to an unsmoothed mesh is above
        reason = "not smoothed"
    return {
        "representation": "polygon, joined by straight segments as PCGrate joins it",
        "reason": reason,
        "vertices": vertices,
    }


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
        corner_grading: int = 4,
        profile_smoothing: "float | Literal['auto'] | None" = None,
        workers: int = 1,
        blas_threads: int | None = 1,
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
            ``p`` from 2 to 6 grades the nodes towards the profile's corners
            with Kress's sigmoidal transformation of order ``p``, which is
            what a cornered profile needs for TM (and, with finite
            conductivity, both polarizations) to converge faster than first
            order. ``4``, the default since 0.2.7, is the recommended order.
            ``0`` places the nodes at equal arc length, the default before
            0.2.7. A profile with no corner
            (a vertex turning by more than
            :data:`gratinglab.profiles.CORNER_TURN_DEGREES`) keeps the
            equal-arc mesh whatever is passed. Both boundary conditions.
        profile_smoothing
            How a measured groove (:class:`~gratinglab.profiles.FromProfileData`)
            is meshed under ``conductivity="tabulated"``. Every point of
            such a polygon that bends by less than
            :data:`gratinglab.profiles.CORNER_TURN_DEGREES` is a slight
            corner of the metal interface, too slight to grade towards and
            too many to resolve, and meshing them as given wrecked the
            independently computed absorption, by 8% of the incident power
            on a measured groove at 80 eV, while leaving the reflected
            orders nearly right (issue #6). ``"auto"`` low-passes the
            polygon along its arc length by a Gaussian of standard
            deviation one sixteenth of the shortest transverse wavelength
            in the scan, either side of the interface, and a number sets
            that width in nm. With no real corner the whole curve is
            smoothed; with real corners each side between them is, with the
            corners kept where they are and graded as before.

            ``None``, the default, meshes the polygon as given, joined by
            straight segments as PCGrate joins it, bit for bit as before.
            Smoothing is opt-in because its width is only negligible where
            it was measured: in the X-ray grazing regime the fault lives in
            it cost R 3.5e-6, but on a finely sampled smooth groove at
            500 nm "auto" moved the efficiencies by up to 4e-3 where the
            polygon as given was already right. Reach for it when the
            energy-balance warning fires on a measured groove, and check
            the result against the C-method. A polygon of corners and
            straight sides, every analytic profile and the perfectly
            conducting condition are untouched whatever is passed.
            Provenance records what was meshed under ``notes["profile"]``.
            See ``docs/theory/integral.md`` section 8.
        workers
            Threads solving wavelengths concurrently. ``1`` (the default) is
            the serial loop. Each thread has its own kernel memo over the
            scan's shared geometry, and LAPACK releases the GIL, so up to the
            core count this divides a scan's time nearly by ``workers``; the
            result is identical to the serial one. See
            :mod:`gratinglab.solvers._execution`.
        blas_threads
            The BLAS thread cap for the length of the call, restored on the
            way out. ``1`` (the default) is the measured fastest at these
            matrix sizes on a busy machine; ``None`` leaves the process's
            setting alone, which can win for one very large solve on an idle
            machine.
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
        if not (
            profile_smoothing is None
            or profile_smoothing == "auto"
            or (
                isinstance(profile_smoothing, (int, float, np.integer, np.floating))
                and not isinstance(profile_smoothing, bool)
                and np.isfinite(profile_smoothing)
                and profile_smoothing > 0
            )
        ):
            raise ValueError(
                'profile_smoothing must be "auto", a positive width in nm, '
                f"or None (mesh the polygon as given), got {profile_smoothing!r}"
            )
        workers = check_workers(workers)
        threads = check_blas_threads(blas_threads)

        started = time.perf_counter()
        sin_alpha = illumination.sin_alpha
        sin_gamma = illumination.sin_gamma
        # Exact zero in-plane (sin_gamma == 1), so the coupled system's
        # cross blocks vanish identically rather than to rounding.
        cos_gamma = float(np.sqrt(1.0 - sin_gamma**2))

        # The field oscillates along the boundary at the *reduced* wavelength.
        reduced_min = float(wavelengths.min()) / sin_gamma

        material: str | None = None
        indices = None
        metal_cone = 0
        metal_reduced_min = np.inf
        if conductivity == "tabulated":
            material, constants = _resolve_index(problem)
            indices = np.asarray(constants.n(wavelengths), dtype=np.complex128)
            if (indices.imag < 0.0).any():
                raise ValueError(
                    f"material {material!r} tabulates a negative absorption "
                    "(gain); the transmission problem loses uniqueness there "
                    "(Goray & Schmidt 2010, section 2.B)"
                )
            # The transverse wavelength inside the material is lambda /
            # |sqrt(n^2 - cos^2 gamma)|, short for visible-light metals.
            factors = np.abs(np.sqrt(indices**2 - cos_gamma**2))
            metal_reduced_min = float((wavelengths / factors).min())
            # The kernel's spectral sum must also clear the metal-side
            # propagating cone |k_t^-| d / 2 pi.
            metal_cone = int(
                np.ceil(
                    (factors / wavelengths).max() * problem.period
                )
            )

        # A measured polygon with no corner is meshed as a smooth curve under
        # the tabulated condition (profile_smoothing, docs/theory/integral.md
        # section 8), at a width set by the shortest transverse wavelength
        # either side unless one is given.
        smoothing = None
        if conductivity == "tabulated" and profile_smoothing is not None:
            smoothing = (
                min(reduced_min, metal_reduced_min) / _SMOOTHING_PER_WAVELENGTH
                if profile_smoothing == "auto"
                else float(profile_smoothing)
            )

        boundary = physical_boundary(
            problem, boundary_points, grading=corner_grading, smoothing=smoothing
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

        # A mesh that cannot resolve the reduced wavelength would not fail
        # loudly on its own. The coarsest spacing must resolve it, which on a graded mesh is
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

        if conductivity == "tabulated":
            # The same resolution floor on the metal side.
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
        # Transmitted orders per row, for a lossless material: {order: value}.
        transmitted_rows: list[dict[int, float]] = [{} for _ in wavelengths]
        propagating = np.zeros_like(efficiency, dtype=bool)
        warnings: list[str] = []
        anomalous: list[float] = []
        conditioning = ConditionLog()
        ill_conditioned: list[float] = []

        with capped_blas(threads):
            # One kernel assembly for the whole scan: the geometry-only half
            # (including the cached Kummer tails) is what every wavelength
            # and both media share; the per-wavelength memo is cleared each
            # row. A worker thread gets a fork of it, sharing the geometry.
            assembler = OperatorAssembler(
                boundary, period=problem.period, terms=terms
            )
            forks = iter([assembler])

            def worker_assembler() -> OperatorAssembler:
                return next(forks, None) or assembler.fork()

            def solve_row(row: int, assembler: OperatorAssembler):
                assembler.clear()
                wavelength = wavelengths[row]
                sines = sin_beta(
                    all_orders, wavelength, problem.period, sin_alpha, sin_gamma
                )
                live = is_propagating(sines)
                if not live.any():
                    return None
                values = np.zeros(int(live.sum()))
                absorbed = 0.0
                transmitted: dict[int, float] = {}
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
                        absorbed += solution.absorption / len(polarizations)
                        for order, value in zip(
                            solution.transmitted_orders.tolist(),
                            solution.transmitted_efficiencies,
                        ):
                            transmitted[order] = transmitted.get(order, 0.0) + float(
                                value
                            ) / len(polarizations)
                    # The boundary integral is the power that crossed into the
                    # material. What a lossless one transmits into orders is
                    # recorded as transmission, so absorption keeps only the
                    # rest, which for a lossless material is the
                    # discretisation error.
                    absorbed -= sum(transmitted.values())

                return live, values, absorbed, transmitted, row_log

            rows = run_rows(
                len(wavelengths),
                solve_row,
                state=worker_assembler,
                workers=workers,
                progress=progress,
            )

        for row, (wavelength, solved) in enumerate(zip(wavelengths, rows)):
            if solved is None:
                continue
            live, values, absorbed, transmitted_rows[row], row_log = solved
            sines = sin_beta(
                all_orders, wavelength, problem.period, sin_alpha, sin_gamma
            )
            if np.abs(cos_beta(sines[live])).min() < _ANOMALY_COS:
                anomalous.append(float(wavelength))
            conditioning.merge(row_log)
            if row_log.worst.get("nystrom", 0.0) > _WARN_KAPPA:
                ill_conditioned.append(float(wavelength))
            efficiency[row, live] = values
            propagating[row, live] = True
            absorption[row] = absorbed

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
        if conductivity == "tabulated" and isinstance(problem.profile, FromProfileData):
            rounded = polyline_rounded_turn(
                problem.profile.t,
                problem.profile.y,
                _ROUNDED_WINDOW_NM / problem.period,
            )
            if rounded > _ROUNDED_TURN_DEGREES:
                warnings.append(
                    f"the measured groove turns {rounded:.0f} deg within "
                    f"{_ROUNDED_WINDOW_NM:g} nm through vertices each under "
                    f"{CORNER_TURN_DEGREES:g} deg: a corner rounded tighter "
                    "than this solver converges on quickly, and smoothing does "
                    "not help it. On such a valley R was still 7.5e-3 from the "
                    "converged C-method at 1024 points; check convergence in "
                    "boundary_points, or compare with the C-method "
                    "(docs/theory/integral.md section 8, issue #6)"
                )
        if conductivity == "perfect" and problem.coating is not None:
            warnings.append(
                f"coating {problem.coating!r} is not consulted: the perfectly "
                "conducting boundary condition has no material in it, and the "
                "sum over orders is relative to a perfect reflector "
                "(docs/conventions.md section 10 item 5)"
            )

        transmission = None
        transmission_total = np.zeros(len(wavelengths))
        if any(transmitted_rows):
            seen = sorted({m for row_orders in transmitted_rows for m in row_orders})
            t_orders = np.arange(seen[0], seen[-1] + 1)
            t_efficiency = np.zeros((len(wavelengths), len(t_orders)))
            t_propagating = np.zeros_like(t_efficiency, dtype=bool)
            for row, row_orders in enumerate(transmitted_rows):
                for order, value in row_orders.items():
                    t_efficiency[row, order - seen[0]] = value
                    t_propagating[row, order - seen[0]] = True
            transmission = Transmission(
                orders=t_orders, efficiency=t_efficiency, propagating=t_propagating
            )
            transmission_total = transmission.total

        # Report, never rescale: the theorem says 1 -- for the sum alone
        # (perfect) or the sum plus absorption (tabulated) -- and the
        # discretisation says how close it got. Corners push the solve off
        # at first order.
        totals = np.where(propagating, efficiency, 0.0).sum(axis=1)
        if conductivity == "tabulated":
            totals = totals + absorption + transmission_total
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
                + (
                    ". This is a measured groove meshed as given, whose "
                    "slight bends the mesh cannot resolve: try "
                    'profile_smoothing="auto" and compare with the C-method '
                    "(docs/theory/integral.md section 8, issue #6)"
                    if conductivity == "tabulated"
                    and boundary.smoothed is None
                    and smooths(problem)
                    else ""
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
        if isinstance(problem.profile, FromProfileData):
            notes["profile"] = _profile_note(
                problem, boundary, conductivity, profile_smoothing
            )

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
            transmission=transmission,
        )


#: The registered singleton.
integral = register(IntegralSolver())
