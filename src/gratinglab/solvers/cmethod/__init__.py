r"""The C-method (Chandezon): rigorous efficiency in curvilinear coordinates.

The translation coordinate ``u = y - a(x)`` maps the groove onto the plane
``u = 0``, so the boundary conditions become a matching of Fourier series
across a flat interface, and in each homogeneous medium the field is a sum
of eigenmodes of a matrix built from the Fourier coefficients of the slope
``a'(x)``. ``docs/theory/cmethod.md`` is the method note: the transform,
the eigenproblem and the factorization it uses, what is refused, and what
has been validated.

Two boundary conditions, named as the integral solver names them, so one
``Problem`` and one option reach either method:

- ``conductivity="perfect"`` (the default): a perfect conductor in any mount.
  The conical problem reduces exactly to the in-plane one at the reduced
  wavelength ``lambda / sin(gamma)``; TE is ``E_z = 0`` and TM is
  ``dH_z/dn = 0`` on the boundary. Efficiencies are relative to a perfect
  reflector, so the sum is 1 in the continuum.
- ``conductivity="tabulated"``: one interface into a material read per
  wavelength from ``coating`` (falling back to ``substrate``). In a conical
  mount ``E_z`` and ``H_z`` couple through ``k_z`` at the boundary.
  Efficiencies are absolute, and the flux into the material is computed from
  the fields on the boundary and recorded as the scan's ``absorption``, so
  ``R + A = 1`` is a check. For a lossless material that flux is the
  transmitted total; the transmitted orders themselves are solved but the
  scan carries reflected orders only.

The accuracy knob is ``truncation``: ``N`` retained harmonics on each side of
the centre of the propagating orders, ``2N + 1`` in all.

What the method cannot describe it refuses, with
:class:`~gratinglab.solvers.base.UnsupportedConfiguration`: an undercut, a
vertical facet or sidewall (``Lamellar``, ``Blazed(antiblaze_angle=90)``, a
measured profile with two heights at one position). A profile with corners (a sawtooth) is solved, and a
provenance warning says that it converges slowly and that its energy balance
overstates its accuracy there (measured; see the method note, section 6).
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
from ._core import Harmonics, solve_finite, solve_perfect
from ._series import profile_series

__all__ = ["CMethodSolver", "cmethod"]

#: Largest ``truncation`` accepted. The eigenproblem is ``2(2N+1)`` square
#: per medium per wavelength, dense, so its cost grows as ``N^3``; at 256 a
#: single wavelength takes seconds, and a ladder stepping past it would take
#: an afternoon. A sweep that has not plateaued by here is the answer.
_MAX_TRUNCATION = 256

#: Deviation of the conserved quantity from unity worth naming in provenance.
#: Reported, never rescaled (the integral solver uses the same figure).
_ENERGY_DEVIATION_WARN = 0.005

#: ``min |cos beta_m|`` below which an order is passing off (Rayleigh anomaly).
_ANOMALY_COS = 1e-3

#: Estimated condition number of the boundary-condition system above which
#: the solve is losing digits worth naming. Measured on sinusoids from 0.3 to
#: 1.0 of the period deep: below 1e12 every case agreed with the integral
#: solver to that solver's own error, and from 1e14 cases on profiles 0.8 P
#: deep or more were wrong by 1e-4 to 1. The estimate is pessimistic on
#: shallow profiles (0.3 P stays correct to 6 digits at 1e17), so it is a
#: warning; the energy balance, which tracks this failure closely on a smooth
#: profile, and the convergence harness are what expose a bad rung.
_CONDITION_WARN = 1e12

#: ``|cos beta_m|`` below which an order is taken to sit exactly on its
#: Rayleigh anomaly. There its upgoing and downgoing plane waves coincide, the
#: plane-wave basis loses a member and the system is singular, so the
#: wavelength is refused by name rather than solved in a degenerate basis.
_GRAZING_COS = 1e-6

#: Slope discontinuity, in degrees of turn, above which a polygon vertex is a
#: corner. The same figure as ``profiles.CORNER_TURN_DEGREES``.
_CORNER_DEGREES = 10.0


def _resolve_index(problem: Problem):
    from ...materials import lookup

    name = problem.coating or problem.substrate
    if name is None:
        raise UnsupportedConfiguration(
            'conductivity="tabulated" needs a material: name one on '
            "Problem.coating (or substrate). The boundary condition reads its "
            "complex index per wavelength."
        )
    return name, lookup(name)


def _has_corners(series) -> bool:
    if series.kind != "polygon":
        return False
    s0, width, a0, slope = series._segments()
    turn = np.degrees(np.abs(np.arctan(slope) - np.arctan(np.roll(slope, 1))))
    return bool((turn > _CORNER_DEGREES).any())


class CMethodSolver:
    """The C-method for single-valued profiles with a finite slope."""

    capabilities = Capabilities(
        name="cmethod",
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
        truncation: int = 32,
        conductivity: Literal["perfect", "tabulated"] = "perfect",
        progress: "Progress | None" = None,
    ) -> EfficiencyScan:
        """Efficiency over a wavelength scan at fixed geometry.

        Parameters
        ----------
        truncation
            ``N``: harmonics ``c - N .. c + N`` are retained, ``c`` the middle
            of the propagating orders at the shortest wavelength. Must leave
            at least one evanescent harmonic beyond each end of them. This is
            the knob the convergence harness sweeps.
        conductivity
            ``"perfect"`` or ``"tabulated"``, as for the integral solver.
            Anything else raises :class:`UnsupportedConfiguration`.
        progress
            Per-wavelength callback under the contract in
            :class:`~gratinglab.solvers.base.Solver`.
        """
        self.capabilities.check(problem, illumination)
        if conductivity not in ("perfect", "tabulated"):
            raise UnsupportedConfiguration(
                f"conductivity={conductivity!r} is not implemented: "
                '"perfect" is the perfectly conducting boundary and '
                '"tabulated" one interface into an absorbing material. A '
                "layer stack or coating is not modelled yet."
            )
        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
        if wavelengths.ndim != 1:
            raise ValueError("wavelengths must be one-dimensional")
        if (wavelengths <= 0).any():
            raise ValueError("wavelengths must be positive")
        if (
            isinstance(truncation, bool)
            or not isinstance(truncation, (int, np.integer))
            or not 1 <= truncation <= _MAX_TRUNCATION
        ):
            raise ValueError(
                f"truncation must be an integer from 1 to {_MAX_TRUNCATION}, "
                f"got {truncation!r}"
            )
        truncation = int(truncation)

        started = time.perf_counter()
        sin_alpha = illumination.sin_alpha
        sin_gamma = illumination.sin_gamma
        cos_gamma = float(np.sqrt(1.0 - sin_gamma**2))

        series = profile_series(problem)
        if series.kind == "sampled":
            # Enough samples for the harmonics and for exp(i beta a), whose
            # local frequency is beta P |a'| / 2 pi cycles per period.
            bandwidth = problem.period * sin_gamma / float(wavelengths.min())
            need = int(8 * (4 * truncation + 1 + bandwidth * (1.0 + series.max_slope)))
            if need > len(series.s):
                series = profile_series(problem, samples=1 << (need - 1).bit_length())

        all_orders = order_range(
            float(wavelengths.min()), problem.period, sin_alpha, sin_gamma
        )
        lo, hi = int(all_orders.min()), int(all_orders.max())
        center = (lo + hi) // 2
        needed = max(hi - center, center - lo) + 1
        if truncation < needed:
            raise ValueError(
                f"truncation={truncation} cannot represent orders {lo} to {hi} "
                f"with an evanescent harmonic beyond each end; need at least "
                f"{needed}"
            )
        harmonics = Harmonics.build(series, center, truncation)

        material: str | None = None
        indices = None
        if conductivity == "tabulated":
            material, constants = _resolve_index(problem)
            indices = np.asarray(constants.n(wavelengths), dtype=np.complex128)
            if (indices.imag < 0.0).any():
                raise ValueError(
                    f"material {material!r} tabulates a negative absorption (gain)"
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
        anomalous: list[float] = []
        condition = 1.0

        for row, wavelength in enumerate(wavelengths):
            if progress is not None:
                progress(row, len(wavelengths))
            sines = sin_beta(all_orders, wavelength, problem.period, sin_alpha, sin_gamma)
            live = is_propagating(sines)
            grazing = np.abs(cos_beta(sines[live])).min()
            if grazing < _GRAZING_COS:
                raise ValueError(
                    f"at {wavelength:.10g} nm order "
                    f"{int(all_orders[live][np.argmin(np.abs(cos_beta(sines[live])))])} "
                    "sits exactly on its Rayleigh anomaly (grazing), where the "
                    "C-method's plane-wave basis is degenerate; move the "
                    "wavelength off it"
                )
            if grazing < _ANOMALY_COS:
                anomalous.append(float(wavelength))

            values = np.zeros(int(live.sum()))
            if conductivity == "perfect":
                solutions = [
                    solve_perfect(
                        harmonics,
                        wavelength=float(wavelength),
                        period=problem.period,
                        sin_alpha=sin_alpha,
                        sin_gamma=sin_gamma,
                        polarization=polarization,
                    )
                    for polarization in polarizations
                ]
            else:
                solutions = solve_finite(
                    harmonics,
                    wavelength=float(wavelength),
                    period=problem.period,
                    sin_alpha=sin_alpha,
                    sin_gamma=sin_gamma,
                    cos_gamma=cos_gamma,
                    index=complex(indices[row]),
                    incidents=[incidents[p] for p in polarizations],
                )
            for solution in solutions:
                columns = np.searchsorted(all_orders[live], solution.orders)
                values[columns] += solution.efficiencies / len(polarizations)
                absorption[row] += solution.absorption / len(polarizations)
                condition = max(condition, solution.condition)

            efficiency[row, live] = values
            propagating[row, live] = True

        if progress is not None:
            progress(len(wavelengths), len(wavelengths))

        warnings: list[str] = []
        if anomalous:
            span = f"{min(anomalous):.6g}-{max(anomalous):.6g} nm"
            warnings.append(
                f"an order passes off within |cos beta| < {_ANOMALY_COS:g} at "
                f"{len(anomalous)} wavelength(s) in {span} (Rayleigh anomaly); "
                "its plane wave's flux carries 1/beta_m and those rows are "
                "numerically delicate"
            )
        if condition > _CONDITION_WARN:
            warnings.append(
                f"the boundary-condition system's condition number reaches "
                f"{condition:.1e}: digits are being lost. On a deep profile the "
                "C-method's evanescent modes span exp(|Im rho| depth), and past "
                "this point a higher truncation can make the answer worse "
                "rather than better. Read the energy balance, sweep truncation, "
                "and prefer a lower one that agrees"
            )
        corners = _has_corners(series)
        if corners:
            warnings.append(
                "the profile has corners (slope discontinuities). The C-method "
                "converges slowly on them, roughly first order in truncation, "
                "and its energy balance closes long before its efficiencies "
                "do, so the balance is not an error estimate here: sweep "
                "truncation, or compare with the integral solver"
            )
        if conductivity == "perfect" and problem.coating is not None:
            warnings.append(
                f"coating {problem.coating!r} is not consulted: the perfectly "
                "conducting boundary condition has no material in it, and the "
                "sum over orders is relative to a perfect reflector"
            )

        totals = np.where(propagating, efficiency, 0.0).sum(axis=1)
        if conductivity == "tabulated":
            totals = totals + absorption
        computed = totals[propagating.any(axis=1)]
        deviation = float(np.abs(computed - 1.0).max()) if computed.size else 0.0
        if deviation > _ENERGY_DEVIATION_WARN:
            warnings.append(
                f"summed efficiency"
                f"{' plus absorption' if conductivity == 'tabulated' else ''} "
                f"strays from unity by up to {100 * deviation:.2f}% "
                "(theorem: exactly 1); raise truncation"
            )

        from ... import __version__

        notes: dict = {
            "harmonics": f"{center - truncation} to {center + truncation}",
            "energy_balance_deviation": deviation,
            "condition": condition,
            "profile_series": (
                "exact polygon integrals" if series.kind == "polygon"
                else f"FFT of {len(series.s)} samples"
            ),
            "factorization": "Laurent throughout, on continuous unknowns (Li)",
        }
        if conductivity == "perfect":
            notes["boundary_condition"] = "perfectly-conducting"
            notes["normalization"] = "relative to a perfect reflector"
        else:
            notes["boundary_condition"] = "finite-conductivity, single interface"
            notes["normalization"] = "absolute"
            notes["material"] = material
            if (indices.imag == 0.0).any():
                notes["absorption"] = (
                    "flux into the material; where it is lossless this is the "
                    "transmitted total"
                )

        return EfficiencyScan(
            wavelengths=wavelengths,
            orders=all_orders,
            efficiency=efficiency,
            propagating=propagating,
            provenance=Provenance(
                method="cmethod",
                version=__version__,
                truncation=truncation,
                converged=None,
                wall_time_s=time.perf_counter() - started,
                warnings=tuple(warnings),
                notes=notes,
            ),
            absorption=absorption if conductivity == "tabulated" else None,
        )


#: The registered singleton.
cmethod = register(CMethodSolver())
