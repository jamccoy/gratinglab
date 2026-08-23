"""The in-plane perfectly-conducting core, against every closed-form answer.

The decisive references, in order of strictness: a flat mirror (exact), the
energy-balance theorem (exact in the continuum), Maystre & Popov's published
Table 4.1 (four digits, at their own P and M), and the Marechal-Stroke
theorem for the echelette (exact, but corner convergence is first-order --
the tolerance says so honestly).
"""

import numpy as np
import pytest

from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, Lamellar, Sinusoidal
from gratinglab.solvers.integral._boundary import PhysicalBoundary, physical_boundary
from gratinglab.solvers.integral._core import solve_transverse


def flat_mirror(period=600.0, points=64):
    x = np.linspace(0.0, period, points, endpoint=False)
    return PhysicalBoundary(
        x=x,
        y=np.zeros(points),
        nx=np.zeros(points),
        ny=np.ones(points),
        arc_length=period,
    )


class TestFlatMirror:
    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_unit_specular(self, polarization):
        solution = solve_transverse(
            flat_mirror(),
            wavelength=550.0,
            period=600.0,
            sin_alpha=0.5,
            polarization=polarization,
            terms=60,
        )
        specular = solution.efficiencies[solution.orders == 0][0]
        assert specular == pytest.approx(1.0, abs=1e-6)
        assert solution.total == pytest.approx(1.0, abs=1e-6)

    def test_amplitude_signs_pin_the_boundary_conditions(self):
        """TE (Dirichlet): the total E field vanishes on the mirror, so
        r0 = -1. TM (Neumann): the total H field doubles, so r0 = +1.
        Wrong global signs cancel in efficiencies; this is the test that
        would catch them."""
        te = solve_transverse(
            flat_mirror(), wavelength=550.0, period=600.0, sin_alpha=0.5,
            polarization="TE", terms=60,
        )
        tm = solve_transverse(
            flat_mirror(), wavelength=550.0, period=600.0, sin_alpha=0.5,
            polarization="TM", terms=60,
        )
        assert te.amplitudes[te.orders == 0][0] == pytest.approx(-1.0, abs=1e-3)
        assert tm.amplitudes[tm.orders == 0][0] == pytest.approx(1.0, abs=1e-3)


class TestMaystrePopovTable41:
    """Sinusoid, period 600 nm, depth 180 nm, alpha = 30 deg, lambda = 600 nm
    (Littrow). M&P's -1 order is this project's m = +1
    (sin beta_m = m lambda / period - sin alpha). Values at their own
    settings P = 110, M = 50."""

    BOUNDARY = physical_boundary(
        Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3)), 110
    )

    @pytest.mark.parametrize(
        "polarization, littrow, specular",
        [("TE", 0.4659, 0.5341), ("TM", 0.9579, 0.0421)],
    )
    def test_efficiencies(self, polarization, littrow, specular):
        solution = solve_transverse(
            self.BOUNDARY,
            wavelength=600.0,
            period=600.0,
            sin_alpha=0.5,
            polarization=polarization,
            terms=50,
        )
        got = dict(zip(solution.orders.tolist(), solution.efficiencies))
        assert got[1] == pytest.approx(littrow, abs=1e-3)
        assert got[0] == pytest.approx(specular, abs=1e-3)
        assert solution.total == pytest.approx(1.0, abs=1e-3)


class TestMarechalStroke:
    """Echelette with blaze 30 deg, apex 90 deg, illuminated normal to the
    blaze facet (alpha = 30 deg, lambda = 600 nm): for TM the -1 order
    (our m = +1) carries exactly unit efficiency (M&P section 4.7.2).

    The corner makes TM convergence first-order (measured: 0.964 at P=100,
    0.982 at P=250, 0.992 at P=500), so the assertion is the trend and a
    floor, not four digits -- M&P needed a graded corner mesh for those."""

    PROBLEM = Problem(
        period=600.0, profile=Blazed(blaze_angle=30.0, antiblaze_angle=60.0)
    )

    def test_tm_blaze_order_dominates_and_improves_with_points(self):
        totals = []
        for points, terms in ((100, 50), (250, 100)):
            solution = solve_transverse(
                physical_boundary(self.PROBLEM, points),
                wavelength=600.0,
                period=600.0,
                sin_alpha=0.5,
                polarization="TM",
                terms=terms,
            )
            got = dict(zip(solution.orders.tolist(), solution.efficiencies))
            assert got[0] == pytest.approx(0.0, abs=1e-3)
            totals.append(got[1])
        assert totals[0] < totals[1]
        assert totals[1] > 0.98

    def test_te_energy_balance_is_untouched_by_the_corner(self):
        solution = solve_transverse(
            physical_boundary(self.PROBLEM, 250),
            wavelength=600.0,
            period=600.0,
            sin_alpha=0.5,
            polarization="TE",
            terms=100,
        )
        assert solution.total == pytest.approx(1.0, abs=1e-4)


class TestEnergyBalance:
    """Sum of efficiencies = 1 is a theorem for a perfect conductor
    (M&P eqs. 4.34/4.42; thesis eq:prop_order_unity). TE holds it to
    quadrature precision on every profile; TM holds it tightly only where
    the boundary is smooth -- the corner deficit is real, measured, and
    documented, so the tolerances differ."""

    SMOOTH = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.5))
    CORNERED = [
        Problem(period=600.0, profile=Blazed(blaze_angle=15.0)),
        Problem(period=600.0, profile=Lamellar(depth_fraction=0.2)),
    ]

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_smooth_profile_is_tight(self, polarization):
        solution = solve_transverse(
            physical_boundary(self.SMOOTH, 200),
            wavelength=500.0,
            period=600.0,
            sin_alpha=0.35,
            polarization=polarization,
            terms=90,
        )
        assert solution.total == pytest.approx(1.0, abs=1e-3)

    @pytest.mark.parametrize("problem", CORNERED, ids=["blazed-vertical", "lamellar"])
    def test_cornered_profiles(self, problem):
        te = solve_transverse(
            physical_boundary(problem, 200),
            wavelength=500.0, period=600.0, sin_alpha=0.35,
            polarization="TE", terms=90,
        )
        tm = solve_transverse(
            physical_boundary(problem, 200),
            wavelength=500.0, period=600.0, sin_alpha=0.35,
            polarization="TM", terms=90,
        )
        assert te.total == pytest.approx(1.0, abs=1e-4)
        assert tm.total == pytest.approx(1.0, abs=0.05)

    def test_xray_regime_reduced_wavelength(self):
        """The off-plane X-ray case after the conical reduction: TASTE-like
        geometry, lambda = 1 nm at gamma = 1.25 deg becomes a 45.8 nm
        transverse wavelength against a 315 nm period."""
        problem = Problem(
            period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
        )
        reduced = 1.0 / np.sin(np.radians(1.25))
        solution = solve_transverse(
            physical_boundary(problem, 200),
            wavelength=reduced,
            period=315.15,
            sin_alpha=float(np.sin(np.radians(19.99))),
            polarization="TE",
            terms=100,
        )
        assert len(solution.orders) > 10
        assert solution.total == pytest.approx(1.0, abs=1e-4)
