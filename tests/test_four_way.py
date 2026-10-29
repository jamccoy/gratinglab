"""All four solvers on one problem, through ``compare.sweep``.

The pairwise files (``test_rcwa.py``, ``test_cmethod.py``,
``test_cross_method.py``) each meet the integral solver on a problem chosen to
suit the method under test. This file does the thing the project exists for:
one ``Problem``, one ``Illumination``, every registered method, one
``Comparison``. The three rigorous methods are held to each other order by
order; scalar theory is held only to what it promises.

The case is a shallow gold sinusoid in grazing incidence, the founding
application's regime: period 40 nm, depth fraction 0.05, alpha = 80 deg, TE,
2 to 6 nm, absolute efficiencies with absorption from the vendored CXRO table
(``materials/data/Au.txt``, in the repository and therefore in CI, never in
the wheel; see ``SOURCES.md``). Forty orders propagate at 2 nm and fourteen at
6 nm.

How the tolerances were set
---------------------------
From a convergence ladder, not from one run. Measured 2026-10-28 on MC1 as
the largest difference over the five wavelengths against the C-method at
truncation 80, which is converged: truncation 24, 40 and 80 agree to 2e-15,
and the integral solver closes on it as its mesh is refined.

==================  ==========  ==========  ==========  ==========
knob                order 0     order 1     orders 2+   absorbed
==================  ==========  ==========  ==========  ==========
integral 600        1.05e-6     4.6e-7      1.6e-7      1.9e-6
integral 800        4.8e-7      2.2e-7      7.7e-8      9.4e-7
integral 1600       8.1e-8      4.2e-8      1.5e-8      1.7e-7
rcwa 45             4.1e-6      4.1e-7      2.2e-7      3.6e-6
rcwa 60             2.6e-6      2.7e-7      1.4e-7      2.3e-6
rcwa 120            9.1e-7      9.7e-8      5.1e-8      7.8e-7
==================  ==========  ==========  ==========  ==========

Both fall as their own knob is refined, so the differences are discretisation
error, not a disagreement about the physics. The test runs integral 800,
RCWA 60 and C-method 40, which is the 2026-10-28 run, and its tolerances sit
about three times above those rows. RCWA's zero order is its weakest number
here because a sinusoid is a staircase to it.

The perfect-conductor variants need no table. RCWA cannot take part: it
expands a permittivity, and a perfect conductor has none, so it refuses, and
that refusal is asserted so that the day it learns the case, this file hears.
"""

import numpy as np
import pytest

from gratinglab.checks import check_energy_balance
from gratinglab.compare import align, sweep
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Sinusoidal
from gratinglab.solvers import UnsupportedConfiguration, get_solver

WAVELENGTHS = np.linspace(2.0, 6.0, 5)
METHODS = ["scalar", "integral", "rcwa", "cmethod"]
RIGOROUS = ("integral", "rcwa", "cmethod")
SINE = Sinusoidal(depth_fraction=0.05)


def worst(comparison, a, b):
    """Largest ``|a - b|`` in order 0, order 1 and the rest, over the grid."""
    delta = np.abs(comparison.difference(a, b))
    zero = int(np.flatnonzero(comparison.orders == 0)[0])
    first = int(np.flatnonzero(comparison.orders == 1)[0])
    rest = np.delete(delta, [zero, first], axis=1)
    return delta[:, zero].max(), delta[:, first].max(), rest.max()


@pytest.fixture(scope="module")
def gold():
    """The 2026-10-28 run: about 30 s, so it is solved once for the module."""
    problem = Problem(period=40.0, profile=SINE, coating="Au")
    illumination = Illumination.classical(alpha=80.0, polarization="TE")
    scans = sweep(
        problem,
        illumination,
        WAVELENGTHS,
        METHODS,
        {
            "integral": {"conductivity": "tabulated", "boundary_points": 800},
            "cmethod": {"conductivity": "tabulated", "truncation": 40},
            "rcwa": {"truncation": 60},
        },
    )
    return problem, illumination, dict(zip(METHODS, scans)), align(scans)


class TestGoldInGrazingIncidence:
    def test_every_method_finds_the_same_propagating_orders(self, gold):
        _, _, scans, comparison = gold
        assert comparison.methods == tuple(METHODS)
        reference = scans["integral"]
        for name in METHODS:
            assert np.array_equal(scans[name].orders, reference.orders), name
            assert np.array_equal(scans[name].propagating, reference.propagating), name
        assert reference.propagating.sum(axis=1).tolist() == [40, 27, 20, 16, 14]

    def test_the_cmethod_is_converged_at_this_truncation(self, gold):
        """What makes the C-method the anchor below: its own change from
        truncation 24 to 40 is far under any tolerance in this file, so the
        C-method-to-integral gap is the integral solver's error."""
        problem, illumination, scans, _ = gold
        (coarse,) = sweep(
            problem, illumination, WAVELENGTHS, ["cmethod"],
            {"cmethod": {"conductivity": "tabulated", "truncation": 24}},
        )
        assert np.abs(coarse.efficiency - scans["cmethod"].efficiency).max() < 1e-9

    def test_the_cmethod_matches_the_integral_solver_order_by_order(self, gold):
        _, _, scans, comparison = gold
        zero, first, rest = worst(comparison, "cmethod", "integral")
        assert zero < 1.5e-6  # measured 4.8e-7
        assert first < 7e-7  # measured 2.2e-7
        assert rest < 3e-7  # measured 7.7e-8
        absorbed = np.abs(scans["cmethod"].absorption - scans["integral"].absorption)
        assert absorbed.max() < 3e-6  # measured 9.4e-7

    def test_rcwa_matches_the_integral_solver_order_by_order(self, gold):
        _, _, scans, comparison = gold
        zero, first, rest = worst(comparison, "rcwa", "integral")
        assert zero < 1e-5  # measured 3.1e-6
        assert first < 1e-6  # measured 2.6e-7
        assert rest < 5e-7  # measured 1.1e-7
        absorbed = np.abs(scans["rcwa"].absorption - scans["integral"].absorption)
        assert absorbed.max() < 1e-5  # measured 2.4e-6

    def test_each_rigorous_method_conserves_energy(self, gold):
        """``R + A = 1`` for a one-interface absorbing grating, with ``A``
        computed independently of ``R`` by all three. The integral solver's
        deviation is its mesh error, measured 1.7e-6 at 800 points; the
        modal methods balance to rounding."""
        _, _, scans, _ = gold
        assert check_energy_balance(scans["integral"], lossless=True, tolerance=5e-6).passed
        assert check_energy_balance(scans["rcwa"], lossless=True, tolerance=1e-10).passed
        assert check_energy_balance(scans["cmethod"], lossless=True, tolerance=1e-10).passed
        for name in RIGOROUS:
            assert scans[name].provenance.notes["normalization"] == "absolute", name

    def test_scalar_theory_is_an_approximation_here_and_says_so(self, gold):
        """Out of its regime: the reduced ratio reaches 0.15 and the facet
        graze exceeds gold's critical angle. It still gets what the grating
        equation decides, the orders, and lands within a factor of two of the
        rigorous first order (measured 0.60 to 0.86 of it), and it warns that
        per-order numbers want a rigorous cross-check. Nothing tighter is
        asserted, because nothing tighter is promised."""
        _, _, scans, comparison = gold
        scalar = scans["scalar"]
        assert check_energy_balance(scalar).passed  # passive: sum <= 1
        ratio = comparison.order(1)["scalar"] / comparison.order(1)["cmethod"]
        assert np.all((ratio > 0.5) & (ratio < 2.0)), ratio
        assert any("reduced ratio" in w for w in scalar.provenance.warnings)


class TestPerfectConductor:
    """The same groove and mount with no material, so no table and no
    absorption: efficiencies relative to a perfect reflector, summing to one.

    Measured 2026-10-28 against the C-method at truncation 80: TE with 400
    boundary points 2.0e-7, falling to 1.1e-8 at 1600; TM with 800 points
    7.1e-6, falling to 1.5e-6 at 1600. TM converges more slowly in the
    integral solver, hence the extra points and the wider tolerance.
    """

    PROBLEM = Problem(period=40.0, profile=SINE)

    @pytest.mark.parametrize(
        "polarization, boundary_points, tolerance, balance",
        [("TE", 400, 1e-6, 1e-9), ("TM", 800, 3e-5, 4e-5)],
    )
    def test_three_methods_agree(self, polarization, boundary_points, tolerance, balance):
        illumination = Illumination.classical(alpha=80.0, polarization=polarization)
        scans = sweep(
            self.PROBLEM,
            illumination,
            WAVELENGTHS,
            ["scalar", "integral", "cmethod"],
            {"integral": {"boundary_points": boundary_points}, "cmethod": {"truncation": 40}},
        )
        _, integral, cmethod = scans
        comparison = align(scans)

        for scan in scans:
            assert np.array_equal(scan.propagating, integral.propagating)
            assert scan.provenance.notes["normalization"].startswith("relative")
        assert comparison.max_abs_difference("cmethod", "integral") < tolerance
        # Measured deviation 4.5e-11 (TE) and 1.2e-5 (TM) for the integral solver.
        assert check_energy_balance(integral, lossless=True, tolerance=balance).passed
        assert check_energy_balance(cmethod, lossless=True, tolerance=1e-10).passed

    def test_scalar_cannot_see_the_polarization_the_rigorous_methods_do(self):
        """Scalar theory is polarization-blind by construction, so its TE and
        TM scans are identical, while the rigorous first order differs by 45%
        to 216% between them here (0.281 against 0.407 at 2 nm, 0.086 against
        0.273 at 6 nm). Which is why no agreement is asserted for it beyond
        the orders."""
        te, tm = (
            sweep(
                self.PROBLEM,
                Illumination.classical(alpha=80.0, polarization=p),
                WAVELENGTHS,
                ["scalar", "cmethod"],
                {"cmethod": {"truncation": 40}},
            )
            for p in ("TE", "TM")
        )
        assert np.array_equal(te[0].efficiency, tm[0].efficiency)
        first = int(np.flatnonzero(te[1].orders == 1)[0])
        rigorous_te, rigorous_tm = te[1].efficiency[:, first], tm[1].efficiency[:, first]
        assert np.all(np.abs(rigorous_tm - rigorous_te) / rigorous_te > 0.2)

    def test_rcwa_refuses_a_perfect_conductor(self):
        with pytest.raises(UnsupportedConfiguration, match="perfect conductor"):
            get_solver("rcwa").solve(
                self.PROBLEM, Illumination.classical(alpha=80.0), WAVELENGTHS
            )
