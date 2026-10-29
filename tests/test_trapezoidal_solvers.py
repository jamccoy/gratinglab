"""``profiles.Trapezoidal`` through every solver, held against each other.

The groove is OGRE-like: period 160 nm, 68 nm deep (depth fraction 0.425), a
crest 0.30 of the period wide, 72 deg sidewalls, CXRO gold (the vendored
``materials/data/Au.txt``), in an extreme off-plane mount at graze 1.5 deg,
1.5 to 4.5 nm. That is the shape the 2021 cross-sections show for the 160 nm
laminar gratings, not a measurement of any one of them (``findings.md``,
"What a cross-section shows that a scan cannot").

Every tolerance here sits about three times above a number measured on
2026-10-29 from a convergence ladder per solver: integral 400 / 800 / 1600
boundary points, C-method truncation 16 to 128, RCWA truncation 10 / 20 / 40
(staircase slices follow it). ``findings.md``, "Trapezoidal through the four
solvers", has the whole map across sidewall angle and corner radius; this
file keeps the cases that are cheap enough to solve on every push.
"""

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Lamellar, Trapezoidal
from gratinglab.solvers import UnsupportedConfiguration, get_solver

WAVELENGTHS = np.array([1.5, 3.0, 4.5])
GROOVE = dict(depth_fraction=0.425, duty_cycle=0.30)
OFFPLANE = Illumination.offplane(graze=1.5, azimuth=0.0, polarization="TE")


def gold(profile):
    return Problem(period=160.0, profile=profile, coating="Au")


def solve(name, problem, illumination=OFFPLANE, **knobs):
    if name in ("integral", "cmethod"):
        knobs.setdefault("conductivity", "tabulated" if problem.coating else "perfect")
    return get_solver(name).solve(problem, illumination, WAVELENGTHS, **knobs)


def gap(a, b):
    return float(np.abs(a.efficiency - b.efficiency).max())


def balance(scan):
    total = np.where(scan.propagating, scan.efficiency, 0.0).sum(axis=1)
    if scan.absorption is not None:
        total = total + scan.absorption
    return float(np.abs(total - 1.0).max())


class TestVerticalAndSharpIsLamellar:
    """``sidewall_angle=90, corner_radius=0`` is ``Lamellar`` translated by
    half a period, which moves every order's phase and no efficiency."""

    LAMELLAR = Lamellar(**GROOVE)
    TRAPEZOID = Trapezoidal(**GROOVE)

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_rcwa_to_rounding(self, polarization):
        """Lamellar is one exact slab; the trapezoid becomes a staircase of
        20 slabs that are all the same slab, so the two agree to rounding.
        Measured 2.0e-13 (TE) and 6.1e-13 (TM)."""
        illumination = Illumination.offplane(graze=1.5, azimuth=0.0, polarization=polarization)
        lamellar = solve("rcwa", gold(self.LAMELLAR), illumination, truncation=10)
        trapezoid = solve("rcwa", gold(self.TRAPEZOID), illumination, truncation=10)
        assert trapezoid.provenance.notes["slices"] == 20
        assert gap(lamellar, trapezoid) < 1e-11
        assert np.abs(lamellar.absorption - trapezoid.absorption).max() < 1e-11

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_integral_perfect_conductor_to_rounding(self, polarization):
        """The same polygon, the same graded mesh relative to its corners:
        measured 4.4e-16 (TE) and 2.5e-13 (TM) at 400 points."""
        illumination = Illumination.offplane(graze=1.5, azimuth=0.0, polarization=polarization)
        lamellar = solve("integral", Problem(period=160.0, profile=self.LAMELLAR),
                         illumination, boundary_points=400)
        trapezoid = solve("integral", Problem(period=160.0, profile=self.TRAPEZOID),
                          illumination, boundary_points=400)
        assert gap(lamellar, trapezoid) < 1e-11

    def test_integral_gold_to_its_mesh_error(self):
        """Gold is not to rounding: 6.0e-6 at 400 points and 1.1e-7 at 800,
        against Lamellar's own change from 400 to 800 of 9.0e-5. The two
        meshes meet the period boundary in different places, so they carry
        different discretisation error that falls with the mesh."""
        lamellar = solve("integral", gold(self.LAMELLAR), boundary_points=400)
        trapezoid = solve("integral", gold(self.TRAPEZOID), boundary_points=400)
        assert gap(lamellar, trapezoid) < 2e-5

    @pytest.mark.parametrize("radius", [0.0, 0.05])
    def test_the_cmethod_refuses_a_vertical_wall(self, radius):
        with pytest.raises(UnsupportedConfiguration, match="vertical sidewalls"):
            solve("cmethod", gold(Trapezoidal(**GROOVE, corner_radius=radius)), truncation=16)


@pytest.fixture(scope="module")
def sharp():
    """72 deg, sharp corners: about 20 s, so solved once for the module."""
    problem = gold(Trapezoidal(**GROOVE, sidewall_angle=72.0))
    return {
        "integral": solve("integral", problem, boundary_points=800),
        "cmethod32": solve("cmethod", problem, truncation=32),
        "cmethod64": solve("cmethod", problem, truncation=64),
        "rcwa10": solve("rcwa", problem, truncation=10),
        "rcwa20": solve("rcwa", problem, truncation=20),
        "scalar": solve("scalar", problem),
    }


class TestSharpCorners:
    """Measured against the integral solver at 1600 points, whose own change
    from 800 is 7.9e-6: integral 800 7.9e-6, C-method 32 / 64 / 128
    3.6e-4 / 5.0e-5 / 9.5e-6, RCWA 10 / 20 / 40 7.2e-3 / 1.3e-3 / 2.4e-4."""

    def test_all_find_the_same_orders(self, sharp):
        orders = sharp["integral"].orders
        for scan in sharp.values():
            assert np.array_equal(scan.orders, orders)
            assert np.array_equal(scan.propagating, sharp["integral"].propagating)

    def test_the_cmethod_closes_on_the_integral_solver(self, sharp):
        coarse = gap(sharp["cmethod32"], sharp["integral"])
        fine = gap(sharp["cmethod64"], sharp["integral"])
        assert fine < 1.5e-4  # measured 5.0e-5
        assert fine < coarse / 3  # measured 7.2x
        assert any("corners" in w for w in sharp["cmethod64"].provenance.warnings)

    def test_rcwa_closes_on_the_integral_solver(self, sharp):
        """Slowly, being a staircase: about 5.5x per doubling."""
        coarse = gap(sharp["rcwa10"], sharp["integral"])
        fine = gap(sharp["rcwa20"], sharp["integral"])
        assert fine < 4e-3  # measured 1.3e-3
        assert fine < coarse / 3

    def test_each_rigorous_method_conserves_energy(self, sharp):
        assert balance(sharp["integral"]) < 5e-6  # measured 9e-7
        assert balance(sharp["rcwa20"]) < 1e-10
        assert balance(sharp["cmethod64"]) < 3e-4  # measured 6.8e-5

    def test_scalar_theory_says_it_is_out_of_its_regime(self, sharp):
        assert any("reduced ratio" in w for w in sharp["scalar"].provenance.warnings)


@pytest.fixture(scope="module")
def rounded():
    """72 deg with 1.6 nm fillets (0.01 of the period)."""
    problem = gold(Trapezoidal(**GROOVE, sidewall_angle=72.0, corner_radius=0.01))
    return {
        "cmethod64": solve("cmethod", problem, truncation=64),
        "cmethod128": solve("cmethod", problem, truncation=128),
        "rcwa10": solve("rcwa", problem, truncation=10),
        "rcwa20": solve("rcwa", problem, truncation=20),
        "integral": solve("integral", problem, boundary_points=400),
    }


class TestRoundedCorners:
    """The C-method is the reference here, and the integral solver is not.

    With the corners rounded the slope is continuous, the C-method converges
    fast (64 to 128 moves it 1.6e-5, balance 3.2e-5 at 64) and RCWA closes on
    it (7.9e-3, 1.9e-3, 4.7e-4 at 10, 20, 40). The integral solver's
    absorption fails on the fillets, which reach it as polylines of
    one-degree vertices: R + A - 1 reads 0.24 at 400 points and 0.04 at 800.
    That is the fault recorded for measured grooves (issue #6), and this
    class pins that it is reported rather than hidden; once it is fixed the
    last test fails and wants turning into an agreement test.
    """

    def test_the_cmethod_is_converged(self, rounded):
        assert gap(rounded["cmethod64"], rounded["cmethod128"]) < 5e-5  # measured 1.6e-5
        assert balance(rounded["cmethod64"]) < 1e-4  # measured 3.2e-5
        assert not any("corners" in w for w in rounded["cmethod64"].provenance.warnings)

    def test_rcwa_closes_on_the_cmethod(self, rounded):
        coarse = gap(rounded["rcwa10"], rounded["cmethod128"])
        fine = gap(rounded["rcwa20"], rounded["cmethod128"])
        assert fine < 6e-3  # measured 1.9e-3
        assert fine < coarse / 3  # measured 4.2x

    def test_the_integral_solver_reports_its_absorption_failing(self, rounded):
        scan = rounded["integral"]
        assert balance(scan) > 0.05  # measured 0.24
        assert any("strays from unity" in w for w in scan.provenance.warnings)

    def test_a_perfect_conductor_has_no_such_failure(self):
        """No absorption to get wrong: integral 400 and C-method 64 agree to
        2.8e-5 (the C-method's own 64-to-128 change is 2.4e-5)."""
        problem = Problem(
            period=160.0,
            profile=Trapezoidal(**GROOVE, sidewall_angle=72.0, corner_radius=0.01),
        )
        integral = solve("integral", problem, boundary_points=400)
        cmethod = solve("cmethod", problem, truncation=64)
        assert gap(integral, cmethod) < 1e-4
        assert balance(integral) < 1e-8  # measured 1.3e-10


class TestTheCmethodOnSteepWalls:
    """Where the C-method stops being usable, measured rather than assumed.

    Sharp corners, gold, against integral 800. Truncation 128 is 1.7e-4 off
    at 85 deg, 6.8e-4 at 87 deg and 6.8e-2 at 89 deg, where 64 to 128 still
    moves it by 0.15. RCWA is indifferent to the angle: 1.2e-4 at 89 deg and
    truncation 40.

    Until #7 the energy balance was the only thing that told you: truncation
    32 answered 0.34 from RCWA with R + A off by 1.14. The C-method now reads
    the wall's slope before solving and refuses a truncation below the
    harmonics it needs (``solvers.cmethod.wall_floor``).
    """

    def test_at_89_deg_it_refuses_too_few_harmonics(self):
        problem = gold(Trapezoidal(**GROOVE, sidewall_angle=89.0))
        with pytest.raises(ValueError, match="too few.*89 deg"):
            solve("cmethod", problem, truncation=32)
        rcwa = solve("rcwa", problem, truncation=10)
        assert balance(rcwa) < 1e-10


class TestTheLargestFillet:
    def test_its_boundary_has_no_phantom_corner(self):
        """At ``max_corner_radius`` the crest's straight run has zero width.
        Its end point used to be appended anyway, a segment of ~1e-17 in an
        arbitrary direction that the graded mesh took for a corner, putting
        two nodes on one point; the integral solver then refused every mesh."""
        base = Trapezoidal(**GROOVE, sidewall_angle=72.0)
        profile = Trapezoidal(**GROOVE, sidewall_angle=72.0, corner_radius=base.max_corner_radius)
        curve = profile.boundary(400, grading=4)
        assert not curve.graded
        assert np.hypot(np.diff(curve.t), np.diff(curve.y)).min() > 1e-3
