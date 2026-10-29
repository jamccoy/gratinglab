"""The C-method, against every answer that does not come from this code.

In order of strictness: the profile's Fourier series against a closed form,
a flat boundary against Fresnel (amplitude and phase), Goray & Schmidt's
Table 3 (which is itself Li's C-method), the Marechal-Stroke theorem, Maystre
& Popov's Table 4.1, the energy theorems, the perfect-conductor limit of the
coupled system, reciprocity, and agreement with the integral solver -- a
different method, sharing nothing below ``Problem``. Then the two failure
regions the method note measures (corners and deep profiles), and the
refusals. ``docs/theory/cmethod.md`` section 6 records the numbers.
"""

import numpy as np
import pytest
from scipy.special import jv

from gratinglab.checks import check_energy_balance, check_reciprocity
from gratinglab.convergence import check_convergence
from gratinglab.illumination import Illumination
from gratinglab.materials import fresnel
from gratinglab.materials.optical import OpticalConstants, write_ari
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, FromProfileData, Lamellar, Sinusoidal, Trapezoidal
from gratinglab.solvers import UnsupportedConfiguration, cmethod, get_solver, integral
from gratinglab.solvers.cmethod import _MAX_TRUNCATION, wall_floor
from gratinglab.solvers.cmethod._core import Harmonics, solve_finite, solve_perfect
from gratinglab.solvers.cmethod._series import ProfileSeries, profile_series
from gratinglab.solvers.integral._boundary import physical_boundary
from gratinglab.solvers.integral._finite import solve_finite_states
from gratinglab.solvers.integral._nystrom import OperatorAssembler

GOLD_VISIBLE = 0.2 + 3.0j
XRAY = 1.0 - 1e-3 + 1e-4j
GLASS = 1.5 + 0.0j

FLAT = Problem(period=600.0, profile=FromProfileData(t=(0.0, 0.5, 1.0), y=(0.0, 0.0, 0.0)))
SINE = Problem(period=500.0, profile=Sinusoidal(depth_fraction=0.3))


def harmonics(problem, truncation, center=0):
    return Harmonics.build(profile_series(problem), center, truncation)


def by_order(orders, values):
    return dict(zip(np.asarray(orders).tolist(), np.asarray(values)))


def integral_finite(problem, *, wavelength, sin_alpha, index, cos_gamma, incident, points):
    boundary = physical_boundary(problem, points)
    terms = points // 2
    solution = solve_finite_states(
        boundary,
        wavelength=wavelength,
        period=problem.period,
        sin_alpha=sin_alpha,
        index=index,
        cos_gamma=cos_gamma,
        incidents=[incident],
        terms=terms,
        assembler=OperatorAssembler(boundary, period=problem.period, terms=terms),
    )[0]
    return solution


class TestRegistration:
    def test_available_by_name(self):
        assert get_solver("cmethod") is cmethod

    def test_capabilities(self):
        caps = cmethod.capabilities
        assert caps.rigorous and caps.conical and caps.reports_progress
        assert caps.accuracy_knob == "truncation"
        assert not caps.handles_undercut


class TestProfileSeries:
    """The two inputs every mode is built from, against closed forms."""

    def test_sinusoid_plane_wave_is_a_bessel_series(self):
        """``exp(i beta a)`` with ``a = (D/2)(1 - cos 2 pi s)`` has
        coefficients ``exp(i beta D/2) (-i)^j J_j(beta D/2)``, from the
        Jacobi-Anger expansion -- an answer the FFT path cannot share."""
        series = profile_series(SINE)
        depth = SINE.depth
        beta = 0.013
        j = np.arange(-12, 13)
        expected = np.exp(0.5j * beta * depth) * (-1j) ** j * jv(j, 0.5 * beta * depth)
        assert np.allclose(series.plane_wave(beta, j), expected, atol=1e-14)

    def test_sawtooth_slope_coefficients_are_exact(self):
        """The polygon path integrates each facet in closed form. The slope
        of a sawtooth is a two-level step, whose coefficients are known."""
        problem = Problem(period=1000.0, profile=Blazed(blaze_angle=20.0, antiblaze_angle=50.0))
        series = profile_series(problem)
        s0, width, _, slope = series._segments()
        j = np.arange(1, 9)
        expected = sum(
            sigma * (np.exp(-2j * np.pi * j * a) - np.exp(-2j * np.pi * j * (a + w)))
            / (2j * np.pi * j)
            for a, w, sigma in zip(s0, width, slope)
        )
        assert np.allclose(series.slope(j), expected, atol=1e-15)
        assert series.slope(np.array([0]))[0] == pytest.approx(0.0, abs=1e-15)
        # Mirrored onto x = (1 - t) P: the blaze facet descends with x.
        assert set(np.round(slope, 12)) == {
            round(-np.tan(np.radians(20.0)), 12),
            round(np.tan(np.radians(50.0)), 12),
        }


class TestTrapezoidalSeries:
    """A trapezoid takes an exact path, never the FFT.

    Until 2026-10-29 it was sampled, and the FFT of a slope that jumps at four
    corners converges at first order in the sample count: on the 72 deg,
    0.425-period groove below the C-method then sat 3.9e-4 from the right
    answer at every truncation, with no corner warning. ``docs/findings.md``,
    "Trapezoidal through the four solvers".
    """

    GROOVE = dict(depth_fraction=0.425, duty_cycle=0.30, sidewall_angle=72.0)

    @staticmethod
    def as_polygon(profile):
        corners = profile._corners()
        return FromProfileData(
            t=(0.0, *(x for x, _ in corners), 1.0), y=(0.0, *(y for _, y in corners), 0.0)
        )

    @staticmethod
    def finely_sampled(profile, samples=1 << 20):
        s = np.arange(samples) / samples
        t = np.mod(1.0 - s, 1.0)
        return ProfileSeries(
            160.0, "sampled", s, profile.height(t) * 160.0, -profile.slope(t)
        )

    def test_a_sharp_trapezoid_is_the_polygon_of_its_corners(self):
        trapezoid = Trapezoidal(**self.GROOVE)
        series = profile_series(Problem(period=160.0, profile=trapezoid))
        polygon = profile_series(Problem(period=160.0, profile=self.as_polygon(trapezoid)))
        assert series.kind == "polygon"
        assert np.array_equal(series.s, polygon.s) and np.array_equal(series.a, polygon.a)

    def test_the_solver_answers_as_for_the_polygon_and_warns_of_corners(self):
        trapezoid = Trapezoidal(**self.GROOVE)
        illumination = Illumination.offplane(graze=1.5, azimuth=0.0, polarization="TE")
        scans = [
            cmethod.solve(
                Problem(period=160.0, profile=profile, coating="Au"),
                illumination, [3.0], conductivity="tabulated", truncation=16,
            )
            for profile in (trapezoid, self.as_polygon(trapezoid))
        ]
        assert np.array_equal(scans[0].efficiency, scans[1].efficiency)
        assert any("corners" in w for w in scans[0].provenance.warnings)

    @pytest.mark.parametrize(
        "sidewall, depth, radius",
        [(72.0, 0.425, 0.005), (72.0, 0.425, 0.02), (45.0, 0.2, 0.01), (30.0, 0.1, 0.05)],
    )
    def test_fillets_match_a_fine_fft(self, sidewall, depth, radius):
        """The fillets are integrated in their own angle, where the integrand
        is entire; a million-sample FFT of the same profile is the reference,
        and its own error (second order, the slope being continuous) is what
        the tolerance allows for: 1e-9 measured at 72 deg, 0.005."""
        trapezoid = Trapezoidal(
            depth_fraction=depth, duty_cycle=0.30, sidewall_angle=sidewall, corner_radius=radius
        )
        series = profile_series(Problem(period=160.0, profile=trapezoid))
        reference = self.finely_sampled(trapezoid)
        j = np.arange(-128, 129)
        assert series.kind == "pieces" and len(series.arcs) == 4
        assert series.max_slope == pytest.approx(reference.max_slope, rel=1e-6)
        assert series.depth == pytest.approx(depth * 160.0, rel=1e-12)
        assert np.abs(series.slope(j) - reference.slope(j)).max() < 1e-8
        for beta in (0.5, 4.0):
            assert np.abs(series.plane_wave(beta, j) - reference.plane_wave(beta, j)).max() < 1e-10

    def test_the_largest_fillet_leaves_no_straight_wall(self):
        """At ``max_corner_radius`` a run is consumed and drops out, and the
        series still closes: the slope's mean over the period is zero."""
        base = Trapezoidal(**self.GROOVE)
        trapezoid = Trapezoidal(**self.GROOVE, corner_radius=base.max_corner_radius)
        series = profile_series(Problem(period=160.0, profile=trapezoid))
        assert len(series.s) < 9 and len(series.arcs) == 4
        assert abs(series.slope(np.array([0]))[0]) < 1e-12
        reference = self.finely_sampled(trapezoid)
        j = np.arange(-64, 65)
        assert np.abs(series.slope(j) - reference.slope(j)).max() < 1e-8

    def test_a_rounded_trapezoid_does_not_warn_of_corners(self):
        trapezoid = Trapezoidal(**self.GROOVE, corner_radius=0.01)
        scan = cmethod.solve(
            Problem(period=160.0, profile=trapezoid, coating="Au"),
            Illumination.offplane(graze=1.5, azimuth=0.0), [3.0],
            conductivity="tabulated", truncation=16,
        )
        assert not any("corners" in w for w in scan.provenance.warnings)
        assert "Gauss-Legendre" in scan.provenance.notes["profile_series"]


class TestFlatInterfaceIsFresnel:
    """Amplitude and phase: this pins every sign in the coupled system and
    the plane-wave flux, TE against r_s and TM against r_p."""

    @pytest.mark.parametrize("index", [GOLD_VISIBLE, XRAY, GLASS], ids=["gold", "xray", "glass"])
    @pytest.mark.parametrize("sin_alpha", [0.0, 0.5])
    def test_amplitudes_and_energy(self, index, sin_alpha):
        te, tm = solve_finite(
            harmonics(FLAT, 6),
            wavelength=550.0,
            period=600.0,
            sin_alpha=sin_alpha,
            sin_gamma=1.0,
            cos_gamma=0.0,
            index=index,
            incidents=[(1.0, 0.0), (0.0, 1.0)],
        )
        r_s, r_p = fresnel.amplitude(index, np.pi / 2.0 - np.arcsin(sin_alpha))
        zero = list(te.orders).index(0)
        assert te.e_amplitudes[zero] == pytest.approx(complex(r_s), abs=1e-12)
        assert tm.h_amplitudes[zero] == pytest.approx(complex(r_p), abs=1e-12)
        assert np.allclose(te.h_amplitudes, 0.0, atol=1e-14)
        for solution in (te, tm):
            assert solution.efficiencies.sum() + solution.absorption == pytest.approx(1.0, abs=1e-12)
        if index.imag == 0.0:
            # Lossless: the flux into the glass is the transmitted order.
            assert te.transmitted_efficiencies.sum() == pytest.approx(te.absorption, abs=1e-12)


class TestGorayEtSchmidtTable3:
    """Dielectric sine, n = 2, conical (their theta = 60, phi = 15 deg),
    pure E_z incidence. Their table is a comparison against Li's
    coordinate-transformation method, so this is the C-method's own
    published result. Measured: within 1.4e-5 of every order, reflected and
    transmitted; the table is quoted to about 1e-4. Their order n is our
    m = -n (as in test_integral_finite.py)."""

    REFLECTED = {0: 0.1033, 1: 0.03873, 2: 0.03741, 3: 0.01121}
    TRANSMITTED = {
        -2: 0.06351, -1: 0.5183, 0: 0.07146, 1: 0.09925,
        2: 0.04922, 3: 0.007396, 4: 0.00002466, 5: 0.0001858,
    }

    def test_every_order(self):
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3))
        (solution,) = solve_finite(
            harmonics(problem, 20, center=1),
            wavelength=300.0,
            period=600.0,
            sin_alpha=float(np.sin(np.radians(60.0))),
            sin_gamma=float(np.cos(np.radians(15.0))),
            cos_gamma=float(np.sin(np.radians(15.0))),
            index=2.0 + 0.0j,
            incidents=[(1.0, 0.0)],
        )
        got_r = by_order(solution.orders, solution.efficiencies)
        got_t = by_order(solution.transmitted_orders, solution.transmitted_efficiencies)
        assert sorted(got_t) == sorted(self.TRANSMITTED)
        for m, expected in self.REFLECTED.items():
            assert got_r[m] == pytest.approx(expected, abs=5e-5), f"R m={m}"
        for m, expected in self.TRANSMITTED.items():
            assert got_t[m] == pytest.approx(expected, abs=5e-5), f"T m={m}"
        total = solution.efficiencies.sum() + solution.transmitted_efficiencies.sum()
        assert total == pytest.approx(1.0, abs=1e-10)
        assert np.abs(solution.h_amplitudes).max() > 0.01  # polarization converts


class TestTransmittedOrdersInTheScan:
    """Into a lossless material the transmitted orders are recorded as the
    scan's ``transmission``, beside the reflected ones, and ``absorption``
    keeps only what is left over -- the truncation error."""

    #: n = 2 exactly, lossless, over the wavelengths these tests use.
    N2 = OpticalConstants(
        name="cmethod-n2", wavelength_nm=(200.0, 700.0),
        decrement=(-1.0, -1.0), absorption=(0.0, 0.0), source="test_cmethod",
    )

    @pytest.fixture(autouse=True)
    def library(self, tmp_path, monkeypatch):
        write_ari(self.N2, tmp_path / f"{self.N2.name}.ari")
        monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))

    def test_goray_schmidt_table_3_through_the_solver(self):
        """The Table 3 case end to end: their theta = 60, phi = 15 deg is
        alpha = 60 deg on a cone of half-angle gamma = 75 deg."""
        problem = Problem(
            period=600.0, profile=Sinusoidal(depth_fraction=0.3), coating=self.N2.name
        )
        illumination = Illumination.offplane(graze=75.0, azimuth=60.0, polarization="TE")
        scan = cmethod.solve(problem, illumination, [300.0], truncation=20, conductivity="tabulated")
        assert scan.transmission is not None
        reflected = by_order(scan.orders, scan.efficiency[0])
        transmitted = by_order(scan.transmission.orders, scan.transmission.efficiency[0])
        live = set(scan.transmission.orders[scan.transmission.propagating[0]].tolist())
        assert live == set(TestGorayEtSchmidtTable3.TRANSMITTED)
        for m, expected in TestGorayEtSchmidtTable3.REFLECTED.items():
            assert reflected[m] == pytest.approx(expected, abs=5e-5), f"R m={m}"
        for m, expected in TestGorayEtSchmidtTable3.TRANSMITTED.items():
            assert transmitted[m] == pytest.approx(expected, abs=5e-5), f"T m={m}"
        assert abs(scan.absorption[0]) < 1e-10  # lossless: the remainder only
        assert check_energy_balance(scan, tolerance=1e-10, lossless=True).passed

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_flat_glass_transmits_the_fresnel_fraction(self, polarization):
        flat = Problem(period=600.0, profile=FLAT.profile, coating=self.N2.name)
        illumination = Illumination.classical(alpha=30.0, polarization=polarization)
        scan = cmethod.solve(flat, illumination, [550.0], truncation=4, conductivity="tabulated")
        r_s, r_p = fresnel.amplitude(2.0 + 0.0j, np.radians(60.0))
        r = complex(r_s if polarization == "TE" else r_p)
        assert scan.transmission.order(0)[0] == pytest.approx(1.0 - abs(r) ** 2, abs=1e-12)
        assert scan.transmission.total[0] == pytest.approx(1.0 - abs(r) ** 2, abs=1e-12)

    def test_an_absorbing_material_records_no_transmission(self):
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.05), coating="Au")
        scan = cmethod.solve(
            problem, Illumination.offplane(graze=1.5, azimuth=25.0), [2.0],
            truncation=24, conductivity="tabulated",
        )
        assert scan.transmission is None
        assert scan.absorption[0] > 0.1

    def test_truncation_must_hold_every_transmitted_order(self):
        """In n = 2 at 300 nm the substrate admits orders -2 to 5 against the
        vacuum's 0 to 3; a truncation that holds the vacuum's but not the
        substrate's would drop transmitted orders, so it is refused."""
        problem = Problem(
            period=600.0, profile=Sinusoidal(depth_fraction=0.3), coating=self.N2.name
        )
        illumination = Illumination.offplane(graze=75.0, azimuth=60.0, polarization="TE")
        cmethod.solve(problem, illumination, [300.0], truncation=4, conductivity="perfect")
        with pytest.raises(ValueError, match="transmitted"):
            cmethod.solve(problem, illumination, [300.0], truncation=4, conductivity="tabulated")


class TestPerfectConductorBenchmarks:
    def test_marechal_stroke(self):
        """A 30/60 deg echelette illuminated normal to its blaze facet sends
        all TM power into the blaze order (M&P 4.7.2). The exact field is two
        plane waves, so the C-method lands on it at any truncation; this pins
        the plane-wave flux and that the blaze is on the documented side."""
        problem = Problem(period=600.0, profile=Blazed(blaze_angle=30.0, antiblaze_angle=60.0))
        solution = solve_perfect(
            harmonics(problem, 8),
            wavelength=600.0, period=600.0, sin_alpha=0.5, sin_gamma=1.0, polarization="TM",
        )
        got = by_order(solution.orders, solution.efficiencies)
        assert got[1] == pytest.approx(1.0, abs=1e-12)
        assert got[0] == pytest.approx(0.0, abs=1e-12)

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_maystre_popov_table_4_1(self, polarization):
        """Sinusoid 600 nm, depth 180 nm, Littrow at 600 nm. TM matches the
        printed 0.9579. TE does not match the printed 0.4659: this solver
        and the integral solver both give 0.466441 (the integral at every
        mesh from 110 to 512 points), so the gap is in the table, and the
        pin is the cross-method value."""
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3))
        solution = solve_perfect(
            harmonics(problem, 16),
            wavelength=600.0, period=600.0, sin_alpha=0.5, sin_gamma=1.0,
            polarization=polarization,
        )
        got = by_order(solution.orders, solution.efficiencies)
        expected = {"TE": 0.466441, "TM": 0.9579}[polarization]
        assert got[1] == pytest.approx(expected, abs=1e-5 if polarization == "TE" else 1e-4)


class TestEnergyTheorem:
    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    @pytest.mark.parametrize("graze", [90.0, 30.0], ids=["in-plane", "conical"])
    def test_perfect_conductor_sums_to_one(self, polarization, graze):
        illumination = (
            Illumination.classical(alpha=20.0, polarization=polarization)
            if graze == 90.0
            else Illumination.offplane(graze=graze, azimuth=20.0, polarization=polarization)
        )
        scan = cmethod.solve(SINE, illumination, [350.0, 400.0, 600.0], truncation=16)
        assert check_energy_balance(scan, tolerance=1e-10, lossless=True).passed
        assert scan.provenance.notes["energy_balance_deviation"] < 1e-10

    def test_perfect_conductor_limit_of_the_coupled_system(self):
        """As |n| grows the conical coupled system must approach the
        decoupled perfect-conductor solve, with cross-polarization and the
        difference both falling as 1/|n| (Leontovich). Measured at 100i and
        300i: differences 1.5e-3 and 4.9e-4 (TE), cross-polarization 8.5e-3
        and 2.8e-3. A real negative permittivity is the case this exercises,
        and the one that once fell through the lossless branch."""
        h = harmonics(SINE, 16)
        args = dict(
            wavelength=280.0, period=500.0, sin_alpha=0.2,
            sin_gamma=float(np.sin(np.radians(30.0))),
        )
        cos_gamma = float(np.cos(np.radians(30.0)))
        perfect = solve_perfect(h, polarization="TE", **args)
        gaps, cross = [], []
        for n in (100j, 300j):
            (finite,) = solve_finite(h, cos_gamma=cos_gamma, index=n, incidents=[(1.0, 0.0)], **args)
            gaps.append(np.abs(finite.efficiencies - perfect.efficiencies).max())
            cross.append(np.abs(finite.h_amplitudes).max())
            assert finite.absorption == pytest.approx(0.0, abs=1e-12)  # lossless
        assert gaps[1] < 2e-3 and cross[1] < 1e-2
        assert gaps[0] / gaps[1] == pytest.approx(3.0, rel=0.2)
        assert cross[0] / cross[1] == pytest.approx(3.0, rel=0.2)


class TestAgreementWithTheIntegralSolver:
    """A different rigorous method, sharing nothing below ``Problem``."""

    @pytest.mark.parametrize("polarization, tolerance", [("TE", 1e-6), ("TM", 3e-5)])
    def test_perfect_conductor_sinusoid(self, polarization, tolerance):
        illumination = Illumination.classical(alpha=20.0, polarization=polarization)
        ours = cmethod.solve(SINE, illumination, [400.0], truncation=16)
        theirs = integral.solve(SINE, illumination, [400.0], boundary_points=400)
        assert np.abs(ours.efficiency - theirs.efficiency).max() < tolerance

    @pytest.mark.parametrize("incident", [(1.0, 0.0), (0.0, 1.0)], ids=["Ez", "Hz"])
    @pytest.mark.parametrize("graze", [90.0, 60.0], ids=["in-plane", "conical"])
    def test_finite_conductivity_gold(self, incident, graze):
        """Measured differences 2e-6 to 1.7e-5, the integral's own error at
        512 points; asserted at 5e-5 with 384."""
        sin_gamma = float(np.sin(np.radians(graze)))
        cos_gamma = 0.0 if graze == 90.0 else float(np.cos(np.radians(graze)))
        (ours,) = solve_finite(
            harmonics(SINE, 16), wavelength=600.0, period=500.0,
            sin_alpha=float(np.sin(np.radians(20.0))), sin_gamma=sin_gamma,
            cos_gamma=cos_gamma, index=GOLD_VISIBLE, incidents=[incident],
        )
        theirs = integral_finite(
            SINE, wavelength=600.0, sin_alpha=float(np.sin(np.radians(20.0))),
            index=GOLD_VISIBLE, cos_gamma=cos_gamma, incident=incident, points=384,
        )
        assert np.array_equal(ours.orders, theirs.orders)
        assert np.abs(ours.efficiencies - theirs.efficiencies).max() < 5e-5
        assert ours.absorption == pytest.approx(theirs.absorption, abs=5e-5)

    def test_off_plane_xray_gold_through_the_solver(self):
        """The founding application's mount, Au from the shipped CXRO table,
        absolute efficiencies with absorption. Measured 8.8e-5 per order
        and 5.4e-5 in absorption against the integral at 512 points. The gap
        is the integral's: this solve is converged (N = 24 and 48 agree to
        1e-14), and the integral closes on it, 2.6e-5 at 768 points, as its
        own energy deviation falls from 2.8e-4 to 8e-5."""
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.05), coating="Au")
        illumination = Illumination.offplane(graze=1.5, azimuth=25.0, polarization="TE")
        ours = cmethod.solve(problem, illumination, [2.0, 5.0], truncation=24, conductivity="tabulated")
        theirs = integral.solve(
            problem, illumination, [2.0, 5.0], boundary_points=512, conductivity="tabulated"
        )
        assert np.abs(ours.efficiency - theirs.efficiency).max() < 2e-4
        assert np.abs(ours.absorption - theirs.absorption).max() < 2e-4
        assert check_energy_balance(ours, tolerance=1e-10, lossless=True).passed
        assert ours.provenance.notes["normalization"] == "absolute"
        assert ours.provenance.notes["material"] == "Au"

    def test_shallow_sawtooth_converges_to_it(self):
        """Corners: the C-method approaches the integral solver's graded-mesh
        answer slowly. Measured on 10/30 deg: 2e-5 at N = 80, and the gap
        shrinks with N."""
        problem = Problem(period=1000.0, profile=Blazed(blaze_angle=10.0, antiblaze_angle=30.0))
        illumination = Illumination.classical(alpha=15.0, polarization="TM")
        reference = integral.solve(
            problem, illumination, [600.0], boundary_points=512, corner_grading=4
        ).efficiency
        gaps = [
            np.abs(cmethod.solve(problem, illumination, [600.0], truncation=n).efficiency
                   - reference).max()
            for n in (20, 64)
        ]
        assert gaps[1] < gaps[0]
        assert gaps[1] < 1e-4


class TestReciprocity:
    def test_perfect_conductor_sawtooth(self):
        """An asymmetric profile, where reciprocity is not a symmetry of the
        geometry. Measured 3.5e-7 (TM) and 1.4e-6 (TE) at N = 40, which is
        the truncation error on a cornered profile."""
        problem = Problem(period=1000.0, profile=Blazed(blaze_angle=10.0, antiblaze_angle=30.0))
        for polarization in ("TE", "TM"):
            report = check_reciprocity(
                cmethod, problem,
                Illumination.classical(alpha=15.0, polarization=polarization),
                [600.0], tolerance=1e-5, truncation=40,
            )
            assert report.passed, report

    def test_absorbing_off_plane(self):
        """Measured 2.1e-7 on Au at 2 nm, gamma = 1.5 deg."""
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.05), coating="Au")
        report = check_reciprocity(
            cmethod, problem,
            Illumination.offplane(graze=1.5, azimuth=25.0, polarization="TE"),
            [2.0], tolerance=1e-6, max_orders=4, truncation=24, conductivity="tabulated",
        )
        assert report.passed, report


class TestConvergence:
    def test_a_smooth_profile_plateaus(self):
        report = check_convergence(
            cmethod,
            Problem(period=500.0, profile=Sinusoidal(depth_fraction=0.2)),
            Illumination.classical(alpha=10.0, polarization="TM"),
            [450.0, 600.0],
            values=(4, 8, 16, 32),
        )
        assert report.passed

    def test_unpolarized_is_the_mean(self):
        illumination = Illumination.offplane(graze=30.0, azimuth=20.0, polarization="TE")
        scans = {
            p: cmethod.solve(SINE, illumination.model_copy(update={"polarization": p}),
                             [400.0], truncation=12)
            for p in ("TE", "TM", "unpolarized")
        }
        mean = 0.5 * (scans["TE"].efficiency + scans["TM"].efficiency)
        assert np.allclose(scans["unpolarized"].efficiency, mean, atol=1e-14)


class TestWhereItIsWeak:
    """Both failure regions the method note measures, each reported."""

    def test_corners_warn_that_the_balance_is_not_an_error_estimate(self):
        """Blazed 20/70, TE: the balance closes to 2e-5 at N = 40 while the
        efficiencies are still 6e-4 from the converged answer. A reader who
        took the balance as the error would be misled by an order of
        magnitude and more, so the solver says so."""
        problem = Problem(period=1000.0, profile=Blazed(blaze_angle=20.0, antiblaze_angle=70.0))
        illumination = Illumination.classical(alpha=15.0, polarization="TE")
        scan = cmethod.solve(problem, illumination, [600.0], truncation=40)
        reference = integral.solve(
            problem, illumination, [600.0], boundary_points=512, corner_grading=4
        )
        balance = scan.provenance.notes["energy_balance_deviation"]
        error = np.abs(scan.efficiency - reference.efficiency).max()
        assert error > 10 * balance
        assert any("corners" in w for w in scan.provenance.warnings)

    def test_a_smooth_profile_does_not_warn_about_corners(self):
        scan = cmethod.solve(SINE, Illumination.classical(alpha=10.0), [450.0], truncation=8)
        assert not any("corners" in w for w in scan.provenance.warnings)

    def test_a_deep_profile_loses_digits_and_says_so(self):
        """Sinusoid one period deep, TM: right at N = 16 and wrong by 1e-2 at
        N = 32, as the evanescent modes exhaust double precision. The rung
        that went wrong carries the condition warning, and its energy balance
        shows the error. The reference is the integral
        solver at 768 points, 0.844609 in the strongest order (measured
        2026-10-27; its own error there is about 1e-5), pinned rather than
        recomputed because that solve takes 40 s."""
        problem = Problem(period=500.0, profile=Sinusoidal(depth_fraction=1.0))
        illumination = Illumination.classical(alpha=10.0, polarization="TM")
        good = cmethod.solve(problem, illumination, [450.0], truncation=16)
        bad = cmethod.solve(problem, illumination, [450.0], truncation=32)
        assert good.efficiency.max() == pytest.approx(0.844609, abs=2e-4)
        assert not good.provenance.warnings
        assert any("condition number" in w for w in bad.provenance.warnings)
        # How wrong a precision-exhausted solve is depends on BLAS rounding
        # (0.833 single-threaded, 0.872 threaded, same machine; 2.8e-4 off on
        # CI's Linux py3.11 runner, run 40626543858, where a fixed 1e-3 bound
        # failed), so the pin is relative: the higher truncation is clearly
        # worse than the lower one, and its energy balance shows it.
        good_error = abs(good.efficiency.max() - 0.844609)
        bad_error = abs(bad.efficiency.max() - 0.844609)
        assert bad_error > 10 * good_error
        good_balance = good.provenance.notes["energy_balance_deviation"]
        bad_balance = bad.provenance.notes["energy_balance_deviation"]
        assert bad_balance > 1e-6 and bad_balance > 1e3 * good_balance


class TestRefusals:
    @pytest.mark.parametrize(
        "profile",
        [
            Lamellar(depth_fraction=0.2),
            Trapezoidal(depth_fraction=0.2),
            Trapezoidal(depth_fraction=0.2, corner_radius=0.05),
            Blazed(blaze_angle=30.0),  # vertical anti-blaze facet
            FromProfileData(t=(0.0, 0.5, 0.5, 1.0), y=(0.0, 0.0, 0.2, 0.0)),
            FromProfileData(t=(0.0, 0.5, 1.0), y=(0.0, 0.2, 0.1)),  # step at the wrap
        ],
        ids=["lamellar", "vertical-trapezoid", "vertical-rounded-trapezoid", "vertical-facet", "vertical-segment", "step-at-the-ends"],
    )
    def test_infinite_slopes(self, profile):
        with pytest.raises(UnsupportedConfiguration):
            cmethod.solve(Problem(period=500.0, profile=profile), Illumination.classical(alpha=10.0), [450.0])

    def test_undercut(self):
        undercut = FromProfileData(t=(0.0, 0.55, 0.45, 1.0), y=(0.0, 0.4, 0.1, 0.0))
        with pytest.raises(UnsupportedConfiguration):
            cmethod.solve(Problem(period=500.0, profile=undercut), Illumination.classical(alpha=10.0), [450.0])

    def test_unknown_boundary_condition(self):
        with pytest.raises(UnsupportedConfiguration):
            cmethod.solve(SINE, Illumination.classical(alpha=10.0), [450.0], conductivity="leontovich")

    def test_tabulated_without_a_material(self):
        with pytest.raises(UnsupportedConfiguration, match="needs a material"):
            cmethod.solve(SINE, Illumination.classical(alpha=10.0), [450.0], conductivity="tabulated")

    @pytest.mark.parametrize("truncation", [0, 1, 257, True, 4.0])
    def test_truncation_outside_its_range(self, truncation):
        """1 is below what orders -1 to 0 need with an evanescent harmonic
        beyond each end; 257 above the cost ceiling."""
        with pytest.raises(ValueError):
            cmethod.solve(SINE, Illumination.classical(alpha=10.0), [450.0], truncation=truncation)

    def test_a_wavelength_exactly_on_a_rayleigh_anomaly(self):
        """sin(alpha) = 0.2 and lambda = 0.6 P puts order -1 at grazing."""
        with pytest.raises(ValueError, match="Rayleigh anomaly"):
            cmethod.solve(
                SINE,
                Illumination.classical(alpha=float(np.degrees(np.arcsin(0.2)))),
                [300.0],
            )

    def test_a_coating_is_not_consulted_in_perfect_mode(self):
        problem = Problem(period=500.0, profile=Sinusoidal(depth_fraction=0.3), coating="Au")
        scan = cmethod.solve(problem, Illumination.classical(alpha=10.0), [450.0], truncation=8)
        assert any("not consulted" in w for w in scan.provenance.warnings)
        assert scan.absorption is None


def asymmetric_groove(steep: float, *, depth=0.3, facet=27.0, land=77 / 400) -> FromProfileData:
    """A 27 deg facet up to a flat land, then a steep wall down to the floor:
    the shape of issue #7's TASTE grooves, whose wall is at 88.8 deg."""
    up = depth / np.tan(np.radians(facet))
    down = depth / np.tan(np.radians(steep))
    floor = 1.0 - up - land - down
    return FromProfileData(
        t=(0.0, floor, floor + up, floor + up + land, 1.0),
        y=(0.0, 0.0, depth, depth, 0.0),
    )


#: Issue #7's mount: X-ray gold, 1.73 deg off-plane.
OFF_PLANE = Illumination.offplane(graze=1.73, azimuth=0.0, polarization="TE")

#: The integral solver on ``asymmetric_groove(80)`` at 1.55 nm, OFF_PLANE,
#: CXRO gold, corner_grading=4: each propagating order, ascending. 1600
#: points; 1200 agree to 6e-6. Pinned because the 1600-point solve took
#: ten minutes (measured 2026-10-29).
STEEP_80_REFERENCE = (
    0.007437, 0.006005, 0.007763, 0.026137, 0.003908, 0.060704, 0.021511,
    0.021711, 0.068545, 0.064224, 0.035783, 0.008218, 0.002726, 0.262197,
    0.049879,
)


class TestSteepWalls:
    """Issue #7. A wall the C-method cannot resolve is refused before any
    solve, rather than answered with efficiencies of 31. The floor and the
    measurements behind it are in ``docs/theory/cmethod.md`` section 6."""

    @pytest.mark.parametrize(
        "profile",
        [
            asymmetric_groove(88.8),
            Trapezoidal(depth_fraction=0.3, duty_cycle=0.3, sidewall_angle=88.8),
            Trapezoidal(
                depth_fraction=0.3, duty_cycle=0.3, sidewall_angle=88.8, corner_radius=0.02
            ),
        ],
        ids=["asymmetric-polygon", "trapezoidal", "trapezoidal-filleted"],
    )
    def test_the_issue_7_wall_is_refused(self, profile):
        """Unrefused, the asymmetric groove returned efficiencies of 1.9e4 at
        N = 16 and 11 at N = 128, where the integral solver gives 0.26."""
        problem = Problem(period=400.0, profile=profile, coating="Au")
        with pytest.raises(UnsupportedConfiguration, match="too steep") as caught:
            cmethod.solve(problem, OFF_PLANE, [1.55, 3.1], truncation=64, conductivity="tabulated")
        assert "88.8 deg" in str(caught.value)
        assert f"beyond the {_MAX_TRUNCATION}" in str(caught.value)

    def test_a_steep_enough_wall_is_refused_at_any_wavelength(self):
        """89.5 deg needs about 3 harmonics per unit of slope however long the
        wavelength; in-plane at 500 nm N = 192 was still 1e-2 off."""
        problem = Problem(
            period=400.0,
            profile=Trapezoidal(depth_fraction=0.3, duty_cycle=0.3, sidewall_angle=89.5),
        )
        with pytest.raises(UnsupportedConfiguration, match="too steep"):
            cmethod.solve(problem, Illumination.classical(alpha=20.0), [500.0])

    def test_too_few_harmonics_names_how_many(self):
        """85 deg at 1.55 nm off-plane: 64 harmonics were 5e-2 off and 32 were
        6.7 off, so both are refused with the number needed."""
        problem = Problem(period=400.0, profile=asymmetric_groove(85.0), coating="Au")
        floor = wall_floor(problem, OFF_PLANE, [1.55])
        assert 64 < floor.needed <= 128
        with pytest.raises(ValueError, match=f"needs about {floor.needed} harmonics"):
            cmethod.solve(problem, OFF_PLANE, [1.55], truncation=64, conductivity="tabulated")

    def test_the_floor_reads_the_boundary_not_the_profile_type(self):
        """A Trapezoidal (sampled) and the same trapezoid as a polygon get the
        same floor, and rounding its corners does not lower it: the fillets
        are tangent to the wall, so the wall's slope is still on the curve."""
        sharp = Trapezoidal(depth_fraction=0.3, duty_cycle=0.3, sidewall_angle=85.0)
        corners = sharp._corners()
        polygon = FromProfileData(
            t=(0.0, *[c[0] for c in corners], 1.0), y=(0.0, *[c[1] for c in corners], 0.0)
        )
        rounded = sharp.model_copy(update={"corner_radius": 0.02})
        floors = [
            wall_floor(Problem(period=400.0, profile=p), OFF_PLANE, [1.55]).needed
            for p in (sharp, polygon, rounded)
        ]
        assert floors[0] == pytest.approx(floors[1], rel=0.05)
        assert floors[2] == pytest.approx(floors[0], rel=0.1)

    def test_the_floor_falls_with_the_reduced_wavelength(self):
        """The same wall at 6 nm off-plane, or in visible light, needs fewer."""
        problem = Problem(period=400.0, profile=asymmetric_groove(85.0))
        short = wall_floor(problem, OFF_PLANE, [1.55]).needed
        longer = wall_floor(problem, OFF_PLANE, [6.0]).needed
        visible = wall_floor(problem, Illumination.classical(alpha=20.0), [500.0]).needed
        assert short > longer > visible

    def test_gentle_profiles_are_untouched(self):
        """Every geometry the rest of this file solves sits far below it."""
        cases = [
            (SINE, Illumination.classical(alpha=10.0), [450.0], 8),
            (Problem(period=500.0, profile=Sinusoidal(depth_fraction=1.0)),
             Illumination.classical(alpha=10.0), [450.0], 16),
            (Problem(period=1000.0, profile=Blazed(blaze_angle=20.0, antiblaze_angle=70.0)),
             Illumination.classical(alpha=15.0), [600.0], 40),
            (Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.05)),
             Illumination.offplane(graze=1.5, azimuth=25.0), [2.0, 5.0], 24),
        ]
        for problem, illumination, wavelengths, truncation in cases:
            assert wall_floor(problem, illumination, wavelengths).needed <= truncation / 2

    def test_inside_the_limit_it_agrees_with_the_integral_solver(self):
        """An 85 deg trapezoid, perfect conductor, TM at 500 nm: the floor is
        26, and at 64 harmonics the C-method was 1.2e-4 from the integral
        solver's graded-mesh answer at 1200 points; asserted at 512. Given
        as a polygon: a sharp Trapezoidal goes down the sampled series,
        whose FFT of the slope's jumps carries an error of its own (3e-3
        here) that is not this test's subject."""
        corners = Trapezoidal(
            depth_fraction=0.3, duty_cycle=0.3, sidewall_angle=85.0
        )._corners()
        problem = Problem(
            period=400.0,
            profile=FromProfileData(
                t=(0.0, *[c[0] for c in corners], 1.0),
                y=(0.0, *[c[1] for c in corners], 0.0),
            ),
        )
        illumination = Illumination.classical(alpha=20.0, polarization="TM")
        assert wall_floor(problem, illumination, [500.0]).needed <= 32
        ours = cmethod.solve(problem, illumination, [500.0], truncation=64)
        theirs = integral.solve(
            problem, illumination, [500.0], boundary_points=512, corner_grading=4
        )
        assert np.abs(ours.efficiency - theirs.efficiency).max() < 5e-4

    def test_off_plane_gold_inside_the_limit(self):
        """The issue's mount and groove with an 80 deg wall, which the method
        can resolve: the floor is 44, and at 128 harmonics every order is
        within 1.3e-4 of the pinned integral answer (measured)."""
        problem = Problem(period=400.0, profile=asymmetric_groove(80.0), coating="Au")
        assert wall_floor(problem, OFF_PLANE, [1.55]).needed <= 64
        scan = cmethod.solve(problem, OFF_PLANE, [1.55], truncation=128, conductivity="tabulated")
        ours = scan.efficiency[0][scan.propagating[0]]
        assert np.abs(ours - np.array(STEEP_80_REFERENCE)).max() < 3e-4

    def test_between_the_floor_and_twice_it_the_solver_warns(self):
        """There the balance under-reports: at 87 deg and N = 128 the integral
        solver was 2.5e-2 away while the balance read 2e-3."""
        problem = Problem(period=400.0, profile=asymmetric_groove(80.0), coating="Au")
        floor = wall_floor(problem, OFF_PLANE, [1.55]).needed
        near = cmethod.solve(problem, OFF_PLANE, [1.55], truncation=floor, conductivity="tabulated")
        assert any("steep wall" in w for w in near.provenance.warnings)
        far = cmethod.solve(
            problem, OFF_PLANE, [1.55], truncation=2 * floor, conductivity="tabulated"
        )
        assert not any("steep wall" in w for w in far.provenance.warnings)
