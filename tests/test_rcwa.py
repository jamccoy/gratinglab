"""The RCWA backend, against references that do not depend on it.

In order of strictness, as in ``docs/theory/rcwa.md`` section 8: a flat
interface is Fresnel, a single film is the Airy formula, a lossless grating
conserves flux, a metallic one is reciprocal, Li's factorization rule is shown
to be in effect by converging where Laurent's does not, and the integral
solver, a different method on the same ``Problem``, agrees order by order.

None of it needs the private corpus. The materials are written into a library
of the test's own, as a user's would be.
"""

import numpy as np
import pytest

from gratinglab.checks import check_energy_balance, check_reciprocity
from gratinglab.convergence import check_convergence
from gratinglab.illumination import Illumination
from gratinglab.materials import fresnel
from gratinglab.materials.optical import OpticalConstants, write_ari
from gratinglab.problem import ConformalLayer, Problem
from gratinglab.profiles import Blazed, FromProfileData, Lamellar, Sinusoidal
from gratinglab.solvers import UnsupportedConfiguration, get_solver
from gratinglab.solvers.rcwa._core import incident_field, solve_stack
from gratinglab.solvers.rcwa._stack import (
    Slab,
    Stack,
    _offset_height,
    _surface,
    groove_geometry,
)

GLASS = 1.5 + 0.0j
#: A visible-light metal, of the size of gold or silver in the red: |eps| ~ 9,
#: strongly negative real part, which is where Laurent's rule fails in TM.
METAL = 0.2 + 3.0j

rcwa = get_solver("rcwa")
integral = get_solver("integral")


def _table(name, index):
    return OpticalConstants(
        name=name,
        wavelength_nm=(200.0, 1000.0),
        decrement=(1.0 - index.real,) * 2,
        absorption=(index.imag,) * 2,
        source="test_rcwa",
    )


@pytest.fixture
def library(tmp_path, monkeypatch):
    for name, index in (("glass", GLASS), ("metal", METAL), ("gain", 1.5 - 0.1j)):
        write_ari(_table(name, index), tmp_path / f"{name}.ari")
    monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))
    return tmp_path


def kvector(alpha_deg, gamma_deg):
    """``(kx0, ky, kz0)`` over ``k0``, as the backend forms them."""
    a, g = np.radians(alpha_deg), np.radians(gamma_deg)
    return -np.sin(a) * np.sin(g), np.cos(g), np.cos(a) * np.sin(g)


def solve(stack, *, alpha=30.0, gamma=90.0, polarizations=("TE", "TM"),
          truncation=4, wavelength=500.0, period=600.0, factorization="li"):
    kx0, ky, kz0 = kvector(alpha, gamma)
    return solve_stack(
        stack,
        wavelength=wavelength,
        period=period,
        kx0=kx0,
        ky=ky,
        incidents=[incident_field(p, kx0, ky, kz0) for p in polarizations],
        truncation=truncation,
        factorization=factorization,
    )


def graze_of(alpha, gamma):
    """The true angle of incidence onto a flat surface, as a graze angle."""
    a, g = np.radians(alpha), np.radians(gamma)
    return np.pi / 2 - np.arccos(np.cos(a) * np.sin(g))


class TestFlatInterfaceIsFresnel:
    """No grating at all: the specular order must carry exactly the Fresnel
    reflectivity of ``materials.fresnel``, which shares no code with this
    backend. This is the test that pins the signs of the external regions,
    the gap medium and the flux normalisation."""

    @pytest.mark.parametrize("index", [GLASS, METAL], ids=["glass", "metal"])
    @pytest.mark.parametrize("alpha", [0.0, 30.0, 70.0])
    def test_in_plane_te_is_s_and_tm_is_p(self, index, alpha):
        solution = solve(Stack(1.0, (), index), alpha=alpha)
        r_s, r_p = fresnel.amplitude(index, graze_of(alpha, 90.0))
        specular = solution.reflected[:, solution.harmonics == 0][:, 0]
        assert specular == pytest.approx([abs(r_s) ** 2, abs(r_p) ** 2], abs=1e-13)
        assert solution.reflected[:, solution.harmonics != 0].max() < 1e-15

    @pytest.mark.parametrize("index", [GLASS, METAL], ids=["glass", "metal"])
    def test_conical_unpolarized_is_the_fresnel_average(self, index):
        """Off-plane, TE and TM by the grooves are not s and p, but they are
        an orthogonal pair, so their average is unpolarized light: the mean of
        Rs and Rp at the true angle of incidence. Each alone differs from it,
        which is the evidence the polarization really is by the grooves."""
        solution = solve(Stack(1.0, (), index), alpha=30.0, gamma=40.0)
        r_s, r_p = fresnel.amplitude(index, graze_of(30.0, 40.0))
        average = 0.5 * (abs(r_s) ** 2 + abs(r_p) ** 2)
        totals = solution.reflected.sum(axis=1)
        assert totals.mean() == pytest.approx(average, abs=1e-13)
        assert abs(totals[0] - totals[1]) > 0.01

    @pytest.mark.parametrize("index", [GLASS, METAL], ids=["glass", "metal"])
    def test_what_is_not_reflected_enters(self, index):
        solution = solve(Stack(1.0, (), index), alpha=30.0, gamma=40.0)
        assert solution.reflected.sum(axis=1) + solution.absorbed == pytest.approx(
            1.0, abs=1e-13
        )


def airy(index_film, thickness, index_substrate, graze, wavelength):
    """Single-film reflectivity ``(Rs, Rp)``, by the textbook Airy sum."""
    cos2 = np.cos(graze) ** 2

    def kz(n):
        return np.sqrt(n**2 - cos2 + 0j)

    k0, k1, k2 = kz(1.0), kz(index_film), kz(index_substrate)
    phase = np.exp(2j * k1 * thickness * 2 * np.pi / wavelength)
    out = []
    for w0, w1, w2 in ((1.0, 1.0, 1.0), (1.0, index_film**2, index_substrate**2)):
        r01 = (k0 / w0 - k1 / w1) / (k0 / w0 + k1 / w1)
        r12 = (k1 / w1 - k2 / w2) / (k1 / w1 + k2 / w2)
        out.append(abs((r01 + r12 * phase) / (1 + r01 * r12 * phase)) ** 2)
    return out


class TestLayerStack:
    """The layer stack is in the data model from the start (roadmap section
    1). ``Problem`` cannot describe one yet, so it is exercised here at the
    level of ``Stack``."""

    FILM = Slab(thickness=137.0, background=1.38 + 0.0j)

    @pytest.mark.parametrize("substrate", [GLASS, METAL], ids=["glass", "metal"])
    @pytest.mark.parametrize("alpha", [0.0, 45.0])
    def test_one_film_is_the_airy_formula(self, substrate, alpha):
        solution = solve(Stack(1.0, (), substrate, above=(self.FILM,)), alpha=alpha)
        expected = airy(1.38, 137.0, substrate, graze_of(alpha, 90.0), 500.0)
        assert solution.reflected.sum(axis=1) == pytest.approx(expected, abs=1e-12)

    def test_one_film_conical_unpolarized(self):
        solution = solve(Stack(1.0, (), METAL, above=(self.FILM,)), alpha=25.0, gamma=50.0)
        expected = airy(1.38, 137.0, METAL, graze_of(25.0, 50.0), 500.0)
        assert solution.reflected.sum(axis=1).mean() == pytest.approx(
            np.mean(expected), abs=1e-12
        )

    def test_above_and_below_are_the_same_layer_without_a_grating(self):
        a = solve(Stack(1.0, (), GLASS, above=(self.FILM,)), gamma=50.0)
        b = solve(Stack(1.0, (), GLASS, below=(self.FILM,)), gamma=50.0)
        assert np.allclose(a.reflected, b.reflected, atol=1e-14)

    def test_a_film_under_a_grating_carries_the_same_flux_at_every_plane(self):
        """Lossless throughout, so the downward flux through every interior
        plane is the same number, and it is what is not reflected. The planes
        are evaluated from the partial S-matrices above and below each one,
        so this checks the interior field recovery, not only the ends."""
        ridge = Slab(200.0, 1.0, ((0.3, 0.7, GLASS),))
        stack = Stack(1.0, (ridge,), GLASS, below=(self.FILM,))
        solution = solve(stack, gamma=50.0, truncation=12)
        assert solution.plane_flux.shape == (2, 3)
        assert np.ptp(solution.plane_flux, axis=1).max() < 1e-12
        assert solution.reflected.sum(axis=1) + solution.absorbed == pytest.approx(
            1.0, abs=1e-12
        )
        assert solution.transmitted.sum(axis=1) == pytest.approx(
            solution.plane_flux[:, -1], abs=1e-12
        )

    def test_absorbing_ridges_take_power_between_the_planes(self):
        ridge = Slab(200.0, 1.0, ((0.3, 0.7, METAL),))
        solution = solve(Stack(1.0, (ridge,), GLASS), truncation=12)
        top, bottom = solution.plane_flux[:, 0], solution.plane_flux[:, 1]
        assert (top - bottom > 1e-3).all()


# ---------------------------------------------------------------------------
# Through the public solver
# ---------------------------------------------------------------------------


def lamellar(material, depth=0.3, duty=0.4, period=600.0):
    return Problem(
        period=period,
        profile=Lamellar(depth_fraction=depth, duty_cycle=duty),
        coating=material,
    )


class TestEnergyBalance:
    """``R + T + A = 1``. For RCWA this is a check of the implementation and
    not of convergence (``rcwa.md`` section 4), so it is asserted at rounding
    level and at a deliberately low truncation."""

    @pytest.mark.parametrize("gamma", [90.0, 50.0])
    @pytest.mark.parametrize("material", ["glass", "metal"])
    def test_reflected_transmitted_and_absorbed_sum_to_one(
        self, library, material, gamma
    ):
        scan = rcwa.solve(
            lamellar(material),
            Illumination(alpha_deg=20.0, gamma_deg=gamma, polarization="unpolarized"),
            [400.0, 550.0, 700.0],
            truncation=6,
        )
        assert check_energy_balance(scan, lossless=True, tolerance=1e-10).passed
        if material == "glass":
            # Lossless: everything that enters is transmitted, into orders.
            assert scan.transmission is not None
            assert (scan.transmission.total > 0.5).all()
            assert np.abs(scan.absorption).max() < 1e-12
        else:
            # Nothing propagates in a metal, so what enters is absorbed.
            assert scan.transmission is None
            assert (scan.absorption > 0.0).all()


class TestTransmission:
    def test_a_flat_interface_transmits_the_fresnel_fraction(self):
        """``Ts = (kz2/kz1)|ts|^2`` and the TM equivalent, formed here from
        the textbook amplitudes and not from ``1 - R``: the transmitted
        amplitudes are a separate output of the solve."""
        solution = solve(Stack(1.0, (), GLASS), alpha=40.0, gamma=60.0)
        kz1 = np.cos(np.radians(40.0)) * np.sin(np.radians(60.0))
        kx = -np.sin(np.radians(40.0)) * np.sin(np.radians(60.0))
        kz2 = np.sqrt(GLASS.real**2 - kx**2 - np.cos(np.radians(60.0)) ** 2)
        eps = GLASS.real**2
        t_s = 2 * kz1 / (kz1 + kz2)  # E-field ratio
        t_p = 2 * eps * kz1 / (eps * kz1 + kz2)  # H-field ratio
        expected = 0.5 * (kz2 / kz1 * t_s**2 + kz2 / (eps * kz1) * t_p**2)
        assert solution.transmitted.sum(axis=1).mean() == pytest.approx(
            expected, abs=1e-13
        )

    def test_the_substrate_admits_more_orders_than_the_vacuum(self, library):
        """At 400 nm on a 600 nm period, alpha = 20: orders -1 and 0 reflect;
        glass at n = 1.5 also transmits order -2 and +1."""
        scan = rcwa.solve(
            lamellar("glass"), Illumination.classical(alpha=20.0), [400.0],
            truncation=8,
        )
        reflected = set(scan.orders[scan.propagating[0]].tolist())
        transmitted = set(
            scan.transmission.orders[scan.transmission.propagating[0]].tolist()
        )
        assert reflected < transmitted
        assert (scan.transmission.efficiency[~scan.transmission.propagating] == 0).all()

    def test_a_truncation_short_of_a_transmitted_order_is_refused(self, library):
        """At 200 nm, orders 4 and 5 do not reflect but do transmit into
        glass, so a window that holds every reflected order (3 is enough) can
        still be too narrow."""
        illumination = Illumination.classical(alpha=10.0)
        with pytest.raises(ValueError, match="order 5 is transmitted into glass"):
            rcwa.solve(lamellar("glass"), illumination, [200.0], truncation=4)
        rcwa.solve(lamellar("glass"), illumination, [200.0], truncation=5)


class TestReciprocity:
    """``E_m(alpha) = E_m(beta_m)``. Not exact at finite truncation, because
    the reversed problem keeps a different window of harmonics, centred on
    its own incident order; the residual must fall as the window grows.

    Measured 2026-10-27 on this grating at 450 and 650 nm, worst pair:

    ==========  =======  =======  =======  =======
    truncation  10       20       30       60
    ==========  =======  =======  =======  =======
    TE          5.8e-5   4.7e-6   9.8e-7   6.4e-8
    TM          8.1e-3   2.2e-3   9.9e-4   2.1e-4
    conical     5.3e-3   1.3e-3   5.4e-4   8.1e-5
    ==========  =======  =======  =======  =======

    (conical: gamma = 60, unpolarized, which is the reciprocal quantity once
    TE and TM convert into each other.) TM carries the metal's corner
    singularities, as in the Li-rule test below, and so it is the slow one.
    """

    PROBLEM_KW = dict(depth=0.3, duty=0.4)

    @pytest.mark.parametrize(
        "polarization, gamma, tolerance",
        [("TE", 90.0, 1e-5), ("TM", 90.0, 2e-3), ("unpolarized", 60.0, 1e-3)],
        ids=["TE", "TM", "conical"],
    )
    def test_metallic_lamellar_at_truncation_30(
        self, library, polarization, gamma, tolerance
    ):
        report = check_reciprocity(
            rcwa,
            lamellar("metal"),
            Illumination(alpha_deg=20.0, gamma_deg=gamma, polarization=polarization),
            [450.0, 650.0],
            tolerance=tolerance,
            truncation=30,
        )
        assert report.passed, report
        assert report.pairs_tested >= 4

    def test_the_residual_falls_with_truncation(self, library):
        def worst(truncation):
            return check_reciprocity(
                rcwa,
                lamellar("metal"),
                Illumination.classical(alpha=20.0, polarization="TM"),
                [650.0],
                truncation=truncation,
            ).max_violation

        assert worst(60) < worst(20) / 5


class TestLisRuleIsInEffect:
    """The only test that can tell Li's rule from Laurent's: metallic TM.
    Everything above passes with either. Specular TM on a metallic lamellar
    grating, against Li's rule at truncation 100.

    Measured 2026-10-27 against a reference at 150 (specular 0.66966), error
    by truncation:

    ====  =======  =======
    M     Li       Laurent
    ====  =======  =======
    5     6.6e-3   1.9e-1
    10    7.9e-3   1.4e-1
    20    4.0e-3   6.7e-2
    40    1.5e-3   9.5e-3
    80    4.1e-4   6.2e-3
    ====  =======  =======

    Li's rule is not fast here either, since the field is singular at the
    metal's corners, but it is 6 to 30 times closer at every truncation, and
    Laurent's is off by a fifth at low truncation. Note Li's error is not
    monotone (5 to 10), which is why the convergence harness asks for three
    agreeing rungs."""

    @staticmethod
    def specular(truncation, factorization, polarization="TM"):
        solution = solve(
            Stack(1.0, (Slab(180.0, 1.0, ((0.6, 1.0, METAL),)),), METAL),
            alpha=20.0,
            polarizations=(polarization,),
            truncation=truncation,
            factorization=factorization,
        )
        return float(solution.reflected[0, truncation])

    @pytest.mark.parametrize("truncation", [5, 10, 20])
    def test_li_converges_where_laurent_does_not(self, metal_tm_reference, truncation):
        li = abs(self.specular(truncation, "li") - metal_tm_reference)
        laurent = abs(self.specular(truncation, "laurent") - metal_tm_reference)
        assert li < 1e-2
        assert laurent > 5e-2
        assert laurent > 5 * li

    def test_te_does_not_care(self):
        """TE involves only ``Eyy``, which is Laurent's under both, so the
        two must agree exactly: the switch reaches TM and nothing else."""
        for truncation in (5, 20):
            assert self.specular(truncation, "li", "TE") == pytest.approx(
                self.specular(truncation, "laurent", "TE"), abs=1e-13
            )


@pytest.fixture(scope="module")
def metal_tm_reference():
    return TestLisRuleIsInEffect.specular(100, "li")


class TestRefusals:
    """``rcwa.md`` section 7: refused rather than approximated."""

    ILLUMINATION = Illumination.classical(alpha=10.0, polarization="TE")

    def test_no_material_is_refused(self, library):
        problem = Problem(period=600.0, profile=Lamellar(depth_fraction=0.2))
        with pytest.raises(UnsupportedConfiguration, match="needs a material"):
            rcwa.solve(problem, self.ILLUMINATION, [500.0])

    def test_an_undercut_is_refused_before_anything_is_built(self, library):
        undercut = FromProfileData(t=(0.0, 0.55, 0.45, 1.0), y=(0.0, 0.4, 0.1, 0.0))
        problem = Problem(period=600.0, profile=undercut, coating="metal")
        with pytest.raises(UnsupportedConfiguration, match="undercut"):
            rcwa.solve(problem, self.ILLUMINATION, [500.0])

    def test_gain_is_refused(self, library):
        with pytest.raises(ValueError, match="gain"):
            rcwa.solve(lamellar("gain"), self.ILLUMINATION, [500.0])

    def test_a_truncation_short_of_a_propagating_order_is_refused(self, library):
        """At 200 nm on a 600 nm period, orders -2..3 propagate; truncation 2
        cannot hold order 3, and quietly dropping it would hand its power to
        nobody."""
        with pytest.raises(ValueError, match="order 3 propagates"):
            rcwa.solve(lamellar("metal"), self.ILLUMINATION, [200.0], truncation=2)
        rcwa.solve(lamellar("metal"), self.ILLUMINATION, [200.0], truncation=3)

    @pytest.mark.parametrize("value", [0, -1, 2.5, True])
    def test_knob_values_are_checked(self, library, value):
        with pytest.raises(ValueError):
            rcwa.solve(lamellar("glass"), self.ILLUMINATION, [500.0], truncation=value)
        with pytest.raises(ValueError):
            rcwa.solve(lamellar("glass"), self.ILLUMINATION, [500.0], slices=value)


class TestProvenance:
    def test_it_records_what_it_did(self, library):
        scan = rcwa.solve(
            lamellar("glass"), Illumination.classical(alpha=10.0), [500.0], truncation=5
        )
        notes = scan.provenance.notes
        assert scan.provenance.method == "rcwa"
        assert scan.provenance.truncation == 5
        assert notes["harmonics"] == 11
        assert notes["profile"] == "exact (lamellar)"
        assert notes["slices"] == 1
        assert notes["material"] == "glass"
        assert "Li" in notes["factorization"]
        assert scan.provenance.warnings == ()

    def test_a_staircase_in_tm_on_a_metal_is_warned(self, library):
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.2),
                          coating="metal")
        tm = rcwa.solve(problem, Illumination.classical(10.0, "TM"), [500.0],
                        truncation=5)
        assert any("staircase" in w for w in tm.provenance.warnings)
        assert tm.provenance.notes["slices"] == 10  # max(8, 2 * 5)
        te = rcwa.solve(problem, Illumination.classical(10.0, "TE"), [500.0],
                        truncation=5)
        assert not te.provenance.warnings

    def test_an_ignored_substrate_is_named(self, library):
        problem = Problem(period=600.0, profile=Lamellar(depth_fraction=0.2),
                          coating="metal", substrate="glass")
        scan = rcwa.solve(problem, Illumination.classical(10.0), [500.0], truncation=4)
        assert any("'glass' is not consulted" in w for w in scan.provenance.warnings)


class TestGrooveGeometry:
    def test_lamellar_is_one_exact_slab(self):
        groove = groove_geometry(Lamellar(depth_fraction=0.25, duty_cycle=0.3), 800.0, 50)
        assert groove.exact
        assert groove.slabs == ((200.0, ((0.7, 1.0, 0),)),)

    def test_blazed_walls_are_the_profiles_own_crossings(self):
        """Each slab's wall sits where the facet crosses the slab's mid-height,
        to rounding, not on the 2048-point grid ``slice_layers`` samples."""
        profile = Blazed(blaze_angle=30.0, antiblaze_angle=70.0)
        groove = groove_geometry(profile, 1000.0, 7)
        assert not groove.exact and len(groove.slabs) == 7
        thickness = profile.depth * 1000.0 / 7
        for k, (slab_thickness, intervals) in enumerate(groove.slabs):
            assert slab_thickness == pytest.approx(thickness)
            level = profile.depth * (1.0 - (k + 0.5) / 7)  # top slab first
            ((start, end, key),) = intervals
            assert key == 0
            # xi = 1 - t, so the interval's ends are 1 - (the two crossings)
            for xi in (start, end):
                assert profile.height(1.0 - xi) == pytest.approx(level, abs=1e-12)

    def test_the_staircase_area_converges_to_the_profiles(self):
        profile = Sinusoidal(depth_fraction=0.3)
        groove = groove_geometry(profile, 1.0, 64)
        area = sum(t * sum(b - a for a, b, _ in iv) for t, iv in groove.slabs)
        # The mean height of a raised cosine is half its depth.
        assert area == pytest.approx(0.15, rel=1e-3)


class TestConvergenceHarness:
    def test_a_dielectric_lamellar_reaches_a_plateau(self, library):
        """The knob the capabilities declare is one the harness can sweep to
        a plateau. TE on glass converges quickly, which is the point: this
        is the harness's contract with the backend, not a hard case."""
        report = check_convergence(
            rcwa,
            lamellar("glass"),
            Illumination.classical(alpha=20.0, polarization="TE"),
            [500.0],
            values=(4, 8, 16, 32, 64),
            tolerance=1e-4,
        )
        assert report.passed, report
        assert report.scan.provenance.converged is True


class TestAgreementWithTheIntegralSolver:
    """A different rigorous method on the same ``Problem``: the integral
    solver's finite-conductivity mode, corner-graded, which shares no code
    with this backend beyond the geometry helpers and the material lookup.

    Measured 2026-10-27 on a 600 nm lamellar grating, depth 0.3, duty 0.4,
    alpha = 20, 500 nm, largest difference over the propagating orders at
    (truncation, boundary_points):

    ===================  ==========  ===========
    case                 (40, 512)   (80, 1024)
    ===================  ==========  ===========
    glass TE in-plane    2.5e-7      4.8e-8
    glass TM in-plane    2.1e-6      5.2e-7
    glass TM gamma=60    1.7e-6      4.2e-7
    metal TE in-plane    4.0e-6      4.6e-7
    metal TE gamma=60    4.6e-4      1.7e-4
    metal TM in-plane    1.9e-3      7.8e-4
    metal TM gamma=60    1.5e-3      6.3e-4
    ===================  ==========  ===========

    The difference falls as both methods refine, which is what makes it a
    comparison of the physics and not of two discretisations. The tests run
    the cheaper column with tolerances set from it.
    """

    @pytest.mark.parametrize(
        "material, polarization, gamma, tolerance",
        [
            ("glass", "TE", 90.0, 1e-6),
            ("glass", "TM", 90.0, 5e-6),
            ("glass", "TM", 60.0, 5e-6),
            ("metal", "TE", 90.0, 1e-5),
            ("metal", "TE", 60.0, 1e-3),
            ("metal", "TM", 90.0, 3e-3),
            ("metal", "TM", 60.0, 3e-3),
        ],
    )
    def test_lamellar(self, library, material, polarization, gamma, tolerance):
        problem = lamellar(material)
        illumination = Illumination(
            alpha_deg=20.0, gamma_deg=gamma, polarization=polarization
        )
        ours = rcwa.solve(problem, illumination, [500.0], truncation=40)
        theirs = integral.solve(
            problem,
            illumination,
            [500.0],
            boundary_points=512,
            conductivity="tabulated",
            corner_grading=4,
        )
        assert np.array_equal(ours.orders, theirs.orders)
        assert np.array_equal(ours.propagating, theirs.propagating)
        assert np.abs(ours.efficiency - theirs.efficiency).max() < tolerance
        # Absorbed power, and the transmitted orders of a lossless substrate,
        # from two independent constructions.
        assert ours.absorption == pytest.approx(theirs.absorption, abs=5 * tolerance)
        if material == "glass":
            mine, other = ours.transmission, theirs.transmission
            assert np.array_equal(mine.orders, other.orders)
            assert np.array_equal(mine.propagating, other.propagating)
            assert np.abs(mine.efficiency - other.efficiency).max() < 5 * tolerance
        else:
            assert ours.transmission is None and theirs.transmission is None


class TestAnAsymmetricGroove:
    """A blazed groove pins which way ``t`` runs (``conventions.md`` section
    3): RCWA must place the facets as the integral solver does, or it solves
    the mirror image. The groove's own mirror is solved too, to show the
    comparison can tell the two apart.

    A staircase, in TE, where staircasing converges well. Measured 2026-10-27
    against the integral solver at 1024 points, 20/70 degree facets, 600 nm
    period, alpha = 10, 500 nm, largest difference over the orders:

    ========================  ======  ======
    (truncation, slices)      glass   metal
    ========================  ======  ======
    (20, 20)                  1.8e-4  1.1e-3
    (20, 40)                  4.6e-5  3.1e-4
    (20, 80)                  1.2e-5  8.9e-5
    (40, 160)                 3.2e-6  2.0e-5
    ========================  ======  ======

    Second order in the slice count, as the staircase's area error is.
    """

    ILLUMINATION = Illumination.classical(alpha=10.0, polarization="TE")

    def groove(self, blaze, antiblaze):
        return Problem(
            period=600.0,
            profile=Blazed(blaze_angle=blaze, antiblaze_angle=antiblaze),
            coating="metal",
        )

    def test_the_facets_face_the_way_the_integral_solver_puts_them(self, library):
        reference = integral.solve(
            self.groove(20.0, 70.0),
            self.ILLUMINATION,
            [500.0],
            boundary_points=512,
            conductivity="tabulated",
            corner_grading=4,
        ).efficiency
        ours = rcwa.solve(
            self.groove(20.0, 70.0), self.ILLUMINATION, [500.0], truncation=20, slices=40
        ).efficiency
        mirror = rcwa.solve(
            self.groove(70.0, 20.0), self.ILLUMINATION, [500.0], truncation=20, slices=40
        ).efficiency
        # Measured: 3.1e-4 from the right groove, 1.35e-2 from its mirror.
        agreement = np.abs(ours - reference).max()
        assert agreement < 1e-3
        assert np.abs(mirror - reference).max() > 10 * agreement


# ---------------------------------------------------------------------------
# Conformal layers (Problem.layers)
# ---------------------------------------------------------------------------


@pytest.fixture
def layer_library(library):
    for name, index in (("vacuum", 1.0 + 0.0j), ("film", 1.38 + 0.0j)):
        write_ari(_table(name, index), library / f"{name}.ari")
    return library


def coated_lamellar(depth, duty, tau, arc=200):
    """The exact true offset of a lamellar groove by ``tau``, all in periods.

    Written out by hand, not through the backend's envelope: ridge top raised
    by ``tau``, walls moved out by ``tau`` and kept vertical, the two convex
    corners rounded into quarter circles of radius ``tau``, the two concave
    ones at the floor left sharp. ``t`` repeats at a vertical wall, which
    ``FromProfileData`` allows.
    """
    theta = np.linspace(0.0, np.pi / 2, arc)
    t = [0.0, duty]
    y = [depth + tau, depth + tau]
    t += list(duty + tau * np.sin(theta[1:]))
    y += list(depth + tau * np.cos(theta[1:]))
    t += [duty + tau, 1.0 - tau, 1.0 - tau]
    y += [tau, tau, depth]
    t += list(1.0 - tau * np.cos(theta[1:]))
    y += list(depth + tau * np.sin(theta[1:]))
    return FromProfileData(t=tuple(t), y=tuple(np.array(y) - tau))


def coated(material, layers, depth=0.3, duty=0.4):
    return Problem(
        period=600.0,
        profile=Lamellar(depth_fraction=depth, duty_cycle=duty),
        coating=material,
        layers=tuple(ConformalLayer(material=m, thickness=t) for m, t in layers),
    )


class TestConformalLayers:
    """``Problem.layers``: a layer following the groove, its outer surface
    the true offset of the one beneath (``ConformalLayer``)."""

    def test_the_offset_is_the_rounded_one(self):
        """On a lamellar ridge: the top raised by T, the floor raised by T,
        and a convex corner rounded into a circle of radius T about it."""
        period, depth, duty, thickness = 600.0, 0.3, 0.4, 40.0
        x, y = _surface(Lamellar(depth_fraction=depth, duty_cycle=duty), period)
        height = _offset_height(x, y, period, thickness)
        # The ridge spans xi = 0.6..1.0 (xi = 1 - t), its corners at 0.6 and 1.
        assert height([0.8]) == pytest.approx([depth * period + thickness])
        assert height([0.3]) == pytest.approx([thickness])
        for dx in (5.0, 20.0, 35.0):
            xi = 0.6 - dx / period
            expected = depth * period + np.sqrt(thickness**2 - dx**2)
            assert height([xi]) == pytest.approx([expected], abs=1e-9)
        # The wall's face moved out by T and kept vertical: a nanometre
        # beyond it the floor's offset, a nanometre inside it the arc.
        face = 0.6 - thickness / period
        assert height([face - 1.0 / period]) == pytest.approx([thickness])
        assert height([face + 1.0 / period]) == pytest.approx(
            [depth * period + np.sqrt(thickness**2 - (thickness - 1.0) ** 2)]
        )

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_a_vacuum_layer_changes_nothing(self, layer_library, polarization):
        """The bookkeeping test: a layer of the superstrate's own index must
        reproduce the bare grating, through a geometry of many more slabs.
        Measured 7e-14."""
        illumination = Illumination.classical(alpha=20.0, polarization=polarization)
        bare = rcwa.solve(lamellar("metal"), illumination, [500.0], truncation=20)
        vacuum = rcwa.solve(
            coated("metal", [("vacuum", 40.0)]), illumination, [500.0], truncation=20
        )
        assert vacuum.provenance.notes["slices"] > 1
        assert np.abs(bare.efficiency - vacuum.efficiency).max() < 1e-12

    def test_a_film_on_a_vanishing_groove_is_the_airy_film(self, layer_library):
        scan = rcwa.solve(
            coated("glass", [("film", 100.0)], depth=1e-4),
            Illumination.classical(alpha=20.0, polarization="TE"),
            [500.0],
            truncation=10,
        )
        rs, _ = airy(1.38, 100.0, GLASS, graze_of(20.0, 90.0), 500.0)
        assert scan.efficiency[0, scan.orders == 0][0] == pytest.approx(rs, abs=1e-4)

    @pytest.mark.parametrize(
        "material, polarization, tolerance",
        [("glass", "TE", 2e-4), ("glass", "TM", 2e-4), ("metal", "TE", 5e-4)],
    )
    def test_a_layer_of_the_gratings_own_material_is_the_offset_groove(
        self, layer_library, material, polarization, tolerance
    ):
        """Coat a grating in its own material and it is a bare grating of the
        offset shape, which the integral solver can solve exactly as a curve.
        Independent of the envelope code: the offset groove is written out by
        hand above.

        Measured 2026-10-27 against the integral solver at 1024 points, 40 nm
        layer, at (truncation, slices) = (20, 40): 2.1e-5 glass TE, 3.1e-5
        glass TM, 1.6e-4 metal TE, against the integral solver's own change
        from 512 to 1024 points of 3e-5 (glass) and 7e-5 (metal TE). So the
        agreement is at the reference's accuracy. Metal TM is the test below.
        """
        illumination = Illumination.classical(alpha=20.0, polarization=polarization)
        ours = rcwa.solve(
            coated(material, [(material, 40.0)]),
            illumination,
            [500.0],
            truncation=20,
            slices=40,
        )
        theirs = integral.solve(
            Problem(period=600.0, profile=coated_lamellar(0.3, 0.4, 40.0 / 600.0),
                    coating=material),
            illumination,
            [500.0],
            boundary_points=512,
            conductivity="tabulated",
            corner_grading=4,
        )
        assert np.abs(ours.efficiency - theirs.efficiency).max() < tolerance

    def test_metallic_tm_on_rounded_corners_converges_slowly(self, layer_library):
        """The staircase weakness (``rcwa.md`` section 5), met where it is
        worst: a metal in TM with rounded corners. Measured 2026-10-27 against
        the integral solver at 1024 points (specular 0.74848): 4.2e-2 at
        (10, 20), 3.3e-2 at (20, 40), 2.8e-2 at (20, 80), 1.1e-2 at (40, 160).
        Recorded so the weakness is visible, and so that it moves the right
        way; this is not a validation of the case."""
        illumination = Illumination.classical(alpha=20.0, polarization="TM")
        reference = 0.74848

        def specular(truncation, slices):
            scan = rcwa.solve(
                coated("metal", [("metal", 40.0)]), illumination, [500.0],
                truncation=truncation, slices=slices,
            )
            return float(scan.efficiency[0, scan.orders == 0][0])

        coarse, finer = abs(specular(10, 20) - reference), abs(specular(20, 80) - reference)
        assert 1e-2 < finer < coarse
        assert any("staircase" in w for w in rcwa.solve(
            coated("metal", [("metal", 40.0)]), illumination, [500.0], truncation=4
        ).provenance.warnings)

    @pytest.mark.parametrize("gamma", [90.0, 50.0])
    def test_lossless_layers_conserve_energy(self, layer_library, gamma):
        scan = rcwa.solve(
            coated("glass", [("film", 30.0), ("glass", 20.0)]),
            Illumination(alpha_deg=20.0, gamma_deg=gamma, polarization="unpolarized"),
            [450.0, 600.0],
            truncation=8,
        )
        assert check_energy_balance(scan, lossless=True, tolerance=1e-10).passed
        assert np.abs(scan.absorption).max() < 1e-12
        assert scan.provenance.notes["layers"] == ["film 30 nm", "glass 20 nm"]

    def test_a_metal_layer_on_glass_absorbs_and_transmits(self, layer_library):
        scan = rcwa.solve(
            coated("glass", [("metal", 10.0)]),
            Illumination.classical(alpha=20.0),
            [500.0],
            truncation=8,
        )
        assert check_energy_balance(scan, lossless=True, tolerance=1e-10).passed
        assert scan.absorption[0] > 1e-3
        assert scan.transmission.total[0] > 1e-3

    def test_a_layer_with_gain_is_refused(self, layer_library):
        with pytest.raises(ValueError, match="'gain'"):
            rcwa.solve(coated("glass", [("gain", 10.0)]), Illumination.classical(10.0),
                       [500.0], truncation=4)

    @pytest.mark.parametrize("name", ["integral", "scalar"])
    def test_solvers_that_cannot_model_layers_refuse_them(self, layer_library, name):
        solver = get_solver(name)
        assert not solver.capabilities.handles_layers
        with pytest.raises(UnsupportedConfiguration, match="conformal layers"):
            solver.capabilities.check(
                coated("metal", [("film", 10.0)]), Illumination.classical(10.0)
            )

    def test_a_problem_with_layers_round_trips_through_json(self):
        problem = coated("metal", [("film", 10.0), ("glass", 2.5)])
        assert Problem.model_validate_json(problem.model_dump_json()) == problem
