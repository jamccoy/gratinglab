"""The finite-conductivity transverse core, against every closed-form answer.

The decisive references, in order of strictness: the flat interface (exact
Fresnel amplitudes, including phase -- the test that pins every sign in the
coupled system, successor to the flat-mirror pins), the two energy theorems
(R + T = 1 for a lossless substrate, R + A = 1 with the absorption computed
by an independent boundary integral), the perfectly conducting limit at
n = 100i against the already-validated milestone-1 solver, and the
shallow-groove limit back to Fresnel.
"""

import numpy as np
import pytest

from gratinglab.materials import fresnel
from gratinglab.problem import Problem
from gratinglab.profiles import Sinusoidal
from gratinglab.solvers.integral._boundary import PhysicalBoundary, physical_boundary
from gratinglab.solvers.integral._core import solve_transverse
from gratinglab.solvers.integral._finite import solve_transverse_finite

GOLD_VISIBLE = 0.2 + 3.0j
XRAY = 1.0 - 1e-3 + 1e-4j
GLASS = 1.5 + 0.0j


def flat_mirror(period=600.0, points=128):
    x = np.linspace(0.0, period, points, endpoint=False)
    return PhysicalBoundary(
        x=x,
        y=np.zeros(points),
        nx=np.zeros(points),
        ny=np.ones(points),
        arc_length=period,
    )


def specular(solution, amplitudes):
    return amplitudes[solution.orders == 0][0]


class TestFlatInterfaceIsFresnel:
    """A flat boundary must reproduce the Fresnel amplitudes of
    ``materials.fresnel`` exactly -- magnitude and phase, TE against r_s and
    TM against r_p (the in-plane identification; the H-form r_p convention
    matches, as its perfectly conducting limit +1 matches the flat-mirror
    TM pin)."""

    @pytest.mark.parametrize(
        "index", [GOLD_VISIBLE, XRAY, GLASS], ids=["gold", "xray", "glass"]
    )
    @pytest.mark.parametrize("sin_alpha", [0.0, 0.5])
    def test_te_amplitude_is_r_s(self, index, sin_alpha):
        solution = solve_transverse_finite(
            flat_mirror(),
            wavelength=550.0,
            period=600.0,
            sin_alpha=sin_alpha,
            index=index,
            incident=(1.0, 0.0),
            terms=64,
        )
        graze = np.pi / 2.0 - np.arcsin(sin_alpha)
        r_s, _ = fresnel.amplitude(index, graze)
        assert specular(solution, solution.e_amplitudes) == pytest.approx(
            complex(r_s), abs=1e-4
        )
        assert np.allclose(solution.b_amplitudes, 0.0, atol=1e-10)

    @pytest.mark.parametrize(
        "index", [GOLD_VISIBLE, XRAY, GLASS], ids=["gold", "xray", "glass"]
    )
    @pytest.mark.parametrize("sin_alpha", [0.0, 0.5])
    def test_tm_amplitude_is_r_p(self, index, sin_alpha):
        solution = solve_transverse_finite(
            flat_mirror(),
            wavelength=550.0,
            period=600.0,
            sin_alpha=sin_alpha,
            index=index,
            incident=(0.0, 1.0),
            terms=64,
        )
        graze = np.pi / 2.0 - np.arcsin(sin_alpha)
        _, r_p = fresnel.amplitude(index, graze)
        assert specular(solution, solution.b_amplitudes) == pytest.approx(
            complex(r_p), abs=1e-4
        )
        assert np.allclose(solution.e_amplitudes, 0.0, atol=1e-10)

    @pytest.mark.parametrize("incident", [(1.0, 0.0), (0.0, 1.0)], ids=["TE", "TM"])
    def test_flat_energy_balance_with_independent_absorption(self, incident):
        """On the flat interface the absorption integral reduces to
        1 - |r|^2 in closed form; the discrete version must agree."""
        solution = solve_transverse_finite(
            flat_mirror(),
            wavelength=550.0,
            period=600.0,
            sin_alpha=0.5,
            index=GOLD_VISIBLE,
            incident=incident,
            terms=64,
        )
        assert solution.total + solution.absorption == pytest.approx(1.0, abs=1e-4)
        assert solution.absorption > 0.0


class TestEnergyTheorems:
    """R + T = 1 (lossless substrate) and R + A = 1 (absorbing substrate)
    are theorems in the continuum; the absorption side is an independent
    boundary integral of the densities, so the check is two-sided."""

    SMOOTH = physical_boundary(
        Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.5)), 200
    )

    @pytest.mark.parametrize("incident", [(1.0, 0.0), (0.0, 1.0)], ids=["TE", "TM"])
    def test_dielectric_reflection_plus_transmission(self, incident):
        solution = solve_transverse_finite(
            self.SMOOTH,
            wavelength=500.0,
            period=600.0,
            sin_alpha=0.35,
            index=GLASS,
            incident=incident,
            terms=90,
        )
        assert len(solution.transmitted_orders) > len(solution.orders)
        assert solution.total + solution.transmitted_total == pytest.approx(
            1.0, abs=5e-4
        )

    @pytest.mark.parametrize("incident", [(1.0, 0.0), (0.0, 1.0)], ids=["TE", "TM"])
    def test_dielectric_boundary_integral_is_the_transmitted_power(self, incident):
        """For a lossless substrate the 'absorption' integral is the power
        crossing the interface, so it must equal the transmitted Rayleigh
        sum -- two entirely different post-processings of one solve."""
        solution = solve_transverse_finite(
            self.SMOOTH,
            wavelength=500.0,
            period=600.0,
            sin_alpha=0.35,
            index=GLASS,
            incident=incident,
            terms=90,
        )
        assert solution.absorption == pytest.approx(
            solution.transmitted_total, abs=2e-4
        )

    @pytest.mark.parametrize("incident", [(1.0, 0.0), (0.0, 1.0)], ids=["TE", "TM"])
    def test_metal_reflection_plus_absorption(self, incident):
        solution = solve_transverse_finite(
            self.SMOOTH,
            wavelength=500.0,
            period=600.0,
            sin_alpha=0.35,
            index=GOLD_VISIBLE,
            incident=incident,
            terms=90,
        )
        assert len(solution.transmitted_orders) == 0
        assert solution.absorption > 0.0
        assert solution.total + solution.absorption == pytest.approx(1.0, abs=1e-3)


class TestPerfectlyConductingLimit:
    """At n = 100i the transmission problem degenerates to the milestone-1
    equations (V^- and c_B vanish like 1/|n|), so the finite-conductivity
    solve must land on the validated perfectly conducting solver -- up to
    the physics it is deliberately adding: the field penetrates
    ~lambda/(2 pi |n|) ~ 1 nm into the metal, which against a 180 nm groove
    reads as a ~0.5% effective-depth change. Measured residual at these
    settings: 1.3e-3 (TE) and 7.4e-4 (TM) worst order, first order in
    1/|n|; a purely imaginary index is purely reactive, so absorption
    stays at numerical zero (measured ~4e-6)."""

    PROBLEM = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3))

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_sinusoid_matches_milestone_one(self, polarization):
        boundary = physical_boundary(self.PROBLEM, 320)
        reference = solve_transverse(
            boundary,
            wavelength=600.0,
            period=600.0,
            sin_alpha=0.5,
            polarization=polarization,
            terms=160,
        )
        incident = (1.0, 0.0) if polarization == "TE" else (0.0, 1.0)
        solution = solve_transverse_finite(
            boundary,
            wavelength=600.0,
            period=600.0,
            sin_alpha=0.5,
            index=100.0j,
            incident=incident,
            terms=160,
        )
        assert np.array_equal(solution.orders, reference.orders)
        assert np.allclose(
            solution.efficiencies, reference.efficiencies, atol=3e-3
        )
        assert abs(solution.absorption) < 1e-4
        assert solution.total + solution.absorption == pytest.approx(
            1.0, abs=1e-4
        )

    def test_absorption_falls_along_the_lossy_path_to_the_conductor(self):
        """Along n = |n| e^{i pi/4} -- a lossy path to the perfect
        conductor -- the absorbed fraction falls monotonically toward zero
        while R + A stays pinned at 1 (measured: 0.603, 0.251, 0.092 with
        R + A - 1 under 1e-6 throughout)."""
        boundary = physical_boundary(self.PROBLEM, 240)
        absorptions = []
        for magnitude in (3.0, 10.0, 30.0):
            solution = solve_transverse_finite(
                boundary,
                wavelength=600.0,
                period=600.0,
                sin_alpha=0.5,
                index=magnitude * np.exp(1j * np.pi / 4.0),
                incident=(1.0, 0.0),
                terms=120,
            )
            absorptions.append(solution.absorption)
            assert solution.total + solution.absorption == pytest.approx(
                1.0, abs=1e-5
            )
        assert all(a > 0.0 for a in absorptions)
        assert absorptions[0] > absorptions[1] > absorptions[2]


class TestConicalFlatInterface:
    """The conical flat interface still has a closed form: rotate into the
    local s/p basis, apply Fresnel, rotate back. With the incident direction
    d = (-sin(gamma) sin(alpha), -sin(gamma) cos(alpha), cos(gamma)) and
    s-hat = (d x y-hat)/|.|, the (E_z, B_z) reflection matrix is
    M_r diag(r_s, r_p) M_i^{-1} with M_i = [[a, b], [-b, a]],
    M_r = [[a, -b], [b, a]], a = d_x/rho, b = d_z (gamma_0/k)/rho. This is
    the test that pins the cross blocks and every sign in the coupled
    system against physics the solver does not share code with
    (measured agreement 6e-5 at these settings, phases included)."""

    SIN_ALPHA = 0.3
    COS_GAMMA = 0.8
    INDEX = GOLD_VISIBLE

    def analytic_matrix(self):
        sin_gamma = np.sqrt(1.0 - self.COS_GAMMA**2)
        cos_alpha = np.sqrt(1.0 - self.SIN_ALPHA**2)
        d_x, d_z = -sin_gamma * self.SIN_ALPHA, self.COS_GAMMA
        gamma0_over_k = sin_gamma * cos_alpha
        rho = np.hypot(d_x, d_z)
        a, b = d_x / rho, d_z * gamma0_over_k / rho
        m_in = np.array([[a, b], [-b, a]])
        m_out = np.array([[a, -b], [b, a]])
        r_s, r_p = fresnel.amplitude(self.INDEX, np.arcsin(gamma0_over_k))
        return m_out @ np.diag([complex(r_s), complex(r_p)]) @ np.linalg.inv(m_in)

    def solver_matrix(self):
        columns = []
        for incident in [(1.0, 0.0), (0.0, 1.0)]:
            solution = solve_transverse_finite(
                flat_mirror(),
                wavelength=550.0,
                period=600.0,
                sin_alpha=self.SIN_ALPHA,
                index=self.INDEX,
                cos_gamma=self.COS_GAMMA,
                incident=incident,
                terms=64,
            )
            assert solution.total + solution.absorption == pytest.approx(
                1.0, abs=1e-4
            )
            columns.append(
                [
                    specular(solution, solution.e_amplitudes),
                    specular(solution, solution.b_amplitudes),
                ]
            )
        return np.array(columns).T

    def test_reflection_matrix_is_rotated_fresnel(self):
        got = self.solver_matrix()
        expected = self.analytic_matrix()
        assert np.allclose(got, expected, rtol=0.0, atol=1e-4)
        # Polarization conversion is real and antisymmetric between the two
        # incident states -- the structure the cross blocks impose.
        assert got[0, 1] == pytest.approx(-got[1, 0], abs=1e-12)
        assert abs(got[0, 1]) > 0.1


class TestNearInPlaneLimit:
    def test_conical_path_joins_the_decoupled_path(self):
        """cos(gamma) = 1e-8 runs the full block machinery; its answer must
        join the decoupled in-plane path continuously."""
        boundary = physical_boundary(
            Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3)), 160
        )
        kwargs = dict(
            wavelength=500.0, period=600.0, sin_alpha=0.35,
            index=GOLD_VISIBLE, incident=(1.0, 0.0), terms=80,
        )
        in_plane = solve_transverse_finite(boundary, cos_gamma=0.0, **kwargs)
        nearly = solve_transverse_finite(boundary, cos_gamma=1e-8, **kwargs)
        assert np.allclose(
            nearly.efficiencies, in_plane.efficiencies, atol=1e-8
        )
        assert nearly.absorption == pytest.approx(
            in_plane.absorption, abs=1e-8
        )


class TestGorayShmidtTable3:
    """The external anchor: G&S Table 3 (dielectric sine, conical, pure
    E_z incidence -- their 'B_z = 0'), itself a cross-method comparison
    against Li's coordinate transformation method. Their theta = 60 deg,
    phi = 15 deg map to sin_alpha = sin 60, cos_gamma = sin 15, and their
    order n is this project's m = -n (opposite Bloch-phase sign). Measured
    agreement at these settings: within 6e-5 of every tabulated order,
    reflected and transmitted (the table itself is quoted to ~1e-4)."""

    # G&S Table 3 / Li's CM, in fractional units, keyed by our order m.
    REFLECTED = {0: 0.1033, 1: 0.03873, 2: 0.03741, 3: 0.01121}
    TRANSMITTED = {
        -2: 0.06351, -1: 0.5183, 0: 0.07146, 1: 0.09925,
        2: 0.04922, 3: 0.007396, 4: 0.00002466, 5: 0.0001858,
    }

    BOUNDARY = physical_boundary(
        Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3)), 200
    )
    KWARGS = dict(
        wavelength=300.0,
        period=600.0,
        sin_alpha=float(np.sin(np.radians(60.0))),
        cos_gamma=float(np.sin(np.radians(15.0))),
        index=2.0 + 0.0j,
        terms=90,
    )

    def test_reflected_and_transmitted_orders(self):
        solution = solve_transverse_finite(
            self.BOUNDARY, incident=(1.0, 0.0), **self.KWARGS
        )
        got_r = dict(zip(solution.orders.tolist(), solution.efficiencies))
        for m, expected in self.REFLECTED.items():
            assert got_r[m] == pytest.approx(expected, abs=2e-3), f"R m={m}"
        got_t = dict(
            zip(
                solution.transmitted_orders.tolist(),
                solution.transmitted_efficiencies,
            )
        )
        assert len(got_t) == len(self.TRANSMITTED)
        for m, expected in self.TRANSMITTED.items():
            assert got_t[m] == pytest.approx(expected, abs=2e-3), f"T m={m}"
        assert solution.total + solution.transmitted_total == pytest.approx(
            1.0, abs=2e-4
        )
        # Conical incidence converts polarization: the reflected field of a
        # pure-E_z wave carries genuine B_z content.
        assert np.abs(solution.b_amplitudes).max() > 0.01

    def test_the_printed_table_4_is_a_duplication(self):
        """G&S Table 4 ('E_z = 0') prints exactly Table 3's efficiencies,
        which our solve shows cannot be right: the true E_z = 0 case
        differs decisively (specular 1.07% against the printed 10.33%)
        while still conserving energy to 2e-5. Recorded in
        docs/findings.md; the assertion here pins the disagreement so the
        finding survives refactors."""
        solution = solve_transverse_finite(
            self.BOUNDARY, incident=(0.0, 1.0), **self.KWARGS
        )
        assert solution.total + solution.transmitted_total == pytest.approx(
            1.0, abs=2e-4
        )
        specular_eff = solution.efficiencies[solution.orders == 0][0]
        assert abs(specular_eff - self.REFLECTED[0]) > 0.05


class TestPerfectlyConductingLimitConical:
    """The conical PC limit: at n = 100i the coupled system must land on
    the milestone-1 conical reduction (one in-plane solve at the reduced
    wavelength), cross-polarization conversion and absorption both dying
    with 1/|n|. Same penetration-depth residual as the in-plane limit."""

    PROBLEM = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3))
    COS_GAMMA = 0.5

    @pytest.mark.parametrize("polarization", ["TE", "TM"])
    def test_sinusoid_matches_the_conical_reduction(self, polarization):
        boundary = physical_boundary(self.PROBLEM, 320)
        sin_gamma = float(np.sqrt(1.0 - self.COS_GAMMA**2))
        reference = solve_transverse(
            boundary,
            wavelength=600.0 / sin_gamma,
            period=600.0,
            sin_alpha=0.5,
            polarization=polarization,
            terms=160,
        )
        incident = (1.0, 0.0) if polarization == "TE" else (0.0, 1.0)
        solution = solve_transverse_finite(
            boundary,
            wavelength=600.0,
            period=600.0,
            sin_alpha=0.5,
            index=100.0j,
            cos_gamma=self.COS_GAMMA,
            incident=incident,
            terms=160,
        )
        assert np.array_equal(solution.orders, reference.orders)
        assert np.allclose(
            solution.efficiencies, reference.efficiencies, atol=3e-3
        )
        assert solution.total + solution.absorption == pytest.approx(
            1.0, abs=1e-4
        )


class TestShallowGrooveLimit:
    """A vanishing groove is a flat interface: the specular order carries
    the Fresnel reflectivity and the diffracted orders carry nothing."""

    def test_specular_approaches_fresnel(self):
        boundary = physical_boundary(
            Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.005)), 200
        )
        solution = solve_transverse_finite(
            boundary,
            wavelength=500.0,
            period=600.0,
            sin_alpha=0.35,
            index=GOLD_VISIBLE,
            incident=(1.0, 0.0),
            terms=90,
        )
        graze = np.pi / 2.0 - np.arcsin(0.35)
        expected = fresnel.reflectivity(GOLD_VISIBLE, graze, polarization="s")
        got = solution.efficiencies[solution.orders == 0][0]
        assert got == pytest.approx(float(expected), abs=1e-3)
        assert solution.total - got < 1e-3
