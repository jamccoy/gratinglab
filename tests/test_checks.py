r"""Physics self-checks.

These validate the model rather than the arithmetic, so the important tests are
the ones showing a check *fails* when the physics is wrong. A check that always
passes is worse than no check.
"""

import numpy as np
import pytest

from gratinglab.checks import check_energy_balance, check_reciprocity
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, Lamellar, Sinusoidal
from gratinglab.result import EfficiencyScan, Provenance
from gratinglab.solvers import scalar

UNPOL = "unpolarized"

GEOMETRIES = [
    pytest.param(315.15, 25.0, 1.5, [3.0], id="off-plane-xray"),
    pytest.param(1400.0, 10.0, 90.0, [600.0], id="in-plane-visible"),
    pytest.param(1400.0, -30.0, 45.0, [550.0], id="general-conical"),
]


class TestReciprocity:
    @pytest.mark.parametrize("period,alpha,gamma,wavelengths", GEOMETRIES)
    @pytest.mark.parametrize(
        "profile",
        [
            Blazed(blaze_angle=29.5, antiblaze_angle=70.5),
            Lamellar(depth_fraction=0.2, duty_cycle=0.4),
            Sinusoidal(depth_fraction=0.15),
        ],
        ids=["blazed", "lamellar", "sinusoid"],
    )
    def test_scalar_solver_is_reciprocal(
        self, period, alpha, gamma, wavelengths, profile
    ):
        r"""E_m(alpha) == E_m(beta_m), to machine precision."""
        report = check_reciprocity(
            scalar,
            Problem(period=period, profile=profile),
            Illumination(alpha_deg=alpha, gamma_deg=gamma, polarization=UNPOL),
            wavelengths,
            quadrature_points=4096,
        )
        assert report.pairs_tested > 0
        assert report.passed, str(report)
        assert report.max_violation < 1e-12

    def test_detects_a_phase_function_that_is_not_symmetric(self, monkeypatch):
        r"""The check must FAIL for wrong physics, or it is worthless.

        Replacing :math:`\cos\alpha + \cos\beta_m` with :math:`2\cos\alpha`
        removes the exit-direction dependence. Every closed-form test in
        ``test_scalar.py`` for a symmetric geometry would still pass, because
        those compare against a formula derived the same way. Reciprocity does
        not, because it constrains structure rather than value.
        """
        original = scalar.solve

        def asymmetric(problem, illumination, wavelengths, **options):
            """Solve with alpha substituted for beta_m in the phase."""
            scan = original(problem, illumination, wavelengths, **options)
            # Emulate a phase that ignores the exit direction by re-solving at
            # normal-ish incidence: the point is only that it breaks symmetry.
            twin = Illumination(
                alpha_deg=0.0,
                gamma_deg=illumination.gamma_deg,
                polarization=illumination.polarization,
            )
            other = original(problem, twin, wavelengths, **options)
            return EfficiencyScan(
                wavelengths=scan.wavelengths,
                orders=scan.orders,
                efficiency=np.where(scan.propagating, other.efficiency, 0.0)
                if other.efficiency.shape == scan.efficiency.shape
                else scan.efficiency,
                propagating=scan.propagating,
                provenance=scan.provenance,
            )

        class Broken:
            capabilities = scalar.capabilities
            solve = staticmethod(asymmetric)

        report = check_reciprocity(
            Broken(),
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)),
            Illumination(alpha_deg=25.0, gamma_deg=1.5, polarization=UNPOL),
            [3.0],
            quadrature_points=2048,
        )
        assert not report.passed, "reciprocity accepted a non-reciprocal solver"

    def test_report_is_usable_when_the_solver_is_exactly_reciprocal(self):
        """A zero violation must still name where it was measured.

        Regression: the running maximum was seeded at 0.0 and only updated on
        a strictly greater violation, so a solver that was *bitwise* reciprocal
        never set worst_order. It passed locally, where floating-point noise
        gave ~1e-18, and failed on CI, where the BLAS returned exact zero --
        the better the solver, the worse the report.
        """

        class PerfectlyReciprocal:
            capabilities = scalar.capabilities

            @staticmethod
            def solve(problem, illumination, wavelengths, **options):
                """Ignore the illumination azimuth, so reciprocity is exact."""
                fixed = Illumination(
                    alpha_deg=0.0,
                    gamma_deg=illumination.gamma_deg,
                    polarization=illumination.polarization,
                )
                return scalar.solve(problem, fixed, wavelengths, **options)

        report = check_reciprocity(
            PerfectlyReciprocal(),
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5)),
            Illumination(alpha_deg=25.0, gamma_deg=1.5, polarization=UNPOL),
            [3.0],
            quadrature_points=2048,
        )
        assert report.max_violation == 0.0
        assert report.passed
        assert report.pairs_tested > 0
        assert report.worst_order is not None, "zero violation left the report unusable"
        assert report.worst_wavelength == 3.0

    def test_reports_where_the_worst_violation_is(self):
        report = check_reciprocity(
            scalar,
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5)),
            Illumination(alpha_deg=25.0, gamma_deg=1.5, polarization=UNPOL),
            [3.0],
            quadrature_points=2048,
        )
        assert report.worst_order is not None
        assert report.worst_wavelength == 3.0
        assert "reciprocity" in str(report)

    def test_max_orders_caps_the_cost(self):
        common = dict(
            problem=Problem(period=1400.0, profile=Blazed(blaze_angle=30.0)),
            illumination=Illumination.classical(alpha=10.0, polarization=UNPOL),
            wavelengths=[600.0],
            quadrature_points=1024,
        )
        few = check_reciprocity(scalar, max_orders=2, **common)
        many = check_reciprocity(scalar, max_orders=6, **common)
        assert few.pairs_tested < many.pairs_tested

    def test_no_testable_pairs_is_not_a_pass(self):
        """An empty check must not report success."""
        report = check_reciprocity(
            scalar,
            Problem(period=1400.0, profile=Blazed(blaze_angle=30.0)),
            Illumination.classical(alpha=10.0, polarization=UNPOL),
            [600.0],
            max_orders=0,
        )
        assert report.pairs_tested == 0
        assert not report.passed


def make_scan(total_per_wavelength):
    """A scan whose orders sum to the given totals."""
    values = np.asarray(total_per_wavelength, dtype=float)
    return EfficiencyScan(
        wavelengths=np.arange(1.0, len(values) + 1.0),
        orders=np.array([0, 1]),
        efficiency=np.column_stack([values * 0.6, values * 0.4]),
        propagating=np.ones((len(values), 2), dtype=bool),
        provenance=Provenance("test"),
    )


class TestEnergyBalance:
    def test_accepts_a_deficit(self):
        """Absorption and evanescent leakage are ordinary."""
        report = check_energy_balance(make_scan([0.4, 0.7, 0.95]))
        assert report.passed
        assert report.max_deficit == pytest.approx(0.6)

    def test_rejects_an_excess(self):
        """No passive grating returns more power than it receives."""
        report = check_energy_balance(make_scan([0.9, 1.5, 0.8]))
        assert not report.passed
        assert report.max_excess == pytest.approx(0.5)
        assert report.unphysical.tolist() == [False, True, False]

    def test_lossless_mode_requires_exact_unity(self):
        near_unity = make_scan([1.0, 0.9999999, 1.0])
        assert check_energy_balance(near_unity, lossless=True).passed
        assert not check_energy_balance(make_scan([1.0, 0.5]), lossless=True).passed

    def test_default_mode_tolerates_what_lossless_mode_does_not(self):
        scan = make_scan([0.5, 0.6])
        assert check_energy_balance(scan).passed
        assert not check_energy_balance(scan, lossless=True).passed

    def test_the_scalar_solver_does_not_conserve_energy(self):
        r"""Recorded, not fixed -- and deliberately so.

        The phase carries :math:`\cos\beta_m`, so the coefficients are not a
        Parseval pair and the sum drifts from unity. That is the price of
        keeping the phase symmetric under :math:`\alpha \leftrightarrow
        \beta_m`, which is what makes the solver *reciprocal*. The
        energy-conserving alternative violates reciprocity instead; no scalar
        formulation satisfies both.

        Asserting the drift *exists* stops it being quietly "fixed" by
        renormalising, which would erase a real limitation of scalar theory.
        """
        scan = scalar.solve(
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5)),
            Illumination(alpha_deg=25.0, gamma_deg=1.5, polarization=UNPOL),
            np.linspace(1.0, 5.0, 15),
            quadrature_points=4096,
        )
        report = check_energy_balance(scan)
        assert not report.passed or abs(report.total - 1.0).max() > 0.01

    def test_the_violation_is_surfaced_as_a_warning(self):
        """A silent violation would be worse than the violation itself."""
        scan = scalar.solve(
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5)),
            Illumination(alpha_deg=25.0, gamma_deg=1.5, polarization=UNPOL),
            np.linspace(1.0, 5.0, 15),
            quadrature_points=4096,
        )
        assert any("summed efficiency" in w for w in scan.provenance.warnings)

    def test_reports_the_range(self):
        report = check_energy_balance(make_scan([0.3, 0.9]))
        assert "energy balance" in str(report)
        assert report.total.tolist() == pytest.approx([0.3, 0.9])
