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


class RecordingSolver:
    """A reciprocal solver whose per-order efficiency we dictate, and which
    remembers every illumination it was asked about.

    The real solver cannot answer "which orders did the check pick?" -- the
    selection happens inside `check_reciprocity` and leaves no trace on the
    report. That made the whole strategy unverifiable, which is how mutation
    testing found `np.argsort(strength)` -> `np.argsort(None)` surviving.

    Reciprocal by construction: efficiency depends on the order index alone,
    never on the incidence azimuth, so a correct check must find zero
    violation and any failure is the check's own.
    """

    def __init__(self, strength_of, polarization_sensitive: bool = False) -> None:
        from gratinglab.solvers.base import Capabilities

        self.capabilities = Capabilities(name="recording", rigorous=False)
        self._strength_of = strength_of
        self._polarization_sensitive = polarization_sensitive
        self.calls: list[tuple[float, str]] = []

    def solve(self, problem, illumination, wavelengths, **options):
        """A physically consistent scan with dictated efficiencies.

        The propagating set comes from the real grating equation, not from a
        fixed range: `check_reciprocity` re-derives it independently and skips
        anything that disagrees, so a stand-in claiming orders that do not
        propagate would silently have every one of them dropped.
        """
        from gratinglab.geometry import is_propagating, order_range, sin_beta

        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=float))
        self.calls.append((illumination.alpha_deg, illumination.polarization))

        orders = order_range(
            float(wavelengths.min()), problem.period,
            illumination.sin_alpha, illumination.sin_gamma,
        )
        values = np.array([self._strength_of(int(m)) for m in orders], dtype=float)
        if self._polarization_sensitive and illumination.polarization != "TM":
            # A backend that resolves polarization answers differently. Scalar
            # does not, which is why this has to be dictated rather than found.
            values = values * 0.5

        efficiency = np.tile(values, (len(wavelengths), 1))
        propagating = np.array([
            is_propagating(sin_beta(orders, float(lam), problem.period,
                                    illumination.sin_alpha, illumination.sin_gamma))
            for lam in wavelengths
        ])
        return EfficiencyScan(
            wavelengths=wavelengths,
            orders=orders,
            efficiency=np.where(propagating, efficiency, 0.0),
            propagating=propagating,
            provenance=Provenance("recording"),
        )

    def reversed_azimuths(self) -> list[float]:
        """Every alpha the check solved at *except* the forward one."""
        return [alpha for alpha, _ in self.calls[1:]]


def orders_behind(azimuths, problem, illumination, wavelength):
    """Which order each reversed azimuth came from.

    `check_reciprocity` sets the reversed incidence to beta_m, so inverting
    that recovers the order it chose to test.
    """
    from gratinglab.geometry import beta as beta_from_sin
    from gratinglab.geometry import order_range, sin_beta

    # Searched over the orders this geometry actually has, not a fixed window:
    # a hardcoded range silently truncates the answer when the geometry moves,
    # which is how this helper first under-reported by five orders.
    candidates = order_range(
        wavelength, problem.period, illumination.sin_alpha, illumination.sin_gamma
    )
    found = []
    for alpha in azimuths:
        for order in candidates:
            sine = sin_beta(
                order, wavelength, problem.period,
                illumination.sin_alpha, illumination.sin_gamma,
            )
            if abs(sine) <= 1.0 and abs(np.degrees(beta_from_sin(sine)) - alpha) < 1e-9:
                found.append(order)
                break
    return found


class TestWhichOrdersGetTested:
    """The strategy `max_orders` implements, which was entirely unverified --
    18 of the 33 survivors in the full mutation sweep were in this function,
    and most of them here.

    A short wavelength on purpose: the cap only does anything when there are
    more orders than it allows, and the visible-light geometry used elsewhere
    in this file propagates four.
    """

    PROBLEM = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
    ILL = Illumination.classical(alpha=10.0, polarization=UNPOL)
    WAVELENGTH = 150.0

    @property
    def live(self):
        """The orders this geometry actually propagates, derived rather than
        counted -- a hardcoded number would be wrong the moment the geometry
        moved."""
        from gratinglab.geometry import order_range

        return order_range(
            self.WAVELENGTH, self.PROBLEM.period,
            self.ILL.sin_alpha, self.ILL.sin_gamma,
        )

    def _run(self, strength_of, **kwargs):
        solver = RecordingSolver(strength_of)
        report = check_reciprocity(
            solver, self.PROBLEM, self.ILL, [self.WAVELENGTH], **kwargs
        )
        chosen = orders_behind(
            solver.reversed_azimuths(), self.PROBLEM, self.ILL, self.WAVELENGTH
        )
        return report, chosen

    def test_the_geometry_has_more_orders_than_the_cap(self):
        """Non-vacuity for the whole class: with fewer, nothing is ever
        trimmed and every selection test below passes trivially."""
        assert len(self.live) > 5

    def test_the_strongest_orders_are_the_ones_tested(self):
        """The docstring says "keep the strongest orders: a violation there
        matters most". Nothing checked it, so the sort could have been
        arbitrary -- reciprocity holds for whichever orders get picked, and
        every existing test would still have passed."""
        strongest = sorted(self.live, key=abs)[:3]
        _, chosen = self._run(lambda m: 1.0 / (1 + abs(m)), max_orders=3)
        assert sorted(chosen) == sorted(strongest)

    def test_and_not_the_weakest(self):
        """Non-vacuity, and the mutant that survived: reversing the sort picks
        the far orders instead, and with this ranking the two sets are
        disjoint."""
        weakest = sorted(self.live, key=abs)[-3:]
        _, chosen = self._run(lambda m: 1.0 / (1 + abs(m)), max_orders=3)
        assert not (set(weakest) & set(chosen))

    def test_the_ranking_follows_the_efficiencies_not_the_order_index(self):
        """A second ranking with a different answer, so the test cannot be
        satisfied by any fixed choice of orders -- here the *outermost* three
        are the strong ones."""
        favoured = set(sorted(self.live, key=abs)[-3:])
        _, chosen = self._run(
            lambda m: 1.0 if m in favoured else 0.01, max_orders=3
        )
        assert sorted(chosen) == sorted(favoured)

    def test_every_order_is_tested_when_the_cap_allows(self):
        report, chosen = self._run(lambda m: 1.0, max_orders=None)
        assert len(chosen) == len(self.live)
        assert report.pairs_tested == len(self.live)

    def test_a_cap_equal_to_the_order_count_trims_nothing(self):
        """The boundary. `>` becoming `>=` here would drop an order whenever
        the cap exactly matched, which is the off-by-one mutation testing keeps
        finding -- `geometry.beta` had the same shape."""
        report, _ = self._run(lambda m: 1.0, max_orders=len(self.live))
        assert report.pairs_tested == len(self.live)

    def test_and_one_below_it_trims_exactly_one(self):
        report, _ = self._run(lambda m: 1.0, max_orders=len(self.live) - 1)
        assert report.pairs_tested == len(self.live) - 1


class TestTheGrazingSkip:
    """Orders diffracting within half a degree of grazing are skipped, because
    the reciprocal illumination is undefined there. Nothing tested that the
    skip happens, nor that it skips only the one order."""

    PROBLEM = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
    WAVELENGTH = 600.0
    #: Chosen, not stumbled on: it solves `sin(beta_1) = sin(89.7 deg)` for
    #: alpha, putting order +1 inside the margin while four others stay well
    #: clear -- which is what lets the `continue`-not-`break` test mean
    #: something.
    GRAZING_ALPHA = -34.8489

    def _near_grazing(self, ill):
        from gratinglab.geometry import beta as beta_from_sin
        from gratinglab.geometry import order_range, sin_beta

        return {
            int(m)
            for m in order_range(
                self.WAVELENGTH, self.PROBLEM.period,
                ill.sin_alpha, ill.sin_gamma,
            )
            if abs(
                np.degrees(
                    beta_from_sin(
                        sin_beta(m, self.WAVELENGTH, self.PROBLEM.period,
                                 ill.sin_alpha, ill.sin_gamma)
                    )
                )
            )
            > 89.5
        }

    def _run(self, alpha, **kwargs):
        ill = Illumination.classical(alpha=alpha, polarization=UNPOL)
        solver = RecordingSolver(lambda m: 1.0)
        report = check_reciprocity(
            solver, self.PROBLEM, ill, [self.WAVELENGTH], **kwargs
        )
        chosen = orders_behind(
            solver.reversed_azimuths(), self.PROBLEM, ill, self.WAVELENGTH
        )
        return report, chosen, ill

    def test_the_chosen_geometry_really_has_one(self):
        """Non-vacuity first: without a near-grazing order the skip has
        nothing to skip and both tests below pass for free."""
        ill = Illumination.classical(alpha=self.GRAZING_ALPHA, polarization=UNPOL)
        assert self._near_grazing(ill) == {1}

    def test_it_is_skipped(self):
        _, chosen, ill = self._run(self.GRAZING_ALPHA, max_orders=None)
        assert not (self._near_grazing(ill) & set(chosen))

    def test_but_the_orders_after_it_are_still_tested(self):
        """`continue`, not `break`. Stopping at the first grazing order would
        silently drop everything past it -- here orders -3 through 0 come
        before +1 and nothing comes after, so the count is what catches it."""
        report, chosen, ill = self._run(self.GRAZING_ALPHA, max_orders=None)
        from gratinglab.geometry import order_range

        live = order_range(
            self.WAVELENGTH, self.PROBLEM.period, ill.sin_alpha, ill.sin_gamma
        )
        assert report.pairs_tested == len(live) - 1
        assert len(chosen) == len(live) - 1

    def test_a_geometry_with_nothing_near_grazing_loses_no_order(self):
        """The margin must not be trimming orders in the ordinary case."""
        report, _, ill = self._run(10.0, max_orders=None)
        from gratinglab.geometry import order_range

        live = order_range(
            self.WAVELENGTH, self.PROBLEM.period, ill.sin_alpha, ill.sin_gamma
        )
        assert self._near_grazing(ill) == set()
        assert report.pairs_tested == len(live)


class TestTheReversedIlluminationMatches:
    def test_it_keeps_the_polarization(self):
        """Reciprocity compares one geometry with its reverse -- everything
        else has to be held fixed. Dropping `polarization` leaves the reverse
        solve on `Illumination`'s default of TE, and the check then compares
        two different physical problems and blames the solver.

        Invisible to the scalar solver, which neglects polarization entirely,
        so it needs a backend that resolves it. That is exactly the class of
        bug the first contributed RCWA would hit.
        """
        solver = RecordingSolver(lambda m: 1.0, polarization_sensitive=True)
        report = check_reciprocity(
            solver,
            Problem(period=1400.0, profile=Blazed(blaze_angle=30.0)),
            Illumination.classical(alpha=10.0, polarization="TM"),
            [600.0],
            max_orders=3,
        )
        assert report.passed
        assert {p for _, p in solver.calls} == {"TM"}

    def test_and_the_stand_in_really_would_notice(self):
        """Non-vacuity for the test above: the same solver answers differently
        for a different polarization, so a mismatch could not go unseen."""
        solver = RecordingSolver(lambda m: 1.0, polarization_sensitive=True)
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        tm = solver.solve(problem, Illumination.classical(alpha=10.0, polarization="TM"), [600.0])
        te = solver.solve(problem, Illumination.classical(alpha=10.0, polarization="TE"), [600.0])
        assert not np.array_equal(tm.efficiency, te.efficiency)


class TestTheEmptyReport:
    """`max_orders=0` reaches the early return. Only `pairs_tested` and
    `passed` were asserted, so the rest of the report could say anything."""

    def _empty(self):
        return check_reciprocity(
            scalar,
            Problem(period=1400.0, profile=Blazed(blaze_angle=30.0)),
            Illumination.classical(alpha=10.0, polarization=UNPOL),
            [600.0],
            max_orders=0,
            tolerance=1e-7,
        )

    def test_it_claims_no_violation_it_did_not_measure(self):
        """A report that never compared anything must not report a violation
        of 1.0 -- that is a measurement it did not make, which is precisely
        what `Provenance` exists to prevent elsewhere."""
        assert self._empty().max_violation == 0.0

    def test_and_carries_the_tolerance_it_was_given(self):
        """Otherwise the report cannot say what it would have accepted."""
        assert self._empty().tolerance == 1e-7

    def test_and_names_no_worst_case(self):
        report = self._empty()
        assert report.worst_order is None
        assert report.worst_wavelength is None

    def test_and_still_does_not_pass(self):
        assert not self._empty().passed
