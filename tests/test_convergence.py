r"""The convergence harness.

Two things are being tested and they are worth keeping apart. Against the real
scalar solver: that a sweep of a real method reaches a real plateau, and that
the measured non-monotonicity which motivated the plateau rule is still there.
Against a scripted stand-in: that the *rule* is what the docstring says it is,
which cannot be shown with a solver whose difference sequence we do not
control.
"""

import dataclasses

import numpy as np
import pytest

from gratinglab.convergence import (
    DEFAULT_PLATEAU,
    check_convergence,
    converged_scan,
    default_ladder,
    doubling_ladder,
)
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, Sinusoidal
from gratinglab.result import EfficiencyScan, Provenance
from gratinglab.solvers import scalar
from gratinglab.solvers.base import Capabilities, UnsupportedConfiguration

UNPOL = "unpolarized"
PERIOD, DELTA, ANTI, ALPHA, GAMMA = 315.15, 29.5, 70.5, 25.0, 1.5

#: Short on purpose. The sweep runs the solver up to ten times, and the
#: plateau is a property of the profile rather than of how finely the scan is
#: sampled -- see `test_the_verdict_does_not_depend_on_scan_length`.
WAVELENGTHS = np.linspace(1.0, 5.0, 9)


def blazed():
    return Problem(period=PERIOD, profile=Blazed(blaze_angle=DELTA, antiblaze_angle=ANTI))


def offplane():
    return Illumination.offplane(graze=GAMMA, azimuth=ALPHA, polarization=UNPOL)


def sweep(**kwargs):
    return check_convergence(scalar, blazed(), offplane(), WAVELENGTHS, **kwargs)


class TestTheLadder:
    def test_it_doubles(self):
        assert doubling_ladder(256, 4) == (256, 512, 1024, 2048)

    def test_one_rung_cannot_be_compared_to_anything(self):
        with pytest.raises(ValueError, match="at least two rungs"):
            doubling_ladder(256, 1)

    def test_a_non_positive_start_is_refused(self):
        with pytest.raises(ValueError, match="must be positive"):
            doubling_ladder(0, 4)

    def test_the_default_straddles_the_solvers_own_default(self):
        """Read from the signature, so a backend that changes its default does
        not leave the sweep quietly exploring the wrong range."""
        import inspect

        declared = inspect.signature(scalar.solve).parameters["quadrature_points"].default
        ladder = default_ladder(scalar)
        assert ladder[0] < declared < ladder[-1]


class TestARealSweepOfARealSolver:
    def test_the_blazed_case_converges(self):
        report = sweep()
        assert report.passed
        assert report.scan.provenance.converged is True

    def test_and_finds_a_setting_cheaper_than_the_one_it_had_to_reach(self):
        """The actionable output. `converged_at` is what a production run
        should use; the finest value swept is only what it took to prove it."""
        report = sweep()
        assert report.converged_at < report.values[-1]

    def test_the_evidence_travels_with_the_scan(self):
        """A result has to be able to defend itself away from the report
        object that produced it."""
        notes = sweep().scan.provenance.notes["convergence"]
        assert notes["knob"] == "quadrature_points"
        assert notes["converged_at"] in notes["values"]
        assert len(notes["differences"]) == len(notes["values"]) - 1
        assert all(d <= notes["tolerance"] for d in notes["differences"][-DEFAULT_PLATEAU:])

    def test_the_returned_scan_is_a_real_solve_at_the_finest_value(self):
        """Not a second source of truth: the harness reports on solves, it
        does not synthesise one."""
        report = sweep()
        direct = scalar.solve(
            blazed(), offplane(), WAVELENGTHS, quadrature_points=report.values[-1]
        )
        assert np.array_equal(report.scan.efficiency, direct.efficiency)

    def test_a_smooth_profile_converges_immediately(self):
        """Contrast, and the reason a single ladder cannot be hardcoded per
        method: a sinusoid is analytic, so the quadrature is spectrally
        accurate and the first rung is already at machine precision."""
        report = check_convergence(
            scalar,
            Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.15)),
            offplane(),
            WAVELENGTHS,
        )
        assert report.passed
        assert report.converged_at == report.values[0]
        assert max(report.differences) < 1e-12

    def test_and_the_blazed_one_does_not(self):
        """Non-vacuity for the contrast above: without this, both cases could
        be converging on the first rung and the ladder would be doing
        nothing."""
        assert sweep().converged_at > 256

    def test_the_verdict_does_not_depend_on_how_finely_the_scan_is_sampled(self):
        """Guards the short WAVELENGTHS grid this file uses for speed. The
        plateau is a property of the profile and the quadrature, not of the
        wavelength count."""
        dense = check_convergence(
            scalar, blazed(), offplane(), np.linspace(1.0, 5.0, 25)
        )
        assert dense.converged_at == sweep().converged_at


class TestConvergenceCanFail:
    """`converged=False` is a result, not an absence of one."""

    def test_an_unreachable_tolerance_is_reported_not_chased(self):
        report = sweep(tolerance=1e-18, values=doubling_ladder(256, 5))
        assert not report.passed
        assert report.converged_at is None
        assert report.scan.provenance.converged is False

    def test_and_says_so_in_a_warning_a_reader_will_see(self):
        report = sweep(tolerance=1e-18, values=doubling_ladder(256, 5))
        text = " ".join(report.scan.provenance.warnings)
        assert "convergence not demonstrated" in text
        assert "not defensible" in text

    def test_a_passing_sweep_adds_no_such_warning(self):
        """Non-vacuity: the warning must be about this case, not about every
        case."""
        text = " ".join(sweep().scan.provenance.warnings)
        assert "convergence not demonstrated" not in text

    def test_the_whole_ladder_runs_when_nothing_plateaus(self):
        """A sweep that gives up early would understate how hard it tried."""
        ladder = doubling_ladder(256, 5)
        assert sweep(tolerance=1e-18, values=ladder).values == ladder


class TestProgressThroughASweep:
    """Two granularities on purpose: called per wavelength, reporting per
    rung. See `check_convergence`'s docstring for why."""

    def _run(self, **kwargs):
        seen = []
        report = sweep(progress=lambda done, total: seen.append((done, total)), **kwargs)
        return seen, report

    def test_the_numbers_are_rungs_not_wavelengths(self):
        seen, report = self._run()
        assert {t for _, t in seen} == {len(default_ladder(scalar))}
        assert max(d for d, _ in seen) == len(report.values) - 1

    def test_they_never_go_backwards(self):
        """The property that rules out a composed per-wavelength counter: it
        would restart on every rung and drive a bar backwards."""
        seen, _ = self._run()
        dones = [d for d, _ in seen]
        assert dones == sorted(dones)

    def test_but_it_is_called_once_per_wavelength(self):
        """The frequency is the point -- every wavelength is a chance for the
        caller to raise, so cancelling does not have to wait out a rung, which
        is up to a quarter of the whole sweep.

        Exact rather than a threshold: rungs x the solver's own (n + 1).
        """
        seen, report = self._run()
        assert len(seen) == len(report.values) * (len(WAVELENGTHS) + 1)
        assert len(seen) > len(report.values), "and far more often than per rung"

    def test_it_does_not_end_at_the_total(self):
        """Finishing early is the *good* outcome -- the sweep stopped at the
        first plateau. Running the counter to the end so a bar looks tidy
        would claim rungs that were never solved."""
        seen, report = self._run()
        assert seen[-1][0] < seen[-1][1]
        assert len(report.values) < len(default_ladder(scalar))

    def test_a_sweep_that_uses_the_whole_ladder_reports_every_rung(self):
        """Non-vacuity for the test above: the early stop is a property of
        this case, not of the reporting."""
        ladder = doubling_ladder(256, 4)
        seen, report = self._run(tolerance=1e-18, values=ladder)
        assert report.values == ladder
        assert max(d for d, _ in seen) == len(ladder) - 1

    def test_raising_aborts_the_sweep_mid_rung(self):
        from gratinglab.solvers.base import SolveCancelled

        calls = {"n": 0}

        def stop(done, total):
            calls["n"] += 1
            if calls["n"] > 3:
                raise SolveCancelled("enough")

        with pytest.raises(SolveCancelled):
            sweep(progress=stop)
        assert calls["n"] == 4, "it stopped where it was told, not at a rung boundary"

    def test_and_no_partial_report_comes_back(self):
        """A sweep that was stopped has not demonstrated anything. Returning a
        report implying otherwise would be the mistake this module exists to
        prevent."""
        from gratinglab.solvers.base import SolveCancelled

        def stop(done, total):
            raise SolveCancelled("immediately")

        with pytest.raises(SolveCancelled):
            sweep(progress=stop)

    def test_not_passing_it_changes_nothing(self):
        assert self._run()[1].converged_at == sweep().converged_at


class TestThePlateauRule:
    """The load-bearing behaviour, tested against a scripted solver.

    A real solver's difference sequence cannot be dictated, and this rule
    exists precisely to handle a sequence that dips below tolerance and comes
    back out. So the sequence is supplied here.
    """

    def test_a_single_accidental_agreement_does_not_count(self):
        """The measured blazed case in miniature: one small step, then a
        larger one. Stopping at the dip would certify a tolerance the very
        next refinement violates."""
        solver = ScriptedSolver([1e-2, 1e-9, 1e-2, 1e-9, 1e-9])
        report = check_convergence(
            solver, blazed(), offplane(), WAVELENGTHS,
            values=doubling_ladder(4, 6), tolerance=1e-6,
        )
        assert report.differences[1] < 1e-6, "the dip must really be a dip"
        # E(32), E(64) and E(128) are the three that actually agree, and 32 is
        # the coarsest of them. The dip at E(8) is passed over.
        assert report.converged_at == 32

    def test_and_plateau_1_would_have_stopped_at_the_dip(self):
        """Non-vacuity, and the whole argument for DEFAULT_PLATEAU > 1: the
        same sequence, judged the tempting way, stops three rungs early on
        evidence the next solve contradicts."""
        solver = ScriptedSolver([1e-2, 1e-9, 1e-2, 1e-9, 1e-9])
        report = check_convergence(
            solver, blazed(), offplane(), WAVELENGTHS,
            values=doubling_ladder(4, 6), tolerance=1e-6, plateau=1,
        )
        assert report.converged_at == 8

    def test_a_stricter_plateau_demands_more_agreements(self):
        solver = ScriptedSolver([1e-9] * 5)
        report = check_convergence(
            solver, blazed(), offplane(), WAVELENGTHS,
            values=doubling_ladder(4, 6), tolerance=1e-6, plateau=3,
        )
        assert len(report.differences) == 3

    def test_requiring_no_agreement_at_all_is_refused(self):
        with pytest.raises(ValueError, match="demonstrates nothing"):
            sweep(plateau=0)


class TestTheMeasurementBehindTheRule:
    """The docstring's table, pinned.

    If the solver ever changes so that quadrature error becomes monotone in
    `n`, the argument for a plateau weakens and the docstring is stale. This
    fails first.
    """

    def test_the_error_still_grows_on_at_least_one_refinement(self):
        report = sweep(tolerance=0.0, values=doubling_ladder(256, 6))
        grew = [
            (a, b) for a, b in zip(report.differences, report.differences[1:]) if b > a
        ]
        assert grew, f"expected non-monotone refinement, got {report.differences}"

    def test_the_kink_that_causes_it_is_between_nodes(self):
        """Why it happens, not merely that it does: the blazed apex sits at an
        awkward fraction of a period, so doubling `n` changes how near a node
        lands to it."""
        apex = 1.0 / (1.0 + np.tan(np.radians(DELTA)) / np.tan(np.radians(ANTI)))
        assert apex == pytest.approx(0.8331, abs=1e-4)
        for n in (256, 512, 1024, 2048):
            assert (apex * n) % 1.0 > 1e-6, f"n={n} happens to land on the kink"


class TestRefusals:
    def test_a_solver_with_no_knob_has_nothing_to_sweep(self):
        with pytest.raises(ValueError, match="no accuracy_knob"):
            check_convergence(
                ScriptedSolver([1e-9] * 3, knob=None), blazed(), offplane(), WAVELENGTHS
            )

    def test_the_swept_knob_cannot_also_be_pinned(self):
        with pytest.raises(ValueError, match="cannot also be pinned"):
            sweep(quadrature_points=2048)

    def test_values_must_ascend(self):
        with pytest.raises(ValueError, match="ascending and distinct"):
            sweep(values=(2048, 512, 1024))

    def test_a_single_value_compares_to_nothing(self):
        with pytest.raises(ValueError, match="at least two"):
            sweep(values=(2048,))

    def test_a_ladder_that_changes_the_order_grid_is_refused(self):
        """Should never happen -- the orders come from the grating equation,
        not from the accuracy knob. It raises rather than silently comparing
        misaligned arrays, which would produce a number that looks like a
        convergence measure and is not."""
        with pytest.raises(ValueError, match="different order grids"):
            check_convergence(
                ScriptedSolver([1e-9] * 3, widen_orders_after=1),
                blazed(), offplane(), WAVELENGTHS, values=doubling_ladder(4, 4),
            )


class TestRungsTheSolverRefuses:
    """A coarse rung below a method's own floor is a fact about the ladder,
    not about convergence."""

    #: 8 is below scalar's hard floor of 16; 16 clears that and then fails the
    #: Nyquist guard, since this geometry reaches order 11 and needs more than
    #: 22 points. Two rungs, two different refusals.
    LOW = (8, 16, 4096, 8192, 16384)

    def test_a_leading_refusal_is_recorded_and_stepped_over(self):
        """Exactly what a ladder starting low will hit, and not a statement
        about convergence."""
        report = sweep(values=self.LOW)
        assert [v for v, _ in report.skipped] == [8, 16]
        assert report.values[0] == 4096

    def test_the_reason_is_kept_not_just_the_fact(self):
        reasons = dict(sweep(values=self.LOW).skipped)
        assert "at least 16" in reasons[8]
        assert "Nyquist" in reasons[16]

    def test_a_refusal_after_a_success_is_a_real_failure_and_propagates(self):
        """Not a floor -- the ladder walked into something. Logging it as a
        skipped rung would leave a sweep that looks complete and is not."""
        with pytest.raises(UnsupportedConfiguration):
            check_convergence(
                ScriptedSolver([1e-9] * 4, refuse_after=1),
                blazed(), offplane(), WAVELENGTHS, values=doubling_ladder(4, 5),
            )

    def test_a_ladder_no_rung_of_which_runs_says_so(self):
        with pytest.raises(RuntimeError, match="no rung"):
            sweep(values=(8, 16))

    def test_and_one_solved_rung_is_still_nothing_to_compare(self):
        """A distinct message, because it is a distinct situation: the sweep
        did run, it just never got a second point."""
        with pytest.raises(RuntimeError, match="only one rung"):
            sweep(values=(16, 4096))


class TestTheConvenienceWrapper:
    def test_it_returns_the_stamped_scan(self):
        scan = converged_scan(scalar, blazed(), offplane(), WAVELENGTHS)
        assert isinstance(scan, EfficiencyScan)
        assert scan.provenance.converged is True
        assert scan.provenance.is_defensible

    def test_which_is_what_is_defensible_meant_all_along(self):
        """`Provenance.is_defensible` has been False for every result this
        project has ever produced. This is the first time it can be True."""
        plain = scalar.solve(blazed(), offplane(), WAVELENGTHS)
        assert plain.provenance.converged is None
        assert not plain.provenance.is_defensible


class ScriptedSolver:
    """A solver whose successive answers differ by a dictated amount.

    Not a mock of the scalar solver -- it computes no physics and claims none.
    It exists so the *rule* can be tested against a difference sequence chosen
    to exercise it, which no real solver can be asked to produce on demand.
    """

    def __init__(
        self,
        steps,
        *,
        knob: str | None = "grid",
        widen_orders_after: int | None = None,
        refuse_after: int | None = None,
    ) -> None:
        self.capabilities = Capabilities(name="scripted", accuracy_knob=knob)
        self._steps = list(steps)
        self._calls = 0
        self._widen_after = widen_orders_after
        self._refuse_after = refuse_after

    def solve(self, problem, illumination, wavelengths, **options):
        index = self._calls
        self._calls += 1
        if self._refuse_after is not None and index > self._refuse_after:
            raise UnsupportedConfiguration("scripted refusal")

        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=float))
        orders = np.array([0, 1] if index <= (self._widen_after or index) else [0, 1, 2])
        # Each answer sits its dictated distance from the one before, so the
        # measured difference sequence is exactly `steps`.
        value = 0.5 + sum(self._steps[:index])
        efficiency = np.full((len(wavelengths), len(orders)), value)
        return EfficiencyScan(
            wavelengths=wavelengths,
            orders=orders,
            efficiency=efficiency,
            propagating=np.ones_like(efficiency, dtype=bool),
            provenance=Provenance(method="scripted"),
        )


def test_the_scripted_solver_really_does_step_by_the_dictated_amounts():
    """Non-vacuity for every test in TestThePlateauRule: if this stand-in did
    not produce the sequence it is told to, those tests would be asserting
    against a fiction."""
    steps = [1e-2, 1e-9, 1e-2]
    solver = ScriptedSolver(steps)
    report = check_convergence(
        solver, blazed(), offplane(), WAVELENGTHS,
        values=doubling_ladder(4, 4), tolerance=0.0,
    )
    assert report.differences == pytest.approx(steps, rel=1e-9)


def test_the_report_is_immutable():
    """It is evidence. Evidence that can be edited after the fact is not."""
    report = sweep()
    with pytest.raises(dataclasses.FrozenInstanceError):
        report.converged_at = 8
