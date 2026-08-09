r"""The progress hook, and cancellation that actually cancels.

Before this, cancelling a solve stopped the *waiting* and not the work: a
Python thread cannot be killed and a NumPy-bound call has no check point. The
callback is the check point. The load-bearing test here is
`test_raising_stops_the_solver_where_it_was_raised`, which measures how far the
solver got rather than taking its word for it -- everything else in this file
is the contract around that one fact.
"""

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.result import EfficiencyScan
from gratinglab.solvers import scalar
from gratinglab.solvers.base import Capabilities, SolveCancelled

UNPOL = "unpolarized"
PERIOD, DELTA, ANTI, ALPHA, GAMMA = 315.15, 29.5, 70.5, 25.0, 1.5

#: Long enough that "it stopped early" is a distinction worth drawing.
WAVELENGTHS = np.linspace(1.0, 5.0, 40)


def blazed():
    return Problem(period=PERIOD, profile=Blazed(blaze_angle=DELTA, antiblaze_angle=ANTI))


def offplane():
    return Illumination.offplane(graze=GAMMA, azimuth=ALPHA, polarization=UNPOL)


def solve(**kwargs):
    return scalar.solve(blazed(), offplane(), WAVELENGTHS, **kwargs)


class Recorder:
    """Collects every call, optionally raising at a chosen one."""

    def __init__(self, raise_at: int | None = None) -> None:
        self.calls: list[tuple[int, int]] = []
        self._raise_at = raise_at

    def __call__(self, done: int, total: int) -> None:
        self.calls.append((done, total))
        if self._raise_at is not None and done >= self._raise_at:
            raise SolveCancelled(f"stopped at {done}")

    @property
    def reached(self) -> int:
        return self.calls[-1][0] if self.calls else -1


class TestCancellationIsMeasured:
    """The point of the milestone, and the only tests here that could not
    have been written before it."""

    def test_raising_stops_the_solver_where_it_was_raised(self):
        """Not "cancellation was requested" -- how far the loop actually got.

        A solver that ignored the callback, or caught the exception, would
        run to 40 and this would fail by 30.
        """
        recorder = Recorder(raise_at=10)
        with pytest.raises(SolveCancelled):
            solve(progress=recorder)
        assert recorder.reached == 10
        assert len(WAVELENGTHS) == 40, "the case must have somewhere left to go"

    def test_and_without_raising_the_same_scan_runs_to_the_end(self):
        """Non-vacuity: without this, the test above would pass on a solver
        that stopped at 10 for some entirely different reason."""
        recorder = Recorder()
        solve(progress=recorder)
        assert recorder.reached == len(WAVELENGTHS)

    def test_a_cancelled_solve_returns_nothing_at_all(self):
        """A partially filled efficiency array -- remaining rows still zero --
        is indistinguishable from a scan whose orders all passed off. Handing
        one back would be exactly the silent wrong number this project is
        arranged against."""
        with pytest.raises(SolveCancelled):
            solve(progress=Recorder(raise_at=5))

    def test_the_solver_does_not_swallow_it(self):
        """Stated in the protocol and worth pinning: an `except Exception`
        anywhere around that loop would quietly restore the old behaviour, in
        which cancelling only stops the waiting."""
        sentinel = SolveCancelled("mine")

        def raise_mine(done, total):
            if done > 0:
                raise sentinel

        with pytest.raises(SolveCancelled) as caught:
            solve(progress=raise_mine)
        assert caught.value is sentinel

    def test_an_ordinary_exception_propagates_too(self):
        """The contract is about *any* exception from the callback, not a
        privileged type. A caller may have its own."""
        with pytest.raises(ZeroDivisionError):
            solve(progress=lambda d, t: 1 / 0)


class TestTheContract:
    def test_it_is_called_once_per_wavelength_plus_one(self):
        recorder = Recorder()
        solve(progress=recorder)
        assert len(recorder.calls) == len(WAVELENGTHS) + 1

    def test_it_starts_at_zero_before_any_work(self):
        """For RCWA a single wavelength can be a minute, so a cancellation
        point ahead of the first one is the difference between stopping now
        and stopping in a minute."""
        recorder = Recorder()
        solve(progress=recorder)
        assert recorder.calls[0] == (0, len(WAVELENGTHS))

    def test_and_ends_exactly_at_the_total(self):
        recorder = Recorder()
        solve(progress=recorder)
        assert recorder.calls[-1] == (len(WAVELENGTHS), len(WAVELENGTHS))

    def test_done_never_goes_backwards(self):
        recorder = Recorder()
        solve(progress=recorder)
        dones = [d for d, _ in recorder.calls]
        assert dones == sorted(dones)

    def test_the_total_never_changes(self):
        recorder = Recorder()
        solve(progress=recorder)
        assert len({t for _, t in recorder.calls}) == 1

    def test_the_count_does_not_depend_on_how_many_orders_propagate(self):
        """Reported at the top of each row rather than the bottom, so the
        loop's `continue` cannot skip a step -- and skipping one would skip a
        cancellation point, which matters more than the number.

        Long wavelengths here: p sin(gamma) = 8.25 nm, so past that only order
        0 survives, against nine orders at 1 nm.
        """
        sparse = np.array([20.0, 30.0, 40.0])
        recorder = Recorder()
        result = scalar.solve(blazed(), offplane(), sparse, progress=recorder)
        assert result.propagating.sum() == len(sparse), "one order each, not nine"
        assert len(recorder.calls) == len(sparse) + 1

    def test_even_though_the_empty_case_cannot_arise(self):
        """Why the top-of-row placement is defensive rather than load-bearing:
        `sin(beta_0) = -sin(alpha)` is real for every legal alpha, so order 0
        always propagates and the `continue` branch is unreachable for any
        valid illumination. Worth knowing before someone deletes the guard."""
        result = scalar.solve(blazed(), offplane(), np.array([1e4]))
        assert result.propagating.any()
        assert result.orders.tolist() == [0]


class TestNotPassingItChangesNothing:
    def test_the_default_is_bit_identical_to_passing_none(self):
        assert np.array_equal(solve().efficiency, solve(progress=None).efficiency)

    def test_and_to_reporting_progress(self):
        """The callback observes; it must not participate. If the reported
        numbers ever came from the same place the efficiency does, this
        drifts."""
        assert np.array_equal(solve().efficiency, solve(progress=Recorder()).efficiency)


class TestTheDeclaration:
    def test_scalar_says_it_reports_progress(self):
        assert scalar.capabilities.reports_progress

    def test_and_a_backend_that_has_not_been_updated_says_otherwise(self):
        """The default is False, so a third-party solver arriving through the
        entry point is never handed a keyword its signature cannot take."""
        assert not Capabilities(name="hypothetical").reports_progress

    def test_a_strict_signature_would_notice_being_handed_one(self):
        """Non-vacuity for the declaration: this is the failure the flag
        exists to prevent, demonstrated rather than asserted."""

        class Strict:
            capabilities = Capabilities(name="strict")

            def solve(self, problem, illumination, wavelengths):
                return None

        with pytest.raises(TypeError):
            Strict().solve(blazed(), offplane(), WAVELENGTHS, progress=Recorder())


class TestCancellationIsAvailableToAnyCaller:
    def test_a_deadline_is_expressible_in_the_same_way(self):
        """Nothing about the mechanism is specific to a GUI Cancel button --
        the callback is just a place to make a decision, and the caller owns
        the policy."""
        budget = {"left": 7}

        def spend(done, total):
            budget["left"] -= 1
            if budget["left"] <= 0:
                raise SolveCancelled("out of budget")

        with pytest.raises(SolveCancelled, match="out of budget"):
            solve(progress=spend)

    def test_a_scan_short_enough_to_finish_inside_the_budget_does(self):
        """Non-vacuity, and it keeps the test above from passing on a solver
        that raised for its own reasons."""
        budget = {"left": 7}

        def spend(done, total):
            budget["left"] -= 1
            if budget["left"] <= 0:
                raise SolveCancelled("out of budget")

        result = scalar.solve(blazed(), offplane(), WAVELENGTHS[:5], progress=spend)
        assert isinstance(result, EfficiencyScan)
