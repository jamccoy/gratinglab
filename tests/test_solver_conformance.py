"""Every registered solver, held to the same contract.

`test_solver_progress.py` pins the contract to the scalar backend by name; this
file states it once and runs it against whatever `available_solvers()` returns,
so a contributed backend -- registered directly or arriving through the
``gratinglab.solvers`` entry point -- is examined the moment it registers,
without anyone writing a test for it. Today that is one solver; the file exists
because the roadmap's RCWA integration makes it two.

Each backend gets a configuration small enough to finish in milliseconds,
chosen per method: the accuracy knob is turned *down*, because conformance is
about the contract, not the physics -- convergence has its own harness.
"""

import inspect
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, FromProfileData, Sinusoidal
from gratinglab.result import EfficiencyScan
from gratinglab.solvers import (
    Capabilities,
    SolveCancelled,
    UnsupportedConfiguration,
    available_solvers,
    get_solver,
)

#: The full polarization vocabulary, for probing what a backend refuses.
ALL_POLARIZATIONS = ("TE", "TM", "unpolarized")

#: A boundary that doubles back in ``t`` -- representable only by a method that
#: parametrises the curve. Everything else must refuse it.
UNDERCUT = FromProfileData(t=(0.0, 0.55, 0.45, 1.0), y=(0.0, 0.4, 0.1, 0.0))


@dataclass(frozen=True)
class Case:
    problem: Problem
    illumination: Illumination
    wavelengths: np.ndarray
    options: dict[str, Any] = field(default_factory=dict)

    def solve(self, solver, **extra):
        return solver.solve(
            self.problem, self.illumination, self.wavelengths, **self.options, **extra
        )


def small_case(name: str, capabilities: Capabilities) -> Case:
    """The smallest legal configuration for a given backend.

    Three wavelengths spanning enough range that the propagating-order count
    changes across the scan -- which is what makes the evanescent-zero test
    below mean something. The final adaptation step exists for a backend this
    table has no row for: it gets the default case bent into the envelope its
    capabilities declare, rather than a spurious refusal.
    """
    if name == "integral":
        case = Case(
            problem=Problem(period=500.0, profile=Sinusoidal(depth_fraction=0.05)),
            illumination=Illumination.classical(alpha=10.0, polarization="TE"),
            wavelengths=np.array([350.0, 500.0, 650.0]),
            options={"boundary_points": 64},
        )
    else:
        # The scalar row, and the default: the off-plane blazed case the rest
        # of the suite uses, cut to three wavelengths.
        case = Case(
            problem=Problem(
                period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
            ),
            illumination=Illumination.offplane(
                graze=1.5, azimuth=25.0, polarization="unpolarized"
            ),
            wavelengths=np.array([1.0, 2.5, 5.0]),
        )

    illumination = case.illumination
    if not capabilities.conical and not illumination.is_in_plane:
        illumination = Illumination.classical(
            alpha=illumination.alpha_deg, polarization=illumination.polarization
        )
    if illumination.polarization not in capabilities.polarizations:
        illumination = illumination.model_copy(
            update={"polarization": capabilities.polarizations[0]}
        )
    return Case(case.problem, illumination, case.wavelengths, case.options)


@pytest.fixture(params=available_solvers())
def name(request):
    return request.param


@pytest.fixture
def solver(name):
    return get_solver(name)


@pytest.fixture
def case(name, solver):
    return small_case(name, solver.capabilities)


class Recorder:
    """Collects every progress call, optionally raising at a chosen one.

    Same shape as the one in `test_solver_progress.py`; the raised instance is
    kept so identity can be asserted -- the contract is that the exception
    propagates *unchanged*.
    """

    def __init__(self, raise_at: int | None = None) -> None:
        self.calls: list[tuple[int, int]] = []
        self.raised: SolveCancelled | None = None
        self._raise_at = raise_at

    def __call__(self, done: int, total: int) -> None:
        self.calls.append((done, total))
        if self._raise_at is not None and done >= self._raise_at:
            self.raised = SolveCancelled(f"stopped at {done}")
            raise self.raised


class TestTheDeclaration:
    def test_the_registry_key_is_the_declared_name(self, name, solver):
        """`register` keys on `capabilities.name`, so these cannot disagree
        today -- this pins that they never start to, e.g. by a registry that
        begins accepting aliases."""
        assert solver.capabilities.name == name

    def test_the_declared_accuracy_knob_is_a_real_keyword(self, solver):
        """The convergence harness sweeps this parameter by name; a knob that
        `solve()` does not actually take would fail there, far from the
        declaration that caused it."""
        knob = solver.capabilities.accuracy_knob
        if knob is None:
            pytest.skip(f"{solver.capabilities.name} declares no accuracy knob")
        assert knob in inspect.signature(solver.solve).parameters


class TestRefusingWhatItCannotDo:
    """`Capabilities.check` exists so a backend refuses loudly instead of
    approximating silently -- the failure mode that makes cross-method plots
    lie. Each test asserts both directions: refused when out of scope,
    accepted when declared, so neither branch can go vacuous."""

    def test_an_undercut_profile_is_refused_unless_handled(self, solver, case):
        caps = solver.capabilities
        undercut = Problem(period=case.problem.period, profile=UNDERCUT)
        assert not UNDERCUT.is_single_valued(), "the probe must actually double back"
        if caps.handles_undercut:
            caps.check(undercut, case.illumination)
        else:
            with pytest.raises(UnsupportedConfiguration):
                caps.check(undercut, case.illumination)

    def test_an_off_plane_mount_is_refused_unless_conical(self, solver, case):
        caps = solver.capabilities
        offplane = Illumination.offplane(
            graze=1.5, azimuth=25.0, polarization=case.illumination.polarization
        )
        if caps.conical:
            caps.check(case.problem, offplane)
        else:
            with pytest.raises(UnsupportedConfiguration):
                caps.check(case.problem, offplane)

    def test_only_declared_polarizations_pass(self, solver, case):
        caps = solver.capabilities
        for polarization in ALL_POLARIZATIONS:
            probe = case.illumination.model_copy(update={"polarization": polarization})
            if polarization in caps.polarizations:
                caps.check(case.problem, probe)
            else:
                with pytest.raises(UnsupportedConfiguration):
                    caps.check(case.problem, probe)

    def test_the_refusal_branch_is_reachable(self):
        """Non-vacuity for the loop above: a backend declaring all three
        polarizations never exercises its `raises` arm, so demonstrate on a
        restricted declaration that the machinery the loop relies on works."""
        restricted = Capabilities(name="restricted", polarizations=("TE",))
        problem = Problem(period=100.0, profile=Sinusoidal(depth_fraction=0.05))
        with pytest.raises(UnsupportedConfiguration):
            restricted.check(problem, Illumination.classical(alpha=5.0, polarization="TM"))


class TestTheScanItReturns:
    def test_it_is_a_valid_scan_attributed_to_the_registry_name(
        self, name, solver, case
    ):
        """`EfficiencyScan.__post_init__` already enforces the shape and
        NaN-free invariants, so a successful construction is itself half the
        test; what it cannot know is *which* solver ran, and the comparison
        harness groups by exactly that string."""
        scan = case.solve(solver)
        assert isinstance(scan, EfficiencyScan)
        assert scan.provenance.method == name
        assert np.array_equal(scan.wavelengths, case.wavelengths)
        assert scan.efficiency.shape == (len(case.wavelengths), len(scan.orders))
        assert np.isfinite(scan.efficiency).all()
        assert (scan.efficiency >= 0.0).all()

    def test_evanescent_orders_are_exactly_zero(self, solver, case):
        """0.0 with ``propagating=False`` -- never NaN, never a stray residual
        (docs/conventions.md 4). Exact equality is the point: any non-zero
        value on a non-propagating order silently inflates the energy
        balance."""
        scan = case.solve(solver)
        assert (~scan.propagating).any(), (
            "the case must span a change in propagating-order count, or this "
            "test passes on an empty selection"
        )
        assert (scan.efficiency[~scan.propagating] == 0.0).all()


class TestTheProgressContract:
    """Only for a backend declaring ``reports_progress``; the declaration is
    precisely the promise that these hold (`Solver` protocol, points 1-4)."""

    @pytest.fixture
    def declared(self, name, solver):
        if not solver.capabilities.reports_progress:
            pytest.skip(f"{name} does not declare reports_progress")

    def test_the_first_call_is_zero_of_the_total(self, declared, solver, case):
        recorder = Recorder()
        case.solve(solver, progress=recorder)
        total = recorder.calls[0][1]
        assert recorder.calls[0] == (0, total)

    def test_done_is_monotone_and_ends_exactly_at_the_total(
        self, declared, solver, case
    ):
        recorder = Recorder()
        case.solve(solver, progress=recorder)
        dones = [done for done, _ in recorder.calls]
        totals = {total for _, total in recorder.calls}
        assert dones == sorted(dones)
        assert len(totals) == 1, "the total must not change mid-scan"
        assert recorder.calls[-1] == (totals.pop(),) * 2

    def test_cancelling_before_any_work_stops_at_the_first_call(
        self, declared, solver, case
    ):
        """The leading ``(0, total)`` is a cancellation point, not decoration:
        for RCWA a single wavelength can be a minute."""
        recorder = Recorder(raise_at=0)
        with pytest.raises(SolveCancelled):
            case.solve(solver, progress=recorder)
        assert len(recorder.calls) == 1

    def test_cancelling_mid_scan_propagates_unchanged_with_no_partial_scan(
        self, declared, solver, case
    ):
        """`pytest.raises` is the assertion that nothing was returned: a
        partially filled scan -- remaining rows still zero -- would be
        indistinguishable from one whose orders all passed off. Identity on
        the exception pins that no ``except`` re-wrapped it on the way out."""
        total = Recorder()
        case.solve(solver, progress=total)
        assert total.calls, "reports_progress is declared, yet no call arrived"
        if total.calls[-1][0] < 2:
            pytest.skip("the scan has no interior in this backend's units")
        recorder = Recorder(raise_at=1)
        with pytest.raises(SolveCancelled) as caught:
            case.solve(solver, progress=recorder)
        assert caught.value is recorder.raised
        assert recorder.calls[-1][0] == 1, "it stopped where the callback raised"
