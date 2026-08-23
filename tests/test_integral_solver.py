"""The public integral solver: registration, physics checks, and the contract.

The physics lives in ``tests/test_integral_core.py``; this file covers what
the *solver* adds on top -- the conical reduction, the polarization average,
the guards, provenance, progress, and the capability nothing else has
(undercut boundaries, whose only cross-check is the physics itself).
"""

import numpy as np
import pytest

from gratinglab.checks import check_energy_balance, check_reciprocity
from gratinglab.convergence import check_convergence
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, FromProfileData, Sinusoidal
from gratinglab.solvers import get_solver, integral
from gratinglab.solvers.base import SolveCancelled, UnsupportedConfiguration

OFFPLANE = Problem(
    period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
)
SINUSOID = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3))
UNDERCUT = Problem(
    period=600.0,
    profile=FromProfileData(
        t=(0.0, 0.30, 0.28, 0.55, 1.0), y=(0.0, 0.45, 0.2, 0.5, 0.0)
    ),
)


class TestRegistration:
    def test_available_by_name(self):
        assert get_solver("integral") is integral

    def test_capabilities(self):
        caps = integral.capabilities
        assert caps.name == "integral"
        assert caps.conical
        assert caps.rigorous
        assert caps.handles_undercut
        assert caps.accuracy_knob == "boundary_points"
        assert caps.reports_progress
        assert set(caps.polarizations) == {"TE", "TM", "unpolarized"}


class TestConicalMount:
    """The reason this solver exists: extreme off-plane at gamma = 1.5 deg."""

    def test_energy_balance_is_a_theorem(self):
        scan = integral.solve(
            OFFPLANE,
            Illumination.offplane(graze=1.5, azimuth=25.0, polarization="TE"),
            [2.0, 3.5],
            boundary_points=96,
        )
        assert check_energy_balance(scan, tolerance=1e-5, lossless=True).passed

    def test_unpolarized_is_the_mean_of_te_and_tm(self):
        kwargs = dict(boundary_points=96)
        wavelengths = [2.5]
        scans = {
            pol: integral.solve(
                OFFPLANE,
                Illumination.offplane(graze=1.5, azimuth=25.0, polarization=pol),
                wavelengths,
                **kwargs,
            )
            for pol in ("TE", "TM", "unpolarized")
        }
        mean = 0.5 * (scans["TE"].efficiency + scans["TM"].efficiency)
        assert np.allclose(scans["unpolarized"].efficiency, mean, atol=1e-12)


class TestReciprocity:
    """Backend-agnostic Lorentz reciprocity, from checks.py unchanged.

    Measured: TE holds to ~1e-8; TM to ~1e-5 (the second-kind equation's
    discretisation error, dominated by the diagonal estimate). Both bounds
    are asserted where measured, not at the harness default of 1e-9."""

    @pytest.mark.parametrize("polarization, bound", [("TE", 1e-6), ("TM", 1e-4)])
    def test_in_plane(self, polarization, bound):
        report = check_reciprocity(
            integral,
            SINUSOID,
            Illumination.classical(alpha=25.0, polarization=polarization),
            [500.0],
            tolerance=bound,
            boundary_points=110,
        )
        assert report.passed


class TestUndercut:
    """The capability no other registered method has. Nothing can
    cross-check it, so the physics invariants *are* the validation."""

    def test_te_conserves_energy_tightly(self):
        scan = integral.solve(
            UNDERCUT,
            Illumination.classical(alpha=20.0, polarization="TE"),
            [450.0],
            boundary_points=128,
        )
        assert check_energy_balance(scan, tolerance=1e-5, lossless=True).passed

    def test_tm_solves_within_corner_tolerance(self):
        scan = integral.solve(
            UNDERCUT,
            Illumination.classical(alpha=20.0, polarization="TM"),
            [450.0],
            boundary_points=128,
        )
        assert check_energy_balance(scan, tolerance=0.05, lossless=True).passed


class TestRefusals:
    def test_finite_conductivity_is_not_quietly_approximated(self):
        with pytest.raises(UnsupportedConfiguration, match="Goray"):
            integral.solve(
                SINUSOID,
                Illumination.classical(alpha=25.0),
                [500.0],
                conductivity="tabulated",
            )

    def test_too_few_points_per_wavelength(self):
        """An X-ray scan on a coarse mesh is refused with the number needed,
        so the convergence ladder skips the rung instead of trusting it."""
        with pytest.raises(ValueError, match="nodes per transverse wavelength"):
            integral.solve(
                OFFPLANE,
                Illumination.offplane(graze=1.5, azimuth=25.0, polarization="TE"),
                [0.5],
                boundary_points=64,
            )

    def test_spectral_terms_must_span_the_orders(self):
        with pytest.raises(ValueError, match="spectral_terms"):
            integral.solve(
                OFFPLANE,
                Illumination.offplane(graze=1.5, azimuth=25.0, polarization="TE"),
                [2.0],
                boundary_points=96,
                spectral_terms=3,
            )


class TestProvenance:
    def test_records_the_boundary_condition_and_knobs(self):
        scan = integral.solve(
            SINUSOID,
            Illumination.classical(alpha=25.0, polarization="TE"),
            [500.0],
            boundary_points=64,
        )
        p = scan.provenance
        assert p.method == "integral"
        assert p.truncation == 64
        assert p.converged is None
        assert p.notes["boundary_condition"] == "perfectly-conducting"
        assert p.notes["spectral_terms"] >= 32
        assert p.notes["energy_balance_deviation"] < 1e-3
        assert p.wall_time_s is not None

    def test_a_named_coating_is_ignored_with_a_warning(self):
        problem = Problem(
            period=600.0, profile=Sinusoidal(depth_fraction=0.3), coating="Au"
        )
        scan = integral.solve(
            problem,
            Illumination.classical(alpha=25.0, polarization="TE"),
            [500.0],
            boundary_points=64,
        )
        assert any("not consulted" in w for w in scan.provenance.warnings)

    def test_corner_deficit_is_reported_never_rescaled(self):
        scan = integral.solve(
            OFFPLANE,
            Illumination.offplane(graze=1.5, azimuth=25.0, polarization="TM"),
            [2.5],
            boundary_points=96,
        )
        deviation = scan.provenance.notes["energy_balance_deviation"]
        assert deviation > 0.001  # the corner deficit is real at this P
        assert any("strays from unity" in w for w in scan.provenance.warnings)
        # and the numbers themselves were not nudged back to 1:
        assert abs(scan.total[0] - 1.0) == pytest.approx(deviation, rel=1e-6)

    def test_evanescent_orders_are_exactly_zero(self):
        scan = integral.solve(
            SINUSOID,
            Illumination.classical(alpha=25.0, polarization="TE"),
            [450.0, 620.0],
            boundary_points=64,
        )
        assert (scan.efficiency[~scan.propagating] == 0.0).all()


class Recorder:
    def __init__(self, raise_at=None):
        self.calls = []
        self._raise_at = raise_at

    def __call__(self, done, total):
        self.calls.append((done, total))
        if self._raise_at is not None and done >= self._raise_at:
            raise SolveCancelled(f"stopped at {done}")


class TestProgressContract:
    WAVELENGTHS = np.linspace(430.0, 560.0, 8)

    def solve(self, **kwargs):
        return integral.solve(
            SINUSOID,
            Illumination.classical(alpha=25.0, polarization="TE"),
            self.WAVELENGTHS,
            boundary_points=64,
            **kwargs,
        )

    def test_first_call_is_zero_of_total_and_last_is_total_of_total(self):
        recorder = Recorder()
        self.solve(progress=recorder)
        n = len(self.WAVELENGTHS)
        assert recorder.calls[0] == (0, n)
        assert recorder.calls[-1] == (n, n)
        done = [d for d, _ in recorder.calls]
        assert done == sorted(done)

    def test_raising_stops_the_solve_where_it_was_raised(self):
        recorder = Recorder(raise_at=3)
        with pytest.raises(SolveCancelled):
            self.solve(progress=recorder)
        assert recorder.calls[-1][0] == 3


class TestConvergenceHarness:
    def test_ladder_reaches_a_plateau(self):
        """An explicit small ladder on the Littrow sinusoid: the harness must
        find a plateau and stamp the scan converged."""
        report = check_convergence(
            integral,
            SINUSOID,
            Illumination.classical(alpha=30.0, polarization="TE"),
            [600.0],
            values=(48, 64, 96, 128, 192),
            tolerance=1e-3,
        )
        assert report.converged_at is not None
        assert report.scan.provenance.converged is True
