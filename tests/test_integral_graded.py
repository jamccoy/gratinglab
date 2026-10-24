"""The corner-graded mesh: geometry, operators, and the solver option.

The thresholds here are the measured ladders of ``docs/theory/integral.md``
section 5 with a margin, so a regression that returns the first-order
corner behaviour fails them by an order of magnitude or more.
"""

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, FromProfileData, Lamellar, Sinusoidal
from gratinglab.solvers import integral
from gratinglab.solvers.base import UnsupportedConfiguration
from gratinglab.solvers.integral._boundary import physical_boundary
from gratinglab.solvers.integral._nystrom import (
    single_layer_matrix,
    tangential_derivative_matrix,
    tangential_layer_matrix,
)

LAMELLAR = Problem(period=1000.0, profile=Lamellar(depth_fraction=0.3))
SAWTOOTH = Problem(period=1000.0, profile=Blazed(blaze_angle=30.0))
OFFPLANE = Problem(
    period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
)


class TestMesh:
    @pytest.mark.parametrize(
        "profile, perimeter",
        [
            (Lamellar(depth_fraction=0.3), 1.6),
            (Blazed(blaze_angle=30.0), 1 / np.cos(np.pi / 6) + np.tan(np.pi / 6)),
        ],
    )
    def test_weights_integrate_the_exact_perimeter(self, profile, perimeter):
        curve = profile.boundary(128, grading=4)
        assert curve.graded
        assert curve.ds.sum() == pytest.approx(perimeter, rel=1e-5)

    def test_points_cluster_at_corners_and_none_sits_on_one(self):
        curve = Lamellar(depth_fraction=0.3).boundary(128, grading=4)
        steps = np.hypot(
            np.diff(np.append(curve.t, curve.t[0] + 1.0)),
            np.diff(np.append(curve.y, curve.y[0])),
        )
        assert steps.min() < 1e-5 < 1e-2 < steps.max()
        assert (curve.ds > 0).all()

    def test_normals_are_exact_unit_facet_normals(self):
        curve = Lamellar(depth_fraction=0.3).boundary(128, grading=4)
        assert np.hypot(curve.nx, curve.ny) == pytest.approx(1.0, abs=1e-15)
        # Every point of a rectangle lies on a horizontal or vertical facet.
        axis = np.minimum(np.abs(curve.nx), np.abs(curve.ny))
        assert axis.max() < 1e-15

    def test_order_runs_from_the_start_of_the_period(self):
        curve = Blazed(blaze_angle=29.5, antiblaze_angle=70.5).boundary(
            96, grading=4
        )
        assert curve.t[0] >= 0.0
        assert curve.t[-1] < curve.t[0] + 1.0

    def test_undercut_is_graded_too(self):
        profile = FromProfileData(
            t=(0.0, 0.30, 0.28, 0.55, 1.0), y=(0.0, 0.45, 0.2, 0.5, 0.0)
        )
        curve = profile.boundary(128, grading=4)
        assert curve.graded
        assert (curve.ny < 0).any()

    def test_no_corner_keeps_the_equal_arc_mesh(self):
        profile = Sinusoidal(depth_fraction=0.3)
        uniform, graded = profile.boundary(128), profile.boundary(128, grading=4)
        assert not graded.graded
        np.testing.assert_array_equal(uniform.t, graded.t)
        np.testing.assert_array_equal(uniform.ds, graded.ds)

    @pytest.mark.parametrize("grading", [1, -2])
    def test_order_one_and_negative_are_refused(self, grading):
        with pytest.raises(ValueError, match="grading"):
            Lamellar(depth_fraction=0.3).boundary(64, grading=grading)

    def test_too_few_points_per_side_is_refused(self):
        with pytest.raises(ValueError, match="sides between corners"):
            Lamellar(depth_fraction=0.3).boundary(12, grading=4)


class TestBoundary:
    def test_equal_arc_weights_are_the_spacing(self):
        boundary = physical_boundary(LAMELLAR, 64)
        assert not boundary.graded
        assert boundary.quadrature == pytest.approx(boundary.spacing)
        assert boundary.speed == pytest.approx(boundary.arc_length)
        assert boundary.resolution_length == boundary.arc_length

    def test_graded_weights_are_carried_in_nm(self):
        boundary = physical_boundary(LAMELLAR, 64, grading=4)
        assert boundary.graded
        assert boundary.quadrature.sum() == pytest.approx(1600.0, rel=1e-4)
        assert boundary.resolution_length > boundary.arc_length


class TestOperators:
    def test_tangential_layer_is_the_derivative_of_the_trace(self):
        """``D_t V sigma = d/ds (V sigma)`` along each facet, at a metal-side
        wavenumber. The equal-arc mesh cannot pass this on a cornered
        profile at all (its spectral derivative of the trace is off by order
        one), which is why the check is on the graded mesh."""
        light = Illumination.offplane(graze=1.5, azimuth=25.0)
        cos_gamma = float(np.sqrt(1.0 - light.sin_gamma**2))
        k = 2.0 * np.pi / 2.0
        km = k * np.sqrt((0.998 + 0.002j) ** 2 - cos_gamma**2)
        alpha0 = -k * light.sin_gamma * light.sin_alpha
        boundary = physical_boundary(OFFPLANE, 256, grading=4)
        kwargs = dict(alpha0=alpha0, period=OFFPLANE.period, terms=128)
        trace = single_layer_matrix(boundary, k=km, **kwargs)
        layer = tangential_layer_matrix(boundary, k=km, **kwargs)
        derivative = tangential_derivative_matrix(boundary, alpha0=alpha0)
        sigma = np.exp(1j * alpha0 * boundary.x) * (
            1.0 + 0.3 * np.cos(2 * np.pi * boundary.x / OFFPLANE.period)
        )
        direct, composed = layer @ sigma, derivative @ (trace @ sigma)
        mid = boundary.quadrature >= np.quantile(boundary.quadrature, 0.5)
        error = np.abs(direct - composed)[mid].max() / np.abs(composed).max()
        assert error < 2e-3


class TestSolver:
    TM = Illumination.classical(alpha=10.0, polarization="TM")

    def test_tm_on_a_lamellar_converges_where_equal_arc_stalls(self):
        uniform = integral.solve(LAMELLAR, self.TM, [600.0], boundary_points=256)
        graded = integral.solve(
            LAMELLAR, self.TM, [600.0], boundary_points=256, corner_grading=4
        )
        # Measured 2.3e-3 and 1.3e-5.
        assert uniform.provenance.notes["energy_balance_deviation"] > 1e-3
        assert graded.provenance.notes["energy_balance_deviation"] < 1e-4

    def test_sawtooth_tm_efficiencies_settle(self):
        light = Illumination.classical(alpha=20.0, polarization="TM")
        coarse, fine = (
            integral.solve(
                SAWTOOTH, light, [500.0], boundary_points=p, corner_grading=4
            )
            for p in (256, 512)
        )
        # Measured 9e-6 apart; the equal-arc mesh moves by 4e-3 over the
        # same step.
        np.testing.assert_allclose(coarse.efficiency, fine.efficiency, atol=1e-4)
        assert fine.provenance.notes["energy_balance_deviation"] < 1e-4

    def test_provenance_names_the_mesh(self):
        graded = integral.solve(
            LAMELLAR, self.TM, [600.0], boundary_points=96, corner_grading=4
        )
        plain = integral.solve(LAMELLAR, self.TM, [600.0], boundary_points=96)
        assert graded.provenance.notes["mesh"] == "corner-graded, order 4"
        assert plain.provenance.notes["mesh"] == "equal arc length"

    def test_smooth_profile_is_unchanged(self):
        problem = Problem(period=1000.0, profile=Sinusoidal(depth_fraction=0.3))
        plain = integral.solve(problem, self.TM, [600.0], boundary_points=96)
        graded = integral.solve(
            problem, self.TM, [600.0], boundary_points=96, corner_grading=4
        )
        np.testing.assert_array_equal(plain.efficiency, graded.efficiency)

    @pytest.mark.parametrize("value", [1, 7, -1, True, 4.0])
    def test_out_of_range_orders_are_refused(self, value):
        with pytest.raises(ValueError, match="corner_grading"):
            integral.solve(
                LAMELLAR, self.TM, [600.0], boundary_points=64,
                corner_grading=value,
            )

    def test_tabulated_conductivity_is_refused(self):
        problem = Problem(
            period=1000.0, profile=Lamellar(depth_fraction=0.3), coating="Au"
        )
        with pytest.raises(UnsupportedConfiguration, match="tabulated"):
            integral.solve(
                problem, self.TM, [600.0], boundary_points=64,
                conductivity="tabulated", corner_grading=4,
            )
