"""The corner-graded mesh: geometry, operators, and the solver option.

The thresholds here are the measured ladders of ``docs/theory/integral.md``
section 5 with a margin, so a regression that returns the first-order
corner behaviour fails them by an order of magnitude or more.
"""

import warnings

import numpy as np
import pytest
from scipy.integrate import IntegrationWarning, quad

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, FromProfileData, Lamellar, Sinusoidal
from gratinglab.solvers import integral
from gratinglab.solvers.integral._boundary import physical_boundary
from gratinglab.solvers.integral._finite import solve_finite_states
from gratinglab.solvers.integral._nystrom import (
    corner_subtraction,
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

    def test_graded_boundary_carries_its_polygon_in_nm(self):
        boundary = physical_boundary(LAMELLAR, 64, grading=4)
        vx, vy = boundary.polygon
        np.testing.assert_allclose(vx, [0.0, 500.0, 500.0, 1000.0])
        np.testing.assert_allclose(vy, [300.0, 300.0, 0.0, 0.0])
        assert physical_boundary(LAMELLAR, 64).polygon is None


def _laplace_double_layer(boundary, i, period):
    """What ``rule + corner_subtraction`` should be: the free-space Laplace
    double layer of a constant at node ``i``, over the polygon and its two
    neighbouring periods, by adaptive quadrature that knows nothing of
    subtended angles. The target's own segment is collinear and adds 0."""
    vx, vy = boundary.polygon
    big_x = np.concatenate((vx - period, vx, vx + period, [vx[0] + 2 * period]))
    big_y = np.concatenate((vy, vy, vy, [vy[0]]))
    x0, y0 = boundary.x[i], boundary.y[i]
    total = 0.0
    for ax, ay, bx, by in zip(big_x[:-1], big_y[:-1], big_x[1:], big_y[1:]):
        length = np.hypot(bx - ax, by - ay)
        tx, ty = (bx - ax) / length, (by - ay) / length
        offset = (x0 - ax) * -ty + (y0 - ay) * tx
        if abs(offset) < 1e-9:
            continue
        foot = float(np.clip((x0 - ax) * tx + (y0 - ay) * ty, 0.0, length))

        def kernel(s, ax=ax, ay=ay, tx=tx, ty=ty):
            dx, dy = ax + s * tx - x0, ay + s * ty - y0
            return (-ty * dx + tx * dy) / (2 * np.pi * (dx * dx + dy * dy))

        # Decade breakpoints either side of the foot: the kernel is a peak
        # of width |offset|, which reaches 1e-9 nm next to a corner.
        steps = abs(offset) * 10.0 ** np.arange(-2, 13)
        breaks = np.unique(np.clip(np.concatenate(
            ([0.0, foot, length], foot - steps, foot + steps)), 0.0, length))
        # quad reports roundoff on panels where the kernel nearly cancels;
        # the assertion on the total is the check, not its error estimate.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            total += sum(
                quad(kernel, lo, hi, epsabs=1e-13, epsrel=1e-10)[0]
                for lo, hi in zip(breaks[:-1], breaks[1:])
            )
    return total


class TestOperators:
    @pytest.mark.parametrize("problem", [LAMELLAR, SAWTOOTH, OFFPLANE])
    def test_corner_subtraction_completes_the_rule(self, problem):
        """At the nodes nearest the corners, where the graded rule is wrong
        by up to a sixth of the identity term, the rule plus the subtraction
        equals an independent quadrature. The blaze is the profile whose
        node coordinates carry rounding, which once made the target's own
        segment count as a whole 1/2."""
        period = problem.period
        boundary = physical_boundary(problem, 128, grading=4)
        correction = corner_subtraction(boundary, period=period)
        near = np.argsort(boundary.quadrature)[:12]
        for i in near:
            rule = 0.0
            for image in (-1, 0, 1):
                dx = boundary.x + image * period - boundary.x[i]
                dy = boundary.y - boundary.y[i]
                with np.errstate(divide="ignore", invalid="ignore"):
                    k0 = (boundary.nx * dx + boundary.ny * dy) / (
                        2 * np.pi * (dx * dx + dy * dy)
                    )
                if image == 0:
                    k0[i] = 0.0
                rule += k0 @ boundary.quadrature
            reference = _laplace_double_layer(boundary, i, period)
            assert rule + correction[i] == pytest.approx(reference, abs=1e-8)
        assert np.abs(correction[near]).max() > 0.05

    def test_corner_subtraction_is_local_and_graded_only(self):
        boundary = physical_boundary(LAMELLAR, 256, grading=4)
        correction = corner_subtraction(boundary, period=LAMELLAR.period)
        # Measured: largest 0.16, median 1.6e-7.
        assert np.median(np.abs(correction)) < 1e-5
        plain = physical_boundary(LAMELLAR, 256)
        assert not corner_subtraction(plain, period=LAMELLAR.period).any()

    def test_tangential_layer_is_the_derivative_of_the_trace(self):
        """``D_t V sigma = d/ds (V sigma)`` along each facet, at a metal-side
        wavenumber. The equal-arc mesh cannot pass this on a cornered
        profile at all (its spectral derivative of the trace is off by order
        one), which is why the check is on the graded mesh. The two agree on
        smooth data; on a solved corner density they part near the corners,
        and the graded solve uses the composition (``_finite``)."""
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

    def test_tabulated_conical_xray_lamellar_converges(self):
        """The case the graded mesh was once refused over: Au at 2 nm, 1.5
        degree graze, a shallow lamellar whose walls are a tenth of its
        floors. Measured 3.7e-5 graded and 0.24 on the equal-arc mesh."""
        problem = Problem(
            period=315.15, profile=Lamellar(depth_fraction=0.05), coating="Au"
        )
        light = Illumination.offplane(graze=1.5, azimuth=25.0)
        kwargs = dict(boundary_points=256, conductivity="tabulated")
        graded = integral.solve(problem, light, [2.0], corner_grading=4, **kwargs)
        plain = integral.solve(problem, light, [2.0], **kwargs)
        assert graded.provenance.notes["energy_balance_deviation"] < 3e-4
        assert plain.provenance.notes["energy_balance_deviation"] > 0.1
        assert graded.provenance.notes["mesh"] == "corner-graded, order 4"


def _finite_deviation(problem, light, *, index, wavelength, points, grading):
    """``|R + A - 1|`` of the coupled system at a fixed index, TM incident."""
    boundary = physical_boundary(problem, points, grading=grading)
    cos_gamma = float(np.sqrt(1.0 - light.sin_gamma**2))
    (state,) = solve_finite_states(
        boundary, wavelength=wavelength, period=problem.period,
        sin_alpha=light.sin_alpha, index=index, cos_gamma=cos_gamma,
        incidents=[(0.0, 1.0)], terms=points // 2,
    )
    return abs(state.total + state.absorption - 1.0)


class TestFiniteConductivity:
    """The two corrections the coupled system needs on a graded mesh, each
    pinned by the case that failed without it."""

    def test_visible_metal_tm_converges(self):
        """Without the double layer's corner subtraction this stalled near
        3e-3 at every refinement. Measured 7.1e-3 equal-arc, 1.9e-5 graded."""
        problem = Problem(period=1000.0, profile=Lamellar(depth_fraction=0.2))
        light = Illumination.classical(alpha=10.0, polarization="TM")
        kwargs = dict(index=0.2 + 3.5j, wavelength=600.0, points=256)
        assert _finite_deviation(problem, light, grading=0, **kwargs) > 5e-3
        assert _finite_deviation(problem, light, grading=4, **kwargs) < 1e-4

    def test_conical_order_five_converges(self):
        """With the single ``D_t V`` operator this stalled near 7e-3 from
        P = 256 to 768. With the composition ``D_t (V w)`` it is 2.1e-6 at
        P = 512 (measured)."""
        light = Illumination.offplane(graze=1.5, azimuth=25.0)
        problem = Problem(period=315.15, profile=Lamellar(depth_fraction=0.05))
        deviation = _finite_deviation(
            problem, light, index=0.998 + 0.002j, wavelength=2.0,
            points=512, grading=5,
        )
        assert deviation < 2e-5
