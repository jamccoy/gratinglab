"""How the rigorous solvers run: the BLAS cap, workers, and RCWA's shortcuts.

None of this may change a number. The BLAS cap must hold inside a solve and
be gone after it, however the solve ends; ``workers`` must give the serial
result exactly and keep the progress contract of ``solvers/base.py`` (calls on
the caller's thread, ``(0, n)`` first, cancellation unchanged, nothing left
running); and RCWA's decoupled in-plane solve and top-plane-only flux must
agree with the coupled solve and the full set of planes they skip.
"""

import threading

import numpy as np
import pytest
from threadpoolctl import threadpool_info

from gratinglab.illumination import Illumination
from gratinglab.materials.optical import OpticalConstants, write_ari
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, Lamellar, Sinusoidal
from gratinglab.solvers import SolveCancelled, get_solver
from gratinglab.solvers._execution import capped_blas
from gratinglab.solvers.rcwa._core import incident_field, solve_stack
from gratinglab.solvers.rcwa._stack import Slab, Stack, groove_geometry

WAVELENGTHS = np.array([420.0, 470.0, 530.0, 590.0, 650.0])


@pytest.fixture
def library(tmp_path, monkeypatch):
    for name, index in (("glass", 1.5 + 0j), ("metal", 0.2 + 3.0j)):
        write_ari(
            OpticalConstants(
                name=name,
                wavelength_nm=(200.0, 1000.0),
                decrement=(1.0 - index.real,) * 2,
                absorption=(index.imag,) * 2,
                source="test_solver_execution",
            ),
            tmp_path / f"{name}.ari",
        )
    monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))


CASES = {
    "integral": (
        Problem(period=600.0, profile=Blazed(blaze_angle=20.0, antiblaze_angle=70.0)),
        Illumination.classical(alpha=20.0, polarization="unpolarized"),
        {"boundary_points": 96},
    ),
    "integral-tabulated": (
        Problem(
            period=600.0,
            profile=Blazed(blaze_angle=20.0, antiblaze_angle=70.0),
            coating="metal",
        ),
        Illumination.offplane(graze=40.0, azimuth=25.0, polarization="unpolarized"),
        {"boundary_points": 128, "conductivity": "tabulated"},
    ),
    "rcwa": (
        Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3), coating="metal"),
        Illumination.offplane(graze=40.0, azimuth=25.0, polarization="unpolarized"),
        {"truncation": 6, "slices": 6},
    ),
    "rcwa-in-plane": (
        Problem(period=600.0, profile=Lamellar(depth_fraction=0.3), coating="metal"),
        Illumination.classical(alpha=20.0, polarization="unpolarized"),
        {"truncation": 8},
    ),
}


@pytest.fixture(params=sorted(CASES))
def case(request, library):
    name = request.param
    problem, illumination, options = CASES[name]
    solver = get_solver(name.split("-")[0])

    def solve(**extra):
        return solver.solve(
            problem, illumination, WAVELENGTHS, **{**options, **extra}
        )

    return solve


def blas_limits():
    return sorted(
        (info["filepath"], info["num_threads"])
        for info in threadpool_info()
        if info["user_api"] == "blas"
    )


class TestTheBlasCap:
    def test_it_holds_during_the_solve_and_is_gone_after(self, case):
        before = blas_limits()
        seen = []
        case(progress=lambda done, total: seen.append(blas_limits()))
        assert seen and all(
            threads == 1 for limits in seen for _, threads in limits
        )
        assert blas_limits() == before

    def test_it_is_restored_when_the_solve_is_cancelled(self, case):
        before = blas_limits()

        def cancel(done, total):
            if done == 1:
                raise SolveCancelled

        with pytest.raises(SolveCancelled):
            case(progress=cancel)
        assert blas_limits() == before

    def test_none_leaves_the_pools_alone(self, case):
        before = blas_limits()
        seen = []
        case(blas_threads=None, progress=lambda d, t: seen.append(blas_limits()))
        assert all(limits == before for limits in seen)

    def test_overlapping_caps_restore_only_when_the_last_one_leaves(self):
        before = blas_limits()
        with capped_blas(1):
            with capped_blas(1):
                pass
            assert all(threads == 1 for _, threads in blas_limits())
        assert blas_limits() == before

    @pytest.mark.parametrize("value", [0, -1, 1.5, True, "2"])
    def test_a_bad_value_is_refused(self, case, value):
        with pytest.raises(ValueError, match="blas_threads"):
            case(blas_threads=value)


class TestWorkers:
    def test_the_scan_is_the_serial_scan_exactly(self, case):
        serial = case()
        threaded = case(workers=3)
        np.testing.assert_array_equal(threaded.efficiency, serial.efficiency)
        np.testing.assert_array_equal(threaded.propagating, serial.propagating)
        if serial.absorption is not None:
            np.testing.assert_array_equal(threaded.absorption, serial.absorption)
        assert threaded.provenance.warnings == serial.provenance.warnings

    def test_progress_keeps_the_contract_on_the_callers_thread(self, case):
        caller = threading.get_ident()
        calls, threads = [], set()

        def record(done, total):
            calls.append((done, total))
            threads.add(threading.get_ident())

        case(workers=3, progress=record)
        n = len(WAVELENGTHS)
        assert calls == [(done, n) for done in range(n + 1)]
        assert threads == {caller}

    def test_cancelling_propagates_unchanged_and_leaves_nothing_running(self, case):
        raised = SolveCancelled("stop")
        calls = []

        def cancel(done, total):
            calls.append(done)
            if done == 1:
                raise raised

        with pytest.raises(SolveCancelled) as caught:
            case(workers=3, progress=cancel)
        assert caught.value is raised
        assert calls == [0, 1]
        assert not [t for t in threading.enumerate() if t.name.startswith("gratinglab-row")]

    def test_cancelling_before_any_work_starts_none(self, case):
        def cancel(done, total):
            raise SolveCancelled

        with pytest.raises(SolveCancelled):
            case(workers=3, progress=cancel)

    @pytest.mark.parametrize("value", [0, -2, 2.0, True, None])
    def test_a_bad_value_is_refused(self, case, value):
        with pytest.raises(ValueError, match="workers"):
            case(workers=value)


def in_plane_stack(index=0.2 + 3.0j):
    groove = groove_geometry(Sinusoidal(depth_fraction=0.3), 600.0, 5)
    film = Slab(thickness=40.0, background=1.3 + 0.01j)
    stack = groove.stack(index)
    return Stack(
        superstrate=stack.superstrate,
        grating=stack.grating,
        substrate=1.5,
        below=(film,),
    )


def solve(stack, *, polarizations, gamma=90.0, **options):
    a, g = np.radians(20.0), np.radians(gamma)
    kx0, ky, kz0 = -np.sin(a) * np.sin(g), np.cos(g), np.cos(a) * np.sin(g)
    return solve_stack(
        stack,
        wavelength=500.0,
        period=600.0,
        kx0=kx0,
        ky=ky,
        incidents=[incident_field(p, kx0, ky, kz0) for p in polarizations],
        truncation=7,
        **options,
    )


class TestRcwaShortcuts:
    @pytest.mark.parametrize("polarizations", [("TE",), ("TM",), ("TE", "TM")])
    @pytest.mark.parametrize("factorization", ["li", "laurent"])
    def test_the_decoupled_in_plane_solve_is_the_coupled_one(
        self, polarizations, factorization
    ):
        stack = in_plane_stack()
        split = solve(stack, polarizations=polarizations, factorization=factorization)
        whole = solve(
            stack,
            polarizations=polarizations,
            factorization=factorization,
            decouple=False,
        )
        for field in ("reflected", "transmitted", "plane_flux"):
            np.testing.assert_allclose(
                getattr(split, field), getattr(whole, field), rtol=0, atol=1e-12
            )

    @pytest.mark.parametrize("gamma", [90.0, 60.0])
    def test_the_top_plane_alone_is_the_first_of_all_planes(self, gamma):
        stack = in_plane_stack()
        every = solve(stack, polarizations=("TE", "TM"), gamma=gamma)
        top = solve(stack, polarizations=("TE", "TM"), gamma=gamma, planes="top")
        assert every.plane_flux.shape == (2, len(stack.slabs) + 1)
        assert top.plane_flux.shape == (2, 1)
        np.testing.assert_allclose(top.absorbed, every.absorbed, rtol=0, atol=1e-13)
        np.testing.assert_allclose(top.reflected, every.reflected, rtol=0, atol=1e-13)
        np.testing.assert_allclose(
            top.transmitted, every.transmitted, rtol=0, atol=1e-13
        )

    def test_decoupling_an_off_plane_mount_is_refused(self):
        with pytest.raises(ValueError, match="in-plane"):
            solve(in_plane_stack(), polarizations=("TE",), gamma=60.0, decouple=True)

    def test_an_unknown_plane_set_is_refused(self):
        with pytest.raises(ValueError, match="planes"):
            solve(in_plane_stack(), polarizations=("TE",), planes="bottom")
