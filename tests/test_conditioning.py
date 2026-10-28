"""Condition estimates of the rigorous solvers' dense systems.

``solvers/_conditioning.py`` replaces ``numpy.linalg.solve`` and ``inv`` in
RCWA and the integral method with a factorisation that also yields LAPACK's
``gecon`` estimate. Checked here: that the solves are unchanged, that the
estimate is the 1-norm condition number to a small factor, and that each
backend records it and warns only past its measured threshold.
"""

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.materials.optical import OpticalConstants, write_ari
from gratinglab.problem import Problem
from gratinglab.profiles import Lamellar
from gratinglab.solvers import get_solver
from gratinglab.solvers._conditioning import ConditionLog, inverse, solve
from gratinglab.solvers.integral import _WARN_KAPPA as INTEGRAL_WARN
from gratinglab.solvers.rcwa import _WARN_KAPPA as RCWA_WARN

METAL = 0.2 + 3.0j
PERIOD = 600.0
ALPHA = 20.0
GRAZING = PERIOD * (1.0 + np.sin(np.radians(ALPHA)))


def hilbert(n):
    i = np.arange(n)
    return 1.0 / (i[:, None] + i[None, :] + 1.0)


@pytest.fixture
def library(tmp_path, monkeypatch):
    write_ari(
        OpticalConstants(
            name="metal",
            wavelength_nm=(200.0, 1000.0),
            decrement=(1.0 - METAL.real,) * 2,
            absorption=(METAL.imag,) * 2,
            source="test_conditioning",
        ),
        tmp_path / "metal.ari",
    )
    monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))
    return tmp_path


def problem():
    return Problem(
        period=PERIOD,
        profile=Lamellar(depth_fraction=0.3, duty_cycle=0.4),
        coating="metal",
    )


class TestTheSolve:
    @pytest.mark.parametrize("dtype", [np.float64, np.complex128])
    def test_it_is_numpys_solve(self, dtype):
        rng = np.random.default_rng(7)
        a = rng.normal(size=(40, 40)).astype(dtype)
        if dtype is np.complex128:
            a = a + 1j * rng.normal(size=(40, 40))
        b = rng.normal(size=(40, 3)).astype(dtype)
        x, _ = solve(a, b)
        np.testing.assert_allclose(x, np.linalg.solve(a, b), rtol=1e-12, atol=0)
        a_inv, _ = inverse(a)
        np.testing.assert_allclose(a_inv, np.linalg.inv(a), rtol=1e-12, atol=1e-14)
        vector, _ = solve(a, b[:, 0])
        assert vector.shape == (40,)

    @pytest.mark.parametrize(
        "matrix",
        [
            np.random.default_rng(3).normal(size=(60, 60)),
            hilbert(8),
            hilbert(12),
            np.diag([1.0, 1e-9, 1.0]),
        ],
        ids=["random", "hilbert8", "hilbert12", "diagonal"],
    )
    def test_the_estimate_is_the_1_norm_condition(self, matrix):
        _, kappa = solve(matrix, np.ones(len(matrix)))
        exact = np.linalg.cond(matrix, 1)
        # gecon is a lower bound, within a small factor in practice
        # (Higham, Accuracy and Stability, section 15.3); measured here it is
        # exact to rounding for all four.
        assert exact / 3.0 <= kappa <= exact * (1.0 + 1e-6)

    def test_a_singular_matrix_raises_as_numpy_does(self):
        with pytest.raises(np.linalg.LinAlgError):
            solve(np.zeros((3, 3)), np.ones(3))


class TestTheLog:
    def test_it_keeps_the_worst_per_group(self):
        log = ConditionLog()
        log.record("a", 10.0)
        log.record("a", 5.0)
        log.merge({"a": 7.0, "b": 3.0})
        assert log.worst == {"a": 10.0, "b": 3.0}

    def test_it_warns_only_past_the_threshold(self):
        log = ConditionLog()
        log.record("modes", 1e3)
        assert log.warning("rcwa", [500.0], 1e10) is None
        log.record("smatrix", 2e12)
        message = log.warning("rcwa", [500.0, 510.0], 1e10)
        assert "smatrix" in message and "modes" not in message
        assert "2 wavelength(s)" in message and "500-510 nm" in message


class TestTheBackendsRecordIt:
    def test_rcwa(self, library):
        scan = get_solver("rcwa").solve(
            problem(), Illumination.classical(alpha=ALPHA, polarization="TM"),
            [500.0, 700.0], truncation=20,
        )
        condition = scan.provenance.notes["condition"]
        assert condition["warn_above"] == RCWA_WARN
        assert "gecon" in condition["estimator"]
        assert 1.0 <= condition["modes"] < 1e5
        assert 1.0 <= condition["smatrix"] < 1e6
        assert not any("condition estimate" in w for w in scan.provenance.warnings)

    @pytest.mark.parametrize("conductivity", ["perfect", "tabulated"])
    def test_integral(self, library, conductivity):
        scan = get_solver("integral").solve(
            problem(), Illumination.classical(alpha=ALPHA, polarization="TE"),
            [500.0], boundary_points=200, conductivity=conductivity,
        )
        condition = scan.provenance.notes["condition"]
        assert condition["warn_above"] == INTEGRAL_WARN
        assert 1.0 <= condition["nystrom"] < 1e6
        assert not any("condition estimate" in w for w in scan.provenance.warnings)


class TestTheyWarnWhereItIsUnusual:
    """Each threshold is crossed by an order at grazing and not by its
    ordinary neighbour (the measurements are in each backend's
    ``_WARN_KAPPA`` comment and the theory pages)."""

    def test_rcwa_exactly_at_grazing(self, library):
        scan = get_solver("rcwa").solve(
            problem(), Illumination.classical(alpha=ALPHA, polarization="TE"),
            [GRAZING - 1.0, GRAZING], truncation=20,
        )
        assert scan.provenance.notes["condition"]["smatrix"] > RCWA_WARN
        flagged = [w for w in scan.provenance.warnings if "condition estimate" in w]
        assert len(flagged) == 1 and "1 wavelength(s)" in flagged[0]
        assert "smatrix" in flagged[0]

    def test_integral_within_1e_12_nm_of_grazing(self, library):
        scan = get_solver("integral").solve(
            problem(), Illumination.classical(alpha=ALPHA, polarization="TE"),
            [GRAZING - 1.0, GRAZING - 1e-12], boundary_points=400,
            corner_grading=4, conductivity="tabulated",
        )
        assert scan.provenance.notes["condition"]["nystrom"] > INTEGRAL_WARN
        flagged = [w for w in scan.provenance.warnings if "condition estimate" in w]
        assert len(flagged) == 1 and "1 wavelength(s)" in flagged[0]
