"""Both rigorous solvers, swept through a Rayleigh anomaly.

Where an order passes off, its ``kz`` goes to zero: RCWA's external region
divides by it (``V = Q kz^-1``) and the integral kernel by ``gamma_m``. A
solver can be fooled there into an infinite, a NaN, or a step that is the
method's and not the physics'. The physics is a square-root branch point:
every efficiency is continuous through the anomaly, with a slope that is
infinite on the side where the order propagates.

So the sweep approaches ``lambda*`` from both sides, to ``1e-8`` nm, and
lands on it. The gap between the rows at ``lambda* - d`` and ``lambda* + d``
should then be ``C sqrt(d / lambda*)`` with ``C`` constant, and so should the
new order's efficiency. A real step would instead make ``C`` grow as
``d^-1/2``, a thousandfold across this sweep. Measured 2026-10-27 on the
grating below, ``C`` for the largest gap over the orders present on both
sides, with the new order's in brackets:

==================  ===========  ===========
backend             TE           TM
==================  ===========  ===========
rcwa (M = 20)       0.96 (0.99)  0.86 (1.09)
integral, tabulated 0.96 (0.99)  0.95 (1.12)
integral, perfect   0.61 (0.61)  12.6 (12.6)
==================  ===========  ===========

each constant to 4% from ``d = 1e-2`` to ``1e-8`` nm. RCWA and the
tabulated integral method agree on the coefficient; the perfect conductor in
TM has a much stronger anomaly, which is the physics (Wood's), not a fault.
Landing exactly on ``lambda*``, a point each backend floors (RCWA's ``kz``,
the kernel's ``gamma_m``), gives a value within ``2.2e-5`` of its
neighbours' mean, half their own spread. Energy balance holds at every row,
at each solver's usual accuracy. The tolerances below are set from these.
"""

import numpy as np
import pytest

from gratinglab.geometry import sin_beta
from gratinglab.illumination import Illumination
from gratinglab.materials.optical import OpticalConstants, write_ari
from gratinglab.problem import Problem
from gratinglab.profiles import Lamellar
from gratinglab.solvers import get_solver

METAL = 0.2 + 3.0j
PERIOD = 600.0
ALPHA = 20.0
#: An order passes off where ``|sin beta_m| = 1``; for this mount the first
#: is order ``+1`` (in ``geometry``'s sign convention) at
#: ``period (1 + sin alpha)``, which the fixture below checks rather than
#: trusts.
GRAZING = PERIOD * (1.0 + np.sin(np.radians(ALPHA)))
OFFSETS = (1e-2, 1e-4, 1e-6, 1e-8)
#: Largest branch coefficient ``C`` accepted; measured 0.6 to 12.6.
BRANCH_LIMIT = 20.0
#: How far ``C`` may drift across the sweep, as max / min; measured 1.04. A
#: step in place of a branch point would drift by 1000.
BRANCH_DRIFT = 1.25


def wavelengths():
    below = [GRAZING - d for d in OFFSETS]
    above = [GRAZING + d for d in reversed(OFFSETS)]
    return np.array(below + [GRAZING] + above)


CASES = {
    "rcwa": dict(solver="rcwa", options=dict(truncation=20), balance=1e-10),
    "integral-perfect": dict(
        solver="integral",
        options=dict(boundary_points=400, corner_grading=4),
        balance=1e-5,
    ),
    "integral-tabulated": dict(
        solver="integral",
        options=dict(boundary_points=400, corner_grading=4, conductivity="tabulated"),
        balance=1e-5,
    ),
}


@pytest.fixture(scope="module", params=list(CASES))
def sweep(request, tmp_path_factory):
    # A module-scoped copy of the library fixture: one solve per backend
    # serves every assertion below.
    case = CASES[request.param]
    directory = tmp_path_factory.mktemp(request.param)
    write_ari(
        OpticalConstants(
            name="metal",
            wavelength_nm=(200.0, 1000.0),
            decrement=(1.0 - METAL.real,) * 2,
            absorption=(METAL.imag,) * 2,
            source="test_rayleigh_anomaly",
        ),
        directory / "metal.ari",
    )
    patch = pytest.MonkeyPatch()
    patch.setenv("GRATINGLAB_MATERIALS", str(directory))
    try:
        scan = get_solver(case["solver"]).solve(
            Problem(
                period=PERIOD,
                profile=Lamellar(depth_fraction=0.3, duty_cycle=0.4),
                coating="metal",
            ),
            Illumination.classical(alpha=ALPHA, polarization="unpolarized"),
            wavelengths(),
            **case["options"],
        )
    finally:
        patch.undo()
    return request.param, case, scan


def test_the_sweep_straddles_an_order_passing_off():
    sines = sin_beta(np.array([1]), GRAZING, PERIOD, np.sin(np.radians(ALPHA)))
    assert abs(abs(float(sines[0])) - 1.0) < 1e-12


def test_every_efficiency_is_finite(sweep):
    _, _, scan = sweep
    assert np.isfinite(scan.efficiency).all()
    assert (scan.efficiency[scan.propagating] >= -1e-12).all()
    if scan.absorption is not None:
        assert np.isfinite(scan.absorption).all()


def test_energy_balance_holds_at_every_row(sweep):
    name, case, scan = sweep
    totals = np.where(scan.propagating, scan.efficiency, 0.0).sum(axis=1)
    if scan.absorption is not None:
        totals = totals + scan.absorption
    assert np.abs(totals - 1.0).max() < case["balance"], name


def _coefficients(values):
    """``value / sqrt(d / lambda*)`` for each offset."""
    return np.asarray(values) / np.sqrt(np.asarray(OFFSETS) / GRAZING)


def _is_a_branch(coefficients):
    return (
        coefficients.max() <= BRANCH_LIMIT
        and coefficients.max() / coefficients.min() <= BRANCH_DRIFT
    )


def test_the_order_comes_on_at_the_anomaly(sweep):
    """The grazing order propagates below ``lambda*`` and not above, and its
    efficiency vanishes as it goes, like ``sqrt(d)``."""
    name, _, scan = sweep
    column = int(np.flatnonzero(scan.orders == 1)[0])
    count = len(OFFSETS)
    assert scan.propagating[:count, column].all()
    assert not scan.propagating[count + 1 :, column].any()
    coefficients = _coefficients(scan.efficiency[:count, column])
    assert _is_a_branch(coefficients), (name, coefficients)


def test_efficiencies_are_continuous_through_it(sweep):
    name, _, scan = sweep
    count = len(OFFSETS)
    both = scan.propagating[0] & scan.propagating[-1]
    assert both.sum() >= 1
    gaps = [
        float(np.abs(scan.efficiency[-1 - step, both] - scan.efficiency[step, both]).max())
        for step in range(count)
    ]
    coefficients = _coefficients(gaps)
    assert _is_a_branch(coefficients), (name, coefficients)
    # Exactly on it: with its neighbours, not a value of its own.
    exact = scan.efficiency[count, both]
    nearest = 0.5 * (scan.efficiency[count - 1, both] + scan.efficiency[count + 1, both])
    assert float(np.abs(exact - nearest).max()) <= gaps[-1], name


def test_the_anomaly_is_flagged(sweep):
    _, _, scan = sweep
    assert any("Rayleigh anomaly" in w for w in scan.provenance.warnings)
