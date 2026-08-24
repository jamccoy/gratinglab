"""Scalar against the integral method: agreement where both are valid, and
measured divergence everywhere else.

This is the committed version of the cross-check the project exists to make
possible. The closed-form tests in ``test_scalar.py`` verify the scalar
integral against answers derived from the same transmittance picture; the
integral solver is validated against PCGrate (``test_integral_corpus.py``,
which needs the external corpus). Here the two live solvers meet with no
external data at all, so this file runs everywhere and skips nothing.

Ground rules that keep the comparison honest:

- **Perfect-conductor terms only.** The integral solver is relative to a
  perfect reflector, so every problem here is uncoated and scalar is relative
  too -- apples to apples, per ``provenance.notes["normalization"]``.
- **Smooth profiles, TE.** Corners limit the integral solver's TM solve to
  first-order convergence, and a comparison must not attribute the integral
  method's mesh error to scalar theory. Sinusoids have no corners; TE is
  clean either way; and every case asserts the integral solver's own
  ``energy_balance_deviation`` is orders of magnitude below the tolerance it
  is being trusted to adjudicate.
- **Divergence is asserted, not just agreement.** Scalar theory has a regime,
  and the tests that show the discrepancy *growing* as the regime is left are
  what make the in-regime agreement information rather than tolerance slack.
  The failure axes measured here: groove depth (dominant), lambda/period, and
  -- the one the naive validity ratio misses -- the half-cone angle, because
  the conical problem decouples at the reduced wavelength ``lambda/sin(gamma)``
  (see ``solvers/integral/__init__.py`` and M&P eq. 4.65).

Measured values quoted in the assertions were taken on this machine at the
stated knobs; tolerances sit a factor of a few above them.
"""

import numpy as np

from gratinglab.compare import align, sweep
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Sinusoidal


def compare_methods(problem, illumination, wavelengths, boundary_points):
    """Both solvers on one problem; returns (comparison, summary, integral_ebd).

    ``integral_ebd`` is the integral solver's own energy-balance deviation --
    its discretisation-error estimate, which every test checks before trusting
    the integral answer as the reference.
    """
    scans = sweep(
        problem,
        illumination,
        wavelengths,
        ["scalar", "integral"],
        {"integral": {"boundary_points": boundary_points}},
    )
    comparison = align(scans)
    return (
        comparison,
        comparison.summary("scalar", "integral"),
        float(scans[1].provenance.notes["energy_balance_deviation"]),
    )


class TestAgreementInTheOverlapRegime:
    """Shallow sinusoid, small lambda/period: both theories are in regime and
    must agree. One case per mount -- an in-plane case is blind to every
    ``sin(gamma)`` in the scalar phase, so the off-plane case is not optional.
    """

    def test_in_plane(self):
        problem = Problem(period=500.0, profile=Sinusoidal(depth_fraction=0.01))
        illumination = Illumination.classical(alpha=10.0, polarization="TE")

        comparison, summary, integral_ebd = compare_methods(
            problem, illumination, [20.0, 25.0], boundary_points=192
        )

        # The integral answer is only a reference if its own discretisation
        # error is far below the disagreement being measured. Measured 5.1e-8.
        assert integral_ebd < 1e-6
        # Measured 5.4e-4 (at lambda=25, order 0); tolerance ~4x above.
        assert summary["max_abs_difference"] < 2e-3

        # Orderwise, not just the headline: every order the integral method
        # puts real power into, scalar reproduces to a few percent relative.
        # Measured worst 0.88% relative on orders above 1e-3.
        scalar_eff, integral_eff = comparison.efficiency
        strong = integral_eff > 1e-3
        relative = np.abs(scalar_eff[strong] - integral_eff[strong])
        assert (relative / integral_eff[strong]).max() < 0.05

    def test_off_plane(self):
        problem = Problem(period=2000.0, profile=Sinusoidal(depth_fraction=0.01))
        illumination = Illumination.offplane(
            graze=2.5, azimuth=20.0, polarization="TE"
        )

        comparison, summary, integral_ebd = compare_methods(
            problem, illumination, [2.0], boundary_points=300
        )

        assert integral_ebd < 1e-6  # measured 2.4e-8
        # Reduced ratio lambda/(p sin gamma) = 0.023 -- in regime. Measured
        # max abs 2.4e-4.
        assert summary["max_abs_difference"] < 1e-3

        scalar_eff, integral_eff = comparison.efficiency
        strong = integral_eff > 1e-3
        relative = np.abs(scalar_eff[strong] - integral_eff[strong])
        assert (relative / integral_eff[strong]).max() < 0.05  # measured 0.49%


class TestDivergenceWithGrooveDepth:
    """Groove depth is scalar theory's dominant failure axis.

    The thin-element phase exp(i k g(t)(cos a + cos b)) is a single-scatter
    picture; as the phase excursion grows the true field develops multiple
    scattering that no transmittance function carries. Fixed lambda/p = 0.04,
    in-plane, so depth is the only thing moving.
    """

    def test_discrepancy_grows_monotonically_with_depth(self):
        illumination = Illumination.classical(alpha=10.0, polarization="TE")

        worst = []
        for depth_fraction in (0.01, 0.02, 0.05, 0.2):
            problem = Problem(
                period=500.0, profile=Sinusoidal(depth_fraction=depth_fraction)
            )
            _, summary, integral_ebd = compare_methods(
                problem, illumination, [20.0], boundary_points=256
            )
            assert integral_ebd < 1e-6  # measured <= 3.0e-7 on every rung
            worst.append(summary["max_abs_difference"])

        # Measured ladder: 3.6e-4, 8.4e-4, 4.4e-3, 1.03e-1.
        assert all(a < b for a, b in zip(worst, worst[1:]))
        # A 20x depth increase costs ~285x in agreement: the failure is real
        # and steep, not tolerance noise.
        assert worst[-1] / worst[0] > 50


class TestDivergenceWithWavelengthOverPeriod:
    """The textbook validity condition, lambda << period, measured directly.

    Shallow profile (depth_fraction = 0.02) so the depth axis stays quiet and
    the wavelength is the only thing moving.
    """

    def test_discrepancy_grows_with_the_ratio(self):
        illumination = Illumination.classical(alpha=10.0, polarization="TE")
        problem = Problem(period=500.0, profile=Sinusoidal(depth_fraction=0.02))

        worst = []
        for ratio in (0.05, 0.1, 0.4):
            _, summary, integral_ebd = compare_methods(
                problem, illumination, [500.0 * ratio], boundary_points=256
            )
            assert integral_ebd < 1e-6
            worst.append(summary["max_abs_difference"])

        # Measured: 8.4e-4, 2.2e-3, 4.1e-3.
        assert all(a <= b for a, b in zip(worst, worst[1:]))
        assert worst[-1] / worst[0] > 3


class TestDivergenceAsTheConeCloses:
    """Scalar validity in a conical mount follows the *reduced* wavelength.

    The perfect-conductor conical problem decouples exactly into an in-plane
    problem at ``lambda / sin(gamma)`` (M&P eq. 4.65; the integral solver is
    built on it), so "wavelength small against the structure" must be judged
    on ``lambda/(p sin gamma)``. This ladder holds lambda/p fixed at 0.001 and
    closes the cone; the discrepancy grows anyway, which is the evidence that
    the naive lambda/p ratio is the wrong guard in an off-plane mount. The
    scalar solver's provenance warning tests the reduced ratio for exactly
    this reason.
    """

    def test_discrepancy_grows_as_gamma_shrinks_at_fixed_lambda_over_p(self):
        problem = Problem(period=1000.0, profile=Sinusoidal(depth_fraction=0.05))

        worst = []
        # boundary_points per rung is the smallest the solver accepts, and it
        # *relaxes* as gamma shrinks: the integral method's cost follows the
        # reduced wavelength, so the extreme off-plane case is its easiest.
        # The scalar solver moves the other way -- that opposition is the
        # finding.
        for gamma, boundary_points in ((5.0, 600), (2.5, 300), (1.25, 192)):
            illumination = Illumination.offplane(
                graze=gamma, azimuth=10.0, polarization="TE"
            )
            _, summary, integral_ebd = compare_methods(
                problem, illumination, [1.0], boundary_points=boundary_points
            )
            assert integral_ebd < 1e-6
            worst.append(summary["max_abs_difference"])

        # Measured: 2.4e-3 (gamma=5), 3.4e-3 (2.5), 4.5e-3 (1.25) -- growing
        # while lambda/p never moved from 0.001.
        assert all(a < b for a, b in zip(worst, worst[1:]))
