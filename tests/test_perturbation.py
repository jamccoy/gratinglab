r"""The scalar solver against first-order Rayleigh perturbation theory.

**Why this file exists.** Every other analytic check in the suite --
``sinc^2``, the binary-grating form, the sinusoid's ``J_m^2`` -- is derived from
the same transmittance-function picture the solver evaluates. They are excellent
tests of the *quadrature* and blind to the *model*: multiply every efficiency by
any function symmetric in :math:`\alpha \leftrightarrow \beta_m` and all three
still pass, and so does ``check_reciprocity``, which constrains the structure of
:math:`\Phi` but cannot see a prefactor.

Perturbation theory comes from somewhere else entirely. Expand the field in
Rayleigh plane waves, impose the boundary condition on the corrugated surface,
and keep terms of first order in the groove height. For a Dirichlet (perfectly
conducting, TE) surface that gives

.. math::
    \eta_m = 4\,k_{z,0}\,k_{z,m}\,\left|\hat{g}_m\right|^2, \qquad
    k_{z} = \frac{2\pi}{\lambda}\cos\theta\,\sin\gamma

with :math:`\hat{g}_m` the Fourier coefficient of the physical groove height. No
Kirchhoff assumption, no transmittance function, no Appendix D.

Kirchhoff theory needs gentle slopes; perturbation theory needs shallow grooves.
A grating that is *both* is described correctly by both, and in that overlap the
two must agree to :math:`O(h)`. That is the only place in this project where an
independent theory pins the scalar solver's **absolute normalisation**, which is
what makes it worth the file.

**What it found.** On first run, against ``E_m = |G_m|^2`` with no flux factor,
every check here failed -- and failed only on the orders furthest from Littrow,
by exactly the ratio in :meth:`TestTheOracleItself
.test_it_disagrees_with_the_bare_fourier_coefficient`:

===============  =========  ======  ==============  =========================
Mount            lambda/p   m       scalar / exact  4c_a c_b / (c_a + c_b)^2
===============  =========  ======  ==============  =========================
in-plane         0.80       +1      1.0137          0.9865
in-plane         0.80       -1      1.6375          0.6107
off-plane        0.0127     +1      1.0003          0.9997
off-plane        0.0127     -1      1.2156          0.8227
===============  =========  ======  ==============  =========================

``geometry.flux_obliquity`` is what closes that gap; ``docs/findings.md`` carries
the narrative.

The Neumann (TM) first-order result is
:math:`4k^2|\hat{g}_m|^2(1 + \sin\alpha\sin\beta_m)^2 / (\cos\alpha\cos\beta_m)`
and diverges from the Dirichlet one away from Littrow. A polarization-blind
model cannot match both; TE is the conventional scalar correspondence, and the
size of the TE/TM spread is itself a measure of where scalar theory stops
meaning anything. See ``docs/theory/scalar.md`` section 3.
"""

import numpy as np
import pytest
from scipy.special import jv

from gratinglab.geometry import cos_beta, is_propagating, sin_beta
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Lamellar, Sinusoidal
from gratinglab.solvers import scalar

UNPOL = "unpolarized"

# Same two-mount discipline as `test_scalar.py`: an in-plane case cannot detect
# an error in the sin(gamma) handling, and off-plane is the primary application.
#
# **These mounts are chosen to disperse, not to be typical.** The quantity this
# file constrains departs from unity as cos(beta_m) departs from cos(alpha), so a
# near-Littrow mount would make the comparison pass no matter what the solver
# did. The first draft of this file used alpha=10 deg at lambda/p = 0.36 and
# alpha=25 deg at lambda/p = 0.01; both agreed to 0.2% and tested nothing.
# `test_the_mounts_actually_discriminate` is what stops that from recurring.
MOUNTS = [
    (1000.0, Illumination(alpha_deg=10.0, gamma_deg=90.0, polarization=UNPOL), 800.0),
    (315.15, Illumination(alpha_deg=20.0, gamma_deg=1.25, polarization=UNPOL), 4.0),
]
MOUNT_IDS = ["in-plane", "off-plane"]


def sinusoid_coefficient(depth_nm, order):
    r""":math:`|\hat{g}_m|` for :math:`g(t) = \frac{h}{2}(1 - \cos 2\pi t)`.

    Only :math:`m = 0, \pm 1` are non-zero -- a sinusoid has exactly one
    harmonic, which is why it is the cleanest profile to test first order on.
    """
    return depth_nm / 4.0 if abs(order) == 1 else 0.0


def lamellar_coefficient(depth_nm, duty, order):
    r""":math:`|\hat{g}_m| = h\,|\sin(\pi m w)| / (\pi |m|)` for a binary groove."""
    if order == 0:
        return depth_nm * duty
    return depth_nm * abs(np.sin(np.pi * order * duty)) / (np.pi * abs(order))


def first_order_efficiency(problem, illumination, wavelength, orders, coefficient):
    r"""Dirichlet first-order result, :math:`4 k_{z,0} k_{z,m} |\hat{g}_m|^2`.

    ``coefficient`` is a callable ``order -> |g_m|`` in nm, supplied analytically
    by the caller rather than transformed numerically here: a numerical
    transform of ``profile.height`` would share code with the thing under test.
    """
    sines = sin_beta(
        orders,
        wavelength,
        problem.period,
        illumination.sin_alpha,
        illumination.sin_gamma,
    )
    live = is_propagating(sines)
    cosines = cos_beta(sines[live])

    k = 2.0 * np.pi / wavelength
    k_z0 = k * illumination.cos_alpha * illumination.sin_gamma
    k_zm = k * cosines * illumination.sin_gamma
    moduli = np.array([coefficient(int(m)) for m in orders[live]])
    return live, 4.0 * k_z0 * k_zm * moduli**2


class TestTheOracleItself:
    """Before trusting it against the solver, show it is not vacuous."""

    def test_it_disagrees_with_the_bare_fourier_coefficient(self):
        r"""Two independent closed forms, and the ratio between them.

        The sinusoid's scalar answer is :math:`J_m^2(\varphi_m/2)`, which in the
        shallow limit tends to :math:`(\varphi_m/4)^2`. Perturbation theory says
        :math:`4k_{z,0}k_{z,m}|\hat{g}_m|^2`. With
        :math:`\varphi_m = k h \sin\gamma(\cos\alpha + \cos\beta_m)` and
        :math:`|\hat{g}_{\pm 1}| = h/4` the ratio is exactly

        .. math::
            \frac{(\cos\alpha + \cos\beta_m)^2}{4\cos\alpha\cos\beta_m}

        which is 1 only when :math:`\cos\alpha = \cos\beta_m`. This is pure
        arithmetic between two closed forms -- no solver in it -- so it is the
        statement of *why* ``geometry.flux_obliquity`` has to exist, and it
        cannot rot as the solver changes.
        """
        period, depth_fraction, alpha_deg, wavelength = 1000.0, 1e-5, 10.0, 800.0
        problem = Problem(
            period=period, profile=Sinusoidal(depth_fraction=depth_fraction)
        )
        ill = Illumination.classical(alpha=alpha_deg, polarization=UNPOL)
        depth_nm = depth_fraction * period

        order = -1
        sine = float(sin_beta(order, wavelength, period, ill.sin_alpha, ill.sin_gamma))
        cosine = float(cos_beta(sine))
        phi = (
            (2.0 * np.pi / wavelength)
            * depth_nm
            * ill.sin_gamma
            * (ill.cos_alpha + cosine)
        )

        bessel = jv(order, phi / 2.0) ** 2
        k = 2.0 * np.pi / wavelength
        perturbation = (
            4.0
            * (k * ill.cos_alpha * ill.sin_gamma)
            * (k * cosine * ill.sin_gamma)
            * sinusoid_coefficient(depth_nm, order) ** 2
        )
        obliquity = (ill.cos_alpha + cosine) ** 2 / (4.0 * ill.cos_alpha * cosine)

        assert bessel / perturbation == pytest.approx(obliquity, rel=1e-6)
        # And at this geometry the two theories are 64% apart, so an oracle that
        # agreed with the solver here would be measuring nothing.
        assert obliquity == pytest.approx(1.637, rel=1e-3)
        del problem  # geometry only; the solver is deliberately absent here

    @pytest.mark.parametrize("mount", MOUNTS, ids=MOUNT_IDS)
    def test_the_mounts_actually_discriminate(self, mount):
        r"""The oracle must have somewhere to disagree, or it proves nothing.

        Every check in this file compares two theories whose ratio is
        :math:`(\cos\alpha + \cos\beta_m)^2 / (4\cos\alpha\cos\beta_m)`. At
        Littrow that ratio is exactly 1 and the comparison is empty. This pins
        that each mount reaches an order where the two theories are at least
        15% apart, so an agreement below is a result rather than a tautology.
        """
        period, ill, wavelength = mount
        orders = np.arange(-8, 9)
        sines = sin_beta(orders, wavelength, period, ill.sin_alpha, ill.sin_gamma)
        live = is_propagating(sines)
        cosines = cos_beta(sines[live])
        obliquity = (
            4.0 * ill.cos_alpha * cosines / (ill.cos_alpha + cosines) ** 2
        )
        assert obliquity.min() < 0.85, (
            f"mount is too near Littrow to test anything: obliquity spans "
            f"{obliquity.min():.4f}-{obliquity.max():.4f}"
        )

    @pytest.mark.parametrize(
        "mount,mount_id", list(zip(MOUNTS, MOUNT_IDS)), ids=MOUNT_IDS
    )
    def test_it_is_an_asymptotic_statement_and_converges_as_the_groove_shallows(
        self, mount, mount_id
    ):
        """The ratio must settle as depth falls, or the comparison is meaningless.

        Guards the choice of ``depth_fraction`` below: if 1e-5 were not already
        deep into the asymptotic regime, an agreement there would be luck.
        """
        period, ill, wavelength = mount
        ratios = []
        for depth_fraction in (1e-3, 1e-4, 1e-5):
            problem = Problem(
                period=period, profile=Sinusoidal(depth_fraction=depth_fraction)
            )
            scan = scalar.solve(problem, ill, [wavelength], quadrature_points=4096)
            live, expected = first_order_efficiency(
                problem,
                ill,
                wavelength,
                scan.orders,
                lambda m: sinusoid_coefficient(depth_fraction * period, m),
            )
            keep = np.abs(scan.orders[live]) == 1
            ratios.append(float(np.mean(scan.efficiency[0][live][keep] / expected[keep])))

        # Successive ratios converge: each step is closer to the last than the
        # step before it, i.e. the sequence is settling on a limit.
        assert abs(ratios[2] - ratios[1]) < abs(ratios[1] - ratios[0]), (
            f"{mount_id}: ratios {ratios} are not converging"
        )


class TestAgainstFirstOrderPerturbationTheory:
    """The check that constrains the absolute normalisation."""

    @pytest.mark.parametrize("mount", MOUNTS, ids=MOUNT_IDS)
    def test_a_shallow_sinusoid_matches_the_exact_first_order_result(self, mount):
        period, ill, wavelength = mount
        depth_fraction = 1e-5
        problem = Problem(
            period=period, profile=Sinusoidal(depth_fraction=depth_fraction)
        )
        scan = scalar.solve(problem, ill, [wavelength], quadrature_points=4096)

        live, expected = first_order_efficiency(
            problem,
            ill,
            wavelength,
            scan.orders,
            lambda m: sinusoid_coefficient(depth_fraction * period, m),
        )
        keep = np.abs(scan.orders[live]) == 1
        assert keep.any(), "mount admits no first order to test"
        # atol=0 is load-bearing. A shallow grating's first-order efficiency is
        # ~1e-9, which is three orders of magnitude below numpy's default
        # atol=1e-8, so the default would pass this comparison against anything
        # at all -- including a solver returning zero.
        assert np.allclose(
            scan.efficiency[0][live][keep], expected[keep], rtol=2e-3, atol=0.0
        ), (
            f"scalar {scan.efficiency[0][live][keep]} vs "
            f"perturbation {expected[keep]}"
        )

    @pytest.mark.parametrize("duty", [0.3, 0.5])
    @pytest.mark.parametrize("mount", MOUNTS, ids=MOUNT_IDS)
    def test_a_shallow_lamellar_matches_the_exact_first_order_result(
        self, mount, duty
    ):
        """A second profile, and one whose spectrum is not a single harmonic.

        The sinusoid can only test :math:`|m| = 1`. A binary groove has
        :math:`\\hat{g}_m \\neq 0` at every order, so this reaches the higher
        ones -- where :math:`\\cos\\beta_m` departs furthest from
        :math:`\\cos\\alpha` and any normalisation error is largest.
        """
        period, ill, wavelength = mount
        depth_fraction = 1e-5
        problem = Problem(
            period=period,
            profile=Lamellar(depth_fraction=depth_fraction, duty_cycle=duty),
        )
        scan = scalar.solve(problem, ill, [wavelength], quadrature_points=8192)

        live, expected = first_order_efficiency(
            problem,
            ill,
            wavelength,
            scan.orders,
            lambda m: lamellar_coefficient(depth_fraction * period, duty, m),
        )
        # Skip m = 0 -- zeroth order is O(1), not O(h), so first order says
        # nothing about it -- and orders the duty cycle extinguishes. The
        # extinction test has to be relative: at w = 0.5, sin(pi m w) for even m
        # is 1.2e-16 rather than 0, so `expected > 0` leaves an order in that
        # compares one numerical zero against another and can fail by any
        # factor at all.
        keep = (scan.orders[live] != 0) & (expected > 1e-9 * expected.max())
        assert keep.any()
        assert np.allclose(  # atol=0: see the sinusoid case above
            scan.efficiency[0][live][keep], expected[keep], rtol=5e-3, atol=0.0
        ), (
            f"orders {scan.orders[live][keep]}: scalar "
            f"{scan.efficiency[0][live][keep]} vs perturbation {expected[keep]}"
        )
