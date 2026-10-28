"""RCWA against published lamellar benchmarks, transcribed from the papers.

Every reference value here was read in a primary text that is freely
available, and the source is named beside it. The canonical papers this was
meant to draw on (Moharam, Grann, Pommet & Gaylord, JOSA A 12, 1068 (1995);
Lalanne & Morris, JOSA A 13, 779 (1996); Li, JOSA A 13, 1870 (1996); Granet &
Guizal, JOSA A 13, 1019 (1996)) are closed access, and none of their numbers
could be confirmed, so none is used. Granet's 1999 values enter only as
quoted by two later papers that print them. What could not be matched is in
``docs/theory/rcwa.md`` section 9, not here.

Each tolerance is set from a truncation sweep, measured 2026-10-27 with
OPENBLAS_NUM_THREADS=1 and recorded in the class docstring: the error at the
truncation used, with margin, and the sweep's direction where it is
monotone, so a regression in convergence fails as well as one in the value.
"""

import numpy as np
import pytest

from gratinglab.convergence import check_convergence
from gratinglab.illumination import Illumination
from gratinglab.materials.optical import OpticalConstants, write_ari
from gratinglab.problem import Problem
from gratinglab.profiles import Lamellar
from gratinglab.solvers import get_solver
from gratinglab.solvers.rcwa._core import incident_field, solve_stack
from gratinglab.solvers.rcwa._stack import Slab, Stack

rcwa = get_solver("rcwa")


def in_plane(stack, *, wavelength, period, alpha, polarization, truncation):
    a = np.radians(alpha)
    kx0, kz0 = -np.sin(a), np.cos(a)
    return solve_stack(
        stack,
        wavelength=wavelength,
        period=period,
        kx0=kx0,
        ky=0.0,
        incidents=[incident_field(polarization, kx0, 0.0, kz0)],
        truncation=truncation,
    )


class TestGranetMetallicTM:
    """Metal ridges on the same metal, air grooves, period = wavelength =
    depth = 1 um, fill 0.5, 30 degrees, n = 0.22 + 6.71i (aluminium near
    10.6 um, scaled). Reflected zeroth order in TM:

    * 84.848%: Lalanne & Hugonin, JOSA A 17, 1033 (2000), Fig. 3 caption,
      p. 1036 ("The exact value of the reflected zeroth order is 84.848% for
      TM"), HAL hal-00867708;
    * 0.8484781: Demesy et al., ch. 5 of Popov (ed.), *Gratings: Theory and
      Numeric Applications*, 2nd ed. (2014), Tab. 5.1, p. 5.14, row "Result
      given by [17]", HAL hal-01084458. The same chapter gives the substrate
      as eps = -44.9757 + 2.9524i, which is (0.22 + 6.71i)^2.

    Both quote G. Granet, JOSA A 16, 2510 (1999). With fill 0.5 it does not
    matter which material the fill names. This is exactly a ``Problem``:
    ridges of ``coating`` in vacuum on a substrate of it.

    Error against 0.8484781 by truncation, Li's rule:

    ====  =========
    M     error
    ====  =========
    10    -6.4e-3
    20    -4.2e-3
    40    -1.7e-3
    80    -6.3e-4
    120   -3.5e-4
    160   -2.5e-4
    200   -2.1e-4
    ====  =========

    from below and monotone, slowly: the metal's corners make the field
    singular (``rcwa.md`` section 3). Lalanne & Hugonin's own RCWA column
    reads 0.84677 at 81 orders, error -1.7e-3, so this backend at M = 40
    (81 orders) matches theirs. Asserted at M = 80 within 1e-3.
    """

    REFERENCE = 0.8484781

    def test_reflected_zeroth_order(self, granet_report):
        scan = granet_report.scan
        assert scan.provenance.truncation == 80
        specular = float(scan.efficiency[0, scan.orders == 0][0])
        assert specular == pytest.approx(self.REFERENCE, abs=1e-3)

    def test_it_converges_towards_it_from_below(self, granet_report):
        assert list(granet_report.values) == [20, 40, 80]
        changes = granet_report.differences
        assert changes[1] < changes[0]
        scan = granet_report.scan
        assert float(scan.efficiency[0, scan.orders == 0][0]) < self.REFERENCE

    def test_one_call_certifies_it(self, granet_report):
        """``convergence.check_convergence`` (and ``converged_scan``, which
        returns its scan) is the opt-in "converge to tolerance" path: the
        stamped scan says it converged, and the published value agrees to
        the tolerance it was certified at."""
        assert granet_report.passed
        provenance = granet_report.scan.provenance
        assert provenance.converged is True
        assert provenance.notes["convergence"]["tolerance"] == 3e-3
        specular = float(granet_report.scan.efficiency[0, granet_report.scan.orders == 0][0])
        assert abs(specular - self.REFERENCE) < 3e-3


@pytest.fixture(scope="module")
def granet_report(tmp_path_factory):
    """The convergence study of :class:`TestGranetMetallicTM`, run once."""
    directory = tmp_path_factory.mktemp("granet")
    index = 0.22 + 6.71j
    write_ari(
        OpticalConstants(
            name="granet-metal",
            wavelength_nm=(500.0, 2000.0),
            decrement=(1.0 - index.real,) * 2,
            absorption=(index.imag,) * 2,
            source="Granet 1999 via Lalanne & Hugonin 2000",
        ),
        directory / "granet-metal.ari",
    )
    patch = pytest.MonkeyPatch()
    patch.setenv("GRATINGLAB_MATERIALS", str(directory))
    try:
        return check_convergence(
            rcwa,
            Problem(
                period=1000.0,
                profile=Lamellar(depth_fraction=1.0, duty_cycle=0.5),
                coating="granet-metal",
            ),
            Illumination.classical(alpha=30.0, polarization="TM"),
            [1000.0],
            values=(20, 40, 80),
            tolerance=3e-3,
        )
    finally:
        patch.undo()


class TestLalanneHugoninDielectricTM:
    """Ridges of n = 2.3 in air on a substrate of n = 1.5, period 2 um,
    wavelength 1 um, depth 1 um, fill 0.234, 30 degrees, TM: the transmitted
    minus-first order is 51.062%. Lalanne & Hugonin, JOSA A 17, 1033 (2000),
    Fig. 6 caption, p. 1038 ("The exact value for the diffraction efficiency
    is 51.062%").

    The paper's Fig. 1, which says what the fill fraction measures and which
    way the orders are numbered, is an image that could not be read. Of the
    four readings, exactly one reaches the value: the fill is the 2.3
    ridge's, and the order is this code's ``-1``. The three others give
    0.0023, 0.020 and 0.034 (asserted below), so the match identifies the
    convention rather than being tuned to it.

    Error against 0.51062 by truncation:

    ====  =========
    M     error
    ====  =========
    10    +1.0e-2
    20    +1.2e-3
    30    +4.1e-4
    40    +2.0e-4
    60    +5.4e-5
    80    +1.4e-5
    120   -1.1e-5
    ====  =========

    The reference is printed to five digits, so 120 sits within its rounding.
    Asserted at M = 60 within 1e-4, and monotone over 20, 30, 40.
    """

    REFERENCE = 0.51062

    @staticmethod
    def transmitted(truncation, ridge_fill=0.234):
        solution = in_plane(
            Stack(1.0, (Slab(1000.0, 1.0, ((0.0, ridge_fill, 2.3),)),), 1.5),
            wavelength=1000.0,
            period=2000.0,
            alpha=30.0,
            polarization="TM",
            truncation=truncation,
        )
        return solution.transmitted[0], truncation

    def test_transmitted_minus_first(self):
        orders, m = self.transmitted(60)
        assert orders[m - 1] == pytest.approx(self.REFERENCE, abs=1e-4)

    def test_it_converges_towards_it(self):
        errors = []
        for truncation in (20, 30, 40):
            orders, m = self.transmitted(truncation)
            errors.append(orders[m - 1] - self.REFERENCE)
        assert errors[0] > errors[1] > errors[2] > 0.0

    def test_no_other_reading_matches(self):
        orders, m = self.transmitted(40)
        other, _ = self.transmitted(40, ridge_fill=1.0 - 0.234)
        for value in (orders[m + 1], other[m - 1], other[m + 1]):
            assert abs(value - self.REFERENCE) > 0.4


class TestGralakSiliconRodsTM:
    """Free-standing rods, eps = 12.96, width 0.28 d, height d / (2 sqrt 2),
    wavelength = period d, 45 degrees, vacuum both sides: reflected zeroth
    order "0.7323 and 0.9487 in s and p polarizations". B. Gralak, ch. 10 of
    Popov (ed.), *Gratings: Theory and Numeric Applications*, 2nd ed.
    (2014), section 10.6, pp. 10.12-10.13, HAL hal-01084458; s is E along
    the grooves and p H along them, as the chapter states, so p is TM here.

    Only p is asserted. s converges cleanly to 0.74230 here (M = 40 to 120
    agree to 6e-5), 0.0100 above the printed 0.7323, and no other order
    is near it; see ``rcwa.md`` section 9.

    Error against 0.9487 in p by truncation:

    ====  =========
    M     error
    ====  =========
    10    -1.7e-2
    20    -2.2e-3
    40    -3.2e-4
    80    -6.6e-5
    100   -4.7e-5
    120   -4.0e-5
    ====  =========

    The reference has four digits, so anything within 5e-5 is agreement.
    Asserted at M = 80 within 1e-4, and monotone over 20, 40, 80.
    """

    REFERENCE = 0.9487

    @staticmethod
    def specular(truncation):
        solution = in_plane(
            Stack(
                1.0,
                (Slab(1000.0 / (2.0 * np.sqrt(2.0)), 1.0, ((0.0, 0.28, np.sqrt(12.96)),)),),
                1.0,
            ),
            wavelength=1000.0,
            period=1000.0,
            alpha=45.0,
            polarization="TM",
            truncation=truncation,
        )
        return float(solution.reflected[0, truncation])

    def test_reflected_zeroth_order_in_p(self):
        assert self.specular(80) == pytest.approx(self.REFERENCE, abs=1e-4)

    def test_it_converges_towards_it_from_below(self):
        errors = [self.specular(m) - self.REFERENCE for m in (20, 40, 80)]
        assert errors[0] < errors[1] < errors[2] < 0.0


def goray_schmidt(eps, depth, theta, phi, delta, psi, truncation, *, flip_p=False):
    """A conical lamellar case of Goray & Schmidt, in their own conventions.

    Their frame (Popov (ed.) 2014, ch. 12, section 12.2.1 and Fig. 12.2) is
    this project's: ``x`` across the grooves, ``y`` the normal pointing into
    the incident medium, ``z`` along the grooves, ``exp(-i omega t)``, and
    the incident wave vector ``k (sin t cos f, -cos t cos f, sin f)``. Orders
    are ``alpha_n = alpha + 2 pi n / d`` (section 12.2.2), the core's ``m``.
    The polarization basis (section 12.8.1) is ``s = k x y / |k x y|``,
    ``p = s x k / |k|``, with ``delta = arctan |E_s / E_p|`` and
    ``psi = -arg(E_s / E_p)``. In the core's frame ``(x, z_gs, -y_gs)``
    (``rcwa.md`` section 2) every vector is the same rotation of theirs, so
    the cross products are taken there directly.

    The grating is one profile between the cover (``eps = 1``) and the
    substrate ``eps``, so the ridges are the substrate's material. With fill
    0.5 a ridge and a groove differ by half a period, which moves only the
    orders' phases, so which one ``c`` measures cannot change an efficiency.
    ``flip_p`` reverses ``p``, the other reading a sign slip would give.
    """
    t, f = np.radians(theta), np.radians(phi)
    k = np.array([np.sin(t) * np.cos(f), np.sin(f), np.cos(t) * np.cos(f)])
    s = np.cross(k, [0.0, 0.0, -1.0])  # y_gs is -z in the core's frame
    s /= np.linalg.norm(s)
    p = np.cross(s, k) * (-1.0 if flip_p else 1.0)
    field = p + np.tan(np.radians(delta)) * np.exp(-1j * np.radians(psi)) * s
    n = np.sqrt(complex(eps))
    return solve_stack(
        Stack(1.0, (Slab(depth, 1.0, ((0.0, 0.5, n),)),), n),
        wavelength=500.0,
        period=1000.0,
        kx0=k[0],
        ky=k[1],
        incidents=[field / np.linalg.norm(field)],
        truncation=truncation,
    )


def order(solution, side, m):
    rows = solution.reflected if side == "R" else solution.transmitted
    return float(rows[0, solution.harmonics == m][0])


class TestGoraySchmidtDielectricConical:
    """Dielectric lamellar, conical, elliptically polarized: ``c/d = 0.5``,
    ``2H/d = 0.5``, cover ``eps = 1``, substrate and ridges ``eps = 2.25``,
    ``lambda/d = 0.5``, ``theta = 35.264``, ``phi = 30``, ``delta = 45``,
    ``psi = 90`` (degrees). L. Goray and G. Schmidt, ch. 12 of Popov (ed.),
    *Gratings: Theory and Numeric Applications*, 2nd ed. (2014), Table 12.2
    and section 12.9.1, p. 12.61, HAL hal-01084458. Their integral method at
    400 points, energy balance ~1e-5. Efficiencies in percent, every
    propagating order, reflected and transmitted. The conventions are
    mapped in :func:`goray_schmidt`.

    Error by truncation, in units of half the last printed digit (1.00 is
    exactly the rounding of the printed value):

    ====  =====  =====  =====  =====  =====  =====  =====  =====
    M     R-2    R-1    R0     T-3    T-2    T-1    T0     T1
    ====  =====  =====  =====  =====  =====  =====  =====  =====
    20    +6.56  -2.93  +0.21  +0.53  +0.41  +0.52  -0.00  -0.82
    40    +1.61  -0.41  -0.26  +0.05  +0.36  +0.27  +0.26  -0.66
    80    +0.82  +0.12  -0.39  -0.16  +0.35  +0.23  +0.33  -0.64
    120   +0.71  +0.21  -0.41  -0.20  +0.34  +0.23  +0.34  -0.64
    ====  =====  =====  =====  =====  =====  =====  =====  =====

    So from M = 80 all eight agree with the table to its printed precision,
    and the remaining residuals are the table's rounding, not truncation.
    Asserted at M = 80 within half a unit of each printed digit. The
    opposite handedness (``p`` reversed, equivalently the conjugate ratio)
    misses every order by more than 100 such units, so the match fixes the
    polarization convention rather than tolerating either. The structure
    is lossless, so ``R + T = 1`` holds to rounding (a check of the
    bookkeeping, not of convergence).
    """

    #: (side, order): (efficiency in percent as printed, decimals printed)
    TABLE = {
        ("R", -2): (0.1612, 4),
        ("R", -1): (0.3807, 4),
        ("R", 0): (1.854, 3),
        ("T", -3): (3.363, 3),
        ("T", -2): (10.35, 2),
        ("T", -1): (31.87, 2),
        ("T", 0): (14.19, 2),
        ("T", 1): (37.83, 2),
    }

    @staticmethod
    def solve(truncation, **kwargs):
        return goray_schmidt(2.25, 500.0, 35.264, 30.0, 45.0, 90.0, truncation, **kwargs)

    @staticmethod
    def half_digit(decimals):
        return 0.5 * 10.0**-decimals / 100.0

    def test_every_order_to_the_printed_digit(self):
        solution = self.solve(80)
        for (side, m), (percent, decimals) in self.TABLE.items():
            value = order(solution, side, m)
            assert abs(value - percent / 100.0) < self.half_digit(decimals), (side, m)

    def test_the_orders_are_all_the_power(self):
        solution = self.solve(80)
        total = solution.reflected[0].sum() + solution.transmitted[0].sum()
        assert total == pytest.approx(1.0, abs=1e-12)
        listed = sum(order(solution, side, m) for side, m in self.TABLE)
        assert listed == pytest.approx(1.0, abs=1e-12)

    def test_the_slowest_order_converges_towards_it(self):
        errors = [
            order(self.solve(m), "R", -2) - 0.1612 / 100.0 for m in (20, 40, 80)
        ]
        assert errors[0] > errors[1] > errors[2] > 0.0

    def test_the_other_handedness_does_not_match(self):
        solution = self.solve(80, flip_p=True)
        for (side, m), (percent, decimals) in self.TABLE.items():
            miss = abs(order(solution, side, m) - percent / 100.0)
            assert miss > 100.0 * self.half_digit(decimals), (side, m)


class TestGoraySchmidtMetallicConical:
    """Metallic lamellar, conical, linearly polarized at 45 degrees between
    s and p: ``c/d = 0.5``, ``2H/d = 1``, cover ``eps = 1``, substrate and
    ridges ``eps = -24.99 + 1i`` (printed as the pair ``(-24.99, 1)``;
    ``exp(-i omega t)``, so ``Im eps > 0`` is loss), ``lambda/d = 0.5``,
    ``theta = 22.208``, ``phi = 20.705``, ``delta = 45``, ``psi = 0``.
    Goray & Schmidt, ch. 12 of Popov (ed.) 2014, Table 12.1 and section
    12.9.1, p. 12.61, HAL hal-01084458: 400 points, energy balance ~1e-6.
    Reflected orders only, in percent: 7.52, 13.25, 44.27, 31.05.

    This backend converges, by truncation, to (percent):

    ====  ========  ========  ========  ========  ======
    M     R-2       R-1       R0        R1        A
    ====  ========  ========  ========  ========  ======
    20    7.32686   13.20862  44.77115  30.76040  3.9330
    40    7.51128   13.25739  44.22653  31.03343  3.9714
    80    7.54694   13.26440  44.15887  31.08989  3.9399
    160   7.55388   13.26487  44.15589  31.10358  3.9218
    240   7.55510   13.26477  44.15716  31.10652  3.9164
    ====  ========  ========  ========  ========  ======

    which is not the table: ``R0`` settles 1.1e-3 below it and ``R1`` and
    ``R-2`` 5.7e-4 and 3.5e-4 above, ten to a hundred times the printed
    rounding. The integral solver (``_finite.solve_finite_states``, a
    different method, fed the same state as ``(E_z, B_z)``, corner grading
    4) agrees with RCWA rather than with the table:

    ======  ========  ========  ========  ========  ======
    points  R-2       R-1       R0        R1        A
    ======  ========  ========  ========  ========  ======
    400     7.53493   13.25233  44.24931  31.05256  3.9093
    600     7.55026   13.26085  44.18537  31.09340  3.9096
    800     7.55387   13.26297  44.17024  31.10321  3.9094
    1200    7.55559   13.26404  44.16283  31.10807  3.9094
    ======  ========  ========  ========  ========  ======

    At 1200 points it is within 6e-5 of RCWA at M = 240 in every order,
    and its ladder moves away from the table towards that answer. At 400
    points, the paper's own count, it sits within 2.2e-4 of the table. So
    the printed values look like an integral solution not yet converged in
    its discretization, which an energy balance of 1e-6 cannot reveal (it
    is a check of the solve, not of convergence; ``rcwa.md`` section 4).
    That is an inference from two of this project's solvers, not something
    the chapter states; see ``rcwa.md`` section 9.

    Asserted at M = 80 within 1.5e-3 of every printed order, which covers
    the measured disagreement and no more, so the case still pins the
    conical kinematics, the metal and the polarization mapping: the other
    linear polarization (``p`` reversed, ``E_s = -E_p``) misses ``R-1`` by
    0.093 and ``R-2`` by 0.051, and that is asserted too.
    """

    TABLE = {-2: 7.52, -1: 13.25, 0: 44.27, 1: 31.05}

    @staticmethod
    def solve(truncation, **kwargs):
        return goray_schmidt(
            -24.99 + 1j, 1000.0, 22.208, 20.705, 45.0, 0.0, truncation, **kwargs
        )

    def test_every_reflected_order_within_the_measured_gap(self):
        solution = self.solve(80)
        for m, percent in self.TABLE.items():
            assert order(solution, "R", m) == pytest.approx(percent / 100.0, abs=1.5e-3), m

    def test_absorption_closes_the_balance(self):
        """The table's reflected orders sum to 96.09%, so 3.91% absorbed;
        this backend's ``A`` at M = 80 is 3.940%."""
        solution = self.solve(80)
        absorbed = float(solution.absorbed[0])
        reflected = float(solution.reflected[0].sum())
        assert reflected + absorbed == pytest.approx(1.0, abs=1e-10)
        assert absorbed == pytest.approx(1.0 - sum(self.TABLE.values()) / 100.0, abs=5e-4)

    def test_the_other_linear_polarization_does_not_match(self):
        solution = self.solve(80, flip_p=True)
        assert abs(order(solution, "R", -1) - self.TABLE[-1] / 100.0) > 0.05
        assert abs(order(solution, "R", -2) - self.TABLE[-2] / 100.0) > 0.03
