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
