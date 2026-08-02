"""The comparison harness -- the end-to-end spine.

If these pass, one Problem really does reach several methods unchanged and
their answers really do line up for comparison.
"""

import numpy as np
import pytest

from gratinglab.compare import align, records, sweep
from gratinglab.illumination import Illumination
from gratinglab.io.efficiency_table import read_scan
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.result import EfficiencyScan, Provenance

from .conftest import reference_dir

UNPOL = "unpolarized"
PROBLEM = Problem(period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5))
ILL = Illumination.offplane(graze=1.25, azimuth=19.99, polarization=UNPOL)


class TestSweep:
    def test_runs_a_single_solver(self):
        scans = sweep(PROBLEM, ILL, np.linspace(1.0, 5.0, 9), ["scalar"])
        assert len(scans) == 1
        assert scans[0].provenance.method == "scalar"

    def test_passes_per_method_options(self):
        scans = sweep(
            PROBLEM,
            ILL,
            [2.4],
            ["scalar"],
            options={"scalar": {"quadrature_points": 4096}},
        )
        assert scans[0].provenance.truncation == 4096

    def test_accepts_a_precomputed_scan_as_a_method(self):
        """An imported table participates exactly like a live solver."""
        precomputed = EfficiencyScan(
            wavelengths=np.array([1.0, 2.0, 3.0]),
            orders=np.array([0, 1]),
            efficiency=np.array([[0.6, 0.4], [0.5, 0.5], [0.7, 0.3]]),
            propagating=np.ones((3, 2), dtype=bool),
            provenance=Provenance("reference"),
        )
        scans = sweep(PROBLEM, ILL, [1.0, 2.0, 3.0], ["scalar", precomputed])
        assert [s.provenance.method for s in scans] == ["scalar", "reference"]

    def test_every_method_lands_on_the_requested_grid(self):
        grid = np.linspace(1.0, 5.0, 7)
        for scan in sweep(PROBLEM, ILL, grid, ["scalar"]):
            assert np.allclose(scan.wavelengths, grid)

    def test_unknown_method_raises(self):
        with pytest.raises(KeyError, match="no solver registered"):
            sweep(PROBLEM, ILL, [2.4], ["nonexistent"])


class TestResampling:
    def test_uses_nearest_neighbour_not_interpolation(self):
        """Interpolating across a passing-off edge would invent a value no
        method produced."""
        source = EfficiencyScan(
            wavelengths=np.array([1.0, 2.0]),
            orders=np.array([0]),
            efficiency=np.array([[0.2], [0.8]]),
            propagating=np.ones((2, 1), dtype=bool),
            provenance=Provenance("reference"),
        )
        resampled = sweep(PROBLEM, ILL, [1.4, 1.6], [source])[0]
        assert resampled.efficiency[0, 0] == 0.2  # nearest 1.0, not 0.44
        assert resampled.efficiency[1, 0] == 0.8  # nearest 2.0, not 0.56

    def test_marks_out_of_range_points_as_non_propagating(self):
        """Better a gap than an extrapolation."""
        source = EfficiencyScan(
            wavelengths=np.array([2.0, 3.0]),
            orders=np.array([0]),
            efficiency=np.array([[0.5], [0.5]]),
            propagating=np.ones((2, 1), dtype=bool),
            provenance=Provenance("reference"),
        )
        resampled = sweep(PROBLEM, ILL, [1.0, 2.5, 9.0], [source])[0]
        assert resampled.propagating[:, 0].tolist() == [False, True, False]
        assert resampled.efficiency[0, 0] == 0.0


class TestAlign:
    def test_unions_the_order_sets(self):
        """A method resolving fewer orders must not truncate the comparison."""
        wide = EfficiencyScan(
            wavelengths=np.array([1.0]),
            orders=np.array([-2, -1, 0, 1, 2]),
            efficiency=np.array([[0.1, 0.2, 0.4, 0.2, 0.1]]),
            propagating=np.ones((1, 5), dtype=bool),
            provenance=Provenance("wide"),
        )
        narrow = EfficiencyScan(
            wavelengths=np.array([1.0]),
            orders=np.array([0, 1]),
            efficiency=np.array([[0.6, 0.4]]),
            propagating=np.ones((1, 2), dtype=bool),
            provenance=Provenance("narrow"),
        )
        comparison = align([wide, narrow])
        assert comparison.orders.tolist() == [-2, -1, 0, 1, 2]
        assert comparison.efficiency[1, 0].tolist() == [0.0, 0.0, 0.6, 0.4, 0.0]

    def test_rejects_mismatched_wavelength_grids(self):
        def scan(wavelengths, name):
            n = len(wavelengths)
            return EfficiencyScan(
                wavelengths=np.asarray(wavelengths, dtype=float),
                orders=np.array([0]),
                efficiency=np.full((n, 1), 0.5),
                propagating=np.ones((n, 1), dtype=bool),
                provenance=Provenance(name),
            )

        with pytest.raises(ValueError, match="different wavelength grids"):
            align([scan([1.0, 2.0], "a"), scan([1.0, 3.0], "b")])

    def test_disambiguates_repeated_method_names(self):
        comparison = align(
            sweep(
                PROBLEM,
                ILL,
                [2.4],
                ["scalar", "scalar"],
                options={"scalar": {"quadrature_points": 2048}},
            )
        )
        assert comparison.methods == ("scalar", "scalar#2")

    def test_rejects_an_empty_comparison(self):
        with pytest.raises(ValueError, match="nothing to compare"):
            align([])


class TestComparisonQueries:
    def setup_method(self):
        self.comparison = align(
            sweep(PROBLEM, ILL, np.linspace(1.0, 5.0, 21), ["scalar"])
        )

    def test_order_returns_a_curve_per_method(self):
        curves = self.comparison.order(1)
        assert set(curves) == {"scalar"}
        assert len(curves["scalar"]) == 21

    def test_identical_methods_have_zero_difference(self):
        pair = align(sweep(PROBLEM, ILL, [2.4], ["scalar", "scalar"]))
        assert pair.max_abs_difference("scalar", "scalar#2") == 0.0

    def test_summary_locates_the_worst_disagreement(self):
        pair = align(
            sweep(
                PROBLEM,
                ILL,
                np.linspace(1.0, 5.0, 11),
                ["scalar", "scalar"],
                options={"scalar": {"quadrature_points": 2048}},
            )
        )
        summary = pair.summary("scalar", "scalar#2")
        assert set(summary) == {
            "max_abs_difference",
            "rms_difference",
            "at_wavelength",
            "at_order",
        }


class TestRecords:
    def test_produces_tidy_rows_across_methods(self):
        scans = sweep(PROBLEM, ILL, np.linspace(1.0, 3.0, 5), ["scalar"])
        rows = records(scans)
        assert len(rows) == 5 * len(scans[0].orders)
        assert {"method", "wavelength", "order", "efficiency"} <= set(rows[0])


REF = reference_dir()


@pytest.mark.skipif(REF is None, reason="reference corpus unavailable")
class TestAgainstImportedReferenceData:
    """The end-to-end evidence that the spine works.

    Scalar and the integral method run on one Problem and land on one grid.
    Their *agreement* is not asserted -- quantifying the disagreement is the
    point of the project, and both are relative here (perfect conductivity, no
    optical constants), so the comparison is apples to apples.
    """

    def test_scalar_and_imported_reference_line_up(self, ref_dir):
        path = ref_dir / "OGRE" / "tastetest_perf_wavescan.txt"
        if not path.exists():
            pytest.skip("TASTE wavescan not in corpus")

        reference = read_scan(path, method="integral")
        wavelengths = reference.wavelengths[::20]
        comparison = align(
            sweep(
                PROBLEM,
                ILL,
                wavelengths,
                ["scalar", reference],
                options={"scalar": {"quadrature_points": 8192}},
            )
        )

        assert comparison.methods == ("scalar", "integral")
        assert comparison.efficiency.shape[0] == 2
        assert np.isfinite(comparison.efficiency).all()
        assert (comparison.efficiency >= 0).all()

        summary = comparison.summary("scalar", "integral")
        assert np.isfinite(summary["max_abs_difference"])
        assert summary["at_wavelength"] in wavelengths

    def test_recovered_geometry_reproduces_the_observed_order_range(self, ref_dir):
        """The strongest check on benchmarks/corpus.toml.

        The geometry was recovered from which orders propagate, so running our
        own grating equation with it must regenerate exactly that order set.
        """
        path = ref_dir / "OGRE" / "tastetest_perf_wavescan.txt"
        if not path.exists():
            pytest.skip("TASTE wavescan not in corpus")

        reference = read_scan(path)
        ours = sweep(PROBLEM, ILL, reference.wavelengths, ["scalar"])[0]

        for row in range(0, len(reference), 25):
            theirs = set(reference.orders[reference.propagating[row]].tolist())
            mine = set(ours.orders[ours.propagating[row]].tolist())
            assert mine == theirs, (
                f"at lambda={reference.wavelengths[row]:.2f} nm we propagate "
                f"{sorted(mine)} but PCGrate propagates {sorted(theirs)}"
            )
