"""Efficiency-table importer.

Split in two: unit tests against a small committed fixture, and corpus tests
that parse the real exported tables. The corpus tests are verification item 1
from the project plan -- column headers carry the ground truth, so the parse is
checkable without running any solver.
"""

import numpy as np
import pytest

from gratinglab.io.efficiency_table import read_table, read_scan, to_scan
from gratinglab.result import EfficiencyScan

from .conftest import reference_dir


class TestSyntheticFixture:
    def test_header_metadata(self, synthetic_wavescan):
        table = read_table(synthetic_wavescan)
        assert table.version == "PCGrate-SX 6.7.1"
        assert table.problem_name == "synthetic_wavescan"
        assert table.solved_at == "1/1/2026 00:00:00"
        assert table.calculating_time == "00:00:01"

    def test_orders_come_from_headers(self, synthetic_wavescan):
        table = read_table(synthetic_wavescan)
        assert table.orders.tolist() == [-1, 0, 1]
        assert table.polarization == "TE"
        assert table.reflected

    def test_wavelength_axis(self, synthetic_wavescan):
        table = read_table(synthetic_wavescan)
        assert table.is_wavelength_scan
        assert table.scan_units == ("nm",)
        assert table.wavelengths.tolist() == [1.0, 2.0, 3.0]

    def test_missing_marker_becomes_zero_not_nan(self, synthetic_wavescan):
        """'--' means evanescent: 0.0 with propagating=False, never NaN."""
        table = read_table(synthetic_wavescan)
        assert table.efficiency[2].tolist() == [0.0, 1.0, 0.0]
        assert table.propagating[2].tolist() == [False, True, False]
        assert not np.isnan(table.efficiency).any()

    def test_converts_to_scan(self, synthetic_wavescan):
        scan = read_scan(synthetic_wavescan)
        assert isinstance(scan, EfficiencyScan)
        assert len(scan) == 3
        assert scan.order(0).tolist() == [0.5, 0.8, 1.0]
        assert np.allclose(scan.total, [1.0, 1.0, 1.0])

    def test_provenance_is_populated(self, synthetic_wavescan):
        scan = read_scan(synthetic_wavescan)
        assert scan.provenance.method == "imported"
        assert scan.provenance.version == "PCGrate-SX 6.7.1"
        assert str(synthetic_wavescan) in scan.provenance.source
        # Imported data has not been convergence-checked by us.
        assert scan.provenance.converged is None
        assert not scan.provenance.is_defensible

    def test_order_lookup_is_by_index_not_position(self, synthetic_wavescan):
        scan = read_scan(synthetic_wavescan)
        assert scan.at(1.0)[-1] == 0.25
        assert scan.at(1.0)[0] == 0.5
        assert scan.at(1.0)[99] == 0.0  # absent order, not an error


class TestMethodLabelling:
    """The caller names the physics; the reader records the program.

    A table says what the efficiencies are, not which numerical method produced
    them, so inferring "integral" from a vendor banner would be an inference
    presented as a fact.
    """

    def test_caller_declares_the_method(self, synthetic_wavescan):
        scan = read_scan(synthetic_wavescan, method="integral")
        assert scan.provenance.method == "integral"

    def test_default_is_neutral(self, synthetic_wavescan):
        assert read_scan(synthetic_wavescan).provenance.method == "imported"

    def test_the_producing_program_is_still_recorded(self, synthetic_wavescan):
        """Neutral labelling must not cost attribution."""
        scan = read_scan(synthetic_wavescan, method="integral")
        assert scan.provenance.version == "PCGrate-SX 6.7.1"
        assert scan.provenance.notes["vendor"] == "PCGrate-SX"

    def test_label_flows_through_to_comparison_records(self, synthetic_wavescan):
        rows = read_scan(synthetic_wavescan, method="integral").to_records()
        assert {row["method"] for row in rows} == {"integral"}


class TestRejection:
    """Bad input must fail loudly. Silently guessing corrupts reference data."""

    def test_rejects_a_file_that_is_not_an_efficiency_table(self, tmp_path):
        path = tmp_path / "notpcgrate.txt"
        # Long enough to get past the length check, so the banner check is what fires.
        path.write_text("\n".join(f"{i}.0 {i}.5 {i}.9" for i in range(12)))
        with pytest.raises(ValueError, match="not a PCGrate banner"):
            read_table(path)

    def test_rejects_file_with_no_efficiency_columns(self, tmp_path):
        path = tmp_path / "empty.txt"
        path.write_text(
            '"PCGrate-SX 6.7.1 (c)1996-2020 I.I.G., Inc."\n\nP\n'
            "Solved at\t1/1/2026\t00:00:00\nCalculating time\t00:00:01\n\n"
            "\tScan Step\t\n P. ang., deg\t20.0\t\n"
        )
        with pytest.raises(ValueError, match="no 'Eff"):
            read_table(path)

    def test_rejects_truncated_file(self, tmp_path):
        path = tmp_path / "short.txt"
        path.write_text('"PCGrate-SX 6.7.1 (c)1996-2020 I.I.G., Inc."\n\nP\n')
        with pytest.raises(ValueError, match="too short"):
            read_table(path)

    def test_angle_scan_refuses_to_become_a_wavelength_scan(self, tmp_path):
        path = tmp_path / "angle.txt"
        path.write_text(
            '"PCGrate-SX 6.7.1 (c)1996-2020 I.I.G., Inc."\n\nP\n'
            "Solved at\t1/1/2026\t00:00:00\nCalculating time\t00:00:01\n\n"
            '\tScan Step\t"Eff.TE(0,R)"\n'
            " P. ang., deg Az. ang., deg\t20.0 89.0\t0.5\n"
        )
        table = read_table(path)
        assert not table.is_wavelength_scan
        assert table.scan_variables == ("P. ang.", "Az. ang.")
        with pytest.raises(ValueError, match="not a\n?\\s*single wavelength axis"):
            to_scan(table)


REF = reference_dir()
requires_corpus = pytest.mark.skipif(
    REF is None, reason="PCGrate reference corpus not available"
)
CORPUS = sorted(REF.rglob("*.txt")) if REF else []

# Files in the corpus that are legitimately not PCGrate efficiency exports.
NOT_EXPORTS = {
    "Au_CXRO_SXR.txt",  # CXRO optical constants
    "Au_optical_constants_n_k.txt",  # converted optical constants
    "pcgratetest3.txt",  # scan geometry only, no efficiency columns
}


@requires_corpus
class TestRealCorpus:
    """Verification item 1: every real export parses, with headers as ground truth."""

    @pytest.mark.parametrize("path", CORPUS, ids=lambda p: p.name)
    def test_parses_or_is_a_known_non_export(self, path):
        if path.name in NOT_EXPORTS:
            with pytest.raises(ValueError):
                read_table(path)
            return

        table = read_table(path)

        assert table.version.startswith("PCGrate")
        assert len(table.orders) > 0
        assert len(np.unique(table.orders)) == len(table.orders)
        assert table.efficiency.shape == (
            len(table.scan_values),
            len(table.orders),
        )
        assert table.propagating.shape == table.efficiency.shape
        assert not np.isnan(table.efficiency).any()
        assert (table.efficiency >= 0).all()
        # Evanescent cells must be exactly zero, never a stray value.
        assert (table.efficiency[~table.propagating] == 0.0).all()

    @pytest.mark.parametrize("path", CORPUS, ids=lambda p: p.name)
    def test_wavelength_scans_convert_and_conserve(self, path):
        if path.name in NOT_EXPORTS:
            pytest.skip("not an efficiency export")
        table = read_table(path)
        if not table.is_wavelength_scan:
            pytest.skip("angle scan")

        scan = to_scan(table)
        assert (scan.wavelengths > 0).all()
        assert np.all(np.diff(scan.wavelengths) > 0), "wavelength axis not ascending"
        assert (scan.efficiency >= 0).all()

    @pytest.mark.parametrize(
        "name", ["panter1_fix_efftable_perf.txt", "OGRE/tastetest_perf_wavescan.txt"]
    )
    def test_perfect_conductivity_conserves_energy(self, ref_dir, name):
        """Chapter-2.tex 2.2.2.3 proves sum(E_n) = 1 for a perfectly conducting
        grating, via Green's theorem. The perfect-conductivity exports must show
        it -- this is a physics theorem, not a tolerance to be relaxed.
        """
        path = ref_dir / name
        if not path.exists():
            pytest.skip(f"{name} not in corpus")

        scan = read_scan(path)
        worst = np.abs(scan.total - 1.0).max()
        assert worst < 0.01, (
            f"{name}: perfect-conductivity energy balance off by {worst:.4f}; "
            "sum over propagating orders must equal unity"
        )

    def test_finite_conductivity_absorbs_where_it_is_converged(self, ref_dir):
        """Where both runs are trustworthy, finite conductivity must lose energy."""
        perfect = ref_dir / "panter1_fix_efftable_perf.txt"
        finite = ref_dir / "panter1_fix_efftable_finite.txt"
        if not (perfect.exists() and finite.exists()):
            pytest.skip("panter1 pair not in corpus")

        p, f = read_scan(perfect), read_scan(finite)
        assert p.orders.tolist() == f.orders.tolist()
        assert np.allclose(p.wavelengths, f.wavelengths)

        trustworthy = f.total <= 1.0
        assert (f.total[trustworthy] < p.total[trustworthy]).all()

    def test_known_defect_finite_conductivity_run_is_unphysical(self, ref_dir):
        """`panter1_fix_efftable_finite.txt` violates energy conservation.

        Total efficiency reaches 3.6 -- impossible for absolute efficiency.
        Verified against the raw file, so this is a defect in the PCGrate run,
        not in this parser. The companion perfect-conductivity run of the same
        problem conserves energy to 4 decimal places, so the solver setup is
        sound; it is the finite-conductivity solve that destabilises.

        79 of 561 points (14%) violate, in five clusters::

            0.60-0.89 nm   30 pts   23-33 orders   peak 2.94
            0.97-1.03 nm    7 pts   19-20 orders   peak 1.23
            1.11-1.12 nm    2 pts      18 orders   peak 3.60
            1.44-1.49 nm    6 pts      14 orders   peak 2.03
            4.16-4.49 nm   34 pts       5 orders   peak 3.56

        The violations are strongly associated with **order passing-off
        (Rayleigh anomalies)**, not with propagating-order count: violating
        points sit a median 0.020 nm from a wavelength where the propagating
        order count changes, against 0.125 nm for clean points. The 4.16-4.49 nm
        cluster carries only five orders and terminates exactly at the 5->4
        passing-off wavelength. Au optical constants are smooth across all five
        bands, so no material edge is involved.

        Recorded as a test so it cannot be silently rediscovered or trusted.
        **Exclude this file from reference use at the listed wavelengths.**
        If the run is regenerated correctly this test fails; delete it then.
        """
        path = ref_dir / "panter1_fix_efftable_finite.txt"
        if not path.exists():
            pytest.skip("panter1 finite run not in corpus")

        scan = read_scan(path)
        bad = scan.total > 1.0

        assert bad.any(), (
            "the finite-conductivity run now conserves energy -- it was "
            "presumably regenerated. Delete this test."
        )
        assert 0.10 < bad.mean() < 0.20, (
            f"violating fraction is now {bad.mean():.2f}, previously 0.14; "
            "the file changed and the exclusion needs revisiting"
        )

        # The association with passing-off is the diagnostic worth locking in.
        n_orders = scan.propagating.sum(axis=1)
        passing_off = scan.wavelengths[:-1][np.diff(n_orders) != 0]
        distance = np.min(
            np.abs(scan.wavelengths[:, None] - passing_off[None, :]), axis=1
        )
        assert np.median(distance[bad]) < np.median(distance[~bad]) / 2, (
            "energy-conservation violations are no longer concentrated near "
            "order passing-off; the failure mode changed"
        )
