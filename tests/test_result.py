"""Result containers.

Every solver and importer flows through these, so the invariants here are what
make an energy-balance check mean anything.
"""

import numpy as np
import pytest

import gratinglab
from gratinglab.result import EfficiencyScan, OrderEfficiency, Provenance


def make_scan(**overrides) -> EfficiencyScan:
    """A 3-wavelength, 3-order scan where the last point has both +/-1 evanescent."""
    kwargs = dict(
        wavelengths=np.array([1.0, 2.0, 3.0]),
        orders=np.array([-1, 0, 1]),
        efficiency=np.array([[0.25, 0.5, 0.25], [0.1, 0.8, 0.1], [0.0, 1.0, 0.0]]),
        propagating=np.array([[True] * 3, [True] * 3, [False, True, False]]),
        provenance=Provenance("test"),
    )
    kwargs.update(overrides)
    return EfficiencyScan(**kwargs)


class TestImmutability:
    def test_arrays_are_read_only(self):
        scan = make_scan()
        for name in ("wavelengths", "orders", "efficiency", "propagating"):
            with pytest.raises(ValueError, match="read-only"):
                getattr(scan, name)[0] = 0

    def test_construction_copies_so_caller_arrays_stay_writable(self):
        """Freezing the caller's own array in place would be a nasty surprise."""
        efficiency = np.array([[0.5, 0.5]])
        EfficiencyScan(
            wavelengths=np.array([1.0]),
            orders=np.array([0, 1]),
            efficiency=efficiency,
            propagating=np.array([[True, True]]),
            provenance=Provenance("test"),
        )
        efficiency[0, 0] = 0.9  # must not raise
        assert efficiency[0, 0] == 0.9

    def test_mutating_the_source_does_not_alter_the_scan(self):
        efficiency = np.array([[0.5, 0.5]])
        scan = EfficiencyScan(
            wavelengths=np.array([1.0]),
            orders=np.array([0, 1]),
            efficiency=efficiency,
            propagating=np.array([[True, True]]),
            provenance=Provenance("test"),
        )
        efficiency[0, 0] = 0.9
        assert scan.efficiency[0, 0] == 0.5

    def test_row_views_are_also_read_only(self):
        with pytest.raises(ValueError, match="read-only"):
            make_scan().at(1.0).efficiency[0] = 0.0


class TestInvariants:
    def test_rejects_efficiency_on_a_non_propagating_order(self):
        """The guard that keeps `total` trustworthy."""
        with pytest.raises(ValueError, match="non-propagating"):
            make_scan(
                efficiency=np.array(
                    [[0.25, 0.5, 0.25], [0.1, 0.8, 0.1], [0.3, 1.0, 0.0]]
                )
            )

    def test_rejects_nan(self):
        with pytest.raises(ValueError, match="NaN"):
            make_scan(
                efficiency=np.array(
                    [[0.25, 0.5, 0.25], [0.1, 0.8, 0.1], [np.nan, 1.0, 0.0]]
                )
            )

    def test_rejects_negative_efficiency(self):
        with pytest.raises(ValueError, match="negative"):
            make_scan(
                efficiency=np.array(
                    [[-0.25, 0.5, 0.25], [0.1, 0.8, 0.1], [0.0, 1.0, 0.0]]
                )
            )

    def test_rejects_duplicate_orders(self):
        with pytest.raises(ValueError, match="duplicates"):
            make_scan(orders=np.array([0, 0, 1]))

    def test_rejects_shape_mismatch(self):
        with pytest.raises(ValueError, match="expected"):
            make_scan(efficiency=np.array([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]))

    def test_rejects_propagating_shape_mismatch(self):
        with pytest.raises(ValueError, match="propagating has shape"):
            make_scan(propagating=np.array([[True, True, True]]))


class TestAccess:
    def test_total_sums_only_propagating_orders(self):
        assert np.allclose(make_scan().total, [1.0, 1.0, 1.0])

    def test_order_returns_the_column_for_that_index(self):
        assert make_scan().order(0).tolist() == [0.5, 0.8, 1.0]

    def test_order_raises_for_an_absent_index(self):
        with pytest.raises(KeyError, match="order 7 not in this scan"):
            make_scan().order(7)

    def test_at_snaps_to_nearest_wavelength(self):
        assert make_scan().at(1.9).wavelength == 2.0

    def test_len_and_iteration(self):
        scan = make_scan()
        assert len(scan) == 3
        rows = list(scan)
        assert len(rows) == 3
        assert all(isinstance(r, OrderEfficiency) for r in rows)
        assert [r.wavelength for r in rows] == [1.0, 2.0, 3.0]

    def test_order_efficiency_getitem_returns_zero_for_absent_order(self):
        row = make_scan().at(1.0)
        assert row[0] == 0.5
        assert row[99] == 0.0

    def test_propagating_orders_excludes_evanescent(self):
        assert make_scan().at(3.0).propagating_orders.tolist() == [0]

    def test_energy_balance_error(self):
        assert make_scan().at(1.0).energy_balance_error() == pytest.approx(0.0)

    def test_to_records_is_complete_and_tidy(self):
        records = make_scan().to_records()
        assert len(records) == 9  # 3 wavelengths x 3 orders
        assert set(records[0]) == {
            "method",
            "wavelength",
            "order",
            "efficiency",
            "propagating",
            "converged",
        }
        # Evanescent cells survive as rows, so column j means order j everywhere.
        evanescent = [r for r in records if not r["propagating"]]
        assert len(evanescent) == 2
        assert all(r["efficiency"] == 0.0 for r in evanescent)


class TestProvenance:
    def test_is_hashable_and_usable_as_a_dict_key(self):
        """The comparison harness groups results by provenance."""
        a, b = Provenance("scalar"), Provenance("rcwa")
        assert len({a, b, Provenance("scalar")}) == 2
        assert {a: 1}[Provenance("scalar")] == 1

    def test_hashable_even_with_notes(self):
        hash(Provenance("rcwa", notes={"factorization": "li", "layers": [1, 2, 3]}))

    def test_notes_do_not_leak_between_instances(self):
        a, b = Provenance("a"), Provenance("b")
        a.notes["x"] = 1
        assert "x" not in b.notes

    def test_is_defensible_requires_demonstrated_convergence(self):
        assert not Provenance("rcwa").is_defensible
        assert not Provenance("rcwa", converged=None).is_defensible
        assert not Provenance("rcwa", converged=False).is_defensible
        assert Provenance("rcwa", converged=True).is_defensible

    def test_with_warning_appends_without_mutating(self):
        original = Provenance("scalar")
        warned = original.with_warning("past total external reflection")
        assert original.warnings == ()
        assert warned.warnings == ("past total external reflection",)
        assert warned.method == "scalar"


class TestPublicAPI:
    def test_documented_names_are_importable(self):
        """The README and docstrings promise these; a bare package would fail."""
        for name in ("Illumination", "EfficiencyScan", "Provenance", "__version__"):
            assert hasattr(gratinglab, name), f"gratinglab.{name} missing"

    def test_all_is_accurate(self):
        for name in gratinglab.__all__:
            assert hasattr(gratinglab, name)
