"""Measured efficiency: the data model, and its place in the comparison harness.

Three parts. The model's own invariants; measured data beside modelled scans
in :func:`align`, :func:`records` and :meth:`Comparison.residuals`; and a
worked comparison against the one fully ported ALS campaign, which needs the
private beamline tree and the ``als_eff_meas`` package and skips cleanly
without either.
"""

import numpy as np
import pytest

from gratinglab import MeasuredEfficiency, MeasurementProvenance
from gratinglab.compare import Residuals, align, records, sweep
from gratinglab.illumination import Illumination
from gratinglab.materials.optical import HC_EV_NM
from gratinglab.problem import Problem
from gratinglab.profiles import Lamellar
from gratinglab.result import EfficiencyScan, Provenance

PROBLEM = Problem(period=158.0, profile=Lamellar(depth_fraction=0.05))
ILL = Illumination.offplane(graze=1.5, azimuth=10.0, polarization="unpolarized")


def _measured(**provenance) -> MeasuredEfficiency:
    return MeasuredEfficiency(
        order=[1, 0, 0, 1],
        wavelength=[2.0, 2.0, 1.5, 1.5],
        value=[0.14, 0.64, 0.43, 0.23],
        sigma=[0.01, 0.02, 0.02, 0.01],
        provenance=MeasurementProvenance(**provenance),
    )


class TestModel:
    def test_points_are_stored_by_order_then_wavelength(self):
        data = _measured()
        assert data.order.tolist() == [0, 0, 1, 1]
        assert data.wavelength.tolist() == [1.5, 2.0, 1.5, 2.0]
        # value and sigma travel with their points, not with their positions
        assert data.value.tolist() == [0.43, 0.64, 0.23, 0.14]
        assert data.sigma.tolist() == [0.02, 0.02, 0.01, 0.01]

    def test_it_is_sparse_per_order(self):
        """The reason it is not an EfficiencyScan: orders need not share
        wavelengths, and nothing is padded with zeros."""
        data = MeasuredEfficiency(
            order=[0, 0, 0, 3], wavelength=[1.0, 2.0, 3.0, 1.0],
            value=[0.5, 0.6, 0.7, 0.01], sigma=[0.01] * 4,
        )
        assert len(data) == 4
        assert len(data.points(3)) == 1
        assert data.points(0).wavelengths.tolist() == [1.0, 2.0, 3.0]

    def test_arrays_are_read_only_copies(self):
        value = np.array([0.5])
        data = MeasuredEfficiency(order=[0], wavelength=[1.0], value=value, sigma=[0.1])
        value[0] = 9.0
        assert data.value[0] == 0.5
        with pytest.raises(ValueError):
            data.value[0] = 1.0

    def test_a_negative_value_is_kept_as_measured(self):
        """Clamping a noisy weak order to zero would bias every residual."""
        data = MeasuredEfficiency(order=[5], wavelength=[1.2], value=[-0.002], sigma=[0.003])
        assert data.value[0] == -0.002

    def test_from_energies_uses_the_materials_constant(self):
        data = MeasuredEfficiency.from_energies([0, 1], [500.0, 1000.0], [0.5, 0.1], [0.01, 0.01])
        assert data.points(0).wavelengths[0] == HC_EV_NM / 500.0
        assert np.allclose(np.sort(data.energy_ev), [500.0, 1000.0])

    def test_wavelength_grid_is_the_distinct_measured_wavelengths(self):
        assert _measured().wavelength_grid().tolist() == [1.5, 2.0]

    @pytest.mark.parametrize(
        "change, message",
        [
            ({"value": [0.1, np.nan]}, "NaN"),
            ({"sigma": [0.1, -0.01]}, "non-negative"),
            ({"wavelength": [1.0, 0.0]}, "positive"),
            ({"order": [0.0, 0.5]}, "non-integer"),
            ({"sigma": [0.1]}, "differ in length"),
            ({"wavelength": [1.0, 1.0], "order": [2, 2]}, "measured twice"),
        ],
    )
    def test_rejects_malformed_points(self, change, message):
        fields = {"order": [0, 1], "wavelength": [1.0, 2.0], "value": [0.1, 0.2], "sigma": [0.1, 0.1]}
        fields.update(change)
        with pytest.raises(ValueError, match=message):
            MeasuredEfficiency(**fields)

    def test_rejects_an_empty_measurement(self):
        with pytest.raises(ValueError, match="at least one point"):
            MeasuredEfficiency(order=[], wavelength=[], value=[], sigma=[])

    def test_an_unmeasured_order_is_a_key_error_not_an_empty_set(self):
        with pytest.raises(KeyError, match="not measured"):
            _measured().points(7)

    def test_provenance_rejects_a_polarization_outside_the_convention(self):
        with pytest.raises(ValueError, match="polarization"):
            MeasurementProvenance(polarization="s")

    def test_provenance_geometry_is_a_copy(self):
        geometry = {"gamma_deg": 1.5, "alpha_deg": 10.0}
        prov = MeasurementProvenance(geometry=geometry)
        geometry["gamma_deg"] = 80.0
        assert prov.geometry["gamma_deg"] == 1.5

    def test_illumination_comes_from_the_recorded_geometry(self):
        data = _measured(geometry={"gamma_deg": 1.5, "alpha_deg": 10.0}, polarization="TE")
        assert data.illumination() == Illumination.offplane(graze=1.5, azimuth=10.0)

    def test_illumination_refuses_to_guess_a_polarization(self):
        data = _measured(geometry={"gamma_deg": 1.5, "alpha_deg": 10.0})
        with pytest.raises(ValueError, match="not recorded"):
            data.illumination()
        assert data.illumination("TM").polarization == "TM"

    def test_illumination_names_the_missing_angle(self):
        with pytest.raises(ValueError, match="alpha_deg"):
            _measured(geometry={"gamma_deg": 1.5}).illumination("TE")

    def test_records_carry_sigma(self):
        rows = _measured(label="ALS").to_records()
        assert len(rows) == 4
        assert rows[0] == {
            "method": "ALS", "wavelength": 1.5, "order": 0, "efficiency": 0.43,
            "sigma": 0.02, "propagating": True, "converged": None,
        }


class TestInTheComparison:
    def setup_method(self):
        self.data = _measured()
        self.scalar = sweep(PROBLEM, ILL, self.data.wavelength_grid(), ["scalar"])[0]

    def test_sweep_refuses_a_measurement(self):
        """Resampling sparse data onto a grid would invent points."""
        with pytest.raises(TypeError, match="align"):
            sweep(PROBLEM, ILL, [1.5], ["scalar", self.data])

    def test_measured_and_modelled_share_one_comparison(self):
        comparison = align([self.data, self.scalar])
        assert comparison.methods == ("scalar",)
        assert comparison.measurements == ("measured",)
        assert comparison.efficiency.shape[0] == 1
        points = comparison.measured_order(1)["measured"]
        assert points.value.tolist() == [0.23, 0.14]
        assert points.sigma.tolist() == [0.01, 0.01]

    def test_measured_orders_join_the_union(self):
        lone = MeasuredEfficiency(order=[9], wavelength=[1.5], value=[0.001], sigma=[0.001])
        comparison = align([self.scalar, lone])
        assert 9 in comparison.orders
        # a model that does not resolve order 9 reports zero there
        assert comparison.order(9)["scalar"].tolist() == [0.0, 0.0]

    def test_a_measurement_alone_still_aligns(self):
        comparison = align([self.data])
        assert comparison.methods == ()
        assert comparison.orders.tolist() == [0, 1]
        assert comparison.efficiency.shape == (0, 0, 2)

    def test_two_measurements_are_told_apart(self):
        comparison = align([self.data, self.data, self.scalar])
        assert comparison.measurements == ("measured", "measured#2")

    def test_measured_order_skips_a_measurement_without_that_order(self):
        other = MeasuredEfficiency(order=[0], wavelength=[1.5], value=[0.4], sigma=[0.02])
        comparison = align([self.data, other, self.scalar])
        assert set(comparison.measured_order(1)) == {"measured"}

    def test_align_refuses_something_that_is_neither(self):
        with pytest.raises(TypeError, match="MeasuredEfficiency"):
            align([self.scalar, {"order": [0]}])

    def test_records_stack_measured_and_modelled_rows(self):
        rows = records([self.scalar, self.data])
        assert len(rows) == self.scalar.efficiency.size + len(self.data)
        measured = [r for r in rows if r["method"] == "measured"]
        assert len(measured) == 4 and all("sigma" in r for r in measured)


class TestResiduals:
    def _synthetic(self, scan: EfficiencyScan, offset: float, sigma: float):
        """Measurements built from a model's own output, so the expected
        residual is known exactly."""
        rows = [
            (m, w, v)
            for w, row, prop in zip(scan.wavelengths, scan.efficiency, scan.propagating)
            for m, v, p in zip(scan.orders, row, prop)
            if p and v > 0.05
        ]
        order, wavelength, value = (np.array(c) for c in zip(*rows))
        return MeasuredEfficiency(
            order=order, wavelength=wavelength, value=value + offset,
            sigma=np.full(len(order), sigma),
        )

    def setup_method(self):
        self.grid = np.linspace(1.4, 2.6, 7)
        self.scalar = sweep(PROBLEM, ILL, self.grid, ["scalar"])[0]

    def test_a_model_against_its_own_output_has_zero_residual(self):
        data = self._synthetic(self.scalar, 0.0, 0.01)
        res = align([data, self.scalar]).residuals("measured", "scalar")
        assert isinstance(res, Residuals)
        assert len(res) == len(data)
        assert np.array_equal(res.residual, np.zeros(len(data)))
        assert res.propagating.all()

    def test_pull_is_the_residual_in_units_of_sigma(self):
        data = self._synthetic(self.scalar, 0.02, 0.01)
        res = align([data, self.scalar]).residuals("measured", "scalar")
        assert np.allclose(res.pull, 2.0)
        summary = res.summary()
        assert summary["chi2"] == pytest.approx(4.0 * len(data))
        assert summary["n_with_sigma"] == len(data)
        assert summary["max_abs_residual"] == pytest.approx(0.02)

    def test_a_zero_sigma_gives_no_pull_and_is_left_out_of_chi2(self):
        data = MeasuredEfficiency(
            order=[0, 0], wavelength=self.grid[:2],
            value=self.scalar.order(0)[:2] + 0.01, sigma=[0.0, 0.01],
        )
        res = align([data, self.scalar]).residuals("measured", "scalar")
        assert np.isnan(res.pull[0]) and res.pull[1] == pytest.approx(1.0)
        assert res.summary()["n_with_sigma"] == 1
        assert res.summary()["chi2"] == pytest.approx(1.0)

    def test_an_order_the_model_calls_evanescent_reads_zero(self):
        """A disagreement worth seeing, so it is reported, not dropped."""
        data = MeasuredEfficiency(order=[-2], wavelength=[self.grid[-1]], value=[0.01], sigma=[0.005])
        res = align([data, self.scalar]).residuals("measured", "scalar")
        assert res.model.tolist() == [0.0]
        assert res.propagating.tolist() == [False]

    def test_off_grid_measurements_are_refused_not_interpolated(self):
        data = MeasuredEfficiency(order=[0], wavelength=[1.45], value=[0.4], sigma=[0.01])
        with pytest.raises(ValueError, match="wavelength_grid"):
            align([data, self.scalar]).residuals("measured", "scalar")

    def test_the_union_grid_serves_both_a_curve_and_the_residuals(self):
        data = MeasuredEfficiency(order=[0], wavelength=[1.45], value=[0.4], sigma=[0.01])
        grid = np.union1d(self.grid, data.wavelength_grid())
        scan = sweep(PROBLEM, ILL, grid, ["scalar"])[0]
        res = align([data, scan]).residuals("measured", "scalar")
        assert res.model[0] == scan.order(0)[np.flatnonzero(grid == 1.45)[0]]

    def test_unknown_names_raise(self):
        comparison = align([_measured(), sweep(PROBLEM, ILL, [1.5, 2.0], ["scalar"])[0]])
        with pytest.raises(ValueError):
            comparison.residuals("nope", "scalar")
        with pytest.raises(ValueError):
            comparison.residuals("measured", "integral")


class TestImportedTableBesideAMeasurement:
    def test_an_imported_scan_and_a_measurement_share_a_comparison(self):
        """The three-way plot the split is for, with a stand-in for the
        imported reference so it runs without the corpus."""
        grid = np.array([1.5, 2.0])
        reference = EfficiencyScan(
            wavelengths=grid, orders=np.array([0, 1]),
            efficiency=np.array([[0.45, 0.22], [0.62, 0.15]]),
            propagating=np.ones((2, 2), dtype=bool),
            provenance=Provenance("integral"),
        )
        scans = sweep(PROBLEM, ILL, grid, ["scalar", reference])
        comparison = align([*scans, _measured()])
        assert comparison.methods == ("scalar", "integral")
        res = comparison.residuals("measured", "integral")
        assert np.allclose(res.residual, [0.43 - 0.45, 0.64 - 0.62, 0.23 - 0.22, 0.14 - 0.15])


# -- the worked comparison ---------------------------------------------------


@pytest.fixture(scope="module")
def ogre_measured():
    als = pytest.importorskip("als_eff_meas")
    campaign = als.load_campaign("ogre_laminar_I_oct21")
    if not campaign.data.is_dir():
        pytest.skip(f"beamline data not available at {campaign.data}; set ALS_DATA_ROOT")
    import warnings

    with warnings.catch_warnings():
        # ODR probes outside the arc while converging; harmless, as
        # als_eff_meas's own conftest notes.
        warnings.simplefilter("ignore", RuntimeWarning)
        geometry = als.fit_geometry(campaign)
    table = als.measure_efficiency(campaign, geometry)

    order, energy, value, sigma = [], [], [], []
    for n in table.orders:
        for e, eff in zip(table.energy[n], table.efficiency[n], strict=True):
            order.append(n)
            energy.append(float(e))
            value.append(eff.nominal_value)
            sigma.append(eff.std_dev)
    angles = geometry.angles.degrees
    return MeasuredEfficiency.from_energies(
        order, energy, value, sigma,
        provenance=MeasurementProvenance(
            facility="ALS",
            campaign=campaign.name,
            polarization=None,
            geometry={
                "gamma_deg": _nominal(angles["gamma"]),
                "alpha_deg": _nominal(angles["alpha"]),
                "eta_deg": _nominal(angles["eta"]),
                "yaw_deg": _nominal(angles["yaw"]),
            },
        ),
    )


class TestOgreLaminarCampaign:
    """OGRE laminar grating I, ALS October 2021, against the scalar model.

    This needs ``als_eff_meas`` importable and its private beamline tree at
    ``$ALS_DATA_ROOT`` (default ``~/Documents/ALS``), and skips without
    either. Checked 2026-10-24: neither MC1 nor MAC1 holds that tree, so it
    has not yet run anywhere.

    Until ``als_eff_meas`` grows its ``to_measured()``, the conversion its
    side will own is written out here, so the test also serves as the
    specification for it.

    Agreement is not asserted. The groove depth and duty cycle of this
    grating are not in the campaign, so the profile below is a placeholder
    and the scalar curve is not a fit; the polarization is not recorded
    either. What is asserted is everything that does not depend on them:
    every reduced point arrives with its sigma, the comparison holds the
    measurement and the model together, and the recovered geometry makes
    every measured order propagate under gratinglab's grating equation.
    """

    def test_every_reduced_point_arrives_with_its_sigma(self, ogre_measured):
        # 249 plotted points across nine orders, plus order 7's two, which
        # are measured and were never plotted (als_eff_meas's baseline test)
        assert len(ogre_measured) == 251
        assert (ogre_measured.sigma > 0).all()
        assert sorted(ogre_measured.orders.tolist()) == list(range(-2, 8))

    def test_measured_and_scalar_on_one_comparison(self, ogre_measured):
        illumination = _oriented(ogre_measured)
        scan = sweep(PROBLEM, illumination, ogre_measured.wavelength_grid(), ["scalar"])[0]
        comparison = align([ogre_measured, scan])
        res = comparison.residuals("measured", "scalar")
        assert len(res) == len(ogre_measured)
        assert np.isfinite(res.residual).all()
        assert set(comparison.measured_order(1)) == {"measured"}

    def test_the_recovered_geometry_propagates_every_measured_order(self, ogre_measured):
        """Depends only on the period and the angles, so it is a real check
        on the geometry whatever the profile."""
        illumination = _oriented(ogre_measured)
        scan = sweep(PROBLEM, illumination, ogre_measured.wavelength_grid(), ["scalar"])[0]
        res = align([ogre_measured, scan]).residuals("measured", "scalar")
        stranded = sorted({(int(m), round(float(w), 3)) for m, w, p in
                           zip(res.order, res.wavelength, res.propagating) if not p})
        assert not stranded, f"measured but evanescent here (order, nm): {stranded[:8]}"


def _nominal(value) -> float:
    return float(getattr(value, "nominal_value", value))


def _oriented(measured: MeasuredEfficiency) -> Illumination:
    """The illumination, with the sign of alpha settled by the data.

    als_eff_meas recovers alpha as ``arcsin(dx_dir / R)`` in its detector
    frame, and whether that frame's positive azimuth is gratinglab's has not
    been established. The two signs propagate different order ranges, so the
    one under which every measured order propagates is taken. If neither
    does, the first is returned and the propagation test reports it.
    """
    base = measured.illumination("unpolarized")
    candidates = [base, base.model_copy(update={"alpha_deg": -base.alpha_deg})]
    for candidate in candidates:
        scan = sweep(PROBLEM, candidate, measured.wavelength_grid(), ["scalar"])[0]
        if align([measured, scan]).residuals("measured", "scalar").propagating.all():
            return candidate
    return base
