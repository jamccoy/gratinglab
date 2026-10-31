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
from gratinglab.profiles import FromProfileData, Lamellar
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

    def test_no_warning_when_every_measured_order_propagates(self):
        data = self._synthetic(self.scalar, 0.0, 0.01)
        res = align([data, self.scalar]).residuals("measured", "scalar")
        assert res.evanescent_warning() is None
        assert res.summary()["n_evanescent"] == 0

    def test_stranded_points_are_counted_per_order_against_that_order(self):
        """The count per order, out of that order's own points, is what
        tells an edge-of-arc point from a whole order lost to a wrong mount."""
        zero = self.scalar.order(0)
        data = MeasuredEfficiency(
            order=[0, 0, -2, -2],
            wavelength=[self.grid[0], self.grid[1], self.grid[-1], self.grid[-2]],
            value=[zero[0], zero[1], 0.01, 0.01],
            sigma=[0.01] * 4,
        )
        res = align([data, self.scalar]).residuals("measured", "scalar")
        # stored by order, then wavelength, so order -2 comes first
        assert res.propagating.tolist() == [False, False, True, True]
        assert res.summary()["n_evanescent"] == 2
        message = res.evanescent_warning()
        assert message.startswith(
            "2 of 4 measured points in 'measured' fall on orders 'scalar' "
            "does not propagate"
        )
        assert "m=-2 (2 of 2)" in message
        assert "m=+0" not in message
        assert "sign of alpha" in message

    def test_the_wrong_sign_of_alpha_strands_whole_orders(self):
        """The SCIL failure in miniature: order +2 measured at an alpha of
        +30 degrees, where it propagates at every wavelength here, and the
        model handed -30 degrees, where it propagates at none of them."""
        mount = Illumination.offplane(graze=1.5, azimuth=30.0, polarization="unpolarized")
        right = sweep(PROBLEM, mount, self.grid, ["scalar"])[0]
        data = MeasuredEfficiency(
            order=[2] * len(self.grid), wavelength=self.grid,
            value=[0.01] * len(self.grid), sigma=[0.005] * len(self.grid),
        )
        assert align([data, right]).residuals("measured", "scalar").evanescent_warning() is None

        flipped = mount.model_copy(update={"alpha_deg": -mount.alpha_deg})
        wrong = sweep(PROBLEM, flipped, self.grid, ["scalar"])[0]
        res = align([data, wrong]).residuals("measured", "scalar")
        assert res.summary()["n_evanescent"] == len(self.grid)
        assert "m=+2 (7 of 7)" in res.evanescent_warning()

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
    ``$ALS_DATA_ROOT`` (default ``~/Documents/Research/ALS``), and skips without
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

    def test_als_eff_meas_recovers_alpha_in_this_projects_sign_convention(
        self, ogre_measured
    ):
        """The roadmap's open question, settled by the data.

        ``als_eff_meas`` recovers alpha in its own detector frame, and whether
        that frame's positive azimuth is :attr:`Illumination.alpha_deg`'s was
        never established. It is: the sign as recovered propagates every
        measured point, and the flip strands 107 of the 251. The margin is the
        whole point -- this is not a tolerance question, so a future change to
        either convention fails here loudly rather than quietly biasing every
        residual. See ``findings.md``, "The measured comparison ran".
        """
        grid = ogre_measured.wavelength_grid()
        as_recovered = ogre_measured.illumination("unpolarized")
        assert as_recovered.alpha_deg > 0.0

        def stranded(illumination):
            scan = sweep(PROBLEM, illumination, grid, ["scalar"])[0]
            res = align([ogre_measured, scan]).residuals("measured", "scalar")
            return int(np.count_nonzero(~res.propagating))

        flipped = as_recovered.model_copy(
            update={"alpha_deg": -as_recovered.alpha_deg}
        )
        assert stranded(as_recovered) == 0
        assert stranded(flipped) == 107


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


# -- the TASTE campaign, where the grating's groove is known -----------------


#: The groove the paper itself fitted: "an acute trapezoid with a near-vertical
#: lateral side opposite a slope that emulates the active blaze facet", plus a
#: flat bottom for the cleared resist, with the best-matching member at blaze
#: 27 deg, depth 120 nm, flat-top width 77 nm and bottom width 85 nm
#: (McCoy, McEntaffer & Miles 2020, section 4). Published geometry, so a test
#: can use it where the AFM scan it agrees with is not in this repository.
TASTE_PERIOD_NM = 400.0


def _taste_trapezoid(
    blaze_deg: float = 27.0, depth_nm: float = 120.0,
    top_nm: float = 77.0, bottom_nm: float = 85.0,
) -> FromProfileData:
    """Valley, bottom land, steep side, apex land, blaze facet, valley.

    That order -- steep side rising, blaze facet descending -- is the one
    ``metrology.boundary`` exports for this grating, i.e. the ``+d-hat`` axis,
    so a constructed groove and a measured one are the same way round. It
    crosses into ``t``, which runs against ``+d-hat``, through
    ``along_dispersion``: build the same points with ``t=`` instead and the
    groove is solved as its mirror image, with the blaze order on the wrong
    side of zero. See ``docs/conventions.md``, "Direction of the profile
    parameter", and ``docs/findings.md``, "A measured groove has no handedness
    of its own".
    """
    facet_nm = depth_nm / np.tan(np.radians(blaze_deg))
    steep_nm = TASTE_PERIOD_NM - bottom_nm - facet_nm - top_nm
    assert steep_nm > 0, f"only {steep_nm:.1f} nm left for the steep side"
    x = np.array([0.0, bottom_nm, bottom_nm + steep_nm,
                  bottom_nm + steep_nm + top_nm, TASTE_PERIOD_NM])
    y = np.array([0.0, 0.0, depth_nm, depth_nm, 0.0])
    grid = np.linspace(0.0, TASTE_PERIOD_NM, 801)
    return FromProfileData.along_dispersion(
        grid / TASTE_PERIOD_NM,
        np.interp(grid, x, y) / TASTE_PERIOD_NM,
        blaze_angle=blaze_deg,
    )


def _taste_problem() -> Problem:
    """The paper's groove at the lithographic period, with no coating.

    No coating on purpose. The campaign reaches 80 eV, i.e. 15.5 nm, and the
    only optical-constant table this repository carries starts at 200 eV, so
    asking for absolute efficiency across this band would raise rather than
    extrapolate -- correctly (``materials/optical.py``). Every assertion in
    this class is about which orders exist and whether the records line up, and
    none of them needs the absolute level. The absolute comparison, which does,
    is in ``docs/findings.md`` and was run against a wider CXRO export supplied
    through ``GRATINGLAB_MATERIALS``.
    """
    return Problem(period=TASTE_PERIOD_NM, profile=_taste_trapezoid())


#: Table 1 of the paper, measured at the beamline: the cone half-angle and the
#: azimuthal incidence angle, with their uncertainties.
TASTE_GAMMA_DEG, TASTE_GAMMA_SIGMA = 1.73, 0.03
TASTE_ALPHA_DEG, TASTE_ALPHA_SIGMA = 23.4, 0.6

TASTE_ROUNDS = ("taste_dec18_round3a", "taste_dec18_round3b")


@pytest.fixture(scope="module")
def taste():
    """The two halves of the December 2018 TASTE campaign, reduced and joined.

    Needs ``als_eff_meas`` with both TASTE campaigns registered and its private
    beamline tree at ``$ALS_DATA_ROOT``, and skips without either. The two
    halves share no (order, wavelength) pair -- 3a runs 440-800 eV and 3b runs
    80-420 eV -- so joining them cannot collide.

    Returns ``(MeasuredEfficiency, {round: recovered angles})``. The angles come
    back because the two rounds disagree about the mount and that disagreement
    is the subject of a test below.
    """
    als = pytest.importorskip("als_eff_meas")
    available = set(als.list_campaigns())
    missing = [name for name in TASTE_ROUNDS if name not in available]
    if missing:
        pytest.skip(
            f"als_eff_meas has no campaign {missing}; available: "
            f"{sorted(available)}"
        )
    import warnings

    order, energy, value, sigma = [], [], [], []
    angles = {}
    for name in TASTE_ROUNDS:
        campaign = als.load_campaign(name)
        if not campaign.data.is_dir():
            pytest.skip(f"beamline data not available at {campaign.data}; "
                        "set ALS_DATA_ROOT")
        with warnings.catch_warnings():
            # ODR probes outside the arc while converging; harmless, as
            # als_eff_meas's own conftest notes.
            warnings.simplefilter("ignore", RuntimeWarning)
            geometry = als.fit_geometry(campaign)
        table = als.measure_efficiency(campaign, geometry)
        for n in table.orders:
            for e, eff in zip(table.energy[n], table.efficiency[n], strict=True):
                order.append(n)
                energy.append(float(e))
                value.append(eff.nominal_value)
                sigma.append(eff.std_dev)
        angles[name] = {k: _nominal(v) for k, v in geometry.angles.degrees.items()}

    measured = MeasuredEfficiency.from_energies(
        order, energy, value, sigma,
        provenance=MeasurementProvenance(
            facility="ALS",
            campaign="taste_dec18",
            polarization=None,
            # Table 1, not either round's own fit: the paper reports one mount
            # for the whole campaign and this is it.
            geometry={"gamma_deg": TASTE_GAMMA_DEG, "alpha_deg": TASTE_ALPHA_DEG},
            source="McCoy, McEntaffer & Miles 2020, ApJ 891, 114",
            notes={"rounds": list(TASTE_ROUNDS), "reduction": "als_eff_meas"},
        ),
    )
    return measured, angles


class TestTasteCampaign:
    """The TASTE blazed grating, ALS December 2018, against its own geometry.

    The campaign the OGRE one could not be: this grating's groove was measured
    by AFM, so the comparison is a prediction rather than a two-parameter fit.
    The AFM scans are not in this repository, so the profile used here is the
    paper's own published trapezoid, which anyone can rebuild from the paper.

    Agreement is not asserted; ``docs/findings.md`` records what it is. What is
    asserted is everything that does not depend on the groove.
    """

    def test_every_reduced_point_arrives_with_its_sigma(self, taste):
        measured, _ = taste
        # 208 from round 3a (orders -4 to 8) plus 88 from round 3b (-2 to 5)
        assert len(measured) == 296
        assert (measured.sigma > 0).all()
        assert measured.orders.tolist() == list(range(-4, 9))
        assert measured.energy_ev.min() == pytest.approx(80.0)
        assert measured.energy_ev.max() == pytest.approx(800.0)

    def test_round_3a_recovers_the_published_table_1(self, taste):
        """The paper's Table 1 is round 3a's arc fit, and this proves it.

        The paper says the arc was mapped "using data gathered at 450 and
        500 eV", which are exactly round 3a's two yaw-scan energies. Both
        angles land inside the published uncertainties, which also pins the
        sign convention: als_eff_meas returns the magnitude the original
        printed with a leading minus.
        """
        _, angles = taste
        a = angles["taste_dec18_round3a"]
        assert a["gamma"] == pytest.approx(TASTE_GAMMA_DEG, abs=TASTE_GAMMA_SIGMA)
        assert a["alpha"] == pytest.approx(TASTE_ALPHA_DEG, abs=TASTE_ALPHA_SIGMA)
        assert a["alpha"] > 0.0

    def test_round_3b_is_a_different_mount_and_says_so(self, taste):
        """The order sorter moved the beam, and the low-energy half knows it.

        The paper allows for this -- the sorter "is expected to slightly shift
        the position of the beam on the grating, and hence the measured
        parameters listed in Table 1" -- and calls the effect small. On gamma it
        is not small: 2.21 deg against 1.73, sixteen times the published sigma.
        Any model of the 80-420 eV half at Table 1's mount carries that.
        """
        _, angles = taste
        b = angles["taste_dec18_round3b"]
        assert b["gamma"] > TASTE_GAMMA_DEG + 10 * TASTE_GAMMA_SIGMA
        assert abs(b["alpha"] - TASTE_ALPHA_DEG) > 4 * TASTE_ALPHA_SIGMA

    def test_measured_and_modelled_on_one_comparison(self, taste):
        measured, _ = taste
        problem = _taste_problem()
        illumination = measured.illumination("TE")
        scan = sweep(problem, illumination, measured.wavelength_grid(),
                     ["scalar"])[0]
        comparison = align([measured, scan])
        res = comparison.residuals("measured", "scalar")
        assert len(res) == len(measured)
        assert np.isfinite(res.residual).all()
        assert comparison.measurements == ("measured",)

    def test_two_measured_points_are_evanescent_at_the_published_mount(self, taste):
        """Named, because they are evidence and not noise.

        At Table 1's mount the grating equation makes order -4 at 680 eV and
        order -2 at 340 eV evanescent, and the campaign has a measured point
        for each. Both are the most negative order their energy recorded, at
        the edge of the arc, and order -2 at 340 eV is the one point of round 3b
        whose reduction is sensitive to the radius correction. The model is not
        wrong to call them evanescent and the reduction is not wrong to have
        found something there -- but one of the two is mistaken, and this test
        fails loudly if their number ever grows.
        """
        measured, _ = taste
        problem = _taste_problem()
        scan = sweep(problem, measured.illumination("TE"),
                     measured.wavelength_grid(), ["scalar"])[0]
        res = align([measured, scan]).residuals("measured", "scalar")
        stranded = sorted(
            (int(m), round(float(HC_EV_NM / w)))
            for m, w, p in zip(res.order, res.wavelength, res.propagating)
            if not p
        )
        assert stranded == [(-4, 680), (-2, 340)]


class TestTheMeasuredGroovesHandedness:
    """Which way the groove faces the beam, and why it is not the file's to say.

    An AFM scan fixes a shape, not an orientation: nothing in the scan records
    which way the sample was mounted at the beamline, or the sign of the fast
    scan axis against the dispersion direction. The convention does fix it --
    a ``.ggp`` and the metrology pipeline's ``x_norm`` run along ``+d-hat``,
    ``t`` runs against it, and ``FromProfileData.along_dispersion`` is the one
    place the two meet (``docs/conventions.md``). Skipping that handover is not
    a cosmetic slip: it sends the blaze order to the wrong side of the arc.

    Needs no measurement and no optical constants, so unlike the campaign tests
    above this one always runs.
    """

    GRID = [2.8178]   # 440 eV, where the campaign's blaze order is n = +4

    def _orders(self, profile):
        problem = Problem(period=TASTE_PERIOD_NM, profile=profile)
        illumination = Illumination(
            gamma_deg=TASTE_GAMMA_DEG, alpha_deg=TASTE_ALPHA_DEG, polarization="TE"
        )
        scan = sweep(problem, illumination, self.GRID, ["integral"],
                     {"integral": {"boundary_points": 256}})[0]
        return {int(m): float(v)
                for m, v, p in zip(scan.orders, scan.efficiency[0],
                                   scan.propagating[0]) if p}

    def _mirrored(self) -> FromProfileData:
        """The same points put straight into ``t``, which is the mistake."""
        correct = _taste_trapezoid()
        x, y = correct.dispersion_axis()
        return FromProfileData(t=tuple(x), y=tuple(y), blaze_angle=27.0)

    def test_the_blaze_order_is_where_the_campaign_found_it(self):
        orders = self._orders(_taste_trapezoid())
        assert max(orders, key=orders.get) == 4
        assert orders[4] > 0.4

    def test_the_mirror_image_moves_the_power_off_the_blaze_order(self):
        """The failure this pins: a plausible curve with the wrong shape.

        Mirrored, order +4 keeps 17% of what it had and order -2 grows sixfold,
        so the model still looks like a blazed grating -- it just blazes
        somewhere the grating does not. A comparison against measured data
        would read as a profile error, and would be a convention error.
        """
        def negative(orders):
            return sum(v for m, v in orders.items() if m < 0)

        right, wrong = self._orders(_taste_trapezoid()), self._orders(self._mirrored())
        assert wrong[4] < 0.25 * right[4]
        assert negative(wrong) > 2.5 * negative(right)

    def test_along_dispersion_round_trips(self):
        """The handover is invertible, so neither direction is special-cased."""
        profile = _taste_trapezoid()
        x, y = profile.dispersion_axis()
        again = FromProfileData.along_dispersion(x, y, blaze_angle=profile.blaze_angle)
        assert np.allclose(again.t, profile.t)
        assert np.allclose(again.y, profile.y)
