"""A measured groove under the finite-conductivity integral solver (issue #6).

A measured profile is a polygon: :class:`FromProfileData` joins its points by
straight segments, and an AFM export has thousands of them. Every vertex is a
real corner of the metal interface, nearly all of them too slight to cross
:data:`~gratinglab.profiles.CORNER_TURN_DEGREES`, so none is graded towards and
the equal-arc mesh never resolves them. The solved boundary data are then wrong
next to every vertex. The reflected orders survive, because they pass through
the Green's representation, which integrates a localised error away, but the
absorption is a pointwise product of trace and normal derivative and does not.
On the TASTE grating's measured groove at 80 eV that put R + A at 1.08.

So the solver meshes such a polygon as a smooth curve by default
(``profile_smoothing="auto"``, ``docs/theory/integral.md`` section 8). These
tests hold the fix to the C-method, which treats the polygon exactly, and hold
everything else to being untouched.

The groove is built here rather than read from a scan, because the measurement
data are not committed (``tests/metrology/scans.py``). It is made the way an AFM
export is made, which is what reproduces the fault: an etched trapezoid (flat
lands, a 66 deg steep side, a shallow facet, the TASTE shape), rounded, sampled
at 102 pixels a period (3.9 nm, the TASTE scan's pixel), given 0.05 nm of
height noise, and linearly resampled to 2000 points. The kinks sit at the
pixels and the largest turns by 8.6 deg, under the threshold, so the polygon
has no corner. Measured 2026-10-29 with ``profile_smoothing=None`` (the polygon
as given), R + A - 1 at 1024 points was 7.3e-3 for gold at 3.1 nm, 2.1e-2 at
6.2 nm and 0.53 at an 80 eV-like index, while R stayed within 4e-4 of the
C-method throughout.
"""

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import (
    Blazed,
    FromProfileData,
    Sinusoidal,
    polyline_has_corners,
    smoothed_polyline,
)
from gratinglab.solvers import get_solver

PERIOD = 400.0
MOUNT = Illumination.offplane(graze=1.73, azimuth=23.4, polarization="TE")

#: An index like gold's at 80 eV (|n| = 0.888), where the fault was largest.
#: The vendored CXRO table starts at 200 eV, so a constant table supplies it.
LOW_ENERGY = "GoldLike80eV"
LOW_ENERGY_DELTA, LOW_ENERGY_BETA = 0.115, 0.075
LOW_ENERGY_NM = 15.5

#: (material, wavelength nm): |n| = 0.994, 0.989 and 0.888.
INDICES = [("Au", 3.1), ("Au", 6.199), (LOW_ENERGY, LOW_ENERGY_NM)]


def measured_groove(pixels=102, points=2000, noise_nm=0.05, round_nm=16.0, seed=6):
    """An etched trapezoid as an AFM export would give it (module docstring)."""
    depth, bottom, top = 100.0, 90.0, 77.5
    steep = depth / np.tan(np.radians(66.0))
    corners_x = np.array([0.0, bottom, bottom + steep, bottom + steep + top, PERIOD])
    corners_y = np.array([0.0, 0.0, depth, depth, 0.0])
    fine = 4096
    grid = np.arange(fine) * PERIOD / fine
    frequency = np.fft.rfftfreq(fine, d=PERIOD / fine)
    rounded = np.fft.irfft(
        np.fft.rfft(np.interp(grid, corners_x, corners_y))
        * np.exp(-2.0 * np.pi**2 * round_nm**2 * frequency**2),
        n=fine,
    )
    pixel_x = np.arange(pixels) * PERIOD / pixels
    pixel_y = np.interp(pixel_x, grid, rounded, period=PERIOD)
    pixel_y += noise_nm * np.random.default_rng(seed).standard_normal(pixels)
    t = np.linspace(0.0, 1.0, points)
    y = np.interp(t * PERIOD, pixel_x, pixel_y, period=PERIOD)
    return FromProfileData.along_dispersion(t, (y - y[0]) / PERIOD, blaze_angle=28.0)


GROOVE = measured_groove()


def reflected(scan):
    return np.nansum(np.where(scan.propagating, scan.efficiency, 0.0), axis=1)


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    """A constant-index CXRO table for the 80 eV-like point, beside the
    vendored gold, for the length of the module."""
    directory = tmp_path_factory.mktemp("materials")
    rows = [f" {LOW_ENERGY} Density=19.32", " Energy(eV), Delta, Beta"]
    rows += [
        f" {energy:.4f} {LOW_ENERGY_DELTA} {LOW_ENERGY_BETA}"
        for energy in np.linspace(60.0, 150.0, 40)
    ]
    (directory / f"{LOW_ENERGY}.txt").write_bytes(("\n".join(rows) + "\n").encode())
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("GRATINGLAB_MATERIALS", str(directory))
        yield directory


def integral_scan(material, wavelength, points, smoothing="auto", profile=GROOVE):
    problem = Problem(period=PERIOD, profile=profile, coating=material)
    return get_solver("integral").solve(
        problem,
        MOUNT,
        [wavelength],
        boundary_points=points,
        conductivity="tabulated",
        profile_smoothing=smoothing,
    )


@pytest.fixture(scope="module")
def cmethod_reflectance(library):
    """The C-method's R on the polygon itself, converged: its exact polygon
    path and the energy balance agree to 1e-7 at truncation 40."""
    out = {}
    for material, wavelength in INDICES:
        problem = Problem(period=PERIOD, profile=GROOVE, coating=material)
        scan = get_solver("cmethod").solve(
            problem, MOUNT, [wavelength], truncation=40, conductivity="tabulated"
        )
        assert abs(reflected(scan)[0] + scan.absorption[0] - 1.0) < 1e-6
        out[material, wavelength] = reflected(scan)[0]
    return out


@pytest.fixture(scope="module")
def ladder(library):
    """The default solve at three meshes, for each index."""
    return {
        key: [integral_scan(*key, points) for points in (256, 512, 1024)]
        for key in INDICES
    }


def test_the_groove_has_no_corner():
    """The case only exists below the threshold; with a corner the graded
    mesh applies and smoothing is not used."""
    assert not polyline_has_corners(GROOVE.t, GROOVE.y)
    assert len(GROOVE.t) == 2000


class TestEnergyBalance:
    @pytest.mark.parametrize("key", INDICES, ids=["3.1nm", "6.2nm", "80eV-like"])
    def test_closes_and_converges_monotonically(self, ladder, key):
        """R + A = 1 to the mesh's own convergence, falling at every rung.
        Measured 2026-10-29: 3.9e-4, 1.8e-5, 3.2e-6 at 3.1 nm; 1.9e-3,
        1.6e-4, 5.3e-6 at 6.2 nm; 1.8e-2, 1.3e-3, 7.5e-5 at the 80 eV-like
        index."""
        deviations = [
            abs(reflected(scan)[0] + scan.absorption[0] - 1.0) for scan in ladder[key]
        ]
        assert deviations[0] > deviations[1] > deviations[2]
        assert deviations[2] < 2e-4

    @pytest.mark.parametrize("key", INDICES, ids=["3.1nm", "6.2nm", "80eV-like"])
    def test_the_polygon_as_given_does_not_close(self, library, key):
        """The fault, kept visible: meshed as straight segments, the same
        groove misses unity by far more than the mesh's error. If this ever
        closes, the case has stopped exercising issue #6 and wants
        replacing, not deleting."""
        scan = integral_scan(*key, 512, smoothing=None)
        assert abs(reflected(scan)[0] + scan.absorption[0] - 1.0) > 5e-3
        assert scan.provenance.notes["profile"]["reason"] == "profile_smoothing=None"


class TestAgainstTheCMethod:
    @pytest.mark.parametrize("key", INDICES, ids=["3.1nm", "6.2nm", "80eV-like"])
    def test_reflectance_closes_on_the_c_method(self, ladder, cmethod_reflectance, key):
        """The smoothed geometry is a different curve from the polygon, and
        the C-method on the polygon is the reference, so the smoothing's
        own cost is a floor this cannot pass. It grows with the width:
        measured 2026-10-29 at 6.2 nm, 1.2e-5 at sigma = 1.46 nm (the width
        a scan reaching 3.1 nm chooses) and 2.6e-5 at 2.3 nm (6.2 nm on its
        own). Against that, the polygon as given misses the energy balance
        by 2.1e-2 at the same point."""
        reference = cmethod_reflectance[key]
        errors = [abs(reflected(scan)[0] - reference) for scan in ladder[key]]
        assert errors[1] < 2e-4
        assert errors[2] < 4e-5

    def test_the_absorption_is_the_c_methods(self, ladder, cmethod_reflectance):
        """The quantity issue #6 was about, at the index where it was
        worst: the C-method closes R + A to 1e-7, so 1 - R is its A."""
        key = (LOW_ENERGY, LOW_ENERGY_NM)
        absorbed = ladder[key][2].absorption[0]
        assert absorbed == pytest.approx(1.0 - cmethod_reflectance[key], abs=1e-4)


class TestProvenance:
    def test_auto_records_the_width_and_the_rule(self, ladder):
        note = ladder["Au", 6.199][0].provenance.notes["profile"]
        assert note["representation"] == "smoothed"
        assert note["rule"].startswith("auto")
        assert note["vertices"] == 2000
        # One sixteenth of gold's metal-side transverse wavelength at 6.2 nm.
        assert 2.0 < note["gaussian_sigma_nm"] < 2.6

    def test_a_given_width_is_used_and_recorded(self, library):
        scan = integral_scan("Au", 6.199, 256, smoothing=3.0)
        note = scan.provenance.notes["profile"]
        assert note["gaussian_sigma_nm"] == 3.0
        assert note["rule"] == "given"

    def test_an_analytic_profile_has_no_profile_note(self, library):
        problem = Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.1), coating="Au")
        scan = get_solver("integral").solve(
            problem, MOUNT, [6.199], boundary_points=256, conductivity="tabulated"
        )
        assert "profile" not in scan.provenance.notes


class TestNothingElseMoves:
    """Smoothing reaches only a measured polygon without a corner, under the
    tabulated condition. Everything else must come out bit for bit the same
    whatever ``profile_smoothing`` says."""

    @staticmethod
    def solve_all_ways(problem, conductivity):
        scans = [
            get_solver("integral").solve(
                problem,
                MOUNT,
                [6.199],
                boundary_points=256,
                conductivity=conductivity,
                profile_smoothing=smoothing,
            )
            for smoothing in ("auto", None, 2.0)
        ]
        for scan in scans[1:]:
            assert np.array_equal(scan.efficiency, scans[0].efficiency, equal_nan=True)
            if conductivity == "tabulated":
                assert np.array_equal(scan.absorption, scans[0].absorption)
        return scans

    @pytest.mark.parametrize(
        "profile",
        [Sinusoidal(depth_fraction=0.1), Blazed(blaze_angle=10.0, antiblaze_angle=30.0)],
        ids=["sinusoid", "blazed"],
    )
    def test_analytic_profiles(self, library, profile):
        self.solve_all_ways(Problem(period=PERIOD, profile=profile, coating="Au"), "tabulated")

    def test_a_measured_polygon_with_a_corner(self, library):
        """A real corner keeps the graded mesh, and the polygon is used as
        given; provenance says why."""
        cornered = FromProfileData(
            t=(0.0, 0.25, 0.30, 0.55, 1.0), y=(0.0, 0.0, 0.08, 0.08, 0.0)
        )
        assert polyline_has_corners(cornered.t, cornered.y)
        scans = self.solve_all_ways(
            Problem(period=PERIOD, profile=cornered, coating="Au"), "tabulated"
        )
        assert "graded" in scans[0].provenance.notes["profile"]["reason"]

    def test_the_perfect_conductor(self, library):
        scans = self.solve_all_ways(Problem(period=PERIOD, profile=GROOVE), "perfect")
        assert "perfectly conducting" in scans[0].provenance.notes["profile"]["reason"]


class TestRefusals:
    @pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), "smooth", True])
    def test_bad_smoothing_is_refused(self, value):
        problem = Problem(period=PERIOD, profile=GROOVE)
        with pytest.raises(ValueError, match="profile_smoothing"):
            get_solver("integral").solve(
                problem, MOUNT, [6.199], boundary_points=256, profile_smoothing=value
            )


class TestSmoothedPolyline:
    def test_a_flat_line_is_unchanged(self):
        x, y = smoothed_polyline([0.0, 0.5], [0.0, 0.0], 0.01, samples=64)
        assert np.allclose(y, 0.0, atol=1e-15)
        assert np.allclose(np.diff(x), 1.0 / 64)

    def test_the_mean_height_is_kept(self):
        """A Gaussian has unit gain at zero frequency, so the curve keeps
        its average height over the arc."""
        t, y = np.asarray(GROOVE.t), np.asarray(GROOVE.y)
        raw_x, raw_y = smoothed_polyline(t, y, 1e-9, samples=8192)
        x, smooth_y = smoothed_polyline(t, y, 0.01, samples=8192)
        assert smooth_y.mean() == pytest.approx(raw_y.mean(), abs=1e-12)
        assert np.abs(smooth_y - raw_y).max() > 1e-4  # and it did something

    def test_width_must_be_positive(self):
        with pytest.raises(ValueError, match="sigma"):
            smoothed_polyline([0.0, 0.5], [0.0, 0.1], 0.0, samples=64)
