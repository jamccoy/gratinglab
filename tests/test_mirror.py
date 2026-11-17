"""The flat mirror: Parratt's recursion, penetration depth and bulk thickness.

The recursion is checked against two things that share none of its code: the
single-interface Fresnel amplitudes in ``materials/fresnel.py`` for a bare
substrate, and a characteristic-matrix (transfer-matrix) calculation written
out here for a real stack. The rest are limits any stack must meet.
"""

import numpy as np
import pytest

from gratinglab import materials
from gratinglab.materials.fresnel import amplitude
from gratinglab.mirror import Mirror, stack_amplitude
from gratinglab.problem import ConformalLayer as Layer


def tmm_s(indices, thicknesses, graze, wavelength):
    """r_s by characteristic matrices (Born and Wolf 1.6), vacuum above."""
    k = 2 * np.pi / wavelength
    q = [np.sqrt(n**2 - np.cos(graze) ** 2 + 0j) for n in [1.0, *indices]]
    m = np.eye(2, dtype=complex)
    for qj, d in zip(q[1:-1], thicknesses):
        delta = k * qj * d
        m = m @ np.array([[np.cos(delta), -1j * np.sin(delta) / qj], [-1j * qj * np.sin(delta), np.cos(delta)]])
    q0, qs = q[0], q[-1]
    num = (m[0, 0] + m[0, 1] * qs) * q0 - (m[1, 0] + m[1, 1] * qs)
    den = (m[0, 0] + m[0, 1] * qs) * q0 + (m[1, 0] + m[1, 1] * qs)
    return num / den


class TestTheRecursion:
    def test_a_bare_substrate_is_fresnel(self):
        au = materials.typical("Au")
        graze = np.array([0.3, 1.0, 3.0, 10.0, 45.0])
        r_s, r_p = Mirror("typical:Au").amplitude(graze, 1.0)
        f_s, f_p = amplitude(au.n(1.0), np.radians(graze))
        assert np.allclose(r_s, f_s, atol=1e-12) and np.allclose(r_p, f_p, atol=1e-12)

    def test_bare_roughness_matches_fresnel_in_magnitude(self):
        # fresnel.py damps by the modulus of the Nevot-Croce exponential; the
        # stack keeps its phase, so the reflectivities agree, not the phases.
        au = materials.typical("Au")
        graze = np.array([0.5, 2.0, 6.0])
        r_s, _ = Mirror("typical:Au", roughness=0.4).amplitude(graze, 1.0)
        f_s, _ = amplitude(au.n(1.0), np.radians(graze), roughness_nm=0.4, wavelength_nm=1.0)
        assert np.allclose(np.abs(r_s), np.abs(f_s), rtol=1e-12)

    @pytest.mark.parametrize("graze_deg", [0.5, 2.0, 7.0, 30.0])
    def test_a_stack_agrees_with_characteristic_matrices(self, graze_deg):
        indices = [1 - 4e-3 + 2e-3j, 1 - 1e-3 + 1e-4j, 1 - 6e-3 + 4e-3j, 1 - 2e-3 + 3e-4j]
        thick = [3.0, 7.5, 2.2]
        graze = np.radians(graze_deg)
        r_s, _ = stack_amplitude(indices, thick, graze, 1.3)
        assert r_s == pytest.approx(tmm_s(indices, thick, graze, 1.3), abs=1e-12)

    def test_a_film_of_the_substrate_material_changes_nothing(self):
        bare = Mirror("typical:Au").reflectivity([0.5, 2, 5], 1.0)
        same = Mirror("typical:Au", (Layer(material="typical:Au", thickness=7.0),)).reflectivity([0.5, 2, 5], 1.0)
        assert np.allclose(bare, same, rtol=1e-12)

    def test_a_thick_absorbing_film_is_its_own_bulk(self):
        film = Mirror("typical:Si", (Layer(material="typical:Au", thickness=500.0),))
        assert np.allclose(film.reflectivity([1, 3, 6], 1.0), Mirror("typical:Au").reflectivity([1, 3, 6], 1.0), rtol=1e-9)

    def test_passive(self):
        m = Mirror("typical:Si", (Layer(material="typical:Ni", thickness=4.0), Layer(material="typical:C", thickness=2.0)))
        r = m.reflectivity(np.linspace(0.1, 89.9, 300)[:, None], np.array([0.5, 1.0, 4.0, 13.5]))
        assert ((r >= 0) & (r <= 1)).all()

    def test_layers_are_innermost_first(self):
        g = np.array([1.0, 4.0])
        ni_on_top = Mirror("typical:Si", (Layer(material="typical:C", thickness=3), Layer(material="typical:Ni", thickness=3)))
        indices = [materials.typical("Ni").n(1.0), materials.typical("C").n(1.0), materials.typical("Si").n(1.0)]
        r_s, _ = stack_amplitude(indices, [3, 3], np.radians(g), 1.0)
        assert np.allclose(ni_on_top.amplitude(g, 1.0)[0], r_s)

    def test_polarization_is_checked(self):
        with pytest.raises(ValueError, match="polarization"):
            Mirror("typical:Au").reflectivity(1.0, 1.0, polarization="TE")


@pytest.fixture
def lossless(tmp_path, monkeypatch):
    """A lossless X-ray-like material, n = 1 - 4e-3, as a plain n,k table."""
    monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
    (tmp_path / "Clear.nk").write_text("wavelength (nm), n, k\n0.5, 0.996, 0\n2.0, 0.996, 0\n", encoding="utf-8")
    return "Clear"


class TestDepth:
    def test_evanescent_below_the_critical_angle(self, lossless):
        graze = np.radians(2.0)  # critical angle is about 5.1 degrees
        expect = 1.0 / (2 * 2 * np.pi / 1.0 * np.sqrt(np.cos(graze) ** 2 - 0.996**2))
        got = Mirror(lossless).penetration_depth(2.0, 1.0)[lossless]
        assert got == pytest.approx(expect, rel=1e-9)

    def test_unbounded_above_it_without_absorption(self, lossless):
        assert np.isinf(Mirror(lossless).penetration_depth(8.0, 1.0)[lossless])

    def test_one_entry_per_material(self):
        m = Mirror("typical:Si", (Layer(material="typical:Au", thickness=2), Layer(material="typical:Si", thickness=2)))
        assert list(m.penetration_depth(1.0, 1.0)) == ["typical:Si", "typical:Au"]


class TestBulkThickness:
    def test_beyond_it_the_film_reads_as_bulk_and_just_short_it_does_not(self):
        m = Mirror("typical:Si", (Layer(material="typical:Au", thickness=1.0),))
        t = m.bulk_thickness(2.0, 1.0, tolerance=1e-3)
        bulk = float(Mirror("typical:Au").reflectivity(2.0, 1.0))

        def r(thickness):
            return float(Mirror("typical:Si", (Layer(material="typical:Au", thickness=thickness),)).reflectivity(2.0, 1.0))

        assert all(abs(r(x) - bulk) <= 1e-3 * bulk for x in t * np.array([1.0, 1.3, 2.0, 5.0]))
        assert abs(r(0.9 * t) - bulk) > 1e-3 * bulk * 0.5  # near the edge, not far inside it

    def test_tighter_tolerance_needs_more(self):
        m = Mirror("typical:Si", (Layer(material="typical:Au", thickness=1.0),))
        assert m.bulk_thickness(2.0, 1.0, 1e-4) > m.bulk_thickness(2.0, 1.0, 1e-2)

    def test_a_lossless_layer_above_its_critical_angle_never_settles(self, lossless):
        m = Mirror("typical:Au", (Layer(material=lossless, thickness=1.0),))
        assert np.isinf(m.bulk_thickness(8.0, 1.0))

    def test_a_bare_substrate_is_refused(self):
        with pytest.raises(ValueError, match="already bulk"):
            Mirror("typical:Au").bulk_thickness(2.0, 1.0)
