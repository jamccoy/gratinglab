"""The integral solver against the PCGrate reference corpus.

This is the project's first rigorous-vs-rigorous comparison: two independent
implementations of the integral method, one commercial and closed, one this
package, on the same measured groove.

The comparison also *identified the profile*. ``benchmarks/corpus.toml``
carried the TASTE run's ``.ggp`` as unconfirmed since the corpus was
recovered; sweeping the candidate profiles against the table found exactly
one that matches -- ``OGRE/AFM_test.ggp`` at azimuth +19.99 deg, our TE
against PCGrate's TE column, orderwise to 8e-4 across the full scan. The
ideal ``Blazed(29.5, 70.5)`` sawtooth disagrees by up to 0.29, so the run was
made from the measured AFM groove, not the design geometry.
"""

import numpy as np
import pytest

from gratinglab.illumination import Illumination
from gratinglab.io.efficiency_table import read_scan
from gratinglab.problem import Problem
from gratinglab.profiles import FromProfileData
from gratinglab.solvers import integral

from corpus import reference_dir

pytestmark = pytest.mark.skipif(
    reference_dir() is None,
    reason="PCGrate reference corpus not found; set GRATINGLAB_REF_DIR to enable",
)


def taste_profile() -> FromProfileData:
    """``OGRE/AFM_test.ggp``, parsed around its defect.

    The file ends with junk after the closing ``1 0`` pair (a lone number and
    a duplicated terminator), which ``io.ggp.read_ggp`` rightly refuses --
    a truncated pair mid-file would silently shift the profile. Here the
    damage is provably *after* the profile: pairs are taken in order and the
    non-monotone stray at the tail is dropped.
    """
    pairs = []
    with (reference_dir() / "OGRE" / "AFM_test.ggp").open() as handle:
        for line in handle:
            parts = line.split()
            if len(parts) == 2:
                try:
                    pairs.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
    t, y = np.array(pairs).T
    keep = np.concatenate(([True], np.diff(t) >= 0))
    return FromProfileData(
        t=tuple(np.clip(t[keep], 0.0, 1.0)), y=tuple(y[keep] - y[keep].min())
    )


class TestTasteWavescan:
    """Geometry from benchmarks/corpus.toml [taste]; the table is PCGrate's
    TE column (one polarization per exported table)."""

    def test_orderwise_agreement_on_the_identified_profile(self, ref_dir):
        path = ref_dir / "OGRE" / "tastetest_perf_wavescan.txt"
        if not path.exists():
            pytest.skip("TASTE wavescan not in corpus")
        reference = read_scan(path, method="integral")

        problem = Problem(period=315.15, profile=taste_profile())
        illumination = Illumination.offplane(
            graze=1.25, azimuth=19.99, polarization="TE"
        )
        wavelengths = reference.wavelengths[::40]
        scan = integral.solve(
            problem, illumination, wavelengths, boundary_points=256
        )

        worst = 0.0
        for w in wavelengths:
            ours, theirs = scan.at(float(w)), reference.at(float(w))
            for m in theirs.propagating_orders:
                worst = max(worst, abs(ours[int(m)] - theirs[int(m)]))
        # Measured 7.8e-4 over the full 541-point scan; PCGrate's own totals
        # sit at 1.0005, so a fair share of the gap is theirs.
        assert worst < 2e-3

    def test_the_design_sawtooth_is_not_what_pcgrate_ran(self, ref_dir):
        """Negative control for the identification: the idealised blaze
        disagrees by two orders of magnitude more, so the match above is
        information, not tolerance slack."""
        from gratinglab.profiles import Blazed

        path = ref_dir / "OGRE" / "tastetest_perf_wavescan.txt"
        if not path.exists():
            pytest.skip("TASTE wavescan not in corpus")
        reference = read_scan(path, method="integral")

        problem = Problem(
            period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
        )
        illumination = Illumination.offplane(
            graze=1.25, azimuth=19.99, polarization="TE"
        )
        wavelength = [2.6]
        scan = integral.solve(
            problem, illumination, wavelength, boundary_points=256
        )
        ours, theirs = scan.at(2.6), reference.at(2.6)
        worst = max(
            abs(ours[int(m)] - theirs[int(m)])
            for m in theirs.propagating_orders
        )
        assert worst > 0.05
