"""
Generate the synthetic AFM scan the metrology tests run on.

Real scans are the research group's data and are not committed (see
`config.DEFAULT_AFM_DIR`). Without a substitute, most of the metrology suite
would skip on any machine that lacks them -- including CI, which would then be
green while testing almost none of the pipeline.

So this writes a *known* grating: an ideal sawtooth of stated period and blaze
angle, in the Gwyddion text export format `core.processing._load_text` already
parses. Two properties matter more than realism:

1. **It is sharp.** Depth and blaze angle are locked together by the geometry,
   which is exactly what a real tip-rounded scan is not -- see `findings.md`,
   "The measured groove is rounded, not faceted". A test can therefore assert
   that the pipeline *recovers* the angle it was given, which no real scan
   allows, because on a real scan the true answer is unknown.
2. **It is deterministic.** Seeded noise, so a failure is a change in the code
   and never a change in the fixture.

Regenerate with:

    python tools/metrology/make_synthetic_scan.py

and commit the result. It is checked in rather than built at test time so that
a fixture change is a reviewable diff.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
TARGET = REPO / "tests" / "metrology" / "fixtures" / "synthetic_blazed_scan.txt"
DILATED_TARGET = (REPO / "tests" / "metrology" / "fixtures"
                  / "synthetic_dilated_scan.txt")

# A 2 um scan at 315 nm gives ~6.3 grooves: enough for the period to be measured
# from spacing, and for edge exclusion to reject one at each end and still leave
# grooves to average.
SCAN_X_UM = 2.0
PERIOD_NM = 315.0
BLAZE_DEG = 30.0
ANTIBLAZE_DEG = 70.0
N_ROWS = 64
N_COLS = 384              # ~61 samples per groove; the real TASTE scan has ~81
NOISE_NM = 0.15          # well under the facet signal; enough to exercise smoothing
TILT_NM_PER_UM = 4.0     # a plausible sample tilt for flattening to remove
SEED = 20260821

# The tip that images the dilated companion fixture. Deliberately *blunter*
# than the nominal sharp probe (R ~ 1 nm): at this fixture's ~5.2 nm pixel
# pitch a 1 nm apex is sub-pixel and dilation would barely move a sample,
# leaving nothing for a recovery test to recover. A 20 nm radius is a worn
# tip, and its rounding spans several pixels -- resolvable, so the test has
# teeth. The 18 deg half angle matches the real probes (flank at 72 deg,
# steeper than both facets below, so the surface is recoverable except the
# trough wedge the sphere cannot enter).
TIP_RADIUS_NM = 20.0
TIP_HALF_ANGLE_DEG = 18.0


def sawtooth(t: np.ndarray, blaze_deg: float, antiblaze_deg: float) -> np.ndarray:
    """One period of an ideal blazed groove, in units of the period.

    Same geometry as `gratinglab.profiles.Blazed`: the apex sits where the two
    facets meet, and depth = 1 / (cot(blaze) + cot(antiblaze)).
    """
    cot_b = 1.0 / np.tan(np.radians(blaze_deg))
    cot_a = 1.0 / np.tan(np.radians(antiblaze_deg))
    apex = cot_b / (cot_b + cot_a)
    depth = 1.0 / (cot_b + cot_a)

    # Mirrored so the *long* facet descends. Real scans of these gratings run
    # that way -- 64% of the TASTE trace has negative slope -- and
    # `config.BLAZE_SIDE = 'negative_slope'` selects the blaze facet by that
    # sign. A fixture with the opposite handedness sends the facet fit onto the
    # anti-blaze and reports a meaningless angle, which is what the first
    # version of this file did.
    t = np.mod(-t, 1.0)
    rising = t / apex * depth
    falling = (1.0 - t) / (1.0 - apex) * depth
    return np.where(t <= apex, rising, falling)


def build(dilated: bool = False) -> np.ndarray:
    """The scan, as heights in metres -- the unit the loader promises.

    ``dilated=True`` images the same sawtooth through the tip above --
    dilation first, instrument noise after, which is the order the physics
    has them in. The RNG is re-seeded per call and drawn in the same order,
    so the two fixtures carry *identical* noise and tilt: any difference
    between them is the tip, nothing else.
    """
    rng = np.random.default_rng(SEED)

    x_nm = np.linspace(0.0, SCAN_X_UM * 1000.0, N_COLS, endpoint=False)
    groove_nm = sawtooth(x_nm / PERIOD_NM, BLAZE_DEG, ANTIBLAZE_DEG) * PERIOD_NM

    rows = np.tile(groove_nm, (N_ROWS, 1))
    if dilated:
        from gratinglab.metrology.core.tip import dilate
        rows = dilate(rows * 1e-9, SCAN_X_UM,
                      radius_nm=TIP_RADIUS_NM,
                      half_angle_deg=TIP_HALF_ANGLE_DEG) * 1e9
    # A per-row offset, so `align_rows` flattening has something to correct and
    # a test asserting it is a no-op on the final angle is not vacuous.
    rows += rng.normal(0.0, 0.8, size=(N_ROWS, 1))
    rows += (x_nm / 1000.0) * TILT_NM_PER_UM
    rows += rng.normal(0.0, NOISE_NM, size=rows.shape)

    return rows * 1e-9


def write(path: Path, data: np.ndarray) -> Path:
    """Four header lines then tab-separated metres, as Gwyddion exports.

    The loader reads the width out of these comments with a regex and falls back
    to its `default_scan_size` if it cannot; writing them keeps the fixture a
    test of that path too.
    """
    height_um = SCAN_X_UM * data.shape[0] / data.shape[1]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("# Channel: ZSensor\n")
        handle.write(f"# Width: {SCAN_X_UM:.3f} µm\n")
        handle.write(f"# Height: {height_um:.3f} µm\n")
        handle.write("# Value units: m\n")
        for row in data:
            handle.write("\t".join(f"{v:.4e}" for v in row) + "\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=TARGET)
    parser.add_argument("--dilated-output", type=Path, default=DILATED_TARGET)
    args = parser.parse_args()
    path = write(args.output, build())
    print(f"wrote {path} ({path.stat().st_size / 1024:.0f} KB, "
          f"{N_ROWS}x{N_COLS}, period {PERIOD_NM:g} nm, blaze {BLAZE_DEG:g} deg)")
    path = write(args.dilated_output, build(dilated=True))
    print(f"wrote {path} (same grating imaged through a "
          f"R = {TIP_RADIUS_NM:g} nm, {TIP_HALF_ANGLE_DEG:g} deg tip)")


if __name__ == "__main__":
    main()
