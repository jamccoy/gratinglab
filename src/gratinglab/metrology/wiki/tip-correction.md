# Tip correction

An AFM image is not the surface. It is the surface **dilated** by the tip: at
every pixel the tip descends until something touches, and the height recorded is
where the *apex* stopped, which is above the true surface everywhere the contact
happened on the tip's flank instead of its point.

The correction (Villarrubia 1997) inverts that where the inversion exists, and
says so where it does not. It is **off by default**, because it changes measured
numbers (groove depth especially), and a corrected depth and an uncorrected
depth are different measurements that must not be filed under one name.

## The tip model

A spherical apex cap on two straight flanks, one either side of the apex along
the scan line, each with its own angle:

| Setting | Means |
|---|---|
| `TIP_RADIUS_NM` | The **apex radius**. Check the datasheet says radius, not width: Bruker quotes the radius. |
| `TIP_HALF_ANGLE_DEG` | The half angle of the flank on the **-x** side, measured **from the tip axis**. |
| `TIP_HALF_ANGLE_RIGHT_DEG` | The **+x** flank. `None` (the default) means the same as -x: a symmetric cone. |
| `TIP_CORRECTION` | `'none'` (default) or `'erosion'`. |

Datasheet angles are **half angles**, measured from the tip axis to one face,
not the full opening angle. A flank of half angle *a* rises at **90 - a**
degrees from the surface, so an 18° tip has 72° flanks. That clears a 70.5°
anti-blaze facet by 1.5° and would not clear a 75° one at all. That
near-parallelism, not the nanometre apex, is what limits blazed-grating
metrology.

## Pyramids, and which face the scan sees

Real tips are pyramids, and in a profile across the grooves the tip meets the
surface with whichever faces point along the fast-scan axis. Bruker's
ScanAsyst-Air, the probe used for the group's scans:

| | Bruker spec | Flank from the surface |
|---|---|---|
| Tip radius | 2 nm nominal, 12 nm maximum | |
| Front angle | 15 ± 2.5° | 75° |
| Back angle | 25 ± 2.5° | **65°** |
| Side angle | 17.5 ± 2.5° | 72.5° |

So the cross-section the correction needs depends on how the cantilever lay
against the scan:

| Scanning | -x flank | +x flank |
|---|---|---|
| across the cantilever | 17.5° (side) | 17.5° (side) |
| along it, front face to -x | 15° (front) | 25° (back) |
| along it, front face to +x | 25° (back) | 15° (front) |

**The orientation matters more than the radius.** A wall rising towards +x is
followed by the tip's +x flank. Scanned along the cantilever with the back face
there, that flank rises at only 65°, and a 70° anti-blaze facet is touched at
its top corner and nowhere below it: the facet is unrecoverable by any
algorithm. Turn the sample or the scan 180° and the 75° front face follows it
instead, and erosion recovers nearly all of it. In the test suite the same
tip, the same 70° wall and the same sampling give at most one certain pixel per
wall one way round and eight or nine of ten the other.

Across the cantilever both flanks are 72.5°, which clears 70° by only 2.5°, the
same as the ±2.5° tolerance on the angle itself. Treat a recovered anti-blaze
facet on such a scan as marginal.

Which way round a scan was taken is a fact about the instrument and the scan
angle, not something the image records reliably, so the choice is yours. If you
do not know it, run both along-cantilever presets: where they disagree is
exactly the part of the groove the orientation decides.

## Setting it

In the metrology window, under **Analysis parameters**: choose `erosion` under
Tip correction, then a **Tip shape**. *Symmetric cone* ties the two flanks
together. The three ScanAsyst-Air presets fill in the radius and both angles,
and editing any number afterwards switches the shape to *Custom*, so the picker
never names a tip the numbers do not describe. The Boundary profile is
corrected straight away, with the fraction of pixels recovered exactly shown
under its metrics; the blaze angles are corrected at the next Run, and the
results say which tip was assumed, both flanks included when they differ.
Outside the window, set the fields in `config.py` or pass an `AnalysisSettings`
with them set. `gratinglab.metrology.core.tip.TIP_PRESETS` holds the presets by
name.

A new tip is near its nominal radius and a worn one is not, so a correction made
with the nominal 2 nm is the optimistic one. If a tip has had a hard life, try
the 12 nm maximum as well.

Expect the exported groove shape to change less than the image does. The
export sets its lowest point to zero, and after erosion that point is the
trough floor the tip never reached, which stays an upper bound. On the
dilated test fixture the image error falls from 9.06 nm to 0.87 nm while the
averaged groove's shape error, with that offset removed, stays about 2 nm,
because the tip's 72 degree flank already followed the 70 degree facet.

## What you get back, and what you do not

Erosion returns the **least upper bound** on the true surface consistent with the
image and the tip. That is two statements:

- Wherever the **apex** made contact, the bound *is* the surface. Recovery there
  is exact, to machine precision on clean data.
- Wherever it did not, the tip physically could not reach, and **no algorithm
  recovers what is there**. The trough wedge between two facets is the standard
  case: a sphere of radius *R* stands off a corner by a distance that grows with
  *R*, and that depth is simply not in the data.

So the number to read is not the corrected profile alone but the **certain
fraction** that comes with it, reported in the metrics, in the `.txt` sidecar,
and in the log line. The uncertain remainder means "at most this high", never
"this high". A corrected profile without that fraction would manufacture exactly
the false confidence the correction exists to remove.

On the test fixture (a known sawtooth imaged through a known 20 nm tip), the
uncorrected surface is wrong by 9.06 nm on average; after erosion, over the
certain pixels, by 0.18 nm. The residual is the noise floor, not the tip.

## What it found on a real scan

The reason this exists. A TASTE grating scan measures a blaze angle of 27.9°
from its facet fit and 20.3° from its depth, a 7.6° disagreement that on a
sharp groove is impossible, since depth and facet angle are locked together.
Tip rounding was the obvious suspect.

Running the correction **acquits the nominal tip**, in both directions:

| Apex radius | Depth/period after erosion | Implied blaze | Pixels certain |
|---|---|---|---|
| uncorrected | 0.3275 | 20.33° | — |
| 1 nm | 0.3275 | 20.33° | 97.7% |
| 2 nm (ScanAsyst-Air nominal) | 0.3275 | 20.33° | 95.9% |
| 5 nm | 0.3275 | 20.33° | 81.8% |
| 10 nm | 0.3274 | 20.32° | 74.8% |

Nothing moves. The measured surface is already reachable by a sharp tip, so
there is no tip-hidden depth to give back.

Forward, dilating an ideal 27.91°/70.5° groove and asking what a tip would
report:

| Apex radius | Image depth/period | Implied blaze |
|---|---|---|
| sharp truth | 0.4460 | 27.91° |
| 1 nm | 0.4439 | 27.77° |
| 10 nm | 0.4316 | 27.00° |
| 40 nm | 0.3909 | 24.40° |
| **80 nm** | **0.3362** | **20.89°** |

A sharp tip costs 0.14° of the 7.6° gap. Reproducing the measurement takes an
apex worn to ~80 nm: forty times the ScanAsyst-Air's 2 nm nominal radius, and
nearly seven times Bruker's 12 nm maximum.

Both tables were made with a symmetric 18° cone, close to the ScanAsyst-Air's
17.5° sides. They predate the asymmetric model, so they say nothing yet about a
scan taken along the cantilever, where the 65° back face would hide more of the
anti-blaze facet than any radius in them.

So the groove is probably rounded, or the tip was badly worn. Which one is
not answerable from the scan: it needs the probe's condition at scan time, a
re-scan with a fresh tip, or a cross-section.

## When to turn it on

- **A tip you have reason to distrust**: late in its life, or after imaging
  something hard. Compare the corrected and uncorrected depth; if they differ,
  the tip was hiding surface.
- **Sharp features near the tip's own scale.** Trenches, undercuts, anything
  where the flank angle is close to the sidewall angle.

And when not to: as a matter of routine on well-behaved scans. If the erosion
returns the image unchanged, it has told you the tip was adequate, which is
worth knowing once, not on every run.
