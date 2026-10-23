# Tip correction

An AFM image is not the surface itself but the surface *dilated* by the tip. At
every pixel the tip descends until it makes contact, and the recorded height is
the position at which the apex stopped. Wherever the contact occurred on the
tip's flank and not at its point, this height lies above the true surface.

The correction (Villarrubia 1997) inverts this dilation where the inversion
exists, and reports where it does not. It is off by default, because it changes
measured numbers (groove depth in particular), and a corrected depth and an
uncorrected depth are different measurements that must not be recorded under the
same name.

## The tip model

The tip is modelled as a spherical apex cap on two straight flanks, one on
either side of the apex along the scan line, each with its own angle:

| Setting | Means |
|---|---|
| `TIP_RADIUS_NM` | The **apex radius**. Check the datasheet says radius, not width: Bruker quotes the radius. |
| `TIP_HALF_ANGLE_DEG` | The half angle of the flank on the **-x** side, measured **from the tip axis**. |
| `TIP_HALF_ANGLE_RIGHT_DEG` | The **+x** flank. `None` (the default) means the same as -x: a symmetric cone. |
| `TIP_CORRECTION` | `'none'` (default) or `'erosion'`. |

Datasheet angles are half angles, measured from the tip axis to one face, and
are not the full opening angle. A flank of half angle *a* rises at 90 - a
degrees from the surface, so an 18° tip has 72° flanks. Such a tip clears a
70.5° anti-blaze facet by 1.5° and would not clear a 75° facet at all. This
near-parallelism of flank and facet, more than the nanometre-scale apex, is what
limits the metrology of blazed gratings.

## Pyramids, and which face the scan sees

Real tips are pyramids, and in a profile across the grooves the tip meets the
surface with whichever of its faces point along the fast-scan axis. The
specifications of Bruker's ScanAsyst-Air, the probe used for the group's scans,
are as follows:

| | Bruker spec | Flank from the surface |
|---|---|---|
| Tip radius | 2 nm nominal, 12 nm maximum | |
| Front angle | 15 ± 2.5° | 75° |
| Back angle | 25 ± 2.5° | **65°** |
| Side angle | 17.5 ± 2.5° | 72.5° |

The cross-section required by the correction therefore depends on the
orientation of the cantilever relative to the scan:

| Scanning | -x flank | +x flank |
|---|---|---|
| across the cantilever | 17.5° (side) | 17.5° (side) |
| along it, front face to -x | 15° (front) | 25° (back) |
| along it, front face to +x | 25° (back) | 15° (front) |

The orientation matters more than the radius. A wall rising towards +x is
followed by the tip's +x flank. When the scan runs along the cantilever with the
back face on that side, the flank rises at only 65°, and a 70° anti-blaze facet
is touched at its top corner and nowhere below it, so the facet cannot be
recovered by any algorithm. If the sample or the scan is rotated by 180°, the
75° front face follows the wall instead, and erosion recovers nearly all of the
facet. In the test suite, the same tip, the same 70° wall and the same sampling
give at most one certain pixel per wall in one orientation and eight or nine of
ten in the other.

Across the cantilever both flanks are at 72.5°, which clears 70° by only 2.5°,
equal to the ±2.5° tolerance on the angle itself. A recovered anti-blaze facet
on such a scan should be regarded as marginal.

The orientation in which a scan was taken is determined by the instrument and
the scan angle, and the image does not record it reliably, so the choice is left
to the user. If the orientation is not known, both along-cantilever presets
should be run: the regions where they disagree are exactly the parts of the
groove that the orientation determines.

## Setting it

In the metrology window, under "Analysis parameters", choose `erosion` under Tip
correction and then select a tip shape. "Symmetric cone" ties the two flanks
together. The three ScanAsyst-Air presets fill in the radius and both angles,
and editing any number afterwards switches the shape to "Custom", so the picker
never names a tip that the numbers do not describe. The Boundary profile is
corrected immediately, and the fraction of pixels recovered exactly is shown
under its metrics. The blaze angles are corrected at the next Run, and the
results state which tip was assumed, including both flanks when they differ.
Outside the window, the fields can be set in `config.py` or passed in an
`AnalysisSettings`. `gratinglab.metrology.core.tip.TIP_PRESETS` holds the
presets by name.

A new tip is close to its nominal radius and a worn one is not, so a correction
made with the nominal 2 nm radius is the optimistic case. For a tip that has seen
heavy use, the 12 nm maximum should be tried as well.

The exported groove shape is expected to change less than the image does. The
export sets its lowest point to zero, and after erosion that point is the trough
floor that the tip never reached, which remains an upper bound. On the dilated
test fixture the image error falls from 9.06 nm to 0.87 nm, while the averaged
groove's shape error, with that offset removed, remains about 2 nm, because the
tip's 72 degree flank already followed the 70 degree facet.

## What you get back, and what you do not

Erosion returns the *least upper bound* on the true surface that is consistent
with the image and the tip. This result has two parts:

- Wherever the apex made contact, the bound is the surface. Recovery there is
  exact, to machine precision on clean data.
- Wherever the apex did not make contact, the tip could not physically reach the
  surface, and no algorithm can recover what is there. The trough wedge between
  two facets is the standard case: a sphere of radius *R* stands off a corner by
  a distance that grows with *R*, and that depth is not contained in the data.

The quantity to read is therefore the corrected profile together with the
*certain fraction* that accompanies it, which is reported in the metrics, in the
`.txt` sidecar, and in the log line. The uncertain remainder means "at most this
high" and never "this high". A corrected profile reported without that fraction
would create the false confidence that the correction is intended to remove.

On the test fixture (a known sawtooth imaged through a known 20 nm tip), the
uncorrected surface is in error by 9.06 nm on average, and after erosion, over
the certain pixels, by 0.18 nm. This residual is the noise floor and is not due
to the tip.

## What it found on a real scan

This scan is the reason the correction was implemented. A TASTE grating scan
gives a blaze angle of 27.9° from its facet fit and 20.3° from its depth. On a
sharp groove a 7.6° disagreement of this kind is impossible, since depth and
facet angle are geometrically linked, and tip rounding was the obvious suspect.

Running the correction shows that the nominal tip is not responsible, in both
directions of the analysis:

| Apex radius | Depth/period after erosion | Implied blaze | Pixels certain |
|---|---|---|---|
| uncorrected | 0.3275 | 20.33° | n/a |
| 1 nm | 0.3275 | 20.33° | 97.7% |
| 2 nm (ScanAsyst-Air nominal) | 0.3275 | 20.33° | 95.9% |
| 5 nm | 0.3275 | 20.33° | 81.8% |
| 10 nm | 0.3274 | 20.32° | 74.8% |

The results do not change. The measured surface is already reachable by a sharp
tip, so there is no depth hidden by the tip to recover.

In the forward direction, an ideal 27.91°/70.5° groove was dilated to determine
what a given tip would report:

| Apex radius | Image depth/period | Implied blaze |
|---|---|---|
| sharp truth | 0.4460 | 27.91° |
| 1 nm | 0.4439 | 27.77° |
| 10 nm | 0.4316 | 27.00° |
| 40 nm | 0.3909 | 24.40° |
| **80 nm** | **0.3362** | **20.89°** |

A sharp tip accounts for 0.14° of the 7.6° discrepancy. Reproducing the
measurement requires an apex worn to ~80 nm, which is forty times the
ScanAsyst-Air's 2 nm nominal radius and nearly seven times Bruker's 12 nm
maximum.

Both tables were produced with a symmetric 18° cone, close to the
ScanAsyst-Air's 17.5° sides. They predate the asymmetric model, and so they do
not yet address a scan taken along the cantilever, where the 65° back face would
hide more of the anti-blaze facet than any radius in the tables.

The results suggest that the groove is rounded, or that the tip was badly worn.
The scan alone cannot distinguish the two: this requires the condition of the
probe at the time of the scan, a re-scan with a fresh tip, or a cross-section.

## When to turn it on

- **A tip you have reason to distrust**: late in its life, or after imaging
  something hard. If the corrected and uncorrected depths differ, the tip was
  hiding part of the surface.
- **Sharp features near the tip's own scale.** Trenches, undercuts, and any
  feature where the flank angle is close to the sidewall angle.

The correction is not needed as a matter of routine on well-behaved scans. If
erosion returns the image unchanged, this shows that the tip was adequate, which
is useful to establish once but not on every run.
