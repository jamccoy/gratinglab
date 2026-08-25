# Tip correction

An AFM image is not the surface. It is the surface **dilated** by the tip: at
every pixel the tip descends until something touches, and the height recorded is
where the *apex* stopped — which is above the true surface everywhere the contact
happened on the tip's flank instead of its point.

The correction (Villarrubia 1997) inverts that where the inversion exists, and
says so where it does not. It is **off by default**, because it changes measured
numbers — groove depth especially — and a corrected depth and an uncorrected
depth are different measurements that must not be filed under one name.

## The tip model

A cone with a spherical cap:

| Setting | Means |
|---|---|
| `TIP_RADIUS_NM` | The **apex radius**. A probe sold as "2 nm wide" is radius 1. |
| `TIP_HALF_ANGLE_DEG` | The cone half-angle, measured **from the tip axis**. 18° is typical. |
| `TIP_CORRECTION` | `'none'` (default) or `'erosion'`. |

The half-angle is the one people misread. Measured from the axis, an 18° tip has
a flank rising at **72°** from the surface — which clears a 70.5° anti-blaze
facet by 1.5°, and would not clear a 75° one at all. That near-parallelism, not
the nanometre apex, is what limits blazed-grating metrology.

Set them in `config.py`, or pass an `AnalysisSettings` with the fields set.
There is no control in the window yet.

## What you get back, and what you do not

Erosion returns the **least upper bound** on the true surface consistent with the
image and the tip. That is two statements:

- Wherever the **apex** made contact, the bound *is* the surface. Recovery there
  is exact — machine precision on clean data.
- Wherever it did not, the tip physically could not reach, and **no algorithm
  recovers what is there**. The trough wedge between two facets is the standard
  case: a sphere of radius *R* stands off a corner by a distance that grows with
  *R*, and that depth is simply not in the data.

So the number to read is not the corrected profile alone but the **certain
fraction** that comes with it — reported in the metrics, in the `.txt` sidecar,
and in the log line. The uncertain remainder means "at most this high", never
"this high". A corrected profile without that fraction would manufacture exactly
the false confidence the correction exists to remove.

On the test fixture — a known sawtooth imaged through a known 20 nm tip — the
uncorrected surface is wrong by 9.06 nm on average; after erosion, over the
certain pixels, by 0.18 nm. The residual is the noise floor, not the tip.

## What it found on a real scan

The reason this exists. A TASTE grating scan measures a blaze angle of 27.9°
from its facet fit and 20.3° from its depth — a 7.6° disagreement that on a
sharp groove is impossible, since depth and facet angle are locked together.
Tip rounding was the obvious suspect.

Running the correction **acquits the nominal tip**, in both directions:

| Apex radius | Depth/period after erosion | Implied blaze | Pixels certain |
|---|---|---|---|
| uncorrected | 0.3275 | 20.33° | — |
| 1 nm (nominal) | 0.3275 | 20.33° | 97.7% |
| 2 nm | 0.3275 | 20.33° | 95.9% |
| 5 nm | 0.3275 | 20.33° | 81.8% |
| 10 nm | 0.3274 | 20.32° | 74.8% |

Nothing moves. The measured surface is already reachable by a sharp tip, so
there is no tip-hidden depth to give back.

Forward, dilating an ideal 27.91°/70.5° groove and asking what a tip would
report:

| Apex radius | Image depth/period | Implied blaze |
|---|---|---|
| sharp truth | 0.4460 | 27.91° |
| 1 nm (nominal) | 0.4439 | 27.77° |
| 10 nm | 0.4316 | 27.00° |
| 40 nm | 0.3909 | 24.40° |
| **80 nm** | **0.3362** | **20.89°** |

A nominal tip costs 0.14° of the 7.6° gap. Reproducing the measurement takes an
apex worn to ~80 nm, forty to eighty times its specification.

**So the groove is probably rounded, or the tip was badly worn.** Which one is
not answerable from the scan: it needs the probe's condition at scan time, a
re-scan with a fresh tip, or a cross-section.

## When to turn it on

- **A tip you have reason to distrust** — late in its life, or after imaging
  something hard. Compare the corrected and uncorrected depth; if they differ,
  the tip was hiding surface.
- **Sharp features near the tip's own scale.** Trenches, undercuts, anything
  where the flank angle is close to the sidewall angle.

And when not to: as a matter of routine on well-behaved scans. If the erosion
returns the image unchanged, it has told you the tip was adequate — which is
worth knowing once, not on every run.
