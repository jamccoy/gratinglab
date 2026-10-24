# Groove orientation

A sample is easily mounted a degree or two off square, and then the grooves
are not perpendicular to the fast scan axis. Every scan line crosses them
obliquely. Nothing in the image looks wrong, and two things are.

| | Effect at angle θ | At 1° | At 3° |
|---|---|---|---|
| **Period** | reads long by 1/cos θ | +0.015% | +0.14% (315 → 315.43 nm) |
| **Blaze angle** | tan reads short by cos θ | 30° reads 29.996° | 30° reads 29.966° |
| **Averaged groove** | smeared across rows | 9% shallower | 32% shallower |

The last row is the one that matters. The grooves sit at a different x in every
row, drifting by tan θ per unit of slow-axis travel, so averaging rows averages
a groove with displaced copies of itself. A row group is short and loses
little, but the boundary profile averages the whole image: on a 2 µm tall scan
at 1° the grooves drift 35 nm from the first row to the last, 11% of a 315 nm
period. The depths in the table are the whole-image average on the synthetic
grating used by the test suite.

## How the angle is measured

The pipeline already finds grooves; this uses them. After image flattening the
scan is cut into up to 32 horizontal bands of at least 3 rows. In each band
the grooves are found with the same detection the analysis uses, and each
position is refined to a fraction of a pixel by a parabola through the minimum.
Positions are linked band to band into one track per groove, and a single slope
dx/dy is fitted to every track at once, each track keeping its own intercept.
Then θ = atan(dx/dy), and its standard error comes from the scatter of the
positions about the fit. One pass of outlier rejection removes a groove that
was linked to its neighbour.

Sign: θ is positive when the grooves move toward larger x as the row index
increases.

dx/dy is a ratio of lengths, so it needs the size of a pixel along the slow
axis as well as the fast one. That comes from the file:

| Format | Slow-axis size from |
|---|---|
| `.spm` | the second `Scan Size` value, times recorded lines over `Lines` |
| `.ibw` | the Y scaling; or `ScanSize` × `SlowRatio` over `ScanLines` in the note |
| text | a `Height:` line in the header, as Gwyddion writes it |

The `.spm` scaling is there because a scan stopped early still carries the size
it was set to cover. The TASTE scan holds 401 lines of 512 samples while its
header still reads 2 µm by 2 µm.

**When the file does not say, the angle is not measured.** A square scan is not
assumed: on a scan twice as wide as it is tall, that assumption reports about
half the real angle, and nothing would show that it had been assumed. The
report says "not measured" and gives this as the reason.

## How it is corrected

Each row is shifted along x so the grooves line up from row to row, and the
fast-axis pixel is then taken as dx · cos θ. For a grating, whose surface does
not change along the grooves, that is exact: a row of the rotated surface is
f(cos θ · (x + y tan θ)), the shift removes the y term and the pixel change
removes the cos θ. Heights are not touched, so the facet slopes come out right
on their own.

This was chosen over rotating the image. A rotation interpolates in two
dimensions, crops four corners from both axes, and needs square pixels, so a
scan with unequal pixels would first be resampled. The shift interpolates
along the row only and loses only the columns it slides off the ends, at most
tan θ of the scan height, split between the two sides. Cropped columns are
reported.

The interpolation is linear and has a measured cost. On an ideal knife-edge
sawtooth it blunts the apex and trough by about one pixel over one period: the
corrected groove reads 3% shallow at 33 pixels per period and under 0.5% at
120. A real groove is rounded by the tip over several nanometres and has no
sub-pixel corner to blunt.

It runs after image flattening and before tip correction. An axisymmetric tip
on a surface that is invariant along the grooves acts as a one-dimensional
dilation along the groove normal, which is the profile the correction
produces, not the oblique one the scan recorded.

## When it is applied

`GROOVE_ORIENTATION` in `config.py`, or the Import tab's *Groove orientation*:

| Mode | Measures | Corrects |
|---|---|---|
| `auto` (default) | yes | when \|θ\| ≥ 0.2° **and** \|θ\| ≥ 3 standard errors |
| `always` | yes | always; refuses if θ cannot be measured |
| `detect` | yes | never |
| `off` | no | never |

The 0.2° threshold is `GROOVE_ORIENTATION_THRESHOLD_DEG`. Below it the
whole-image drift on a square scan is under 0.4% of the scan width, smaller
than the tip's own effect on the groove, and the interpolation would cost more
than it returns. The significance condition stops a noisy scan from being
corrected for a rotation it cannot tell from zero.

## Where the angle is reported

Whether or not it was applied:

- the analysis result carries `groove_orientation`, with the angle, its
  standard error, how many grooves and positions it came from, whether it was
  applied and why or why not;
- the results panel prints a *Groove rotation* line;
- the Import tab shows it under the mode, and the Boundary tab's metrics;
- a boundary profile's metrics carry `groove_rotation_deg`,
  `groove_rotation_stderr_deg` and `groove_rotation_corrected`, so the `.json`
  export and the metrics sidecar both say whether its period is along the scan
  line or along the groove normal.

The swath statistics measure the same thing from the other side: a steady
lateral shift down the image is this rotation (or thermal drift, which a
linear fit cannot tell apart from it). See the `swath-statistics` page.
