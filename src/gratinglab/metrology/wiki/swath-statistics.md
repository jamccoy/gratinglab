# Swath statistics

The boundary profile averages every row of the scan into one trace, and every
groove along that trace into one shape. This is the appropriate quantity to
export, but it conceals how much the groove varies. *Swath statistics* builds
the same averaged groove separately for each horizontal band of rows (the same
bands that the blaze-angle analysis calls row groups) and compares the results.

To use it in the metrology window, load a scan, open the Boundary tab and press
"Swath statistics…". The dialog uses the Boundary tab's settings and the tip
correction, so that each swath is built in exactly the same way as the exported
profile. The swath count is initially the row-group count and can be changed in
the dialog.

## The two spreads

| Number | Means |
|---|---|
| Shape spread between swaths | How much the averaged groove changes **down the image**, along the length of the grooves. RMS over the period, nm. |
| Groove-to-groove within one | How much individual grooves differ **side by side**, from their swath's average. |
| Whole scan vs swath mean | How far the exported whole-scan groove sits from the swaths' mean. |

A between-swath spread well above the within-swath spread indicates that the
shape changes along the grooves, which may be caused by a ruling or etch
gradient or by a scan that changed as it progressed. The reverse indicates that
the grooves differ from one another more than any one groove changes along its
length.

## Drift, and why the exported groove can be rounder than every swath

Each swath's groove is aligned to the others before they are compared, so a
groove of the same shape that sits slightly to one side is counted as the same
shape. The sideways position is reported separately as the *lateral shift*,
measured from the positions at which each swath finds its grooves.

A steady shift down the image means either that the grooves are not parallel to
the slow scan axis (the sample is slightly rotated) or that the scan drifted
thermally. This affects the export, because the whole-scan profile averages rows
whose grooves sit in different places and therefore blurs the shape. On the
synthetic grating, 20 nm of drift across the scan causes the whole-scan groove
depth to read about 5% low, while every swath reads correctly. The "Whole scan
vs swath mean" figure detects this case. If it is large and the shift plot shows
a slope, the image or scan should be rotated before the exported depth and
angles are relied upon.

## Not independent samples

Swaths re-measure the same grooves, so their spread describes how the groove
varies and does not represent the uncertainty of a mean. It should not be
divided by the square root of the swath count. The blaze-angle row groups have
the same property and are corrected for it, as described on the ICC correction
page.

## Export

"Export per-swath CSV…" writes one row per swath, containing its rows, lateral
shift, grooves used, period, groove-to-groove spread and the profile metrics,
with the three spreads above written as `#` comment lines under the header. The
same output is available from Python:

```python
from gratinglab.metrology.boundary import swath_statistics
stats = swath_statistics(data, scan_x_um, settings, n_swaths=8)
print(stats.summary())
stats.write_csv("swaths.csv")
```
