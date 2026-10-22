# Swath statistics

The boundary profile averages every row of the scan into one trace, and every
groove along that trace into one shape. That is the right thing to export, and
it hides how much the groove varies. **Swath statistics** builds the same
averaged groove separately for each horizontal band of rows (the same bands the
blaze-angle analysis calls row groups) and compares them.

In the metrology window: load a scan, open the **Boundary** tab and press
**Swath statistics…**. The dialog uses the Boundary tab's settings and the tip
correction, so each swath is built exactly as the exported profile is. The
swath count starts at the row-group count and can be changed there.

## The two spreads

| Number | Means |
|---|---|
| Shape spread between swaths | How much the averaged groove changes **down the image**, along the length of the grooves. RMS over the period, nm. |
| Groove-to-groove within one | How much individual grooves differ **side by side**, from their swath's average. |
| Whole scan vs swath mean | How far the exported whole-scan groove sits from the swaths' mean. |

A between-swath spread well above the within-swath one says the shape changes
along the grooves: a ruling or etch gradient, or a scan that changed as it went.
The reverse says the grooves differ from one another more than any one groove
changes along its length.

## Drift, and why the exported groove can be rounder than every swath

Each swath's groove is aligned to the others before they are compared, so a
groove that has the same shape but sits a little to one side counts as the same
shape. The sideways position is reported separately as **lateral shift**,
measured from where each swath finds its grooves.

A steady shift down the image means the grooves are not parallel to the slow
scan axis (the sample is rotated slightly) or the scan drifted thermally. It
matters for the export: the whole-scan profile averages rows whose grooves sit
in different places, so it **blurs** the shape. On the synthetic grating, 20 nm
of drift across the scan makes the whole-scan groove depth read about 5% low,
while every swath reads correctly. That is what **Whole scan vs swath mean**
catches. If it is large and the shift plot is a slope, rotate the image or scan
before trusting the exported depth and angles.

## Not independent samples

Swaths re-measure the same grooves, so their spread describes how the groove
varies, not the uncertainty of a mean. Do not divide it by the square root of
the swath count. The blaze-angle row groups have the same problem and correct
for it; see the ICC correction page.

## Export

**Export per-swath CSV…** writes one row per swath: its rows, lateral shift,
grooves used, period, groove-to-groove spread and the profile metrics, with the
three spreads above as `#` comment lines under the header. From Python:

```python
from gratinglab.metrology.boundary import swath_statistics
stats = swath_statistics(data, scan_x_um, settings, n_swaths=8)
print(stats.summary())
stats.write_csv("swaths.csv")
```
