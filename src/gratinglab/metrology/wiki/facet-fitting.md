# Facet fitting and trim

Each detected groove is windowed and split into two facets, and one of the two
is fitted with a straight line. The slope of that line is the blaze angle.

## Choosing the facet

`BLAZE_SIDE` selects which face is measured. The selection is made by the sign
of the slope and not by position:

- `negative_slope`: the down-sloping facet (the default)
- `positive_slope`: the up-sloping facet
- `longer`: whichever facet is wider

Selecting by sign makes the setting independent of how the sample was oriented
under the tip, where "left" or "right" would depend on it. A scan taken with the
grating rotated 180° requires no change.

Any other value raises a `ValueError`, so that the wrong face is never measured
without warning. Older documentation for this project referred to `'left'` and
`'right'`, which were never valid values.

## The split

The window is divided at its lowest point, which is the groove trough. All
samples before that index form one facet, and all samples after it form the
other.

This is also the step that fails when a groove is clipped by a scan edge. If the
window is truncated, the minimum can fall on the very first sample, and one
"facet" is then only a few pixels wide. For this reason edge grooves are
rejected before they reach this step, as described in "Row groups and scan
edges".

## Trimming

Fitting the full facet would include the rounded groove top and the flattened
trough, both of which curve away from the straight section and bias the angle
low. `FACET_TRIM` removes a fraction from each end of the facet before fitting.

The trim is asymmetric: 2.5× more is removed from the trough side than from the
land side, because the bottom of the groove is the more strongly rounded of the
two. The total fraction removed is therefore

    FACET_TRIM × 3.5

### The hard limit

That multiplier has a direct consequence. At `FACET_TRIM = 1/3.5 ≈ 0.286` the
trim consumes the entire facet. Beyond that value every groove fit fails and the
analysis returns no measurements at all, so the failure is abrupt.

For this reason the setting is capped at 0.28, in `settings.py` as
`MAX_FACET_TRIM`, and the GUI rejects any larger value with a message that names
the control, so the limit is reported directly and is not first encountered as
an empty result.

### It moves the answer more than anything else

On the master sample:

| FACET_TRIM | Mean blaze angle |
|---|---|
| 0.05 | 32.89° |
| 0.10 | 33.13° |
| 0.15 | 32.71° |
| 0.20 | 31.49° |
| 0.25 | 29.94° |

This is a 3.2° change across a plausible range of a single parameter, larger
than most of the between-sample differences that this software is used to
detect.

This is a property of the measurement more than a flaw in the software: the
portion of a curved facet that is taken to be "the facet" determines the angle
obtained. What matters is that the same trim is used for every sample being
compared and that the value is reported alongside the result. Both conditions
are met, since `FACET_TRIM` is written into every summary file.

The GUI is designed in part to make this easy to explore: the trim spinner can
be changed and the analysis re-run to observe the effect on the result.

## Before any of this: flattening

The facet is fitted to the flattened profile, so the flattening method affects
what is being measured. Across the four profile-flattening methods the mean
blaze angle varies by about 0.49° on a typical scan, which is comparable to the
effect of facet trim and to the differences between samples. This is described
in "Flattening".

## Fit quality

Every groove fit records an R². The R² is reported but is not used as a
criterion: no measurement is dropped for a poor fit.

In practice the edge-exclusion rule removes the poor fits, and after it is
applied no fit in the sample dataset falls below R² 0.95. The worst R² is shown
in the GUI results panel and written to the per-groove CSV, so a degrading
dataset would be visible.
