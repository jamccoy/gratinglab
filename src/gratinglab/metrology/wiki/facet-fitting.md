# Facet fitting and trim

Each detected groove is windowed, split into two facets, and one of them is fitted
with a straight line. The slope of that line is the blaze angle.

## Choosing the facet

`BLAZE_SIDE` selects which face to measure, **by slope sign rather than by
position**:

- `negative_slope` — the down-sloping facet (the default)
- `positive_slope` — the up-sloping facet
- `longer` — whichever facet is wider

Choosing by sign rather than by "left" or "right" makes the setting independent of
how the sample happened to be oriented under the tip. A scan taken with the
grating rotated 180° needs no change.

Anything else raises a `ValueError` rather than silently measuring the wrong face.
Older documentation for this project referred to `'left'` and `'right'`, which
were never valid values.

## The split

The window is divided at its **lowest point** — the groove trough. Everything
before that index is one facet, everything after is the other.

This is also the step that fails when a groove is clipped by a scan edge: if the
window is truncated, the minimum can land on the very first sample and one
"facet" becomes a few pixels wide. That is why edge grooves are rejected before
they reach this point. See **Row groups and scan edges**.

## Trimming

Fitting the full facet would include the rounded groove top and the flattened
trough, both of which curve away from the straight section and bias the angle
low. `FACET_TRIM` removes a fraction from each end before fitting.

The trim is **asymmetric**: 2.5× more is removed from the trough side than from
the land side, because the bottom of the groove is the more strongly rounded of
the two. The total removed is therefore

    FACET_TRIM × 3.5

### The hard limit

That multiplier has a consequence worth stating plainly. At
`FACET_TRIM = 1/3.5 ≈ 0.286`, the trim consumes the entire facet. Past that point
**every** groove fit fails and the analysis returns no measurements at all —
it does not degrade gracefully.

The setting is capped at **0.28** for this reason, in `settings.py` as
`MAX_FACET_TRIM`, and the GUI rejects anything above it with a message naming the
control rather than letting you discover it as an empty result.

### It moves the answer more than anything else

On the master sample:

| FACET_TRIM | Mean blaze angle |
|---|---|
| 0.05 | 32.94° |
| 0.10 | 33.23° |
| 0.15 | 32.84° |
| 0.20 | 31.58° |
| 0.25 | 29.99° |

That is a **3.2° swing** across a plausible range of a single parameter — larger
than most of the between-sample differences this software is used to detect.

It is not a flaw in the software so much as a fact about the measurement: how much
of a curved facet you call "the facet" changes the angle you get. What matters is
that the same trim is used across every sample being compared, and that the value
is reported alongside the result. Both are true — `FACET_TRIM` is written into
every summary file.

The GUI exists partly so this is easy to explore. Drag the trim spinner, re-run,
and watch the number move.

## Before any of this: flattening

The facet is fitted to the *flattened* profile, so the flattening method shapes
what is being measured. Across the four profile-flattening methods the mean blaze
angle moves by about 0.49° on a typical scan — comparable to facet trim, and
comparable to the differences between samples. See **Flattening**.

## Fit quality

Every groove fit records an R². It is **reported, not enforced** — no measurement
is dropped for a poor fit.

In practice the edge-exclusion rule removes the bad fits anyway: after it, no fit
in the sample dataset falls below R² 0.95. The worst R² is shown in the GUI
results panel and written to the per-groove CSV, so a degrading dataset would be
visible.
