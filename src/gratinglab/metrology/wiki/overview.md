# Overview

This software measures *blaze angles* on diffraction gratings from AFM
topography scans, and it exports the averaged groove shape as a *boundary
profile* for efficiency modelling.

Both outputs are produced by a single pipeline that divides only at its last
step.

## The pipeline

```
AFM scan file (.txt)
      ↓  load, detect the scan width from the header
      ↓  optionally mask out a region
      ↓  level the image, and optionally repair scan-line scars
      ↓  correct the groove orientation when it matters
      ↓  optionally undo the tip's dilation
      ↓  average rows into a height profile
      ↓  flatten (remove the background tilt or curvature)
      ↓  detect groove centres, reject those clipped by a scan edge
      ↓  extract a window around each groove
      ├──────────────────────────┬───────────────────────────
      ↓                          ↓
  fit each facet             average the grooves
  → blaze angles             → one normalised groove
  → statistics, CSV          → .json, .csv or .ggp
```

Every step above the fork is shared code, and this is deliberate. The two
outputs were once produced by separate scripts, each with its own copy of the
loading and detection logic, and the copies drifted apart until a bug fixed in
one remained in the other.

## What a "measurement" is

A single scan is divided into *row groups*, which are horizontal bands (20 by
default). Each band is averaged into its own profile and analysed separately, so
that one image yields 80–100 blaze-angle measurements, where averaging the whole
image into a single profile would yield 4–10.

The statistical treatment of these measurements is where most of the care in
this project has gone. The measurements are not independent of one another,
because every band re-measures the same physical grooves; this is discussed in
"The correlation correction".

## Where to go next

- **Flattening**: the two stages, which of them changes the answer, and by how
  much
- **Scars and masking**: repairing scan lines the tip got wrong, and leaving a
  particle or a bad patch out of the analysis
- **Tip correction**: undoing the tip's effect on the measured surface where
  that is possible, and stating where it is not
- **Row groups and scan edges**: what row-group analysis provides, and why
  grooves at the edge of a scan are discarded
- **The correlation correction**: why the error bars are wider than √N suggests
- **Facet fitting and trim**: how the blaze facet is chosen and measured, and
  the parameter that has the largest effect on the answer
- **Reading the outputs**: the meaning of every column and file
