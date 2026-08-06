# Overview

This software measures **blaze angles** on diffraction gratings from AFM
topography scans, and exports the averaged groove shape as a **PCGrate boundary
profile** for efficiency modelling.

Both outputs come from one pipeline that forks only at the last step.

## The pipeline

```
AFM scan file (.txt)
      ↓  load, detect the scan width from the header
      ↓  average rows into a height profile
      ↓  flatten (remove the background tilt or curvature)
      ↓  detect groove centres, reject those clipped by a scan edge
      ↓  extract a window around each groove
      ├──────────────────────────┬───────────────────────────
      ↓                          ↓
  fit each facet             average the grooves
  → blaze angles             → one normalised groove
  → statistics, CSV          → .ggp for PCGrate
```

Everything above the fork is shared code. That is deliberate: these were once two
separate scripts, each with its own copy of the loading and detection logic, and
they drifted apart until a bug fixed in one quietly survived in the other.

## What a "measurement" is

A single scan is divided into **row groups** — horizontal bands, 20 by default.
Each band is averaged into its own profile and analysed separately, so one image
yields 80–100 blaze-angle measurements rather than the 4–10 you would get by
averaging the whole image into a single profile.

This is where the statistics get interesting, and where most of the care in this
project has gone. Those measurements are *not* independent of one another: every
band re-measures the same physical grooves. See **The correlation correction**.

## Where to go next

- **Row groups and scan edges** — what row-group analysis buys, and why grooves
  at the edge of a scan are thrown away
- **The correlation correction** — why the error bars are wider than √N suggests
- **Facet fitting and trim** — how the blaze facet is chosen and measured, and
  the one parameter that moves the answer most
- **Reading the outputs** — what every column and file means
