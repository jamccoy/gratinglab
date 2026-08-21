# Row groups and scan edges

## Row-group analysis

A traditional AFM profile analysis averages every row of the image into one
height profile, finds the grooves in it, and measures each one. A 2 µm scan of a
315 nm grating gives about six grooves — six measurements per image.

Row-group analysis divides the image into `N_ROW_GROUPS` horizontal bands (20 by
default), averages each band separately, and analyses all of them. The same scan
now yields 80–100 measurements.

**What that buys:** a view of how the groove shape varies down the image, and
enough measurements to see the distribution rather than six points.

**What it costs:** the measurements are not independent — see
**The correlation correction**. More bands over the same grooves does not
straightforwardly mean more information.

Each band is flattened independently, because the background tilt is not
necessarily the same across the image.

## Why grooves at a scan edge are rejected

`EDGE_EXCLUSION_PERIODS` (0.6 by default) discards any groove detected within
that many periods of either end of the scan line.

The groove is real. The problem is that the scan starts or stops part-way
through it, so **its facet does not exist in the data**. A fit to what remains is
not a measurement of anything.

### The failure this prevents

Four measurements in `20250905_280C_00004.txt` came back as blaze angles of
**2.40°, 2.64°, 2.69° and 3.60°** — in a population averaging 30°. A blazed
facet cannot be 2.4°.

Tracing them: the pixel size in that scan is 3.70 nm, and a healthy blaze facet
spans about 50 pixels. Those four spanned **5–6 pixels**. The groove sat at pixel
5 of a 512-pixel profile, the extraction window was clipped at the array
boundary, and the facet-splitting step landed on the very first sample — leaving
a sliver. A straight line fitted through five noisy points returns an arbitrary
shallow slope.

They were not rare noise. They were the visible tail of a broader problem: about
13% of measurements came from windows clipped by a scan edge, most of which
degraded quietly instead of producing an obviously absurd number.

### What rejecting them fixed

- Minimum angle across all measurements: **2.40° → 25.58°**, inside the physical
  range for these gratings.
- No fit now falls below R² 0.95 — every poor fit in the dataset was an edge
  groove.
- The 280 °C standard deviation fell from 3.04° to 1.65°, in line with the other
  samples rather than looking anomalously noisy.

It also puts the samples on equal terms. How many edge grooves a scan contains is
pure luck of where the scan happened to start relative to the grating phase — the
master sample had none, every other sample had 18–26. Leaving them in weighted
the samples differently for no physical reason.

The cost is about 13% of measurements, which barely moves the standard error, and
those were the least trustworthy points in the set.

### It matters for the boundary profile too

The exported PCGrate profile averages grooves onto a shared window, sized by the
*most restrictive* groove in the set. An edge groove narrows that shared window,
and the normalisation step then stretches the result to span a full period. On
the TASTE scan this stretched the exported profile horizontally by **5.5%**.
Excluding it brings that to 0.31%.

Set `EDGE_EXCLUSION_PERIODS = 0` to disable the check and reproduce the old
behaviour.
