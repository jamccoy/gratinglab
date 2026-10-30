# Row groups and scan edges

## Row-group analysis

A traditional AFM profile analysis averages every row of the image into one
height profile, finds the grooves in it, and measures each groove. A 2 µm scan
of a 315 nm grating contains about six grooves, and so gives six measurements
per image.

Row-group analysis divides the image into `N_ROW_GROUPS` horizontal bands (20 by
default), averages each band separately, and analyses all of them. The same scan
then yields 80–100 measurements.

This provides a view of how the groove shape varies down the image, and enough
measurements for the distribution to be examined, where the traditional analysis
gives only six points.

The cost is that the measurements are not independent, as described in "The
correlation correction". Adding more bands over the same grooves does not
straightforwardly add information.

Each band is flattened independently, because the background tilt is not
necessarily uniform across the image.

## Where the leftover lines go

A line count rarely divides by the number of bands. The lines left over are
left out, from the centre line and from the two ends, so that the bands sit
symmetrically about the middle of the image. The count is reported as
`rows_unbanded`, separately from the lines dropped as scars or masked away.

The placement matters because the same scan can arrive either way up: a
Nanoscope `.spm` and the Gwyddion text export of that scan disagree about
which end of the slow axis is stored first. Until 0.3.2 the leftover lines
were given to the last band, which made every boundary move when the rows were
reversed. On an image that `groove_orientation` has sheared, each line is
shifted by an amount that depends on its distance from the centre, so the
reshuffle changed what each band averaged: on one 401-line scan it moved the
mean blaze angle by 0.067°, against the 0.002° the two readings of that scan
otherwise agree to. Centred bands are the same set of bands read either way,
so the measurement no longer depends on the direction of the slow axis.

## Why grooves at a scan edge are rejected

`EDGE_EXCLUSION_PERIODS` (0.6 by default) discards any groove detected within
that many periods of either end of the scan line.

Such a groove is real, but the scan starts or stops part-way through it, so its
facet is not present in the data. A fit to the remaining points does not measure
any physical quantity.

### The failure this prevents

Four measurements in `20250905_280C_00004.txt` returned blaze angles of 2.40°,
2.64°, 2.69° and 3.60°, in a population averaging 30°. A blazed facet cannot
have an angle of 2.4°.

These measurements were traced as follows. The pixel size in that scan is
3.70 nm, and a typical blaze facet spans about 50 pixels, whereas these four
spanned 5–6 pixels. The groove sat at pixel 5 of a 512-pixel profile, the
extraction window was clipped at the array boundary, and the facet-splitting
step placed the split on the very first sample, leaving a narrow sliver. A
straight line fitted through five noisy points returns an arbitrary shallow
slope.

These were not isolated noise. They were the visible part of a broader problem:
about 13% of measurements came from windows clipped by a scan edge, and most of
these were degraded without producing an obviously implausible number.

### What rejecting them fixed

- Minimum angle across all measurements: 2.40° → 25.58°, which lies inside the
  physical range for these gratings.
- No fit now falls below R² 0.95, since every poor fit in the dataset was an
  edge groove.
- The 280 °C standard deviation fell from 3.04° to 1.65°, in line with the
  other samples, where previously it appeared anomalously noisy.

The rejection also places the samples on equal terms. The number of edge grooves
in a scan depends only on where the scan happened to start relative to the
grating phase: the master sample had none, and every other sample had 18–26.
Including them weighted the samples differently for no physical reason.

The cost is about 13% of measurements, which changes the standard error very
little, and those measurements were the least reliable points in the set.

### It matters for the boundary profile too

The exported PCGrate profile averages grooves onto a shared window whose size is
set by the most restrictive groove in the set. An edge groove narrows that shared
window, and the normalisation step then stretches the result to span a full
period. On the TASTE scan this stretched the exported profile horizontally by
5.5%. Excluding edge grooves reduces this to 0.31%.

Setting `EDGE_EXCLUSION_PERIODS = 0` disables the check and reproduces the old
behaviour.
