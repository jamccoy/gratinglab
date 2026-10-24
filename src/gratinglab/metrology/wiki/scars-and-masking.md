# Scars and masking

Two things used to send a scan back to Gwyddion before it could be analysed
here: scan lines the tip got wrong, and a patch of the surface that is not
grating at all. Both are now settings, set in the Import tab and carried in
`AnalysisSettings`, so a batch or command-line run reproduces exactly what was
chosen in the window, and both are recorded in every summary file they affect.

## Where they sit in the pipeline

```
load
  ↓  region mask
  ↓  image flattening (none / plane / align_rows)
  ↓  scar repair (none / interpolate / drop)
  ↓  groove orientation (measured; corrected when it matters)
  ↓  tip correction
  ↓  rows averaged into row groups
```

Scar repair comes after levelling because an un-levelled row offset looks
exactly like a scar, and every line would be marked. It comes before tip
correction because erosion spreads each pixel sideways by the tip's reach, which
would widen a scar into its neighbours before anything could find it.

The mask comes first of all, for the reason given under "What happens to
excluded pixels" below.

## Scar repair

A *scar* is a run of one to four adjacent scan lines that stands above (or
below) the line on each side of it: the tip skipping over a contaminant,
feedback ringing after a step, a line where the Z piezo briefly lost track. This
is the same idea as Gwyddion's *Mark scars* and *Remove scars*.

It is a stage with its own setting rather than one more image-flattening method.
Image flattening is a choice of one, and scar repair is used *together* with a
levelling method, normally `align_rows`; as a method it would have had to be
offered once per levelling method.

### How scars are found

1. The lines are levelled against each other for this purpose only, by the
   median *step* from each line to the next. `align_rows` takes each line's
   median height, and a grating's heights spread over its whole depth, so a scar
   along part of a line moves that median and leaves the rest of the line offset
   by tens of nanometres. The median step is not moved by a scar covering under
   half the line.
2. The scale of the line-to-line noise is measured robustly, as 1.4826 times the
   median absolute deviation of the vertical differences, so neither the scars
   nor a steep facet can inflate it.
3. A pixel is a candidate when every line of its band lies above both the line
   just above the band and the line just below it by more than
   **threshold** × noise (default 4), or below both by as much.
4. Candidates that are only the shadow of a stronger scar are set aside (below).
5. A run along a line is kept when its strong core, at the full threshold, is at
   least **min length** pixels long (default 16), and it is then grown along the
   line while the test holds at half the threshold.

The length rule is what keeps noise out: a 4σ excursion on both sides is rare,
and sixteen of them in a row along a line essentially never happens by chance.
On the committed synthetic scan, which carries only noise, nothing is marked, and
neither is a ten-pixel spike.

A vertical *maximum* is required, not merely a large vertical difference, so a
groove running at a slight angle to the slow axis, which shifts from line to line
monotonically, is not mistaken for a scar.

**Shadows.** That same angle is what makes step 4 necessary. On a rotated
grating neighbouring lines differ by several nanometres at every steep facet, so
the lines on the far side of a raised scar sit below both their neighbours there
and look like a lowered scar about one facet long at every groove. On a grating
rotated 3° with three scarred lines, eight more lines were marked before this
step. A vertical run of candidates that touches a run of the opposite sign at
least twice as strong is taken for its shadow. Twice, not merely stronger,
because two raised lines with an untouched one between them look locally exactly
like one lowered line, with the same margin, and only the length rule can tell
them apart. Scars are repaired before the groove orientation is corrected, which
shifts lines sideways and would carry each scar with its line.

### What is done with them

| Setting | What happens to a scar line |
|---|---|
| `none` | nothing; the default |
| `interpolate` | the line is first re-levelled against its clean neighbours on its unmarked pixels, removing whatever offset the scar gave its median, then each marked pixel is refilled down its own column, linearly between the nearest unmarked pixels above and below |
| `drop` | the whole line leaves the analysis and joins no row group |

Interpolation runs down the column because that is the direction the grooves
run, so the lines either side are the best estimate of what the tip should have
seen. `drop` invents nothing and costs the data on those lines.

The count of lines repaired is shown under the control, ticked in orange at the
right edge of the image, and written to the results panel, the summary file and
the boundary metrics file.

### How much it matters

Measured on the synthetic scan (period 315 nm, blaze 30°, 64 lines) with seven
lines scarred by +60 nm over a third of their length:

| | Mean blaze angle | Averaged groove, worst point |
|---|---|---|
| clean scan | 29.8955° | reference |
| scarred, no repair | 29.7309° | 1.01 nm off |
| scarred, `interpolate` | 29.8951° | 0.006 nm off |
| scarred, `drop` | 29.8937° | 0.007 nm off |

A third of a line, not a half: over half, which part of the line is the scar is a
matter of opinion, and the median step takes the majority's.

On the same grating rotated 3° with three lines scarred over a third of their
length, exactly those three lines are marked, and with `interpolate` and the
orientation correction the blaze angle comes back within 0.02° of 30°.

## Region mask

Draw on the image in the Import tab to leave part of the scan out. The drag mode
chooses what a drag adds:

| Mode | Setting | Effect |
|---|---|---|
| rectangle | `mask_rect = (row0, row1, col0, col1)` | that patch, such as a particle; a new rectangle replaces the old one |
| rows | `mask_rows = ((row0, row1), ...)` | whole scan lines; they join no row group |
| columns | `mask_cols = ((col0, col1), ...)` | whole columns; no groove is measured across them |

The values are **pixel indices of the image as loaded**, half-open like a Python
slice, with line 0 at the top of the displayed image (which is how the y axis is
labelled). Pixels rather than micrometres because the slow axis has no reliable
physical size in every format, and because a mask belongs to one scan. The
panel shows the values to copy into `config.py` (`MASK_RECT`, `MASK_ROWS`,
`MASK_COLS`). Loading a different file clears the mask.

### What happens to excluded pixels

They become NaN straight after loading, before image flattening, because
levelling reads whole lines: `align_rows` takes each line's median and `plane`
fits every pixel, so a particle left in would move every line it sits on before
any later step could hide it. Both levelling methods skip NaN. Tip erosion cannot
take a NaN, so it is handed a copy with each gap bridged along its line, and the
NaN is put back afterwards. From there:

- A line that is excluded end to end is removed before the row groups are cut,
  so the groups are formed from what is left and stay equal in size.
- Within a row group, each column is averaged over the lines not excluded there.
  A rectangle therefore removes only its own pixels; the rest of each column
  still counts.
- A column excluded across a whole row group is a gap in that group's profile.
  Any groove whose analysis window, 0.65 of a period either side of its centre,
  reaches into a gap is discarded, and the number discarded is reported.
- The two grooves either side of a gap are not neighbours. Their spacing is
  several periods, and taken as a period it widened both grooves' windows into
  the gap (they read 7.6° and 25.6° on a 30° grating) and inflated the measured
  period. A gap is treated as a scan edge: that spacing counts towards neither
  the period nor either groove's window.
- To flatten the profile and find grooves at all, a gap is bridged by a straight
  line. The bridge is never measured: every groove that would read it is the one
  discarded above, profile flattening fits its background to measured points
  only and takes no peak or trough within two pixels of a gap, and the boundary
  export fits its residual tilt to measured columns only. With the bridge left
  in the fits, a gap of 2.5 periods moved the blaze angle by up to 1.3°.

### How much it matters

The same synthetic scan with a 150 nm Gaussian particle on lines 20–44, columns
100–160:

| | Mean blaze angle | Measurements |
|---|---|---|
| clean scan | 29.8955° | 100 |
| particle, unmasked | 30.3866° | 98 |
| particle, rectangle masked | 29.9042° | 79 |
| particle, those lines masked | 29.8947° | 65 |

Masked, the particle scan gives exactly the result of the clean scan with the
same mask: the particle has no influence left. What remains against the
unmasked clean scan is the sampling difference of measuring fewer grooves. With
whole columns masked only 40 measurements remain, which is the cost of that
choice on a scan two micrometres wide.

## Where it is recorded

| Output | Scar repair | Region mask |
|---|---|---|
| results panel | method, lines | the mask, grooves rejected |
| `analysis_summary_*.txt` | method, lines, pixels | the mask, pixels, lines out of grouping, grooves rejected |
| `groove_analysis_metrics_*.txt` | method, lines | the mask, pixels, grooves rejected |

Each line appears only when the stage did something, so a run without them
writes exactly what it always wrote.
