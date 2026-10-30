# Reading the outputs

Every run writes timestamped files into `results/`, and no existing file is
overwritten. This page describes each output file in turn, followed by the GUI
results panel.

## `analysis_data_<timestamp>.csv`

This is the summary table, with one row per sample.

| Column | Meaning |
|---|---|
| `Sample` | Label from `SAMPLES_TO_COMPARE` |
| `File` | Source scan(s); several are joined with `;` when scans were combined |
| `N_scans` | How many scans went into this row |
| `N_grooves` | Total measurements: bands × grooves, **not** independent samples |
| `Mean_angle_deg` | The measurement |
| `Std_angle_deg` | Spread of individual measurements (σ) |
| `Min_angle_deg`, `Max_angle_deg` | Range. A minimum far below the rest is the classic sign of a bad fit |
| `Mean_slope` | dy/dx of the fitted facet |
| `Period_nm`, `Period_std_nm` | Measured groove spacing |
| `Mean_depth_nm` | Groove depth |
| `Mean_facet_width_nm` | Blaze facet width. Much smaller than the rest means clipped grooves |
| `Local_angle_std_deg`, `Local_angle_range_deg` | Within-facet variation (camber) |
| `Temperature_C` | From the sample definition |
| `ICC` | Intraclass correlation: how correlated the row groups are |
| `N_eff` | Effective sample size |
| `SEM_deg` | Naive standard error, σ/√N |
| `SEM_corrected_deg` | **The one to quote**: σ/√N_eff |

The last four columns are appended at the end of the table, so that results
recorded before the correction still diff cleanly on columns 1–16.

### SEM vs SEM_corrected

`SEM_deg` assumes that the measurements are independent, which they are not.
`SEM_corrected_deg` accounts for the correlation and is 1.17–1.65× larger on
this project's data. Every significance test uses the corrected value.

`SEM_corrected_deg` is the value to quote. `SEM_deg` is retained so that the
size of the correction remains visible.

## `per_groove_data_<timestamp>.csv`

This file contains one row per individual measurement, and is the raw material
behind the summary. Its columns are:

`Sample`, `Scan_file`, `Row_group`, `Groove_number`, `Blaze_angle_deg`,
`Groove_depth_nm`, `Blaze_width_nm`, `Steep_width_nm`, `R2`, `Local_period_nm`.

`Row_group` identifies the band from which each measurement came, and it is
this column that allows the correlation structure to be analysed outside this
software. Without it the ICC cannot be computed.

Two quick sanity checks can be made on this file:

- **Sort by `Blaze_angle_deg`.** A value far from the rest of the population is
  a failed fit and should not be interpreted as a physical result.
- **Sort by `Blaze_width_nm`.** A facet much narrower than the rest was clipped,
  and its angle has no meaning however good its R² appears.

## `analysis_summary_<timestamp>.txt`

This is a prose version of the results, together with the settings the run used
(`BLAZE_SIDE`, `FACET_TRIM`, `FLATTEN_METHOD`), so that a result can be traced
back to the parameters that produced it. For each sample it reports the
uncertainty decomposition: the measurement uncertainty from the fits, the
physical groove-to-groove variation, the intraclass correlation, the effective
sample size, and both standard errors.

It also contains the pairwise comparisons: differences, Welch t-tests on the
effective sample sizes, and Cohen's d.

## `icc_report_<timestamp>.txt`

This report is produced by `ANALYSIS_MODE = 'icc'`. It gives the per-scan ICC,
the effective sample size, and the factor by which the naive SEM is understated.
It is diagnostic only: it changes no result, and it shows how large an effect
the correction has.

## The boundary profile: `.ggp`, `.json`, `.csv` and its metrics

The boundary profile is produced by the Boundary tab or by
`ANALYSIS_MODE = 'ggp'`. Both routes share one computation, so their results
cannot drift apart. "Export…" writes any of three formats, chosen by the file
type in the save dialog:

| Format | Holds | For |
|---|---|---|
| `.ggp` | the shape only, normalised to one period | PCGrate |
| `.json` | the shape at full precision, the measured period in nm, the groove counts, the metrics, and provenance (gratinglab version, source scan, analysis settings) | reloading into gratinglab's efficiency window or `read_boundary_json`, with nothing to retype |
| `.csv` | one row per point: `t`, `y` as fractions of the period and `x_nm`, `y_nm` in nm, with the period and groove counts in `#` lines after the header | a spreadsheet, MATLAB, or another code |

A `.ggp` or `.csv` export is accompanied by the metrics file described below,
while the `.json` contains the metrics internally. "Send to efficiency window"
opens the groove in the main gratinglab window with the measured period already
filled in, by way of a `.json` written to the temporary folder. Only the `.json`
and that button preserve the period: a `.ggp` cannot hold it, so the period must
be entered again for a profile loaded from one.

### `averaged_groove_profile_<sample>.ggp`

This is the file PCGrate reads. It is one `x y` pair per line, with no header,
which is what PCGrate's border-profile import takes:

```
0.000000 0.000000
0.000500 0.000001
```

Files written before 0.3.0 may start with the legacy
`3 0 - Polygonal type` / `Period: 1 PSC: 1` lines; they still load. A script
that needs those lines can call `write_ggp(..., header=True)`.

Both columns are normalised: x runs from 0 to 1 across exactly one period, and
y is the height as a fraction of that period. A groove 0.33 deep in this file is
therefore a third as deep as the grating's period is wide, whatever that period
is in nanometres.

The header must not be commented. Writing this file with
`np.savetxt(header=...)` prepends `# `, which PCGrate rejects. Three different
header variants were found among the hand-made files before this step was
automated: the correct one, a `#`-commented one, and one carrying only the numpy
header with no PCGrate lines at all.

The endpoints are forced to exactly zero so that successive periods tile without
a step. A seam at the period boundary is a defect that makes an efficiency curve
wrong while the file still appears entirely reasonable, and for this reason the
panel plots the profile in addition to reporting numbers.

### `<sample>_metrics.txt`

| Field | Meaning |
|---|---|
| Period | Measured groove spacing, nm |
| Grooves averaged | How many went into the average |
| Groove depth | Peak height, as a fraction of the period |
| Peak-to-valley | Full range, as a fraction of the period |
| RMS slope | Root-mean-square gradient of the normalised profile |
| Max slope magnitude | Steepest gradient |
| Max curvature | Sharpest bend, a proxy for how rounded the groove is |

### Partial grooves are skipped

`average_grooves` averages every groove over the narrowest window in the set,
and the normalisation step then stretches that window to span a full period. A
single groove cut into by the scan edge would therefore narrow the window for
all of the grooves and distort the exported shape. For this reason the function
skips any groove whose full period does not fit inside the scan, independently
of the edge-exclusion setting. The summary's "n of m grooves averaged" counts
the skipped grooves.

Before this check was added, edge exclusion was the only protection, and edge
exclusion is a user setting. On the TASTE scan the last groove sits 38 samples
from the edge against a half-width of 40, which stretched the exported profile
horizontally by 5.5%. On the synthetic sawtooth, one partial groove narrowed the
window to 0.66 of a period, turned a 70° anti-blaze facet into 58° and reduced
the depth by 13%. A change of this size in the groove's shape is a real change
in the efficiency it predicts, and nothing in the file would reveal it.

### Controls worth understanding

- **Points in the profile**: how finely the averaged groove is resampled. This
  is the number of lines in the file.
- **Smoothing**: light smoothing that wraps at the period boundary, removing
  kinks left by interpolating each groove onto a common axis. It does change
  the exported curve.
- **Minimum half-width**: grooves whose symmetric extent is at or below this
  value are skipped. On a typical scan the exclusion is all-or-nothing, because
  the surviving grooves all share one extent.

## The GUI results panel

The panel shows the corrected SEM with the ICC and effective N beside it, the
spread and range, the measurement count and mode, the period, and the worst fit
R² in the set. The worst R² serves as an at-a-glance quality check: if it falls
well below 0.95, the detection or trim settings should be examined.

## Regenerating and comparing

`results/` is git-ignored, since these files are outputs and not source. To
confirm that a change to the code did not alter the scientific results, a fresh
run should be compared by diff against a stored baseline; `docs/BASELINES.md`
records which baseline to use and why.
