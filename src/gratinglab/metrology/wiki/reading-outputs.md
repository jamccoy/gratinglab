# Reading the outputs

Every run writes into `results/`, timestamped. Nothing is overwritten.

## `analysis_data_<timestamp>.csv`

One row per sample. The summary table.

| Column | Meaning |
|---|---|
| `Sample` | Label from `SAMPLES_TO_COMPARE` |
| `File` | Source scan(s); several are joined with `;` when scans were combined |
| `N_scans` | How many scans went into this row |
| `N_grooves` | Total measurements — bands × grooves, **not** independent samples |
| `Mean_angle_deg` | The measurement |
| `Std_angle_deg` | Spread of individual measurements (σ) |
| `Min_angle_deg`, `Max_angle_deg` | Range. A minimum far below the rest is the classic sign of a bad fit |
| `Mean_slope` | dy/dx of the fitted facet |
| `Period_nm`, `Period_std_nm` | Measured groove spacing |
| `Mean_depth_nm` | Groove depth |
| `Mean_facet_width_nm` | Blaze facet width. Much smaller than the rest means clipped grooves |
| `Local_angle_std_deg`, `Local_angle_range_deg` | Within-facet variation (camber) |
| `Temperature_C` | From the sample definition |
| `ICC` | Intraclass correlation — how correlated the row groups are |
| `N_eff` | Effective sample size |
| `SEM_deg` | Naive standard error, σ/√N |
| `SEM_corrected_deg` | **The one to quote**: σ/√N_eff |

The last four columns are appended rather than inserted, so results recorded
before the correction still diff cleanly on columns 1–16.

### SEM vs SEM_corrected

`SEM_deg` assumes the measurements are independent. They are not.
`SEM_corrected_deg` accounts for that and is 1.17–1.65× larger on this project's
data. Every significance test uses the corrected value.

Quote `SEM_corrected_deg`. `SEM_deg` is retained so the size of the correction is
visible rather than hidden.

## `per_groove_data_<timestamp>.csv`

One row per individual measurement — the raw material behind the summary.

`Sample`, `Scan_file`, `Row_group`, `Groove_number`, `Blaze_angle_deg`,
`Groove_depth_nm`, `Blaze_width_nm`, `Steep_width_nm`, `R2`, `Local_period_nm`.

`Row_group` is what makes the correlation structure analysable outside this
software — it identifies which band each measurement came from. Without it the
ICC cannot be computed at all.

Two quick sanity checks on this file:

- **Sort by `Blaze_angle_deg`.** Anything far from the population is a failed fit,
  not a discovery.
- **Sort by `Blaze_width_nm`.** A facet much narrower than the rest was clipped;
  its angle is meaningless regardless of how good its R² looks.

## `analysis_summary_<timestamp>.txt`

Prose version, plus the settings the run used — `BLAZE_SIDE`, `FACET_TRIM`,
`FLATTEN_METHOD` — so a result can be traced back to the parameters that produced
it. For each sample it reports the uncertainty decomposition: measurement
uncertainty from the fits, physical groove-to-groove variation, the intraclass
correlation, the effective sample size, and both standard errors.

It also carries the pairwise comparisons: differences, Welch t-tests on the
effective sample sizes, and Cohen's d.

## `icc_report_<timestamp>.txt`

Produced by `ANALYSIS_MODE = 'icc'`. Per-scan ICC, effective sample size, and the
factor by which the naive SEM is understated. Diagnostic only — it changes
nothing, it just tells you how much the correction is doing.

## The GUI results panel

Shows the corrected SEM with the ICC and effective N beside it, the spread and
range, the measurement count and mode, the period, and the worst fit R² in the
set. The worst R² is there as an at-a-glance quality check: if it drops well below
0.95, something in the detection or trim settings needs looking at.

## Regenerating and comparing

`results/` is git-ignored — these files are outputs, not source. To check that a
change to the code did not alter the science, diff a fresh run against a stored
baseline; `docs/BASELINES.md` records which baseline to use and why.
