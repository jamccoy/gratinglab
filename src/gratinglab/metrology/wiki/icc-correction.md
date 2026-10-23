# The correlation correction

Row-group analysis produces around 100 blaze-angle measurements per scan.
Treating these as 100 independent samples would understate the uncertainty on
the mean, because they are not independent: every band re-measures the same
physical grooves.

This page describes how the software accounts for that correlation.

## The problem

The standard error of a mean assumes independent observations:

    SEM = σ / √N

Dividing by √100 for 100 correlated measurements claims ten times the precision
of a single measurement. If all 20 bands are essentially re-imaging the same
five grooves, however, the real information content is closer to that of five
carefully measured grooves than to that of 100 independent draws.

The measurements are clearly correlated, so the relevant question is the degree
of correlation, since that determines whether it matters.

## Measuring it: the intraclass correlation

The *intraclass correlation coefficient* (ICC) is the fraction of the total
variance that lies between row groups, as opposed to within them:

    ICC = σ²_between / (σ²_between + σ²_within)

- **ICC = 0**: groups are indistinguishable. The measurements behave
  independently, and dividing by √N is approximately correct.
- **ICC = 1**: all the variance is between groups. Measurements within a group
  carry no new information, and the effective sample size is the number of
  groups instead of the number of measurements.

The ICC is computed in `gratinglab/metrology/stats/icc.py`. Running with
`ANALYSIS_MODE = 'icc'` produces a per-scan report.

## What was measured

On the eight bundled sample scans:

| Sample | Scan | N | ICC | N_eff | SEM × |
|---|---|---|---|---|---|
| Master | ALD_master_1p5um | 80 | 0.244 | 46.2 | 1.32 |
| 150°C | 20250820_150C_00003 | 100 | 0.109 | 69.7 | 1.20 |
| 215°C | 20250820_215C_00001 | 100 | 0.377 | 39.9 | 1.58 |
| 280°C | 20250820_280C_00004 | 102 | 0.128 | 67.0 | 1.23 |
| 280°C | 20250905_280C_00005 | 94 | 0.097 | 69.2 | 1.17 |
| 280°C | 20250905_280C_00004 | 100 | 0.429 | 36.8 | 1.65 |
| 280°C | 20250905_280C_00000 | 80 | 0.357 | 38.6 | 1.44 |
| 500°C | 500C_N2_flatten | 78 | 0.201 | 49.3 | 1.26 |

The ICC ranges over 0.097–0.429, with a median of 0.244, so roughly a quarter of
the variance lies between bands. This is large enough to matter, but not so
large that the row groups are without value.

## The correction

The *design effect* converts an ICC into an effective sample size:

    N_eff = N / (1 + (m − 1) · ICC)

where `m` is the mean number of measurements per group. Standard errors then use
`N_eff` in place of `N`:

    SEM = σ / √N_eff

At ICC = 0 this reduces to `N_eff = N` and the correction vanishes, as it
should. As the ICC rises, `N_eff` falls toward the number of groups.

Every standard error, confidence interval and p-value in the software uses the
effective sample size. As a result, reported SEMs are 1.17–1.65× wider than the
naive value, and effective sizes are 37–70 where the raw sizes are 78–102.

## What it changed

No conclusion changed. On the bundled samples every comparison kept its verdict:

- Master vs each treated sample: highly significant before and after the
  correction. These are ~3° effects against standard errors near 0.25°, and a
  1.6× inflation does not affect them.
- Adjacent temperature steps: null before and after, and slightly more clearly
  null afterwards. The 215 °C vs 280 °C comparison went from p = 0.61 to
  p = 0.74.

This is the expected outcome for a correction that is real but moderate: it
widens the error bars without altering the findings. Had it reversed a result,
that reversal would itself have been reported as a finding.

## Two details worth knowing

Multiple scans of one sample add their effective sizes. Separate scans are
independent images, and the correlation being corrected for exists within a
scan. The 280 °C sample, which combines four scans, therefore has
67.0 + 69.2 + 36.8 + 38.6 = 211.7 effective measurements from 376 raw ones.
Deriving a single design effect over the pooled measurements would incorrectly
treat four separate images as one correlated cluster.

Cohen's d uses the raw counts. It is a standardised mean difference (a
description of how far apart two samples are in units of their own spread) and
is not an inference statistic. The correlation is therefore accounted for in the
standard error and not in the pooled variance.

## If you want fewer correlated measurements

Reducing `N_ROW_GROUPS` does not help much, because fewer, thicker bands are
each less noisy but there are fewer of them. The effective sample size is raised
by measuring more grooves (a wider scan) or by taking more independent scans of
the same sample. Both of these add new information, whereas more bands over the
same grooves mostly do not.
