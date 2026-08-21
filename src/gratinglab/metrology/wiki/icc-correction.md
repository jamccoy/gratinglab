# The correlation correction

Row-group analysis produces around 100 blaze-angle measurements per scan. Treating
those as 100 independent samples would understate the uncertainty on the mean,
because they are not independent — every band re-measures the **same physical
grooves**.

This page explains what the software does about that.

## The problem

The standard error of a mean assumes independent observations:

    SEM = σ / √N

Divide 100 correlated measurements by √100 and you claim ten times the precision
of a single measurement. But if all 20 bands are essentially re-photographing the
same five grooves, the real information content is closer to five grooves
measured carefully than to 100 independent draws.

The question is not whether the measurements are correlated — they obviously are
— but **how much**, because that determines whether it matters.

## Measuring it: the intraclass correlation

The intraclass correlation coefficient (ICC) is the fraction of the total
variance that sits *between* row groups rather than *within* them:

    ICC = σ²_between / (σ²_between + σ²_within)

- **ICC = 0** — groups are indistinguishable. The measurements behave
  independently and dividing by √N is roughly right.
- **ICC = 1** — all the variance is between groups. Measurements within a group
  carry no new information; the effective sample size is the number of *groups*,
  not the number of measurements.

Computed in `gratinglab/metrology/stats/icc.py`. Run `ANALYSIS_MODE = 'icc'` to produce a
per-scan report.

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

Range 0.097–0.429, median 0.244. Roughly a quarter of the variance sits between
bands. Substantial enough to matter, not so large that the row groups are
worthless.

## The correction

The **design effect** converts an ICC into an effective sample size:

    N_eff = N / (1 + (m − 1) · ICC)

where `m` is the mean number of measurements per group. Standard errors then use
`N_eff` in place of `N`:

    SEM = σ / √N_eff

At ICC = 0 this reduces to `N_eff = N` and the correction vanishes, as it should.
As ICC rises, `N_eff` falls toward the number of groups.

Every standard error, confidence interval and p-value in the software uses the
effective sample size. Reported SEMs are 1.17–1.65× wider than the naive value as
a result, and effective sizes are 37–70 rather than 78–102.

## What it changed

**No conclusion.** On the bundled samples every comparison kept its verdict:

- Master vs each treated sample — highly significant before and after. These are
  ~3° effects against standard errors near 0.25°; a 1.6× inflation does not
  threaten them.
- Adjacent temperature steps — null before and after, and slightly more clearly
  null. The 215 °C vs 280 °C comparison went from p = 0.61 to p = 0.74.

That is the expected outcome when a correction is real but moderate: it changes
the error bars without rewriting the findings. Had it flipped a result, that
would have been the finding.

## Two details worth knowing

**Multiple scans of one sample add their effective sizes.** Separate scans are
genuinely independent images — the correlation being corrected for is *within* a
scan. So the 280 °C sample, combining four scans, has
67.0 + 69.2 + 36.8 + 38.6 = 211.7 effective measurements from 376 raw ones.
Re-deriving a single design effect over the pooled measurements would wrongly
treat four separate images as one correlated cluster.

**Cohen's d deliberately keeps the raw counts.** It is a standardised mean
difference — a description of how far apart two samples are in units of their own
spread — not an inference statistic. Correlation belongs in the standard error,
not in the pooled variance.

## If you want fewer correlated measurements

Reducing `N_ROW_GROUPS` does not help much: fewer, thicker bands are each less
noisy but there are fewer of them. The honest way to raise the effective sample
size is more *grooves* — a wider scan — or more independent scans of the same
sample. Both add genuinely new information; more bands over the same grooves
mostly do not.
