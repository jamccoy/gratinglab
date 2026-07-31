# Hierarchical statistics — NOT INTEGRATED

**Status: staged, unused.** Nothing in `afm_analysis/` imports from this directory.
This was previously the `_new/` folder.

> **The ICC has now been measured (2026-07-31), and it is substantial.**
>
> Run `ANALYSIS_MODE = 'icc'` to reproduce. Across the eight sample scans the ICC
> ranges **0.097 to 0.429, median 0.244** — i.e. roughly a quarter of the variance
> sits between row groups rather than within them. The measurements are clearly
> not independent.
>
> | Sample | Scan | N | ICC | N_eff | SEM × |
> |---|---|---|---|---|---|
> | Master | ALD_master_1p5um_flatten | 80 | 0.244 | 46.2 | 1.32 |
> | 150°C | 20250820_150C_00003 | 100 | 0.109 | 69.7 | 1.20 |
> | 215°C | 20250820_215C_00001 | 100 | 0.377 | 39.9 | 1.58 |
> | 280°C | 20250820_280C_00004 | 102 | 0.128 | 67.0 | 1.23 |
> | 280°C | 20250905_280C_00005 | 94 | 0.097 | 69.2 | 1.17 |
> | 280°C | 20250905_280C_00004 | 100 | 0.429 | 36.8 | 1.65 |
> | 280°C | 20250905_280C_00000 | 80 | 0.357 | 38.6 | 1.44 |
> | 500°C | 500C_N2_flatten | 78 | 0.201 | 49.3 | 1.26 |
>
> **What this means.** Reported SEMs are understated by **1.17× to 1.65×**, so
> confidence intervals are that much too narrow and p-values correspondingly too
> small. Effective sample sizes are 37–70, not 78–102.
>
> It does not overturn the large results: Master → 150°C is −3.15° against a SEM
> of ~0.25, which survives a 1.6× inflation comfortably. It does bear on the small
> ones — the 215°C → 280°C step of +0.08° was already inside the noise and is more
> clearly so now.
>
> So integrating this directory is warranted. The decision it was waiting on has
> been made.

## The problem it addresses

With `USE_ROW_GROUPS = True` and `N_ROW_GROUPS = 20`, each AFM image yields roughly
20 × 5 ≈ 100 blaze-angle measurements. But those measurements re-measure the **same
physical grooves** from different horizontal bands of the same image — they are not
independent samples.

The live code treats them as if they were:

- `afm_analysis/analyzer.py` → `_calculate_statistics` computes `sem = total_std / sqrt(N)`
  where `N` is the total measurement count.
- `afm_analysis/stats/analysis.py` runs Welch t-tests on that same `N`.

The consequence is that reported SEMs are too small and p-values too significant,
potentially by a large factor. Effect sizes and mean angles are unaffected — this is
purely about the uncertainty on those means.

## What's here

| File | Purpose |
|---|---|
| `correlation_source_analysis.py` | ICC calculation — measures how much of the variance is between row groups vs. within. This is the diagnostic that tells you whether hierarchical stats are actually needed. |
| `improved_statistics.py` | Hierarchical / variance-components replacement for `stats/analysis.py`. |
| `diagnostic_plots.py` | Visualizations for the above. |
| `statistical_decision_guide.md` | Decision tree: simple vs. hierarchical, based on ICC. |
| `statistical_guide.md` | Conceptual background (ICC, hierarchical models). |
| `integration_guide.md` | Step-by-step wiring instructions. |
| `INSTALLATION_README.md` | Overview of the above. |

## Suggested next step

Run the ICC check on an existing result set first. If ICC < 0.1 the row groups are
effectively independent and the current statistics are close enough; if higher, the
hierarchical path in `integration_guide.md` is warranted. Don't integrate blind —
the ICC number decides whether any of this is necessary.
