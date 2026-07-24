# Hierarchical statistics — NOT INTEGRATED

**Status: staged, unused.** Nothing in `afm_analysis/` imports from this directory.
This was previously the `_new/` folder.

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
