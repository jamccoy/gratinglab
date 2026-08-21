# Historical documents

These describe migrations and feature installations that have **already been completed**.
They are kept for reference on how the code reached its current shape, but none of them
are instructions to follow today — the steps in them are done.

| File | What it documented |
|---|---|
| `RESTRUCTURING_GUIDE.md` | Migrating from a flat script layout to the `afm_analysis/` package. Done. |
| `REFACTORING_SUMMARY.md` | Summary of that same restructuring. Done. |
| `README_ROW_GROUPS.md` | Installing row-group analysis. Done — it is live and enabled (`USE_ROW_GROUPS = True`). |
| `IMPLEMENTATION_GUIDE.md` | Step-by-step install for the row-group feature. Done. |
| `EXPECTED_RESULTS.md` | Before/after comparison used to validate the row-group rollout. Done. |
| `stat_analysis.md` | January 2026 review of the statistics. It opens by naming "3 major statistical issues that need addressing" — all since resolved. It predates the scan-edge fix and the entire ICC correction, so anyone reading it today would conclude the statistics are still broken. Superseded by the **Wiki** tab in the application (`src/afm_analysis/wiki/`). |

For the current state of the project, see the top-level `README.md`.
For where things stand and what is still open, see `docs/PROGRESS_SUMMARY.md`.
