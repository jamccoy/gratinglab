# Verification baselines

There is no test suite for the blaze-angle path. The way to prove a change did
what you intended — and nothing else — is to diff its output against a stored
reference run. This file records which reference to use and what each one means.

`results/` is git-ignored, so **these files exist only on disk. Don't delete them.**

## The baselines

| File | What it is |
|---|---|
| `results/analysis_data_20260318_162808.csv` | Last known-good run before the project was revived (18 Mar 2026). Predates every change made since. |
| `results/baseline_prefix_analysis.csv` | Immediately **before** the scan-edge fix. Still contains the four broken 2.4-3.6° measurements. |
| `results/baseline_prefix_per_groove.csv` | Per-groove detail for the same run. |
| `tests/fixtures/rev3_reference.ggp` | Output of `afm_scan_avg_profile_rev3.py`, the standalone script the `boundary/` code replaced. Checked into git. |

## Which one to use

**Refactoring — no behaviour change intended.** Diff against the *most recent*
run before your change. Any difference is a bug in the refactor.

```bash
MPLBACKEND=Agg .venv/bin/python main.py
diff <(cut -d, -f1,3- results/<before>.csv) <(cut -d, -f1,3- results/<after>.csv)
```

Field 2 is the filename column; the timestamp lives in the filename rather than
the data, so everything else must match exactly. This held across the jump from
Python 3.9 to 3.12 with numpy 2.5 / scipy 1.18, so genuine behaviour-preserving
changes really do reproduce bit-for-bit.

**Changing the science deliberately.** Capture a baseline first, then diff and
read every difference to confirm each one is the change you meant:

```bash
cp results/analysis_data_<latest>.csv results/baseline_<what>.csv
```

**Boundary/PCGrate changes.** Run the test suite; it compares against the fixture
automatically.

```bash
.venv/bin/python tests/test_ggp_equivalence.py
```

## What the current numbers mean

The 18 Mar baseline and `baseline_prefix_*` both **predate the scan-edge fix**, so
they contain measurements taken from grooves clipped by the edge of the scan —
including four with angles of 2.4-3.6°, physically impossible for a ~30° blazed
facet. Current output will not match them, and should not.

Headline effect of that fix: N 843 → 734, minimum angle 2.40° → 25.58°, no fit
below R² 0.95, and the 280°C standard deviation falling from 3.04 to 1.65. The
Master sample is unchanged to four decimal places, because it had no edge grooves
to reject.

## What these baselines do and don't prove

A diff proves output *changed* or *didn't*. It says nothing about whether the
statistics are right.

The row-group independence problem that this warning used to describe is now
**fixed**: measurements are clustered by row group, the ICC is measured, and every
standard error and p-value uses the effective sample size. The correction was
applied on 2026-08-06.

Note this when diffing across that date: `analysis_data_*.csv` gained four
columns (`ICC`, `N_eff`, `SEM_deg`, `SEM_corrected_deg`) appended after
`Temperature_C`. Columns 1–16 are unchanged and still diff clean against every
earlier baseline — cut to those columns when comparing across the change.

## Verified dependency versions

`pyproject.toml` gives lower bounds so the package installs broadly. The stored
baselines were produced with these exact versions, on Python 3.12.9:

```
numpy==2.5.1  scipy==1.18.0  matplotlib==3.11.1  PySide6-Essentials==6.11.1
```

They replaced a pinned `requirements.txt`, which was retired when the project
became pip-installable. If a baseline diff fails and nothing in the code changed,
check these first.
