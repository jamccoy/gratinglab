# GEBL dose / write-uniformity metrology (2023) — FROZEN

**Status: not maintained, not imported by anything.** Moved here on 2026-07-24 from
`~/Desktop/nanofab/AFM/`, where these files sat loose among the raw scan folders.
Preserved because it isn't obvious whether the results are still wanted; nothing depends
on it, and it can be deleted or revived without touching the rest of the project.

## What it measured

A different quantity from the rest of this repo. Not blaze angle, not a PCGrate boundary
profile — **groove depth and spacing**, and how they vary with two things:

- **e-beam dose** — `GEBL_array1.py`, over dose steps 40–52, across several exposure runs
- **position along a long grating write** — `TASTE_22longwrite.py`, `TASTE_22longwrite_.py`
  and `afm_panter1.py`, sampling at 1/5/10/15/20/25 mm from the top or bottom of the write

The point was write-quality calibration: does the groove come out the same depth
everywhere on the grating, and how does dose change it.

## Files

| File | Role |
|---|---|
| `AFM_GEBL_functions.py` | Shared library. `raw_data`, `linear_flat_corr`, `find_maxima`, `groov_meas`, `space_meas`, `line_fit`. Imported by the others with `from AFM_GEBL_functions import *`. |
| `GEBL_array1.py` | Dose arrays → `meas_space_array{1,2}.txt`, `meas_space_arrayold{1,2}.txt` |
| `TASTE_22longwrite.py` | Long-write uniformity → `meas_space.txt` |
| `TASTE_22longwrite_.py` | Extended version of the above → `meas_space_longwrite.txt` |
| `afm_panter1.py` | Same measurement for the PANTER sample → `meas_space_panter1.txt` |

The `meas_space*.txt` files here are copies of the **results**, so the numbers survive
even if the scripts are never run again. Format is three rows: dose (or position),
measured value in nm, and uncertainty.

## Why it wasn't ported to the shared core

Two reasons, both of which a future port would have to deal with:

1. **Different stack.** These use `uncertainties` (`ufloat`, `unumpy`) and `scipy.odr` for
   orthogonal-distance regression. Neither is in `requirements.txt`, so this code will not
   run in the project venv as it stands — `pip install uncertainties` first.
2. **Hard-coded relative paths into the Desktop scan tree.** Every script opens files like
   `'2022-10-17_GEBLbeamtest/3300_40nA_nomdose.0_00003.txt'`, relative to the working
   directory. Those folders are still at `~/Desktop/nanofab/AFM/`. Running these from the
   repo requires either `cd`-ing there or repointing the paths.

## If it comes back

It carries its own copy of the front-end — `raw_data` and `linear_flat_corr` duplicate
what `afm_analysis/core/processing.py` does, with a cruder flattening (two-point linear
level vs. the package's four flattening methods). That duplication is exactly the drift
that caused the scan-edge bug elsewhere in this project.

So a port means: delete `raw_data` / `linear_flat_corr` / `find_maxima` here, take profile
loading, flattening and groove detection from `afm_analysis.core`, and keep only
`groov_meas` / `space_meas` as a third back-end alongside `blaze/` and `boundary/`.
