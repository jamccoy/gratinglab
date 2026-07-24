# AFM Statistical Improvements - Progress Summary

## Overview

**Goal:** Improve statistical analysis of AFM blaze angle measurements to properly quantify uncertainties and leverage all available data.

**Key Issue Identified:** Currently only getting 4 measurements per AFM image (one per groove from averaged profile), when we could be getting hundreds by using row-level data.

---

## What We've Completed

### ✅ Priority 1: Fit Uncertainty Tracking (DONE)

**What changed:**
- Added covariance matrix calculation to `np.polyfit` calls
- Propagated slope uncertainty → angle uncertainty using calculus
- Decomposed total uncertainty into measurement + physical components
- Added new metrics: `blaze_angle_stderr`, `total_std`, `sem`

**Files modified:**
- `afm_analysis/core/analysis.py` 
- `afm_analysis/analyzer.py`

**What you see:**
```
Uncertainty analysis:
  Average measurement uncertainty per groove: 0.167 deg
  Physical variation (groove-to-groove): 0.909 deg
  Total uncertainty (combined): 0.927 deg
  Standard error of mean: 0.463 deg
  95% confidence interval on mean: ±0.908 deg
```

**Status:** ✅ Working correctly, but limited by analysis approach (see below)

---

### ✅ Priority 2a: Bar Chart Error Bars (DONE)

**What changed:**
- Bar charts now show SEM instead of standard deviation
- Clear labels: "Error bars: ±1 SEM - Uncertainty in Mean"
- Bar labels show both SEM and σ: "17.15° ± 0.46° SEM (σ=0.91°)"

**Files modified:**
- `afm_analysis/visualization/statistics.py`

**Status:** ✅ Working correctly

---

## The Fundamental Issue We Discovered

### Current Analysis Approach:
```
AFM Image (512 rows × N columns)
    ↓
Average all rows → 1D profile
    ↓
Find ~4-10 grooves
    ↓
Measure each groove once
    ↓
Result: 4-10 measurements per image
```

### The Problem:
- **Only using 4-10 measurements** per AFM image
- **Throwing away row-level information** by averaging first
- **Cannot assess within-image variation** properly
- SEM is large because n is small, NOT because we lack data

### What We SHOULD Be Doing:
```
AFM Image (512 rows × N columns)
    ↓
Process each row (or row groups) independently
    ↓
Find grooves in each row
    ↓
Measure angles from many rows
    ↓
Result: 100s-1000s of measurements per image
```

---

## Why This Matters

### Current Situation (Your Data):
- **n = 4 grooves** per image (from averaged profile)
- **SEM = 0.463°** (large because n is small)
- **Measurement uncertainty = 0.167°** per groove (from one averaged profile)

### What We're Missing:
- **512 rows of data** in each AFM image
- Each row could give independent measurements
- Could assess **row-to-row variation** within same image
- Could dramatically reduce SEM by having n=100+ instead of n=4

### Example:
If you analyzed grooves from 10 different row groups:
- Instead of 4 measurements → 40 measurements per image
- SEM would decrease by factor of √10 ≈ 3.2×
- 0.463° → 0.146° SEM

If you analyzed each groove from multiple rows:
- Could get 50+ measurements per groove
- Better assess manufacturing consistency
- Much tighter confidence intervals

---

## What Needs To Be Done (The Real Priority 1!)

### Option A: Row-Group Analysis (Simpler)

**Approach:**
1. Divide AFM image into N row groups (e.g., 10 groups of 50 rows each)
2. Average within each group → N independent profiles
3. Find grooves in each profile
4. Measure angles from each profile
5. Result: N× more measurements

**Advantages:**
- Moderate code changes
- Still benefits from averaging to reduce noise
- N can be tuned (10-100 groups)

**Implementation:**
- Modify `raw_data()` to return multiple profiles instead of one
- Process each profile through existing pipeline
- Aggregate results

### Option B: Per-Row Analysis (More Complex)

**Approach:**
1. Process each row (or small groups) independently
2. Find groove positions for each row
3. Measure angles from each row
4. Handle missing grooves, alignment issues

**Advantages:**
- Maximum data utilization
- Can assess row-to-row variation
- Most measurements possible

**Challenges:**
- Individual rows are noisy
- Groove detection may be unreliable
- Need robust alignment/matching
- More complex implementation

### Option C: Bootstrap/Resampling (Statistical)

**Approach:**
1. Keep current averaged profile
2. Use row-level data to estimate uncertainty via bootstrap
3. Resample rows, recompute profile, remeasure angles
4. Build distribution of angle measurements

**Advantages:**
- Proper uncertainty quantification
- Doesn't require detecting grooves in noisy rows
- Statistically rigorous

**Challenges:**
- More sophisticated statistics
- Computational cost
- Interpretation complexity

---

## Recommended Path Forward

### Immediate Next Step: **Option A - Row-Group Analysis**

This gives you the best balance of:
- More measurements (10-50× increase)
- Reasonable implementation complexity
- Maintains existing quality controls
- Easy to understand/interpret

### Implementation Plan:

**Phase 1: Core Function Changes**
1. Modify `raw_data()` to accept parameter `n_groups` 
2. Return list of profiles instead of single profile
3. Update analyzer to loop over profiles

**Phase 2: Aggregation**
4. Collect all measurements across groups
5. Calculate within-image variation
6. Calculate between-image variation (for multiple scans)

**Phase 3: Uncertainty Refinement**
7. Add "within-image variance" to statistics
8. Update SEM calculation to reflect true sample size
9. Update visualizations

### Expected Improvements:

**With 10 row groups per image:**
- 4 grooves × 10 groups = 40 measurements per image
- SEM would decrease: 0.463° → ~0.15° 
- Much better statistical power

**With 20 row groups per image:**
- 4 grooves × 20 groups = 80 measurements per image  
- SEM would decrease: 0.463° → ~0.10°
- Excellent confidence intervals

---

## Current State of Files

### Modified Files (Ready to Use):
1. ✅ `afm_analysis/core/analysis.py` - Has fit uncertainty tracking
2. ✅ `afm_analysis/analyzer.py` - Has uncertainty decomposition
3. ✅ `afm_analysis/visualization/statistics.py` - Has SEM error bars

### Files Still Needed for Complete Solution:
4. ⏳ `afm_analysis/core/processing.py` - Needs row-group extraction
5. ⏳ `afm_analysis/analyzer.py` - Needs multi-profile handling  
6. ⏳ `afm_analysis/data/aggregation.py` - Needs within-image variance
7. ⏳ `afm_analysis/stats/analysis.py` - Already good, but could add more tests

---

## Summary: What We Thought vs. Reality

### What We Thought:
- Problem: Not tracking fit uncertainties
- Solution: Add covariance matrices and error propagation
- Result: Better uncertainty estimates

### What We Actually Have:
- **Real Problem:** Only using 4 measurements per image when we have 512 rows of data!
- **Our Solution:** Added uncertainty tracking (good, but incomplete)
- **What's Needed:** Analyze data at row-group level to get many more measurements

### The Good News:
- Everything we did (Priority 1 & 2a) is correct and useful
- It will work perfectly once we fix the fundamental issue
- The infrastructure is in place for proper uncertainty tracking
- We just need to feed it more measurements!

---

## Questions for You

1. **How many rows** are in your typical AFM image? (512? 1024?)

2. **Would you prefer:**
   - Option A: 10-50 row groups (moderate increase in measurements)
   - Option B: Per-row analysis (maximum measurements, more complex)
   - Option C: Bootstrap resampling (statistical approach)

3. **Do you want to:**
   - Fix this first (row-group analysis) before continuing Priority 2/3?
   - Continue with visualization/comparison improvements?
   - See a proof-of-concept implementation?

4. **For your immediate needs:**
   - Are you comparing multiple samples? (Then current approach may be "good enough" temporarily)
   - Do you need precise absolute angles? (Then row-group analysis is critical)
   - Publishing soon? (Need to decide quickly)

---

## My Recommendation

**Do Option A (Row-Group Analysis) next.** This is the real "Priority 1" we should have done first! It will:

1. Give you 10-100× more measurements per image
2. Dramatically reduce SEM (error bars will shrink appropriately) 
3. Properly use all your AFM data
4. Make all our Priority 1 & 2 improvements actually meaningful

Would you like me to implement this? I can modify the code to do row-group analysis.
