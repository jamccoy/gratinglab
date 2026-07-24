# Statistical Decision Guide: Which Approach for Your AFM Data?

## The Key Question You Raised

**"Grooves can vary independently, even adjacent ones, due to manufacturing."**

This is **critically important** and changes everything!

## Two Scenarios

### Scenario A: Independent Grooves (Your Case?)

**Physical reality:**
- Each groove's blaze angle is determined independently during fabrication
- Adjacent grooves can differ significantly
- No spatial correlation in the actual sample

**What creates ICC in this case:**
- ONLY measurement/processing artifacts
- Shared flattening algorithm per row-group
- Common noise in that image region

**ICC should be:** LOW (< 0.1)

**Correct analysis:** 
- If ICC < 0.1 → Use SIMPLE statistics!
- All ~100 measurements are essentially independent
- No penalty needed

### Scenario B: Spatially Correlated (Might not be your case)

**Physical reality:**
- Blaze angle varies smoothly across sample
- Spatial gradients (e.g., temperature during fabrication)
- Adjacent grooves are more similar

**ICC should be:** MODERATE to HIGH (> 0.2)

**Correct analysis:**
- Must use hierarchical statistics
- Effective n < 100

## How to Find Out Which You Have

### Quick Test (5 minutes):

Run this on your data:

```python
from correlation_source_analysis import quick_icc_check

# For each of your samples:
need_hierarchical = quick_icc_check(
    angles=your_angles,
    row_groups=your_row_groups,
    label="Master Sample"
)

if not need_hierarchical:
    print("Good news! Use simple statistics.")
    sem = np.std(your_angles) / np.sqrt(len(your_angles))
else:
    print("Use hierarchical statistics.")
    sem = sem_best  # from hierarchical analysis
```

### Full Analysis (20 minutes):

```python
from correlation_source_analysis import analyze_correlation_source

analysis = analyze_correlation_source(
    angles=your_angles,
    row_group_labels=your_row_groups,
    spatial_positions=your_positions,  # If available
    label="Sample Name"
)

# This will:
# 1. Calculate ICC
# 2. Test for spatial variation (ANOVA)
# 3. Check autocorrelation between adjacent grooves
# 4. Create diagnostic plots
# 5. Recommend approach
```

## Interpreting Results

### If ICC < 0.1:

**Interpretation:**
- ✓ Grooves are essentially independent
- ✓ Manufacturing variation dominates
- ✓ Row-grouping is not creating artificial correlation

**Statistical Approach:**
```python
# SIMPLE - All measurements independent
n = len(all_angles)
mean = np.mean(all_angles)
std = np.std(all_angles, ddof=1)
sem = std / np.sqrt(n)  # Use this!

print(f"Mean: {mean:.3f}° ± {sem:.3f}° (SEM, N={n})")
```

**For comparisons:**
```python
# Standard t-test (equal or unequal variance)
sem1 = std1 / np.sqrt(n1)
sem2 = std2 / np.sqrt(n2)
se_diff = np.sqrt(sem1**2 + sem2**2)

t_stat = (mean2 - mean1) / se_diff
# ... standard t-test ...
```

**Report in paper:**
> "Each groove represents an independent measurement (ICC = 0.03, indicating negligible spatial correlation). We obtained N = 100 independent measurements per sample."

### If 0.1 < ICC < 0.2:

**Interpretation:**
- ≈ Weak but non-negligible correlation
- Manufacturing variation dominates but slight spatial trends exist

**Statistical Approach:**
```python
# HYBRID - Report both for transparency
sem_simple = std / np.sqrt(n_total)  # Optimistic
sem_hierarchical = ...  # From hierarchical analysis (conservative)

print(f"Mean: {mean:.3f}°")
print(f"  SEM (simple): {sem_simple:.3f}°")
print(f"  SEM (hierarchical): {sem_hierarchical:.3f}°")
print(f"  (ICC = {icc:.2f})")
```

**Report in paper:**
> "Accounting for weak spatial correlation within images (ICC = 0.15), we obtained n_eff = 85 effective independent measurements from 100 total measurements per sample."

### If ICC > 0.2:

**Interpretation:**
- ! Significant spatial correlation
- May indicate fabrication gradients or measurement artifacts

**Statistical Approach:**
```python
# HIERARCHICAL - Essential
stats = calculate_hierarchical_statistics(angles, row_groups, meas_errors)
sem = stats['sem_best']  # Must use this

print(f"Mean: {mean:.3f}° ± {sem:.3f}° (SEM)")
print(f"  Effective N = {stats['n_effective']:.1f}")
print(f"  (ICC = {icc:.2f})")
```

**Report in paper:**
> "Measurements exhibited spatial clustering (ICC = 0.35). Using hierarchical statistics to account for this correlation, we obtained n_eff = 45 effective independent measurements from 100 total measurements per sample."

## What I Expect for Your Data

Based on your description ("grooves can vary independently"), I predict:

**ICC ≈ 0.05 to 0.15**
- Mostly independent
- Slight correlation from measurement process
- n_effective ≈ 80-95 (out of 100)

This means:
- ✓ Your row-grouping strategy is working well
- ✓ You're getting nearly all the benefit of 100 measurements
- ≈ Small correction for measurement correlation
- ✓ Can use "best estimate" SEM confidently

## Practical Workflow

### Step 1: Run ICC Check on All Samples

```python
# Quick check
for sample in your_samples:
    quick_icc_check(sample['angles'], sample['row_groups'], sample['label'])
```

Expected output:
```
Master: ICC = 0.08 → Simple statistics OK ✓
150°C:  ICC = 0.12 → Hierarchical stats recommended ≈
280°C:  ICC = 0.09 → Simple statistics OK ✓
```

### Step 2: Choose Strategy

**If all ICC < 0.1:**
- Use simple statistics for everything
- Report: "Measurements were independent (ICC < 0.1 for all samples)"

**If all ICC < 0.2:**
- Use hierarchical for comparisons (conservative)
- Report both SEMs for transparency

**If any ICC > 0.2:**
- Must use hierarchical throughout
- Investigate why (real gradient vs artifact?)

### Step 3: Update Your Code

**Option A: ICC < 0.1 (Simple Stats)**

```python
# In your analyzer - minimal changes
result = {
    'mean_angle': np.mean(all_angles),
    'std_angle': np.std(all_angles, ddof=1),
    'sem': np.std(all_angles, ddof=1) / np.sqrt(len(all_angles)),
    'n_grooves': len(all_angles),
    'icc': icc,  # Report for transparency
}

# In comparisons - standard t-test
sem_combined = np.sqrt(sem1**2 + sem2**2)
t_stat = diff / sem_combined
# ... proceed as normal ...
```

**Option B: ICC > 0.1 (Hierarchical)**

```python
# Use the full hierarchical statistics
# (As in the integration guide I provided earlier)
```

## The Bottom Line

**Your question revealed an important insight:**

If your grooves are truly independent physically, then:
1. ICC should be low
2. You CAN use (nearly) all 100 measurements as independent
3. The hierarchical approach will correctly identify this (via low ICC)
4. You'll get credit for all your measurements!

The hierarchical statistics framework is **flexible** - it will:
- Give you ~100 effective measurements if grooves are independent (ICC ≈ 0)
- Give you ~20-50 effective measurements if they're correlated (ICC ≈ 0.5)
- **Automatically determine which** based on your data

## Next Steps

1. **Run the ICC check** on one representative sample
2. **Look at the diagnostic plots** - are group means all similar?
3. **Check ANOVA** - are groups significantly different?
4. Based on results, **choose simple or hierarchical**
5. Apply consistently to all samples

Want me to help you run this analysis on your actual data to see which scenario you're in?
