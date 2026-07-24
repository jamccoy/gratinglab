# Understanding Your AFM Statistics: A Practical Guide

## The Core Issue: Not All Measurements Are Equal

### Your Current Approach
✗ 20 groups × 5 grooves = 100 **independent** measurements
✓ 20 groups × 5 grooves = 100 **correlated** measurements

**Why they're correlated:**
- Same flattening operation per group
- Spatially adjacent within group
- Share local sample properties

## Real-World Analogy

Imagine measuring student heights:

**Scenario A:** Measure 100 students from 20 different schools (5 per school)
- Students within same school are more similar (nutrition, genetics, demographics)
- **Effective sample size** < 100

**Scenario B:** Measure 100 randomly selected students from across a country
- All independent
- **Effective sample size** = 100

Your row-groups are like schools: measurements within groups are more similar.

## The Solution: Intraclass Correlation (ICC)

**ICC measures**: How much correlation exists within groups

- **ICC = 0**: No correlation → all 100 measurements independent
- **ICC = 0.5**: Moderate correlation → ~50 effective measurements
- **ICC = 1.0**: Perfect correlation → only 20 effective measurements (group means)

**My code calculates ICC automatically** from your data!

## Three Estimates of Uncertainty

The improved code gives you three SEM estimates:

### 1. **Conservative** (Most skeptical)
```
SEM = SD(group_means) / sqrt(n_groups)
```
- Treats only group means as independent
- n_groups = 20
- **When to use:** Publication, grant proposals, being cautious

### 2. **Best Estimate** (ICC-adjusted)
```
n_effective = n_total / design_effect
SEM = SD(all data) / sqrt(n_effective)
```
- Accounts for actual correlation
- n_effective ≈ 50-70 (typically)
- **When to use:** Normal reporting, comparisons

### 3. **Liberal** (Most optimistic)
```
SEM = SD(all data) / sqrt(n_total)
```
- Ignores correlation completely
- n_total = 100
- **When to use:** Never! (But shown for comparison)

## Variance Decomposition: Where Does Uncertainty Come From?

```
Total Variance = Measurement² + Physical² + Spatial²
```

**Measurement variance**: From your linear fits
- σ_measurement ≈ 0.05° (small, good!)

**Physical variance**: Real groove-to-groove differences
- σ_physical ≈ 0.2° (within row-group variation)

**Spatial variance**: Row-group to row-group
- σ_spatial ≈ 0.3° (between row-group variation)

**Total SD**: ~0.4° (what you currently report)

### What This Tells You

If **spatial variance >> measurement variance**:
→ Your sample has real spatial non-uniformity (interesting physics!)

If **measurement variance >> physical variance**:
→ Need better fitting or more measurements

## Addressing Your Questions

### Q1: Can I use all 100 measurements?

**Yes, but properly weighted!**

The improved code:
1. Calculates ICC from your data
2. Determines n_effective (usually 50-70, not 100)
3. Uses this for correct SEM and p-values

You still **report** "100 measurements" but **analyze** with n_effective.

### Q2: Multiple scans at same temperature

Since they're the same sample, **combine them**:
- Pool all measurements
- Track which scan each came from
- Calculate both within-scan and between-scan variance

The code handles this with the row_group_labels approach.

### Q3: What precision is acceptable?

With your current setup (σ_total ≈ 0.4°, n_effective ≈ 60):
- **SEM ≈ 0.05°** (uncertainty in mean)
- **95% CI ≈ ±0.1°**

**Can you detect 0.1° difference?**
- No, too small vs measurement noise
- Power analysis suggests need ~500 measurements for 80% power

**Can you detect 0.5° difference?**
- Yes! This is ~1.25× your SD
- High probability of detection

**Recommendation**: Report differences > 0.3° as potentially meaningful

## How to Use the New Code

### Step 1: Calculate hierarchical statistics

```python
from improved_statistics import calculate_hierarchical_statistics

# Your measurement data
angles = [17.2, 17.5, 17.3, ...]  # All measurements
row_groups = [0, 0, 0, 0, 0,  # Group 0: 5 measurements
              1, 1, 1, 1, 1,  # Group 1: 5 measurements
              ...]            # etc.
meas_errors = [0.05, 0.06, ...]  # From blaze_angle_stderr

stats = calculate_hierarchical_statistics(angles, row_groups, meas_errors)
```

### Step 2: Print diagnostics

```python
from improved_statistics import print_hierarchical_statistics

print_hierarchical_statistics(stats, "150°C Sample")
```

This shows:
- Variance decomposition
- ICC and design effect  
- All three SEM estimates
- Recommended confidence intervals

### Step 3: Compare samples

```python
from improved_statistics import compare_samples_hierarchical

comparison = compare_samples_hierarchical(
    stats1, stats2, 
    "Master", "150°C"
)

# Shows: difference, proper p-value, effect size
```

## What Changes in Your Workflow?

### Minimal! Just add row-group tracking:

**Before:**
```python
all_angles = []
for groove in grooves:
    angle = extract_blaze_angle(...)
    all_angles.append(angle)

mean = np.mean(all_angles)
sem = np.std(all_angles) / np.sqrt(len(all_angles))  # WRONG!
```

**After:**
```python
all_angles = []
all_groups = []
all_errors = []

for group_id in range(n_groups):
    for groove in grooves_in_group:
        angle, _, _, quality = extract_blaze_angle(...)
        all_angles.append(angle)
        all_groups.append(group_id)
        all_errors.append(quality['blaze_angle_stderr'])

stats = calculate_hierarchical_statistics(
    all_angles, all_groups, all_errors
)
# Now stats['sem_best'] is CORRECT!
```

## Diagnostic Plots to Add

I can also create visualization functions for:

1. **Variance component pie chart**
   - Shows measurement/physical/spatial contributions

2. **Row-group consistency plot**
   - Mean ± SE for each row-group
   - Shows spatial variation across image

3. **Q-Q plot for normality**
   - Validates t-test assumptions

4. **Measurement uncertainty vs residuals**
   - Tests if fit uncertainties are realistic

Would you like me to implement these?

## Bottom Line

**Your intuition is right**: You do get more information from 100 measurements than 20.

**The statistics correction**: They're not 100 *independent* measurements, more like ~60 effective ones.

**The good news**: This is still great! And now your p-values will be trustworthy.

**The implementation**: Just track which row-group each measurement comes from, and use the new statistical functions. That's it!
