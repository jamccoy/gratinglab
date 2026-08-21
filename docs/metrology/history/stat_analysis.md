# AFM Statistical Analysis - Current Implementation Review

## EXECUTIVE SUMMARY

After analyzing your code, I've identified **3 major statistical issues** that need addressing:

1. **Row-averaging is NOT tracked** - Pixel averaging reduces noise but this isn't propagated
2. **Error bars show different things** - Bar charts show σ (spread), histograms show distribution
3. **Multiple scan combination may double-count uncertainty** - When combining scans, both within-scan and between-scan variation contribute

---

## DETAILED STATISTICAL FLOW ANALYSIS

### Level 0: Raw AFM Data → Row Average (MISSING ERROR TRACKING)

**Location**: `analyzer.py` → `raw_data()` function (not shown, but inferred)

**What happens:**
```python
# Inferred from your code structure:
raw_y = np.mean(afm_data, axis=1)  # Average each row of pixels
```

**Statistical issue:**
- Each row has ~512 pixels that get averaged
- This averaging REDUCES noise by factor of √512 ≈ 22.6×
- **BUT this uncertainty reduction is NOT tracked anywhere**

**Impact:**
- Your downstream uncertainties (fit errors, angle uncertainties) are too small
- They only reflect row-average noise, not original pixel noise

**Recommendation:**
```python
def raw_data_with_uncertainty(data, scan_x_size):
    """Extract profile with proper uncertainty tracking"""
    # Average across rows
    raw_y = np.mean(data, axis=1)
    
    # Calculate standard error of the mean for each row
    raw_y_stderr = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])
    
    # Create x-axis
    raw_x = np.linspace(0, scan_x_size, len(raw_y))
    
    return raw_x, raw_y, raw_y_stderr
```

---

### Level 1: Linear Fit to Facet (PARTIAL ERROR TRACKING)

**Location**: `analysis.py` → `extract_blaze_angle()` lines 88-119

**What happens:**
```python
blaze_coeffs = np.polyfit(blaze_x_trim, blaze_y_trim, 1)
blaze_slope = blaze_coeffs[0]
blaze_angle = abs(np.arctan(blaze_slope) * 180 / np.pi)
```

**Statistical issue:**
- `np.polyfit` can return covariance matrix but you're not using it
- Fit uncertainty exists but is NOT calculated
- R² is calculated (good!) but doesn't give you angle uncertainty

**Current quality metrics:**
```python
quality_metrics = {
    'blaze_r2': blaze_r2,  # ✓ Good
    'steep_r2': steep_r2,
    'blaze_width_nm': blaze_width,
    'steep_width_nm': steep_width,
    'groove_depth_nm': groove_depth,
    # MISSING: angle_uncertainty from fit
    # MISSING: propagated pixel uncertainty
}
```

**Recommendation:**
```python
# Get fit with full covariance
blaze_coeffs, blaze_cov = np.polyfit(blaze_x_trim, blaze_y_trim, 1, cov=True)
blaze_slope = blaze_coeffs[0]
blaze_slope_stderr = np.sqrt(blaze_cov[0, 0])  # Standard error of slope

# Calculate angle with error propagation
blaze_angle = abs(np.arctan(blaze_slope) * 180 / np.pi)

# Error propagation: dθ/dm = 1/(1+m²) for θ = arctan(m)
blaze_angle_stderr = abs(blaze_slope_stderr / (1 + blaze_slope**2)) * 180 / np.pi

quality_metrics['angle_uncertainty_deg'] = blaze_angle_stderr
```

---

### Level 2: Per-Groove Statistics (CURRENTLY JUST A POINT ESTIMATE)

**Location**: `analyzer.py` → `_extract_all_angles()` lines 158-179

**What happens:**
```python
for i, (center, local_period) in enumerate(zip(groove_centers, local_periods)):
    angle, steep, slope, qual = extract_blaze_angle(...)
    
    if angle is not None:
        blaze_angles.append(angle)  # ← Just the point estimate
        quality.append(qual)
```

**Statistical issue:**
- Each groove has a fit uncertainty (from Level 1) but it's stored in `qual` and never used
- `blaze_angles` is just a list of numbers with no associated uncertainties
- Local angles (within-facet variation) are tracked but separate from the main angle

**Current structure:**
```
blaze_angles = [17.1, 17.3, 16.9, ...]  # Just numbers
quality = [
    {'blaze_r2': 0.998, 'local_angles': [...], ...},
    {'blaze_r2': 0.997, 'local_angles': [...], ...},
    ...
]
```

**Recommendation:**
```python
# Store angles WITH uncertainties
groove_measurements = []
for i, (center, local_period) in enumerate(zip(groove_centers, local_periods)):
    angle, steep, slope, qual = extract_blaze_angle(...)
    
    if angle is not None:
        groove_measurements.append({
            'angle': angle,
            'angle_uncertainty': qual.get('angle_uncertainty_deg', 0),
            'r_squared': qual['blaze_r2'],
            'local_angles': qual.get('local_angles', []),
            'groove_id': i
        })

# Then extract arrays
blaze_angles = [g['angle'] for g in groove_measurements]
angle_uncertainties = [g['angle_uncertainty'] for g in groove_measurements]
```

---

### Level 3: Per-Scan Statistics (CORRECT BUT INCOMPLETE)

**Location**: `analyzer.py` → `_calculate_statistics()` lines 182-202

**What happens:**
```python
mean_angle = np.mean(blaze_angles)
std_angle = np.std(blaze_angles)
```

**This is CORRECT for groove-to-groove variation!**

**But there are TWO sources of uncertainty:**

1. **Measurement uncertainty**: How well can we measure each groove?
   - Comes from fit uncertainty (Level 1)
   - Average per-groove uncertainty ≈ 0.02° (estimated from your R² values)

2. **Physical variation**: How much do grooves actually differ?
   - This is your `std_angle` ≈ 0.18°
   - This is REAL variation, not measurement error

**Current calculation ONLY gives you #2 (physical variation)**

**Standard error calculation:**
```python
# In statistics.py, lines 35-36
se_combined = np.sqrt(r1['std_angle']**2 / r1['n_grooves'] + 
                     r2['std_angle']**2 / r2['n_grooves'])
```

This is **correct** but it assumes:
- Measurement uncertainty is negligible compared to physical variation
- OR it's already included in `std_angle`

**BUT it's NOT included because std_angle comes from variation of point estimates!**

**Proper calculation should be:**
```python
# Separate the two sources
measurement_variance = np.mean([unc**2 for unc in angle_uncertainties])
physical_variance = np.var(blaze_angles, ddof=1)

# Total variance is the sum
total_variance = physical_variance + measurement_variance

# Standard deviation
total_std = np.sqrt(total_variance)

# Standard error of the mean
sem = total_std / np.sqrt(n_grooves)
```

---

### Level 4: Within-Facet (Local) Angles (TRACKED BUT SEPARATE)

**Location**: `analysis.py` → `extract_blaze_angle()` lines 193-215

**What happens:**
```python
# Local angles calculated by sliding window
window_size = max(5, len(blaze_x_trim) // 10)
for i in range(window_size, len(blaze_x_trim) - window_size):
    x_window = blaze_x_trim[i-window_size:i+window_size]
    y_window = blaze_y_trim[i-window_size:i+window_size]
    local_fit = np.polyfit(x_window, y_window, 1)
    local_slopes.append(local_fit[0])

local_angles = np.abs(np.arctan(np.array(local_slopes)) * 180 / np.pi)

quality_metrics['local_angles'] = local_angles
quality_metrics['angle_std'] = np.std(local_angles)
```

**This gives you within-facet curvature/variation!**

**Statistical interpretation:**
- This measures how much the angle varies WITHIN one facet
- It's different from groove-to-groove variation
- Average `angle_std` ≈ 0.15° (from your typical results)

**Usage in visualization:**
- Histograms can show these (statistics.py, line 267-270)
- This gives you MANY more data points (hundreds vs tens)

**Issue: Two different things called "std"**

In your histograms (statistics.py, lines 290-294):
```python
# Title includes:
title += f"{mean_val:.2f}° ± {std_val:.2f}° "
# Where std_val = np.std(all_local_angles)
```

This shows the **spread of local angle measurements**, not uncertainty in the mean!

---

### Level 5: Multiple Scan Combination (CORRECT METHOD)

**Location**: `aggregation.py` → `combine_scans()` lines 10-77

**What happens:**
```python
# Collect ALL angles from ALL scans
all_angles = []
for scan in scan_results:
    all_angles.extend(scan['all_angles'])

# Calculate combined statistics
mean_angle = np.mean(all_angles)
std_angle = np.std(all_angles)
```

**This is a POOLING approach - you're treating all grooves as one big sample.**

**Is this correct?**

**YES, if:** The scans are just different spatial regions of the same sample
**MAYBE NOT, if:** The scans were taken at different times/conditions

**What about scan-to-scan variation?**

You're NOT explicitly calculating it, but it's included in the pooled `std_angle`.

**Alternative approach (more robust):**
```python
def combine_scans_with_variance_decomposition(scan_results):
    """
    Properly separate within-scan and between-scan variation
    """
    # Get mean angle from each scan
    scan_means = [s['mean_angle'] for s in scan_results]
    scan_stds = [s['std_angle'] for s in scan_results]
    scan_n = [s['n_grooves'] for s in scan_results]
    
    # Overall mean (weighted by sample size)
    total_grooves = sum(scan_n)
    overall_mean = sum(m * n for m, n in zip(scan_means, scan_n)) / total_grooves
    
    # Between-scan variance (how much do scan means differ?)
    between_variance = np.var(scan_means, ddof=1) if len(scan_means) > 1 else 0
    
    # Within-scan variance (average of individual scan variances)
    within_variance = sum(s**2 * (n-1) for s, n in zip(scan_stds, scan_n)) / (total_grooves - len(scan_results))
    
    # Total variance
    total_variance = between_variance + within_variance
    total_std = np.sqrt(total_variance)
    
    # Standard error (accounts for all sources)
    sem = total_std / np.sqrt(total_grooves)
    
    return {
        'mean_angle': overall_mean,
        'std_angle': total_std,
        'sem': sem,
        'between_scan_std': np.sqrt(between_variance),
        'within_scan_std': np.sqrt(within_variance),
        'n_scans': len(scan_results),
        'n_grooves': total_grooves,
        # ... other fields
    }
```

**Your current method:**
- Pools all grooves: `std_angle` = total variation (within + between)
- This is actually FINE for most purposes
- But you're not explicitly separating the components

---

## VISUALIZATION ISSUES

### Issue 1: Bar Chart Error Bars (statistics.py, lines 152-156)

```python
bars = ax_bar.bar(x_pos, means, yerr=stds, ...)
```

**What this shows:**
- `yerr=stds` means error bars are ± one standard deviation
- This is the **spread of groove angles**, not uncertainty in the mean

**What users might think:**
- "These are error bars on the mean"
- "Overlapping error bars means not significantly different"

**But actually:**
- Error bars represent groove-to-groove variation
- To test significance, you need **standard error**, not standard deviation

**The fix:**
```python
# Calculate standard errors
sems = [r['std_angle'] / np.sqrt(r['n_grooves']) for r in results]

# Plot with SEM
bars = ax_bar.bar(x_pos, means, yerr=sems, ...)

# Update title
ax_bar.set_title('Mean Blaze Angle Comparison\n(Error bars: ±1 SEM)', ...)
```

**OR keep std but clarify:**
```python
ax_bar.set_title('Mean Blaze Angle Comparison\n(Error bars: ±1σ groove variation, NOT uncertainty in mean)', ...)
```

### Issue 2: Histogram Uses Local Angles (statistics.py, lines 265-271)

```python
# Use all local angles (within-facet measurements) if available
if 'all_local_angles' in r and len(r['all_local_angles']) > 0:
    angles = r['all_local_angles']
    angle_type = "local"
else:
    angles = r['all_angles']
    angle_type = "per-groove"
```

**This is mixing two different things!**

- **Local angles**: Within-facet variation (curvature/camber)
- **Per-groove angles**: Groove-to-groove variation

**In the comparison plot:**
- Bar chart shows per-groove means with per-groove variation
- Histograms show local angles with within-facet variation

**These are NOT the same statistic!**

**Recommendation:**

Either:
1. Use per-groove angles for BOTH (consistent but limited data)
2. Use local angles for BOTH (more data but different meaning)
3. Show BOTH but make it clear they're different

**Best option: #3**
```python
# Create two histogram subplots per sample
# Top: Per-groove distribution (groove-to-groove variation)
# Bottom: Local angle distribution (within-facet variation)
```

---

## COMPARISON STATISTICS (analysis.py from stats folder)

### Pairwise Comparisons (lines 35-48)

```python
diff = r2['mean_angle'] - r1['mean_angle']

se_combined = np.sqrt(r1['std_angle']**2 / r1['n_grooves'] + 
                     r2['std_angle']**2 / r2['n_grooves'])
```

**This is the standard error of the difference in means.**

**Formula is correct!**

**BUT** it assumes:
- Independent samples ✓
- Normal distribution (reasonable with n>10)
- Equal measurement uncertainty (NOT explicitly tracked)

**Issue:** Should this be standard ERROR or standard DEVIATION?

Currently you're reporting: `Difference: +0.15° (±0.057° SE)`

This means:
- The difference is 0.15°
- The uncertainty in this difference is 0.057°
- 95% CI would be roughly ±0.11° (2×SEM)

**This is CORRECT usage of SEM!**

But you should consider adding:
```python
# Calculate t-statistic for significance test
t_stat = diff / se_combined
df = r1['n_grooves'] + r2['n_grooves'] - 2  # Approximate df
from scipy import stats
p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df))

print(f"  t-statistic: {t_stat:.2f}")
print(f"  p-value: {p_value:.4f}")
if p_value < 0.05:
    print(f"  Difference is statistically significant (p < 0.05)")
```

---

## SUMMARY OF ISSUES

### Critical Issues (Fix These!)

1. **Missing: Pixel-level uncertainty propagation**
   - Location: `raw_data()` function
   - Impact: All downstream uncertainties underestimated by √512 ≈ 23×
   - Fix: Track row-averaging uncertainty

2. **Missing: Fit uncertainty in angle calculation**
   - Location: `extract_blaze_angle()`, line 88
   - Impact: No per-groove measurement uncertainty
   - Fix: Use `np.polyfit(..., cov=True)` and propagate

3. **Unclear: Bar chart error bars**
   - Location: `plot_multi_file_comparison()`, line 152
   - Impact: Users may misinterpret as SEM when it's σ
   - Fix: Either change to SEM or clarify label

### Important Improvements

4. **Inconsistent: Histogram shows local angles, bar chart shows per-groove**
   - Location: `plot_multi_file_comparison()`, lines 265-271
   - Impact: Comparing different statistics
   - Fix: Make consistent or show both explicitly

5. **Missing: Variance decomposition for multiple scans**
   - Location: `combine_scans()`, line 33
   - Impact: Can't separate within-scan vs between-scan variation
   - Fix: Calculate and report both components

6. **Missing: Significance tests**
   - Location: `print_pairwise_comparisons()`, line 42
   - Impact: No way to know if differences are significant
   - Fix: Add t-tests and p-values

---

## RECOMMENDED STATISTICAL WORKFLOW

### Complete Uncertainty Tracking

```python
# Level 0: Row averaging
raw_x, raw_y, raw_y_stderr = raw_data_with_uncertainty(afm_data, scan_x_size)

# Level 1: Facet fitting with uncertainty
angle, angle_unc, quality = fit_facet_with_uncertainty(raw_x, raw_y, raw_y_stderr)

# Level 2: Per-scan statistics with separated variance
scan_stats = {
    'mean_angle': mean(angles),
    'measurement_variance': mean(angle_unc**2),  # From fits
    'physical_variance': var(angles),  # Groove-to-groove
    'total_variance': measurement_var + physical_var,
    'std_angle': sqrt(total_variance),
    'sem': sqrt(total_variance / n_grooves)
}

# Level 3: Multiple scan combination
combined = combine_with_variance_decomposition(scan_list)
```

### Visualization Best Practices

```python
# Bar chart: Show SEM for mean comparison
bar(means, yerr=sems, label='±1 SEM (uncertainty in mean)')

# Histogram: Show per-groove distribution
hist(per_groove_angles, label=f'σ={std_grooves:.2f}° (groove variation)')

# Secondary histogram: Show local angle distribution
hist(local_angles, label=f'σ={std_local:.2f}° (within-facet variation)')

# Make clear distinction in titles!
```

---

## NEXT STEPS

### Priority 1: Add Uncertainty Tracking
1. Modify `raw_data()` to return uncertainty
2. Modify `extract_blaze_angle()` to calculate fit uncertainty
3. Store uncertainties alongside angles in results

### Priority 2: Fix Visualizations
4. Clarify what bar chart error bars represent
5. Make histogram and bar chart use consistent statistics
6. Add labels explaining the difference between σ and SEM

### Priority 3: Add Statistical Tests
7. Implement t-tests for pairwise comparisons
8. Add p-values to comparison output
9. Calculate proper confidence intervals

### Priority 4: Decompose Variance (Optional)
10. Separate within-scan and between-scan variation
11. Report both components in output
12. Use for more robust multi-scan statistics

---

## QUESTIONS FOR YOU

1. **What is your primary goal?**
   - Comparing different samples?
   - Measuring absolute angles precisely?
   - Understanding temperature effects?

2. **What error bars do you want on bar charts?**
   - Standard deviation (shows spread of grooves)
   - Standard error (shows uncertainty in mean)
   - 95% confidence interval (shows statistical uncertainty)

3. **How should histograms work?**
   - Show per-groove angles only?
   - Show local angles only?
   - Show both side-by-side?

4. **Do you care about pixel-level uncertainty?**
   - If yes, we need to track it through the whole pipeline
   - If no, we can assume it's negligible (but it's probably not!)

5. **For publications, do you need:**
   - Significance tests (p-values)?
   - Confidence intervals?
   - Formal statistical reporting?

Let me know your priorities and I can create:
- Modified code with proper uncertainty tracking
- Improved visualization with clear labels
- Statistical test implementations
- Documentation of what each error bar means