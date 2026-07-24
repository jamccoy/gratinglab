# Integration Guide: Adding Hierarchical Statistics to Your AFM Analysis

## Overview

This guide shows you exactly how to modify your existing code to use the improved statistical analysis.

## Step 1: Copy the New Modules

Place these files in your project:
- `improved_statistics.py` → in your `stats/` directory (or rename existing `analysis.py`)
- `diagnostic_plots.py` → in your `visualization/` directory

## Step 2: Modify Your Analyzer

### In `analyzer.py` (or wherever you collect measurements):

**BEFORE:**
```python
def analyze_single_file(filename, ...):
    # ... load data ...
    
    all_angles = []
    quality_list = []
    
    for groove_center in groove_centers:
        angle, steep, slope, quality = extract_blaze_angle(...)
        if angle is not None:
            all_angles.append(angle)
            quality_list.append(quality)
    
    # Simple statistics
    result = {
        'filename': filename,
        'mean_angle': np.mean(all_angles),
        'std_angle': np.std(all_angles, ddof=1),
        'n_grooves': len(all_angles),
        # ...
    }
    
    return result
```

**AFTER:**
```python
from stats.improved_statistics import calculate_hierarchical_statistics

def analyze_single_file(filename, use_row_groups=True, n_groups=20, ...):
    # ... load data ...
    
    all_angles = []
    all_row_groups = []  # NEW: Track which row-group
    all_meas_errors = []  # NEW: Track measurement uncertainty
    quality_list = []
    
    if use_row_groups:
        # Extract multiple row-group profiles
        disp_um, profiles_nm, group_info = raw_data_multi_group(data, scan_x_size, n_groups)
        
        for group_id, profile in enumerate(profiles_nm):
            # Flatten this group's profile
            x, y_flat, background = flatten_profile(disp_um, profile, ...)
            
            # Find grooves in this group
            groove_centers = find_groove_positions(x, y_flat, ...)
            
            # Extract each groove
            for groove_center in groove_centers:
                angle, steep, slope, quality = extract_blaze_angle(
                    x, y_flat, groove_center, ...
                )
                
                if angle is not None:
                    all_angles.append(angle)
                    all_row_groups.append(group_id)  # NEW: Track group
                    all_meas_errors.append(quality['blaze_angle_stderr'])  # NEW
                    quality_list.append(quality)
    
    else:
        # Original single-profile analysis
        # ... existing code ...
        all_row_groups = None  # No hierarchy
    
    # NEW: Calculate hierarchical statistics
    stats = calculate_hierarchical_statistics(
        angles=all_angles,
        row_group_labels=all_row_groups,
        measurement_errors=all_meas_errors
    )
    
    # Build result dictionary with hierarchical stats
    result = {
        'filename': filename,
        'mean_angle': stats['mean_angle'],
        'std_angle': stats['std_angle'],
        'n_grooves': stats['n_measurements'],
        
        # NEW: Add hierarchical info
        'sem_best': stats['sem_best'],  # Use this for comparisons
        'sem_conservative': stats['sem_conservative'],
        'sem_liberal': stats['sem_liberal'],
        'n_effective': stats.get('n_effective', stats['n_measurements']),
        
        # Variance components
        'measurement_variance': stats.get('measurement_variance'),
        'physical_variance': stats.get('physical_variance'),
        'between_group_variance': stats.get('between_group_variance'),
        'intraclass_correlation': stats.get('intraclass_correlation'),
        
        # Keep existing fields
        'angles': all_angles,
        'row_groups': all_row_groups,  # NEW
        'meas_errors': all_meas_errors,  # NEW
        'quality': quality_list,
        # ... other fields ...
    }
    
    return result
```

## Step 3: Update Comparison Code

### In `stats/analysis.py` (comparison functions):

**BEFORE:**
```python
def print_pairwise_comparisons(results):
    for i, r1 in enumerate(results):
        for j, r2 in enumerate(results):
            if i < j:
                diff = r2['mean_angle'] - r1['mean_angle']
                
                # OLD: Simple SEM
                sem1 = r1['std_angle'] / np.sqrt(r1['n_grooves'])
                sem2 = r2['std_angle'] / np.sqrt(r2['n_grooves'])
                se_combined = np.sqrt(sem1**2 + sem2**2)
                
                # ... t-test ...
```

**AFTER:**
```python
from stats.improved_statistics import compare_samples_hierarchical

def print_pairwise_comparisons_hierarchical(results, labels):
    """Enhanced pairwise comparisons using hierarchical statistics."""
    
    comparisons = []
    
    for i, r1 in enumerate(results):
        for j, r2 in enumerate(results):
            if i < j:
                # NEW: Use hierarchical comparison
                comp = compare_samples_hierarchical(
                    r1['hierarchical_stats'],  # Stored from analyzer
                    r2['hierarchical_stats'],
                    label1=labels[i],
                    label2=labels[j]
                )
                comparisons.append(comp)
    
    # Print results
    print(f"\n{'='*80}")
    print(f"PAIRWISE COMPARISONS (Hierarchical Statistics)")
    print(f"{'='*80}\n")
    
    for comp in comparisons:
        print(f"{comp['label1']} vs {comp['label2']}:")
        print(f"  Difference: {comp['difference']:+.3f}° ± {comp['se_difference']:.3f}°")
        print(f"  95% CI: [{comp['ci_95_lower']:+.3f}°, {comp['ci_95_upper']:+.3f}°]")
        print(f"  t({comp['df']:.1f}) = {comp['t_statistic']:.2f}, p = {comp['p_value']:.4f}", end="")
        
        if comp['significant_001']:
            print(" ***")
        elif comp['significant_01']:
            print(" **")
        elif comp['significant_05']:
            print(" *")
        else:
            print()
        
        print(f"  Cohen's d = {comp['cohens_d']:.2f}\n")
    
    return comparisons
```

## Step 4: Add Diagnostic Plots

### In your workflow (e.g., `workflows.py`):

```python
from visualization.diagnostic_plots import (
    create_diagnostic_report,
    plot_comparison_forest
)

def run_comparison_analysis():
    # ... existing analysis code ...
    
    # NEW: Create diagnostic report for each sample
    for result, label in zip(results, labels):
        if 'row_groups' in result and result['row_groups'] is not None:
            create_diagnostic_report(
                stats_dict=result['hierarchical_stats'],
                angles=result['angles'],
                row_group_labels=result['row_groups'],
                measurement_errors=result['meas_errors'],
                label=label
            )
    
    # NEW: Forest plot of all comparisons
    if len(comparisons) > 0:
        plot_comparison_forest(comparisons)
    
    plt.show()
```

## Step 5: Update Config File

### In `config.py`:

```python
# ============ ROW-GROUP ANALYSIS ============
USE_ROW_GROUPS = True
N_ROW_GROUPS = 20  # Adjust based on image size

# NEW: Statistical options
STATISTICAL_METHOD = 'hierarchical'  # 'simple' or 'hierarchical'
USE_CONSERVATIVE_SEM = False  # True for publications, False for exploratory
SHOW_DIAGNOSTIC_PLOTS = True  # Show variance decomposition, Q-Q plots, etc.

# ============ REPORTING OPTIONS ============
REPORT_SEM_TYPE = 'best'  # 'conservative', 'best', or 'liberal'
SIGNIFICANCE_LEVEL = 0.05
```

## Step 6: Modify Print Functions

### Update your summary printing:

```python
def print_comparison_summary(results, labels):
    print(f"\n{'='*80}")
    print(f"COMPARISON SUMMARY")
    print(f"{'='*80}")
    
    # Header
    print(f"{'Sample':<20} {'N_meas':<7} {'N_eff':<7} {'Mean ± SEM':<20} {'ICC':<6}")
    print(f"{'-'*80}")
    
    for r, label in zip(results, labels):
        n_meas = r['n_grooves']
        n_eff = r.get('n_effective', n_meas)
        mean = r['mean_angle']
        sem = r['sem_best']
        icc = r.get('intraclass_correlation', 0)
        
        print(f"{label:<20} {n_meas:<7} {n_eff:<7.1f} "
              f"{mean:.2f} ± {sem:.3f}°    {icc:<6.3f}")
    
    print(f"\nNote: N_eff = Effective sample size (accounts for clustering)")
    print(f"      ICC = Intraclass correlation (0=independent, 1=fully clustered)")
    print(f"      SEM uses '{REPORT_SEM_TYPE}' estimate")
```

## Example: Complete Workflow Integration

Here's how a complete analysis might look:

```python
# main.py or workflows.py

def run_comparison_analysis_v2():
    """
    Enhanced comparison analysis with hierarchical statistics.
    """
    print("Starting AFM Comparison Analysis (Hierarchical Statistics)")
    print("="*80)
    
    results = []
    labels = []
    
    # Analyze each sample
    for filename, label, temp in SAMPLES_TO_COMPARE:
        print(f"\nAnalyzing: {label}")
        print(f"  File: {filename}")
        
        result = analyze_single_file(
            filename,
            use_row_groups=USE_ROW_GROUPS,
            n_groups=N_ROW_GROUPS,
            # ... other params ...
        )
        
        results.append(result)
        labels.append(label)
        
        # Print individual statistics
        if SHOW_DIAGNOSTIC_PLOTS:
            from stats.improved_statistics import print_hierarchical_statistics
            print_hierarchical_statistics(result['hierarchical_stats'], label)
    
    # Compare all samples
    print("\n" + "="*80)
    print("STATISTICAL COMPARISONS")
    print("="*80)
    
    comparisons = print_pairwise_comparisons_hierarchical(results, labels)
    
    # Temperature analysis (if applicable)
    temperatures = [temp for _, _, temp in SAMPLES_TO_COMPARE]
    if any(t is not None for t in temperatures):
        print_temperature_analysis(results, labels, temperatures)
    
    # Diagnostic plots
    if SHOW_DIAGNOSTIC_PLOTS:
        print("\nGenerating diagnostic plots...")
        
        for result, label in zip(results, labels):
            if result.get('row_groups') is not None:
                create_diagnostic_report(
                    result['hierarchical_stats'],
                    result['angles'],
                    result['row_groups'],
                    result['meas_errors'],
                    label
                )
        
        # Forest plot
        plot_comparison_forest(comparisons)
    
    # Profiles plot
    plot_sample_profiles_by_temperature(...)
    
    plt.show()
    print("\n✓ Analysis complete!")
```

## Migration Checklist

- [ ] Copy `improved_statistics.py` and `diagnostic_plots.py` to project
- [ ] Modify analyzer to track row_groups and measurement_errors
- [ ] Update result dictionary to include hierarchical stats
- [ ] Replace comparison functions with hierarchical versions
- [ ] Add diagnostic plot calls to workflow
- [ ] Update config file with new options
- [ ] Test on one sample first
- [ ] Verify SEM values are reasonable
- [ ] Check diagnostic plots make sense
- [ ] Update all samples
- [ ] Review and interpret variance components

## Testing Your Integration

Run this quick test:

```python
# test_hierarchical_stats.py

import numpy as np
from stats.improved_statistics import calculate_hierarchical_statistics, print_hierarchical_statistics

# Simulate data like yours
angles = np.random.normal(17.5, 0.4, 100)
groups = np.repeat(range(20), 5)  # 20 groups, 5 per group
errors = np.random.uniform(0.05, 0.15, 100)

stats = calculate_hierarchical_statistics(angles, groups, errors)
print_hierarchical_statistics(stats, "Test Sample")

# Check that:
# 1. n_effective < n_measurements (usually)
# 2. sem_conservative > sem_best > sem_liberal
# 3. ICC is between 0 and 1
# 4. Variance components are non-negative

print("\nTest passed!" if stats['n_effective'] < 100 else "WARNING: Check calculations")
```

## Questions or Issues?

Common problems and solutions:

**Q: My ICC is negative or > 1**
A: This can happen with very small samples. The code clamps it to [0,1], but check your data.

**Q: n_effective is larger than n_measurements**
A: This shouldn't happen. Check that row_group_labels are correctly assigned.

**Q: SEM values seem too small**
A: Make sure you're using 'sem_best' or 'sem_conservative', not 'sem_liberal'.

**Q: Diagnostic plots show non-normality**
A: This is OK for moderate deviations. t-tests are robust. For severe non-normality, consider bootstrap CIs.

**Q: Variance components don't sum correctly**
A: They're not supposed to sum! Variances add, standard deviations don't.
