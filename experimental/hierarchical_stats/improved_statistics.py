"""
Improved Statistical Analysis for AFM Blaze Angle Measurements
Properly handles hierarchical data structure (row-groups within images)

Key improvements:
1. Variance component decomposition
2. Proper handling of correlated measurements
3. Conservative and liberal uncertainty estimates
4. Diagnostic statistics
"""
import numpy as np
from scipy import stats as scipy_stats
import warnings


def calculate_hierarchical_statistics(angles, row_group_labels=None, measurement_errors=None):
    """
    Calculate statistics accounting for hierarchical data structure.
    
    Parameters:
        angles: array of blaze angle measurements
        row_group_labels: array indicating which row-group each measurement belongs to
                         (None if no row-grouping)
        measurement_errors: array of measurement uncertainties (stderr from fits)
    
    Returns:
        stats_dict: Dictionary with comprehensive statistics
    """
    angles = np.array(angles)
    n_total = len(angles)
    
    if n_total == 0:
        return None
    
    # Basic statistics (assuming independence)
    mean_angle = np.mean(angles)
    std_angle = np.std(angles, ddof=1)
    
    stats_dict = {
        'mean_angle': mean_angle,
        'std_angle': std_angle,
        'n_measurements': n_total,
    }
    
    # Measurement uncertainty component
    if measurement_errors is not None:
        measurement_errors = np.array(measurement_errors)
        mean_meas_error = np.mean(measurement_errors)
        measurement_variance = np.mean(measurement_errors**2)  # Average of variances
        
        stats_dict['mean_measurement_error'] = mean_meas_error
        stats_dict['measurement_variance'] = measurement_variance
    else:
        measurement_variance = 0
        stats_dict['mean_measurement_error'] = None
        stats_dict['measurement_variance'] = None
    
    # If we have row-group information, do hierarchical analysis
    if row_group_labels is not None:
        row_group_labels = np.array(row_group_labels)
        unique_groups = np.unique(row_group_labels)
        n_groups = len(unique_groups)
        
        # Calculate mean for each group
        group_means = []
        group_sizes = []
        group_stds = []
        
        for group in unique_groups:
            mask = row_group_labels == group
            group_angles = angles[mask]
            group_means.append(np.mean(group_angles))
            group_sizes.append(len(group_angles))
            if len(group_angles) > 1:
                group_stds.append(np.std(group_angles, ddof=1))
            else:
                group_stds.append(0)
        
        group_means = np.array(group_means)
        group_sizes = np.array(group_sizes)
        group_stds = np.array(group_stds)
        
        # Between-group variance (how much do group means vary?)
        between_group_variance = np.var(group_means, ddof=1)
        
        # Within-group variance (average variance within groups)
        # Weighted by group size
        if n_groups > 1 and np.sum(group_sizes > 1) > 0:
            within_group_variance = np.average(
                group_stds**2, 
                weights=group_sizes
            )
        else:
            within_group_variance = 0
        
        # Physical variance (remove measurement component)
        physical_variance = max(0, within_group_variance - measurement_variance)
        
        # Total variance decomposition
        # Var(total) = Var(between groups) + Var(within groups) + Var(measurement)
        total_variance_hierarchical = (
            between_group_variance + 
            physical_variance + 
            measurement_variance
        )
        
        # Standard errors
        # CONSERVATIVE: Use group means as independent units
        sem_conservative = np.std(group_means, ddof=1) / np.sqrt(n_groups)
        
        # LIBERAL: Assume all measurements independent (ignore hierarchy)
        sem_liberal = std_angle / np.sqrt(n_total)
        
        # BEST ESTIMATE: Account for correlation within groups
        # Effective sample size depends on intraclass correlation
        if within_group_variance > 0:
            icc = between_group_variance / (between_group_variance + within_group_variance)
        else:
            icc = 0
        
        # Design effect: how much correlation reduces effective sample size
        avg_group_size = n_total / n_groups
        design_effect = 1 + (avg_group_size - 1) * icc
        n_effective = n_total / design_effect
        sem_best = std_angle / np.sqrt(n_effective)
        
        # Add hierarchical statistics
        stats_dict.update({
            'n_groups': n_groups,
            'group_means': group_means,
            'group_sizes': group_sizes,
            'group_stds': group_stds,
            'between_group_variance': between_group_variance,
            'within_group_variance': within_group_variance,
            'physical_variance': physical_variance,
            'total_variance_hierarchical': total_variance_hierarchical,
            'intraclass_correlation': icc,
            'design_effect': design_effect,
            'n_effective': n_effective,
            'sem_conservative': sem_conservative,
            'sem_liberal': sem_liberal,
            'sem_best': sem_best,
            'use_hierarchical': True
        })
        
    else:
        # No row-group information: use simple statistics
        sem_simple = std_angle / np.sqrt(n_total)
        
        stats_dict.update({
            'n_groups': None,
            'sem_conservative': sem_simple,
            'sem_liberal': sem_simple,
            'sem_best': sem_simple,
            'use_hierarchical': False
        })
    
    return stats_dict


def compare_samples_hierarchical(stats1, stats2, label1="Sample 1", label2="Sample 2"):
    """
    Compare two samples with proper hierarchical statistics.
    
    Returns:
        comparison_dict: Dictionary with comparison statistics
    """
    mean1 = stats1['mean_angle']
    mean2 = stats2['mean_angle']
    diff = mean2 - mean1
    
    # Use best SEM estimates
    sem1 = stats1['sem_best']
    sem2 = stats2['sem_best']
    
    # Combined standard error of the difference
    se_diff = np.sqrt(sem1**2 + sem2**2)
    
    # Degrees of freedom
    # If hierarchical, use group-level df
    if stats1['use_hierarchical'] and stats2['use_hierarchical']:
        # Welch-Satterthwaite with group-level statistics
        var1 = stats1['between_group_variance'] + stats1['within_group_variance']/stats1['n_groups']
        var2 = stats2['between_group_variance'] + stats2['within_group_variance']/stats2['n_groups']
        n1 = stats1['n_groups']
        n2 = stats2['n_groups']
        
        df = welch_satterthwaite_df(var1, n1, var2, n2)
    else:
        # Use effective sample sizes
        n1_eff = stats1.get('n_effective', stats1['n_measurements'])
        n2_eff = stats2.get('n_effective', stats2['n_measurements'])
        var1 = stats1['std_angle']**2
        var2 = stats2['std_angle']**2
        
        df = welch_satterthwaite_df(var1, n1_eff, var2, n2_eff)
    
    # t-statistic and p-value
    t_stat = diff / se_diff if se_diff > 0 else 0
    p_value = 2 * (1 - scipy_stats.t.cdf(abs(t_stat), df))
    
    # Confidence interval
    t_critical = scipy_stats.t.ppf(0.975, df)
    ci_lower = diff - t_critical * se_diff
    ci_upper = diff + t_critical * se_diff
    
    # Effect size (Cohen's d)
    # Use total SD for effect size calculation
    pooled_std = np.sqrt(
        ((stats1['n_measurements']-1) * stats1['std_angle']**2 + 
         (stats2['n_measurements']-1) * stats2['std_angle']**2) / 
        (stats1['n_measurements'] + stats2['n_measurements'] - 2)
    )
    cohens_d = diff / pooled_std if pooled_std > 0 else 0
    
    comparison = {
        'label1': label1,
        'label2': label2,
        'difference': diff,
        'se_difference': se_diff,
        'ci_95_lower': ci_lower,
        'ci_95_upper': ci_upper,
        't_statistic': t_stat,
        'df': df,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'significant_05': p_value < 0.05,
        'significant_01': p_value < 0.01,
        'significant_001': p_value < 0.001
    }
    
    return comparison


def welch_satterthwaite_df(var1, n1, var2, n2):
    """Calculate Welch-Satterthwaite degrees of freedom."""
    if n1 < 2 or n2 < 2:
        return min(n1, n2) - 1
    
    numerator = (var1/n1 + var2/n2)**2
    denominator = (var1/n1)**2/(n1-1) + (var2/n2)**2/(n2-1)
    
    if denominator > 0:
        df = numerator / denominator
    else:
        df = min(n1, n2) - 1
    
    return max(1, df)


def print_hierarchical_statistics(stats, label="Sample"):
    """
    Print hierarchical statistics in a clear, interpretable format.
    """
    print(f"\n{'='*80}")
    print(f"HIERARCHICAL STATISTICS: {label}")
    print(f"{'='*80}")
    
    print(f"\nMean angle: {stats['mean_angle']:.3f}°")
    print(f"Total measurements: N = {stats['n_measurements']}")
    
    if stats['use_hierarchical']:
        print(f"Row groups: {stats['n_groups']}")
        print(f"Avg measurements per group: {stats['n_measurements']/stats['n_groups']:.1f}")
        
        print(f"\n--- VARIANCE DECOMPOSITION ---")
        total_var = stats['std_angle']**2
        
        if stats['measurement_variance'] is not None:
            meas_pct = 100 * stats['measurement_variance'] / total_var if total_var > 0 else 0
            print(f"Measurement error:     {np.sqrt(stats['measurement_variance']):.3f}° "
                  f"({meas_pct:.1f}% of variance)")
        
        phys_pct = 100 * stats['physical_variance'] / total_var if total_var > 0 else 0
        between_pct = 100 * stats['between_group_variance'] / total_var if total_var > 0 else 0
        
        print(f"Within-group (physical): {np.sqrt(stats['physical_variance']):.3f}° "
              f"({phys_pct:.1f}% of variance)")
        print(f"Between-group (spatial): {np.sqrt(stats['between_group_variance']):.3f}° "
              f"({between_pct:.1f}% of variance)")
        print(f"Total SD:                {stats['std_angle']:.3f}°")
        
        print(f"\n--- CORRELATION & EFFECTIVE SAMPLE SIZE ---")
        print(f"Intraclass correlation (ICC): {stats['intraclass_correlation']:.3f}")
        print(f"  (How correlated are measurements within same row-group?)")
        print(f"  ICC=0: No correlation (independent)")
        print(f"  ICC=1: Perfect correlation (all info in group means)")
        
        print(f"Design effect: {stats['design_effect']:.2f}")
        print(f"  (How much does clustering reduce effective sample size)")
        
        print(f"Effective sample size: {stats['n_effective']:.1f}")
        print(f"  (Equivalent number of independent measurements)")
        
        print(f"\n--- STANDARD ERROR ESTIMATES ---")
        print(f"Conservative (group means only): {stats['sem_conservative']:.4f}°")
        print(f"Best estimate (ICC-adjusted):    {stats['sem_best']:.4f}°")
        print(f"Liberal (ignore clustering):     {stats['sem_liberal']:.4f}°")
        print(f"\nRecommended: Use 'Best estimate' for comparisons")
        
        print(f"\n--- CONFIDENCE INTERVALS (95%) ---")
        ci_cons = 1.96 * stats['sem_conservative']
        ci_best = 1.96 * stats['sem_best']
        ci_lib = 1.96 * stats['sem_liberal']
        
        print(f"Conservative: {stats['mean_angle']:.3f}° ± {ci_cons:.3f}°  "
              f"[{stats['mean_angle']-ci_cons:.3f}°, {stats['mean_angle']+ci_cons:.3f}°]")
        print(f"Best estimate: {stats['mean_angle']:.3f}° ± {ci_best:.3f}°  "
              f"[{stats['mean_angle']-ci_best:.3f}°, {stats['mean_angle']+ci_best:.3f}°]")
        print(f"Liberal:      {stats['mean_angle']:.3f}° ± {ci_lib:.3f}°  "
              f"[{stats['mean_angle']-ci_lib:.3f}°, {stats['mean_angle']+ci_lib:.3f}°]")
        
    else:
        print(f"\nNo hierarchical structure (single-level analysis)")
        print(f"Standard deviation: {stats['std_angle']:.3f}°")
        print(f"Standard error: {stats['sem_best']:.4f}°")
        print(f"95% CI: {stats['mean_angle']:.3f}° ± {1.96*stats['sem_best']:.3f}°")


def print_comparison_summary_hierarchical(comparisons, stats_list, labels):
    """
    Print summary of all comparisons with hierarchical statistics.
    """
    print(f"\n{'='*80}")
    print(f"PAIRWISE COMPARISONS (with hierarchical statistics)")
    print(f"{'='*80}\n")
    
    for comp in comparisons:
        print(f"{comp['label1']} vs {comp['label2']}:")
        print(f"  Difference: {comp['difference']:+.3f}° ± {comp['se_difference']:.3f}° (SE)")
        print(f"  95% CI: [{comp['ci_95_lower']:+.3f}°, {comp['ci_95_upper']:+.3f}°]")
        print(f"  t({comp['df']:.1f}) = {comp['t_statistic']:.2f}, p = {comp['p_value']:.4f}", end="")
        
        if comp['significant_001']:
            print(" ***  (highly significant)")
        elif comp['significant_01']:
            print(" **   (very significant)")
        elif comp['significant_05']:
            print(" *    (significant)")
        else:
            print("      (not significant)")
        
        print(f"  Effect size: Cohen's d = {comp['cohens_d']:.2f}", end="")
        if abs(comp['cohens_d']) < 0.2:
            print(" (negligible)")
        elif abs(comp['cohens_d']) < 0.5:
            print(" (small)")
        elif abs(comp['cohens_d']) < 0.8:
            print(" (medium)")
        else:
            print(" (large)")
        print()


def test_normality(angles, label="Sample"):
    """
    Test if measurements are normally distributed.
    Important for validating t-test assumptions.
    """
    from scipy import stats
    
    angles = np.array(angles)
    n = len(angles)
    
    print(f"\n{'='*60}")
    print(f"NORMALITY TESTS: {label}")
    print(f"{'='*60}")
    
    if n < 3:
        print("Not enough data for normality tests (N < 3)")
        return
    
    # Shapiro-Wilk test (good for small-medium samples)
    if n <= 5000:
        stat_sw, p_sw = stats.shapiro(angles)
        print(f"Shapiro-Wilk test: W = {stat_sw:.4f}, p = {p_sw:.4f}")
        if p_sw < 0.05:
            print("  → Data deviates from normality (p < 0.05)")
        else:
            print("  → Data consistent with normality")
    
    # Anderson-Darling test
    result_ad = stats.anderson(angles, dist='norm')
    print(f"\nAnderson-Darling test: A² = {result_ad.statistic:.4f}")
    print(f"  Critical values: {result_ad.critical_values}")
    print(f"  Significance levels: {result_ad.significance_level}%")
    
    # D'Agostino-Pearson test (combines skewness and kurtosis)
    if n >= 8:
        stat_dp, p_dp = stats.normaltest(angles)
        print(f"\nD'Agostino-Pearson test: χ² = {stat_dp:.4f}, p = {p_dp:.4f}")
        if p_dp < 0.05:
            print("  → Data deviates from normality")
        else:
            print("  → Data consistent with normality")
    
    # Descriptive statistics
    skewness = stats.skew(angles)
    kurtosis = stats.kurtosis(angles)
    print(f"\nSkewness: {skewness:.3f} (0 = symmetric)")
    print(f"Kurtosis: {kurtosis:.3f} (0 = normal tails)")
    
    print(f"\nNote: For small samples (N < 30), normality tests have low power.")
    print(f"      Visual inspection (Q-Q plot) is often more informative.")


# Example usage
if __name__ == "__main__":
    # Simulate some hierarchical data
    np.random.seed(42)
    
    # 20 row-groups, 5 measurements each
    n_groups = 20
    n_per_group = 5
    
    group_means = np.random.normal(17.5, 0.3, n_groups)  # Group-to-group variation
    
    angles = []
    groups = []
    errors = []
    
    for i, gm in enumerate(group_means):
        # Within-group variation
        group_angles = np.random.normal(gm, 0.15, n_per_group)
        angles.extend(group_angles)
        groups.extend([i] * n_per_group)
        # Measurement errors
        errors.extend(np.random.uniform(0.05, 0.15, n_per_group))
    
    angles = np.array(angles)
    groups = np.array(groups)
    errors = np.array(errors)
    
    # Calculate statistics
    stats_result = calculate_hierarchical_statistics(angles, groups, errors)
    print_hierarchical_statistics(stats_result, "Example Sample")
