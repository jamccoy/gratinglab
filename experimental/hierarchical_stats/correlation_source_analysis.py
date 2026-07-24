"""
Physical vs Measurement Correlation Analysis

Helps determine: Are grooves within row-groups correlated due to
(A) Measurement/processing artifacts, OR
(B) Real physical spatial correlation?

Key insight: If grooves are physically independent (true manufacturing variation),
then ICC should be LOW even with row-grouping.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats


def analyze_correlation_source(angles, row_group_labels, spatial_positions=None, label="Sample"):
    """
    Distinguish between measurement artifacts and real physical correlation.
    
    Parameters:
        angles: array of blaze angle measurements
        row_group_labels: which row-group each measurement belongs to
        spatial_positions: actual physical position of each groove (if available)
        label: sample name
    
    Returns:
        analysis_dict: Dictionary with diagnostic information
    """
    angles = np.array(angles)
    row_group_labels = np.array(row_group_labels)
    
    unique_groups = np.unique(row_group_labels)
    n_groups = len(unique_groups)
    
    print(f"\n{'='*80}")
    print(f"CORRELATION SOURCE ANALYSIS: {label}")
    print(f"{'='*80}\n")
    
    # Calculate ICC
    group_means = []
    group_sizes = []
    group_vars = []
    
    for g in unique_groups:
        mask = row_group_labels == g
        group_angles = angles[mask]
        group_means.append(np.mean(group_angles))
        group_sizes.append(len(group_angles))
        if len(group_angles) > 1:
            group_vars.append(np.var(group_angles, ddof=1))
        else:
            group_vars.append(0)
    
    group_means = np.array(group_means)
    group_vars = np.array(group_vars)
    
    # Between-group variance
    between_var = np.var(group_means, ddof=1)
    
    # Within-group variance (average)
    within_var = np.mean(group_vars)
    
    # ICC
    total_var = between_var + within_var
    icc = between_var / total_var if total_var > 0 else 0
    
    print(f"VARIANCE BREAKDOWN:")
    print(f"  Between-group variance: {between_var:.6f}° (σ = {np.sqrt(between_var):.3f}°)")
    print(f"  Within-group variance:  {within_var:.6f}° (σ = {np.sqrt(within_var):.3f}°)")
    print(f"  Total variance:         {total_var:.6f}° (σ = {np.sqrt(total_var):.3f}°)")
    print(f"\n  ICC = {icc:.4f}")
    
    # Interpret ICC
    print(f"\nINTERPRETATION:")
    if icc < 0.05:
        print(f"  ✓ ICC < 0.05: Grooves within row-groups are nearly INDEPENDENT")
        print(f"    → Manufacturing variation dominates")
        print(f"    → Effective N ≈ {len(angles)} (nearly all measurements count!)")
        print(f"    → Your row-grouping is NOT creating artificial correlation")
        conclusion = "physical_independent"
    elif icc < 0.15:
        print(f"  ≈ ICC = {icc:.3f}: Weak correlation within row-groups")
        print(f"    → Mostly independent, with slight spatial trends")
        print(f"    → Effective N ≈ {len(angles) * (1-icc):.0f}")
        conclusion = "weak_correlation"
    elif icc < 0.40:
        print(f"  ≈ ICC = {icc:.3f}: Moderate correlation within row-groups")
        print(f"    → Spatial gradients present but substantial local variation")
        print(f"    → Effective N ≈ {len(angles) * (1-icc):.0f}")
        conclusion = "moderate_correlation"
    else:
        print(f"  ! ICC = {icc:.3f}: Strong correlation within row-groups")
        print(f"    → Group means are more important than individual measurements")
        print(f"    → Effective N ≈ {n_groups * (1 + (1-icc)*len(angles)/n_groups):.0f}")
        print(f"    → Consider spatial gradients or measurement artifacts")
        conclusion = "strong_correlation"
    
    # Test for spatial autocorrelation (if positions available)
    if spatial_positions is not None:
        print(f"\n{'='*60}")
        print(f"SPATIAL AUTOCORRELATION TEST")
        print(f"{'='*60}")
        
        spatial_positions = np.array(spatial_positions)
        
        # Sort by position
        sort_idx = np.argsort(spatial_positions)
        sorted_angles = angles[sort_idx]
        sorted_positions = spatial_positions[sort_idx]
        
        # Calculate lag-1 autocorrelation (adjacent grooves)
        if len(sorted_angles) > 1:
            diffs = np.diff(sorted_angles)
            mean_diff = np.mean(sorted_angles)
            
            # Pearson correlation between adjacent pairs
            autocorr_1 = np.corrcoef(sorted_angles[:-1], sorted_angles[1:])[0, 1]
            
            print(f"\nLag-1 autocorrelation (adjacent grooves): {autocorr_1:.4f}")
            
            if abs(autocorr_1) < 0.1:
                print(f"  → Adjacent grooves are nearly UNCORRELATED")
                print(f"  → Supports groove-to-groove independence")
            elif abs(autocorr_1) < 0.3:
                print(f"  → Weak correlation between adjacent grooves")
            else:
                print(f"  → Significant correlation between adjacent grooves")
                print(f"  → May indicate spatial gradients")
    
    # Statistical test: Are group means significantly different?
    print(f"\n{'='*60}")
    print(f"ONE-WAY ANOVA: Do row-groups have different means?")
    print(f"{'='*60}")
    
    # ANOVA
    groups_list = [angles[row_group_labels == g] for g in unique_groups]
    f_stat, p_value = scipy_stats.f_oneway(*groups_list)
    
    print(f"\nF({n_groups-1}, {len(angles)-n_groups}) = {f_stat:.3f}, p = {p_value:.4f}")
    
    if p_value > 0.05:
        print(f"  → Group means are NOT significantly different (p > 0.05)")
        print(f"  → No evidence of spatial variation")
        print(f"  → Treating all grooves as independent is justified!")
    else:
        print(f"  → Group means ARE significantly different (p < 0.05)")
        print(f"  → Spatial variation is real")
        print(f"  → Hierarchical analysis is appropriate")
    
    # Create diagnostic plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Top-left: Group means
    ax = axes[0, 0]
    ax.bar(range(n_groups), group_means, color='steelblue', edgecolor='black')
    overall_mean = np.mean(angles)
    ax.axhline(overall_mean, color='red', linestyle='--', linewidth=2, label='Overall mean')
    ax.set_xlabel('Row Group', fontsize=11)
    ax.set_ylabel('Mean Angle (degrees)', fontsize=11)
    ax.set_title(f'Group Means (ICC = {icc:.3f})', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Top-right: Within-group scatter
    ax = axes[0, 1]
    for g in unique_groups:
        mask = row_group_labels == g
        group_angles = angles[mask]
        x_jitter = np.random.normal(g, 0.1, len(group_angles))
        ax.scatter(x_jitter, group_angles, alpha=0.6, s=30)
    ax.set_xlabel('Row Group', fontsize=11)
    ax.set_ylabel('Angle (degrees)', fontsize=11)
    ax.set_title('Within-Group Scatter', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Bottom-left: Within vs between variance
    ax = axes[1, 0]
    categories = ['Within-Group\n(groove-to-groove)', 'Between-Group\n(spatial)']
    variances = [within_var, between_var]
    colors = ['#66b3ff', '#99ff99']
    ax.bar(categories, variances, color=colors, edgecolor='black', linewidth=2)
    ax.set_ylabel('Variance (degrees²)', fontsize=11)
    ax.set_title('Variance Components', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Bottom-right: Autocorrelation (if spatial data available)
    ax = axes[1, 1]
    if spatial_positions is not None:
        # Plot angles vs position
        sort_idx = np.argsort(spatial_positions)
        ax.plot(spatial_positions[sort_idx], angles[sort_idx], 'o-', 
               markersize=4, alpha=0.6, linewidth=1)
        ax.set_xlabel('Spatial Position (µm)', fontsize=11)
        ax.set_ylabel('Angle (degrees)', fontsize=11)
        ax.set_title(f'Spatial Trend (r_lag1 = {autocorr_1:.3f})', fontsize=12, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'Spatial position\ndata not available', 
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Summary
    analysis = {
        'icc': icc,
        'between_var': between_var,
        'within_var': within_var,
        'f_statistic': f_stat,
        'anova_p_value': p_value,
        'conclusion': conclusion,
        'n_effective': len(angles) / (1 + (len(angles)/n_groups - 1) * icc)
    }
    
    if spatial_positions is not None:
        analysis['autocorr_lag1'] = autocorr_1
    
    return analysis


def recommend_statistical_approach(icc, anova_p_value):
    """
    Based on ICC and ANOVA results, recommend best statistical approach.
    """
    print(f"\n{'='*80}")
    print(f"RECOMMENDED STATISTICAL APPROACH")
    print(f"{'='*80}\n")
    
    if icc < 0.1 and anova_p_value > 0.05:
        print("✓ SIMPLE STATISTICS (treating all grooves as independent)")
        print("\nReason:")
        print(f"  - ICC = {icc:.3f} < 0.1 (very weak clustering)")
        print(f"  - ANOVA p = {anova_p_value:.3f} > 0.05 (no spatial variation)")
        print("\nUse:")
        print("  sem = std(all_angles) / sqrt(n_measurements)")
        print("\nYou can confidently report all measurements as independent!")
        
    elif icc < 0.2:
        print("≈ HIERARCHICAL STATISTICS (with minimal correction)")
        print("\nReason:")
        print(f"  - ICC = {icc:.3f} is low but non-negligible")
        print("\nUse:")
        print("  - Report both simple and hierarchical SEM for transparency")
        print("  - Use hierarchical for comparisons (slightly conservative)")
        print("  - Mention ICC in methods: 'weak spatial correlation (ICC={icc:.2f})'")
        
    else:
        print("! HIERARCHICAL STATISTICS (essential)")
        print("\nReason:")
        print(f"  - ICC = {icc:.3f} indicates meaningful clustering")
        print(f"  - ANOVA p = {anova_p_value:.3f} < 0.05 (significant spatial variation)")
        print("\nUse:")
        print("  - MUST use hierarchical statistics")
        print("  - Report: 'Accounting for spatial clustering (ICC={icc:.2f})'")
        print("  - Consider modeling spatial trends explicitly")


def quick_icc_check(angles, row_groups, label="Sample"):
    """
    Quick check to see if hierarchical statistics are needed.
    
    Returns True if hierarchical stats needed, False if simple stats OK.
    """
    angles = np.array(angles)
    row_groups = np.array(row_groups)
    
    unique_groups = np.unique(row_groups)
    group_means = [np.mean(angles[row_groups == g]) for g in unique_groups]
    group_vars = [np.var(angles[row_groups == g], ddof=1) if np.sum(row_groups == g) > 1 else 0 
                  for g in unique_groups]
    
    between_var = np.var(group_means, ddof=1)
    within_var = np.mean(group_vars)
    icc = between_var / (between_var + within_var) if (between_var + within_var) > 0 else 0
    
    print(f"\n{label}: ICC = {icc:.4f}", end="")
    
    if icc < 0.1:
        print(" → Simple statistics OK ✓")
        return False
    elif icc < 0.2:
        print(" → Hierarchical stats recommended ≈")
        return True
    else:
        print(" → Hierarchical stats essential !")
        return True


# Example usage
if __name__ == "__main__":
    import numpy as np
    
    # Test Case 1: Truly independent grooves (your scenario?)
    print("\n" + "="*80)
    print("TEST CASE 1: Independent Grooves (Real Manufacturing Variation)")
    print("="*80)
    
    np.random.seed(42)
    
    # Each groove is independent, drawn from same distribution
    # NO spatial correlation
    angles_indep = np.random.normal(17.5, 0.4, 100)
    groups_indep = np.repeat(range(20), 5)
    positions_indep = np.linspace(0, 2, 100)  # 2 µm scan
    
    analysis_indep = analyze_correlation_source(
        angles_indep, groups_indep, positions_indep, 
        "Independent Grooves"
    )
    
    recommend_statistical_approach(
        analysis_indep['icc'], 
        analysis_indep['anova_p_value']
    )
    
    # Test Case 2: Spatially correlated (smooth gradient)
    print("\n\n" + "="*80)
    print("TEST CASE 2: Spatial Gradient (Correlated)")
    print("="*80)
    
    # Grooves vary smoothly across image
    gradient = np.linspace(-0.3, 0.3, 20)  # Smooth gradient across groups
    angles_corr = []
    groups_corr = []
    
    for i, offset in enumerate(gradient):
        group_angles = np.random.normal(17.5 + offset, 0.15, 5)
        angles_corr.extend(group_angles)
        groups_corr.extend([i] * 5)
    
    angles_corr = np.array(angles_corr)
    groups_corr = np.array(groups_corr)
    positions_corr = np.linspace(0, 2, 100)
    
    analysis_corr = analyze_correlation_source(
        angles_corr, groups_corr, positions_corr,
        "Spatial Gradient"
    )
    
    recommend_statistical_approach(
        analysis_corr['icc'],
        analysis_corr['anova_p_value']
    )
    
    plt.show()
