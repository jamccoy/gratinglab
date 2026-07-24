"""
Diagnostic Visualization Functions for Hierarchical AFM Statistics

Creates plots to help understand:
1. Where your uncertainty comes from
2. Whether measurements are normally distributed
3. Consistency across row-groups
4. Quality of uncertainty estimates
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats


def plot_variance_components(stats_dict, label="Sample"):
    """
    Pie chart showing variance decomposition.
    
    Helps answer: Where does my uncertainty come from?
    """
    if not stats_dict['use_hierarchical']:
        print(f"No hierarchical structure for {label}, skipping variance plot")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left plot: Variance components
    components = []
    labels = []
    colors = []
    
    if stats_dict['measurement_variance'] is not None:
        components.append(stats_dict['measurement_variance'])
        labels.append(f"Measurement\n{np.sqrt(stats_dict['measurement_variance']):.3f}°")
        colors.append('#ff9999')
    
    components.append(stats_dict['physical_variance'])
    labels.append(f"Physical\n(within-group)\n{np.sqrt(stats_dict['physical_variance']):.3f}°")
    colors.append('#66b3ff')
    
    components.append(stats_dict['between_group_variance'])
    labels.append(f"Spatial\n(between-group)\n{np.sqrt(stats_dict['between_group_variance']):.3f}°")
    colors.append('#99ff99')
    
    # Calculate percentages
    total_var = sum(components)
    percentages = [100 * c / total_var for c in components]
    
    wedges, texts, autotexts = ax1.pie(
        components, 
        labels=labels,
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 11}
    )
    
    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_weight('bold')
        autotext.set_fontsize(12)
    
    ax1.set_title(f'{label}\nVariance Decomposition', fontsize=13, fontweight='bold')
    
    # Right plot: Standard deviation contributions
    std_components = [np.sqrt(c) for c in components]
    total_std = stats_dict['std_angle']
    
    ax2.barh(range(len(std_components)), std_components, color=colors, edgecolor='black', linewidth=1.5)
    ax2.axvline(total_std, color='red', linestyle='--', linewidth=2, label=f'Total SD = {total_std:.3f}°')
    ax2.set_yticks(range(len(labels)))
    ax2.set_yticklabels([l.replace('\n', ' ') for l in labels], fontsize=10)
    ax2.set_xlabel('Standard Deviation (degrees)', fontsize=11)
    ax2.set_title('SD Components (not additive!)', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(axis='x', alpha=0.3)
    
    # Add interpretation text
    fig.text(0.5, 0.02, 
             f'ICC = {stats_dict["intraclass_correlation"]:.3f} | '
             f'Design Effect = {stats_dict["design_effect"]:.2f} | '
             f'Effective N = {stats_dict["n_effective"]:.1f} / {stats_dict["n_measurements"]}',
             ha='center', fontsize=11, style='italic',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])


def plot_row_group_consistency(stats_dict, label="Sample"):
    """
    Plot mean and uncertainty for each row-group.
    
    Shows spatial variation across the image.
    """
    if not stats_dict['use_hierarchical']:
        print(f"No hierarchical structure for {label}, skipping consistency plot")
        return
    
    group_means = stats_dict['group_means']
    group_stds = stats_dict['group_stds']
    group_sizes = stats_dict['group_sizes']
    n_groups = len(group_means)
    
    # Calculate SE for each group
    group_sems = group_stds / np.sqrt(group_sizes)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Top plot: Group means with error bars
    x = np.arange(n_groups)
    ax1.errorbar(x, group_means, yerr=group_sems, 
                 fmt='o', markersize=8, capsize=5, capthick=2,
                 color='steelblue', ecolor='gray', elinewidth=2,
                 label='Group mean ± SEM')
    
    # Overall mean line
    overall_mean = stats_dict['mean_angle']
    ax1.axhline(overall_mean, color='red', linestyle='--', linewidth=2, 
               label=f'Overall mean = {overall_mean:.3f}°')
    
    # Confidence band (±1.96 SEM conservative)
    sem_cons = stats_dict['sem_conservative']
    ax1.fill_between([-0.5, n_groups-0.5], 
                     overall_mean - 1.96*sem_cons,
                     overall_mean + 1.96*sem_cons,
                     alpha=0.2, color='red',
                     label=f'95% CI on mean (±{1.96*sem_cons:.3f}°)')
    
    ax1.set_xlabel('Row Group Number', fontsize=12)
    ax1.set_ylabel('Blaze Angle (degrees)', fontsize=12)
    ax1.set_title(f'{label}: Spatial Consistency Across Row Groups', 
                 fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-0.5, n_groups-0.5)
    
    # Bottom plot: Range within each group
    ax2.bar(x, group_stds, color='lightcoral', edgecolor='black', linewidth=1.5)
    ax2.axhline(np.sqrt(stats_dict['physical_variance']), 
               color='blue', linestyle='--', linewidth=2,
               label=f'Average within-group SD = {np.sqrt(stats_dict["physical_variance"]):.3f}°')
    ax2.set_xlabel('Row Group Number', fontsize=12)
    ax2.set_ylabel('Within-Group SD (degrees)', fontsize=12)
    ax2.set_title('Groove-to-Groove Variation Within Each Group', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xlim(-0.5, n_groups-0.5)
    
    plt.tight_layout()


def plot_qq_normality(angles, row_group_labels=None, label="Sample"):
    """
    Q-Q plot to assess normality of measurements.
    
    If data is normal, points should follow the red line.
    """
    angles = np.array(angles)
    
    fig, axes = plt.subplots(1, 2 if row_group_labels is not None else 1, 
                            figsize=(12 if row_group_labels is not None else 7, 5))
    
    if row_group_labels is None:
        axes = [axes]
    
    # Left plot: All measurements
    ax = axes[0]
    scipy_stats.probplot(angles, dist="norm", plot=ax)
    ax.set_title(f'{label}: Q-Q Plot (All Measurements)\nN = {len(angles)}', 
                fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add Shapiro-Wilk test result
    if len(angles) >= 3 and len(angles) <= 5000:
        stat, p_value = scipy_stats.shapiro(angles)
        ax.text(0.05, 0.95, 
               f'Shapiro-Wilk:\nW = {stat:.4f}\np = {p_value:.4f}',
               transform=ax.transAxes, fontsize=10,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Right plot: Group means (if hierarchical)
    if row_group_labels is not None:
        ax = axes[1]
        unique_groups = np.unique(row_group_labels)
        group_means = [np.mean(angles[row_group_labels == g]) for g in unique_groups]
        
        scipy_stats.probplot(group_means, dist="norm", plot=ax)
        ax.set_title(f'{label}: Q-Q Plot (Group Means)\nN = {len(group_means)}', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Shapiro-Wilk for group means
        if len(group_means) >= 3 and len(group_means) <= 5000:
            stat, p_value = scipy_stats.shapiro(group_means)
            ax.text(0.05, 0.95, 
                   f'Shapiro-Wilk:\nW = {stat:.4f}\np = {p_value:.4f}',
                   transform=ax.transAxes, fontsize=10,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle('Normality Assessment: Points should follow red line', 
                fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()


def plot_measurement_uncertainty_validation(angles, measurement_errors, 
                                           row_group_labels, label="Sample"):
    """
    Check if measurement uncertainties (from fits) match actual scatter.
    
    If fit uncertainties are realistic:
    - Points should scatter around zero
    - Scatter should match predicted uncertainty
    """
    angles = np.array(angles)
    measurement_errors = np.array(measurement_errors)
    row_group_labels = np.array(row_group_labels)
    
    unique_groups = np.unique(row_group_labels)
    
    # Calculate residuals (deviation from group mean)
    residuals = np.zeros_like(angles)
    for g in unique_groups:
        mask = row_group_labels == g
        group_mean = np.mean(angles[mask])
        residuals[mask] = angles[mask] - group_mean
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left plot: Residuals vs predicted uncertainty
    ax1.scatter(measurement_errors, residuals, alpha=0.6, s=40, edgecolors='black', linewidth=0.5)
    ax1.axhline(0, color='red', linestyle='--', linewidth=2)
    
    # Add ±1 sigma bounds (predicted)
    x_range = np.array([0, max(measurement_errors)])
    ax1.fill_between(x_range, -x_range, x_range, alpha=0.2, color='gray',
                     label='±1σ expected')
    ax1.fill_between(x_range, -2*x_range, 2*x_range, alpha=0.1, color='gray',
                     label='±2σ expected')
    
    ax1.set_xlabel('Predicted Uncertainty (degrees)', fontsize=12)
    ax1.set_ylabel('Residual from Group Mean (degrees)', fontsize=12)
    ax1.set_title('Measurement Uncertainty Validation', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Right plot: Histogram of normalized residuals
    # Normalized residual = residual / predicted_uncertainty
    # Should be ~N(0,1) if uncertainties are correct
    normalized_residuals = residuals / measurement_errors
    
    ax2.hist(normalized_residuals, bins=20, density=True, 
            alpha=0.7, color='steelblue', edgecolor='black')
    
    # Overlay theoretical N(0,1)
    x_theory = np.linspace(-4, 4, 100)
    y_theory = scipy_stats.norm.pdf(x_theory, 0, 1)
    ax2.plot(x_theory, y_theory, 'r-', linewidth=2, label='N(0,1) expected')
    
    ax2.set_xlabel('Normalized Residual (σ units)', fontsize=12)
    ax2.set_ylabel('Density', fontsize=12)
    ax2.set_title('Distribution of Normalized Residuals', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add statistics
    mean_norm = np.mean(normalized_residuals)
    std_norm = np.std(normalized_residuals)
    ax2.text(0.05, 0.95,
            f'Mean: {mean_norm:.2f}\nSD: {std_norm:.2f}\n(Expect: 0, 1)',
            transform=ax2.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle(f'{label}: Are Fit Uncertainties Realistic?', 
                fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()


def plot_comparison_forest(comparisons, reference_label=None):
    """
    Forest plot showing all pairwise comparisons.
    
    Visualizes effect sizes and confidence intervals.
    """
    n_comp = len(comparisons)
    
    fig, ax = plt.subplots(figsize=(10, max(6, n_comp * 0.4)))
    
    y_positions = np.arange(n_comp)
    
    for i, comp in enumerate(comparisons):
        diff = comp['difference']
        ci_lower = comp['ci_95_lower']
        ci_upper = comp['ci_95_upper']
        
        # Color by significance
        if comp['significant_001']:
            color = 'darkgreen'
            marker = 'D'
        elif comp['significant_01']:
            color = 'green'
            marker = 'o'
        elif comp['significant_05']:
            color = 'orange'
            marker = 'o'
        else:
            color = 'gray'
            marker = 'o'
        
        # Plot point and CI
        ax.plot([ci_lower, ci_upper], [i, i], color=color, linewidth=2)
        ax.plot(diff, i, marker=marker, color=color, markersize=10,
               markeredgecolor='black', markeredgewidth=1)
        
        # Add label
        label_text = f"{comp['label2']} - {comp['label1']}"
        ax.text(-0.02, i, label_text, ha='right', va='center', fontsize=10,
               transform=ax.get_yaxis_transform())
    
    # Zero line
    ax.axvline(0, color='black', linestyle='--', linewidth=1.5)
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels([''] * n_comp)
    ax.set_xlabel('Difference in Blaze Angle (degrees)', fontsize=12)
    ax.set_title('Pairwise Comparisons: Forest Plot\n(Error bars = 95% CI)', 
                fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='D', color='w', markerfacecolor='darkgreen',
              markersize=10, label='p < 0.001 ***', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green',
              markersize=10, label='p < 0.01 **', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange',
              markersize=10, label='p < 0.05 *', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
              markersize=10, label='n.s.', markeredgecolor='black')
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=10)
    
    plt.tight_layout()


def create_diagnostic_report(stats_dict, angles, row_group_labels, 
                            measurement_errors, label="Sample"):
    """
    Create a comprehensive diagnostic report with all plots.
    """
    print(f"\n{'='*80}")
    print(f"CREATING DIAGNOSTIC PLOTS FOR: {label}")
    print(f"{'='*80}\n")
    
    # 1. Variance decomposition
    print("1. Variance components...")
    plot_variance_components(stats_dict, label)
    
    # 2. Row-group consistency
    print("2. Row-group consistency...")
    plot_row_group_consistency(stats_dict, label)
    
    # 3. Normality assessment
    print("3. Normality tests...")
    plot_qq_normality(angles, row_group_labels, label)
    
    # 4. Measurement uncertainty validation
    print("4. Measurement uncertainty validation...")
    plot_measurement_uncertainty_validation(angles, measurement_errors,
                                           row_group_labels, label)
    
    print("\n✓ All diagnostic plots created!")
    print(f"{'='*80}\n")


# Example usage
if __name__ == "__main__":
    # Generate example data
    np.random.seed(42)
    n_groups = 20
    n_per_group = 5
    
    group_means = np.random.normal(17.5, 0.3, n_groups)
    
    angles = []
    groups = []
    errors = []
    
    for i, gm in enumerate(group_means):
        group_angles = np.random.normal(gm, 0.15, n_per_group)
        angles.extend(group_angles)
        groups.extend([i] * n_per_group)
        errors.extend(np.random.uniform(0.05, 0.15, n_per_group))
    
    angles = np.array(angles)
    groups = np.array(groups)
    errors = np.array(errors)
    
    # Import the statistics module (sits alongside this file)
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from improved_statistics import calculate_hierarchical_statistics
    
    stats = calculate_hierarchical_statistics(angles, groups, errors)
    
    # Create all diagnostic plots
    create_diagnostic_report(stats, angles, groups, errors, "Example Sample")
    
    plt.show()
