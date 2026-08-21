"""
Statistical analysis functions
Handles comparisons between samples

UPDATED: Proper uncertainty propagation in all comparisons
"""
import numpy as np
from scipy import stats as scipy_stats


def _effective_n(result):
    """
    Sample size to use for inference.

    Row-group analysis re-measures the same physical grooves in every band, so
    the raw count overstates how much independent information there is. Falls
    back to the raw count for results produced before n_effective existed, and
    for traditional (non row-group) analysis where the two are the same thing.
    """
    n_eff = result.get('n_effective')
    if n_eff is None or not np.isfinite(n_eff) or n_eff <= 0:
        return float(result['n_grooves'])
    return float(n_eff)


def print_comparison_summary(results, labels):
    """
    Print a summary table comparing all samples
    
    Parameters:
        results: list of result dictionaries
        labels: list of sample labels
    """
    print(f"\n{'='*80}")
    print("COMPARISON SUMMARY")
    print(f"{'='*80}")
    
    # Determine if we have combined scans
    has_multiple_scans = any('n_scans' in r and r['n_scans'] > 1 for r in results)
    
    if has_multiple_scans:
        print(f"{'Sample':<20} {'File(s)':<30} {'Scans':<6} {'N':<5} {'Mean ± SEM':<18} {'σ_total':<10} {'Spacing':<12}")
        print(f"{'-'*100}")
        for r in results:
            n_scans = r.get('n_scans', 1)
            if n_scans > 1:
                file_str = f"{n_scans} files"
            else:
                file_str = r['filename'].split('/')[-1][:30]
            
            # Get SEM
            sem = r.get('sem_corrected', r.get('sem', r['std_angle'] / np.sqrt(r['n_grooves'])))
            total_std = r.get('total_std', r['std_angle'])
            
            print(f"{r.get('label', 'N/A'):<20} "
                  f"{file_str:<30} "
                  f"{n_scans:<6} "
                  f"{r['n_grooves']:<5} "
                  f"{r['mean_angle']:.2f} ± {sem:.3f}°    "
                  f"{total_std:<10.2f} "
                  f"{r['period_nm']:<10.2f} nm")
    else:
        print(f"{'Sample':<20} {'File':<30} {'N':<5} {'Mean ± SEM':<18} {'σ_total':<10} {'Spacing':<12}")
        print(f"{'-'*100}")
        for r in results:
            file_str = r['filename'].split('/')[-1][:30]
            
            # Get SEM
            sem = r.get('sem_corrected', r.get('sem', r['std_angle'] / np.sqrt(r['n_grooves'])))
            total_std = r.get('total_std', r['std_angle'])
            
            print(f"{r.get('label', 'N/A'):<20} "
                  f"{file_str:<30} "
                  f"{r['n_grooves']:<5} "
                  f"{r['mean_angle']:.2f} ± {sem:.3f}°    "
                  f"{total_std:<10.2f} "
                  f"{r['period_nm']:<10.2f} nm")
    
    print("\nNote: SEM = Standard Error of Mean (uncertainty in mean estimate)")
    print("      σ_total = Total standard deviation (includes measurement + physical variation)")


def print_pairwise_comparisons(results):
    """
    Print pairwise statistical comparisons between samples
    NOW WITH PROPER UNCERTAINTY PROPAGATION AND SIGNIFICANCE TESTS
    
    Parameters:
        results: list of result dictionaries
    """
    print(f"\n{'='*80}")
    print("STATISTICAL COMPARISONS")
    print(f"{'='*80}")
    print("Using proper uncertainty propagation (measurement + physical variation)")
    print("")
    
    for i, r1 in enumerate(results):
        for j, r2 in enumerate(results):
            if i < j:  # Only compare each pair once
                label1 = r1.get('label', f'Sample {i+1}')
                label2 = r2.get('label', f'Sample {j+1}')
                
                diff = r2['mean_angle'] - r1['mean_angle']
                
                # Use total_std if available (includes measurement uncertainty)
                # Otherwise fall back to std_angle (physical variation only)
                std1 = r1.get('total_std', r1['std_angle'])
                std2 = r2.get('total_std', r2['std_angle'])
                
                # Combined standard error of the difference, on the effective
                # sample sizes rather than the raw measurement counts.
                n1, n2 = _effective_n(r1), _effective_n(r2)
                se_combined = np.sqrt(std1**2 / n1 + std2**2 / n2)

                # Calculate t-statistic for significance test
                t_stat = diff / se_combined if se_combined > 0 else 0

                # Degrees of freedom (Welch-Satterthwaite approximation)
                df = _calculate_welch_df(std1, n1, std2, n2)
                
                # Two-tailed p-value
                p_value = 2 * (1 - scipy_stats.t.cdf(abs(t_stat), df))
                
                # 95% confidence interval on the difference
                t_critical = scipy_stats.t.ppf(0.975, df)
                ci_lower = diff - t_critical * se_combined
                ci_upper = diff + t_critical * se_combined
                
                print(f"{label1} vs {label2}:")
                print(f"  Difference: {diff:+.3f}° ± {se_combined:.3f}° (SE)")
                print(f"  95% CI: [{ci_lower:+.3f}°, {ci_upper:+.3f}°]")
                print(f"  {label1}: {r1['mean_angle']:.2f}° ± {std1:.2f}° "
                      f"(σ_total, N={r1['n_grooves']}, N_eff={n1:.1f})")
                print(f"  {label2}: {r2['mean_angle']:.2f}° ± {std2:.2f}° "
                      f"(σ_total, N={r2['n_grooves']}, N_eff={n2:.1f})")
                print(f"  t-statistic: {t_stat:.2f} (df={df:.1f})")
                print(f"  p-value: {p_value:.4f}", end="")
                
                # Interpret significance
                if p_value < 0.001:
                    print(" ***  (highly significant)")
                elif p_value < 0.01:
                    print(" **   (very significant)")
                elif p_value < 0.05:
                    print(" *    (significant)")
                else:
                    print("      (not significant)")
                
                # Effect size (Cohen's d). Deliberately on the raw counts, not
                # the effective ones: this is a standardised mean difference - a
                # description of how far apart the samples are in units of their
                # own spread - not an inference statistic. Correlation belongs in
                # the standard error above, not in the pooled variance.
                pooled_std = np.sqrt(((r1['n_grooves']-1)*std1**2 + (r2['n_grooves']-1)*std2**2) /
                                    (r1['n_grooves'] + r2['n_grooves'] - 2))
                cohens_d = diff / pooled_std if pooled_std > 0 else 0
                print(f"  Effect size (Cohen's d): {cohens_d:.2f}")
                print()


def print_temperature_analysis(results, labels, temperatures):
    """
    Print temperature-dependent analysis with proper uncertainties
    
    Parameters:
        results: list of result dictionaries
        labels: list of sample labels
        temperatures: list of temperatures (None for master)
    """
    # Separate temperature and master samples
    temp_samples = [(i, t, labels[i], results[i]) for i, t in enumerate(temperatures) if t is not None]
    master_samples = [(i, labels[i], results[i]) for i, t in enumerate(temperatures) if t is None]
    
    if len(temp_samples) >= 2:
        print(f"\n{'='*80}")
        print("TEMPERATURE-DEPENDENT CHANGES")
        print(f"{'='*80}")
        
        # Sort by temperature
        temp_samples.sort(key=lambda x: x[1])
        
        # Consecutive temperature steps
        print("\nConsecutive temperature steps:")
        print(f"{'-'*80}")
        
        for idx in range(len(temp_samples) - 1):
            i, temp_i, label_i, result_i = temp_samples[idx]
            j, temp_j, label_j, result_j = temp_samples[idx + 1]
            
            angle_diff = result_j['mean_angle'] - result_i['mean_angle']
            temp_diff = temp_j - temp_i
            
            if temp_diff != 0:
                rate = angle_diff / temp_diff
                
                # Use total_std for uncertainty
                std_i = result_i.get('total_std', result_i['std_angle'])
                std_j = result_j.get('total_std', result_j['std_angle'])
                
                n_i, n_j = _effective_n(result_i), _effective_n(result_j)
                se_combined = np.sqrt(std_i**2 / n_i + std_j**2 / n_j)
                se_rate = se_combined / abs(temp_diff)

                # Significance test
                t_stat = angle_diff / se_combined if se_combined > 0 else 0
                df = _calculate_welch_df(std_i, n_i, std_j, n_j)
                p_value = 2 * (1 - scipy_stats.t.cdf(abs(t_stat), df))
                
                print(f"\n{label_i} ({temp_i}°C) → {label_j} ({temp_j}°C):")
                print(f"  Temperature change: {temp_diff:+.0f}°C")
                print(f"  Angle change: {angle_diff:+.3f}° ± {se_combined:.3f}°")
                print(f"  Rate: {rate:+.6f} ± {se_rate:.6f} °/°C")
                print(f"  Significance: p = {p_value:.4f}", end="")
                if p_value < 0.05:
                    print(" *")
                else:
                    print()
    
    # Master vs treated comparisons
    if len(master_samples) > 0 and len(temp_samples) > 0:
        print(f"\n{'='*80}")
        print("MASTER vs TREATED COMPARISONS")
        print(f"{'='*80}\n")
        
        for i, label_i, result_i in master_samples:
            for j, temp_j, label_j, result_j in temp_samples:
                diff = result_j['mean_angle'] - result_i['mean_angle']
                
                # Use total_std
                std_i = result_i.get('total_std', result_i['std_angle'])
                std_j = result_j.get('total_std', result_j['std_angle'])
                
                n_i, n_j = _effective_n(result_i), _effective_n(result_j)
                se_combined = np.sqrt(std_i**2 / n_i + std_j**2 / n_j)

                # Significance test
                t_stat = diff / se_combined if se_combined > 0 else 0
                df = _calculate_welch_df(std_i, n_i, std_j, n_j)
                p_value = 2 * (1 - scipy_stats.t.cdf(abs(t_stat), df))
                
                print(f"{label_i} → {label_j} ({temp_j}°C):")
                print(f"  Angle change: {diff:+.3f}° ± {se_combined:.3f}°")
                print(f"  {label_i}: {result_i['mean_angle']:.2f}° ± {std_i:.2f}°")
                print(f"  {label_j}: {result_j['mean_angle']:.2f}° ± {std_j:.2f}°")
                print(f"  Significance: p = {p_value:.4f}", end="")
                if p_value < 0.05:
                    print(" *")
                else:
                    print()
                print()


def _calculate_welch_df(std1, n1, std2, n2):
    """
    Calculate Welch-Satterthwaite degrees of freedom for unequal variances
    
    This is the proper df for t-tests when variances are different
    """
    var1 = std1**2
    var2 = std2**2
    
    numerator = (var1/n1 + var2/n2)**2
    denominator = (var1/n1)**2/(n1-1) + (var2/n2)**2/(n2-1)
    
    if denominator > 0:
        df = numerator / denominator
    else:
        df = min(n1, n2) - 1
    
    return df


def print_uncertainty_breakdown(results, labels):
    """
    Print detailed breakdown of uncertainty components
    
    This helps understand where your uncertainty comes from:
    - Measurement uncertainty (from fits)
    - Physical variation (groove-to-groove)
    - Within-image variation (between row groups, if applicable)
    """
    print(f"\n{'='*80}")
    print("UNCERTAINTY BREAKDOWN")
    print(f"{'='*80}")
    print("Understanding your error sources:\n")
    
    for r, label in zip(results, labels):
        print(f"{label}:")
        print(f"  Mean angle: {r['mean_angle']:.2f}°")
        
        # Measurement uncertainty
        if 'mean_measurement_uncertainty' in r:
            print(f"  Measurement uncertainty: {r['mean_measurement_uncertainty']:.3f}° (from fits)")
        
        # Within-image variation (row-group analysis)
        if 'within_image_std' in r and r['within_image_std'] > 0:
            print(f"  Within-image variation: {r['within_image_std']:.3f}° (between row groups)")
        
        # Physical variation
        print(f"  Physical variation: {r['std_angle']:.3f}° (groove-to-groove)")
        
        # Total uncertainty
        total_std = r.get('total_std', r['std_angle'])
        print(f"  Total uncertainty: {total_std:.3f}° (combined)")
        
        # SEM
        sem = r.get('sem', total_std / np.sqrt(r['n_grooves']))
        print(f"  Standard error of mean: {sem:.3f}°")
        print(f"  95% CI on mean: ±{1.96*sem:.3f}°")
        
        # Sample size info
        print(f"  Sample size: N={r['n_grooves']} measurements", end="")
        if 'n_groups' in r and r['n_groups'] > 1:
            print(f" from {r['n_groups']} row groups")
        else:
            print()
        
        # Variance decomposition
        if 'measurement_variance' in r and 'physical_variance' in r:
            meas_var = r['measurement_variance']
            phys_var = r['physical_variance']
            total_var = meas_var + phys_var
            
            if total_var > 0:
                meas_pct = 100 * meas_var / total_var
                phys_pct = 100 * phys_var / total_var
                print(f"  Variance components: {meas_pct:.1f}% measurement, {phys_pct:.1f}% physical")
        
        print()
