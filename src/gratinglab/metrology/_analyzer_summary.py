"""
Statistics, the printed summary and the result dict of one analysed file.

:mod:`.analyzer` measures the grooves; this turns its lists of blaze angles,
quality reports and periods into the statistics, prints them, and packages
them as the dict ``analyze_single_file`` returns. Each step has a traditional
and a row-group version, because a row-group run's measurements are not
independent of each other and its statistics have to say so.
"""
import numpy as np
import matplotlib.pyplot as plt
import os

from .stats.icc import compute_icc, effective_sample_size


# ============ STATISTICS ============

def _calculate_statistics(blaze_angles, quality, all_local_angles, angle_uncertainties):
    """
    Calculate summary statistics with an uncertainty decomposition.

    Two sources of uncertainty are separated:
    1. Measurement uncertainty: how well each groove can be measured.
    2. Physical variation: how much the grooves differ from one another.
    """
    mean_angle = np.mean(blaze_angles)
    
    # Physical variation (groove-to-groove differences)
    physical_variance = np.var(blaze_angles, ddof=1)
    std_angle = np.sqrt(physical_variance)
    
    # Measurement uncertainty (average uncertainty per groove from fits)
    measurement_variance = np.mean([unc**2 for unc in angle_uncertainties]) if len(angle_uncertainties) > 0 else 0
    mean_measurement_uncertainty = np.mean(angle_uncertainties) if len(angle_uncertainties) > 0 else 0
    
    # Total variance = physical variation + measurement uncertainty
    total_variance = physical_variance + measurement_variance
    total_std = np.sqrt(total_variance)
    
    # Standard error of the mean (uncertainty in the mean angle estimate)
    n_grooves = len(blaze_angles)
    sem = total_std / np.sqrt(n_grooves) if n_grooves > 0 else 0
    
    mean_slope = np.mean([q.get('blaze_slope', 0) for q in quality])
    
    # Collect steep angles if available
    steep_angles = [q.get('steep_angle') for q in quality if q.get('steep_angle') is not None]
    mean_steep = np.mean(steep_angles) if len(steep_angles) > 0 else None
    
    # Local angle statistics
    local_angle_std = np.std(all_local_angles) if len(all_local_angles) > 0 else 0
    local_angle_range = (np.max(all_local_angles) - np.min(all_local_angles)) if len(all_local_angles) > 1 else 0
    
    return {
        'mean_angle': mean_angle,
        'std_angle': std_angle,  # This is now ONLY physical variation
        'total_std': total_std,  # Physical + measurement
        'sem': sem,  # Standard error of the mean
        'measurement_variance': measurement_variance,
        'physical_variance': physical_variance,
        'mean_measurement_uncertainty': mean_measurement_uncertainty,
        'mean_slope': mean_slope,
        'mean_steep': mean_steep,
        'local_angle_std': local_angle_std,
        'local_angle_range': local_angle_range
    }


def _calculate_statistics_row_groups(blaze_angles, quality, all_local_angles,
                                     angle_uncertainties, groove_row_groups):
    """
    Calculate statistics for row-group analysis.

    The variance is decomposed into three components:
    1. Measurement uncertainty (how well each groove is fitted)
    2. Within-image variation (variation between row groups)
    3. Groove-to-groove variation (real physical differences)

    Parameters:
        groove_row_groups: row group index of each measurement, same length as
            blaze_angles. This used to be reconstructed by walking the angle list
            in slices the size of each group's detected-centre count, but angles
            are only recorded for grooves whose fit succeeded, so one failed fit
            shifted every later group and a guard then silently dropped the
            tail. Selection by label is not subject to that misalignment.
    """
    stats = _calculate_statistics(blaze_angles, quality, all_local_angles,
                                  angle_uncertainties)

    angles = np.asarray(blaze_angles)
    labels = np.asarray(groove_row_groups)

    if len(labels) != len(angles):
        raise ValueError(f"groove_row_groups and blaze_angles must be the same "
                         f"length, got {len(labels)} and {len(angles)}")

    unique_groups = np.unique(labels)
    stats['n_groups'] = len(unique_groups)

    group_means = [np.mean(angles[labels == g]) for g in unique_groups]

    # Between-group spread, i.e. spatial variation across the image
    stats['within_image_std'] = (float(np.std(group_means, ddof=1))
                                 if len(group_means) > 1 else 0.0)

    # Correct the standard error for the fact that these measurements are not
    # independent: every row group re-measures the same physical grooves. The
    # plain SEM above divides by sqrt(N) over all of them, which is only valid at
    # ICC = 0. Measured ICC on this project's samples is 0.10-0.43.
    #
    # stats/icc.py is the single source of this arithmetic. An earlier prototype
    # (improved_statistics.py, removed 2026-08-09) computed its own ICC with a
    # size-weighted within-group variance; keeping two implementations was
    # exactly the duplication that let the scan-edge bug survive in one code
    # path and not another.
    icc_stats = compute_icc(angles, labels)
    stats['icc'] = icc_stats['icc']
    stats['design_effect'] = float('nan')
    stats['n_effective'] = float(len(angles))
    stats['sem_corrected'] = stats['sem']

    if np.isfinite(icc_stats['icc']):
        n_eff = effective_sample_size(icc_stats['n_measurements'],
                                      icc_stats['mean_group_size'],
                                      icc_stats['icc'])
        stats['design_effect'] = float(len(angles) / n_eff) if n_eff else float('nan')
        stats['n_effective'] = float(n_eff)
        # Same total_std as the uncorrected SEM, over the effective sample size.
        # At ICC = 0, n_eff == N and this equals stats['sem'] exactly.
        stats['sem_corrected'] = (float(stats['total_std'] / np.sqrt(n_eff))
                                  if n_eff > 0 else float('nan'))

    return stats


def _plot_row_group_variation(group_results, all_blaze_angles, all_quality):
    """Plot variation across row groups"""
    # Calculate mean angle for each group
    group_means = []
    group_stds = []
    group_labels = []
    start_idx = 0
    
    for group in group_results:
        n_measurements = group['n_grooves']
        if start_idx + n_measurements <= len(all_blaze_angles):
            group_angles = all_blaze_angles[start_idx:start_idx + n_measurements]
            if len(group_angles) > 0:
                group_means.append(np.mean(group_angles))
                group_stds.append(np.std(group_angles) if len(group_angles) > 1 else 0)
                group_labels.append(f"Group {group['group_idx']+1}")
            start_idx += n_measurements
    
    if len(group_means) < 2:
        return
    
    plt.figure(figsize=(12, 5))
    
    # Plot 1: Mean angle per group
    plt.subplot(1, 2, 1)
    x_pos = np.arange(len(group_means))
    plt.bar(x_pos, group_means, yerr=group_stds, capsize=5, alpha=0.7, color='steelblue')
    plt.axhline(np.mean(group_means), color='r', linestyle='--', linewidth=2, 
                label=f'Overall mean: {np.mean(group_means):.2f}°')
    plt.xlabel('Row Group')
    plt.ylabel('Mean Blaze Angle (degrees)')
    plt.title('Variation Across Row Groups')
    plt.xticks(x_pos, [f"{i+1}" for i in range(len(group_means))], rotation=45)
    plt.legend()
    plt.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: All measurements colored by group
    plt.subplot(1, 2, 2)
    colors = plt.cm.viridis(np.linspace(0, 1, len(group_results)))
    start_idx = 0
    
    for idx, group in enumerate(group_results):
        n_measurements = group['n_grooves']
        if start_idx + n_measurements <= len(all_blaze_angles):
            group_angles = all_blaze_angles[start_idx:start_idx + n_measurements]
            x_vals = np.arange(start_idx, start_idx + len(group_angles))
            plt.scatter(x_vals, group_angles, color=colors[idx], alpha=0.6,
                       label=f"Group {group['group_idx']+1}" if idx < 5 else "")
            start_idx += n_measurements
    
    plt.axhline(np.mean(all_blaze_angles), color='r', linestyle='--', linewidth=2,
                label=f'Overall: {np.mean(all_blaze_angles):.2f}°')
    plt.xlabel('Measurement Number')
    plt.ylabel('Blaze Angle (degrees)')
    plt.title('All Measurements (colored by row group)')
    if len(group_results) <= 5:
        plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()


# ============ PRINTED SUMMARY AND RESULT DICTS ============

def _print_summary(filename, blaze_angles, stats, period_nm, period_std, quality,
                   settings):
    """Print analysis summary to console"""
    print(f"\n{'='*60}")
    print(f"RESULTS FOR {os.path.basename(filename)}")
    print(f"{'='*60}")
    print(f"Analysis side: {settings.blaze_side}")
    print(f"Grooves analyzed: {len(blaze_angles)}")
    print("\nPer-groove statistics (groove-to-groove variation):")
    print(f"  Mean blaze angle: {stats['mean_angle']:.2f} deg ± {stats['std_angle']:.2f} deg (physical variation)")
    print(f"  Min/Max blaze angle: {np.min(blaze_angles):.2f} deg / {np.max(blaze_angles):.2f} deg")
    print(f"  Mean slope (dy/dx): {stats['mean_slope']:.4f}")
    
    # Report uncertainty decomposition
    print("\nUncertainty analysis:")
    print(f"  Average measurement uncertainty per groove: {stats['mean_measurement_uncertainty']:.3f} deg")
    print(f"  Physical variation (groove-to-groove): {stats['std_angle']:.3f} deg")
    print(f"  Total uncertainty (combined): {stats['total_std']:.3f} deg")
    print(f"  Standard error of mean: {stats['sem']:.3f} deg")
    print(f"  95% confidence interval on mean: ±{1.96 * stats['sem']:.3f} deg")
    
    if stats['local_angle_std'] > 0:
        print("\nWithin-facet statistics (camber/curvature):")
        print(f"  Local angle std: {stats['local_angle_std']:.2f} deg")
        print(f"  Local angle range: {stats['local_angle_range']:.2f} deg")
        print(f"  Mean within-facet variation: {np.mean([q.get('angle_std', 0) for q in quality]):.2f} deg")
    
    print("\nGroove geometry:")
    print(f"  Measured groove spacing: {period_nm:.2f} nm ± {period_std:.2f} nm")
    print(f"  Mean groove depth: {np.mean([q['groove_depth_nm'] for q in quality]):.2f} nm")
    print(f"  Mean blaze facet width: {np.mean([q['blaze_width_nm'] for q in quality]):.2f} nm")


def _print_summary_row_groups(filename, blaze_angles, stats, period_nm,
                              period_std, quality, group_info, settings):
    """Print summary for row-group analysis"""
    print(f"\n{'='*60}")
    print(f"RESULTS FOR {os.path.basename(filename)}")
    print(f"{'='*60}")
    print("Analysis mode: ROW-GROUP ANALYSIS")
    print(f"Analysis side: {settings.blaze_side}")
    print(f"Row groups: {stats.get('n_groups', group_info['n_groups'])}")
    print(f"Total measurements: {len(blaze_angles)}")
    
    print("\nOverall statistics:")
    print(f"  Mean blaze angle: {stats['mean_angle']:.2f} deg")
    print(f"  Min/Max blaze angle: {np.min(blaze_angles):.2f} deg / {np.max(blaze_angles):.2f} deg")
    print(f"  Mean slope (dy/dx): {stats['mean_slope']:.4f}")
    
    # Three-way uncertainty decomposition
    print("\nUncertainty analysis (with row-group decomposition):")
    print(f"  Average measurement uncertainty per groove: {stats['mean_measurement_uncertainty']:.3f} deg")
    
    if stats.get('within_image_std', 0) > 0:
        print(f"  Within-image variation (between row groups): {stats['within_image_std']:.3f} deg")
    
    print(f"  Physical variation (groove-to-groove): {stats['std_angle']:.3f} deg")
    print(f"  Total uncertainty (combined): {stats['total_std']:.3f} deg")
    print(f"  Standard error of mean: {stats['sem']:.3f} deg")
    print(f"  95% confidence interval on mean: ±{1.96 * stats['sem']:.3f} deg")
    
    if stats['local_angle_std'] > 0:
        print("\nWithin-facet statistics (camber/curvature):")
        print(f"  Local angle std: {stats['local_angle_std']:.2f} deg")
        print(f"  Local angle range: {stats['local_angle_range']:.2f} deg")
    
    print("\nGroove geometry:")
    print(f"  Measured groove spacing: {period_nm:.2f} nm ± {period_std:.2f} nm")
    print(f"  Mean groove depth: {np.mean([q['groove_depth_nm'] for q in quality]):.2f} nm")
    print(f"  Mean blaze facet width: {np.mean([q['blaze_width_nm'] for q in quality]):.2f} nm")


def _package_results(filename, blaze_angles, stats, period_nm, period_std,
                     groove_periods, quality, all_local_angles,
                     raw_x, flat_y, groove_centers, flat_report=None):
    """Package all results into a dictionary"""
    return {
        'filename': filename,
        'n_grooves': len(blaze_angles),
        'mean_angle': stats['mean_angle'],
        'std_angle': stats['std_angle'],
        'min_angle': np.min(blaze_angles),
        'max_angle': np.max(blaze_angles),
        'mean_slope': stats['mean_slope'],
        'mean_steep': stats['mean_steep'],
        'period_nm': period_nm,
        'period_std': period_std,
        'groove_periods': groove_periods,
        'all_angles': blaze_angles,
        'quality': quality,
        'local_angle_std': stats['local_angle_std'],
        'local_angle_range': stats['local_angle_range'],
        'all_local_angles': all_local_angles,
        **_aggregate_flattening([flat_report]),
        # Uncertainty metrics
        'total_std': stats['total_std'],
        'sem': stats['sem'],
        'measurement_variance': stats['measurement_variance'],
        'physical_variance': stats['physical_variance'],
        'mean_measurement_uncertainty': stats['mean_measurement_uncertainty'],
        # Store profile data for visualization
        'raw_x': raw_x,
        'flat_y': flat_y,
        'groove_centers': groove_centers[:len(quality)]
    }


def _aggregate_flattening(reports):
    """
    One sample's flattening, across however many row groups it had.

    Three things a reader of the summary needs: which method was asked for, how
    many groups did not get it, and how flat the lands actually were. The
    residual is the worst across the groups rather than the mean, because it is
    a question of trust and the weakest reference is the one that matters.

    Groups whose fit had no spare features are left out of that worst case: such
    a fit passes through every point and reports a residual of zero, which would
    otherwise pull the figure towards a flatness nothing measured.
    """
    reports = [r for r in reports if r is not None]
    if not reports:
        return {}

    informative = [r for r in reports
                   if r.residual_rms_nm is not None and r.dof > 0]
    return {
        'flatten_method': reports[0].method_requested,
        'flatten_fallback_groups': sum(1 for r in reports if r.fell_back),
        'flatten_residual_rms_nm': (max(r.residual_rms_nm for r in informative)
                                    if informative else None),
        'flatten_residual_groups': len(informative),
    }


def _package_results_row_groups(filename, blaze_angles, stats, period_nm, period_std,
                                groove_periods, quality, all_local_angles,
                                raw_x, group_results, groove_row_groups):
    """Package row-group analysis results"""
    # Use first group's data for visualization compatibility
    first_group = group_results[0] if len(group_results) > 0 else None
    
    result = {
        'filename': filename,
        'n_grooves': len(blaze_angles),  # Total measurements
        'n_groups': stats.get('n_groups', len(group_results)),
        'mean_angle': stats['mean_angle'],
        'std_angle': stats['std_angle'],
        'min_angle': np.min(blaze_angles),
        'max_angle': np.max(blaze_angles),
        'mean_slope': stats['mean_slope'],
        'mean_steep': stats['mean_steep'],
        'period_nm': period_nm,
        'period_std': period_std,
        'groove_periods': groove_periods,
        'all_angles': blaze_angles,
        # Row group each measurement came from, aligned 1:1 with all_angles.
        # Measurements sharing a label are not independent of each other.
        'groove_row_groups': groove_row_groups,
        'quality': quality,
        'local_angle_std': stats['local_angle_std'],
        'local_angle_range': stats['local_angle_range'],
        'all_local_angles': all_local_angles,
        **_aggregate_flattening([g.get('flatten_report')
                                 for g in group_results]),
        # Uncertainty metrics
        'total_std': stats['total_std'],
        'sem': stats['sem'],
        'measurement_variance': stats['measurement_variance'],
        'physical_variance': stats['physical_variance'],
        'mean_measurement_uncertainty': stats['mean_measurement_uncertainty'],
        'within_image_std': stats.get('within_image_std', 0),
        # Correlation-corrected uncertainty. Row groups re-measure the same
        # physical grooves, so 'sem' above is optimistic; see stats/icc.py.
        'icc': stats.get('icc'),
        'design_effect': stats.get('design_effect'),
        'n_effective': stats.get('n_effective'),
        'sem_corrected': stats.get('sem_corrected'),
        # Store profile data for visualization (use first group)
        'raw_x': raw_x,
        'flat_y': first_group['flat_y'] if first_group else None,
        'groove_centers': first_group['groove_centers'] if first_group else []
    }
    
    return result
