"""
Core AFM file analyzer
Handles single file analysis with all the detailed logic
Now supports row-group analysis for many more measurements per image!
"""
import numpy as np
import matplotlib.pyplot as plt
import os

from .settings import AnalysisSettings
from .stats.icc import compute_icc, effective_sample_size
from .core.processing import (
    raw_data, raw_data_multi_group,  # NEW: added raw_data_multi_group
    flatten_profile, find_groove_positions, load_afm_data
)
from .core.analysis import extract_blaze_angle
from .visualization.diagnostics import plot_analyzed_regions_overlay, plot_flattening_diagnostic
from .visualization.statistics import plot_summary_statistics


def analyze_single_file(filename, show_plots=True, settings=None):
    """
    Analyze a single AFM file and extract blaze angles
    
    Supports two modes:
    1. Traditional: Average all rows → single profile → few measurements
    2. Row-group: Multiple row-group profiles → many measurements

    Set USE_ROW_GROUPS = True in config.py to enable row-group analysis, or pass
    settings with use_row_groups set.

    Parameters:
        filename: path to AFM data file
        show_plots: whether to show diagnostic plots
        settings: AnalysisSettings to use. Defaults to the values in config.py,
            which is what every config-driven workflow relies on.

    Returns:
        dict with analysis results, or None if analysis failed
    """
    if settings is None:
        settings = AnalysisSettings.from_config()

    # Check if row-group analysis is enabled
    if settings.use_row_groups:
        return analyze_single_file_with_row_groups(filename, show_plots, settings)
    else:
        return _analyze_single_file_traditional(filename, show_plots, settings)


def analyze_single_file_with_row_groups(filename, show_plots=True, settings=None):
    """
    Analyze a single AFM file using row-group analysis
    
    This extracts multiple profiles from different regions of the image,
    giving you many more measurements and allowing assessment of within-image variation.
    
    Parameters:
        filename: path to AFM data file
        show_plots: whether to show diagnostic plots
        
    Returns:
        dict with analysis results, or None if analysis failed
    """
    
    if settings is None:
        settings = AnalysisSettings.from_config()

    print(f"\n{'='*60}")
    print(f"Analyzing: {filename}")
    print(f"Mode: ROW-GROUP ANALYSIS (n_groups={settings.n_row_groups})")
    print(f"{'='*60}")
    
    # Load and validate data
    data, scan_x_size = _load_and_validate(filename, settings)
    if data is None:
        return None
    
    # Show 2D image if requested
    if settings.show_2d_image and show_plots:
        _plot_2d_image(data, scan_x_size, filename)
    
    # Extract multiple profiles from row groups
    raw_x, profiles_list, group_info = raw_data_multi_group(data, scan_x_size, settings.n_row_groups)
    
    # Calculate scan parameters (same for all profiles)
    scan_width_nm = scan_x_size * 1000
    estimated_grooves = max(2, int(scan_width_nm / settings.period_est))
    period_nm_est = scan_width_nm / estimated_grooves
    
    print(f"Scan width: {scan_width_nm:.1f} nm")
    print(f"Estimated period: {period_nm_est:.2f} nm")
    
    # Process each row-group profile
    all_blaze_angles = []
    all_quality = []
    all_local_angles = []
    all_angle_uncertainties = []
    all_groove_periods = []
    group_results = []  # Store results per group for diagnostics
    n_edge_rejected = 0  # Grooves dropped for sitting on a scan edge
    all_groove_positions = []  # x position (um) of each accepted measurement
    all_groove_groups = []     # row group index of each accepted measurement

    for group_idx, raw_y in enumerate(profiles_list):
        print(f"\n  Processing row group {group_idx + 1}/{settings.n_row_groups}...")
        
        # Flatten this profile
        flat_y, background = flatten_profile(raw_x, raw_y,
                                             method=settings.flatten_method,
                                             poly_order=settings.flatten_poly_order,
                                             exclude_edges=settings.flatten_exclude_edges,
                                             feature=settings.flatten_feature,
                                             period_nm=period_nm_est)
        
        # Show flattening diagnostic for first group only
        if show_plots and settings.show_flattening_diagnostic and group_idx == 0:
            plot_flattening_diagnostic(raw_x, raw_y, flat_y, background,
                                      settings.flatten_method, settings.flatten_feature, period_nm_est)
        
        # Find grooves in this profile
        groove_centers, n_edge = find_groove_positions(raw_x, flat_y, period_nm_est,
                                              prominence_factor=settings.prominence_factor,
                                              distance_factor=settings.distance_factor,
                                              edge_exclusion=settings.edge_exclusion_periods,
                                              return_n_edge_rejected=True)
        n_edge_rejected += n_edge

        edge_note = f" ({n_edge} rejected: clipped by scan edge)" if n_edge else ""
        print(f"    Found {len(groove_centers)} grooves{edge_note}")

        if len(groove_centers) == 0:
            print(f"    No grooves detected in group {group_idx + 1}, skipping")
            continue
        
        # Calculate local periods for this group
        if len(groove_centers) > 1:
            groove_spacings = np.diff(groove_centers) * (raw_x[1] - raw_x[0]) * 1000
            local_period_nm = np.mean(groove_spacings)
            groove_periods_this_group = groove_spacings.tolist()
        else:
            local_period_nm = period_nm_est
            groove_periods_this_group = []
        
        # Calculate local periods for each groove
        local_periods = _calculate_local_periods_for_grooves(
            groove_centers, groove_periods_this_group, local_period_nm
        )
        
        # Extract blaze angles from all grooves in this group
        show_individual = show_plots and settings.show_individual_grooves and group_idx == 0
        
        for groove_idx, (center, local_period) in enumerate(zip(groove_centers, local_periods)):
            angle, steep, slope, qual = extract_blaze_angle(
                raw_x, flat_y, center, local_period_nm,
                trim_fraction=settings.facet_trim,
                side=settings.blaze_side,
                show_plot=show_individual,
                groove_num=groove_idx + 1,
                return_local_angles=True,
                return_regions=show_plots and settings.show_analyzed_regions and group_idx == 0,
                local_period_nm=local_period
            )
            
            if angle is not None:
                all_blaze_angles.append(angle)
                all_quality.append(qual)
                all_angle_uncertainties.append(qual.get('blaze_angle_stderr', 0))
                # Position of this specific measurement, so plots can pair angles
                # with positions 1:1 even though each groove is measured once per
                # row group. Recorded here rather than reconstructed later,
                # because grooves whose fit failed never reach this branch.
                all_groove_positions.append(raw_x[center])
                # Which row group this measurement came from. Required to compute
                # the intraclass correlation - measurements from one group are not
                # independent of each other - and recorded here for the same
                # reason as the position above: a groove whose fit failed never
                # reaches this branch, so counting detected centres afterwards
                # would misattribute every later measurement.
                all_groove_groups.append(group_idx)

                if 'local_angles' in qual:
                    all_local_angles.extend(qual['local_angles'])
        
        # Store groove periods from this group
        all_groove_periods.extend(groove_periods_this_group)
        
        # Store group result
        group_results.append({
            'group_idx': group_idx,
            'n_grooves': len(groove_centers),
            'raw_y': raw_y,
            'flat_y': flat_y,
            'groove_centers': groove_centers
        })
    
    # Check if we got any measurements
    if len(all_blaze_angles) == 0:
        print("Could not extract any blaze angles from any row group!")
        return None
    
    print(f"\n{'='*60}")
    print(f"ROW-GROUP ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total measurements: {len(all_blaze_angles)}")
    print(f"  (compared to ~{len(group_results) * 4} with traditional averaging)")
    print(f"Row groups processed: {len(group_results)}/{settings.n_row_groups}")
    if n_edge_rejected:
        print(f"Edge-clipped grooves rejected: {n_edge_rejected} "
              f"(within {settings.edge_exclusion_periods}x period of a scan edge)")
    
    # Calculate overall period statistics
    if len(all_groove_periods) > 0:
        period_nm = np.mean(all_groove_periods)
        period_std = np.std(all_groove_periods)
        print(f"Measured period: {period_nm:.2f} ± {period_std:.2f} nm")
    else:
        period_nm = period_nm_est
        period_std = 0
    
    # Show analyzed regions for first group only
    if show_plots and settings.show_analyzed_regions and len(group_results) > 0:
        first_group = group_results[0]
        n_qual_first = min(len(first_group['groove_centers']), 
                          sum(1 for q in all_quality if 'regions' in q))
        
        if n_qual_first > 0:
            qual_with_regions = [q for q in all_quality if 'regions' in q][:n_qual_first]
            plot_analyzed_regions_overlay(
                raw_x, first_group['flat_y'], 
                first_group['groove_centers'][:n_qual_first],
                period_nm, qual_with_regions, settings.facet_trim
            )
    
    # Calculate statistics with row-group data
    stats = _calculate_statistics_row_groups(
        all_blaze_angles, all_quality, all_local_angles,
        all_angle_uncertainties, all_groove_groups
    )
    
    # Summary plots
    if show_plots:
        first_group = group_results[0] if len(group_results) > 0 else None

        if first_group is not None:
            plot_summary_statistics(
                all_blaze_angles, all_groove_positions,
                all_local_angles, stats['mean_angle'], stats['std_angle'],
                settings.show_local_angle_distribution
            )
        
        # Additional row-group specific plot
        if len(group_results) > 1:
            _plot_row_group_variation(group_results, all_blaze_angles, all_quality)
    
    # Print summary
    _print_summary_row_groups(filename, all_blaze_angles, stats,
                             period_nm, period_std, all_quality, group_info,
                             settings)
    
    # Package results
    result = _package_results_row_groups(
        filename, all_blaze_angles, stats, period_nm, period_std,
        all_groove_periods, all_quality, all_local_angles,
        raw_x, group_results, all_groove_groups
    )
    
    return result


def _analyze_single_file_traditional(filename, show_plots=True, settings=None):
    """
    Traditional analysis: average all rows into single profile
    
    This is the original analysis method. Kept for compatibility and comparison.
    """
    
    if settings is None:
        settings = AnalysisSettings.from_config()

    print(f"\n{'='*60}")
    print(f"Analyzing: {filename}")
    print(f"Mode: TRADITIONAL (single averaged profile)")
    print(f"{'='*60}")
    
    # Load and validate data
    data, scan_x_size = _load_and_validate(filename, settings)
    if data is None:
        return None
    
    # Show 2D image if requested
    if settings.show_2d_image and show_plots:
        _plot_2d_image(data, scan_x_size, filename)
    
    # Extract and flatten profile
    raw_x, raw_y, flat_y, background, period_nm_est = _process_profile(data, scan_x_size, show_plots, settings)
    
    # Find grooves
    groove_centers, n_edge_rejected = find_groove_positions(raw_x, flat_y, period_nm_est,
                                          prominence_factor=settings.prominence_factor,
                                          distance_factor=settings.distance_factor,
                                          edge_exclusion=settings.edge_exclusion_periods,
                                          return_n_edge_rejected=True)

    edge_note = f" ({n_edge_rejected} rejected: clipped by scan edge)" if n_edge_rejected else ""
    print(f"Found {len(groove_centers)} grooves{edge_note}")

    if len(groove_centers) == 0:
        print("No grooves detected!")
        return None
    
    # Calculate groove spacing
    period_nm, period_std, groove_periods = _calculate_periods(raw_x, groove_centers, period_nm_est)
    local_periods = _calculate_local_periods(groove_centers, groove_periods, period_nm)
    
    # Plot full profile if requested
    if show_plots and settings.show_full_profile:
        _plot_full_profile(raw_x, flat_y, groove_centers, filename)
    
    # Extract blaze angles
    blaze_angles, quality, all_local_angles, angle_uncertainties = _extract_all_angles(
        raw_x, flat_y, groove_centers, period_nm, local_periods, show_plots, settings
    )
    
    if len(blaze_angles) == 0:
        print("Could not extract blaze angles!")
        return None
    
    # Show analyzed regions if requested
    if show_plots and settings.show_analyzed_regions:
        plot_analyzed_regions_overlay(raw_x, flat_y, groove_centers[:len(quality)],
                                     period_nm, quality, settings.facet_trim)
    
    # Calculate and display statistics
    stats = _calculate_statistics(blaze_angles, quality, all_local_angles, angle_uncertainties)
    
    # Summary plots
    if show_plots:
        plot_summary_statistics(blaze_angles,
                              raw_x[groove_centers[:len(blaze_angles)]],
                              all_local_angles,
                              stats['mean_angle'], stats['std_angle'], 
                              settings.show_local_angle_distribution)
    
    # Print summary
    _print_summary(filename, blaze_angles, stats, period_nm, period_std, quality,
                   settings)
    
    # Combine all results
    return _package_results(filename, blaze_angles, stats, period_nm, period_std,
                           groove_periods, quality, all_local_angles,
                           raw_x, flat_y, groove_centers)


# ============ HELPER FUNCTIONS ============

def _load_and_validate(filename, settings):
    """Load AFM data and validate it"""
    try:
        data, scan_x_size = load_afm_data(filename,
                                          default_scan_size=settings.scan_x_size,
                                          settings=settings)
        return data, scan_x_size
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None, None


def _plot_2d_image(data, scan_x_size, filename):
    """Plot 2D AFM topography image"""
    plt.figure(figsize=(10, 6))
    plt.imshow(data, aspect='auto', cmap='viridis',
              extent=[0, scan_x_size, 0, len(data)])
    plt.title(f'2D AFM Topography: {os.path.basename(filename).replace("_", " ")}')
    plt.xlabel('X (µm)')
    plt.ylabel('Y (scan lines)')
    plt.colorbar(label='Height (m)')


def _process_profile(data, scan_x_size, show_plots, settings):
    """Extract profile, flatten it, and show diagnostics"""
    # Extract profile
    raw_x, raw_y = raw_data(data, scan_x_size)
    
    # Calculate scan parameters
    scan_width_nm = scan_x_size * 1000
    estimated_grooves = max(2, int(scan_width_nm / settings.period_est))
    period_nm_est = scan_width_nm / estimated_grooves
    
    print(f"Scan width: {scan_width_nm:.1f} nm")
    print(f"Estimated period: {period_nm_est:.2f} nm")
    
    # Flatten profile
    flat_y, background = flatten_profile(raw_x, raw_y,
                                         method=settings.flatten_method,
                                         poly_order=settings.flatten_poly_order,
                                         exclude_edges=settings.flatten_exclude_edges,
                                         feature=settings.flatten_feature,
                                         period_nm=period_nm_est)
    
    # Show flattening diagnostic if requested
    if show_plots and settings.show_flattening_diagnostic:
        plot_flattening_diagnostic(raw_x, raw_y, flat_y, background,
                                  settings.flatten_method, settings.flatten_feature, period_nm_est)
    
    # Print flattening info
    _print_flattening_info(raw_x, background, settings)
    
    return raw_x, raw_y, flat_y, background, period_nm_est


def _print_flattening_info(raw_x, background, settings):
    """Print information about the flattening process"""
    tilt_angle = np.arctan(np.polyfit(raw_x, background, 1)[0] / 1000) * 180 / np.pi
    print(f"Flattening method: {settings.flatten_method}")
    if settings.flatten_method == 'polynomial':
        print(f"  Polynomial order: {settings.flatten_poly_order}")
    elif settings.flatten_method == 'level_grooves':
        print(f"  Polynomial order: {settings.flatten_poly_order}")
        print(f"  Feature: {settings.flatten_feature}")
    print(f"  Approximate tilt angle: {tilt_angle:.4f} degrees")
    print(f"  Background range: {np.max(background) - np.min(background):.2f} nm")


def _calculate_periods(raw_x, groove_centers, period_nm_est):
    """Calculate groove periods from detected positions"""
    groove_periods = []
    
    if len(groove_centers) > 1:
        groove_spacings = np.diff(groove_centers) * (raw_x[1] - raw_x[0]) * 1000
        groove_periods = groove_spacings.tolist()
        actual_period = np.mean(groove_spacings)
        period_std = np.std(groove_spacings)
        print(f"Measured period: {actual_period:.2f} ± {period_std:.2f} nm")
        print(f"  Individual spacings: min={np.min(groove_spacings):.2f} nm, "
              f"max={np.max(groove_spacings):.2f} nm")
        period_nm = actual_period
    else:
        period_nm = period_nm_est
        period_std = 0
    
    return period_nm, period_std, groove_periods


def _calculate_local_periods(groove_centers, groove_periods, period_nm):
    """Calculate local period for each groove (average of neighbors)"""
    local_periods = []
    
    for i in range(len(groove_centers)):
        if len(groove_centers) == 1:
            # Only one groove - use global estimate
            local_periods.append(period_nm)
        elif i == 0:
            # First groove - use spacing to next groove
            local_periods.append(groove_periods[0])
        elif i == len(groove_centers) - 1:
            # Last groove - use spacing from previous groove
            local_periods.append(groove_periods[-1])
        else:
            # Middle grooves - average of left and right spacing
            local_periods.append((groove_periods[i-1] + groove_periods[i]) / 2)
    
    return local_periods


def _calculate_local_periods_for_grooves(groove_centers, groove_periods, period_nm):
    """Helper for row-group analysis: calculate local periods"""
    return _calculate_local_periods(groove_centers, groove_periods, period_nm)


def _plot_full_profile(raw_x, flat_y, groove_centers, filename):
    """Plot the full profile with detected grooves"""
    plt.figure(figsize=(12, 4))
    plt.plot(raw_x, flat_y, 'k-', linewidth=0.5, label='Full profile')
    plt.plot(raw_x[groove_centers], flat_y[groove_centers], 'ro',
            markersize=8, label=f'Detected grooves (N={len(groove_centers)})')
    plt.xlabel('Displacement (µm)')
    plt.ylabel('Height (nm)')
    plt.title(f'Full Profile: {os.path.basename(filename).replace("_", " ")}')
    plt.legend()
    plt.grid(True, alpha=0.3)


def _extract_all_angles(raw_x, flat_y, groove_centers, period_nm, local_periods,
                        show_plots, settings):
    """Extract blaze angles from all grooves"""
    blaze_angles = []
    steep_angles = []
    slopes = []
    quality = []
    all_local_angles = []
    angle_uncertainties = []
    
    for i, (center, local_period) in enumerate(zip(groove_centers, local_periods)):
        show_individual = show_plots and settings.show_individual_grooves
        
        angle, steep, slope, qual = extract_blaze_angle(
            raw_x, flat_y, center, period_nm,
            trim_fraction=settings.facet_trim,
            side=settings.blaze_side,
            show_plot=show_individual,
            groove_num=i+1,
            return_local_angles=True,
            return_regions=True,
            local_period_nm=local_period
        )
        
        if angle is not None:
            blaze_angles.append(angle)
            slopes.append(slope)
            quality.append(qual)
            angle_uncertainties.append(qual.get('blaze_angle_stderr', 0))
            if steep is not None:
                steep_angles.append(steep)
            
            if 'local_angles' in qual:
                all_local_angles.extend(qual['local_angles'])
    
    return blaze_angles, quality, all_local_angles, angle_uncertainties


def _calculate_statistics(blaze_angles, quality, all_local_angles, angle_uncertainties):
    """
    Calculate summary statistics with proper uncertainty decomposition
    
    Separates two sources of uncertainty:
    1. Measurement uncertainty: How well can we measure each groove?
    2. Physical variation: How much do grooves actually differ?
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
    Calculate statistics for row-group analysis

    Decomposes variance into three components:
    1. Measurement uncertainty (how well we fit each groove)
    2. Within-image variation (variation between row groups)
    3. Groove-to-groove variation (real physical differences)

    Parameters:
        groove_row_groups: row group index of each measurement, same length as
            blaze_angles. This used to be reconstructed by walking the angle list
            in slices the size of each group's detected-centre count, but angles
            are only recorded for grooves whose fit succeeded - so one failed fit
            shifted every later group, and a guard then silently dropped the
            tail. Selecting by label cannot drift.
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
    # stats/icc.py is the single source of this arithmetic. improved_statistics.py
    # in experimental/ computes its own ICC with a size-weighted within-group
    # variance; keeping two implementations is exactly the duplication that let
    # the scan-edge bug survive in one code path and not another.
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


def _print_summary(filename, blaze_angles, stats, period_nm, period_std, quality,
                   settings):
    """Print analysis summary to console"""
    print(f"\n{'='*60}")
    print(f"RESULTS FOR {os.path.basename(filename)}")
    print(f"{'='*60}")
    print(f"Analysis side: {settings.blaze_side}")
    print(f"Grooves analyzed: {len(blaze_angles)}")
    print(f"\nPer-groove statistics (groove-to-groove variation):")
    print(f"  Mean blaze angle: {stats['mean_angle']:.2f} deg ± {stats['std_angle']:.2f} deg (physical variation)")
    print(f"  Min/Max blaze angle: {np.min(blaze_angles):.2f} deg / {np.max(blaze_angles):.2f} deg")
    print(f"  Mean slope (dy/dx): {stats['mean_slope']:.4f}")
    
    # Report uncertainty decomposition
    print(f"\nUncertainty analysis:")
    print(f"  Average measurement uncertainty per groove: {stats['mean_measurement_uncertainty']:.3f} deg")
    print(f"  Physical variation (groove-to-groove): {stats['std_angle']:.3f} deg")
    print(f"  Total uncertainty (combined): {stats['total_std']:.3f} deg")
    print(f"  Standard error of mean: {stats['sem']:.3f} deg")
    print(f"  95% confidence interval on mean: ±{1.96 * stats['sem']:.3f} deg")
    
    if stats['local_angle_std'] > 0:
        print(f"\nWithin-facet statistics (camber/curvature):")
        print(f"  Local angle std: {stats['local_angle_std']:.2f} deg")
        print(f"  Local angle range: {stats['local_angle_range']:.2f} deg")
        print(f"  Mean within-facet variation: {np.mean([q.get('angle_std', 0) for q in quality]):.2f} deg")
    
    print(f"\nGroove geometry:")
    print(f"  Measured groove spacing: {period_nm:.2f} nm ± {period_std:.2f} nm")
    print(f"  Mean groove depth: {np.mean([q['groove_depth_nm'] for q in quality]):.2f} nm")
    print(f"  Mean blaze facet width: {np.mean([q['blaze_width_nm'] for q in quality]):.2f} nm")


def _print_summary_row_groups(filename, blaze_angles, stats, period_nm,
                              period_std, quality, group_info, settings):
    """Print summary for row-group analysis"""
    print(f"\n{'='*60}")
    print(f"RESULTS FOR {os.path.basename(filename)}")
    print(f"{'='*60}")
    print(f"Analysis mode: ROW-GROUP ANALYSIS")
    print(f"Analysis side: {settings.blaze_side}")
    print(f"Row groups: {stats.get('n_groups', group_info['n_groups'])}")
    print(f"Total measurements: {len(blaze_angles)}")
    
    print(f"\nOverall statistics:")
    print(f"  Mean blaze angle: {stats['mean_angle']:.2f} deg")
    print(f"  Min/Max blaze angle: {np.min(blaze_angles):.2f} deg / {np.max(blaze_angles):.2f} deg")
    print(f"  Mean slope (dy/dx): {stats['mean_slope']:.4f}")
    
    # Three-way uncertainty decomposition
    print(f"\nUncertainty analysis (with row-group decomposition):")
    print(f"  Average measurement uncertainty per groove: {stats['mean_measurement_uncertainty']:.3f} deg")
    
    if stats.get('within_image_std', 0) > 0:
        print(f"  Within-image variation (between row groups): {stats['within_image_std']:.3f} deg")
    
    print(f"  Physical variation (groove-to-groove): {stats['std_angle']:.3f} deg")
    print(f"  Total uncertainty (combined): {stats['total_std']:.3f} deg")
    print(f"  Standard error of mean: {stats['sem']:.3f} deg")
    print(f"  95% confidence interval on mean: ±{1.96 * stats['sem']:.3f} deg")
    
    if stats['local_angle_std'] > 0:
        print(f"\nWithin-facet statistics (camber/curvature):")
        print(f"  Local angle std: {stats['local_angle_std']:.2f} deg")
        print(f"  Local angle range: {stats['local_angle_range']:.2f} deg")
    
    print(f"\nGroove geometry:")
    print(f"  Measured groove spacing: {period_nm:.2f} nm ± {period_std:.2f} nm")
    print(f"  Mean groove depth: {np.mean([q['groove_depth_nm'] for q in quality]):.2f} nm")
    print(f"  Mean blaze facet width: {np.mean([q['blaze_width_nm'] for q in quality]):.2f} nm")


def _package_results(filename, blaze_angles, stats, period_nm, period_std,
                     groove_periods, quality, all_local_angles,
                     raw_x, flat_y, groove_centers):
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
