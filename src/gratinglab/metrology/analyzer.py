"""
Core AFM file analyzer.

This module handles the analysis of a single file in full detail. It also
supports row-group analysis, which yields many more measurements per image.
"""
import numpy as np
import matplotlib.pyplot as plt
import os

from .settings import AnalysisSettings
from .core.image_flatten import flatten_image, row_offset_spread
from .core.orientation import apply_groove_orientation
from .core.region_mask import (GROOVE_REACH_PERIODS, apply_region_mask,
                               bridge_nan, describe_mask, groove_spacings_px,
                               reject_grooves_near_nan,
                               tip_correct_around_nan)
from .core.scars import repair_scars_for
from .core.processing import (
    raw_data, raw_data_multi_group,  # NEW: added raw_data_multi_group
    flatten_profile, find_groove_positions, load_afm_data, load_scan_y_size
)
from .core.analysis import extract_blaze_angle
from .visualization.diagnostics import plot_analyzed_regions_overlay, plot_flattening_diagnostic
from .visualization.statistics import plot_summary_statistics
from ._analyzer_summary import (
    _calculate_statistics,
    _calculate_statistics_row_groups,
    _plot_row_group_variation,
    _print_summary,
    _print_summary_row_groups,
    _package_results,
    _package_results_row_groups,
)


# ============ ENTRY POINTS: ONE FILE, IN EITHER MODE ============

def analyze_single_file(filename, show_plots=True, settings=None):
    """
    Analyze a single AFM file and extract blaze angles
    
    Two modes are supported:
    1. Traditional: all rows are averaged into a single profile, which gives
       few measurements.
    2. Row-group: several row-group profiles are extracted, which gives many
       measurements.

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
    
    This extracts multiple profiles from different regions of the image, which
    gives many more measurements and allows within-image variation to be
    assessed.
    
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
    data, scan_x_size, prep = _load_and_validate(filename, settings)
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
    n_mask_rejected = 0  # Grooves dropped for reaching into a masked column
    all_groove_positions = []  # x position (um) of each accepted measurement
    all_groove_groups = []     # row group index of each accepted measurement

    # ------------ each row group: flatten, find grooves, fit facets ------------
    for group_idx, raw_y in enumerate(profiles_list):
        print(f"\n  Processing row group {group_idx + 1}/{settings.n_row_groups}...")

        # Columns masked across this whole group are NaN. Bridged for
        # flattening and detection only; grooves reading them are dropped below.
        measured_y = raw_y
        raw_y = bridge_nan(raw_y)

        # Flatten this profile
        flat_y, background, flat_report = flatten_profile(raw_x, raw_y,
                                             method=settings.flatten_method,
                                             poly_order=settings.flatten_poly_order,
                                             exclude_edges=settings.flatten_exclude_edges,
                                             feature=settings.flatten_feature,
                                             period_nm=period_nm_est,
                                             return_report=True,
                                             exclude=~np.isfinite(measured_y))
        if flat_report.fell_back:
            print(f"    Flattening: {flat_report.summary()}")
        
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
        groove_centers, n_masked = reject_grooves_near_nan(
            groove_centers, measured_y,
            GROOVE_REACH_PERIODS * period_nm_est / ((raw_x[1] - raw_x[0]) * 1000))
        n_mask_rejected += n_masked

        edge_note = f" ({n_edge} rejected: clipped by scan edge)" if n_edge else ""
        if n_masked:
            edge_note += f" ({n_masked} rejected: masked region)"
        print(f"    Found {len(groove_centers)} grooves{edge_note}")

        if len(groove_centers) == 0:
            print(f"    No grooves detected in group {group_idx + 1}, skipping")
            continue
        
        # Calculate local periods for this group
        # A spacing across a masked gap is NaN: the grooves either side of it
        # are not neighbours (core/region_mask.groove_spacings_px).
        groove_spacings = (groove_spacings_px(groove_centers, measured_y)
                           * (raw_x[1] - raw_x[0]) * 1000)
        measured_spacings = groove_spacings[np.isfinite(groove_spacings)]
        if len(measured_spacings):
            local_period_nm = np.mean(measured_spacings)
            groove_periods_this_group = measured_spacings.tolist()
        else:
            local_period_nm = period_nm_est
            groove_periods_this_group = []

        # Calculate local periods for each groove
        local_periods = _calculate_local_periods(
            groove_centers, groove_spacings.tolist(), local_period_nm
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
            'groove_centers': groove_centers,
            # How this group's background was actually fitted. A group that
            # fell back is not comparable with one that did not.
            'flatten_report': flat_report,
        })
    
    # ------------ across all groups: periods, statistics, report ------------
    # Check if we got any measurements
    if len(all_blaze_angles) == 0:
        print("Could not extract any blaze angles from any row group!")
        return None
    
    print(f"\n{'='*60}")
    print("ROW-GROUP ANALYSIS SUMMARY")
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
    result['groove_orientation'] = _orientation_dict(prep.pop('orientation'))
    result.update(prep)
    result['mask_rejected_grooves'] = n_mask_rejected
    result['rows_excluded'] = group_info.get('rows_excluded', 0)

    return result


def _analyze_single_file_traditional(filename, show_plots=True, settings=None):
    """
    Traditional analysis: average all rows into a single profile.

    This is the original analysis method, kept for compatibility and comparison.
    """
    
    if settings is None:
        settings = AnalysisSettings.from_config()

    print(f"\n{'='*60}")
    print(f"Analyzing: {filename}")
    print("Mode: TRADITIONAL (single averaged profile)")
    print(f"{'='*60}")
    
    # Load and validate data
    data, scan_x_size, prep = _load_and_validate(filename, settings)
    if data is None:
        return None
    
    # Show 2D image if requested
    if settings.show_2d_image and show_plots:
        _plot_2d_image(data, scan_x_size, filename)
    
    # Extract and flatten profile
    (raw_x, raw_y, flat_y, background, period_nm_est,
     flat_report, measured_y) = _process_profile(data, scan_x_size, show_plots,
                                                 settings)
    
    # Find grooves
    groove_centers, n_edge_rejected = find_groove_positions(raw_x, flat_y, period_nm_est,
                                          prominence_factor=settings.prominence_factor,
                                          distance_factor=settings.distance_factor,
                                          edge_exclusion=settings.edge_exclusion_periods,
                                          return_n_edge_rejected=True)
    groove_centers, n_mask_rejected = reject_grooves_near_nan(
        groove_centers, measured_y,
        GROOVE_REACH_PERIODS * period_nm_est / ((raw_x[1] - raw_x[0]) * 1000))

    edge_note = f" ({n_edge_rejected} rejected: clipped by scan edge)" if n_edge_rejected else ""
    if n_mask_rejected:
        edge_note += f" ({n_mask_rejected} rejected: masked region)"
    print(f"Found {len(groove_centers)} grooves{edge_note}")

    if len(groove_centers) == 0:
        print("No grooves detected!")
        return None
    
    # Calculate groove spacing
    spacings = (groove_spacings_px(groove_centers, measured_y)
                * (raw_x[1] - raw_x[0]) * 1000)
    period_nm, period_std, groove_periods = _calculate_periods(
        raw_x, groove_centers, period_nm_est, profile=measured_y)
    local_periods = _calculate_local_periods(groove_centers, spacings.tolist(),
                                             period_nm)
    
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
    result = _package_results(filename, blaze_angles, stats, period_nm, period_std,
                              groove_periods, quality, all_local_angles,
                              raw_x, flat_y, groove_centers,
                              flat_report=flat_report)
    result['groove_orientation'] = _orientation_dict(prep.pop('orientation'))
    result.update(prep)
    result['mask_rejected_grooves'] = n_mask_rejected
    return result


# ============ HELPER FUNCTIONS ============

def _orientation_dict(report):
    """The groove-orientation report as plain values, or None."""
    return None if report is None else report.as_dict()


# ============ LOADING AND FLATTENING ============

def _load_and_validate(filename, settings):
    """
    Load AFM data, flatten the image, and validate it.

    Image flattening happens here because this is the only point where the 2-D
    array exists and every workflow passes through it; the next step averages
    the rows away. Groove orientation is measured, and corrected when the
    settings say so, here for the same reason: it needs the rows apart.

    Returns:
        (data, scan_x_size, prep), where prep records what was done to the
        scan: the OrientationReport under 'orientation', and the scar repair
        and region mask. (None, None, None) when the file cannot be read.
    """
    try:
        data, scan_x_size = load_afm_data(filename,
                                          default_scan_size=settings.scan_x_size,
                                          settings=settings)

        # The region mask first: levelling reads every pixel of a line, so a
        # particle left in would move every line it sits on. Excluded pixels
        # are NaN from here on; see core/region_mask.py.
        data, n_masked = apply_region_mask(data, settings)
        if n_masked:
            print(f"Region mask: {describe_mask(settings)} "
                  f"({n_masked} pixels excluded)")

        method = getattr(settings, 'image_flatten_method', 'none')
        if method != 'none':
            before = row_offset_spread(data)
            data = flatten_image(data, method)
            print(f"Image flattening: {method} "
                  f"(row-offset spread {before*1e9:.2f} -> "
                  f"{row_offset_spread(data)*1e9:.2f} nm)")

        # After levelling, which would otherwise read every row offset as a
        # scar, and before the orientation correction and tip correction,
        # which would carry each scar along with its row and widen it
        # sideways. See core/scars.py.
        data, scars = repair_scars_for(data, settings)
        if scars.method != 'none':
            print(f"Scar repair: {scars.summary}")

        # Before tip correction: see core/orientation.py for why the tip is
        # undone along the groove normal rather than along the scan line.
        data, scan_x_size, orientation = apply_groove_orientation(
            data, scan_x_size, load_scan_y_size(filename, settings), settings)
        print(f"Groove orientation: {orientation.summary()}")

        data, correction = tip_correct_around_nan(data, scan_x_size, settings)
        if correction is not None:
            # The uncertain pixels are upper bounds, not measurements. Facet
            # fits mostly live on open facets where the apex tracks, but the
            # reader of a blaze angle from a corrected scan should know the
            # scan was corrected and by what tip.
            print(f"Tip correction: {correction.summary}")

        # Recorded with the result: a repaired, masked or rotated scan is a
        # different measurement of the same file, and the summary has to say
        # which. The caller takes 'orientation' out before merging the rest.
        prep = {
            'orientation': orientation,
            'scar_repair': scars.method,
            'scar_lines': scars.n_lines,
            'scar_pixels': scars.n_pixels,
            'region_mask': describe_mask(settings),
            'masked_pixels': n_masked,
        }
        return data, scan_x_size, prep
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None, None, None


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
    # Extract profile. Masked columns are NaN in `measured_y`; the bridged
    # `raw_y` is only for flattening and detection.
    raw_x, measured_y = raw_data(data, scan_x_size)
    raw_y = bridge_nan(measured_y)
    
    # Calculate scan parameters
    scan_width_nm = scan_x_size * 1000
    estimated_grooves = max(2, int(scan_width_nm / settings.period_est))
    period_nm_est = scan_width_nm / estimated_grooves
    
    print(f"Scan width: {scan_width_nm:.1f} nm")
    print(f"Estimated period: {period_nm_est:.2f} nm")
    
    # Flatten profile
    flat_y, background, flat_report = flatten_profile(raw_x, raw_y,
                                         method=settings.flatten_method,
                                         poly_order=settings.flatten_poly_order,
                                         exclude_edges=settings.flatten_exclude_edges,
                                         feature=settings.flatten_feature,
                                         period_nm=period_nm_est,
                                         return_report=True,
                                         exclude=~np.isfinite(measured_y))
    
    # Show flattening diagnostic if requested
    if show_plots and settings.show_flattening_diagnostic:
        plot_flattening_diagnostic(raw_x, raw_y, flat_y, background,
                                  settings.flatten_method, settings.flatten_feature, period_nm_est)
    
    # Print flattening info
    _print_flattening_info(raw_x, background, settings)
    
    return raw_x, raw_y, flat_y, background, period_nm_est, flat_report, measured_y


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


# ============ PERIODS AND BLAZE ANGLES ============

def _calculate_periods(raw_x, groove_centers, period_nm_est, profile=None):
    """
    Calculate groove periods from detected positions

    profile: the averaged profile with masked columns NaN. A spacing across a
        gap is left out, since it is several periods and not one
        (core/region_mask.groove_spacings_px).
    """
    groove_periods = []
    if profile is None:
        spacings = np.diff(groove_centers).astype(float)
    else:
        spacings = groove_spacings_px(groove_centers, profile)
    groove_spacings = (spacings[np.isfinite(spacings)]
                       * (raw_x[1] - raw_x[0]) * 1000)
    
    if len(groove_spacings) > 0:
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
    
    # The spacing on each side, where there is one: none before the first
    # groove, none after the last, and none across a masked gap, which the
    # caller marks NaN. Averaged when there are two; the global value when
    # there are none. With no NaN this is first/last/average as it always was.
    n = len(groove_centers)
    for i in range(n):
        sides = []
        if i > 0 and np.isfinite(groove_periods[i - 1]):
            sides.append(groove_periods[i - 1])
        if i < n - 1 and np.isfinite(groove_periods[i]):
            sides.append(groove_periods[i])
        local_periods.append(sum(sides) / len(sides) if sides else period_nm)
        
    return local_periods


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
