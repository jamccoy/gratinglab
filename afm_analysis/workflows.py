"""
Analysis workflow orchestration
Handles the three main analysis modes: single, multiple, compare
"""
import glob
import os

import numpy as np

from .config import (
    SINGLE_FILE, FILE_PATTERN, SAMPLES_TO_COMPARE, DATA_DIR, RESULTS_DIR,
    resolve_path, SCAN_X_SIZE, PERIOD_EST, PROMINENCE_FACTOR, DISTANCE_FACTOR,
    EDGE_EXCLUSION_PERIODS, GGP_SOURCE_FILE, GGP_N_POINTS, GGP_APPLY_SMOOTHING,
    GGP_SMOOTHING_WINDOW, GGP_MIN_HALF_WIDTH,
)
from .analyzer import analyze_single_file
from .io.file_io import save_results_to_file
from .data.aggregation import group_by_temperature, extract_temperatures_for_output
from .stats.analysis import (
    print_comparison_summary, print_pairwise_comparisons,
    print_temperature_analysis
)
from .visualization.statistics import plot_multi_file_comparison
from .visualization.profiles import plot_sample_profiles_by_temperature


def run_single_file_analysis():
    """
    Analyze a single AFM file with full diagnostic plots
    """
    print(f"Mode: Single file analysis")
    result = analyze_single_file(resolve_path(SINGLE_FILE), show_plots=True)
    
    if result is not None:
        save_results_to_file([result])
    
    return result


def run_multiple_file_analysis():
    """
    Analyze multiple files matching a pattern
    """
    print(f"Mode: Multiple file analysis (pattern matching)")
    # Patterns are matched inside the data directory, not the working directory
    files = sorted(glob.glob(os.path.join(DATA_DIR, FILE_PATTERN)))

    if len(files) == 0:
        print(f"No files found in {DATA_DIR} matching pattern: {FILE_PATTERN}")
        return []
    
    print(f"Found {len(files)} files to analyze")
    
    # Analyze all files
    results = []
    for f in files:
        result = analyze_single_file(f, show_plots=False)
        if result is not None:
            results.append(result)
    
    if len(results) == 0:
        return []
    
    # Generate comparison plots
    plot_multi_file_comparison(results)
    
    # Print summary table
    _print_multiple_file_summary(results)
    
    # Save results
    save_results_to_file(results)
    
    return results


def run_comparison_analysis():
    """
    Compare specific samples with optional temperature grouping
    """
    print(f"Mode: Compare specific samples")
    print(f"Number of sample entries: {len(SAMPLES_TO_COMPARE)}")
    
    # Parse sample information and analyze files
    results, labels, temperatures = _analyze_comparison_samples()
    
    if len(results) == 0:
        print("No samples were successfully analyzed!")
        return []
    
    # Group by temperature and combine multiple scans
    grouped_results, grouped_labels, temp_order = group_by_temperature(
        results, labels, temperatures
    )
    
    # Convert to lists for compatibility
    combined_results = [grouped_results[key] for key in temp_order]
    combined_labels = [grouped_labels[key] for key in temp_order]
    combined_temps = extract_temperatures_for_output(grouped_results, temp_order)
    
    # Update labels in results
    for result, label in zip(combined_results, combined_labels):
        result['label'] = label
    
    # Generate all plots and statistics
    _generate_comparison_outputs(combined_results, combined_labels, combined_temps,
                                 grouped_results, grouped_labels, temp_order)
    
    return combined_results


def run_boundary_profile_export():
    """
    Average the grooves of one scan into a PCGrate .ggp boundary profile.

    Shares the profile front-end with the blaze-angle workflows, then diverges:
    instead of fitting facets, the grooves are averaged, normalised to one
    period and exported for grating-efficiency modelling.
    """
    from .core.processing import load_afm_data, raw_data, find_groove_positions
    from .boundary import (flatten_endpoints, average_grooves,
                           normalize_profile, profile_metrics)
    from .io.ggp import write_ggp, write_profile_metrics

    filename = resolve_path(GGP_SOURCE_FILE)
    print(f"Mode: PCGrate boundary profile export")
    print(f"Source: {filename}")

    data, scan_x_size = load_afm_data(filename, default_scan_size=SCAN_X_SIZE)
    raw_x, raw_y = raw_data(data, scan_x_size)

    # Boundary-profile flattening: endpoints to zero, then remove residual tilt.
    # Deliberately not core.flatten_profile - see boundary/average.py.
    flat_y = flatten_endpoints(raw_x, raw_y)
    flat_y = flat_y - np.polyval(np.polyfit(raw_x, flat_y, 1), raw_x)

    scan_width_nm = scan_x_size * 1000
    period_nm = scan_width_nm / max(2, int(scan_width_nm / PERIOD_EST))

    groove_centers, n_edge = find_groove_positions(
        raw_x, flat_y, period_nm,
        prominence_factor=PROMINENCE_FACTOR,
        distance_factor=DISTANCE_FACTOR,
        edge_exclusion=EDGE_EXCLUSION_PERIODS,
        return_n_edge_rejected=True)

    if len(groove_centers) == 0:
        print("No grooves detected!")
        return None

    edge_note = f" ({n_edge} rejected: clipped by scan edge)" if n_edge else ""
    print(f"Found {len(groove_centers)} grooves{edge_note}")

    if len(groove_centers) > 1:
        period_nm = np.mean(np.diff(groove_centers) * (raw_x[1] - raw_x[0]) * 1000)
    print(f"Measured period: {period_nm:.2f} nm")

    x_avg, y_avg, y_std, n_used = average_grooves(
        raw_x, flat_y, groove_centers, period_nm,
        margin=0.0, n_points=GGP_N_POINTS, min_half_width=GGP_MIN_HALF_WIDTH)
    print(f"Averaged {n_used} of {len(groove_centers)} grooves")

    x_norm, y_norm, edge_height = normalize_profile(
        x_avg, y_avg, period_nm,
        apply_smoothing=GGP_APPLY_SMOOTHING,
        smoothing_window=GGP_SMOOTHING_WINDOW)

    metrics = profile_metrics(x_norm, y_norm, period_nm, n_used)

    stem = os.path.splitext(os.path.basename(filename))[0]
    ggp_path = os.path.join(RESULTS_DIR, f'averaged_groove_profile_{stem}.ggp')
    met_path = os.path.join(RESULTS_DIR, f'groove_analysis_metrics_{stem}.txt')
    write_ggp(ggp_path, x_norm, y_norm)
    write_profile_metrics(met_path, metrics)

    print(f"\n  Groove depth: {metrics['groove_depth']:.4f} of period")
    print(f"  Max sidewall angle: {metrics['max_angle_deg']:.2f}deg")
    print(f"\n✓ Boundary profile saved to: {ggp_path}")
    print(f"✓ Profile metrics saved to: {met_path}")

    return {'x': x_norm, 'y': y_norm, 'metrics': metrics, 'ggp_path': ggp_path}


def _analyze_comparison_samples():
    """
    Parse SAMPLES_TO_COMPARE and analyze each file
    
    Returns:
        results: list of analysis results
        labels: list of sample labels
        temperatures: list of temperatures (None for master)
    """
    results = []
    labels = []
    temperatures = []
    
    for sample_info in SAMPLES_TO_COMPARE:
        # Handle both (filename, label) and (filename, label, temperature) formats
        if len(sample_info) == 2:
            filename, label = sample_info
            temperature = None
        elif len(sample_info) == 3:
            filename, label, temperature = sample_info
        else:
            print(f"Error: Invalid sample format: {sample_info}")
            continue
        
        print(f"\nProcessing: {label} ({filename})")
        result = analyze_single_file(resolve_path(filename), show_plots=False)
        
        if result is not None:
            results.append(result)
            labels.append(label)
            temperatures.append(temperature)
            # Store metadata in result
            result['label'] = label
            result['temperature'] = temperature
    
    return results, labels, temperatures


def _print_multiple_file_summary(results):
    """Print summary table for multiple file analysis"""
    print(f"\n{'='*80}")
    print(f"SUMMARY TABLE")
    print(f"{'='*80}")
    print(f"{'File':<40} {'N':<5} {'Mean':<10} {'Std':<10} {'Period(nm)':<12}")
    print(f"{'-'*80}")
    
    for r in results:
        print(f"{os.path.basename(r['filename']):<40} "
              f"{r['n_grooves']:<5} "
              f"{r['mean_angle']:<10.2f} "
              f"{r['std_angle']:<10.2f} "
              f"{r['period_nm']:<12.2f}")


def _generate_comparison_outputs(combined_results, combined_labels, combined_temps,
                                 grouped_results, grouped_labels, temp_order):
    """
    Generate all comparison plots and statistics
    
    Parameters:
        combined_results: list of combined results
        combined_labels: list of labels
        combined_temps: list of temperatures
        grouped_results: dict of results by temperature key
        grouped_labels: dict of labels by temperature key
        temp_order: ordered list of temperature keys
    """
    # Comparison plots with custom labels
    plot_multi_file_comparison(combined_results, labels=combined_labels,
                              temperatures=combined_temps)
    
    # AFM profiles organized by temperature
    plot_sample_profiles_by_temperature(grouped_results, grouped_labels, temp_order)
    
    # Print comparison summary table
    print_comparison_summary(combined_results, combined_labels)
    
    # Statistical comparison
    print_pairwise_comparisons(combined_results)
    
    # Temperature analysis if applicable
    if any(t is not None for t in combined_temps):
        print_temperature_analysis(combined_results, combined_labels, combined_temps)
    
    # Save all results to files
    save_results_to_file(combined_results, labels=combined_labels,
                        temperatures=combined_temps)