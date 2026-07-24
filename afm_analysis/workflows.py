"""
Analysis workflow orchestration
Handles the three main analysis modes: single, multiple, compare
"""
import glob
import os

from .config import SINGLE_FILE, FILE_PATTERN, SAMPLES_TO_COMPARE, DATA_DIR, resolve_path
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