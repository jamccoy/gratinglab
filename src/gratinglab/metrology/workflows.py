"""
Analysis workflow orchestration
Handles the three main analysis modes: single, multiple, compare
"""
import glob
import os
from datetime import datetime

from .config import (
    SINGLE_FILE, FILE_PATTERN, SAMPLES_TO_COMPARE, DATA_DIR, RESULTS_DIR,
    resolve_path, GGP_SOURCE_FILE,
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
    print("Mode: Single file analysis")
    result = analyze_single_file(resolve_path(SINGLE_FILE), show_plots=True)
    
    if result is not None:
        save_results_to_file([result])
    
    return result


def run_multiple_file_analysis():
    """
    Analyze multiple files matching a pattern
    """
    print("Mode: Multiple file analysis (pattern matching)")
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
    print("Mode: Compare specific samples")
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


def run_boundary_profile_export(settings=None, filename=None):
    """
    Average the grooves of one scan into a PCGrate .ggp boundary profile.

    Shares the front-end with the blaze-angle workflows - loading, image
    flattening, groove detection - then diverges: instead of fitting facets, the
    grooves are averaged, normalised to one period and exported for
    grating-efficiency modelling.

    Parameters:
        settings: AnalysisSettings; defaults to config.py, so the CLI behaves as
            it always has. The GUI passes its own.
        filename: overrides GGP_SOURCE_FILE, for the GUI which already has a file
            open.

    The computation lives in boundary.build_boundary_profile so the GUI can
    preview exactly what this writes.
    """
    from .core.processing import load_afm_data
    from .core.image_flatten import flatten_image
    from .boundary import build_boundary_profile
    from gratinglab.io.ggp import write_ggp

    from .io.metrics import write_profile_metrics
    from .settings import AnalysisSettings

    if settings is None:
        settings = AnalysisSettings.from_config()
    filename = resolve_path(filename or GGP_SOURCE_FILE)

    print("Mode: PCGrate boundary profile export")
    print(f"Source: {filename}")

    data, scan_x_size = load_afm_data(filename,
                                      default_scan_size=settings.scan_x_size,
                                      settings=settings)

    # Applied here too, so both outputs treat a scan the same way. Provably a
    # no-op for the affine methods - the endpoint flattening downstream removes
    # constant and linear terms again - and verified identical at the written
    # precision, so the stored .ggp fixture still matches.
    if settings.image_flatten_method != 'none':
        data = flatten_image(data, settings.image_flatten_method)

    try:
        profile = build_boundary_profile(data, scan_x_size, settings)
    except ValueError as exc:
        print(f"{exc}")
        return None

    if profile.metrics.get('tip_correction', 'none') != 'none':
        print(f"Tip correction: erosion "
              f"(R = {profile.metrics['tip_radius_nm']:g} nm, "
              f"half angle = {profile.metrics['tip_half_angle_deg']:g} deg), "
              f"{100.0 * profile.metrics['tip_certain_fraction']:.1f}% of "
              f"pixels certain; the rest are upper bounds")

    edge_note = (f" ({profile.n_edge_rejected} rejected: clipped by scan edge)"
                 if profile.n_edge_rejected else "")
    print(f"Found {profile.n_grooves} grooves{edge_note}")
    print(f"Measured period: {profile.period_nm:.2f} nm")
    print(f"Averaged {profile.n_used} of {profile.n_grooves} grooves")

    stem = os.path.splitext(os.path.basename(filename))[0]
    ggp_path = os.path.join(RESULTS_DIR, f'averaged_groove_profile_{stem}.ggp')
    met_path = os.path.join(RESULTS_DIR, f'groove_analysis_metrics_{stem}.txt')
    write_ggp(ggp_path, t=profile.x_norm, y=profile.y_norm)
    write_profile_metrics(met_path, profile.metrics)

    print(f"\n  Groove depth: {profile.metrics['groove_depth']:.4f} of period")
    print(f"  Max sidewall angle: {profile.metrics['max_angle_deg']:.2f}deg")
    print(f"\n✓ Boundary profile saved to: {ggp_path}")
    print(f"✓ Profile metrics saved to: {met_path}")

    return {'x': profile.x_norm, 'y': profile.y_norm,
            'metrics': profile.metrics, 'ggp_path': ggp_path,
            'profile': profile}


def run_icc_report():
    """
    Measure how correlated the row-group measurements actually are.

    Diagnostic only - reports numbers, changes nothing. Row-group analysis
    re-measures the same grooves in each band of a scan, so the measurements are
    not independent, but the reported SEMs and p-values assume they are. The ICC
    quantifies that, per scan, and the design effect turns it into the factor by
    which the current SEMs are understated.
    """
    from .stats.icc import (compute_icc, effective_sample_size,
                            sem_inflation_factor)

    print("Mode: ICC diagnostic (row-group correlation)")
    print(f"Scans to check: {len(SAMPLES_TO_COMPARE)}\n")

    rows = []
    for sample_info in SAMPLES_TO_COMPARE:
        filename, label = sample_info[0], sample_info[1]
        result = analyze_single_file(resolve_path(filename), show_plots=False)
        if result is None:
            print(f"  {label}: analysis produced no measurements, skipping")
            continue

        groups = result.get('groove_row_groups')
        if not groups:
            print(f"  {label}: no row-group labels "
                  f"(USE_ROW_GROUPS is off?), skipping")
            continue

        stats = compute_icc(result['all_angles'], groups)
        stats['n_eff'] = effective_sample_size(
            stats['n_measurements'], stats['mean_group_size'], stats['icc'])
        stats['sem_factor'] = sem_inflation_factor(
            stats['mean_group_size'], stats['icc'])
        stats['label'] = label
        stats['file'] = os.path.basename(filename)
        stats['sem_reported'] = result.get('sem', float('nan'))
        rows.append(stats)

    if not rows:
        print("No scans could be assessed.")
        return []

    report = _format_icc_report(rows)
    print(report)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(RESULTS_DIR, f'icc_report_{timestamp}.txt')
    with open(path, 'w') as f:
        f.write(report)
    print(f"\n✓ ICC report saved to: {path}")

    return rows


def _format_icc_report(rows):
    """Render the ICC table and its interpretation"""
    from .stats.icc import interpret_icc

    lines = []
    lines.append("=" * 100)
    lines.append("ROW-GROUP CORRELATION (ICC) DIAGNOSTIC")
    lines.append("=" * 100)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("Row-group analysis measures the same physical grooves once per")
    lines.append("horizontal band, so those measurements are not independent. The")
    lines.append("reported SEM divides by sqrt(N) over all of them, which is only")
    lines.append("valid at ICC = 0.")
    lines.append("")
    lines.append(f"{'Sample':<9}{'Scan file':<34}{'N':>5}{'grp':>5}{'ICC':>8}"
                 f"{'N_eff':>8}{'SEM x':>7}  interpretation")
    lines.append("-" * 100)

    for r in rows:
        lines.append(
            f"{r['label']:<9}{r['file']:<34}{r['n_measurements']:>5}"
            f"{r['n_groups']:>5}{r['icc']:>8.3f}{r['n_eff']:>8.1f}"
            f"{r['sem_factor']:>7.2f}  {r['interpretation']}")

    lines.append("-" * 100)

    iccs = [r['icc'] for r in rows]
    factors = [r['sem_factor'] for r in rows]
    lines.append("")
    lines.append(f"ICC range      : {min(iccs):.3f} to {max(iccs):.3f}   "
                 f"(median {sorted(iccs)[len(iccs)//2]:.3f})")
    lines.append(f"SEM inflation  : {min(factors):.2f}x to {max(factors):.2f}x")
    lines.append("")
    lines.append("Reading:")
    lines.append("  ICC    fraction of variance sitting between row groups rather")
    lines.append("         than within them. 0 = independent, 1 = each group adds")
    lines.append("         only one measurement's worth of information.")
    lines.append("  N_eff  measurements the data is actually worth, n / design effect.")
    lines.append("  SEM x  factor the reported standard error should be multiplied")
    lines.append("         by. Confidence intervals scale the same way.")
    lines.append("")

    worst = max(rows, key=lambda r: r['icc'])
    lines.append(f"Overall: {interpret_icc(max(iccs))} "
                 f"(driven by {worst['label']} / {worst['file']}).")
    lines.append("")
    lines.append("This report changes nothing. It sizes the correction that")
    lines.append("row-group analysis applies via stats/icc.py.")
    lines.append("=" * 100)
    return "\n".join(lines)


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
    print("SUMMARY TABLE")
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