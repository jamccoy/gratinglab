"""
File I/O operations for AFM analysis
Handles saving results to text and CSV files
"""
import numpy as np
import os
from datetime import datetime


def save_results_to_file(results, labels=None, temperatures=None, output_dir=None):
    """
    Save analysis results to text and CSV files

    Parameters:
        results: list of result dictionaries
        labels: list of sample labels
        temperatures: list of treatment temperatures (or None)
        output_dir: directory to save results (defaults to the project's results/)
    """
    from ..config import ANALYSIS_MODE, BLAZE_SIDE, FACET_TRIM, FLATTEN_METHOD, RESULTS_DIR

    if output_dir is None:
        output_dir = RESULTS_DIR

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate timestamp for filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    if labels is None:
        labels = [os.path.basename(r['filename']) for r in results]
    
    # ===== DETAILED TEXT SUMMARY =====
    txt_filename = os.path.join(output_dir, f'analysis_summary_{timestamp}.txt')
    
    with open(txt_filename, 'w') as f:
        f.write("="*80 + "\n")
        f.write("AFM BLAZE ANGLE ANALYSIS SUMMARY\n")
        f.write("="*80 + "\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Number of samples: {len(results)}\n")
        f.write(f"Analysis mode: {ANALYSIS_MODE}\n")
        f.write(f"Blaze side: {BLAZE_SIDE}\n")
        f.write(f"Facet trim: {FACET_TRIM*100:.0f}%\n")
        f.write(f"Flattening method: {FLATTEN_METHOD}\n")
        f.write("\n")
        
        # Individual sample details
        f.write("="*80 + "\n")
        f.write("INDIVIDUAL SAMPLE RESULTS\n")
        f.write("="*80 + "\n\n")
        
        for r, label in zip(results, labels):
            f.write(f"Sample: {label}\n")
            
            # Check if this is a combined result
            if 'individual_scans' in r:
                f.write(f"Combined from {r['n_scans']} scans:\n")
                for scan in r['individual_scans']:
                    f.write(f"  - {os.path.basename(scan['filename'])}\n")
            else:
                f.write(f"File: {r['filename']}\n")
            
            f.write(f"{'-'*80}\n")
            f.write(f"  Number of grooves analyzed: {r['n_grooves']}")
            if 'n_scans' in r and r['n_scans'] > 1:
                f.write(f" ({r['n_scans']} scans)")
            f.write("\n")
            
            f.write(f"  Mean blaze angle: {r['mean_angle']:.3f}° ± {r['std_angle']:.3f}°\n")
            f.write(f"  Angle range: {r['min_angle']:.3f}° to {r['max_angle']:.3f}°\n")
            f.write(f"  Mean slope (dy/dx): {r['mean_slope']:.6f}\n")
            if r['mean_steep'] is not None:
                f.write(f"  Mean steep facet angle: {r['mean_steep']:.3f}°\n")
            f.write(f"\n  Groove geometry:\n")
            f.write(f"    Measured period: {r['period_nm']:.3f} nm ± {r.get('period_std', 0):.3f} nm\n")
            f.write(f"    Mean groove depth: {np.mean([q['groove_depth_nm'] for q in r['quality']]):.3f} nm\n")
            f.write(f"    Mean blaze facet width: {np.mean([q['blaze_width_nm'] for q in r['quality']]):.3f} nm\n")
            f.write(f"\n  Within-facet variation (camber/curvature):\n")
            f.write(f"    Local angle std: {r.get('local_angle_std', 0):.3f}°\n")
            f.write(f"    Local angle range: {r.get('local_angle_range', 0):.3f}°\n")
            
            # Individual groove periods if available
            if 'groove_periods' in r and len(r['groove_periods']) > 0:
                f.write(f"\n  Individual groove spacings (nm):\n")
                for i, period in enumerate(r['groove_periods']):
                    f.write(f"    Groove {i+1}-{i+2}: {period:.3f} nm\n")
            
            f.write("\n")
        
        # Comparison statistics
        if len(results) > 1:
            f.write("="*80 + "\n")
            f.write("COMPARISON STATISTICS\n")
            f.write("="*80 + "\n\n")
            
            # Sample-to-sample comparisons
            f.write("Pairwise comparisons:\n")
            f.write(f"{'-'*80}\n")
            for i, r1 in enumerate(results):
                for j, r2 in enumerate(results):
                    if i < j:
                        diff = r2['mean_angle'] - r1['mean_angle']
                        se_combined = np.sqrt(r1['std_angle']**2 / r1['n_grooves'] + 
                                             r2['std_angle']**2 / r2['n_grooves'])
                        
                        f.write(f"\n{labels[i]} vs {labels[j]}:\n")
                        f.write(f"  {labels[i]}: {r1['mean_angle']:.3f}° ± {r1['std_angle']:.3f}°\n")
                        f.write(f"  {labels[j]}: {r2['mean_angle']:.3f}° ± {r2['std_angle']:.3f}°\n")
                        f.write(f"  Difference: {diff:+.3f}° ± {se_combined:.3f}° (SE)\n")
            
            # Temperature-dependent analysis
            if temperatures is not None:
                temp_samples = [(i, t, labels[i], results[i]) for i, t in enumerate(temperatures) if t is not None]
                
                if len(temp_samples) >= 2:
                    f.write(f"\n{'-'*80}\n")
                    f.write("Temperature-dependent changes (consecutive samples):\n")
                    f.write(f"{'-'*80}\n")
                    
                    for idx in range(len(temp_samples) - 1):
                        i, temp_i, label_i, result_i = temp_samples[idx]
                        j, temp_j, label_j, result_j = temp_samples[idx + 1]
                        
                        angle_diff = result_j['mean_angle'] - result_i['mean_angle']
                        temp_diff = temp_j - temp_i
                        
                        if temp_diff != 0:
                            rate = angle_diff / temp_diff
                            se_combined = np.sqrt(result_i['std_angle']**2 / result_i['n_grooves'] + 
                                                 result_j['std_angle']**2 / result_j['n_grooves'])
                            se_rate = se_combined / abs(temp_diff)
                            
                            f.write(f"\n{label_i} ({temp_i}°C) → {label_j} ({temp_j}°C):\n")
                            f.write(f"  Temperature change: {temp_diff:+.0f}°C\n")
                            f.write(f"  Angle change: {angle_diff:+.3f}°\n")
                            f.write(f"  Rate: {rate:+.6f} ± {se_rate:.6f} °/°C\n")
                
                # Master vs treated comparisons
                no_temp_samples = [(i, labels[i], results[i]) for i, t in enumerate(temperatures) if t is None]
                if len(no_temp_samples) > 0 and len(temp_samples) > 0:
                    f.write(f"\n{'-'*80}\n")
                    f.write("Master vs treated samples:\n")
                    f.write(f"{'-'*80}\n")
                    
                    for i, label_i, result_i in no_temp_samples:
                        for j, temp_j, label_j, result_j in temp_samples:
                            diff = result_j['mean_angle'] - result_i['mean_angle']
                            se_combined = np.sqrt(result_i['std_angle']**2 / result_i['n_grooves'] + 
                                                 result_j['std_angle']**2 / result_j['n_grooves'])
                            f.write(f"\n{label_i} → {label_j} ({temp_j}°C):\n")
                            f.write(f"  Difference: {diff:+.3f}° ± {se_combined:.3f}°\n")
    
    print(f"\n✓ Detailed summary saved to: {txt_filename}")
    
    # ===== CSV FILE FOR DATA ANALYSIS =====
    csv_filename = os.path.join(output_dir, f'analysis_data_{timestamp}.csv')
    
    with open(csv_filename, 'w') as f:
        # Header
        f.write("Sample,File,N_scans,N_grooves,Mean_angle_deg,Std_angle_deg,Min_angle_deg,Max_angle_deg,")
        f.write("Mean_slope,Period_nm,Period_std_nm,Mean_depth_nm,Mean_facet_width_nm,")
        f.write("Local_angle_std_deg,Local_angle_range_deg,Temperature_C\n")
        
        # Data rows
        for i, (r, label) in enumerate(zip(results, labels)):
            temp = temperatures[i] if temperatures is not None and i < len(temperatures) else ''
            
            # Get filename(s)
            if 'individual_scans' in r:
                filenames = '; '.join([os.path.basename(s['filename']) for s in r['individual_scans']])
            else:
                filenames = os.path.basename(r['filename'])
            
            n_scans = r.get('n_scans', 1)
            
            f.write(f"{label},{filenames},{n_scans},")
            f.write(f"{r['n_grooves']},{r['mean_angle']:.4f},{r['std_angle']:.4f},")
            f.write(f"{r['min_angle']:.4f},{r['max_angle']:.4f},")
            f.write(f"{r['mean_slope']:.6f},{r['period_nm']:.4f},{r.get('period_std', 0):.4f},")
            f.write(f"{np.mean([q['groove_depth_nm'] for q in r['quality']]):.4f},")
            f.write(f"{np.mean([q['blaze_width_nm'] for q in r['quality']]):.4f},")
            f.write(f"{r.get('local_angle_std', 0):.4f},{r.get('local_angle_range', 0):.4f},")
            f.write(f"{temp}\n")
    
    print(f"✓ CSV data saved to: {csv_filename}")
    
    # ===== PER-GROOVE CSV (detailed data) =====
    detailed_csv = os.path.join(output_dir, f'per_groove_data_{timestamp}.csv')
    
    def _row_group_labels(scan):
        """
        Row group of each measurement, or blanks when not applicable.

        Traditional (non row-group) analysis has no such label, and older result
        dicts predate the field, so fall back to empty strings rather than
        guessing - a wrong group label would silently corrupt any ICC computed
        from this file.
        """
        groups = scan.get('groove_row_groups')
        if groups is None or len(groups) != len(scan['all_angles']):
            return [''] * len(scan['all_angles'])
        return list(groups)

    with open(detailed_csv, 'w') as f:
        f.write("Sample,Scan_file,Row_group,Groove_number,Blaze_angle_deg,Groove_depth_nm,")
        f.write("Blaze_width_nm,Steep_width_nm,R2,Local_period_nm\n")

        for r, label in zip(results, labels):
            # Handle both single and combined results
            scans = r['individual_scans'] if 'individual_scans' in r else [r]
            for scan in scans:
                scan_file = os.path.basename(scan['filename'])
                groups = _row_group_labels(scan)
                for i, (angle, qual, grp) in enumerate(
                        zip(scan['all_angles'], scan['quality'], groups)):
                    f.write(f"{label},{scan_file},{grp},{i+1},{angle:.4f},")
                    f.write(f"{qual['groove_depth_nm']:.4f},{qual['blaze_width_nm']:.4f},")
                    f.write(f"{qual['steep_width_nm']:.4f},{qual['blaze_r2']:.4f},")
                    f.write(f"{qual.get('local_period_used_nm', scan['period_nm']):.4f}\n")
    
    print(f"✓ Per-groove data saved to: {detailed_csv}")
    
    return txt_filename, csv_filename, detailed_csv