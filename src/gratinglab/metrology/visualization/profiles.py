"""
AFM profile plotting functions
Complex multi-panel layouts organized by temperature
"""
import matplotlib.pyplot as plt
import os


def plot_sample_profiles_by_temperature(grouped_results, grouped_labels, temp_order):
    """
    Plot AFM profiles organized by temperature (one figure per temperature)
    
    Parameters:
        grouped_results: dict mapping temperature key to result
        grouped_labels: dict mapping temperature key to label
        temp_order: list of temperature keys in order
    """
    
    for temp_key in temp_order:
        result = grouped_results[temp_key]
        label = grouped_labels[temp_key]
        
        # Check if this is a combined result with multiple scans
        if 'individual_scans' in result:
            # Multiple scans at this temperature
            _plot_multiple_scans(result, label)
        else:
            # Single scan at this temperature
            _plot_single_scan(result, label)


def _plot_multiple_scans(result, label):
    """
    Plot multiple AFM scans on separate subplots (for same temperature)
    
    Parameters:
        result: combined result dictionary with 'individual_scans' key
        label: sample label
    """
    scans = result['individual_scans']
    n_scans = len(scans)
    
    # Create figure with one subplot per scan
    subplot_height = 3
    total_height = min(subplot_height * n_scans, 11)
    fig, axes = plt.subplots(n_scans, 1, figsize=(14, total_height))
    
    if n_scans == 1:
        axes = [axes]
    
    # Overall figure title
    fig.suptitle(f'{label} - Combined: {result["mean_angle"]:.2f}° ± {result["std_angle"]:.2f}° '
                f'(N={result["n_grooves"]} grooves from {n_scans} scans)', 
                fontsize=12, fontweight='bold', y=0.995)
    
    for i, scan in enumerate(scans):
        ax = axes[i]
        
        if 'raw_x' in scan and 'flat_y' in scan:
            # Plot profile
            ax.plot(scan['raw_x'], scan['flat_y'], 'k-', linewidth=0.8, alpha=0.5, 
                   label='Flattened profile')
            
            # Plot groove centers
            if 'groove_centers' in scan:
                ax.plot(scan['raw_x'][scan['groove_centers']], 
                       scan['flat_y'][scan['groove_centers']], 
                       'ro', markersize=4, alpha=0.7, label='Groove centers')
            
            # Plot analyzed regions and fits
            if 'quality' in scan:
                for j, qual in enumerate(scan['quality']):
                    if 'regions' in qual:
                        reg = qual['regions']
                        blaze_x_um = reg['blaze_x_trim'] / 1000
                        blaze_y = reg['blaze_y_trim']
                        fit_line = reg['blaze_fit']
                        
                        if j == 0:
                            ax.plot(blaze_x_um, blaze_y, 'b-', linewidth=2.5, alpha=0.8,
                                   label='Analyzed regions', zorder=3)
                            ax.plot(blaze_x_um, fit_line, 'cyan', linewidth=1.5, 
                                   linestyle='--', alpha=0.7, label='Linear fits', zorder=4)
                        else:
                            ax.plot(blaze_x_um, blaze_y, 'b-', linewidth=2.5, alpha=0.8, zorder=3)
                            ax.plot(blaze_x_um, fit_line, 'cyan', linewidth=1.5, 
                                   linestyle='--', alpha=0.7, zorder=4)
            
            ax.set_ylabel('Height (nm)', fontsize=9)
            
            # Scan-specific title
            scan_title = f'Scan {i+1}: {os.path.basename(scan["filename"])}\n'
            scan_title += f'Mean = {scan["mean_angle"]:.2f}° ± {scan["std_angle"]:.2f}° '
            scan_title += f'(N={scan["n_grooves"]} grooves) | Spacing: {scan["period_nm"]:.2f} nm'
            ax.set_title(scan_title, fontsize=9)
            ax.legend(loc='best', fontsize=7)
            ax.grid(True, alpha=0.3)
            
            if i == n_scans - 1:
                ax.set_xlabel('Position (µm)', fontsize=9)
        else:
            ax.text(0.5, 0.5, f'Scan {i+1}\nProfile data not available', 
                   ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])


def _plot_single_scan(result, label):
    """
    Plot a single AFM scan
    
    Parameters:
        result: result dictionary
        label: sample label
    """
    fig, ax = plt.subplots(1, 1, figsize=(14, 4))
    
    if 'raw_x' in result and 'flat_y' in result:
        # Plot profile
        ax.plot(result['raw_x'], result['flat_y'], 'k-', linewidth=0.8, alpha=0.5, 
               label='Flattened profile')
        
        # Plot groove centers
        if 'groove_centers' in result:
            ax.plot(result['raw_x'][result['groove_centers']], 
                   result['flat_y'][result['groove_centers']], 
                   'ro', markersize=4, alpha=0.7, label='Groove centers')
        
        # Plot analyzed regions and fits
        if 'quality' in result:
            for j, qual in enumerate(result['quality']):
                if 'regions' in qual:
                    reg = qual['regions']
                    blaze_x_um = reg['blaze_x_trim'] / 1000
                    blaze_y = reg['blaze_y_trim']
                    fit_line = reg['blaze_fit']
                    
                    if j == 0:
                        ax.plot(blaze_x_um, blaze_y, 'b-', linewidth=2.5, alpha=0.8,
                               label='Analyzed regions', zorder=3)
                        ax.plot(blaze_x_um, fit_line, 'cyan', linewidth=1.5, 
                               linestyle='--', alpha=0.7, label='Linear fits', zorder=4)
                    else:
                        ax.plot(blaze_x_um, blaze_y, 'b-', linewidth=2.5, alpha=0.8, zorder=3)
                        ax.plot(blaze_x_um, fit_line, 'cyan', linewidth=1.5, 
                               linestyle='--', alpha=0.7, zorder=4)
        
        ax.set_ylabel('Height (nm)', fontsize=10)
        ax.set_xlabel('Position (µm)', fontsize=10)
        
        # Title
        title = f'{label}: Mean = {result["mean_angle"]:.2f}° ± {result["std_angle"]:.2f}° '
        title += f'(N={result["n_grooves"]} grooves) | Spacing: {result["period_nm"]:.2f} nm'
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, f'{label}\nProfile data not available', 
               ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()


# Legacy function for backwards compatibility
def plot_sample_profiles_comparison(results, labels):
    """
    DEPRECATED: Use plot_sample_profiles_by_temperature instead
    
    Plot AFM profiles for all samples in a single multi-panel figure
    This is kept for backwards compatibility but not recommended for many samples
    """
    import warnings
    warnings.warn("plot_sample_profiles_comparison is deprecated. "
                 "Use plot_sample_profiles_by_temperature instead.", 
                 DeprecationWarning)
    
    n_samples = len(results)
    
    if n_samples <= 3:
        subplot_height = 3
    elif n_samples == 4:
        subplot_height = 2.5
    else:
        subplot_height = 2.0
    
    total_height = min(subplot_height * n_samples, 11)
    fig, axes = plt.subplots(n_samples, 1, figsize=(14, total_height))
    
    if n_samples == 1:
        axes = [axes]
    
    for i, (r, label) in enumerate(zip(results, labels)):
        if 'raw_x' in r and 'flat_y' in r:
            axes[i].plot(r['raw_x'], r['flat_y'], 'k-', linewidth=0.8, alpha=0.5, 
                        label='Flattened profile')
            
            if 'groove_centers' in r:
                axes[i].plot(r['raw_x'][r['groove_centers']], r['flat_y'][r['groove_centers']], 
                           'ro', markersize=4, alpha=0.7, label='Groove centers')
            
            if 'quality' in r:
                for j, qual in enumerate(r['quality']):
                    if 'regions' in qual:
                        reg = qual['regions']
                        blaze_x_um = reg['blaze_x_trim'] / 1000
                        blaze_y = reg['blaze_y_trim']
                        fit_line = reg['blaze_fit']
                        
                        if j == 0:
                            axes[i].plot(blaze_x_um, blaze_y, 'b-', linewidth=2.5, alpha=0.8,
                                       label='Analyzed regions', zorder=3)
                            axes[i].plot(blaze_x_um, fit_line, 'cyan', linewidth=1.5, 
                                       linestyle='--', alpha=0.7, label='Linear fits', zorder=4)
                        else:
                            axes[i].plot(blaze_x_um, blaze_y, 'b-', linewidth=2.5, alpha=0.8, zorder=3)
                            axes[i].plot(blaze_x_um, fit_line, 'cyan', linewidth=1.5, 
                                       linestyle='--', alpha=0.7, zorder=4)
            
            axes[i].set_ylabel('Height (nm)', fontsize=9)
            
            # Include spacing in the title
            spacing_str = f"Spacing: {r['period_nm']:.2f} nm"
            axes[i].set_title(f'{label}: Mean = {r["mean_angle"]:.2f}° ± {r["std_angle"]:.2f}° (N={r["n_grooves"]}) | {spacing_str}', 
                            fontsize=10)
            axes[i].legend(loc='best', fontsize=7)
            axes[i].grid(True, alpha=0.3)
            
            if i == n_samples - 1:
                axes[i].set_xlabel('Position (µm)', fontsize=9)
        else:
            axes[i].text(0.5, 0.5, f'{label}\nProfile data not available', 
                       ha='center', va='center', transform=axes[i].transAxes)
    
    plt.tight_layout()