"""
Diagnostic plotting functions
Shows analysis details like flattening and analyzed regions
"""
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib.lines import Line2D


def plot_analyzed_regions_overlay(raw_x, flat_y, groove_centers, period_nm, quality_list, facet_trim):
    """
    Create an overlay plot showing which portions of each groove are actually analyzed.
    
    Parameters:
        raw_x: x-axis positions (µm)
        flat_y: flattened height profile (nm)
        groove_centers: indices of detected groove centers
        period_nm: groove period (nm)
        quality_list: list of quality dictionaries with region info
        facet_trim: fraction of facet trimmed from edges
    """
    fig, ax = plt.subplots(figsize=(14, 6))
    
    ax.plot(raw_x, flat_y, 'k-', linewidth=1, alpha=0.4, label='Flattened profile', zorder=1)
    ax.plot(raw_x[groove_centers], flat_y[groove_centers], 'ro', 
            markersize=6, label='Groove centers', zorder=3, alpha=0.7)
    
    for i, qual in enumerate(quality_list):
        if 'regions' in qual:
            reg = qual['regions']
            blaze_x_um = reg['blaze_x_trim'] / 1000
            blaze_y = reg['blaze_y_trim']
            ax.plot(blaze_x_um, blaze_y, 'b-', linewidth=3, alpha=0.8, zorder=2)
            
            blaze_x_full_um = reg['blaze_x_full'] / 1000
            blaze_y_full = reg['blaze_y_full']
            ax.plot(blaze_x_full_um, blaze_y_full, 'b-', linewidth=1.5, alpha=0.3, zorder=1)
            
            if len(blaze_x_um) > 0:
                fit_line = reg['blaze_fit']
                ax.plot(blaze_x_um, fit_line, 'cyan', linewidth=2, 
                       linestyle='--', alpha=0.7, zorder=2)
    
    custom_lines = [
        Line2D([0], [0], color='k', linewidth=1, alpha=0.4),
        Line2D([0], [0], color='b', linewidth=1.5, alpha=0.3),
        Line2D([0], [0], color='b', linewidth=3, alpha=0.8),
        Line2D([0], [0], color='cyan', linewidth=2, linestyle='--', alpha=0.7),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='r', markersize=6)
    ]
    ax.legend(custom_lines, 
             ['Full profile (flattened)', 'Full blaze facets', 
              'Analyzed regions (trimmed)', 'Linear fits', 'Groove centers'],
             loc='best', fontsize=10)
    
    ax.set_xlabel('Position (µm)', fontsize=12)
    ax.set_ylabel('Height (nm)', fontsize=12)
    ax.set_title(f'Analyzed Regions Overlay - Trim fraction: {facet_trim:.0%} (2.5x extra trough trimming)', fontsize=13)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()


def plot_flattening_diagnostic(raw_x, raw_y, flat_y, background, flatten_method, 
                               flatten_feature, period_nm_est):
    """
    Show before/after flattening comparison
    
    Parameters:
        raw_x: x-axis positions (µm)
        raw_y: raw height profile (nm)
        flat_y: flattened height profile (nm)
        background: fitted background curve (nm)
        flatten_method: method used for flattening
        flatten_feature: feature used (for level_grooves method)
        period_nm_est: estimated period (nm)
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))
    
    ax1.plot(raw_x, raw_y, 'b-', linewidth=1, label='Raw profile', alpha=0.7)
    ax1.plot(raw_x, background, 'r--', linewidth=2, label=f'Background ({flatten_method})')
    
    if flatten_method == 'level_grooves':
        if flatten_feature in ['peaks', 'both']:
            dx_nm = (raw_x[1] - raw_x[0]) * 1000
            min_distance = int(0.3 * period_nm_est / dx_nm)
            peaks, _ = find_peaks(raw_y, distance=min_distance)
            ax1.plot(raw_x[peaks], raw_y[peaks], 'go', markersize=8, 
                    label=f'Peaks used ({len(peaks)})', zorder=5)
        if flatten_feature in ['troughs', 'both']:
            y_inv = -raw_y
            dx_nm = (raw_x[1] - raw_x[0]) * 1000
            min_distance = int(0.3 * period_nm_est / dx_nm)
            troughs, _ = find_peaks(y_inv, distance=min_distance)
            ax1.plot(raw_x[troughs], raw_y[troughs], 'mo', markersize=8, 
                    label=f'Troughs used ({len(troughs)})', zorder=5)
    
    ax1.set_xlabel('Displacement (µm)', fontsize=11)
    ax1.set_ylabel('Height (nm)', fontsize=11)
    ax1.set_title('Before Flattening: Raw Profile + Background Fit', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(raw_x, flat_y, 'g-', linewidth=1, label='Flattened profile')
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    ax2.set_xlabel('Displacement (µm)', fontsize=11)
    ax2.set_ylabel('Height (nm)', fontsize=11)
    ax2.set_title('After Flattening: Background Removed', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()