"""
Statistical plotting functions
Histograms, distributions, and multi-file comparisons

UPDATED: Proper error bar representation with clear labels
"""
import numpy as np
import matplotlib.pyplot as plt
import os


def plot_summary_statistics(blaze_angles, positions, all_local_angles,
                           mean_angle, std_angle, show_local_distribution):
    """
    Plot histogram, position, and local angle distribution for a single sample

    Parameters:
        blaze_angles: list of blaze angles (degrees)
        positions: x position (µm) of each measurement in blaze_angles, same
            length. Under row-group analysis many measurements share a position,
            because each row group re-measures the same grooves - this used to
            take raw_x and groove_centers and slice them to match, which cannot
            work when there are n_groups x n_grooves angles but only n_grooves
            centers.
        all_local_angles: list of all local angle measurements
        mean_angle: mean blaze angle (degrees)
        std_angle: standard deviation of blaze angles (degrees)
        show_local_distribution: whether to show within-facet distribution
    """
    blaze_angles = np.asarray(blaze_angles)
    positions = np.asarray(positions)
    if len(positions) != len(blaze_angles):
        raise ValueError(f"positions and blaze_angles must be the same length, "
                         f"got {len(positions)} and {len(blaze_angles)}")

    n_plots = 3 if show_local_distribution and len(all_local_angles) > 0 else 2
    fig, axes = plt.subplots(1, n_plots, figsize=(6*n_plots, 4))
    
    if n_plots == 2:
        ax1, ax2 = axes
    else:
        ax1, ax2, ax3 = axes
    
    # Per-groove angle histogram
    ax1.hist(blaze_angles, bins=15, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(mean_angle, color='r', linestyle='--', linewidth=2,
               label=f'Mean: {mean_angle:.2f} deg')
    ax1.set_xlabel('Blaze Angle (degrees)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(f'Per-Groove Angle Distribution\n(σ = {std_angle:.2f}° - Physical Groove Variation)', 
                 fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Angle vs position. With row groups the same groove is measured once per
    # group, so positions repeat - connecting them with a line would draw
    # meaningless zig-zags between groups. Scatter in that case.
    repeated = len(np.unique(positions)) < len(positions)
    if repeated:
        ax2.plot(positions, blaze_angles, 'o', markersize=5, alpha=0.4,
                 color='steelblue')
    else:
        ax2.plot(positions, blaze_angles, 'o-', markersize=8, linewidth=1.5,
                 color='steelblue')
    ax2.axhline(mean_angle, color='r', linestyle='--', linewidth=2, label='Mean')
    ax2.fill_between([np.min(positions), np.max(positions)],
                    mean_angle - std_angle, mean_angle + std_angle,
                    alpha=0.2, color='red', label=f'±1σ ({std_angle:.2f}°)')
    ax2.set_xlabel('Position (µm)', fontsize=12)
    ax2.set_ylabel('Blaze Angle (degrees)', fontsize=12)
    ax2.set_title('Blaze Angle vs Position', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Within-facet angle distribution
    if show_local_distribution and len(all_local_angles) > 0:
        mean_local = np.mean(all_local_angles)
        std_local = np.std(all_local_angles)
        range_local = np.max(all_local_angles) - np.min(all_local_angles)
        
        ax3.hist(all_local_angles, bins=50, edgecolor='black', alpha=0.7, color='seagreen')
        ax3.axvline(mean_local, color='r', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_local:.2f} deg')
        ax3.set_xlabel('Local Angle (degrees)', fontsize=12)
        ax3.set_ylabel('Count', fontsize=12)
        ax3.set_title(f'Within-Facet Angle Distribution\n(σ = {std_local:.2f}° - Facet Curvature, range = {range_local:.2f}°)', 
                     fontsize=13)
        ax3.legend(fontsize=11)
        ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()


def plot_multi_file_comparison(results, labels=None, temperatures=None):
    """
    Comparison plots for multiple files - bar chart and histograms
    
    NOW WITH IMPROVED ERROR BAR HANDLING:
    - Bar charts show ±1 SEM (uncertainty in the mean)
    - Histograms show per-groove distributions (consistent with bar chart)
    - Clear labels distinguish between σ (spread) and SEM (uncertainty)
    
    Parameters:
        results: list of result dictionaries
        labels: list of sample labels (optional)
        temperatures: list of temperatures (optional)
    """
    
    if labels is None:
        labels = [os.path.basename(r['filename']) for r in results]
    
    means = [r['mean_angle'] for r in results]
    stds = [r['std_angle'] for r in results]
    
    # Calculate SEM (standard error of the mean) for error bars
    # SEM = (total_std) / sqrt(n_grooves)
    # Use total_std if available (includes measurement uncertainty), otherwise fall back to std_angle
    sems = []
    total_stds = []  # Store for labels
    
    for r in results:
        if 'sem' in r:
            # Use pre-calculated SEM that includes both measurement and physical uncertainty
            sems.append(r['sem'])
            total_stds.append(r.get('total_std', r['std_angle']))
        elif 'total_std' in r:
            # Calculate SEM from total_std
            total_std = r['total_std']
            sems.append(total_std / np.sqrt(r['n_grooves']))
            total_stds.append(total_std)
        else:
            # Fall back to std/sqrt(n) if new metrics not available
            sems.append(r['std_angle'] / np.sqrt(r['n_grooves']))
            total_stds.append(r['std_angle'])
    
    x_pos = np.arange(len(results))
    
    # Determine subplot layout based on number of samples
    n_samples = len(results)
    
    if n_samples <= 3:
        # 1 row for bar chart + 1 row with histograms side by side
        fig = plt.figure(figsize=(14, 8))
        gs = fig.add_gridspec(2, n_samples, height_ratios=[1, 1])
        
        # Bar chart spans full width of top row
        ax_bar = fig.add_subplot(gs[0, :])
        
        # Histograms in bottom row
        ax_hists = [fig.add_subplot(gs[1, i]) for i in range(n_samples)]
        
    elif n_samples <= 6:
        # 1 row for bar chart + 2 rows of histograms (3 per row)
        n_cols = 3
        n_rows_hist = int(np.ceil(n_samples / n_cols))
        
        fig = plt.figure(figsize=(14, 4 + 3*n_rows_hist))
        gs = fig.add_gridspec(1 + n_rows_hist, n_cols, height_ratios=[1.2] + [1]*n_rows_hist)
        
        # Bar chart spans full width of top row
        ax_bar = fig.add_subplot(gs[0, :])
        
        # Histograms in subsequent rows
        ax_hists = []
        for i in range(n_samples):
            row = 1 + i // n_cols
            col = i % n_cols
            ax_hists.append(fig.add_subplot(gs[row, col]))
    else:
        # Many samples - use 4 columns for histograms
        n_cols = 4
        n_rows_hist = int(np.ceil(n_samples / n_cols))
        
        fig = plt.figure(figsize=(16, 4 + 2.5*n_rows_hist))
        gs = fig.add_gridspec(1 + n_rows_hist, n_cols, height_ratios=[1.2] + [1]*n_rows_hist)
        
        ax_bar = fig.add_subplot(gs[0, :])
        
        ax_hists = []
        for i in range(n_samples):
            row = 1 + i // n_cols
            col = i % n_cols
            ax_hists.append(fig.add_subplot(gs[row, col]))
    
    # Determine colors: master samples get orange, others get steelblue
    bar_colors = []
    for i, (label, temp) in enumerate(zip(labels, temperatures if temperatures else [None]*len(labels))):
        if temp is None or 'master' in label.lower():
            bar_colors.append('darkorange')  # Master samples
        else:
            bar_colors.append('steelblue')   # Treated samples
    
    # Bar chart with SEM error bars
    bars = ax_bar.bar(x_pos, means, yerr=sems, capsize=5, alpha=0.7, 
                      edgecolor='black', color=bar_colors, width=0.6)
    
    # Add headroom for labels (15% extra space at top)
    y_max = max([m + s for m, s in zip(means, sems)])
    ax_bar.set_ylim(0, y_max * 1.15)
    
    ax_bar.set_ylabel('Blaze Angle (degrees)', fontsize=12)
    ax_bar.set_title('Mean Blaze Angle Comparison\n(Error bars: ±1 SEM - Uncertainty in the Mean)', 
                     fontsize=13, fontweight='bold')
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(labels, rotation=45, ha='right', fontsize=10)
    ax_bar.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    # Show: Mean ± SEM (with total_std for reference)
    for i, (bar, mean, sem, std, total_std) in enumerate(zip(bars, means, sems, stds, total_stds)):
        height = bar.get_height()
        # Clarify what each number represents
        ax_bar.text(bar.get_x() + bar.get_width()/2., height + sem + 0.2,
                    f'{mean:.2f}°\n±{sem:.3f}° SEM\n(σ_total={total_std:.2f}°)', 
                    ha='center', va='bottom', 
                    fontsize=8, fontweight='bold')
    
    # Individual histograms for each sample
    # NOW USING PER-GROOVE ANGLES (consistent with bar chart)
    colors = plt.cm.Set3(np.linspace(0, 1, n_samples))
    
    for i, (r, label, ax_hist) in enumerate(zip(results, labels, ax_hists)):
        # IMPORTANT CHANGE: Use per-groove angles (same as bar chart)
        # This makes histograms and bar chart statistically consistent
        angles = r['all_angles']  # These are the per-groove measurements
        
        # But we can show local angles as well if desired (see below)
        has_local = 'all_local_angles' in r and len(r['all_local_angles']) > 0
        
        # Create histogram of per-groove angles
        n_bins = min(30, max(8, len(angles) // 3))
        ax_hist.hist(angles, bins=n_bins, edgecolor='black', alpha=0.7, 
                     color=colors[i], linewidth=1.0, label='Per-groove')
        
        # Add mean line
        mean_val = np.mean(angles)
        std_val = np.std(angles)
        sem_val = sems[i]
        
        ax_hist.axvline(mean_val, color='red', linestyle='--', 
                       linewidth=2, label=f"Mean: {mean_val:.2f}°")
        
        # Add ±1 SEM shading (uncertainty in the mean)
        ax_hist.axvspan(mean_val - sem_val, 
                       mean_val + sem_val,
                       alpha=0.3, color='red', label=f'±1 SEM ({sem_val:.3f}°)')
        
        # Add ±1σ for reference (shows groove-to-groove variation)
        ax_hist.axvspan(mean_val - std_val, 
                       mean_val + std_val,
                       alpha=0.15, color='blue', label=f'±1σ ({std_val:.2f}°)')
        
        # Labels and title
        ax_hist.set_xlabel('Blaze Angle (°)', fontsize=10)
        ax_hist.set_ylabel('Count', fontsize=10)
        
        # Title with sample info
        title = f"{label}\n"
        title += f"{mean_val:.2f}° ± {sem_val:.3f}° SEM "
        
        # Include number of scans if available
        if 'n_scans' in r and r['n_scans'] > 1:
            title += f"({r['n_scans']} scans)"
        
        # Show sample size
        title += f"\nN={len(angles)} measurements"
        
        # Add info about row-group analysis if used
        if 'n_groups' in r and r['n_groups'] > 1:
            title += f" from {r['n_groups']} regions"
        
        ax_hist.set_title(title, fontsize=9, fontweight='bold')
        
        ax_hist.legend(fontsize=7, loc='upper right')
        ax_hist.grid(True, alpha=0.3, axis='y')
    
    # Hide any unused subplot axes
    if n_samples < len(ax_hists):
        for ax in ax_hists[n_samples:]:
            ax.set_visible(False)
    
    plt.tight_layout()


def plot_multi_file_comparison_with_local_angles(results, labels=None, temperatures=None):
    """
    Alternative comparison plot that shows BOTH per-groove AND local angle distributions
    
    This gives you the full picture:
    - Per-groove angles: Shows groove-to-groove variation
    - Local angles: Shows within-facet variation (curvature/camber)
    
    Use this if you want to see both types of variation!
    """
    
    if labels is None:
        labels = [os.path.basename(r['filename']) for r in results]
    
    means = [r['mean_angle'] for r in results]
    stds = [r['std_angle'] for r in results]
    
    # Calculate SEM for error bars
    sems = []
    for r in results:
        if 'sem' in r:
            sems.append(r['sem'])
        elif 'total_std' in r:
            sems.append(r['total_std'] / np.sqrt(r['n_grooves']))
        else:
            sems.append(r['std_angle'] / np.sqrt(r['n_grooves']))
    
    n_samples = len(results)
    
    # Layout: Bar chart + 2 histograms per sample (per-groove + local)
    if n_samples <= 2:
        n_cols = n_samples
        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(3, n_cols, height_ratios=[1, 1, 1])
        ax_bar = fig.add_subplot(gs[0, :])
    else:
        # For many samples, use grid layout
        n_cols = min(3, n_samples)
        n_rows_hist = int(np.ceil(n_samples / n_cols))
        fig = plt.figure(figsize=(16, 4 + 4*n_rows_hist))
        gs = fig.add_gridspec(1 + 2*n_rows_hist, n_cols, 
                            height_ratios=[1.2] + [1, 1]*n_rows_hist)
        ax_bar = fig.add_subplot(gs[0, :])
    
    x_pos = np.arange(len(results))
    
    # Determine colors
    bar_colors = []
    for label, temp in zip(labels, temperatures if temperatures else [None]*len(labels)):
        if temp is None or 'master' in label.lower():
            bar_colors.append('darkorange')
        else:
            bar_colors.append('steelblue')
    
    # Bar chart
    bars = ax_bar.bar(x_pos, means, yerr=sems, capsize=5, alpha=0.7, 
                      edgecolor='black', color=bar_colors, width=0.6)
    
    y_max = max([m + s for m, s in zip(means, sems)])
    ax_bar.set_ylim(0, y_max * 1.15)
    
    ax_bar.set_ylabel('Blaze Angle (degrees)', fontsize=12)
    ax_bar.set_title('Mean Blaze Angle Comparison (Error bars: ±1 SEM)', 
                     fontsize=13, fontweight='bold')
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(labels, rotation=45, ha='right', fontsize=10)
    ax_bar.grid(True, alpha=0.3, axis='y')
    
    for bar, mean, sem in zip(bars, means, sems):
        height = bar.get_height()
        ax_bar.text(bar.get_x() + bar.get_width()/2., height + sem + 0.2,
                    f'{mean:.2f}° ± {sem:.3f}°', 
                    ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    # Dual histograms for each sample
    colors = plt.cm.Set3(np.linspace(0, 1, n_samples))
    
    for i, (r, label) in enumerate(zip(results, labels)):
        if n_samples <= 2:
            row_per_groove = 1
            row_local = 2
            col = i
        else:
            row_base = 1 + 2 * (i // n_cols)
            row_per_groove = row_base
            row_local = row_base + 1
            col = i % n_cols
        
        ax_per_groove = fig.add_subplot(gs[row_per_groove, col])
        ax_local = fig.add_subplot(gs[row_local, col])
        
        # Per-groove histogram
        angles = r['all_angles']
        mean_val = np.mean(angles)
        std_val = np.std(angles)
        sem_val = sems[i]
        
        n_bins = min(20, max(8, len(angles) // 3))
        ax_per_groove.hist(angles, bins=n_bins, edgecolor='black', alpha=0.7, 
                          color=colors[i], linewidth=1.0)
        ax_per_groove.axvline(mean_val, color='red', linestyle='--', linewidth=2)
        ax_per_groove.axvspan(mean_val - sem_val, mean_val + sem_val,
                             alpha=0.3, color='red')
        
        ax_per_groove.set_ylabel('Count', fontsize=9)
        ax_per_groove.set_title(f'{label}\nPer-Groove: {mean_val:.2f}° ± {sem_val:.3f}° SEM (N={len(angles)})',
                               fontsize=9, fontweight='bold')
        ax_per_groove.grid(True, alpha=0.3)
        
        # Local angle histogram
        if 'all_local_angles' in r and len(r['all_local_angles']) > 0:
            local_angles = r['all_local_angles']
            mean_local = np.mean(local_angles)
            std_local = np.std(local_angles)
            
            n_bins_local = min(40, max(15, len(local_angles) // 20))
            ax_local.hist(local_angles, bins=n_bins_local, edgecolor='black', 
                         alpha=0.7, color='seagreen', linewidth=1.0)
            ax_local.axvline(mean_local, color='red', linestyle='--', linewidth=2)
            ax_local.axvspan(mean_local - std_local, mean_local + std_local,
                           alpha=0.2, color='red')
            
            ax_local.set_xlabel('Blaze Angle (°)', fontsize=9)
            ax_local.set_ylabel('Count', fontsize=9)
            ax_local.set_title(f'Within-Facet: {mean_local:.2f}° ± {std_local:.2f}° σ (N={len(local_angles)})',
                              fontsize=9)
            ax_local.grid(True, alpha=0.3)
        else:
            ax_local.text(0.5, 0.5, 'No local angle data', 
                         ha='center', va='center', transform=ax_local.transAxes)
            ax_local.set_xlabel('Blaze Angle (°)', fontsize=9)
    
    plt.tight_layout()
