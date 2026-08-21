"""
Blaze angle extraction and analysis functions
"""
import numpy as np
import matplotlib.pyplot as plt
from .processing import extract_single_groove


def extract_blaze_angle(x, y, groove_center, period_nm, trim_fraction=0.2, 
                       side='negative_slope', show_plot=False, groove_num=None,
                       return_local_angles=False, return_regions=False,
                       local_period_nm=None):
    """
    Extract blaze angle from a single groove.
    
    Parameters:
        side: 'negative_slope', 'positive_slope', or 'longer'
        return_local_angles: if True, also return array of local angles
        return_regions: if True, return dict with x,y coordinates of analyzed regions
        local_period_nm: if provided, use this local period instead of global period_nm
    
    Returns:
        blaze_angle, steep_angle, blaze_slope, quality_metrics
    """
    
    # Use local period if provided, otherwise fall back to global estimate
    period_to_use = local_period_nm if local_period_nm is not None else period_nm
    
    # Extract groove region (allow asymmetric for edge grooves)
    x_groove, y_groove = extract_single_groove(x, y, groove_center, period_to_use, 
                                               margin=0.2, allow_asymmetric=True)
    
    # Convert to nm
    x_nm = x_groove * 1000
    y_nm = y_groove
    
    # Find the groove minimum (trough)
    imin = np.argmin(y_nm)
    
    # Split into left and right facets
    left_x = x_nm[:imin+1]
    left_y = y_nm[:imin+1]
    right_x = x_nm[imin:]
    right_y = y_nm[imin:]
    
    # Calculate slopes
    if len(left_x) > 3:
        left_slope = np.polyfit(left_x, left_y, 1)[0]
    else:
        left_slope = 0
    
    if len(right_x) > 3:
        right_slope = np.polyfit(right_x, right_y, 1)[0]
    else:
        right_slope = 0
    
    # Calculate widths
    left_width = left_x[-1] - left_x[0] if len(left_x) > 1 else 0
    right_width = right_x[-1] - right_x[0] if len(right_x) > 1 else 0
    
    # Determine which side to analyze
    if side == 'negative_slope':
        if left_slope < 0 and abs(left_slope) > abs(right_slope):
            blaze_x, blaze_y = left_x, left_y
            steep_x, steep_y = right_x, right_y
            blaze_side_actual = f'left (slope={left_slope:.4f})'
            blaze_width = left_width
            steep_width = right_width
        elif right_slope < 0:
            blaze_x, blaze_y = right_x, right_y
            steep_x, steep_y = left_x, left_y
            blaze_side_actual = f'right (slope={right_slope:.4f})'
            blaze_width = right_width
            steep_width = left_width
        else:
            if left_slope < right_slope:
                blaze_x, blaze_y = left_x, left_y
                steep_x, steep_y = right_x, right_y
                blaze_side_actual = f'left (slope={left_slope:.4f})'
                blaze_width = left_width
                steep_width = right_width
            else:
                blaze_x, blaze_y = right_x, right_y
                steep_x, steep_y = left_x, left_y
                blaze_side_actual = f'right (slope={right_slope:.4f})'
                blaze_width = right_width
                steep_width = left_width
                
    elif side == 'positive_slope':
        if left_slope > 0 and left_slope > right_slope:
            blaze_x, blaze_y = left_x, left_y
            steep_x, steep_y = right_x, right_y
            blaze_side_actual = f'left (slope={left_slope:.4f})'
            blaze_width = left_width
            steep_width = right_width
        elif right_slope > 0:
            blaze_x, blaze_y = right_x, right_y
            steep_x, steep_y = left_x, left_y
            blaze_side_actual = f'right (slope={right_slope:.4f})'
            blaze_width = right_width
            steep_width = left_width
        else:
            if left_slope > right_slope:
                blaze_x, blaze_y = left_x, left_y
                steep_x, steep_y = right_x, right_y
                blaze_side_actual = f'left (slope={left_slope:.4f})'
                blaze_width = left_width
                steep_width = right_width
            else:
                blaze_x, blaze_y = right_x, right_y
                steep_x, steep_y = left_x, left_y
                blaze_side_actual = f'right (slope={right_slope:.4f})'
                blaze_width = right_width
                steep_width = left_width
                
    elif side == 'longer':
        if left_width > right_width:
            blaze_x, blaze_y = left_x, left_y
            steep_x, steep_y = right_x, right_y
            blaze_side_actual = f'left (longer, slope={left_slope:.4f})'
            blaze_width = left_width
            steep_width = right_width
        else:
            blaze_x, blaze_y = right_x, right_y
            steep_x, steep_y = left_x, left_y
            blaze_side_actual = f'right (longer, slope={right_slope:.4f})'
            blaze_width = right_width
            steep_width = left_width
    else:
        raise ValueError(f"side must be 'negative_slope', 'positive_slope', or 'longer', got '{side}'")
    
    # Trim edges to avoid rounded corners AND trough region
    def trim_facet(x, y, trim_frac, is_blaze=False):
        if len(x) < 10:
            return x, y
        n_trim = int(trim_frac * len(x))
        if n_trim < 1:
            return x, y
        
        # For blaze facet, trim more aggressively on the trough side
        # to avoid the flattened bottom region
        if is_blaze:
            # Determine which end is near the trough (lower height)
            if y[0] < y[-1]:  # Trough is at start
                # Trim 2.5x more from the start (trough side)
                n_trim_start = int(n_trim * 2.5)
                n_trim_end = n_trim
            else:  # Trough is at end
                # Trim 2.5x more from the end (trough side)
                n_trim_start = n_trim
                n_trim_end = int(n_trim * 2.5)
            
            return x[n_trim_start:-n_trim_end], y[n_trim_start:-n_trim_end]
        else:
            return x[n_trim:-n_trim], y[n_trim:-n_trim]
    
    blaze_x_trim, blaze_y_trim = trim_facet(blaze_x, blaze_y, trim_fraction, is_blaze=True)
    steep_x_trim, steep_y_trim = trim_facet(steep_x, steep_y, trim_fraction, is_blaze=False)
    
    # Fit linear slope to blaze facet with covariance
    if len(blaze_x_trim) < 3:
        return None, None, None, None
    
    # Get fit coefficients AND covariance matrix for uncertainty estimation
    blaze_coeffs, blaze_cov = np.polyfit(blaze_x_trim, blaze_y_trim, 1, cov=True)
    blaze_slope = blaze_coeffs[0]

    # Calculate uncertainty in slope from covariance matrix
    # Diagonal element [0,0] is variance of slope parameter
    blaze_slope_stderr = np.sqrt(blaze_cov[0, 0])
    
    # Calculate angle and propagate uncertainty
    blaze_angle = abs(np.arctan(blaze_slope) * 180 / np.pi)  # Absolute value
    
    # Error propagation: for θ = arctan(m), dθ/dm = 1/(1+m²)
    # Therefore: σ_θ = |σ_m / (1 + m²)| * (180/π)
    blaze_angle_stderr = abs(blaze_slope_stderr / (1 + blaze_slope**2)) * 180 / np.pi
    
    # Sanity check: warn if angle seems too large (might be measuring wrong facet)
    if blaze_angle > 45:
        import warnings
        warnings.warn(
            f"Groove {groove_num}: Blaze angle of {blaze_angle:.1f}° is unusually large! "
            f"You may be measuring the steep facet instead of the blaze facet. "
            f"Check BLAZE_SIDE setting in config.py",
            UserWarning
        )
    
    # Calculate R² for blaze fit
    blaze_fit = np.polyval(blaze_coeffs, blaze_x_trim)
    blaze_residuals = blaze_y_trim - blaze_fit
    blaze_ss_res = np.sum(blaze_residuals**2)
    blaze_ss_tot = np.sum((blaze_y_trim - np.mean(blaze_y_trim))**2)
    blaze_r2 = 1 - (blaze_ss_res / blaze_ss_tot) if blaze_ss_tot > 0 else 0
    
    # Fit steep facet
    steep_angle = None
    steep_r2 = None
    steep_angle_stderr = None
    if len(steep_x_trim) >= 3:
        steep_coeffs, steep_cov = np.polyfit(steep_x_trim, steep_y_trim, 1, cov=True)
        steep_slope = steep_coeffs[0]
        steep_slope_stderr = np.sqrt(steep_cov[0, 0])
        steep_angle = abs(np.arctan(steep_slope) * 180 / np.pi)  # Absolute value
        steep_angle_stderr = abs(steep_slope_stderr / (1 + steep_slope**2)) * 180 / np.pi
        
        steep_fit = np.polyval(steep_coeffs, steep_x_trim)
        steep_residuals = steep_y_trim - steep_fit
        steep_ss_res = np.sum(steep_residuals**2)
        steep_ss_tot = np.sum((steep_y_trim - np.mean(steep_y_trim))**2)
        steep_r2 = 1 - (steep_ss_res / steep_ss_tot) if steep_ss_tot > 0 else 0
    
    # Quality metrics
    quality_metrics = {
        'blaze_r2': blaze_r2,
        'steep_r2': steep_r2,
        'blaze_width_nm': blaze_width,
        'steep_width_nm': steep_width,
        'groove_depth_nm': np.max(y_nm) - np.min(y_nm),
        'local_period_used_nm': period_to_use,  # Store which period was used
        # NEW: Uncertainty estimates from linear fits
        'blaze_angle_stderr': blaze_angle_stderr,  # Uncertainty in blaze angle (degrees)
        'blaze_slope': blaze_slope,  # Store slope for reference
        'blaze_slope_stderr': blaze_slope_stderr,  # Uncertainty in slope
        'steep_angle_stderr': steep_angle_stderr  # Uncertainty in steep angle (degrees), or None
    }
    
    # Store region information if requested
    if return_regions:
        regions = {
            'groove_x_nm': x_nm,
            'groove_y_nm': y_nm,
            'blaze_x_full': blaze_x,
            'blaze_y_full': blaze_y,
            'blaze_x_trim': blaze_x_trim,
            'blaze_y_trim': blaze_y_trim,
            'steep_x_full': steep_x,
            'steep_y_full': steep_y,
            'steep_x_trim': steep_x_trim,
            'steep_y_trim': steep_y_trim,
            'blaze_fit': blaze_fit,
            'blaze_side': blaze_side_actual
        }
        quality_metrics['regions'] = regions
    
    # Calculate local angles
    if return_local_angles and len(blaze_x_trim) > 10:
        window_size = max(5, len(blaze_x_trim) // 10)
        local_slopes = []
        local_positions = []
        
        for i in range(window_size, len(blaze_x_trim) - window_size):
            x_window = blaze_x_trim[i-window_size:i+window_size]
            y_window = blaze_y_trim[i-window_size:i+window_size]
            local_fit = np.polyfit(x_window, y_window, 1)
            local_slopes.append(local_fit[0])
            local_positions.append(blaze_x_trim[i])
        
        local_angles = np.abs(np.arctan(np.array(local_slopes)) * 180 / np.pi)  # Absolute value
        
        quality_metrics['local_angles'] = local_angles
        quality_metrics['local_positions'] = np.array(local_positions)
        quality_metrics['angle_std'] = np.std(local_angles) if len(local_angles) > 0 else 0
        quality_metrics['angle_range'] = (np.max(local_angles) - np.min(local_angles)) if len(local_angles) > 1 else 0
    
    # Individual groove plot if requested
    if show_plot:
        plt.figure(figsize=(14, 4))
        
        ax1 = plt.subplot(1, 3, 1)
        ax1.plot(x_nm, y_nm, 'k-', linewidth=2, alpha=0.5, label='Full groove')
        ax1.plot(blaze_x, blaze_y, 'b-', linewidth=2, label=f'Blaze facet ({blaze_side_actual})')
        ax1.plot(steep_x, steep_y, 'r-', linewidth=2, label='Steep facet')
        ax1.plot(blaze_x_trim, blaze_y_trim, 'bo', markersize=4, label='Fitted region')
        ax1.plot(blaze_x_trim, blaze_fit, 'b--', linewidth=2, 
                label=f'Fit: {blaze_angle:.2f} deg (R2={blaze_r2:.3f})')
        ax1.set_xlabel('Position (nm)')
        ax1.set_ylabel('Height (nm)')
        title = f'Groove {groove_num}: Blaze = {blaze_angle:.2f} deg'
        if steep_angle is not None:
            title += f', Steep = {steep_angle:.2f} deg'
        if local_period_nm is not None:
            title += f'\nLocal period: {local_period_nm:.2f} nm'
        ax1.set_title(title)
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
        
        ax2 = plt.subplot(1, 3, 2)
        ax2.plot(blaze_x_trim, blaze_residuals, 'bo-', markersize=3)
        ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)
        ax2.set_xlabel('Position (nm)')
        ax2.set_ylabel('Residual (nm)')
        ax2.set_title(f'Fit Quality (RMS = {np.sqrt(np.mean(blaze_residuals**2)):.2f} nm)')
        ax2.grid(True, alpha=0.3)
        
        if return_local_angles and 'local_angles' in quality_metrics:
            ax3 = plt.subplot(1, 3, 3)
            local_ang = quality_metrics['local_angles']
            local_pos = quality_metrics['local_positions']
            ax3.plot(local_pos, local_ang, 'go-', markersize=3, linewidth=1.5)
            ax3.axhline(blaze_angle, color='b', linestyle='--', linewidth=2, 
                       label=f'Mean: {blaze_angle:.2f}°')
            ax3.fill_between([local_pos[0], local_pos[-1]], 
                            blaze_angle - quality_metrics['angle_std'], 
                            blaze_angle + quality_metrics['angle_std'],
                            alpha=0.2, color='blue')
            ax3.set_xlabel('Position (nm)')
            ax3.set_ylabel('Local Angle (degrees)')
            ax3.set_title(f'Local Variation (std = {quality_metrics["angle_std"]:.2f}°)')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
        
        plt.tight_layout()
    
    return blaze_angle, steep_angle, blaze_slope, quality_metrics