"""
Core data processing functions for AFM analysis
"""
import numpy as np
from scipy.signal import find_peaks


def load_afm_data(filename, default_scan_size=2.0):
    """
    Load AFM data and try to extract scan size from header
    
    Returns:
        data: 2D array of height data
        scan_x_size: scan width in microns
    """
    import re
    
    scan_x_size = default_scan_size
    
    # Try to read scan size from header
    try:
        with open(filename, 'r') as f:
            # Read first 10 lines to look for scan size info
            header_lines = [f.readline() for _ in range(10)]
            
            for line in header_lines:
                # Common AFM file formats include scan size info
                # Look for patterns like "Width: 2.000 µm"
                if 'width' in line.lower():
                    match = re.search(r'(\d+\.?\d*)\s*(um|µm|micron)', line, re.IGNORECASE)
                    if match:
                        scan_x_size = float(match.group(1))
                        print(f"  Detected scan width from header: {scan_x_size} µm")
                        break
                    
                    match = re.search(r'(\d+\.?\d*)\s*nm', line, re.IGNORECASE)
                    if match:
                        scan_x_size = float(match.group(1)) / 1000  # Convert nm to µm
                        print(f"  Detected scan width from header: {scan_x_size} µm ({match.group(1)} nm)")
                        break
                        
                elif 'scan' in line.lower() and 'size' in line.lower():
                    # Extract number
                    match = re.search(r'(\d+\.?\d*)\s*(um|µm|micron)', line, re.IGNORECASE)
                    if match:
                        scan_x_size = float(match.group(1))
                        print(f"  Detected scan size from header: {scan_x_size} µm")
                        break
    except:
        pass
    
    # Load the data
    try:
        data = np.genfromtxt(filename, skip_header=4, skip_footer=0)
    except Exception as e:
        raise Exception(f"Error loading {filename}: {e}")
    
    if scan_x_size == default_scan_size:
        print(f"  Using default scan size: {scan_x_size} µm")
    
    return data, scan_x_size


def raw_data(data, scan_x_size):
    """Extract displacement and height profile from 2D AFM data (averages all rows)"""
    disp_um = scan_x_size * np.arange(0, len(data[0])) / (len(data[0]) - 1)
    groov_nm = (10**9) * (np.mean(data, axis=0) - np.min(np.mean(data, axis=0)))
    return disp_um, groov_nm


def raw_data_multi_group(data, scan_x_size, n_groups=10):
    """
    Extract multiple height profiles from different row groups in 2D AFM data
    
    This function divides the AFM image into n_groups horizontal bands and extracts
    an averaged profile from each band. This allows for many more measurements per
    image and enables assessment of within-image variation.
    
    Parameters:
        data: 2D array of AFM height data (rows × columns)
        scan_x_size: scan width in microns
        n_groups: number of row groups to create (default: 10)
        
    Returns:
        disp_um: 1D array of x-positions in microns (same for all groups)
        profiles_nm: list of n_groups 1D arrays, each containing a height profile in nm
        group_info: dict with metadata about the grouping
    """
    n_rows = data.shape[0]
    n_cols = data.shape[1]
    
    # Calculate x-displacement (same for all profiles)
    disp_um = scan_x_size * np.arange(0, n_cols) / (n_cols - 1)
    
    # Determine row ranges for each group
    rows_per_group = n_rows // n_groups
    
    if rows_per_group < 3:
        print(f"Warning: Only {rows_per_group} rows per group with n_groups={n_groups}")
        print(f"  Reducing to n_groups={n_rows // 3} to ensure at least 3 rows per group")
        n_groups = max(1, n_rows // 3)
        rows_per_group = n_rows // n_groups
    
    profiles_nm = []
    group_ranges = []
    
    for i in range(n_groups):
        # Calculate row range for this group
        start_row = i * rows_per_group
        
        # For the last group, include any remaining rows
        if i == n_groups - 1:
            end_row = n_rows
        else:
            end_row = (i + 1) * rows_per_group
        
        # Extract and average rows in this group
        group_data = data[start_row:end_row, :]
        profile_avg = np.mean(group_data, axis=0)
        
        # Convert to nm and subtract minimum
        profile_nm = (10**9) * (profile_avg - np.min(profile_avg))
        
        profiles_nm.append(profile_nm)
        group_ranges.append((start_row, end_row))
    
    # Package metadata
    group_info = {
        'n_groups': n_groups,
        'n_rows': n_rows,
        'rows_per_group': rows_per_group,
        'group_ranges': group_ranges
    }
    
    print(f"  Extracted {n_groups} row-group profiles")
    print(f"    Rows per group: {rows_per_group}")
    print(f"    Total rows used: {n_groups * rows_per_group} of {n_rows}")
    
    return disp_um, profiles_nm, group_info


def flatten_profile(x, y, method='linear', poly_order=2, exclude_edges=0.0, 
                   feature='both', period_nm=None):
    """
    Remove background tilt/curvature from AFM profile
    
    Parameters:
        method: 'linear', 'polynomial', 'groove_peaks', or 'level_grooves'
        poly_order: polynomial order (for 'polynomial' and 'level_grooves' methods)
        exclude_edges: fraction of data to exclude from edges (0.0 to 0.5)
        feature: for 'level_grooves': 'peaks', 'troughs', or 'both'
        period_nm: estimated period for finding grooves (for 'level_grooves')
    
    Returns:
        y_flat: flattened profile
        background: the removed background
    """
    
    if method == 'linear':
        poly_coeffs = np.polyfit(x, y, 1)
        background = np.polyval(poly_coeffs, x)
        
    elif method == 'polynomial':
        if exclude_edges > 0:
            n_exclude = int(len(x) * exclude_edges)
            x_fit = x[n_exclude:-n_exclude]
            y_fit = y[n_exclude:-n_exclude]
            poly_coeffs = np.polyfit(x_fit, y_fit, poly_order)
        else:
            poly_coeffs = np.polyfit(x, y, poly_order)
        
        background = np.polyval(poly_coeffs, x)
        
    elif method == 'groove_peaks':
        peaks, _ = find_peaks(y, distance=len(y)//20)
        
        if len(peaks) < 3:
            print("Warning: Not enough peaks found for groove_peaks method, using polynomial")
            poly_coeffs = np.polyfit(x, y, poly_order)
            background = np.polyval(poly_coeffs, x)
        else:
            poly_coeffs = np.polyfit(x[peaks], y[peaks], poly_order)
            background = np.polyval(poly_coeffs, x)
            
    elif method == 'level_grooves':
        if period_nm is None:
            raise ValueError("period_nm must be provided for level_grooves method")
        
        feature_x = []
        feature_y = []
        
        if feature in ['peaks', 'both']:
            dx_nm = (x[1] - x[0]) * 1000
            min_distance = int(0.5 * period_nm / dx_nm)  # Changed from 0.3 to 0.5 for more reliable detection
            # Use a lower prominence to catch more peaks
            height_range = np.max(y) - np.min(y)
            min_prominence = 0.2 * height_range  # Require peaks to be at least 20% of height range
            peaks, _ = find_peaks(y, distance=min_distance, prominence=min_prominence)
            if len(peaks) > 0:
                feature_x.extend(x[peaks])
                feature_y.extend(y[peaks])
        
        if feature in ['troughs', 'both']:
            y_inv = -y
            dx_nm = (x[1] - x[0]) * 1000
            min_distance = int(0.5 * period_nm / dx_nm)  # Changed from 0.3 to 0.5
            height_range = np.max(y) - np.min(y)
            min_prominence = 0.2 * height_range
            troughs, _ = find_peaks(y_inv, distance=min_distance, prominence=min_prominence)
            if len(troughs) > 0:
                feature_x.extend(x[troughs])
                feature_y.extend(y[troughs])
        
        if len(feature_x) < 2:
            print("    Warning: Not enough features found for level_grooves, using linear")
            poly_coeffs = np.polyfit(x, y, 1)
            background = np.polyval(poly_coeffs, x)
        else:
            feature_x = np.array(feature_x)
            feature_y = np.array(feature_y)
            
            if exclude_edges > 0:
                x_min = x[0] + exclude_edges * (x[-1] - x[0])
                x_max = x[-1] - exclude_edges * (x[-1] - x[0])
                mask = (feature_x >= x_min) & (feature_x <= x_max)
                feature_x = feature_x[mask]
                feature_y = feature_y[mask]
            
            if len(feature_x) < 2:
                print("    Warning: Not enough features after edge exclusion, using all")
            
            poly_coeffs = np.polyfit(feature_x, feature_y, poly_order)
            background = np.polyval(poly_coeffs, x)
    else:
        raise ValueError(f"Unknown flattening method: {method}")
    
    y_flat = y - background
    
    return y_flat, background


def find_groove_positions(x, y, period_nm, prominence_factor=0.3, distance_factor=0.7):
    """Find positions of groove centers (minima)"""
    y_inv = -y  # Invert so grooves become peaks
    dx_nm = (x[1] - x[0]) * 1000
    period_idx = int(period_nm / dx_nm)
    height_range = np.max(y) - np.min(y)
    min_prominence = prominence_factor * height_range
    
    peaks, properties = find_peaks(y_inv, 
                                   distance=int(distance_factor * period_idx), 
                                   prominence=min_prominence)
    
    return peaks


def extract_single_groove(x, y, center_idx, period_nm, margin=0.1, allow_asymmetric=True):
    """
    Extract a single groove around center_idx
    
    Parameters:
        allow_asymmetric: If True, allows asymmetric extraction near edges
    """
    dx_nm = (x[1] - x[0]) * 1000
    half_width = int(round((period_nm / 2) * (1 + margin) / dx_nm))
    
    if allow_asymmetric:
        # Allow asymmetric extraction for edge grooves
        left_extent = min(half_width, center_idx)
        right_extent = min(half_width, len(x) - 1 - center_idx)
        
        start = center_idx - left_extent
        end = center_idx + right_extent + 1
    else:
        # Original symmetric extraction
        L = min(half_width, center_idx, len(x) - 1 - center_idx)
        start = center_idx - L
        end = center_idx + L + 1
    
    x_groove = x[start:end]
    y_groove = y[start:end]
    
    return x_groove, y_groove
