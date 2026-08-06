"""
Data aggregation functions
Handles combining multiple scans and grouping by temperature
"""
import numpy as np
from collections import defaultdict


def combine_scans(scan_results):
    """
    Combine multiple scan results into a single aggregated result
    
    Parameters:
        scan_results: list of individual scan result dictionaries
        
    Returns:
        combined_result: aggregated result dictionary
    """
    if len(scan_results) == 1:
        # Single scan - return as is but mark it
        result = scan_results[0].copy()
        result['n_scans'] = 1
        return result
    
    # Collect all angles and quality metrics across scans
    all_angles = []
    all_slopes = []
    all_steep = []
    all_quality = []
    all_local_angles = []
    all_groove_periods = []
    
    for scan in scan_results:
        all_angles.extend(scan['all_angles'])
        all_slopes.extend([scan['mean_slope']] * len(scan['all_angles']))
        if scan['mean_steep'] is not None:
            all_steep.append(scan['mean_steep'])
        all_quality.extend(scan['quality'])
        if 'all_local_angles' in scan:
            all_local_angles.extend(scan['all_local_angles'])
        if 'groove_periods' in scan:
            all_groove_periods.extend(scan['groove_periods'])
    
    # Calculate combined statistics
    mean_angle = np.mean(all_angles)
    std_angle = np.std(all_angles)
    mean_slope = np.mean(all_slopes)
    mean_steep = np.mean(all_steep) if len(all_steep) > 0 else None
    
    # Calculate combined period statistics
    if len(all_groove_periods) > 0:
        period_nm = np.mean(all_groove_periods)
        period_std = np.std(all_groove_periods)
    else:
        # Fall back to average of individual scan periods
        period_nm = np.mean([s['period_nm'] for s in scan_results])
        period_std = np.std([s['period_nm'] for s in scan_results])
    
    # Local angle statistics
    local_angle_std = np.std(all_local_angles) if len(all_local_angles) > 0 else 0
    local_angle_range = (np.max(all_local_angles) - np.min(all_local_angles)) if len(all_local_angles) > 1 else 0
    
    # Create combined result
    combined_result = {
        'filename': scan_results[0]['filename'],  # Keep first filename as reference
        'n_scans': len(scan_results),
        'n_grooves': len(all_angles),
        'mean_angle': mean_angle,
        'std_angle': std_angle,
        'min_angle': np.min(all_angles),
        'max_angle': np.max(all_angles),
        'mean_slope': mean_slope,
        'mean_steep': mean_steep,
        'period_nm': period_nm,
        'period_std': period_std,
        'groove_periods': all_groove_periods,
        'all_angles': all_angles,
        'quality': all_quality,
        'local_angle_std': local_angle_std,
        'local_angle_range': local_angle_range,
        'all_local_angles': all_local_angles,
        'individual_scans': scan_results  # Store individual scans for detailed analysis
    }
    
    # Copy label and temperature if present
    if 'label' in scan_results[0]:
        combined_result['label'] = scan_results[0]['label']
    if 'temperature' in scan_results[0]:
        combined_result['temperature'] = scan_results[0]['temperature']
    
    return combined_result


def group_by_temperature(results, labels, temperatures):
    """
    Group results by temperature, combining multiple scans at same temperature
    
    Parameters:
        results: list of result dictionaries
        labels: list of sample labels
        temperatures: list of temperatures (None for master)
        
    Returns:
        grouped_results: dict mapping temperature key to combined result
        grouped_labels: dict mapping temperature key to label
        temp_order: list of temperature keys in order (master first, then sorted temps)
    """
    # Group results by temperature
    temp_groups = defaultdict(list)
    temp_labels = {}
    
    for result, label, temp in zip(results, labels, temperatures):
        # Create a key for grouping
        if temp is None:
            key = 'master'
            temp_labels[key] = label
        else:
            key = f'{temp}C'
            # Use a generic label for the temperature if multiple scans
            if key in temp_labels:
                temp_labels[key] = f'{temp}°C'
            else:
                temp_labels[key] = label
        
        temp_groups[key].append(result)
    
    # Combine scans within each temperature group
    grouped_results = {}
    for key, scans in temp_groups.items():
        grouped_results[key] = combine_scans(scans)
    
    # Create ordered list: master first, then temperatures in ascending order
    temp_order = []
    if 'master' in grouped_results:
        temp_order.append('master')
    
    # Sort numeric temperature keys
    numeric_keys = [k for k in grouped_results.keys() if k != 'master']
    numeric_keys.sort(key=lambda x: float(x.replace('C', '')))
    temp_order.extend(numeric_keys)
    
    return grouped_results, temp_labels, temp_order


def extract_temperatures_for_output(grouped_results, temp_order):
    """
    Extract temperature values for CSV output
    
    Parameters:
        grouped_results: dict mapping temperature key to result
        temp_order: list of temperature keys in order
        
    Returns:
        temperatures: list of temperature values (None for master)
    """
    temperatures = []
    for key in temp_order:
        if key == 'master':
            temperatures.append(None)
        else:
            temperatures.append(float(key.replace('C', '')))
    
    return temperatures