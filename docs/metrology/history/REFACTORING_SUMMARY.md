# AFM Analysis Code Refactoring Summary

## New File Structure

The code has been reorganized into focused, maintainable modules:

### Core Analysis Files (unchanged)
- **`config.py`** - Configuration parameters
- **`processing.py`** - Data loading and processing
- **`analysis.py`** - Blaze angle extraction

### New Modular Structure

#### Data Management
- **`data_aggregation.py`** (~90 lines)
  - `combine_scans()` - Combines multiple AFM scans into single result
  - `group_by_temperature()` - Groups results by temperature
  - `extract_temperatures_for_output()` - Helper for export

- **`file_io.py`** (~200 lines)
  - `save_results_to_file()` - Saves text summaries, CSV, and per-groove data
  - Handles both single and combined scan results

#### Statistical Analysis
- **`statistical_analysis.py`** (~120 lines)
  - `print_comparison_summary()` - Summary tables
  - `print_pairwise_comparisons()` - Statistical pairwise tests
  - `print_temperature_analysis()` - Temperature-dependent analysis

#### Visualization (split from 360+ lines)
- **`plot_diagnostics.py`** (~90 lines)
  - `plot_analyzed_regions_overlay()` - Shows analyzed regions
  - `plot_flattening_diagnostic()` - Before/after flattening

- **`plot_statistics.py`** (~170 lines)
  - `plot_summary_statistics()` - Histograms and distributions
  - `plot_multi_file_comparison()` - Bar charts and multi-sample histograms

- **`plot_profiles.py`** (~180 lines)
  - `plot_sample_profiles_by_temperature()` - Temperature-organized profiles
  - `_plot_multiple_scans()` - Multiple scans at same temperature
  - `_plot_single_scan()` - Single scan plotting
  - `plot_sample_profiles_comparison()` - Legacy function (deprecated)

#### Main Script
- **`main.py`** (~320 lines, down from ~470)
  - `analyze_file()` - Single file analysis
  - Main execution logic for three modes: single, multiple, compare

## Key Improvements

### 1. Better Organization
- Each file has a clear, single responsibility
- Related functions are grouped together
- Easier to find and modify specific functionality

### 2. Multiple Scans Support
- Scans at the same temperature are automatically combined
- Statistics calculated across all grooves from all scans
- Individual scan data preserved for detailed visualization

### 3. Temperature-Organized Plotting
- AFM profiles organized by temperature (one figure per temperature)
- Multiple scans at same temperature shown as subplots
- Cleaner, more logical organization

### 4. Improved Maintainability
- Smaller, focused files (~90-200 lines each vs 360+ lines)
- Clear function names and responsibilities
- Better separation of concerns

## Migration Guide

### No changes needed to:
- `config.py`
- `processing.py`
- `analysis.py`
- Your analysis workflow

### New files to add:
1. `data_aggregation.py`
2. `file_io.py`
3. `statistical_analysis.py`
4. `plot_diagnostics.py`
5. `plot_statistics.py`
6. `plot_profiles.py`

### Files to replace:
1. `main.py` - Updated with new imports
2. `visualization.py` - Can be deleted (replaced by 3 plot_*.py files)

## Usage Examples

### Multiple Scans at Same Temperature
```python
SAMPLES_TO_COMPARE = [
    ('data/master.txt', 'Master', None),
    ('data/150C_scan1.txt', '150°C', 150),
    ('data/150C_scan2.txt', '150°C', 150),  # Automatically combined
    ('data/175C.txt', '175°C', 175),
]
```

The two 150°C scans will be:
- Combined into single statistics
- Shown as separate subplots on the 150°C profile figure
- Tracked in CSV with `N_scans=2`

### Profile Organization
- **Master grating**: One figure
- **150°C**: One figure with 2 subplots (if 2 scans)
- **175°C**: One figure
- Much easier to navigate than one giant multi-panel figure!

## File Size Comparison

| Old Structure | Lines | New Structure | Lines |
|---------------|-------|---------------|-------|
| visualization.py | ~360 | plot_diagnostics.py | ~90 |
| | | plot_statistics.py | ~170 |
| | | plot_profiles.py | ~180 |
| main.py | ~470 | main.py | ~320 |
| | | data_aggregation.py | ~90 |
| | | file_io.py | ~200 |
| | | statistical_analysis.py | ~120 |

**Result**: Longest file went from 470 lines → 320 lines, visualization split into 3 focused modules