# AFM Analysis Code Restructuring Guide

## New Package Structure

Your code is now organized as a proper Python package:

```
afm_analysis/                    # Project root
├── main.py                      # Entry point (40 lines!)
├── README.md                    # Documentation
├── requirements.txt             # Dependencies
├── results/                     # Output directory (auto-created)
└── afm_analysis/               # Main package
    ├── __init__.py             # Package initialization
    ├── config.py               # Configuration (unchanged)
    ├── analyzer.py             # Core single-file analyzer (~230 lines)
    ├── workflows.py            # Workflow orchestration (~140 lines)
    │
    ├── core/                   # Core analysis algorithms
    │   ├── __init__.py
    │   ├── processing.py       # Data loading & processing
    │   └── analysis.py         # Blaze angle extraction
    │
    ├── io/                     # Input/output operations
    │   ├── __init__.py
    │   └── file_io.py          # Save results to files
    │
    ├── data/                   # Data management
    │   ├── __init__.py
    │   └── aggregation.py      # Combine & group scans
    │
    ├── stats/                  # Statistical analysis
    │   ├── __init__.py
    │   └── analysis.py         # Comparisons & statistics
    │
    └── visualization/          # All plotting functions
        ├── __init__.py
        ├── diagnostics.py      # Diagnostic plots
        ├── statistics.py       # Statistical plots
        └── profiles.py         # AFM profile plots
```

## File Size Breakdown

| File | Lines | Purpose |
|------|-------|---------|
| **main.py** | **40** | Simple entry point |
| analyzer.py | 230 | Single file analysis logic |
| workflows.py | 140 | Orchestrate the 3 modes |
| core/processing.py | ~200 | Data processing |
| core/analysis.py | ~150 | Angle extraction |
| io/file_io.py | ~200 | File saving |
| data/aggregation.py | ~90 | Data combining |
| stats/analysis.py | ~120 | Statistics |
| visualization/diagnostics.py | ~90 | Diagnostic plots |
| visualization/statistics.py | ~170 | Statistical plots |
| visualization/profiles.py | ~180 | Profile plots |

**Result**: Longest file is now 230 lines (was 470)!

## Migration Steps

### 1. Create the new folder structure

```bash
mkdir -p afm_analysis/core
mkdir -p afm_analysis/io
mkdir -p afm_analysis/data
mkdir -p afm_analysis/stats
mkdir -p afm_analysis/visualization
```

### 2. Move and rename files

**Old files → New locations:**

```bash
# Create __init__.py files
touch afm_analysis/__init__.py
touch afm_analysis/core/__init__.py
touch afm_analysis/io/__init__.py
touch afm_analysis/data/__init__.py
touch afm_analysis/stats/__init__.py
touch afm_analysis/visualization/__init__.py

# Move existing files
mv config.py afm_analysis/config.py
mv processing.py afm_analysis/core/processing.py
mv analysis.py afm_analysis/core/analysis.py
mv file_io.py afm_analysis/io/file_io.py
mv data_aggregation.py afm_analysis/data/aggregation.py
mv statistical_analysis.py afm_analysis/stats/analysis.py
mv plot_diagnostics.py afm_analysis/visualization/diagnostics.py
mv plot_statistics.py afm_analysis/visualization/statistics.py
mv plot_profiles.py afm_analysis/visualization/profiles.py

# Add new files (from artifacts)
# - afm_analysis/__init__.py
# - afm_analysis/analyzer.py
# - afm_analysis/workflows.py
# - main.py (replace existing)
```

### 3. Update imports in moved files

When you move files into subdirectories, update their imports:

**Example for `afm_analysis/core/processing.py`:**
```python
# Old
from config import SOME_SETTING

# New
from ..config import SOME_SETTING
```

**Example for `afm_analysis/visualization/diagnostics.py`:**
```python
# Old
import numpy as np

# New (unchanged for external imports)
import numpy as np
```

### 4. Create simple `__init__.py` files

**`afm_analysis/core/__init__.py`:**
```python
from .processing import load_afm_data, raw_data, flatten_profile, find_groove_positions
from .analysis import extract_blaze_angle

__all__ = ['load_afm_data', 'raw_data', 'flatten_profile', 
           'find_groove_positions', 'extract_blaze_angle']
```

**`afm_analysis/io/__init__.py`:**
```python
from .file_io import save_results_to_file

__all__ = ['save_results_to_file']
```

**`afm_analysis/data/__init__.py`:**
```python
from .aggregation import (group_by_temperature, extract_temperatures_for_output,
                         combine_scans)

__all__ = ['group_by_temperature', 'extract_temperatures_for_output', 'combine_scans']
```

**`afm_analysis/stats/__init__.py`:**
```python
from .analysis import (print_comparison_summary, print_pairwise_comparisons,
                      print_temperature_analysis)

__all__ = ['print_comparison_summary', 'print_pairwise_comparisons',
          'print_temperature_analysis']
```

**`afm_analysis/visualization/__init__.py`:**
```python
from .diagnostics import plot_analyzed_regions_overlay, plot_flattening_diagnostic
from .statistics import plot_summary_statistics, plot_multi_file_comparison
from .profiles import plot_sample_profiles_by_temperature

__all__ = ['plot_analyzed_regions_overlay', 'plot_flattening_diagnostic',
          'plot_summary_statistics', 'plot_multi_file_comparison',
          'plot_sample_profiles_by_temperature']
```

## Usage Examples

### Running analysis (unchanged!)

```bash
python main.py
```

### Advanced usage in Python

```python
from afm_analysis import analyze_single_file, run_comparison_analysis

# Analyze a single file programmatically
result = analyze_single_file('data/my_sample.txt', show_plots=False)
print(f"Mean angle: {result['mean_angle']:.2f}°")

# Or run full comparison workflow
results = run_comparison_analysis()
```

## Benefits of New Structure

### 1. **Clarity** 
- Each file has a single, clear purpose
- Related code grouped in subpackages
- Easy to find what you need

### 2. **Maintainability**
- No file over 230 lines
- Changes isolated to specific modules
- Easy to add new features

### 3. **Testability**
- Can import and test individual components
- Clear dependencies between modules

### 4. **Scalability**
- Easy to add new analysis modes
- Simple to add new plot types
- Clean extension points

### 5. **Professional**
- Standard Python package structure
- Can be pip-installed if needed
- Easy for others to use

## Key Improvements

### main.py: 470 lines → 40 lines
```python
# Before: 470 lines of mixed logic
# After: Simple, clear entry point
def main():
    if ANALYSIS_MODE == 'single':
        run_single_file_analysis()
    elif ANALYSIS_MODE == 'multiple':
        run_multiple_file_analysis()
    elif ANALYSIS_MODE == 'compare':
        run_comparison_analysis()
```

### analyzer.py: Extract single-file logic
- All single-file analysis in one place
- Broken into small, focused functions
- Each function < 50 lines

### workflows.py: Orchestrate analysis modes
- One function per mode
- Clean separation of concerns
- Easy to modify behavior

## Next Steps

### Optional enhancements:

1. **Add tests/**
   ```
   tests/
   ├── test_processing.py
   ├── test_analysis.py
   └── test_aggregation.py
   ```

2. **Add setup.py for installation**
   ```python
   from setuptools import setup, find_packages
   
   setup(
       name='afm_analysis',
       version='1.0.0',
       packages=find_packages(),
       install_requires=['numpy', 'matplotlib', 'scipy'],
   )
   ```

3. **Add examples/**
   ```
   examples/
   ├── basic_analysis.py
   ├── batch_processing.py
   └── custom_workflow.py
   ```

## Checklist

- [ ] Create folder structure
- [ ] Move files to new locations
- [ ] Create all `__init__.py` files
- [ ] Update imports (use relative imports like `from ..config import X`)
- [ ] Test that `python main.py` still works
- [ ] Update your documentation
- [ ] Celebrate having clean, organized code! 🎉