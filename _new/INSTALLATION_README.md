# AFM Statistical Analysis - Installation Guide

## Files Overview

You have 6 files to integrate into your project:

### 📚 Documentation (Read First!)
1. **statistical_guide.md** - Explains the concepts (ICC, hierarchical stats, etc.)
2. **statistical_decision_guide.md** - Helps you decide which approach to use
3. **integration_guide.md** - Step-by-step code integration instructions

### 🔧 Code Files (Add to Your Project)
4. **improved_statistics.py** - Core statistical functions (REPLACES your current stats/analysis.py)
5. **diagnostic_plots.py** - Visualization functions (ADD to visualization/)
6. **correlation_source_analysis.py** - Tools to check ICC and determine approach (ADD to stats/)

---

## Quick Start (5 Steps)

### Step 1: Read the Decision Guide (5 min)
```bash
# Read this first to understand which approach you need
cat statistical_decision_guide.md
```

### Step 2: Check Your Data (5 min)

Run this quick test on ONE sample to see if you need hierarchical stats:

```python
# test_my_data.py
import numpy as np
from correlation_source_analysis import quick_icc_check

# Load one of your existing result files
# Replace with your actual data loading:
angles = [17.2, 17.5, 17.3, ...]  # Your measurements
row_groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, ...]  # Which group each belongs to

# Check ICC
quick_icc_check(angles, row_groups, "Test Sample")
```

**If output says:** 
- `ICC < 0.1 → Simple statistics OK ✓` → **You're done!** Use simple stats (see option A below)
- `ICC > 0.1 → Hierarchical stats recommended` → Use full hierarchical approach (option B below)

### Step 3A: If ICC < 0.1 (Simple Statistics)

**GOOD NEWS:** Your grooves are independent! Minimal changes needed.

**Files to add:**
- `correlation_source_analysis.py` → Copy to `stats/` directory (for reporting ICC)

**Your existing code is mostly fine!** Just add ICC reporting:

```python
# In your analyzer, add this:
from stats.correlation_source_analysis import quick_icc_check

# After collecting measurements:
quick_icc_check(all_angles, all_row_groups, label)

# Continue with your existing simple statistics:
mean = np.mean(all_angles)
sem = np.std(all_angles, ddof=1) / np.sqrt(len(all_angles))
```

**In your paper, report:**
> "Measurements from individual grooves were statistically independent (ICC = 0.08), yielding N = 100 independent measurements per sample."

### Step 3B: If ICC > 0.1 (Hierarchical Statistics)

**Need full hierarchical approach.**

**Files to add:**
1. `improved_statistics.py` → **REPLACE** `stats/analysis.py` 
   - Or rename your current file to `analysis_old.py`
   - Then copy this as `analysis.py`

2. `diagnostic_plots.py` → **ADD** to `visualization/` directory
   - Copy as-is, imports will work

3. `correlation_source_analysis.py` → **ADD** to `stats/` directory
   - For diagnostics and ICC checking

**Then follow the integration_guide.md** for detailed code changes.

---

## File Installation Details

### Directory Structure (After Installation)

```
your_afm_project/
├── afm_analysis/
│   ├── core/
│   │   ├── processing.py         (your existing file)
│   │   ├── analysis.py           (your existing file)
│   │   └── __init__.py
│   ├── stats/
│   │   ├── analysis.py           ← REPLACE with improved_statistics.py
│   │   ├── correlation_source_analysis.py  ← NEW (add this)
│   │   └── __init__.py
│   ├── visualization/
│   │   ├── diagnostics.py        (your existing file)
│   │   ├── profiles.py           (your existing file)
│   │   ├── statistics.py         (your existing file)
│   │   ├── diagnostic_plots.py   ← NEW (add this)
│   │   └── __init__.py
│   └── workflows.py
├── config.py
├── main.py
└── README.md
```

### Specific File Actions

#### 1. improved_statistics.py

**Location:** `stats/improved_statistics.py` (or rename to `analysis.py`)

**Purpose:** Core hierarchical statistical functions

**What it replaces:** Your current `stats/analysis.py` comparison functions

**Key functions:**
- `calculate_hierarchical_statistics()` - Main statistics calculator
- `compare_samples_hierarchical()` - Proper t-tests with ICC
- `print_hierarchical_statistics()` - Enhanced output
- `test_normality()` - Distribution checks

**Import as:**
```python
from stats.improved_statistics import calculate_hierarchical_statistics
# or if renamed:
from stats.analysis import calculate_hierarchical_statistics
```

#### 2. diagnostic_plots.py

**Location:** `visualization/diagnostic_plots.py`

**Purpose:** Diagnostic visualization for hierarchical data

**What it adds:** New plot types (doesn't replace anything)

**Key functions:**
- `plot_variance_components()` - Pie chart showing uncertainty sources
- `plot_row_group_consistency()` - Spatial variation across image
- `plot_qq_normality()` - Test normality assumptions
- `plot_measurement_uncertainty_validation()` - Check if fit errors are realistic
- `create_diagnostic_report()` - Generate all plots at once

**Import as:**
```python
from visualization.diagnostic_plots import create_diagnostic_report
```

#### 3. correlation_source_analysis.py

**Location:** `stats/correlation_source_analysis.py`

**Purpose:** Determine if you need hierarchical stats

**What it adds:** Diagnostic tools (doesn't replace anything)

**Key functions:**
- `quick_icc_check()` - Fast check for hierarchical vs simple
- `analyze_correlation_source()` - Full diagnostic analysis
- `recommend_statistical_approach()` - Tells you what to use

**Import as:**
```python
from stats.correlation_source_analysis import quick_icc_check
```

---

## Testing Your Installation

### Test 1: Check Imports

```python
# test_imports.py
try:
    from stats.improved_statistics import calculate_hierarchical_statistics
    print("✓ improved_statistics imported successfully")
except ImportError as e:
    print(f"✗ Error importing improved_statistics: {e}")

try:
    from visualization.diagnostic_plots import create_diagnostic_report
    print("✓ diagnostic_plots imported successfully")
except ImportError as e:
    print(f"✗ Error importing diagnostic_plots: {e}")

try:
    from stats.correlation_source_analysis import quick_icc_check
    print("✓ correlation_source_analysis imported successfully")
except ImportError as e:
    print(f"✗ Error importing correlation_source_analysis: {e}")

print("\n✓ All imports successful!")
```

### Test 2: Run Example

```python
# test_statistics.py
import numpy as np
from stats.improved_statistics import calculate_hierarchical_statistics, print_hierarchical_statistics

# Simulate data like yours
np.random.seed(42)
angles = np.random.normal(17.5, 0.4, 100)
row_groups = np.repeat(range(20), 5)
meas_errors = np.random.uniform(0.05, 0.15, 100)

# Calculate statistics
stats = calculate_hierarchical_statistics(angles, row_groups, meas_errors)

# Print results
print_hierarchical_statistics(stats, "Test Sample")

# Verify
assert 0 <= stats['intraclass_correlation'] <= 1, "ICC out of range!"
assert stats['sem_conservative'] >= stats['sem_best'] >= stats['sem_liberal'], "SEM ordering wrong!"
print("\n✓ All tests passed!")
```

### Test 3: Check ICC on Your Data

```python
# test_my_icc.py
from stats.correlation_source_analysis import quick_icc_check

# Load YOUR actual data here
# angles = load_your_data(...)
# row_groups = ...

quick_icc_check(angles, row_groups, "My Sample")

# This tells you whether to use simple or hierarchical stats
```

---

## Integration Roadmap

### Phase 1: ICC Check (Day 1)
- [ ] Copy `correlation_source_analysis.py` to `stats/`
- [ ] Run `quick_icc_check()` on all your samples
- [ ] Determine: Simple or Hierarchical?

### Phase 2: Statistics Update (Day 1-2)
**If ICC < 0.1 (Simple):**
- [ ] Add ICC reporting to your existing code
- [ ] Done! ✓

**If ICC > 0.1 (Hierarchical):**
- [ ] Copy `improved_statistics.py` to `stats/`
- [ ] Copy `diagnostic_plots.py` to `visualization/`
- [ ] Follow `integration_guide.md` step-by-step
- [ ] Update analyzer to track row_groups and meas_errors
- [ ] Update comparison functions

### Phase 3: Testing (Day 2)
- [ ] Run analysis on one sample
- [ ] Check that ICC values make sense
- [ ] Verify SEM values are reasonable
- [ ] Compare old vs new results

### Phase 4: Diagnostics (Day 3)
- [ ] Generate diagnostic plots for all samples
- [ ] Review variance decomposition
- [ ] Check normality assumptions
- [ ] Validate measurement uncertainties

### Phase 5: Production (Day 3+)
- [ ] Process all samples with new statistics
- [ ] Generate final plots
- [ ] Update manuscript/report
- [ ] Archive old results

---

## Minimal Changes Example

If you just want to **check ICC without changing anything else**:

```python
# Add this ONE function to your existing code:

def check_and_report_icc(all_angles, all_row_groups, label):
    """Quick ICC check - add to your analyzer"""
    from stats.correlation_source_analysis import quick_icc_check
    
    need_hierarchical = quick_icc_check(all_angles, all_row_groups, label)
    
    if not need_hierarchical:
        print(f"  → Using simple statistics (all measurements independent)")
        return 'simple'
    else:
        print(f"  → Should use hierarchical statistics")
        return 'hierarchical'

# Then in your workflow:
approach = check_and_report_icc(result['angles'], result['row_groups'], label)

if approach == 'simple':
    # Your existing simple statistics are fine!
    sem = result['std_angle'] / np.sqrt(result['n_grooves'])
else:
    # Need to implement hierarchical
    print("  WARNING: Hierarchical stats recommended but not implemented")
    print("  Current SEM may be too optimistic")
```

---

## Support & Troubleshooting

### Common Issues

**Import errors:**
- Check file paths match your directory structure
- Make sure `__init__.py` files exist in all directories
- Update imports if you renamed files

**ICC values seem wrong:**
- Check that row_groups are correctly assigned (0, 1, 2, ... not scattered)
- Verify each group has multiple measurements
- Make sure angles and row_groups arrays have same length

**SEM values very different:**
- This is expected if ICC > 0.2
- Review diagnostic plots to understand why
- Check if spatial gradients are real vs artifacts

**Plots don't generate:**
- Check matplotlib backend is working
- Try `plt.show()` explicitly
- Verify all data fields are present

### Questions?

If you have questions about:
- **Concepts:** Read `statistical_guide.md`
- **Which approach:** Read `statistical_decision_guide.md`  
- **Implementation:** Read `integration_guide.md`
- **Specific functions:** Check docstrings in the code

---

## Summary

**Minimum to get started:**
1. Copy `correlation_source_analysis.py` to `stats/`
2. Run `quick_icc_check()` on your data
3. If ICC < 0.1: You're done! Keep using simple stats
4. If ICC > 0.1: Follow integration guide for full hierarchical approach

**Files you MUST copy:**
- `correlation_source_analysis.py` (for ICC checking)

**Files you MIGHT need:**
- `improved_statistics.py` (only if ICC > 0.1)
- `diagnostic_plots.py` (only if you want diagnostic plots)

**Files that are just documentation:**
- `statistical_guide.md` (read for understanding)
- `statistical_decision_guide.md` (read to decide approach)
- `integration_guide.md` (read if doing full integration)

Good luck! Start with the ICC check and go from there.
