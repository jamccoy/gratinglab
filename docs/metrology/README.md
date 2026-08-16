# AFM Blaze Angle Analysis (in progress)

A Python toolkit for analyzing Atomic Force Microscopy (AFM) data of diffraction
gratings. One pipeline, two outputs:

- **Blaze angles** — facet angles measured per groove, with row-group analysis for
  many measurements per image and statistical comparison across samples
- **PCGrate boundary profiles** — the averaged groove shape, normalised to one
  period and exported as `.ggp` for grating-efficiency modelling

Both share the same front-end (loading, flattening, groove detection, groove
windowing) and diverge only at the final step. That sharing is deliberate: these
were once separate scripts with their own copies of that code, and they drifted
apart until a bug fixed in one silently persisted in the other.

## Features

- 📊 **High-precision blaze angle extraction** from AFM topography data
- 🔬 **Row-group analysis**: Extract 50-100+ measurements per image (vs. 4-10 with traditional averaging)
- 🌡️ **Temperature studies**: Automatic scan aggregation and temperature-dependent analysis
- 📈 **Comprehensive statistics**: Proper uncertainty quantification with SEM, confidence intervals, and significance tests
- 📉 **Facet curvature analysis**: Quantify within-facet angle variation (camber) across each blaze facet
- 🎯 **Quality control**: Automatic filtering to avoid measuring groove tops and ensure accurate facet detection
- 💾 **Rich outputs**: Detailed reports, CSV data, publication-ready visualizations

## Key Capabilities

### 1. Row-Group Analysis for More Measurements
Traditional AFM analysis averages all rows into a single profile, extracting only 4-10 measurements per image. This package uses **row-group analysis** to divide the image into multiple horizontal bands, giving you:

- **50-100+ measurements** per image instead of 4-10
- **Spatial variation mapping** across the sample
- **Smaller error bars** — but not by √N, see below

**The bands are not independent**, and the statistics account for that. Each band
re-measures the same physical grooves, so SEM does *not* scale as 1/√N. Measured
intraclass correlation is 0.097–0.429 (median 0.244), which puts the effective
sample size at 37–70 rather than 78–102.

Every standard error, confidence interval and p-value is computed on the
**effective** sample size, `N_eff = N / (1 + (m−1)·ICC)`. Reported SEMs are
1.17–1.65× wider than the naive value as a result. `ANALYSIS_MODE = 'icc'`
reproduces the diagnostic; `ICC`, `N_eff`, `SEM_deg` and `SEM_corrected_deg` are
columns in every `analysis_data_*.csv`.

Applying it changed no conclusion on the bundled samples: master-vs-treated and
the 500°C step stay highly significant, and the adjacent-temperature steps stay
null — slightly more clearly so. Effect sizes (Cohen's d) deliberately keep the
raw counts, since a standardised mean difference describes spread rather than
inferring from it.

### 2. Facet Curvature (Camber) Quantification
Beyond measuring the overall blaze angle of each groove, the software quantifies **within-facet angle variation** by analyzing local slopes along each facet. This reveals:

- Facet curvature/camber (how much the angle varies along the facet)
- Manufacturing quality (consistent vs. variable facets)
- Thermal effects on facet shape
- Shown in dedicated histograms and statistics

### 3. Robust Groove Detection with Top Avoidance
Ensures accurate measurements by:

- Automatically detecting groove centers via prominence-based peak finding
- Rejecting grooves clipped by the scan edge, whose facet is only partly present
- Trimming facet edges to avoid measuring groove tops/bottoms
- Recording a fit R² for every groove, in the console output and the per-groove CSV
- Visual diagnostics to verify proper flattening and groove detection

Note that R² is **reported, not enforced** — no measurement is currently dropped for
a poor fit. In practice the edge-exclusion rule removes them anyway: after it, no fit
in the sample dataset falls below R² 0.95.

### 4. Proper Statistical Analysis
All comparisons include:

- Standard Error of Mean (SEM) - uncertainty in mean estimates
- 95% confidence intervals
- Significance tests (p-values) with Welch's t-test
- Effect sizes (Cohen's d) for practical interpretation
- Variance decomposition (measurement vs. physical variation)

---

## Quick Start

### Installation

Requires Python 3.12 (pinned in `.python-version`).

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

That installs the analysis, the PCGrate export and the ICC diagnostic, and pulls
in **no GUI toolkit**. The Qt window is an optional extra:

```bash
.venv/bin/pip install -e '.[dev,gui]'
```

PySide6 (LGPL) rather than PyQt (GPL), so this BSD-3 package stays free of a GPL
dependency. Earlier versions listed PyQt5 as a hard requirement, which put GPL
code in the install path of anyone who only wanted the analysis.

Verified against numpy 2.5.1 / scipy 1.18.0 / matplotlib 3.11.1 — see
`docs/BASELINES.md` if you need to reproduce a stored result exactly.

### Basic Usage

1. **Configure your analysis** in `afm_analysis/config.py`:
   ```python
   # Analysis mode
   ANALYSIS_MODE = 'single'  # or 'multiple' or 'compare'
   SINGLE_FILE = 'data/your_sample.txt'

   # Grating specifications (UPDATE THIS FOR YOUR SAMPLE!)
   PERIOD_EST = 315.0  # Groove spacing in nm (NOT 833 nm!)

   # Row-group analysis (RECOMMENDED)
   USE_ROW_GROUPS = True  # Extract many measurements per image
   N_ROW_GROUPS = 20      # Number of spatial regions to analyze

   # Analysis parameters
   BLAZE_SIDE = 'negative_slope'  # 'negative_slope', 'positive_slope', or 'longer'
   FACET_TRIM = 0.1               # Trim edges to avoid groove tops (0.05-0.28)
   ```

2. **Run the analysis**:
   ```bash
   .venv/bin/python main.py     # or: afm-analysis
   ```

   The Qt window is `afm-gui`. It has two tabs: **Analysis**, and a **Wiki**
   explaining how the statistics work — the correlation correction, why edge
   grooves are rejected, what each output column means.

   Paths in `config.py` resolve against the project root, so this works from any
   working directory.

3. **Review results**:
   - **Console output**: Statistical summary with SEM, p-values, significance
   - **Plots**: Bar charts (with SEM error bars), histograms (per-groove and facet curvature)
   - **Files**: `results/` folder contains detailed reports and CSV data

---

## Configuration Guide

### Critical Parameters (Set These First!)

#### Grating Specifications
```python
# Groove spacing - MUST MATCH YOUR SAMPLE
PERIOD_EST = 315.0   # nm (e.g., ~315 nm for typical echelle gratings)
                     # NOT 833 nm! That's for a different grating type

SCAN_X_SIZE = 2.0    # um - fallback scan width, used only when the width
                     # cannot be read from the data file header
```

#### Row-Group Analysis (Highly Recommended)
```python
USE_ROW_GROUPS = True   # Enable multi-region analysis
N_ROW_GROUPS = 20       # Number of spatial regions (10-50 typical)
                        # More groups = more measurements = smaller SEM
                        # But requires good S/N ratio in data
```

### Data Processing Parameters

#### Flattening (Background Removal)
```python
FLATTEN_METHOD = 'level_grooves'  # 'linear', 'polynomial', 'groove_peaks',
                                  # or 'level_grooves'
FLATTEN_POLY_ORDER = 2            # For 'polynomial' and 'level_grooves':
                                  # 1=linear, 2=quadratic, 3=cubic
FLATTEN_EXCLUDE_EDGES = 0.05      # Exclude edge fraction to avoid boundary effects
FLATTEN_FEATURE = 'peaks'         # For 'level_grooves' only:
                                  # 'peaks' (lands), 'troughs' (groove bottoms), 'both'
```

**Key Point**: The flattening diagnostic plot shows the **averaged profile** to verify proper background removal, even when using row-group analysis. This ensures all subsequent measurements are accurate.

#### Groove Detection
```python
PROMINENCE_FACTOR = 0.01  # Peak detection sensitivity
                          # Lower = detect more grooves (including noise)
                          # Higher = only strong grooves

DISTANCE_FACTOR = 0.3     # Minimum separation between grooves,
                          # as a fraction of the period
                          # Prevents detecting multiple peaks per groove

EDGE_EXCLUSION_PERIODS = 0.6  # Reject grooves within this many periods of either
                              # end of the scan line. Such a groove is real, but
                              # the scan starts or stops part-way through it, so
                              # its facet is clipped and the fitted angle is
                              # meaningless. Set to 0 to disable.
```

**Why edge exclusion matters**: without it, a groove sitting at the start of a scan line
is fitted over a 5-6 pixel sliver instead of a full ~50 pixel facet, returning angles
around 2-3° that are then averaged in with real ~30° measurements. It also equalises
samples: how many edge grooves a scan has is pure luck of where the scan started
relative to the grating phase, so leaving them in weights samples differently for no
physical reason.

### Analysis Parameters

#### Facet Selection and Trimming
```python
BLAZE_SIDE = 'negative_slope'  # Which facet to measure. Selected by slope sign,
                               # not by position, so it is orientation-independent:
                               #   'negative_slope' - the down-sloping facet
                               #   'positive_slope' - the up-sloping facet
                               #   'longer'         - whichever facet is wider
                               # Any other value raises ValueError.

FACET_TRIM = 0.1         # Fraction to trim from facet edges (0.0-0.28)
                         # CRITICAL: Prevents measuring groove tops!
                         # 0.15 = use middle 70% of facet
                         # Increase if getting very low angles (measuring tops)
                         # Decrease if fits look poor (too little data)
```

**Understanding FACET_TRIM**: This is your main defense against measuring groove tops:
- Too small (0.05): Risk measuring rounded groove tops → low angles
- Just right (0.10-0.20): Measures clean facet regions → correct angles
- Too large (0.25+): Few points left to fit → noisy angles

**Hard limit ~0.286.** The blaze facet is trimmed 2.5× harder on the trough side
to avoid the flattened groove bottom, so `FACET_TRIM × 3.5` is removed in total.
Above ~0.286 nothing is left, every groove fit fails, and the analysis returns no
measurements at all rather than degrading gracefully.

### Visualization Options

```python
# Diagnostic plots (recommended for first-time analysis)
SHOW_2D_IMAGE = True                  # 2D topography map
SHOW_FLATTENING_DIAGNOSTIC = True     # Before/after flattening (uses averaged profile)
SHOW_FULL_PROFILE = True              # Full profile with detected grooves
SHOW_INDIVIDUAL_GROOVES = False       # Per-groove fit plots (verbose!)
SHOW_ANALYZED_REGIONS = True          # Overlay of analyzed facet regions

# Statistical plots (always useful)
SHOW_LOCAL_ANGLE_DISTRIBUTION = True  # Within-facet curvature histograms
```

---

## Analysis Modes

### Single File Analysis
Analyze one AFM scan with full diagnostic plots - **recommended for first time**:

```python
ANALYSIS_MODE = 'single'
SINGLE_FILE = 'data/master.txt'
```

**Outputs:**
- Console: Detailed statistics with uncertainty breakdown
- Plots: Diagnostics, histograms (per-groove + facet curvature), spatial variation
- Files: `results/analysis_summary_*.txt`, `analysis_data_*.csv`, `per_groove_data_*.csv`

### Multiple File Analysis
Batch process all files matching a pattern:

```python
ANALYSIS_MODE = 'multiple'
FILE_PATTERN = '*.txt'   # matched inside data/, so no 'data/' prefix
```

**Outputs:**
- Individual analysis for each file
- Combined comparison plots with significance tests
- Batch CSV with all samples

### Comparison Analysis (Temperature Studies)
Compare specific samples with automatic temperature grouping:

```python
ANALYSIS_MODE = 'compare'
SAMPLES_TO_COMPARE = [
    ('data/master.txt', 'Master', None),
    ('data/150C_scan1.txt', '150°C', 150),
    ('data/150C_scan2.txt', '150°C', 150),  # Auto-combined with above
    ('data/175C_scan1.txt', '175°C', 175),
    ('data/175C_scan2.txt', '175°C', 175),
    ('data/200C.txt', '200°C', 200),
]
```

**Features:**
- Multiple scans at same temperature are automatically combined
- Temperature-dependent rate calculations with significance tests
- Master vs. treated comparisons
- Publication-ready comparison plots

### ICC Diagnostic (row-group correlation)

Measures how correlated the row-group measurements actually are, and reports how far
that puts the current standard errors off. Diagnostic only — it changes no reported
number.

```python
ANALYSIS_MODE = 'icc'
```

It runs every entry in `SAMPLES_TO_COMPARE` and prints a per-scan table of ICC,
effective sample size and the SEM inflation factor, saving it to
`results/icc_report_<timestamp>.txt`.

**Output on the bundled sample data**: ICC 0.097–0.429, median 0.244; SEM understated
by 1.17–1.65×.

### PCGrate Boundary Profile Export

A different output from the same measurement. Instead of fitting facet angles,
this averages the detected grooves of one scan into a single representative
groove, normalises it to one period, and writes a PCGrate `.ggp` boundary profile
for grating-efficiency modelling.

```python
ANALYSIS_MODE = 'ggp'
GGP_SOURCE_FILE = 'data/TASTE_ALS_A205_Ti_Pt_flatten.txt'
GGP_N_POINTS = 2000          # Points in the exported profile
GGP_APPLY_SMOOTHING = True   # Light smoothing to remove interpolation kinks
GGP_SMOOTHING_WINDOW = 5
```

**Outputs:**
- `results/averaged_groove_profile_<sample>.ggp` — the boundary profile
- `results/groove_analysis_metrics_<sample>.txt` — depth, RMS slope, max sidewall
  angle and curvature of the normalised profile

**The `.ggp` format** is two header lines then one `x y` pair per line, both
normalised — x runs 0 to 1 across exactly one period, y is height as a fraction of
that period:

```
3 0 - Polygonal type
Period: 1 PSC: 1
0.000000 0.000000
0.000500 0.000001
```

The header must **not** be commented. Writing these with `np.savetxt(header=...)`
prepends `# `, which PCGrate rejects — the writer in `afm_analysis/io/ggp.py`
always emits the correct form.

---

## Understanding Your Results

### Console Output

```
Sample           File         N    Mean ± SEM         σ_total    
Master           master.txt   78   30.15 ± 0.021°     0.184      
150°C            150C.txt     82   30.32 ± 0.019°     0.172      

Note: SEM = Standard Error of Mean (uncertainty in mean estimate)
      σ_total = Total standard deviation (measurement + physical variation)

STATISTICAL COMPARISONS:
Master vs 150°C:
  Difference: +0.170° ± 0.028° (SE)
  95% CI: [+0.114°, +0.226°]
  t-statistic: 6.07 (df=157.8)
  p-value: 0.0000 ***  (highly significant)
  Effect size (Cohen's d): 0.95 (large)
```

**Interpretation:**
- **Mean ± SEM**: Best estimate of blaze angle with uncertainty
- **σ_total**: How much individual measurements vary (groove-to-groove + measurement noise)
- **p < 0.001 (\*\*\*)**: Very confident the difference is real, not random
- **95% CI**: True difference very likely between 0.114° and 0.226°
- **Cohen's d > 0.8**: Large practical effect

### Visualizations

#### 1. Bar Chart (Mean Comparison)
- **Error bars**: ±1 SEM (uncertainty in the mean)
- **Bar labels**: Show mean, SEM, and σ_total
- **Use for**: Comparing mean angles between samples

#### 2. Per-Groove Histogram
- **Distribution**: Individual groove measurements
- **Red shading**: ±1 SEM (narrow - uncertainty in mean)
- **Blue shading**: ±1σ (wide - groove-to-groove variation)
- **Use for**: Understanding measurement spread and precision

#### 3. Facet Curvature Histogram
- **Distribution**: Local angles along each facet (within-facet variation)
- **Shows**: How much blaze angle varies along the facet (camber)
- **Interpretation**: 
  - Narrow distribution (σ < 0.5°): Flat, well-made facets
  - Wide distribution (σ > 1°): Curved/cambered facets
- **Use for**: Assessing manufacturing quality and thermal effects on facet shape

#### 4. Flattening Diagnostic
- **Shows**: Averaged AFM profile before and after background removal
- **Purpose**: Verify that flattening preserves groove structure
- **Critical**: Even with row-group analysis, check this plot to ensure proper flattening!

#### 5. Spatial Variation (Row-Group Analysis)
- **Left panel**: Mean angle per row group (bars)
- **Right panel**: All measurements colored by row group
- **Use for**: Identifying spatial trends or non-uniformity across sample

---

## Troubleshooting

### Problem: Blaze Angles Too Low (~17° instead of ~30°)

**Causes:**
1. ❌ **Wrong `PERIOD_EST`**: Using 833 nm instead of actual spacing (e.g., 315 nm)
2. ❌ **Measuring groove tops**: `FACET_TRIM` too small, measuring rounded peaks
3. ❌ **Wrong facet**: Measuring steep side instead of blaze side

**Solutions:**
```python
# Fix 1: Set correct groove spacing
PERIOD_EST = 315.15  # Match your actual grating!

# Fix 2: Increase trim to avoid groove tops
FACET_TRIM = 0.20    # Use middle 60% of facet

# Fix 3: Switch facet
BLAZE_SIDE = 'positive_slope'  # Try the other facet
```

### Problem: No Grooves Detected

**Causes:**
1. ❌ Poor flattening - background not removed properly
2. ❌ `PROMINENCE_FACTOR` too high
3. ❌ Wrong `PERIOD_EST` - expecting grooves too far apart

**Solutions:**
```python
# Fix 1: Adjust flattening
FLATTEN_METHOD = 'level_grooves'  # Try alternative method
FLATTEN_POLY_ORDER = 2            # Lower order = gentler removal

# Fix 2: Lower detection threshold
PROMINENCE_FACTOR = 0.2

# Fix 3: Check expected spacing
PERIOD_EST = 315.15  # Verify this matches your grating
```

**Diagnostic**: Check `SHOW_FLATTENING_DIAGNOSTIC = True` to verify grooves are visible after flattening.

### Problem: High Variability in Measurements

**Causes:**
1. Noisy AFM data
2. Poor flattening
3. Too many row groups for data quality

**Solutions:**
```python
# Reduce number of row groups
N_ROW_GROUPS = 10  # Instead of 20-50

# Increase facet trim (use cleaner central regions)
FACET_TRIM = 0.20

# Check flattening
SHOW_FLATTENING_DIAGNOSTIC = True
```

### Problem: Measuring Wrong Facet

**Symptom**: Angles near 60° instead of 30° (measuring steep side)

**Solution:**
```python
BLAZE_SIDE = 'positive_slope'  # Switch to other facet
```

### Problem: Error Bars Seem Too Small

**This is usually CORRECT!** Row-group analysis gives proper small error bars because:
- SEM = σ / √N  
- With N=80 instead of N=4, SEM is ~4.5× smaller
- This is real increased precision from using all your data

**Verify**:
- σ_total should be similar between traditional and row-group (~0.2°)
- SEM should scale as 1/√N
- If σ_total is unrealistically small (<0.05°), check data quality

---

## Statistical Concepts

### Standard Error of Mean (SEM)
**Definition**: Uncertainty in the mean value estimate  
**Formula**: `SEM = σ_total / √N`  
**Use**: Comparing mean values between samples  
**Shown**: Error bars on bar charts

### Standard Deviation (σ)
**Definition**: Spread of individual measurements  
**Types**:
- `σ_total`: Includes measurement uncertainty + physical variation
- `std_angle`: Physical variation only (groove-to-groove)
**Use**: Describing consistency of measurements  
**Shown**: Shaded regions in histograms

### Within-Facet Variation (Facet Curvature/Camber)
**Definition**: How much the blaze angle varies along a single facet  
**Calculated**: Local slopes in sliding windows along facet  
**Interpretation**:
- Small (< 0.5°): Flat facets (good manufacturing)
- Large (> 1°): Curved/cambered facets (thermal stress, manufacturing variation)  
**Shown**: Dedicated "Facet Curvature" histogram

### Significance Testing
- **p < 0.05 (*)**: Significant
- **p < 0.01 (**)**: Very significant
- **p < 0.001 (***)**: Highly significant

**95% Confidence Interval**: Range containing true difference with 95% probability

**Effect Size (Cohen's d)**:
- Small: d < 0.2
- Medium: d ≈ 0.5  
- Large: d > 0.8

---

## Output Files

### Text Summary (`results/analysis_summary_*.txt`)
- Configuration parameters used
- Individual sample results with uncertainty breakdown
- Pairwise comparisons with p-values and confidence intervals
- Temperature-dependent analysis
- Variance decomposition

### CSV Data (`results/analysis_data_*.csv`)
Columns include:
- `filename`, `label`, `temperature`
- `mean_angle`, `sem`, `std_angle`, `total_std`
- `n_grooves`, `n_groups` (for row-group analysis)
- `period_nm`, `period_std`
- `mean_measurement_uncertainty`, `within_image_std`

### Per-Groove Data (`results/per_groove_data_*.csv`)
Individual groove measurements:
- `groove_id`, `blaze_angle`, `angle_uncertainty`
- `blaze_r2`, `groove_depth_nm`
- `blaze_width_nm`, `steep_width_nm`
- `local_period_nm`

---

## Project Structure

```
afm_blaze_meas/
├── pyproject.toml                 # Packaging; Qt is an optional [gui] extra
├── main.py                        # Entry point (same dispatch as afm-analysis)
├── .python-version                # Python 3.12.9 (pyenv)
├── data/                          # Your AFM data files
├── results/                       # Analysis outputs (git-ignored, regenerable)
├── docs/
│   ├── BASELINES.md               # How to prove a change altered nothing
│   ├── PROGRESS_SUMMARY.md        # Status and open work
│   └── history/                   # Superseded guides (reference only)
├── tests/                         # pytest; GUI tests skip without PySide6
├── examples/
└── src/afm_analysis/
    ├── wiki/                      # In-app documentation (Wiki tab, package data)
    ├── config.py                  # ⚙️ CONFIGURE HERE
    ├── settings.py                # AnalysisSettings - config as a passed value
    ├── cli.py                     # Mode dispatch
    ├── analyzer.py  workflows.py
    ├── core/                      # SHARED FRONT-END - used by both outputs
    │   ├── processing.py          # Loading, flattening, groove detection
    │   └── analysis.py            # Blaze angle extraction with uncertainty
    ├── boundary/                  # Back-end: PCGrate boundary profiles
    ├── io/                        # Result CSV/TXT, and the .ggp writer
    ├── data/  stats/              # Scan aggregation; comparisons, p-values, ICC
    ├── visualization/
    └── gui/                       # PySide6 front-end
        ├── state.py               # PURE - form validation, no toolkit import
        └── qt/                    # The ONLY modules importing PySide6
            ├── app.py  main_window.py  canvas.py  worker.py
```

The `gui/` vs `gui/qt/` split is enforced by `tests/test_gui_boundary.py`, which
fails if a Qt binding is imported outside `gui/qt`. It caught a violation on its
first run.

### Running the tests

```bash
.venv/bin/python tests/test_ggp_equivalence.py
```

This is the only test suite in the project. It pins the boundary-profile port to
the output of the standalone script it replaced, so a refactor cannot silently
change the exported profile.

For the blaze path there is no test suite; the regression check is a diff of
`results/analysis_data_*.csv` against a known-good run. See `docs/BASELINES.md`.

### Known limitation

Row-group analysis produces many measurements of the *same* physical grooves, but the
statistics treat them as independent samples, so reported SEMs and p-values are
optimistic. Mean angles are unaffected.

This has now been **measured** rather than merely suspected. Run
`ANALYSIS_MODE = 'icc'` to reproduce: across the eight sample scans the intraclass
correlation is 0.097–0.429 (median 0.244), meaning reported standard errors are
understated by **1.17× to 1.65×** and effective sample sizes are 37–70 rather than
78–102. Large effects survive this comfortably; differences below a few tenths of a
degree do not.

The correction is applied: `stats/icc.py` computes the ICC and effective sample
size, and row-group results carry `icc`, `n_effective`, `design_effect`, and
`sem_corrected` alongside the uncorrected `sem`.

---

## Advanced Usage

### Programmatic Access

```python
from afm_analysis import analyze_single_file

# Analyze with row-group analysis
result = analyze_single_file('data/sample.txt', show_plots=False)

print(f"Mean: {result['mean_angle']:.2f}° ± {result['sem']:.3f}° SEM")
print(f"Total std: {result['total_std']:.2f}°")
print(f"Measurements: N={result['n_grooves']} from {result['n_groups']} regions")
print(f"Facet curvature: {result['local_angle_std']:.2f}° (within-facet variation)")
```

### Custom Comparison with Detailed Statistics

```python
from afm_analysis.workflows import run_comparison_analysis
from afm_analysis.stats.analysis import print_uncertainty_breakdown

# Run comparison
results = run_comparison_analysis()

# Get detailed uncertainty breakdown
print_uncertainty_breakdown(results, labels)

# Access individual components
for r in results:
    print(f"{r['label']}:")
    print(f"  Measurement uncertainty: {r['mean_measurement_uncertainty']:.3f}°")
    print(f"  Physical variation: {r['std_angle']:.3f}°")
    print(f"  Within-image variation: {r.get('within_image_std', 0):.3f}°")
    print(f"  Facet curvature (mean): {r['local_angle_std']:.3f}°")
```

### Dual Histogram Plots (Per-Groove + Facet Curvature)

```python
from afm_analysis.visualization.statistics import plot_multi_file_comparison_with_local_angles

# Show both distributions side-by-side
plot_multi_file_comparison_with_local_angles(results, labels, temperatures)
```

---

## Best Practices

### 1. First-Time Analysis Workflow

```python
# config.py - Start with full diagnostics
ANALYSIS_MODE = 'single'
SINGLE_FILE = 'data/test_sample.txt'

# Set YOUR grating specs
PERIOD_EST = 315.15  # YOUR groove spacing!

# Enable all diagnostics
SHOW_2D_IMAGE = True
SHOW_FLATTENING_DIAGNOSTIC = True  # CHECK THIS FIRST!
SHOW_FULL_PROFILE = True
SHOW_ANALYZED_REGIONS = True

# Start with row-group analysis
USE_ROW_GROUPS = True
N_ROW_GROUPS = 20

# Conservative facet trim
FACET_TRIM = 0.15
```

**Check:**
1. Flattening diagnostic - grooves should be clear after background removal
2. Full profile - grooves properly detected
3. Mean angle - should be near expected (~30°, not ~17°)
4. Facet regions overlay - trimmed regions should be clean facet portions

### 2. Temperature Study Workflow

```python
# config.py
ANALYSIS_MODE = 'compare'
SAMPLES_TO_COMPARE = [
    ('data/master_1.txt', 'Master', None),
    ('data/master_2.txt', 'Master', None),  # Multiple masters combined
    ('data/150C_1.txt', '150°C', 150),
    ('data/150C_2.txt', '150°C', 150),
    # ... more temperatures
]

# Optimize for statistics
USE_ROW_GROUPS = True
N_ROW_GROUPS = 20  # More measurements = better p-values

# Consistent analysis
FACET_TRIM = 0.15
BLAZE_SIDE = 'negative_slope'
```

### 3. Publication-Ready Analysis

- Use row-group analysis for maximum precision
- Report mean ± SEM for comparisons
- Include p-values and confidence intervals
- Show both per-groove and facet curvature distributions
- Document all parameters in config.py
- Save all plots and CSV files

**Example Statement**:
> "Blaze angles were measured using AFM topography with row-group analysis (N=20 spatial regions per sample). The master grating showed a mean blaze angle of 30.15° ± 0.02° (SEM, n=78 measurements). Heat treatment at 150°C resulted in a significant increase of 0.17° ± 0.03° (p < 0.001, 95% CI [0.11°, 0.23°], Cohen's d = 0.95). Within-facet angle variation (facet curvature) was 0.18° ± 0.01° for both samples, indicating no change in facet shape."

---

## Requirements

- Python 3.7+
- NumPy - Numerical computations
- Matplotlib - Visualization
- SciPy - Signal processing, peak detection, statistical tests

---

## Citation

If you use this software in your research, please cite:

```
[Your citation information here]
```

---

## License

[Your license here - e.g., MIT, GPL, etc.]

---

## Support

For questions or issues:
- Check troubleshooting section above
- Review diagnostic plots (especially flattening)
- Verify PERIOD_EST matches your sample
- Open an issue on GitHub with diagnostic plots attached

---

## Version History

### v2.0.0 (2025-02)
- **Major**: Row-group analysis implementation for 10-100× more measurements per image
- **Major**: Proper statistical analysis (SEM, p-values, confidence intervals, effect sizes)
- **Major**: Facet curvature (camber) quantification with dedicated histograms
- **Improved**: Flattening diagnostic uses averaged profile for verification
- **Fixed**: Updated default PERIOD_EST to 315.15 nm (was incorrectly 833 nm)
- **Fixed**: Updated expected angles to ~30° (was ~17°)
- **Added**: Automatic groove top avoidance via FACET_TRIM
- **Added**: Variance decomposition and uncertainty breakdown
- **Added**: Publication-ready visualization with clear error bar labels

### v1.0.0 (2025-01)
- Initial release
- Single file, multiple file, and comparison analysis modes
- Basic temperature grouping and scan combination
- Traditional row-averaging analysis

---

**Remember**: Always verify that `PERIOD_EST` (groove spacing) matches your actual grating before running analysis! The defaults are examples only.
