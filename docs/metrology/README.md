# AFM Blaze Angle Analysis (in progress)

A Python toolkit for analyzing Atomic Force Microscopy (AFM) data of diffraction
gratings. A single pipeline produces two outputs:

- **Blaze angles**: facet angles measured per groove, with row-group analysis to
  obtain many measurements per image and statistical comparison across samples.
- **PCGrate boundary profiles**: the averaged groove shape, normalised to one
  period and exported as `.ggp` for grating-efficiency modelling.

Both outputs share the same front-end (loading, flattening, groove detection and
groove windowing) and diverge only at the final step. The sharing is deliberate.
The two outputs were once produced by separate scripts, each with its own copy of
that code, and the copies drifted apart until a bug fixed in one persisted
unnoticed in the other.

## Features

- **High-precision blaze angle extraction** from AFM topography data.
- **Row-group analysis**: 50-100+ measurements are extracted per image, compared
  with 4-10 from traditional averaging.
- **Temperature studies**: scans are aggregated automatically and analyzed as a
  function of temperature.
- **Comprehensive statistics**: uncertainty is quantified with SEM, confidence
  intervals and significance tests.
- **Facet curvature analysis**: the variation of angle within each blaze facet
  (camber) is quantified.
- **Quality control**: automatic filtering avoids measuring groove tops and
  ensures accurate facet detection.
- **Rich outputs**: detailed reports, CSV data and publication-ready
  visualizations.

## Key Capabilities

### 1. Row-Group Analysis for More Measurements
Traditional AFM analysis averages all rows into a single profile and so extracts only 4-10 measurements per image. This package uses *row-group analysis*, in which the image is divided into multiple horizontal bands. This provides:

- **50-100+ measurements** per image instead of 4-10
- **Spatial variation mapping** across the sample
- **Smaller error bars**, although the reduction is not by √N (see below)

The bands are not independent, and the statistics account for this. Each band
re-measures the same physical grooves, so the SEM does not scale as 1/√N. The
measured intraclass correlation is 0.097–0.429 (median 0.244), which places the
effective sample size at 37–70 instead of 78–102.

Every standard error, confidence interval and p-value is computed on the
effective sample size, `N_eff = N / (1 + (m−1)·ICC)`. As a result, the reported
SEMs are 1.17–1.65× wider than the naive value. `ANALYSIS_MODE = 'icc'`
reproduces the diagnostic, and `ICC`, `N_eff`, `SEM_deg` and `SEM_corrected_deg`
are columns in every `analysis_data_*.csv`.

Applying the correction changed no conclusion on the bundled samples. The
master-versus-treated comparison and the 500°C step remain highly significant,
and the steps between adjacent temperatures remain null, slightly more clearly
than before. Effect sizes (Cohen's d) deliberately retain the raw counts, since a
standardised mean difference describes spread and is not used for inference.

### 2. Facet Curvature (Camber) Quantification
In addition to the overall blaze angle of each groove, the software quantifies the *within-facet angle variation* by analyzing local slopes along each facet. This reveals:

- Facet curvature/camber (how much the angle varies along the facet)
- Manufacturing quality (consistent vs. variable facets)
- Thermal effects on facet shape
- Shown in dedicated histograms and statistics

### 3. Robust Groove Detection with Top Avoidance
Accurate measurements are ensured by:

- Automatically detecting groove centers via prominence-based peak finding
- Rejecting grooves clipped by the scan edge, whose facet is only partly present
- Trimming facet edges to avoid measuring groove tops/bottoms
- Recording a fit R² for every groove, in the console output and the per-groove CSV
- Visual diagnostics to verify proper flattening and groove detection

The R² value is reported but not enforced: no measurement is currently dropped
for a poor fit. In practice the edge-exclusion rule removes such fits anyway;
after it is applied, no fit in the sample dataset falls below R² 0.95.

### 4. Statistical Analysis
All comparisons include:

- Standard Error of Mean (SEM): the uncertainty in the mean estimate
- 95% confidence intervals
- Significance tests (p-values) with Welch's t-test
- Effect sizes (Cohen's d) for practical interpretation
- Variance decomposition (measurement vs. physical variation)

---

## Quick Start

### Installation

Python 3.12 is required (pinned in `.python-version`).

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

This installs the analysis, the PCGrate export and the ICC diagnostic, and it
installs no GUI toolkit. The Qt window is an optional extra:

```bash
.venv/bin/pip install -e '.[dev,gui]'
```

The GUI uses PySide6 (LGPL) in place of PyQt (GPL), so that this BSD-3 package
has no GPL dependency. Earlier versions listed PyQt5 as a hard requirement, which
placed GPL code in the install path of anyone who wanted only the analysis.

The package has been verified against numpy 2.5.1 / scipy 1.18.0 / matplotlib
3.11.1; see `docs/BASELINES.md` for how to reproduce a stored result exactly.

### Basic Usage

1. **Configure the analysis** in `afm_analysis/config.py`:
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

   The Qt window is `afm-gui`. It has four tabs: **Import** (load and flatten a
   scan), **Analysis** (the blaze-angle results), **Boundary** (PCGrate `.ggp`
   export), and a **Wiki** that explains how the statistics work, including the
   correlation correction, the reason edge grooves are rejected, and the meaning
   of each output column.

   Paths in `config.py` resolve against the project root, so the analysis can be
   run from any working directory.

3. **Review the results**:
   - **Console output**: statistical summary with SEM, p-values and significance
   - **Plots**: bar charts (with SEM error bars) and histograms (per-groove and facet curvature)
   - **Files**: the `results/` folder contains detailed reports and CSV data

---

## Configuration Guide

### Parameters to Set First

#### Grating Specifications
```python
# Groove spacing - MUST MATCH YOUR SAMPLE
PERIOD_EST = 315.0   # nm (e.g., ~315 nm for typical echelle gratings)
                     # NOT 833 nm! That's for a different grating type

SCAN_X_SIZE = 2.0    # um - fallback scan width, used only when the width
                     # cannot be read from the data file header
```

#### Row-Group Analysis (Recommended)
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

The flattening diagnostic plot shows the averaged profile, so that the
background removal can be verified even when row-group analysis is used. This
check ensures that all subsequent measurements are accurate.

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

Edge exclusion serves two purposes. Without it, a groove at the start of a scan
line is fitted over a 5-6 pixel sliver instead of a full ~50 pixel facet, which
returns angles of around 2-3° that are then averaged together with real ~30°
measurements. It also places samples on an equal footing: the number of edge
grooves in a scan depends only on where the scan started relative to the grating
phase, so retaining them weights samples differently for no physical reason.

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

FACET_TRIM is the main safeguard against measuring groove tops:
- Too small (0.05): rounded groove tops may be measured, giving low angles
- Appropriate (0.10-0.20): clean facet regions are measured, giving correct angles
- Too large (0.25+): few points remain for the fit, giving noisy angles

The parameter has a hard limit of ~0.286. The blaze facet is trimmed 2.5× harder
on the trough side to avoid the flattened groove bottom, so `FACET_TRIM × 3.5` is
removed in total. Above ~0.286 nothing remains, every groove fit fails, and the
analysis returns no measurements at all instead of degrading gradually.

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
This mode analyzes one AFM scan with full diagnostic plots and is recommended for a first analysis:

```python
ANALYSIS_MODE = 'single'
SINGLE_FILE = 'data/master.txt'
```

**Outputs:**
- Console: detailed statistics with an uncertainty breakdown
- Plots: diagnostics, histograms (per-groove and facet curvature), spatial variation
- Files: `results/analysis_summary_*.txt`, `analysis_data_*.csv`, `per_groove_data_*.csv`

### Multiple File Analysis
This mode processes all files matching a pattern as a batch:

```python
ANALYSIS_MODE = 'multiple'
FILE_PATTERN = '*.txt'   # matched inside data/, so no 'data/' prefix
```

**Outputs:**
- Individual analysis for each file
- Combined comparison plots with significance tests
- Batch CSV with all samples

### Comparison Analysis (Temperature Studies)
This mode compares specified samples and groups them by temperature automatically:

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
- Multiple scans at the same temperature are combined automatically
- Temperature-dependent rate calculations with significance tests
- Master vs. treated comparisons
- Publication-ready comparison plots

### ICC Diagnostic (row-group correlation)

This mode measures how strongly the row-group measurements are correlated and
reports how far that correlation places the current standard errors from their
correct values. It is a diagnostic only and changes no reported number.

```python
ANALYSIS_MODE = 'icc'
```

It runs every entry in `SAMPLES_TO_COMPARE`, prints a per-scan table of ICC,
effective sample size and SEM inflation factor, and saves the table to
`results/icc_report_<timestamp>.txt`.

On the bundled sample data the output is an ICC of 0.097–0.429 (median 0.244),
with the SEM understated by 1.17–1.65×.

### PCGrate Boundary Profile Export

This mode produces a different output from the same measurement. In place of
fitting facet angles, it averages the detected grooves of one scan into a single
representative groove, normalises it to one period, and writes a PCGrate `.ggp`
boundary profile for grating-efficiency modelling.

```python
ANALYSIS_MODE = 'ggp'
GGP_SOURCE_FILE = 'data/TASTE_ALS_A205_Ti_Pt_flatten.txt'
GGP_N_POINTS = 2000          # Points in the exported profile
GGP_APPLY_SMOOTHING = True   # Light smoothing to remove interpolation kinks
GGP_SMOOTHING_WINDOW = 5
```

**Outputs:**
- `results/averaged_groove_profile_<sample>.ggp`: the boundary profile
- `results/groove_analysis_metrics_<sample>.txt`: depth, RMS slope, maximum
  sidewall angle and curvature of the normalised profile

The `.ggp` format is one `x y` pair per line, with no header and no line
break after the last point. Both coordinates are normalised: x runs from 0 to 1
across exactly one period, and y is the height as a fraction of that period.

```
0.000000 0.000000
0.000500 0.000001
```

PCGrate's border-profile import takes the points alone (Aditi, 2026-10-29).
A script that needs the legacy `3 0 - Polygonal type` / `Period: 1 PSC: 1`
lines can call `gratinglab.io.ggp.write_ggp(..., header=True)`; the GUI and
the CLI never write them. Never write these files with
`np.savetxt(header=...)`, which prepends `# `.

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
- **Mean ± SEM**: the best estimate of the blaze angle, with its uncertainty
- **σ_total**: the variation of individual measurements (groove-to-groove variation plus measurement noise)
- **p < 0.001 (\*\*\*)**: high confidence that the difference is real and not due to chance
- **95% CI**: the true difference very likely lies between 0.114° and 0.226°
- **Cohen's d > 0.8**: a large practical effect

### Visualizations

#### 1. Bar Chart (Mean Comparison)
- **Error bars**: ±1 SEM (uncertainty in the mean)
- **Bar labels**: mean, SEM and σ_total
- **Use for**: comparing mean angles between samples

#### 2. Per-Groove Histogram
- **Distribution**: individual groove measurements
- **Red shading**: ±1 SEM (narrow; uncertainty in the mean)
- **Blue shading**: ±1σ (wide; groove-to-groove variation)
- **Use for**: understanding the spread and precision of the measurements

#### 3. Facet Curvature Histogram
- **Distribution**: local angles along each facet (within-facet variation)
- **Shows**: how much the blaze angle varies along the facet (camber)
- **Interpretation**: 
  - Narrow distribution (σ < 0.5°): flat, well-made facets
  - Wide distribution (σ > 1°): curved or cambered facets
- **Use for**: assessing manufacturing quality and thermal effects on facet shape

#### 4. Flattening Diagnostic
- **Shows**: the averaged AFM profile before and after background removal
- **Purpose**: to verify that flattening preserves the groove structure
- **Note**: this plot should be checked to confirm proper flattening even when row-group analysis is used

#### 5. Spatial Variation (Row-Group Analysis)
- **Left panel**: mean angle per row group (bars)
- **Right panel**: all measurements colored by row group
- **Use for**: identifying spatial trends or non-uniformity across the sample

---

## Troubleshooting

### Problem: Blaze Angles Too Low (~17° instead of ~30°)

**Causes:**
1. **Wrong `PERIOD_EST`**: 833 nm is used instead of the actual spacing (e.g., 315 nm)
2. **Measuring groove tops**: `FACET_TRIM` is too small, so rounded peaks are measured
3. **Wrong facet**: the steep side is measured instead of the blaze side

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
1. Poor flattening: the background has not been removed properly
2. `PROMINENCE_FACTOR` is too high
3. Wrong `PERIOD_EST`: grooves are expected too far apart

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

**Diagnostic**: set `SHOW_FLATTENING_DIAGNOSTIC = True` to verify that the grooves are visible after flattening.

### Problem: High Variability in Measurements

**Causes:**
1. Noisy AFM data
2. Poor flattening
3. Too many row groups for the quality of the data

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

**Symptom**: angles near 60° instead of 30°, indicating that the steep side is being measured

**Solution:**
```python
BLAZE_SIDE = 'positive_slope'  # Switch to other facet
```

### Problem: Error Bars Seem Too Small

Small error bars are usually correct. Row-group analysis gives smaller error bars
than the traditional analysis because it uses all of the data. The naive SEM,
SEM = σ / √N, would be ~4.5× smaller with N=80 than with N=4, but the row groups
re-measure the same grooves and are therefore not independent, so the reported
error does not shrink as fast as that. `SEM_corrected_deg` uses the effective
sample size from the ICC correction (see "Known limitation" below), and
comparisons between samples use it together with the effective sample size.

**Verify** that:
- σ_total is similar between traditional and row-group analysis (~0.2°)
- `SEM_corrected_deg` is wider than `SEM_deg`, by more when the row groups are strongly correlated
- σ_total is not unrealistically small (<0.05°); if it is, the data quality should be checked

---

## Statistical Concepts

### Standard Error of Mean (SEM)
**Definition**: the uncertainty in the estimate of the mean value  
**Formula**: `SEM = σ_total / √N`  
**Use**: comparing mean values between samples  
**Shown**: error bars on bar charts

### Standard Deviation (σ)
**Definition**: the spread of individual measurements  
**Types**:
- `σ_total`: includes measurement uncertainty and physical variation
- `std_angle`: physical variation only (groove-to-groove)
**Use**: describing the consistency of measurements  
**Shown**: shaded regions in histograms

### Within-Facet Variation (Facet Curvature/Camber)
**Definition**: how much the blaze angle varies along a single facet  
**Calculated**: from local slopes in sliding windows along the facet  
**Interpretation**:
- Small (< 0.5°): flat facets (good manufacturing)
- Large (> 1°): curved or cambered facets (thermal stress, manufacturing variation)  
**Shown**: a dedicated "Facet Curvature" histogram

### Significance Testing
- **p < 0.05 (*)**: Significant
- **p < 0.01 (**)**: Very significant
- **p < 0.001 (***)**: Highly significant

**95% Confidence Interval**: the range containing the true difference with 95% probability

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
The columns include:
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

The separation between `gui/` and `gui/qt/` is enforced by
`tests/test_gui_boundary.py`, which fails if a Qt binding is imported outside
`gui/qt`. The test found a violation on its first run.

### Running the tests

```bash
.venv/bin/python -m pytest tests/ -v
```

The suite contains 135 tests, covering the blaze-angle pipeline
(`tests/test_analyzer.py`), the wiring of the ICC correction
(`tests/test_icc_correction.py`), the boundary/PCGrate export, and the Qt GUI
(skipped automatically if PySide6 is not installed). CI runs the same suite, plus
`ruff check .`, on every push and pull request.

`tests/test_ggp_equivalence.py` is the oldest of these tests. It pins the
boundary-profile port to the output of the standalone script it replaced, so that
a refactor cannot change the exported profile without detection.

The blaze-angle path also has a manual baseline-diffing procedure
(`docs/BASELINES.md`) that predates the test suite above. It remains useful for
proving that a refactor changed nothing, since a diff against a stored CSV is a
stronger guarantee than any set of assertions, but it is no longer the only
regression check on that path.

### Known limitation

Row-group analysis produces many measurements of the same physical grooves, so
statistics that treated them as independent samples would report optimistic SEMs
and p-values. Mean angles are unaffected.

This effect has now been measured. It can be reproduced with
`ANALYSIS_MODE = 'icc'`: across the eight sample scans the intraclass correlation
is 0.097–0.429 (median 0.244), which means that the reported standard errors are
understated by 1.17× to 1.65× and that the effective sample sizes are 37–70
instead of 78–102. Large effects remain significant under this correction, while
differences below a few tenths of a degree do not.

The correction is applied: `stats/icc.py` computes the ICC and the effective
sample size, and row-group results carry `icc`, `n_effective`, `design_effect`
and `sem_corrected` alongside the uncorrected `sem`.

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
1. Flattening diagnostic: the grooves should be clear after background removal
2. Full profile: the grooves should be properly detected
3. Mean angle: it should be near the expected value (~30°, not ~17°)
4. Facet regions overlay: the trimmed regions should be clean portions of the facet

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
- Show both the per-groove and the facet curvature distributions
- Document all parameters in config.py
- Save all plots and CSV files

**Example Statement**:
> "Blaze angles were measured using AFM topography with row-group analysis (N=20 spatial regions per sample). The master grating showed a mean blaze angle of 30.15° ± 0.02° (SEM, n=78 measurements). Heat treatment at 150°C resulted in a significant increase of 0.17° ± 0.03° (p < 0.001, 95% CI [0.11°, 0.23°], Cohen's d = 0.95). Within-facet angle variation (facet curvature) was 0.18° ± 0.01° for both samples, indicating no change in facet shape."

---

## Requirements

- Python 3.7+
- NumPy: numerical computations
- Matplotlib: visualization
- SciPy: signal processing, peak detection and statistical tests

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
- Consult the troubleshooting section above
- Review the diagnostic plots (especially flattening)
- Verify that PERIOD_EST matches the sample
- Open an issue on GitHub with the diagnostic plots attached

---

## Version History

### v2.0.0 (2025-02)
- **Major**: row-group analysis implemented, for 10-100× more measurements per image
- **Major**: statistical analysis (SEM, p-values, confidence intervals, effect sizes)
- **Major**: facet curvature (camber) quantification with dedicated histograms
- **Improved**: the flattening diagnostic uses the averaged profile for verification
- **Fixed**: the default PERIOD_EST was updated to 315.15 nm (it was incorrectly 833 nm)
- **Fixed**: the expected angles were updated to ~30° (they were ~17°)
- **Added**: automatic groove top avoidance via FACET_TRIM
- **Added**: variance decomposition and uncertainty breakdown
- **Added**: publication-ready visualization with clear error bar labels

### v1.0.0 (2025-01)
- Initial release
- Single file, multiple file, and comparison analysis modes
- Basic temperature grouping and scan combination
- Traditional row-averaging analysis

---

Before an analysis is run, `PERIOD_EST` (the groove spacing) should be checked
against the actual grating. The default values are examples only.
