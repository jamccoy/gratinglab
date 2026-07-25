"""
Configuration file for AFM blaze angle analysis
"""
import os

# ============ PROJECT PATHS ============
# All paths below are resolved relative to the project root, not the current
# working directory, so the analysis runs the same from anywhere.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')


def resolve_path(name):
    """Resolve a configured path against the project root, leaving absolute paths alone"""
    return name if os.path.isabs(name) else os.path.join(PROJECT_ROOT, name)


# ============ SCAN PARAMETERS ============
# Adjust these for your data
SCAN_X_SIZE = 2.0  # um (actual scan width)
PERIOD_EST = 315.0  # nm (estimated grating period)

# ============ ANALYSIS PARAMETERS ============
PROMINENCE_FACTOR = 0.01  # Lower = more sensitive groove detection
DISTANCE_FACTOR = 0.3     # Minimum spacing between grooves (fraction of period)
EDGE_EXCLUSION_PERIODS = 0.6  # Reject grooves within this many periods of either
                              # end of the scan line. Such a groove is real, but
                              # the scan starts/stops part-way through it, so its
                              # facet is clipped and any fitted angle is garbage.
                              # 0 disables (restores the old behaviour).
FACET_TRIM = 0.1          # Trim 10% from each end of facet (reduced from 20% to capture more data)
BLAZE_SIDE = 'negative_slope'  # 'negative_slope', 'positive_slope', or 'longer'

# ============ ROW-GROUP ANALYSIS (NEW!) ============
# Enable this to get many more measurements per AFM image
USE_ROW_GROUPS = True  # Set to True to enable row-group analysis
N_ROW_GROUPS = 20      # Number of row groups to extract from each image
                       # More groups = more measurements (10-50 recommended)
                       # Each group gets averaged separately and analyzed
                       # This gives you N_ROW_GROUPS × N_GROOVES measurements per image!

# ============ FLATTENING OPTIONS ============
FLATTEN_METHOD = 'level_grooves'  # 'linear', 'polynomial', 'groove_peaks', or 'level_grooves'
FLATTEN_POLY_ORDER = 2  # For polynomial and level_grooves: 1=linear, 2=quadratic, 3=cubic
FLATTEN_EXCLUDE_EDGES = 0.05  # Exclude this fraction from each edge when fitting
FLATTEN_FEATURE = 'peaks'  # For level_grooves: 'peaks' (lands), 'troughs' (groove bottoms), or 'both'

# ============ DISPLAY OPTIONS ============
SHOW_2D_IMAGE = False     # Show 2D AFM topography
SHOW_INDIVIDUAL_GROOVES = False  # Show each groove's blaze angle
SHOW_FULL_PROFILE = True  # Show full profile with all grooves overlaid
SHOW_FLATTENING_DIAGNOSTIC = False  # Show before/after flattening comparison (now per-group, can be verbose)
SHOW_LOCAL_ANGLE_DISTRIBUTION = True  # Show distribution of local angles within facets
SHOW_ANALYZED_REGIONS = False  # Show which portions of facets are being analyzed (now per-group, can be verbose)

# ============ FILE SELECTION ============
# Choose analysis mode: 'single', 'multiple', or 'compare'
ANALYSIS_MODE = 'compare'  # 'single' for one file, 'multiple' for pattern matching, 'compare' for specific samples

# Single file mode
SINGLE_FILE = 'data/150C_2um_flatten.txt'

# Multiple file mode (pattern matching)
FILE_PATTERN = 'TASTE_*.txt'

# Compare mode - define your samples here
# Each entry: (filename, label for plots, temperature)
SAMPLES_TO_COMPARE = [
    ('data/ALD_master_1p5um_flatten.txt', 'Master', None),
    #('data/20250820_150C_00002.txt', '150°C', 150), #only a partial scan, same as below
    ('data/20250820_150C_00003.txt', '150°C', 150),
    ('data/20250820_215C_00001.txt', '215°C', 215),
    ('data/20250820_280C_00004.txt', '280°C', 280),
    ('data/20250905_280C_00005.txt', '280°C', 280),
    ('data/20250905_280C_00004.txt', '280°C', 280),
    #('data/20250905_280C_00003.txt', '280°C', 280), #accidentally measuring some of the groove bottom
    ('data/20250905_280C_00000.txt', '280°C', 280),
    ('data/500C_N2_flatten.txt', '500°C', 500)
    # Add more samples as needed:
    # ('file3.txt', '200°C'),
    # ('file4.txt', '225°C'),
]
