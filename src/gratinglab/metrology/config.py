"""
Configuration file for AFM blaze angle analysis
"""
import os

# ============ PROJECT PATHS ============
# All paths below are resolved relative to the project root, not the current
# working directory, so the analysis runs the same from anywhere.
#
# The package lives at <root>/src/gratinglab/metrology. The depth has changed
# before and will again; found by looking for the marker rather than counting parents,
# because counting broke silently when the package moved under src/ - data/ and
# results/ resolved to src/data and src/results, and the failure looked like a
# stale output file rather than a path bug.
#
# An installed (non-editable) copy has no project root above it; fall back to the
# working directory so `pip install gratinglab[metrology]` still runs against a user's
# own data and results folders.
def _find_project_root(start):
    # pyproject.toml only. Searching for a "data" directory looked reasonable and
    # was wrong: this package contains its own data subpackage
    # (gratinglab/metrology/data/aggregation.py), so the walk stopped immediately and
    # DATA_DIR pointed inside the source tree.
    path = os.path.dirname(os.path.abspath(start))
    while True:
        if os.path.isfile(os.path.join(path, 'pyproject.toml')):
            return path
        parent = os.path.dirname(path)
        if parent == path:
            return os.getcwd()
        path = parent


PROJECT_ROOT = _find_project_root(__file__)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')

# Scans live *outside* the checkout, and are not committed.
#
# They are the research group's measurement data, and the same rule already
# applies to the PCGrate reference corpus (see tests/corpus.py and
# GRATINGLAB_REF_DIR). The metrology package arrived from a repository that
# committed its scans instead; keeping both conventions in one project would
# have meant the answer to "is measurement data in the repo?" depended on which
# half you were standing in.
#
# Point GRATINGLAB_AFM_DIR somewhere else to override.
DEFAULT_AFM_DIR = os.path.join(os.path.expanduser('~'), 'Documents', 'afm_scans')
DATA_DIR = os.environ.get('GRATINGLAB_AFM_DIR') or DEFAULT_AFM_DIR


def resolve_path(name):
    """Resolve a configured scan name against DATA_DIR, leaving absolute paths alone.

    Every caller passes a scan filename, and scans no longer live in the
    checkout, so this resolves against DATA_DIR rather than PROJECT_ROOT. An
    absolute path is still honoured untouched, which is how the GUI passes a
    file the user picked from anywhere.
    """
    return name if os.path.isabs(name) else os.path.join(DATA_DIR, name)


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

# ============ IMAGE FLATTENING (2-D, before rows are averaged) ============
# Corrects scan lines relative to one another, on the image as it came off the
# instrument. 'none', 'plane', or 'align_rows'.
#
# Affine methods cannot change a blaze angle in this pipeline: the profile
# flattening below removes any constant or linear term again. Measured 0.0000 deg
# across all eight compare-mode samples. It matters for viewing the image, and as
# the place to add a method that is not affine.
IMAGE_FLATTEN_METHOD = 'align_rows'

# ============ TIP CORRECTION (2-D, after image flattening) ============
# Undo the tip's dilation of the surface where that is possible (Villarrubia
# 1997 erosion; see core/tip.py). Off by default: the correction changes
# measured numbers -- groove depth especially -- and must be asked for. The
# reconstruction is an upper bound on the true surface, exact where the tip
# apex made contact; the certainty fraction is reported wherever the results
# land, and a facet steeper than the tip flank (90 - half angle from the
# surface) is unrecoverable no matter the algorithm.
TIP_CORRECTION = 'none'      # 'none' or 'erosion'
TIP_RADIUS_NM = 1.0          # Apex radius; a "2 nm wide" tip is radius 1
TIP_HALF_ANGLE_DEG = 18.0    # Cone half angle from the tip axis

# ============ PROFILE FLATTENING (1-D, after row averaging) ============
# This is the one that moves the answer: about 0.49 deg between methods.
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

# ============ RAW NANOSCOPE (.spm) INPUT ============
# Used only when the input is a raw .spm rather than a text export. A Nanoscope
# file holds several planes; these choose which one to analyse.
SPM_CHANNEL = 'Height Sensor'   # 'Height Sensor' is the topography channel.
                                # 'Peak Force Error' is not a height map and
                                # cannot yield a blaze angle.
SPM_DIRECTION = 'Retrace'       # 'Retrace' or 'Trace'. Retrace matches the
                                # existing Gwyddion exports in data/.

# ============ PCGRATE BOUNDARY PROFILE EXPORT ============
# Used by ANALYSIS_MODE = 'ggp'. Averages the grooves of one scan into a single
# representative groove and writes it as a PCGrate .ggp boundary profile.
GGP_SOURCE_FILE = 'TASTE_ALS_A205_Ti_Pt_flatten.txt'
GGP_N_POINTS = 2000          # Points in the exported profile
GGP_APPLY_SMOOTHING = True   # Light smoothing to remove interpolation kinks
GGP_SMOOTHING_WINDOW = 5     # Larger = smoother
GGP_MIN_HALF_WIDTH = 10      # Skip grooves whose symmetric half-extent is <= this
                             # many samples (i.e. clipped by the scan edge)

# ============ FILE SELECTION ============
# Choose analysis mode: 'single', 'multiple', 'compare', 'ggp', or 'icc'
ANALYSIS_MODE = 'compare'  # 'single' for one file, 'multiple' for pattern matching,
                           # 'compare' for specific samples, 'ggp' for a PCGrate
                           # profile, 'icc' for the row-group correlation diagnostic

# Single file mode
SINGLE_FILE = '150C_2um_flatten.txt'

# Multiple file mode (pattern matching)
FILE_PATTERN = 'TASTE_*.txt'

# Compare mode - define your samples here
# Each entry: (filename, label for plots, temperature)
SAMPLES_TO_COMPARE = [
    ('ALD_master_1p5um_flatten.txt', 'Master', None),
    #('20250820_150C_00002.txt', '150°C', 150), #only a partial scan, same as below
    ('20250820_150C_00003.txt', '150°C', 150),
    ('20250820_215C_00001.txt', '215°C', 215),
    ('20250820_280C_00004.txt', '280°C', 280),
    ('20250905_280C_00005.txt', '280°C', 280),
    ('20250905_280C_00004.txt', '280°C', 280),
    #('20250905_280C_00003.txt', '280°C', 280), #accidentally measuring some of the groove bottom
    ('20250905_280C_00000.txt', '280°C', 280),
    ('500C_N2_flatten.txt', '500°C', 500)
    # Add more samples as needed:
    # ('file3.txt', '200°C'),
    # ('file4.txt', '225°C'),
]
