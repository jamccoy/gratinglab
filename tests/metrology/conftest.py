"""
Shared test setup.

Puts src/ on the path so the suite runs from a checkout without an install.
An editable install (`pip install -e '.[dev,gui]'`) makes this redundant but not
harmful - the installed package resolves to the same files.
"""
import os
import sys

import matplotlib

# Every test must run without a display. Set before any pyplot import.
matplotlib.use("Agg")

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src')
if SRC not in sys.path:
    sys.path.insert(0, SRC)
