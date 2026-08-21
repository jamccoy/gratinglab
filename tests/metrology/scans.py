"""Where AFM scans come from, for the tests that need one.

Two kinds, and the distinction is the point of this module.

**The synthetic scan** is committed, and is a *known* grating: an ideal sawtooth
of stated period and blaze angle, written by
``tools/metrology/make_synthetic_scan.py``. Most of the suite uses it, because
most of the suite is testing the pipeline rather than a measurement, and a
fixture whose true answer is known lets a test assert recovery instead of merely
absence of crashes.

**Real scans** are the research group's measurement data and are deliberately
not committed -- the same rule the PCGrate corpus already follows in
``corpus.py``. Only tests that genuinely need a real measurement should reach for
them: the Nanoscope binary reader, the byte-level ``rev3`` equivalence pin, and
the ICC statistics, which need real groove-to-groove correlation to mean
anything at all.

Separate from ``conftest.py`` for the same reason ``corpus.py`` is: pytest
imports conftest under a name of its own choosing, so importing *from* it breaks
under a tool that drives pytest in-process from a copied tree. mutmut is one.
"""

import os
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"

#: A known grating: period 315 nm, blaze 30 deg, anti-blaze 70 deg, sharp.
#: Regenerate with ``python tools/metrology/make_synthetic_scan.py``.
SYNTHETIC = FIXTURES / "synthetic_blazed_scan.txt"
SYNTHETIC_PERIOD_NM = 315.0
SYNTHETIC_BLAZE_DEG = 30.0
SYNTHETIC_ANTIBLAZE_DEG = 70.0

#: Real scans, not committed. Override with ``GRATINGLAB_AFM_DIR``. Kept in step
#: with ``gratinglab.metrology.config.DEFAULT_AFM_DIR``.
DEFAULT_SCAN_DIR = Path.home() / "Documents" / "afm_scans"


def scan_dir() -> Path | None:
    """The directory of real scans, or ``None`` if unavailable."""
    override = os.environ.get("GRATINGLAB_AFM_DIR")
    candidate = Path(override) if override else DEFAULT_SCAN_DIR
    return candidate if candidate.is_dir() else None


def real_scan(name: str) -> Path:
    """A named real scan, skipping cleanly when it is not on this machine.

    Skips rather than fails: a contributor without the group's measurement data
    must still get a green suite, exactly as with the PCGrate corpus. The
    message names the environment variable, because "file not found" does not
    tell anyone what to do about it.
    """
    directory = scan_dir()
    if directory is None:
        pytest.skip(
            "AFM scans not found; set GRATINGLAB_AFM_DIR to enable "
            f"(looked in {DEFAULT_SCAN_DIR})"
        )
    path = directory / name
    if not path.exists():
        pytest.skip(f"{name} not in {directory}")
    return path
