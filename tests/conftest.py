"""Shared fixtures.

The real PCGrate corpus is reference data belonging to the research group and is
deliberately **not** committed. Tests that need it skip cleanly when it is
absent, so the suite runs anywhere. Point ``GRATINGLAB_REF_DIR`` elsewhere to
override the default location.
"""

import os
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"

DEFAULT_REF_DIR = Path.home() / "Documents" / "diffraction_efficiency"


def reference_dir() -> Path | None:
    """The PCGrate reference corpus, or ``None`` if unavailable."""
    override = os.environ.get("GRATINGLAB_REF_DIR")
    candidate = Path(override) if override else DEFAULT_REF_DIR
    return candidate if candidate.is_dir() else None


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    directory = reference_dir()
    if directory is None:
        pytest.skip(
            "PCGrate reference corpus not found; set GRATINGLAB_REF_DIR to enable"
        )
    return directory


@pytest.fixture(scope="session")
def synthetic_wavescan() -> Path:
    return DATA / "synthetic_wavescan.txt"
