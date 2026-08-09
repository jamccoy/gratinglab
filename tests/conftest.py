"""Shared fixtures.

The real PCGrate corpus is reference data belonging to the research group and is
deliberately **not** committed. Tests that need it skip cleanly when it is
absent, so the suite runs anywhere. Point ``GRATINGLAB_REF_DIR`` elsewhere to
override the default location.

This module is also where the Qt environment is settled, for the same reason:
so that plain ``pytest`` works everywhere and CI and a local checkout agree on
what "run the tests" means.
"""

import os
from pathlib import Path

import pytest

# Qt must be told there is no display before anything imports it, or the widget
# tests fail on a headless runner. `setdefault`, not assignment, so that
# `QT_QPA_PLATFORM=cocoa pytest ...` still lets you watch a test drive a real
# window -- which is genuinely the fastest way to debug a layout problem.
#
# QT_API matters because matplotlib's shim (backends/qt_compat.py) picks a
# binding from sys.modules first and this variable second; without it, a
# stray PyQt5 in the environment would win.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_API", "PySide6")

# Re-exported so existing fixtures keep working. The definitions live in
# `corpus.py` because a module-level `from .conftest import ...` breaks under
# a tool that drives pytest in-process from a copied tree -- see that file.
from corpus import DATA, DEFAULT_REF_DIR, reference_dir  # noqa: F401


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
