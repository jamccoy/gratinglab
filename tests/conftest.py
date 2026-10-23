"""Shared fixtures.

The real PCGrate corpus is reference data belonging to the research group and is
deliberately **not** committed. Tests that need it skip cleanly when it is
absent, so the suite runs anywhere. Point ``GRATINGLAB_REF_DIR`` elsewhere to
override the default location.

The Qt and matplotlib environment used to be settled here. It now comes from
``endstation``'s pytest plugin, which is loaded through an entry point and so
runs *before* any conftest is imported -- which is the timing this actually
needs, since both the Qt platform plugin and the matplotlib backend bind at
first import and a root conftest is imported before ``pytest_configure``.

It still uses ``setdefault`` throughout, so ``QT_QPA_PLATFORM=cocoa pytest ...``
still lets you watch a test drive a real window.
"""

import os
import tempfile
from pathlib import Path

import pytest

# Re-exported so existing fixtures keep working. The definitions live in
# `corpus.py` because a module-level `from .conftest import ...` breaks under
# a tool that drives pytest in-process from a copied tree -- see that file.
from corpus import DATA, DEFAULT_REF_DIR, reference_dir  # noqa: F401

# Your own material library (``~/.gratinglab/materials``, or wherever
# GRATINGLAB_MATERIALS points) must not leak into the suite: a test that passes
# because of a table on one machine fails on the next. An empty directory, set
# at import rather than in a fixture, because `parametrize(materials.available())`
# runs at collection, before any fixture would. A test that wants a library
# sets the variable itself with monkeypatch, which restores this afterwards.
os.environ["GRATINGLAB_MATERIALS"] = tempfile.mkdtemp(prefix="gratinglab-no-library-")


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


def pytest_sessionfinish(session, exitstatus):
    """Free the Qt objects the suite leaves behind while Qt still exists.

    Python's final garbage collection runs after Qt has started to shut down,
    and if it finds Qt objects then it can abort the interpreter after every
    test has passed: ``Fatal Python error: bool_dealloc`` with
    "Garbage-collecting" on the stack, which the ubuntu / py3.11 CI job hit
    intermittently. ``gui.qt.gc_guard`` collects at exit for the same reason,
    but it cannot reach windows that were closed and never deleted. Closing and
    deleting every top-level widget here, then collecting, leaves the final
    collection nothing of Qt's to free.

    Nothing is imported for a run that never loaded Qt.
    """
    import gc
    import sys

    if "PySide6.QtWidgets" not in sys.modules:
        return
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        return
    app.closeAllWindows()
    for widget in app.topLevelWidgets():
        widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()
    gc.collect()
