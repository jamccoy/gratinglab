"""The changelog has a section for the version being released.

A release commit moves the "Unreleased" entries under a heading named for the
new version. This test fails on the bump commit if that step was forgotten,
which is the point at which it is cheapest to fix, the same arrangement as the
CITATION.cff test.
"""
import re
from pathlib import Path

import gratinglab

ROOT = Path(__file__).resolve().parent.parent
CHANGELOG = ROOT / "CHANGELOG.md"


def _sections():
    text = CHANGELOG.read_text(encoding="utf-8")
    return re.findall(r"^## (\S+)", text, re.M)


def test_the_current_version_has_a_section():
    assert gratinglab.__version__ in _sections(), (
        f"CHANGELOG.md has no '## {gratinglab.__version__}' section; move the "
        "Unreleased entries under it in the release commit")


def test_unreleased_comes_first_and_versions_descend():
    sections = _sections()
    versions = [s for s in sections if s != "Unreleased"]
    if "Unreleased" in sections:
        assert sections[0] == "Unreleased", "Unreleased must be the first section"
    key = [tuple(int(p) for p in v.split(".")) for v in versions]
    assert key == sorted(key, reverse=True), "versions must be listed newest first"
