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


def test_no_version_has_two_sections():
    """
    One heading per version, including Unreleased.

    Two edits that each anchored where their own entry belonged, without
    looking at what the other had already written, produced a second
    `## Unreleased` above the first. The tests above did not notice:
    Unreleased was still first, and the version ordering ignores it. The
    duplicate would have been released as a version section with half its
    entries missing.
    """
    sections = _sections()
    duplicated = {s for s in sections if sections.count(s) > 1}
    assert not duplicated, f"CHANGELOG.md has two sections for {duplicated}"


def _subsections():
    """The `###` groups under each `##` heading, in order."""
    groups, current = {}, None
    for line in CHANGELOG.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            current = line[3:].split(" ")[0]
            groups[current] = []
        elif line.startswith("### ") and current is not None:
            groups[current].append(line[4:].strip())
    return groups


def test_no_section_lists_a_group_twice():
    """
    Added, Changed and Fixed appear at most once per version.

    The same failure one level down: entries written at different times, each
    under a heading of its own, left a section reading Fixed, Changed, Fixed.
    """
    for version, groups in _subsections().items():
        duplicated = {g for g in groups if groups.count(g) > 1}
        assert not duplicated, f"{version} lists {duplicated} twice"


def test_groups_follow_the_keep_a_changelog_order():
    """
    Added before Changed before Fixed, as every released section here does.

    Not a matter of taste: a reader scanning for what broke should find it in
    the same place in every section.
    """
    order = ["Added", "Changed", "Deprecated", "Removed", "Fixed", "Security"]
    for version, groups in _subsections().items():
        known = [g for g in groups if g in order]
        assert known == sorted(known, key=order.index), (
            f"{version} has {known}, expected Keep a Changelog order")
