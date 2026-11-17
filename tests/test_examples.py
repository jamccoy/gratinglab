"""Every script in ``examples/`` runs, and says what it claims to.

An example nobody runs goes stale the first time an interface moves, and the
person who finds out is the one who copied it. Each is run here in-process,
start to finish, and a line or two of its output is checked so that "it ran"
also means "it computed something".
"""

import runpy
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

#: Each numbered example, and text its output must contain.
EXPECTED = {
    "01_efficiency_scan.py": ("solved 41 wavelengths", "order +2: 0.7342", "first record:"),
    "02_compare_solvers.py": ("order 1 efficiency", "largest difference over all orders", "absorbed"),
    "03_your_own_material.py": ("Glass: covers 400 to 700 nm", "index at 500 nm: 1.5220", "refused", "typical:Au at 1 nm"),
    "04_mirror_and_coating.py": ("gold needed to read as bulk", "20 nm is enough", "20 nm is not enough", "r_s ="),
    "05_parameter_sweep.py": ("blaze angle   order +1 at 2 nm", "best of these: 15.0 deg", "typical:Ni"),
}


def test_every_numbered_example_is_checked():
    """A new example joins the table, or this says so."""
    numbered = sorted(p.name for p in EXAMPLES.glob("[0-9]*.py"))
    assert numbered == sorted(EXPECTED)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_example_runs_and_reports(name, capsys, monkeypatch):
    # 03 points GRATINGLAB_MATERIALS at a folder of its own; setting it here
    # first means monkeypatch puts back what the suite had.
    monkeypatch.setenv("GRATINGLAB_MATERIALS", "")
    runpy.run_path(str(EXAMPLES / name), run_name="__main__")
    out = capsys.readouterr().out
    for text in EXPECTED[name]:
        assert text in out, f"{name} did not print {text!r}"


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_examples_print_plain_ascii(name):
    """A degree sign in a print raises UnicodeEncodeError on a console whose
    code page lacks it, which is the first thing a new user would see."""
    text = (EXAMPLES / name).read_text(encoding="utf-8")
    assert all(ord(c) < 128 for c in text), name
