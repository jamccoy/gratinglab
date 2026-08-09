r"""Optical constants: looking them up, and reading your own.

```python
from gratinglab import materials

materials.available()          # ('Au',)
gold = materials.lookup("Au")  # a vendored table
mine = materials.from_cxro_file("Ni_CXRO.txt")   # your own download
```

**A vendored table is a convenience, never a requirement.** Everything here
works on a user-supplied CXRO export, and `data/SOURCES.md` records what is
bundled and why that directory could be deleted outright without changing any
code. See that file before the first public release.

Lookup is by name, deliberately, because that is what a
:class:`~gratinglab.problem.Problem` carries: `coating="Au"` keeps a benchmark
case a small serialisable file instead of embedding a 500-row table in it, and
keeps `Problem` free of any dependency on which tables happen to be installed.

An unknown name **raises and lists what is available**. That matters more here
than it looks: `Problem.coating` was a free-form string that nothing read, so
setting it to anything at all silently relabelled a result as "absolute" while
changing no number. Resolution is what turns that string into a claim something
has to honour.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from .optical import (
    HC_EV_NM,
    OpticalConstants,
    read_ari,
    read_cxro,
    write_ari,
)

__all__ = [
    "HC_EV_NM",
    "OpticalConstants",
    "UnknownMaterial",
    "available",
    "data_dir",
    "from_ari_file",
    "from_cxro_file",
    "lookup",
    "read_ari",
    "read_cxro",
    "write_ari",
]

#: Vendored tables live here. One CXRO export per material, named for it.
_DATA = Path(__file__).parent / "data"


class UnknownMaterial(KeyError):
    """A material name with no table behind it.

    A ``KeyError`` because that is what a failed lookup is, but its own type so
    a caller can distinguish "you asked for a material I do not have" from any
    other missing key -- the GUI wants to say the first differently.
    """

    def __str__(self) -> str:  # pragma: no cover - KeyError quotes its arg
        return self.args[0] if self.args else ""


def data_dir() -> Path:
    """Where vendored tables are read from. Exposed so a test can say that the
    directory being empty is a supported state rather than a broken install."""
    return _DATA


def available() -> tuple[str, ...]:
    """Vendored material names, sorted.

    Reads the directory rather than a hardcoded list, so adding a table is a
    file drop -- see ``data/SOURCES.md``. Returns empty if the directory has
    been removed, which is a supported configuration, not an error.
    """
    if not _DATA.is_dir():
        return ()
    return tuple(sorted(p.stem for p in _DATA.glob("*.txt")))


@lru_cache(maxsize=None)
def lookup(name: str) -> OpticalConstants:
    """A vendored table by material name, e.g. ``"Au"``.

    Cached: the tables are immutable (`OpticalConstants` freezes its arrays),
    a solve may resolve the same coating once per wavelength grid, and parsing
    500 rows repeatedly is pure waste.
    """
    path = _DATA / f"{name}.txt"
    if not path.is_file():
        known = ", ".join(available()) or "none installed"
        raise UnknownMaterial(
            f"no optical constants for {name!r}; available: {known}. "
            "Supply your own with materials.from_cxro_file(path), or see "
            "src/gratinglab/materials/data/SOURCES.md to add a table."
        )
    return read_cxro(path, name=name)


def from_cxro_file(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a CXRO/Henke ``(energy, delta, beta)`` export you downloaded."""
    return read_cxro(path, name=name)


def from_ari_file(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a PCGrate ``.ari`` table.

    Note the columns are ``(wavelength_nm, decrement, absorption)`` despite the
    ``_n_k`` in the corpus filenames -- see :func:`~.optical.read_ari`.
    """
    return read_ari(path, name=name)
