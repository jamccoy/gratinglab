r"""Optical constants: looking them up, and reading your own.

```python
from gratinglab import materials

materials.available()          # every material in every library
pt = materials.lookup("Pt")    # ~/.gratinglab/materials/Pt.txt, say
mine = materials.from_cxro_file("Ni_CXRO.txt")   # any file, no library
```

A library is a directory of tables, one material per file, named for it:
`Pt.txt` (a CXRO export, X-ray), `Pt.yml` (a refractiveindex.info file,
ultraviolet to infrared) or `Pt.ari` (PCGrate) is `coating="Pt"`. Yours is
`~/.gratinglab/materials`, or whatever `GRATINGLAB_MATERIALS` names, and it is
searched before the tables bundled with a checkout, so your own replaces a
bundled one of the same name. The release packages bundle none, because CXRO's
terms for redistributing its tables are not settled; `data/SOURCES.md` says so.

Lookup is by name, because that is what a
:class:`~gratinglab.problem.Problem` carries: `coating="Au"` keeps a benchmark
case a small serialisable file instead of embedding a 500-row table in it, and
keeps `Problem` free of any dependency on which tables happen to be installed.

An unknown name raises an error that lists what is available. This is more
important than it may appear. `Problem.coating` was formerly a free-form string
that nothing read, so setting it to any value silently relabelled a result as
"absolute" while changing no number. Resolving the name turns that string into
a claim that the code has to honour.
"""

from __future__ import annotations

import os
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

from .optical import (
    HC_EV_NM,
    OpticalConstants,
    read_ari,
    read_cxro,
    write_ari,
)
from .refractiveindex import read_refractiveindex

__all__ = [
    "HC_EV_NM",
    "LIBRARY_ENV",
    "OpticalConstants",
    "UnknownMaterial",
    "available",
    "data_dir",
    "from_ari_file",
    "from_cxro_file",
    "from_refractiveindex_file",
    "library_dirs",
    "lookup",
    "read_ari",
    "read_cxro",
    "read_refractiveindex",
    "user_dir",
    "write_ari",
]

#: Tables bundled with a checkout live here. The release packages ship none;
#: see ``data/SOURCES.md``.
_DATA = Path(__file__).parent / "data"

#: Names your own library directories, separated by ``os.pathsep``. When set it
#: replaces :func:`user_dir` rather than adding to it, the same rule
#: ``GRATINGLAB_AFM_DIR`` follows.
LIBRARY_ENV = "GRATINGLAB_MATERIALS"

#: One reader per file suffix. A library is a directory of these, one material
#: per file, named for it: ``Pt.txt`` is ``coating="Pt"``.
_READERS = {
    ".txt": read_cxro,
    ".ari": read_ari,
    ".yml": read_refractiveindex,
    ".yaml": read_refractiveindex,
}


class UnknownMaterial(KeyError):
    """A material name with no table behind it.

    It is a ``KeyError`` because a failed lookup is one, and it has its own
    type so that a caller can distinguish "you asked for a material I do not
    have" from any other missing key. The GUI reports the first case
    differently.
    """

    def __str__(self) -> str:  # pragma: no cover - KeyError quotes its arg
        return self.args[0] if self.args else ""


def data_dir() -> Path:
    """Where bundled tables are read from. Exposed so a test can say that the
    directory being empty is a supported state rather than a broken install."""
    return _DATA


def user_dir() -> Path:
    """Your default library, used when ``GRATINGLAB_MATERIALS`` is unset.

    It sits under the home directory and not beside the package, because an
    install's own directory is site-packages: a table placed there is lost on
    the next upgrade, and the directory may not be writable at all.
    """
    return Path.home() / ".gratinglab" / "materials"


def library_dirs() -> tuple[Path, ...]:
    """Every directory searched for a material, in the order searched.

    Your own directories come first, then the bundled one, so a table of your
    own replaces a bundled one of the same name. The table's ``source`` names
    the file it was read from, so a result's provenance says which one was
    used. Directories that do not exist are listed anyway, since they are
    where a lookup looked.
    """
    configured = os.environ.get(LIBRARY_ENV)
    if configured:
        yours = [Path(p).expanduser() for p in configured.split(os.pathsep) if p]
    else:
        yours = [user_dir()]
    return (*yours, _DATA)


def _tables(directory: Path) -> dict[str, Path]:
    """The tables in one library directory, by material name.

    Two files for one name in the same directory (``Pt.txt`` and ``Pt.ari``)
    are refused instead of one being picked, because neither suffix outranks
    the other and any choice would be a guess about which the user meant.
    """
    if not directory.is_dir():
        return {}
    found: dict[str, Path] = {}
    for path in sorted(directory.iterdir()):
        if path.suffix.lower() not in _READERS or not path.is_file():
            continue
        if path.stem in found:
            raise ValueError(
                f"{directory} holds two tables for {path.stem!r}, "
                f"{found[path.stem].name} and {path.name}; keep one"
            )
        found[path.stem] = path
    return found


def available() -> tuple[str, ...]:
    """Material names across every library, sorted.

    The directories are read in place of a hardcoded list, so adding a
    material only requires adding a file. An empty result is a supported
    configuration and not an error, since it is the state an install starts
    in.
    """
    names: set[str] = set()
    for directory in library_dirs():
        names.update(_tables(directory))
    return tuple(sorted(names))


def lookup(name: str) -> OpticalConstants:
    """A material by name, e.g. ``"Au"``, from the first library that has it."""
    for directory in library_dirs():
        path = _tables(directory).get(name)
        if path is not None:
            return _read(path, path.stat().st_mtime_ns, name)
    known = ", ".join(available()) or "none installed"
    searched = "; ".join(str(d) for d in library_dirs())
    raise UnknownMaterial(
        f"no optical constants for {name!r}; available: {known}. Searched: "
        f"{searched}. Add a CXRO export as {name}.txt, a refractiveindex.info "
        f"file as {name}.yml, or a PCGrate {name}.ari, to {user_dir()} or to a "
        f"directory named in {LIBRARY_ENV}; or read one directly with "
        "materials.from_cxro_file(path)."
    )


@lru_cache(maxsize=None)
def _read(path: Path, mtime_ns: int, name: str) -> OpticalConstants:
    """Cached by file and modification time instead of by name.

    The tables are immutable (`OpticalConstants` freezes its arrays) and a
    solve may resolve the same coating once per wavelength grid, so parsing
    500 rows each time would be wasteful. Keying on the name alone would keep serving a
    table after its file was edited or a library was added in front of it.
    """
    table = _READERS[path.suffix.lower()](path, name=name)
    # The reader records the file name; the directory is what tells two
    # libraries' Au.txt apart, and `source` is what provenance carries.
    return replace(table, source=f"{table.source}, in {path.parent}")


def from_cxro_file(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a CXRO/Henke ``(energy, delta, beta)`` export you downloaded."""
    return read_cxro(path, name=name)


def from_refractiveindex_file(
    path: "str | Path", name: str | None = None
) -> OpticalConstants:
    """Read a refractiveindex.info YAML file, ultraviolet to infrared.

    Needs the ``refractiveindex`` extra. See :mod:`.refractiveindex` for what
    each kind of entry becomes, and which are refused.
    """
    return read_refractiveindex(path, name=name)


def from_ari_file(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a PCGrate ``.ari`` table.

    Note the columns are ``(wavelength_nm, decrement, absorption)`` despite the
    ``_n_k`` in the corpus filenames; see :func:`~.optical.read_ari`.
    """
    return read_ari(path, name=name)
