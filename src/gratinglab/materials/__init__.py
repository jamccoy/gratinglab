r"""Optical constants: looking them up, and reading your own.

```python
from gratinglab import materials

materials.available()          # every material in every library
pt = materials.lookup("Pt")    # ~/.gratinglab/materials/Pt.txt, say
mine = materials.from_cxro_file("Ni_CXRO.txt")   # any file, no library
film = materials.from_nk_file("Ni_film.csv")     # wavelength, n, k columns
au = materials.typical("Au")   # a model, not a table: see typical.py
sio2 = materials.lookup("factors:SiO2@2.2")   # from Si.nff and O.nff
film = materials.lookup("Au@17.5")            # Au.txt rescaled to 17.5 g/cm^3
```

A library is a directory of tables, one material per file, named for it:
`Pt.txt` (a CXRO export, X-ray), `Pt.yml` (a refractiveindex.info file,
ultraviolet to infrared), `Pt.ari` (PCGrate) or `Pt.nk` / `Pt.csv` (plain
wavelength, n, k columns, see :mod:`.nk`) is `coating="Pt"`. Yours is
`~/.gratinglab/materials`, or whatever `GRATINGLAB_MATERIALS` names, and it is
searched before the tables bundled with a checkout, so your own replaces a
bundled one of the same name. Neither the repository nor a release bundles
any, because no published table is clearly public domain; `data/SOURCES.md`
says so.

A name beginning `typical:` (`coating="typical:Au"`) is not looked up in any
library: it is the X-ray model in :mod:`.typical`, built from physical facts,
approximate, and named so that a result says it used a model. A plain `"Au"`
with no table still refuses rather than falling back to it.

A name beginning `factors:` (`coating="factors:SiO2@2.2"`) is computed from a
formula, a density and one tabulated scattering-factor file per element,
`Si.nff` and `O.nff`, found in the libraries like any table: see
:mod:`.factors`. The density after `@` is in g/cm^3 and may be left off where
:data:`.typical.DENSITIES` has one.

A plain name with a density, `coating="Au@17.5"`, is the `Au` table rescaled
to that density, since the decrement and absorption are both proportional to
it. Only a CXRO export states the density it was made at, so only those can.

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
import re
import shutil
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

from . import factors as _factors
from .factors import ScatteringFactors, read_nff
from .optical import (
    HC_EV_NM,
    OpticalConstants,
    cxro_density,
    read_ari,
    read_cxro,
    write_ari,
)
from .nk import read_nk
from .refractiveindex import read_refractiveindex
from .typical import DENSITIES
from .typical import PREFIX as TYPICAL_PREFIX
from .typical import TypicalConstants, from_name, typical

__all__ = [
    "HC_EV_NM",
    "LIBRARY_ENV",
    "OpticalConstants",
    "ScatteringFactors",
    "TypicalConstants",
    "UnknownMaterial",
    "available",
    "data_dir",
    "factors",
    "from_ari_file",
    "from_cxro_file",
    "from_nk_file",
    "from_refractiveindex_file",
    "install",
    "library_dirs",
    "lookup",
    "read_ari",
    "read_cxro",
    "read_nff",
    "read_nk",
    "read_refractiveindex",
    "typical",
    "user_dir",
    "write_ari",
]

#: Bundled tables would live here, and there are none, in a checkout or a
#: release; see ``data/SOURCES.md``. The test suite points this elsewhere.
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
    ".nk": read_nk,
    ".csv": read_nk,
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
    """A material by name, e.g. ``"Au"``, from the first library that has it.

    ``typical:<formula>[@<density>]`` is the model in :mod:`.typical` instead,
    ``factors:<formula>[@<density>]`` is computed from tabulated scattering
    factors (:func:`factors`), and ``<name>@<density>`` is a CXRO table
    rescaled to that density.
    """
    if name.startswith(TYPICAL_PREFIX):
        return from_name(name)
    if name.startswith(_factors.PREFIX):
        return factors(*_split_density(name[len(_factors.PREFIX):], name))
    if "@" in name:
        return _rescaled(*_split_density(name, name), name)
    for directory in library_dirs():
        path = _tables(directory).get(name)
        if path is not None:
            return _read(path, path.stat().st_mtime_ns, name)
    known = ", ".join(available()) or "none installed"
    searched = "; ".join(str(d) for d in library_dirs())
    raise UnknownMaterial(
        f"no optical constants for {name!r}; available: {known}. Searched: "
        f"{searched}. Add a CXRO export as {name}.txt, a refractiveindex.info "
        f"file as {name}.yml, a PCGrate {name}.ari, or wavelength, n, k columns as "
        f"{name}.nk, to {user_dir()} or to a "
        f"directory named in {LIBRARY_ENV}; or read one directly with "
        "materials.from_cxro_file(path)."
    )


def _split_density(body: str, name: str) -> "tuple[str, float | None]":
    """``"SiO2@2.2"`` as ``("SiO2", 2.2)``, and ``"Au"`` as ``("Au", None)``."""
    base, _, rho = body.partition("@")
    if not rho:
        return base, None
    try:
        return base, float(rho)
    except ValueError:
        raise ValueError(f"{name!r}: the density after @ must be a number in g/cm^3") from None


def factors(formula: str, density: float | None = None) -> OpticalConstants:
    """``formula`` at ``density`` g/cm^3 from the ``<Symbol>.nff`` files in
    your libraries; see :mod:`.factors`. ``density`` defaults as
    :func:`typical`'s does, and is required for anything else."""
    if density is None:
        if formula not in DENSITIES:
            raise ValueError(
                f"no default density for {formula!r}; write "
                f"{_factors.PREFIX}{formula}@<g/cm^3> (defaults exist for "
                f"{', '.join(DENSITIES)})"
            )
        density = DENSITIES[formula]
    tables = {}
    missing = []
    for symbol, _ in _factors.parse_formula(formula):
        path = _factor_file(symbol)
        if path is None:
            missing.append(symbol)
        else:
            tables[symbol] = _read_nff(path, path.stat().st_mtime_ns, symbol)
    if missing:
        searched = "; ".join(str(d) for d in library_dirs())
        raise UnknownMaterial(
            f"no scattering factors for {', '.join(missing)} in {formula!r}. "
            "Download each element's table from "
            "https://henke.lbl.gov/optical_constants/asf.html and save it as "
            f"{missing[0]}{_factors.SUFFIX} in {user_dir()} or a directory "
            f"named in {LIBRARY_ENV}. Searched: {searched}."
        )
    suffix = "" if density == DENSITIES.get(formula) else f"@{density:g}"
    name = f"{_factors.PREFIX}{formula}{suffix}"
    return _factors.factors(formula, density, tables, name=name)


def _factor_file(symbol: str) -> "Path | None":
    """The first library's ``<symbol>.nff``, matched without case, since CXRO
    names its downloads ``au.nff``."""
    want = f"{symbol}{_factors.SUFFIX}".lower()
    for directory in library_dirs():
        if directory.is_dir():
            for path in sorted(directory.iterdir()):
                if path.name.lower() == want and path.is_file():
                    return path
    return None


@lru_cache(maxsize=None)
def _read_nff(path: Path, mtime_ns: int, symbol: str) -> ScatteringFactors:
    """Cached by file and modification time, as :func:`_read` is."""
    table = read_nff(path, symbol)
    return replace(table, source=f"{table.source} in {path.parent}")


def _rescaled(base: str, density: "float | None", name: str) -> OpticalConstants:
    """The CXRO table ``base`` at ``density`` instead of the one it states."""
    if density is None or not density > 0:
        raise ValueError(f"{name!r}: the density after @ must be positive")
    table = lookup(base)  # raises the usual UnknownMaterial
    path = next(t[base] for d in library_dirs() if base in (t := _tables(d)))
    stated = cxro_density(path) if path.suffix.lower() == ".txt" else None
    if stated is None:
        raise ValueError(
            f"{name!r}: {path.name} does not state the density it was made "
            "at, so it cannot be rescaled. Only a CXRO export's header "
            "(Au Density=19.32) does"
        )
    k = density / stated
    return replace(
        table,
        name=name,
        decrement=table.decrement * k,
        absorption=table.absorption * k,
        source=f"{table.source}, rescaled from {stated:g} to {density:g} g/cm^3",
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


#: A material name is a file stem in a library and a ``coating`` string on a
#: ``Problem``, so it has to be safe as both: no separators, no leading dot.
_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._+-]*")


def install(
    path: "str | Path", name: str | None = None, *, overwrite: bool = False
) -> Path:
    """Copy a table into your library, after proving that it reads.

    The table goes into the first directory :func:`library_dirs` searches,
    which is :func:`user_dir` unless ``GRATINGLAB_MATERIALS`` names others,
    as ``<name><suffix>``; ``name`` defaults to the file's stem. A
    refractiveindex.info file is named for its author (``Johnson.yml``), so
    that is the case for passing one.

    The file is read with the reader for its suffix first, and whatever that
    raises is raised here, so an unreadable table never reaches the library
    to fail later inside a solve. A table already there under that name is
    refused unless ``overwrite`` is set (``FileExistsError``), and one under
    the same name with a different suffix is refused regardless
    (``ValueError``), because :func:`lookup` would
    then refuse the pair and replacing it means deleting a file of yours.
    Returns the path written.
    """
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix == _factors.SUFFIX:
        return _install_factors(source, name, overwrite)
    if suffix not in _READERS:
        raise ValueError(
            f"{source.name} is not a table this library reads; use one of "
            f"{', '.join(sorted([*_READERS, _factors.SUFFIX]))}"
        )
    name = source.stem if name is None else name
    if not _NAME.fullmatch(name):
        raise ValueError(
            f"{name!r} cannot be a material name: use letters, digits and "
            "._+- only, starting with a letter or digit (e.g. Pt, SiO2, B4C)"
        )
    _READERS[suffix](source, name=name)

    directory = library_dirs()[0]
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / f"{name}{suffix}"
    existing = _tables(directory).get(name)
    if existing is not None:
        if existing.suffix.lower() != suffix:
            raise ValueError(
                f"{directory} already holds {existing.name} for {name!r}; "
                f"remove it to install {source.name} in its place, since a "
                "library cannot hold two tables for one name"
            )
        if existing.resolve() == source.resolve():
            return existing
        if not overwrite:
            raise FileExistsError(f"{existing} already exists")
        destination = existing
    shutil.copyfile(source, destination)
    return destination


def _install_factors(source: Path, symbol: "str | None", overwrite: bool) -> Path:
    """Copy one element's scattering factors in as ``<Symbol>.nff``."""
    table = read_nff(source, symbol.capitalize() if symbol else None)
    directory = library_dirs()[0]
    directory.mkdir(parents=True, exist_ok=True)
    existing = _factor_file(table.symbol)
    destination = directory / f"{table.symbol}{_factors.SUFFIX}"
    if existing is not None and existing.parent == directory:
        if existing.resolve() == source.resolve():
            return existing
        if not overwrite:
            raise FileExistsError(f"{existing} already exists")
        destination = existing
    shutil.copyfile(source, destination)
    return destination


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


def from_nk_file(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a plain ``wavelength, n[, k]`` table, such as an ellipsometry
    export; the header must name the first column's unit. See :mod:`.nk`."""
    return read_nk(path, name=name)


def from_ari_file(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a PCGrate ``.ari`` table.

    Note the columns are ``(wavelength_nm, decrement, absorption)`` despite the
    ``_n_k`` in the corpus filenames; see :func:`~.optical.read_ari`.
    """
    return read_ari(path, name=name)
