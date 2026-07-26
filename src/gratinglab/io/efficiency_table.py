"""Reader for efficiency tables exported by external grating codes.

An exported efficiency table is just data, so reference results from a
commercial solver can be compared against without needing a licence to run it.
This module reads the tabular dialect written by PCGrate-SX, which implements
the integral method and is the de-facto standard for grazing-incidence X-ray
gratings.

**The reader does not label the physics.** It cannot know from a file which
numerical method produced it, so :func:`read_scan` takes a ``method`` argument
that the caller supplies -- ``method="integral"`` for an integral-method code.
That keeps comparison legends describing physics ("scalar vs integral") rather
than products, while :attr:`Provenance.version` still records exactly which
program and version generated the reference data.

Dialects are distinguished by their banner line. Only one is implemented today;
:func:`read_table` is where a GSolver or RETICOLO dialect would be added, rather
than as a separate vendor module.

Format, confirmed by inspection of PCGrate-SX 6.7.1 output::

    line 1   "PCGrate-SX 6.7.1 (c)1996-2020 I.I.G., Inc."
    line 2   (blank)
    line 3   problem name
    line 4   Solved at <TAB> <date> <TAB> <time>
    line 5   Calculating time <TAB> hh:mm:ss
    line 6   (blank)
    line 7   <empty> <TAB> Scan Step <TAB> "Eff.TE(-7,R)" <TAB> ...
    line 8+  " Wav., nm" <TAB> 0.6 <TAB> 0.0196 <TAB> ...

Every line is tab-separated and right-padded with empty fields.

**Order indices come from the column headers, never from a positional offset.**
The prototype scripts in ``~/Documents/diffraction_efficiency`` read columns as
``data[i+11]``, which works only because ``np.genfromtxt`` splits the quoted
label ``" Wav., nm"`` into two whitespace-delimited fields. That offset is
different for every file and silently wrong if the order range changes.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..result import EfficiencyScan, Provenance

__all__ = ["EfficiencyTable", "read_table", "MISSING"]

#: PCGrate writes this for an order that does not propagate at a given scan point.
MISSING = "--"

_HEADER_LINES = 7

# "Eff.TE(-7,R)" -> polarization TE, order -7, reflected
_COLUMN_RE = re.compile(r'Eff\.(?P<pol>TE|TM)\((?P<order>-?\d+),(?P<rt>[RT])\)')

_VERSION_RE = re.compile(r"PCGrate-(?P<edition>\w+)\s+(?P<version>[\d.]+)")

# " Wav., nm"                        -> [("Wav.", "nm")]
# " P. ang., deg Az. ang., deg"      -> [("P. ang.", "deg"), ("Az. ang.", "deg")]
# PCGrate concatenates "<name>, <unit>" pairs for a multi-variable scan.
_SCAN_LABEL_RE = re.compile(r"([^,]+),\s*(\S+)")


@dataclass(frozen=True, slots=True)
class EfficiencyTable:
    """A parsed PCGrate export, before conversion to an :class:`EfficiencyScan`.

    ``scan_values`` is always 2-D, shape ``(n_points, n_scan_variables)``.
    PCGrate supports multi-variable scans -- an angle scan sweeping polar and
    azimuthal angle together writes both into one tab-delimited field -- so a
    1-D assumption would silently mis-read those files.
    """

    problem_name: str
    version: str
    solved_at: str
    calculating_time: str
    scan_variables: tuple[str, ...]
    scan_units: tuple[str, ...]
    scan_values: np.ndarray
    orders: np.ndarray
    polarization: str
    reflected: bool
    efficiency: np.ndarray
    propagating: np.ndarray
    source: Path

    @property
    def is_wavelength_scan(self) -> bool:
        """True for a single-variable scan over wavelength."""
        return len(self.scan_variables) == 1 and self.scan_variables[
            0
        ].lower().startswith("wav")

    @property
    def wavelengths(self) -> np.ndarray:
        """The wavelength axis. Raises if this is not a wavelength scan."""
        if not self.is_wavelength_scan:
            raise ValueError(
                f"{self.source}: scan variables are {self.scan_variables}, "
                "not a single wavelength axis"
            )
        return self.scan_values[:, 0]


def _split(line: str) -> list[str]:
    """Tab-split and strip surrounding quotes/whitespace from each field."""
    return [f.strip().strip('"').strip() for f in line.rstrip("\n").split("\t")]


def _parse_scan_label(label: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split a scan label into parallel (names, units) tuples.

    ``" Wav., nm"`` -> ``(("Wav.",), ("nm",))``
    ``" P. ang., deg Az. ang., deg"`` -> ``(("P. ang.", "Az. ang."), ("deg", "deg"))``
    """
    pairs = _SCAN_LABEL_RE.findall(label)
    if not pairs:
        return ((label.strip(),), ("",)) if label.strip() else ((), ())
    names = tuple(name.strip() for name, _ in pairs)
    units = tuple(unit.strip() for _, unit in pairs)
    return names, units


def read_table(path: str | Path, *, encoding: str = "cp1252") -> EfficiencyTable:
    """Parse a PCGrate-SX exported efficiency table.

    Parameters
    ----------
    path
        Path to the exported ``.txt``.
    encoding
        PCGrate runs on Windows and writes cp1252. Falls back to UTF-8 and
        latin-1 automatically, so this rarely needs setting.

    Raises
    ------
    ValueError
        If the file is not a PCGrate export, or if the header declares columns
        the data rows do not supply. Both indicate the file was edited or
        truncated, and silently guessing would corrupt reference data.
    """
    path = Path(path)
    raw = _read_text(path, encoding)
    lines = raw.splitlines()

    if len(lines) <= _HEADER_LINES:
        raise ValueError(f"{path}: too short to be a PCGrate export ({len(lines)} lines)")

    banner = _split(lines[0])[0]
    match = _VERSION_RE.search(banner)
    if match is None:
        raise ValueError(
            f"{path}: first line is not a PCGrate banner, got {banner!r}. "
            "Refusing to guess the layout."
        )
    version = f"PCGrate-{match['edition']} {match['version']}"

    problem_name = _split(lines[2])[0]
    solved_at = " ".join(f for f in _split(lines[3])[1:] if f)
    calculating_time = next((f for f in _split(lines[4])[1:] if f), "")

    orders, polarizations, rt_flags, columns = _parse_column_headers(
        _split(lines[_HEADER_LINES - 1]), path
    )

    scan_values, efficiency, propagating, scan_labels = _parse_rows(
        lines[_HEADER_LINES:], columns, path
    )

    scan_variables, scan_units = _parse_scan_label(scan_labels[0] if scan_labels else "")

    if scan_values.shape[1] != len(scan_variables):
        raise ValueError(
            f"{path}: scan label {scan_labels[0]!r} names {len(scan_variables)} "
            f"variable(s) but rows carry {scan_values.shape[1]} value(s)"
        )

    # PCGrate exports one polarization and one R/T sense per table.
    if len(set(polarizations)) != 1 or len(set(rt_flags)) != 1:
        raise ValueError(
            f"{path}: mixed polarizations {set(polarizations)} or senses "
            f"{set(rt_flags)} in one table; this reader expects one of each"
        )

    return EfficiencyTable(
        problem_name=problem_name,
        version=version,
        solved_at=solved_at,
        calculating_time=calculating_time,
        scan_variables=scan_variables,
        scan_units=scan_units,
        scan_values=scan_values,
        orders=orders,
        polarization=polarizations[0],
        reflected=rt_flags[0] == "R",
        efficiency=efficiency,
        propagating=propagating,
        source=path,
    )


def _read_text(path: Path, encoding: str) -> str:
    for enc in (encoding, "utf-8", "latin-1"):
        try:
            return path.read_text(encoding=enc)
        except UnicodeDecodeError:
            continue
    raise ValueError(f"{path}: could not decode as cp1252, utf-8 or latin-1")


def _parse_column_headers(
    fields: list[str], path: Path
) -> tuple[np.ndarray, list[str], list[str], list[int]]:
    """Map header fields to order indices. Returns (orders, pols, rt, col_indices)."""
    orders: list[int] = []
    polarizations: list[str] = []
    rt_flags: list[str] = []
    columns: list[int] = []

    for index, field in enumerate(fields):
        match = _COLUMN_RE.fullmatch(field)
        if match is None:
            continue
        orders.append(int(match["order"]))
        polarizations.append(match["pol"])
        rt_flags.append(match["rt"])
        columns.append(index)

    if not orders:
        raise ValueError(
            f"{path}: no 'Eff.TE(m,R)'-style columns found in the header row. "
            f"Got fields: {[f for f in fields if f][:6]}"
        )

    order_array = np.asarray(orders, dtype=np.int64)
    if len(np.unique(order_array)) != len(order_array):
        raise ValueError(f"{path}: duplicate order columns in header")

    return order_array, polarizations, rt_flags, columns


def _parse_rows(
    lines: list[str], columns: list[int], path: Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Read data rows into (scan_values, efficiency, propagating, labels)."""
    scan_values: list[list[float]] = []
    scan_labels: list[str] = []
    efficiency_rows: list[list[float]] = []
    propagating_rows: list[list[bool]] = []

    for line_number, line in enumerate(lines, start=_HEADER_LINES + 1):
        fields = _split(line)
        if not any(fields):
            continue  # trailing blank / tab-padded line

        if len(fields) <= max(columns):
            raise ValueError(
                f"{path}:{line_number}: row has {len(fields)} fields but the "
                f"header declares a column at index {max(columns)}. "
                "The file looks truncated."
            )

        # A multi-variable scan packs its values space-separated into one field.
        try:
            scan_values.append([float(v) for v in fields[1].split()])
        except ValueError as exc:
            raise ValueError(
                f"{path}:{line_number}: scan value {fields[1]!r} is not a number "
                "(or space-separated numbers)"
            ) from exc
        scan_labels.append(fields[0])

        if len(scan_values[-1]) != len(scan_values[0]):
            raise ValueError(
                f"{path}:{line_number}: row has {len(scan_values[-1])} scan "
                f"value(s), but the first row has {len(scan_values[0])}"
            )

        values: list[float] = []
        flags: list[bool] = []
        for column in columns:
            raw = fields[column]
            if raw in ("", MISSING):
                # Not propagating at this scan point. Recorded as 0.0, never NaN.
                values.append(0.0)
                flags.append(False)
            else:
                try:
                    values.append(float(raw))
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: efficiency {raw!r} is neither a "
                        f"number nor {MISSING!r}"
                    ) from exc
                flags.append(True)

        efficiency_rows.append(values)
        propagating_rows.append(flags)

    if not efficiency_rows:
        raise ValueError(f"{path}: header parsed but no data rows found")

    return (
        np.asarray(scan_values, dtype=np.float64),
        np.asarray(efficiency_rows, dtype=np.float64),
        np.asarray(propagating_rows, dtype=np.bool_),
        scan_labels,
    )


def to_scan(table: EfficiencyTable, *, method: str = "imported") -> EfficiencyScan:
    """Convert a parsed table to the common :class:`EfficiencyScan` container.

    Only wavelength scans convert: an angle scan is a different object and
    silently relabelling its axis would corrupt a comparison.

    Parameters
    ----------
    method
        Label for this data in comparisons -- ``"integral"`` for output from an
        integral-method code. **The caller supplies it because the file does
        not say.** A table records efficiencies, not which numerical method
        produced them, so guessing from the vendor banner would be an
        inference dressed up as a fact. The exact program and version stay on
        ``Provenance.version`` regardless.
    """
    if not table.is_wavelength_scan:
        raise ValueError(
            f"{table.source}: scan variables are {table.scan_variables}, not a "
            "single wavelength axis. EfficiencyScan is wavelength-indexed; an "
            "angle scan is a different object and relabelling its axis would "
            "corrupt any comparison built on it."
        )
    unit = table.scan_units[0]
    if unit and unit.lower() != "nm":
        raise ValueError(
            f"{table.source}: wavelength unit is {unit!r}; "
            "gratinglab works in nm (docs/conventions.md 1)"
        )

    provenance = Provenance(
        method=method,
        version=table.version,
        source=str(table.source),
        # The exporting program's own truncation is not in the table; it lives
        # in its project file. Left None rather than guessed.
        truncation=None,
        converged=None,
        notes={
            "vendor": table.version.split()[0] if table.version else "unknown",
            "problem_name": table.problem_name,
            "solved_at": table.solved_at,
            "calculating_time": table.calculating_time,
            "polarization": table.polarization,
            "sense": "R" if table.reflected else "T",
        },
    )

    return EfficiencyScan(
        wavelengths=table.wavelengths,
        orders=table.orders,
        efficiency=table.efficiency,
        propagating=table.propagating,
        provenance=provenance,
    )


def read_scan(
    path: str | Path, *, method: str = "imported", **kwargs
) -> EfficiencyScan:
    """Read an exported efficiency table straight into an :class:`EfficiencyScan`.

    >>> read_scan("run.txt", method="integral")   # doctest: +SKIP

    See :func:`to_scan` for why ``method`` is the caller's to declare.
    """
    return to_scan(read_table(path, **kwargs), method=method)
