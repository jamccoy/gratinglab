"""Saving a measured groove profile with its scale and where it came from.

A ``.ggp`` holds a groove's shape and nothing else: it is normalised to one
period and has no room for the period itself, so a profile written to one and
read back has lost the number the solver's phase term scales by. That is the
gap ``docs/roadmap.md`` item 5 names. Two formats close it here:

``.json``, gratinglab's own
    The normalised profile at full precision, the measured period in nm, the
    groove count, the profile metrics, and provenance: the gratinglab version,
    the scan it came from, how many grooves were detected and averaged, and
    the analysis settings. :func:`read_boundary_json` reads it back, and its
    :meth:`BoundaryFile.to_problem` gives the same :class:`~gratinglab.problem.Problem`
    the in-process ``BoundaryProfile.to_problem`` does, so saving loses nothing.
    The main window's "Load boundary profile" reads it and fills in the period.

``.csv``, for everything else
    One row per point, ``t`` and ``y`` as fractions of the period plus the same
    points in nm, under ``#`` comment lines carrying the period, groove counts
    and source. It is intended for a spreadsheet, MATLAB, or another code. It
    is write-only; the JSON is the format to read back.

Both writers take any object with the fields of
``gratinglab.metrology.boundary.BoundaryProfile`` (``x_norm``, ``y_norm``,
``period_nm``, ``n_grooves``, ``n_used``, ``n_edge_rejected``, ``metrics``)
rather than importing it, because the metrology package needs the plotting
extra and this module must not.
"""

from __future__ import annotations

import dataclasses
import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "FORMAT",
    "VERSION",
    "BoundaryFile",
    "read_boundary_json",
    "read_profile_file",
    "write_boundary_csv",
    "write_boundary_json",
]

#: Identifies the file, so a reader can refuse an unrelated JSON by name
#: rather than failing on a missing key three lines later.
FORMAT = "gratinglab-boundary"

#: Bumped when a reader of the previous version could misread a new file.
#: Adding an optional key is not such a change.
VERSION = 1


@dataclass(frozen=True)
class BoundaryFile:
    """A groove profile read back from a gratinglab ``.json``."""

    period_nm: float
    #: Fractions of one period along ``+d-hat``, as ``x_norm`` and a ``.ggp``
    #: hold them; :meth:`profile` reflects them into the profile parameter.
    t: tuple[float, ...]
    y: tuple[float, ...]
    #: Grooves averaged into this profile, carried onto the problem.
    n_grooves: int | None = None
    #: A fitted facet angle in degrees, if the writer had one.
    blaze_angle: float | None = None
    metrics: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)

    def profile(self):
        """The shape, as the solvers take it."""
        from ..profiles import FromProfileData

        return FromProfileData.along_dispersion(
            self.t, self.y, blaze_angle=self.blaze_angle
        )

    def to_problem(self, *, coating=None, substrate=None, roughness=0.0):
        """This groove as an efficiency problem, with the period it was measured at.

        This is the same problem ``BoundaryProfile.to_problem`` builds in
        process, which is the purpose of the format; see
        ``tests/test_io_boundary.py``.
        Roughness is not applied, for the reason that method's docstring
        gives. The measured groove residual travels in
        ``metrics['groove_residual_rms_nm']`` and is passed only on request,
        as ``roughness=``.
        """
        from ..problem import Problem

        return Problem(
            period=self.period_nm,
            profile=self.profile(),
            coating=coating,
            substrate=substrate,
            roughness=roughness,
            n_grooves=self.n_grooves,
        )


def write_boundary_json(
    path: "str | Path",
    profile: Any,
    *,
    source: "str | Path | None" = None,
    settings: Any = None,
    blaze_angle: float | None = None,
) -> Path:
    """Write ``profile`` (a ``BoundaryProfile``) as a gratinglab ``.json``.

    ``source`` names the scan it came from and ``settings`` is the
    ``AnalysisSettings`` it was built with; both land in the provenance and
    neither is required. Floats are written at full precision, so reading the
    file back reproduces the profile exactly, where a ``.ggp`` keeps six
    decimals.
    """
    path = Path(path)
    t, y = _points(profile)
    document = {
        "format": FORMAT,
        "version": VERSION,
        "period_nm": float(profile.period_nm),
        "profile": {
            "units": "t and y are fractions of one period",
            "t": t,
            "y": y,
        },
        "n_grooves": int(profile.n_used),
        "blaze_angle_deg": None if blaze_angle is None else float(blaze_angle),
        "metrics": _plain(profile.metrics),
        "provenance": _provenance(profile, source, settings),
    }
    text = json.dumps(document, indent=1, allow_nan=False) + "\n"
    path.write_bytes(text.encode("utf-8"))
    return path


def read_boundary_json(path: "str | Path") -> BoundaryFile:
    """Read a gratinglab ``.json`` boundary profile.

    Refuses, by name, a JSON that is not one, a newer version than this reader
    knows, and a profile whose two columns differ in length or whose period is
    not a positive number. The profile's own checks (at least two points, and
    so on) run when :meth:`BoundaryFile.profile` builds it.
    """
    path = Path(path)
    try:
        document = json.loads(path.read_bytes().decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path}: not a JSON file ({exc})") from exc
    if not isinstance(document, dict) or document.get("format") != FORMAT:
        raise ValueError(f"{path}: not a gratinglab boundary profile (no format: {FORMAT!r})")
    version = document.get("version")
    if not isinstance(version, int) or version > VERSION:
        raise ValueError(
            f"{path}: boundary format version {version!r}; this gratinglab reads "
            f"up to {VERSION}, so upgrade it to read this file")
    period = document.get("period_nm")
    if not isinstance(period, (int, float)) or not math.isfinite(period) or period <= 0:
        raise ValueError(f"{path}: period_nm must be a positive number, got {period!r}")
    points = document.get("profile") or {}
    t, y = points.get("t"), points.get("y")
    if not isinstance(t, list) or not isinstance(y, list) or len(t) != len(y):
        raise ValueError(f"{path}: profile t and y must be two lists of equal length")
    blaze = document.get("blaze_angle_deg")
    n_grooves = document.get("n_grooves")
    return BoundaryFile(
        period_nm=float(period),
        t=tuple(float(v) for v in t),
        y=tuple(float(v) for v in y),
        n_grooves=None if n_grooves is None else int(n_grooves),
        blaze_angle=None if blaze is None else float(blaze),
        metrics=document.get("metrics") or {},
        provenance=document.get("provenance") or {},
    )


def write_boundary_csv(
    path: "str | Path",
    profile: Any,
    *,
    source: "str | Path | None" = None,
    settings: Any = None,
) -> Path:
    """Write ``profile`` as a CSV table under ``#`` comment lines.

    Columns ``t, y`` are fractions of one period and ``x_nm, y_nm`` the same
    points scaled by the measured period, so a reader needs no other number to
    plot it at true size. ``settings`` is accepted for symmetry with
    :func:`write_boundary_json`; only the JSON records them in full.
    """
    path = Path(path)
    t, y = _points(profile)
    period = float(profile.period_nm)
    prov = _provenance(profile, source, None)
    # The column header comes FIRST and the comments after it: numpy's
    # genfromtxt(names=True) takes names from line one even when it is a
    # comment, so a comment there containing a comma became two junk column
    # names. Header-first reads cleanly in numpy, pandas (comment="#") and a
    # spreadsheet, whose first row is then the column names.
    lines = [
        "t,y,x_nm,y_nm",
        f"# gratinglab boundary profile written by gratinglab {prov['gratinglab_version']}",
        f"# period_nm: {period!r}",
        f"# grooves: {prov['grooves_averaged']} averaged of {prov['grooves_detected']} detected",
    ]
    if prov.get("source"):
        lines.append(f"# source: {prov['source']}")
    lines.append(
        "# t and y are fractions of one period; x_nm and y_nm are the same points times the period")
    lines += [f"{ti!r},{yi!r},{ti * period!r},{yi * period!r}" for ti, yi in zip(t, y)]
    path.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
    return path


def read_profile_file(path: "str | Path"):
    """A profile from either a gratinglab ``.json`` or a PCGrate ``.ggp``.

    The format is chosen by suffix, because the two are unrelated formats and
    sniffing the content would only turn a clear refusal into a confusing one.
    Only the JSON carries a period; read it with :func:`read_boundary_json`
    to get that too.
    """
    path = Path(path)
    if path.suffix.lower() == ".json":
        return read_boundary_json(path).profile()
    from .ggp import read_profile

    return read_profile(path)


def _points(profile: Any) -> tuple[list[float], list[float]]:
    t = [float(v) for v in np.asarray(profile.x_norm, dtype=float)]
    y = [float(v) for v in np.asarray(profile.y_norm, dtype=float)]
    if len(t) != len(y):
        raise ValueError("profile x_norm and y_norm differ in length")
    return t, y


def _provenance(profile: Any, source: "str | Path | None", settings: Any) -> dict:
    from .. import __version__

    record = {
        "gratinglab_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": None if source is None else str(source),
        "grooves_detected": int(profile.n_grooves),
        "grooves_averaged": int(profile.n_used),
        "grooves_edge_rejected": int(profile.n_edge_rejected),
    }
    if settings is not None:
        record["settings"] = _plain(
            dataclasses.asdict(settings) if dataclasses.is_dataclass(settings) else dict(settings)
        )
    return record


def _plain(value: Any) -> Any:
    """JSON-safe: numpy scalars and arrays to Python, tuples to lists, NaN refused."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value
