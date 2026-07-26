"""Solver output, and the provenance that makes it citable.

A bare efficiency number is not a result. Comparison plots across methods are
only credible if every point carries the method, the version, the truncation it
was computed at, and the evidence that it converged -- so :class:`Provenance` is
mandatory, not optional.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Any, Iterator

import numpy as np
from numpy.typing import NDArray

__all__ = ["Provenance", "OrderEfficiency", "EfficiencyScan"]


def _frozen_array(values: Any, dtype: Any) -> np.ndarray:
    """Return an independent, read-only copy of ``values``.

    Copying matters: freezing the caller's array in place would silently make
    *their* object read-only, a nasty surprise if they meant to reuse the
    buffer. Freezing matters because reference data that can be mutated after
    construction is reference data that cannot be trusted.
    """
    array = np.array(values, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class Provenance:
    """Where a number came from and whether it can be defended.

    Attributes
    ----------
    method
        Backend identifier, e.g. ``"rcwa"``, ``"scalar.harvey"``, ``"pcgrate:file"``.
    version
        Version of the code that produced it. For imported data, the version
        string of the foreign code.
    source
        Path or URL for imported data; ``None`` for computed results.
    truncation
        Fourier truncation order, boundary discretisation count, or whatever
        governs accuracy for the method. ``None`` where not applicable.
    converged
        ``True``/``False`` if a convergence study was run, ``None`` if not.
        **A result with ``converged=None`` has not been shown to be correct**
        and must not be presented as a rigorous value.
    wall_time_s
        Seconds of wall clock, for method-cost comparisons.
    warnings
        Validity-guard breaches (e.g. scalar theory used past its regime).
        Recorded rather than suppressed -- mapping where a model breaks down is
        a deliverable.
    notes
        Free-form extras: solver options, boundary conditions, anything needed
        to reproduce.
    """

    method: str
    version: str = "unknown"
    source: str | None = None
    truncation: int | None = None
    converged: bool | None = None
    wall_time_s: float | None = None
    warnings: tuple[str, ...] = ()
    notes: dict[str, Any] = field(default_factory=dict)

    @property
    def is_defensible(self) -> bool:
        """True only if convergence was actually demonstrated."""
        return self.converged is True

    def with_warning(self, message: str) -> "Provenance":
        """Return a copy carrying an additional warning."""
        return dataclasses.replace(self, warnings=self.warnings + (message,))

    def __hash__(self) -> int:
        """Hash over the identifying fields only, skipping ``notes``.

        ``notes`` is a free-form dict and therefore unhashable, but provenance
        needs to work as a dict key -- the comparison harness groups results by
        it. Two records differing only in ``notes`` collide on hash and are
        still distinguished by ``__eq__``, which is exactly what the hash
        contract permits.
        """
        return hash(
            (
                self.method,
                self.version,
                self.source,
                self.truncation,
                self.converged,
                self.warnings,
            )
        )


@dataclass(frozen=True, slots=True)
class OrderEfficiency:
    """Efficiencies for every order at a single wavelength."""

    wavelength: float
    orders: NDArray[np.int64]
    efficiency: NDArray[np.float64]
    propagating: NDArray[np.bool_]
    provenance: Provenance

    def __getitem__(self, order: int) -> float:
        """Efficiency in a given order. Evanescent and absent orders give 0.0."""
        hit = np.flatnonzero(self.orders == order)
        return float(self.efficiency[hit[0]]) if hit.size else 0.0

    @property
    def total(self) -> float:
        """Sum over propagating orders.

        Equals 1.0 for a lossless structure; equals the reflectivity of an
        equivalent surface for a blazed grating at grazing incidence.
        """
        return float(self.efficiency[self.propagating].sum())

    @property
    def propagating_orders(self) -> NDArray[np.int64]:
        return self.orders[self.propagating]

    def energy_balance_error(self) -> float:
        """``|total - 1|`` -- the lossless self-check. Meaningless if absorbing."""
        return abs(self.total - 1.0)


@dataclass(frozen=True, slots=True)
class EfficiencyScan:
    """Efficiency over a wavelength scan at fixed geometry.

    This is the shape both solvers and importers produce, and the unit the
    comparison harness works in.

    ``efficiency`` is ``(n_wavelengths, n_orders)``. Orders that do not
    propagate at a given wavelength carry ``0.0`` with ``propagating=False`` --
    never ``NaN``, and never dropped, so that column *j* means order
    ``orders[j]`` at every wavelength (``docs/conventions.md`` 4).
    """

    wavelengths: NDArray[np.float64]
    orders: NDArray[np.int64]
    efficiency: NDArray[np.float64]
    propagating: NDArray[np.bool_]
    provenance: Provenance

    def __post_init__(self) -> None:
        # Replace every field with an independent read-only copy, so a scan can
        # never be altered after construction. object.__setattr__ is how a
        # frozen dataclass legitimately initialises itself.
        for name, dtype in (
            ("wavelengths", np.float64),
            ("orders", np.int64),
            ("efficiency", np.float64),
            ("propagating", np.bool_),
        ):
            object.__setattr__(self, name, _frozen_array(getattr(self, name), dtype))

        n, m = len(self.wavelengths), len(self.orders)
        if self.efficiency.shape != (n, m):
            raise ValueError(
                f"efficiency has shape {self.efficiency.shape}, "
                f"expected ({n}, {m}) from {n} wavelengths and {m} orders"
            )
        if self.propagating.shape != (n, m):
            raise ValueError(
                f"propagating has shape {self.propagating.shape}, expected ({n}, {m})"
            )
        if np.isnan(self.efficiency).any():
            raise ValueError(
                "efficiency contains NaN; evanescent orders must be recorded as "
                "0.0 with propagating=False (docs/conventions.md 4)"
            )
        if len(np.unique(self.orders)) != m:
            raise ValueError("orders contains duplicates")

        # The invariant that keeps `total` trustworthy. Without it, a stray
        # value on a non-propagating order silently inflates the energy
        # balance -- the check the whole validation strategy rests on.
        stray = self.efficiency[~self.propagating]
        if stray.size and (stray != 0.0).any():
            bad = int(np.count_nonzero(stray))
            raise ValueError(
                f"{bad} non-propagating order(s) carry non-zero efficiency "
                f"(largest {np.abs(stray).max():.6g}). Evanescent orders must be "
                "exactly 0.0 (docs/conventions.md 4)"
            )
        if (self.efficiency < 0.0).any():
            raise ValueError("efficiency contains negative values")

    def __len__(self) -> int:
        return len(self.wavelengths)

    def __iter__(self) -> Iterator[OrderEfficiency]:
        for i in range(len(self)):
            yield self._row(i)

    def _row(self, i: int) -> OrderEfficiency:
        return OrderEfficiency(
            wavelength=float(self.wavelengths[i]),
            orders=self.orders,
            efficiency=self.efficiency[i],
            propagating=self.propagating[i],
            provenance=self.provenance,
        )

    def at(self, wavelength: float) -> OrderEfficiency:
        """Nearest sampled wavelength. Does not interpolate."""
        return self._row(int(np.argmin(np.abs(self.wavelengths - wavelength))))

    def order(self, m: int) -> NDArray[np.float64]:
        """Efficiency of order ``m`` across the scan."""
        hit = np.flatnonzero(self.orders == m)
        if not hit.size:
            raise KeyError(
                f"order {m} not in this scan; available: "
                f"{self.orders.min()}..{self.orders.max()}"
            )
        return self.efficiency[:, hit[0]]

    @property
    def total(self) -> NDArray[np.float64]:
        """Summed efficiency over propagating orders, per wavelength."""
        return np.where(self.propagating, self.efficiency, 0.0).sum(axis=1)

    def to_records(self) -> list[dict[str, Any]]:
        """Tidy long-form records: one row per (wavelength, order).

        The interchange shape for the comparison harness, and what a
        ``DataFrame`` gets built from without making pandas a hard dependency.
        """
        return [
            {
                "method": self.provenance.method,
                "wavelength": float(w),
                "order": int(m),
                "efficiency": float(self.efficiency[i, j]),
                "propagating": bool(self.propagating[i, j]),
                "converged": self.provenance.converged,
            }
            for i, w in enumerate(self.wavelengths)
            for j, m in enumerate(self.orders)
        ]
