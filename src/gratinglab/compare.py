"""Run one problem through several methods and line the answers up.

This is the point of the project. A comparison is only honest if every method
saw the same physics, so :func:`sweep` takes one :class:`~gratinglab.problem.Problem`
and one :class:`~gratinglab.illumination.Illumination` and hands them to each
backend unchanged.

Precomputed results are first-class participants. An imported PCGrate table is
passed in exactly like a live solver, so scalar-versus-integral-method plots are
available without a PCGrate licence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .illumination import Illumination
from .problem import Problem
from .result import EfficiencyScan
from .solvers.base import get_solver

__all__ = ["sweep", "align", "Comparison", "records"]

#: A method is either a registered solver name or an already-computed scan.
Method = str | EfficiencyScan


def sweep(
    problem: Problem,
    illumination: Illumination,
    wavelengths: ArrayLike,
    methods: Sequence[Method],
    options: dict[str, dict[str, Any]] | None = None,
) -> list[EfficiencyScan]:
    """Solve ``problem`` with each method on a common wavelength grid.

    Parameters
    ----------
    methods
        Registered solver names (``"scalar"``), or :class:`EfficiencyScan`
        objects for precomputed data such as an imported PCGrate table. A
        precomputed scan is **resampled onto** ``wavelengths`` -- see
        :func:`align` for how, and for the caveat.
    options
        Per-method keyword arguments, keyed by method name, e.g.
        ``{"scalar": {"quadrature_points": 8192}}``.
    """
    wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
    options = options or {}

    scans: list[EfficiencyScan] = []
    for method in methods:
        if isinstance(method, EfficiencyScan):
            scans.append(_resample(method, wavelengths))
            continue
        solver = get_solver(method)
        scans.append(
            solver.solve(
                problem, illumination, wavelengths, **options.get(method, {})
            )
        )
    return scans


def _resample(scan: EfficiencyScan, wavelengths: NDArray[np.float64]) -> EfficiencyScan:
    """Put a precomputed scan on ``wavelengths`` by nearest-neighbour lookup.

    Nearest-neighbour rather than interpolation, deliberately: efficiency curves
    have sharp features at order passing-off, and interpolating across one
    invents a value that no method produced. Points outside the source range
    are marked non-propagating rather than extrapolated.
    """
    source = scan.wavelengths
    index = np.abs(wavelengths[:, None] - source[None, :]).argmin(axis=1)
    inside = (wavelengths >= source.min()) & (wavelengths <= source.max())

    efficiency = scan.efficiency[index]
    propagating = scan.propagating[index]
    efficiency = np.where(inside[:, None], efficiency, 0.0)
    propagating = propagating & inside[:, None]

    return EfficiencyScan(
        wavelengths=wavelengths,
        orders=scan.orders,
        efficiency=efficiency,
        propagating=propagating,
        provenance=scan.provenance,
    )


@dataclass(frozen=True, slots=True)
class Comparison:
    """Several methods on one wavelength grid and one order set."""

    wavelengths: NDArray[np.float64]
    orders: NDArray[np.int64]
    methods: tuple[str, ...]
    #: ``(n_methods, n_wavelengths, n_orders)``
    efficiency: NDArray[np.float64]
    scans: tuple[EfficiencyScan, ...]

    def order(self, m: int) -> dict[str, NDArray[np.float64]]:
        """Efficiency in order ``m`` versus wavelength, per method."""
        column = int(np.flatnonzero(self.orders == m)[0])
        return {
            name: self.efficiency[i, :, column]
            for i, name in enumerate(self.methods)
        }

    def difference(self, a: str, b: str) -> NDArray[np.float64]:
        """Signed ``a - b`` over the full grid."""
        return self.efficiency[self.methods.index(a)] - self.efficiency[
            self.methods.index(b)
        ]

    def max_abs_difference(self, a: str, b: str) -> float:
        return float(np.abs(self.difference(a, b)).max())

    def summary(self, a: str, b: str) -> dict[str, float]:
        """Where and how much two methods disagree.

        The headline number of the whole project: not "do they agree" but
        "by how much, and where".
        """
        delta = self.difference(a, b)
        flat = int(np.abs(delta).argmax())
        row, column = np.unravel_index(flat, delta.shape)
        return {
            "max_abs_difference": float(np.abs(delta).max()),
            "rms_difference": float(np.sqrt(np.mean(delta**2))),
            "at_wavelength": float(self.wavelengths[row]),
            "at_order": int(self.orders[column]),
        }


def align(scans: Iterable[EfficiencyScan]) -> Comparison:
    """Put several scans on a shared order set.

    All scans must already share a wavelength grid -- :func:`sweep` arranges
    that. Orders are unioned, so a method that resolves fewer orders reports
    zero for the rest rather than being silently truncated to the smallest
    common set.
    """
    scans = tuple(scans)
    if not scans:
        raise ValueError("nothing to compare")

    reference = scans[0].wavelengths
    for scan in scans[1:]:
        if not np.allclose(scan.wavelengths, reference):
            raise ValueError(
                "scans are on different wavelength grids; use sweep() so every "
                "method is evaluated on the same one"
            )

    orders = np.unique(np.concatenate([scan.orders for scan in scans]))
    stacked = np.zeros((len(scans), len(reference), len(orders)))
    for i, scan in enumerate(scans):
        column = np.searchsorted(orders, scan.orders)
        stacked[i][:, column] = scan.efficiency

    names = _unique_names(scans)
    return Comparison(
        wavelengths=reference,
        orders=orders,
        methods=names,
        efficiency=stacked,
        scans=scans,
    )


def _unique_names(scans: Sequence[EfficiencyScan]) -> tuple[str, ...]:
    """Method labels, disambiguated if two scans share a method name."""
    names, seen = [], {}
    for scan in scans:
        base = scan.provenance.method
        seen[base] = seen.get(base, 0) + 1
        names.append(base if seen[base] == 1 else f"{base}#{seen[base]}")
    return tuple(names)


def records(scans: Iterable[EfficiencyScan]) -> list[dict[str, Any]]:
    """Tidy long-form rows across several scans, ready for a DataFrame."""
    return [row for scan in scans for row in scan.to_records()]
