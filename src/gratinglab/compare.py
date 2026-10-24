"""Run one problem through several methods and line the answers up.

This is the central purpose of the project. A comparison is meaningful only if
every method saw the same physics, so :func:`sweep` takes one
:class:`~gratinglab.problem.Problem` and one
:class:`~gratinglab.illumination.Illumination` and hands them to each backend
unchanged.

Precomputed results are first-class participants. An imported PCGrate table is
passed in exactly like a live solver, so scalar-versus-integral-method plots are
available without a PCGrate licence.

Measured data joins at :func:`align`, not at :func:`sweep`. A
:class:`~gratinglab.measured.MeasuredEfficiency` is sparse, with a sigma on
every point, so resampling it onto a model grid the way a precomputed table is
resampled would invent points that were never measured. :func:`align` carries
measurements beside the dense stack instead, so measured, scalar and integral
results land on one :class:`Comparison`, and :meth:`Comparison.residuals`
evaluates a model at each measured point. To make that exact, sweep the
models on :meth:`MeasuredEfficiency.wavelength_grid`, or on its union with a
denser grid when a smooth curve is wanted as well.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .illumination import Illumination
from .measured import MeasuredEfficiency, OrderPoints
from .problem import Problem
from .result import EfficiencyScan
from .solvers.base import Progress, get_solver

__all__ = ["sweep", "align", "Comparison", "Residuals", "records"]

#: A method is either a registered solver name or an already-computed scan.
Method = str | EfficiencyScan

#: Anything :func:`align` and :func:`records` take: a dense scan or a
#: sparse measurement.
Participant = EfficiencyScan | MeasuredEfficiency


def sweep(
    problem: Problem,
    illumination: Illumination,
    wavelengths: ArrayLike,
    methods: Sequence[Method],
    options: dict[str, dict[str, Any]] | None = None,
    *,
    progress: "Progress | None" = None,
) -> list[EfficiencyScan]:
    """Solve ``problem`` with each method on a common wavelength grid.

    Parameters
    ----------
    methods
        Registered solver names (``"scalar"``), or :class:`EfficiencyScan`
        objects for precomputed data such as an imported PCGrate table. A
        precomputed scan is resampled onto ``wavelengths``; see
        :func:`align` for how, and for the caveat.
    options
        Per-method keyword arguments, keyed by method name, e.g.
        ``{"scalar": {"quadrature_points": 8192}}``.
    progress
        Called ``(done, total)`` over the whole sweep instead of per method,
        so a bar driven by it never restarts at a method boundary. Units are
        wavelengths, and ``total`` is ``len(methods) * len(wavelengths)``.

        Raising from it aborts the sweep (see
        :class:`~gratinglab.solvers.base.SolveCancelled`). It is a keyword of
        its own instead of an entry in ``options``, so that it can never
        collide with a solver's own parameter name.
    """
    wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
    options = options or {}
    per_method = len(wavelengths)
    total = len(methods) * per_method

    scans: list[EfficiencyScan] = []
    for index, method in enumerate(methods):
        finished = index * per_method

        if isinstance(method, MeasuredEfficiency):
            raise TypeError(
                "a MeasuredEfficiency is not a method to sweep: resampling "
                "sparse measured points onto a model grid would invent points "
                "nobody measured. Pass it to align() beside the swept scans, "
                "and sweep on measured.wavelength_grid() to compare at the "
                "measured wavelengths."
            )

        if isinstance(method, EfficiencyScan):
            # Resampling is instant, so this block completes at once. Stepping
            # the counter by the whole block is honest; pretending a
            # precomputed table has wavelengths still to compute is not.
            scans.append(_resample(method, wavelengths))
            if progress is not None:
                progress(finished + per_method, total)
            continue

        solver = get_solver(method)
        forwarded = _progress_option(solver, progress, finished, total)
        scans.append(
            solver.solve(
                problem, illumination, wavelengths,
                **forwarded, **options.get(method, {}),
            )
        )
        if progress is not None and not forwarded:
            # Only for a solver that did *not* report: this is what stops the
            # bar stalling at a method boundary. One that did report already
            # ended its block exactly here, and repeating the call would emit
            # a duplicate for no reason.
            progress(finished + per_method, total)
    return scans


def _progress_option(
    solver, progress: "Progress | None", finished: int, total: int
) -> dict[str, Any]:
    """``{"progress": ...}`` for a solver that declared it, else ``{}``.

    A backend that has not been updated, including a third-party one arriving
    through the ``gratinglab.solvers`` entry point, has an explicit signature
    that would raise ``TypeError`` on an unexpected keyword. For that reason
    :attr:`~gratinglab.solvers.base.Capabilities.reports_progress` is declared
    by each solver and is not assumed.
    """
    if progress is None or not solver.capabilities.reports_progress:
        return {}
    return {"progress": lambda done, _per_method: progress(finished + done, total)}


def _resample(scan: EfficiencyScan, wavelengths: NDArray[np.float64]) -> EfficiencyScan:
    """Put a precomputed scan on ``wavelengths`` by nearest-neighbour lookup.

    Nearest-neighbour lookup is used in place of interpolation because
    efficiency curves have sharp features at order passing-off, and
    interpolating across one invents a value that no method produced. Points
    outside the source range are marked non-propagating and are not
    extrapolated.
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
class Residuals:
    """Measured minus modelled efficiency at every measured point.

    Parallel 1-D arrays in the measurement's own order (by order, then
    wavelength). ``pull`` is ``residual / sigma``, the residual in units of
    its own uncertainty, and is NaN where ``sigma`` is zero, since a point
    with no stated uncertainty cannot say how far off it is in those units.
    ``propagating`` is the *model's* flag: a measured order that the model
    calls evanescent is a disagreement worth seeing, and its ``model`` is 0.
    """

    measurement: str
    method: str
    order: NDArray[np.int64]
    wavelength: NDArray[np.float64]
    measured: NDArray[np.float64]
    sigma: NDArray[np.float64]
    model: NDArray[np.float64]
    propagating: NDArray[np.bool_]

    @property
    def residual(self) -> NDArray[np.float64]:
        return self.measured - self.model

    @property
    def pull(self) -> NDArray[np.float64]:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.sigma > 0.0, self.residual / self.sigma, np.nan)

    def __len__(self) -> int:
        return len(self.order)

    def summary(self) -> dict[str, float]:
        """How far the model is from the data, and where it is furthest.

        ``chi2`` sums ``pull**2`` over the points that have a sigma, and
        ``n_with_sigma`` says how many that was, so the number is never read
        without its denominator. ``at_wavelength`` and ``at_order`` locate the
        largest absolute residual.
        """
        residual = self.residual
        pull = self.pull
        finite = np.isfinite(pull)
        worst = int(np.abs(residual).argmax())
        return {
            "max_abs_residual": float(np.abs(residual).max()),
            "rms_residual": float(np.sqrt(np.mean(residual**2))),
            "chi2": float(np.sum(pull[finite] ** 2)),
            "n_points": float(len(self)),
            "n_with_sigma": float(np.count_nonzero(finite)),
            "at_wavelength": float(self.wavelength[worst]),
            "at_order": int(self.order[worst]),
        }


@dataclass(frozen=True, slots=True)
class Comparison:
    """Several methods on one wavelength grid and one order set, with any
    measurements beside them.

    ``methods`` and ``efficiency`` describe the dense scans only.
    Measurements are sparse and keep their own points, so they are named in
    ``measurements`` and held in ``measured`` instead of being forced onto
    the grid. ``orders`` is the union over both.
    """

    wavelengths: NDArray[np.float64]
    orders: NDArray[np.int64]
    methods: tuple[str, ...]
    #: ``(n_methods, n_wavelengths, n_orders)``
    efficiency: NDArray[np.float64]
    scans: tuple[EfficiencyScan, ...]
    measurements: tuple[str, ...] = ()
    measured: tuple[MeasuredEfficiency, ...] = ()

    def measured_order(self, m: int) -> dict[str, OrderPoints]:
        """Measured points in order ``m``, with error bars, per measurement.

        A measurement with no point in order ``m`` is left out rather than
        reported as an empty set. The partner of :meth:`order` for plotting.
        """
        return {
            name: data.points(m)
            for name, data in zip(self.measurements, self.measured)
            if m in data.orders
        }

    def residuals(self, measurement: str, method: str) -> Residuals:
        """The named method evaluated at every point of the named measurement.

        The model is read off the comparison grid at each measured
        wavelength, which must be on that grid. Nothing here interpolates,
        because a value between two model points is one no method produced;
        sweep on the measurement's ``wavelength_grid()`` to arrange it.
        """
        data = self.measured[self.measurements.index(measurement)]
        scan = self.scans[self.methods.index(method)]

        row = np.abs(self.wavelengths[None, :] - data.wavelength[:, None]).argmin(axis=1)
        matched = np.isclose(self.wavelengths[row], data.wavelength, rtol=1e-9, atol=0.0)
        if not matched.all():
            missing = np.unique(data.wavelength[~matched])
            raise ValueError(
                f"{len(missing)} measured wavelength(s) of {measurement!r} are not "
                f"on the comparison grid (first {missing[0]:.6g} nm). Sweep on "
                "measured.wavelength_grid(), or its union with your grid, so the "
                "model is evaluated where the data is."
            )

        model = np.zeros(len(data))
        propagating = np.zeros(len(data), dtype=bool)
        for k, (i, m) in enumerate(zip(row, data.order)):
            hit = np.flatnonzero(scan.orders == m)
            if hit.size:
                model[k] = scan.efficiency[i, hit[0]]
                propagating[k] = scan.propagating[i, hit[0]]

        return Residuals(
            measurement=measurement,
            method=method,
            order=data.order,
            wavelength=data.wavelength,
            measured=data.value,
            sigma=data.sigma,
            model=model,
            propagating=propagating,
        )

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

        This is the headline number of the project. It answers the question
        "by how much, and where" in place of "do they agree".
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


def align(scans: Iterable[Participant]) -> Comparison:
    """Put several scans, and any measurements, on a shared order set.

    All scans must already share a wavelength grid, which :func:`sweep`
    arranges. Orders are unioned, so a method that resolves fewer orders
    reports zero for the rest; the comparison is not truncated to the
    smallest common set.

    A :class:`~gratinglab.measured.MeasuredEfficiency` may appear anywhere in
    the list. It is not put on the grid but carried beside it (see
    :class:`Comparison`), and its orders join the union. Labels are made
    unique across scans and measurements together, so two measurements both
    called ``"measured"`` become ``"measured"`` and ``"measured#2"``.
    """
    participants = tuple(scans)
    if not participants:
        raise ValueError("nothing to compare")
    for item in participants:
        if not isinstance(item, (EfficiencyScan, MeasuredEfficiency)):
            raise TypeError(
                "align() takes EfficiencyScan and MeasuredEfficiency, got "
                f"{type(item).__name__}"
            )

    names = _unique_names(participants)
    scans = tuple(p for p in participants if isinstance(p, EfficiencyScan))
    measured = tuple(p for p in participants if isinstance(p, MeasuredEfficiency))
    method_names = tuple(
        n for n, p in zip(names, participants) if isinstance(p, EfficiencyScan)
    )
    measurement_names = tuple(
        n for n, p in zip(names, participants) if isinstance(p, MeasuredEfficiency)
    )

    if not scans:
        orders = np.unique(np.concatenate([m.orders for m in measured]))
        return Comparison(
            wavelengths=np.zeros(0),
            orders=orders,
            methods=(),
            efficiency=np.zeros((0, 0, len(orders))),
            scans=(),
            measurements=measurement_names,
            measured=measured,
        )

    reference = scans[0].wavelengths
    for scan in scans[1:]:
        if not np.allclose(scan.wavelengths, reference):
            raise ValueError(
                "scans are on different wavelength grids; use sweep() so every "
                "method is evaluated on the same one"
            )

    orders = np.unique(
        np.concatenate([scan.orders for scan in scans] + [m.orders for m in measured])
    )
    stacked = np.zeros((len(scans), len(reference), len(orders)))
    for i, scan in enumerate(scans):
        column = np.searchsorted(orders, scan.orders)
        stacked[i][:, column] = scan.efficiency

    return Comparison(
        wavelengths=reference,
        orders=orders,
        methods=method_names,
        efficiency=stacked,
        scans=scans,
        measurements=measurement_names,
        measured=measured,
    )


def _unique_names(scans: Sequence[Participant]) -> tuple[str, ...]:
    """Labels, disambiguated if two participants share a name."""
    names, seen = [], {}
    for scan in scans:
        base = scan.provenance.method
        seen[base] = seen.get(base, 0) + 1
        names.append(base if seen[base] == 1 else f"{base}#{seen[base]}")
    return tuple(names)


def records(scans: Iterable[Participant]) -> list[dict[str, Any]]:
    """Tidy long-form rows across several scans, ready for a DataFrame.

    A measurement contributes one row per measured point, with a ``sigma``
    key that model rows do not have, so a DataFrame built from a mixed list
    has ``sigma`` as NaN on the model rows.
    """
    return [row for scan in scans for row in scan.to_records()]
