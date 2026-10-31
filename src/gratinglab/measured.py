"""Measured diffraction efficiency, with an uncertainty on every point.

A measurement is not a model, and this module keeps the two apart.
:class:`~gratinglab.result.EfficiencyScan` is what solvers produce: a dense
wavelength grid, every order in every row, and the invariant that an
evanescent order is exactly 0.0. Beamline data has none of that shape. Each
order is measured only where it propagates and is bright enough to centroid,
points are dropped after inspection, and each one carries a standard
uncertainty. So :class:`MeasuredEfficiency` is sparse and long-form: four
parallel 1-D arrays, one entry per measured ``(order, wavelength)`` point.

The uncertainty is a plain float array, ``sigma``, beside ``value``. The
reduction in ``als_eff_meas`` propagates errors with the ``uncertainties``
package, and converting its ufloats into two float arrays is that package's
job (a ``to_measured()`` there); ``uncertainties`` never becomes a dependency
here. ``sigma`` is a standard uncertainty (one standard deviation) in
absolute efficiency, the same units as ``value``.

Values are kept as measured. A weak order after dark subtraction can come out
slightly negative, and that is allowed: clamping it to zero would bias every
residual taken against it, which is "report, never rescale" applied to data
instead of to a model.

The design, and which parts of ``als_eff_meas`` move here and which stay, is
in ``docs/roadmap.md`` under "Measured efficiency and the als_eff_meas split".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .illumination import Illumination
from .materials.optical import HC_EV_NM
from .result import _frozen_array

__all__ = ["MeasurementProvenance", "MeasuredEfficiency", "OrderPoints"]

_POLARIZATIONS = ("TE", "TM", "unpolarized")


@dataclass(frozen=True, slots=True)
class MeasurementProvenance:
    """Where a measurement was taken, and under what conditions.

    Attributes
    ----------
    label
        The name this data carries in a comparison: ``"measured"`` unless
        two measurements need telling apart.
    facility
        For example ``"ALS"``. ``None`` if not recorded.
    campaign
        The run it belongs to, for example ``"ogre_laminar_I_oct21"``.
    polarization
        ``"TE"``, ``"TM"`` or ``"unpolarized"`` in gratinglab's convention
        (referenced to the grooves, ``docs/conventions.md`` §7), or ``None``
        when it is not known. Unknown is a legitimate answer and is never
        filled in: a synchrotron beam is strongly linearly polarized, but
        whether that is TE or TM at the grating depends on the mount.
    geometry
        The mount, in degrees. ``alpha_deg`` and ``gamma_deg`` mean what
        they mean on :class:`~gratinglab.illumination.Illumination`, which
        :meth:`MeasuredEfficiency.illumination` relies on. Anything else
        (``yaw_deg``, ``eta_deg``, their uncertainties) is carried for the
        record.

        A beamline's own record may not use that convention, and may not
        agree with itself: the June 2017 SCIL campaign's
        ``results/arc_angles.txt`` records α negative where this convention
        needs it positive, and its ``arc_angles_check.txt`` gives a γ a full
        σ away. So take the sign from the propagation count rather than from
        the file — the wrong sign strands measured orders on evanescent
        branches, 59 of 123 points in that campaign — and treat the recorded
        angles as approximate (``docs/findings.md``, "The June 2017 SCIL
        campaign").
    source
        Path or URL of the reduced data.
    notes
        Free-form extras needed to reproduce, such as the reduction's version
        or the points it excluded.
    """

    label: str = "measured"
    facility: str | None = None
    campaign: str | None = None
    polarization: str | None = None
    geometry: Mapping[str, float] = field(default_factory=dict)
    source: str | None = None
    notes: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.polarization is not None and self.polarization not in _POLARIZATIONS:
            raise ValueError(
                f"polarization must be one of {_POLARIZATIONS} or None, got "
                f"{self.polarization!r}"
            )
        # A copy, so the caller's dict cannot change a frozen record later.
        object.__setattr__(
            self, "geometry", {str(k): float(v) for k, v in dict(self.geometry).items()}
        )

    @property
    def method(self) -> str:
        """The comparison label, under the name a scan's provenance uses."""
        return self.label


@dataclass(frozen=True, slots=True)
class OrderPoints:
    """The measured points of one order, ascending in wavelength."""

    order: int
    wavelengths: NDArray[np.float64]
    value: NDArray[np.float64]
    sigma: NDArray[np.float64]

    def __len__(self) -> int:
        return len(self.wavelengths)


@dataclass(frozen=True, slots=True)
class MeasuredEfficiency:
    """Absolute efficiency measured per order, sparse, with a σ per point.

    ``order``, ``wavelength``, ``value`` and ``sigma`` are parallel 1-D
    arrays with one entry per measured point. They are stored sorted by order
    and then wavelength, whatever order they were passed in, and each
    ``(order, wavelength)`` pair appears at most once.

    Wavelengths are in nm (``docs/conventions.md`` §1). Beamline data is
    usually labelled by photon energy; :meth:`from_energies` converts.
    """

    order: NDArray[np.int64]
    wavelength: NDArray[np.float64]
    value: NDArray[np.float64]
    sigma: NDArray[np.float64]
    provenance: MeasurementProvenance = field(default_factory=MeasurementProvenance)

    def __post_init__(self) -> None:
        arrays = {
            "order": np.asarray(self.order),
            "wavelength": np.asarray(self.wavelength, dtype=np.float64),
            "value": np.asarray(self.value, dtype=np.float64),
            "sigma": np.asarray(self.sigma, dtype=np.float64),
        }
        for name, array in arrays.items():
            if array.ndim != 1:
                raise ValueError(f"{name} must be 1-D, got shape {array.shape}")
        lengths = {name: len(a) for name, a in arrays.items()}
        if len(set(lengths.values())) != 1:
            raise ValueError(
                f"order, wavelength, value and sigma differ in length: {lengths}"
            )
        if not lengths["order"]:
            raise ValueError("a measurement needs at least one point")

        order = arrays["order"]
        if not np.issubdtype(order.dtype, np.integer):
            as_float = order.astype(np.float64)
            if not (np.isfinite(as_float).all() and np.array_equal(as_float, np.round(as_float))):
                raise ValueError("order contains non-integer values")
        order = order.astype(np.int64)

        wavelength, value, sigma = arrays["wavelength"], arrays["value"], arrays["sigma"]
        for name, array in (("wavelength", wavelength), ("value", value), ("sigma", sigma)):
            if not np.isfinite(array).all():
                raise ValueError(
                    f"{name} contains NaN or infinity. A point that was not "
                    "measured is left out, not recorded as NaN."
                )
        if (wavelength <= 0.0).any():
            raise ValueError("wavelength must be positive (nm)")
        if (sigma < 0.0).any():
            raise ValueError("sigma must be non-negative; it is a standard uncertainty")

        sort = np.lexsort((wavelength, order))
        order, wavelength = order[sort], wavelength[sort]
        repeated = (np.diff(order) == 0) & (np.diff(wavelength) == 0)
        if repeated.any():
            i = int(np.flatnonzero(repeated)[0])
            raise ValueError(
                f"order {order[i]} is measured twice at {wavelength[i]} nm; "
                "combine repeated measurements before building the record"
            )

        object.__setattr__(self, "order", _frozen_array(order, np.int64))
        object.__setattr__(self, "wavelength", _frozen_array(wavelength, np.float64))
        object.__setattr__(self, "value", _frozen_array(value[sort], np.float64))
        object.__setattr__(self, "sigma", _frozen_array(sigma[sort], np.float64))

    @classmethod
    def from_energies(
        cls,
        order: ArrayLike,
        energy_ev: ArrayLike,
        value: ArrayLike,
        sigma: ArrayLike,
        provenance: MeasurementProvenance | None = None,
    ) -> "MeasuredEfficiency":
        """Build from photon energies in eV, as beamline data is labelled.

        Converted with ``wavelength_nm = HC_EV_NM / energy_ev``, the constant
        :mod:`gratinglab.materials` uses, so a measured point and the optical
        constants looked up for it agree on its wavelength.
        """
        energy = np.asarray(energy_ev, dtype=np.float64)
        if not (np.isfinite(energy).all() and (energy > 0.0).all()):
            raise ValueError("energy must be finite and positive (eV)")
        return cls(
            order=np.asarray(order),
            wavelength=HC_EV_NM / energy,
            value=value,
            sigma=sigma,
            provenance=provenance or MeasurementProvenance(),
        )

    def __len__(self) -> int:
        return len(self.order)

    @property
    def method(self) -> str:
        """The comparison label, from the provenance."""
        return self.provenance.label

    @property
    def orders(self) -> NDArray[np.int64]:
        """The orders with at least one measured point, ascending."""
        return np.unique(self.order)

    @property
    def energy_ev(self) -> NDArray[np.float64]:
        """Photon energy of each point in eV, parallel to :attr:`wavelength`."""
        return HC_EV_NM / self.wavelength

    def points(self, m: int) -> OrderPoints:
        """The measured points of order ``m``.

        Raises ``KeyError`` for an order with no points, rather than returning
        an empty set that a plot would silently draw as nothing.
        """
        hit = self.order == m
        if not hit.any():
            raise KeyError(f"order {m} was not measured; measured: {self.orders.tolist()}")
        return OrderPoints(
            order=int(m),
            wavelengths=self.wavelength[hit],
            value=self.value[hit],
            sigma=self.sigma[hit],
        )

    def wavelength_grid(self) -> NDArray[np.float64]:
        """Every distinct measured wavelength, ascending.

        Sweep a model on this grid, or on its ``np.union1d`` with a denser
        one for a smooth curve, and
        :meth:`gratinglab.compare.Comparison.residuals` then evaluates the
        model exactly at every measured point.
        """
        return np.unique(self.wavelength)

    def illumination(self, polarization: str | None = None) -> Illumination:
        """The :class:`~gratinglab.illumination.Illumination` this was measured in.

        Needs ``alpha_deg`` and ``gamma_deg`` in the geometry. The
        polarization comes from the provenance unless one is given, and if
        neither has it this raises rather than default to TE.
        """
        geometry = self.provenance.geometry
        missing = [k for k in ("alpha_deg", "gamma_deg") if k not in geometry]
        if missing:
            raise ValueError(
                f"geometry has no {', '.join(missing)}; it holds "
                f"{sorted(geometry) or 'nothing'}"
            )
        polarization = polarization or self.provenance.polarization
        if polarization is None:
            raise ValueError(
                "the polarization of this measurement is not recorded; pass one "
                "explicitly rather than let it default"
            )
        return Illumination(
            alpha_deg=geometry["alpha_deg"],
            gamma_deg=geometry["gamma_deg"],
            polarization=polarization,
        )

    def to_records(self) -> list[dict[str, Any]]:
        """Tidy long-form records, one per measured point.

        The keys of :meth:`EfficiencyScan.to_records` plus ``sigma``, so
        modelled and measured rows stack into one table. ``propagating`` is
        True by construction, since an order is measured only where it
        propagates; ``converged`` is None, since it describes a numerical
        method and a measurement is not one.
        """
        return [
            {
                "method": self.method,
                "wavelength": float(w),
                "order": int(m),
                "efficiency": float(v),
                "sigma": float(s),
                "propagating": True,
                "converged": None,
            }
            for m, w, v, s in zip(self.order, self.wavelength, self.value, self.sigma)
        ]
