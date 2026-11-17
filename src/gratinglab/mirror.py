r"""A flat mirror: a coating stack's reflectivity, and how deep the light goes.

```python
from gratinglab.mirror import Mirror
from gratinglab.problem import ConformalLayer

m = Mirror(substrate="Si", layers=(ConformalLayer(material="Au", thickness=20),))
m.reflectivity(graze_deg=[1, 2, 3], wavelength_nm=1.0)
m.penetration_depth(graze_deg=2, wavelength_nm=1.0)   # nm, per material
m.bulk_thickness(graze_deg=2, wavelength_nm=1.0)      # when Au might as well be bulk
```

No grooves: this answers what a coating does before a grating is cut into it.
How much it reflects, how far the field reaches into each material, and
whether the outermost layer is thick enough that what lies beneath it no longer
matters, in which case a grating calculation may treat it as semi-infinite
(``Problem.coating``) instead of modelling the stack.

The stack
=========

``substrate`` and ``layers`` mean what they mean on a
:class:`~gratinglab.problem.Problem`: material names, resolved by
:func:`gratinglab.materials.lookup`, and layers **innermost first**, so the
last one is the surface the light meets. Above it is vacuum. ``roughness`` is
one RMS value applied at every interface, as Névot–Croce on the amplitude.

Reflection is Parratt's recursion, from the substrate up:

.. math:: X_j = \frac{r_{j,j+1} + X_{j+1}\, e^{2 i k k_{z,j+1} d_{j+1}}}
                     {1 + r_{j,j+1} X_{j+1}\, e^{2 i k k_{z,j+1} d_{j+1}}},

with :math:`X = 0` in the substrate, :math:`k_{z,j} = \sqrt{n_j^2 - \cos^2\theta}`
(the decaying branch, as in :mod:`.materials.fresnel`), and :math:`r_{j,j+1}`
the single-interface Fresnel coefficient, multiplied by
:math:`e^{-2 k^2 k_{z,j} k_{z,j+1}\sigma^2}` for roughness. With no layers it
is exactly :func:`gratinglab.materials.fresnel.amplitude`; a test holds it to
that. Angles are grazing, from the surface, and ``s`` and ``p`` are local to
the mirror, as in :mod:`.materials.fresnel`.

How deep
========

:meth:`Mirror.penetration_depth` is the depth at which the intensity of the
transmitted wave in each material, taken as semi-infinite, falls to
:math:`1/e`: :math:`1 / (2 k\, \mathrm{Im}\, k_z)`. Below the critical angle
it is short and set by the angle (the wave is evanescent); above it, it is set
by absorption and is far longer.

:meth:`Mirror.bulk_thickness` asks the question that matters for modelling:
how thick must the outermost layer be before the reflectivity is within
``tolerance`` of that layer's own bulk value, at this angle and wavelength,
**with the real stack underneath**. It sweeps the thickness rather than
assuming the exponential, because interference fringes from the layer below
make the approach to bulk oscillate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .problem import ConformalLayer

__all__ = ["Mirror", "stack_amplitude", "Polarization"]

Polarization = Literal["s", "p", "unpolarized"]


def _kz(n: NDArray[np.complex128], graze: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Normal wavevector over k; NumPy's principal root has Im >= 0, which is
    the wave decaying into the material."""
    return np.sqrt(n**2 - np.cos(graze) ** 2)


def stack_amplitude(
    indices: "list[ArrayLike]",
    thicknesses: "list[float]",
    graze_rad: ArrayLike,
    wavelength_nm: ArrayLike,
    roughness_nm: float = 0.0,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """``(r_s, r_p)`` of a stack under vacuum.

    ``indices`` runs from the outermost medium down to the substrate, each an
    index or an array of them broadcast against ``graze_rad`` and
    ``wavelength_nm``; ``thicknesses`` are the films', in nm, in the same
    order, so there is one fewer of them. A thickness may itself be an array,
    which is how :meth:`Mirror.bulk_thickness` sweeps one in a single call.
    """
    if len(thicknesses) != len(indices) - 1:
        raise ValueError("one thickness per film: len(indices) - 1")
    graze = np.asarray(graze_rad, dtype=np.float64)
    k = 2 * np.pi / np.asarray(wavelength_nm, dtype=np.float64)
    media = [np.complex128(1.0) + 0 * graze, *(np.asarray(n, dtype=np.complex128) for n in indices)]
    kz = [_kz(n, graze) for n in media]
    out = []
    for pol in ("s", "p"):
        x = np.complex128(0.0)
        for j in range(len(media) - 2, -1, -1):
            a, b = kz[j], kz[j + 1]
            if pol == "s":
                r = (a - b) / (a + b)
            else:
                r = (media[j + 1] ** 2 * a - media[j] ** 2 * b) / (media[j + 1] ** 2 * a + media[j] ** 2 * b)
            if roughness_nm:
                r = r * np.exp(-2 * k**2 * a * b * roughness_nm**2)
            # Film j+1 lies between interface j and the one below it; the
            # substrate (the last medium) has no thickness and X = 0 beneath.
            if j + 1 < len(media) - 1:
                phase = np.exp(2j * k * b * np.asarray(thicknesses[j], dtype=np.float64))
                x = (r + x * phase) / (1 + r * x * phase)
            else:
                x = r
        out.append(np.asarray(x, dtype=np.complex128))
    return out[0], out[1]


@dataclass(frozen=True)
class Mirror:
    """A flat ``substrate`` under ``layers`` (innermost first), under vacuum."""

    substrate: str
    layers: tuple[ConformalLayer, ...] = ()
    roughness: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "layers", tuple(self.layers))
        if self.roughness < 0:
            raise ValueError(f"roughness must be non-negative, got {self.roughness}")

    # -- materials ---------------------------------------------------------

    @staticmethod
    def _table(name: str):
        from .materials import lookup  # cached by file there

        return lookup(name)

    def _media(self) -> list[str]:
        """Material names from the surface down to the substrate."""
        return [layer.material for layer in reversed(self.layers)] + [self.substrate]

    def _indices(self, wavelength_nm: NDArray[np.float64]) -> list[NDArray[np.complex128]]:
        return [self._table(name).n(wavelength_nm) for name in self._media()]

    # -- reflection --------------------------------------------------------

    def amplitude(self, graze_deg: ArrayLike, wavelength_nm: ArrayLike):
        """``(r_s, r_p)``, broadcast over ``graze_deg`` and ``wavelength_nm``."""
        graze, wl = np.broadcast_arrays(
            np.radians(np.asarray(graze_deg, dtype=np.float64)),
            np.asarray(wavelength_nm, dtype=np.float64),
        )
        thick = [layer.thickness for layer in reversed(self.layers)]
        return stack_amplitude(self._indices(wl), thick, graze, wl, self.roughness)

    def reflectivity(
        self, graze_deg: ArrayLike, wavelength_nm: ArrayLike, polarization: Polarization = "unpolarized"
    ) -> NDArray[np.float64]:
        """Specular reflectivity, in [0, 1]."""
        return _intensity(*self.amplitude(graze_deg, wavelength_nm), polarization)

    # -- depth -------------------------------------------------------------


    def penetration_depth(self, graze_deg: ArrayLike, wavelength_nm: ArrayLike) -> dict[str, NDArray[np.float64]]:
        """1/e intensity depth in nm in each material, each taken as bulk,
        from the surface material down to the substrate."""
        graze, wl = np.broadcast_arrays(
            np.radians(np.asarray(graze_deg, dtype=np.float64)),
            np.asarray(wavelength_nm, dtype=np.float64),
        )
        k = 2 * np.pi / wl
        out = {}
        with np.errstate(divide="ignore"):
            for name, n in zip(self._media(), self._indices(wl)):
                out.setdefault(name, 1.0 / (2 * k * np.imag(_kz(n, graze))))
        return out

    def bulk_thickness(
        self,
        graze_deg: float,
        wavelength_nm: float,
        tolerance: float = 1e-3,
        polarization: Polarization = "unpolarized",
    ) -> float:
        """Thinnest outermost layer, in nm, beyond which the reflectivity stays
        within ``tolerance`` (relative) of that layer as bulk, with the rest of
        the stack beneath it. ``inf`` when it never settles within forty 1/e
        depths, as for a material with no absorption above its critical angle.
        """
        if not self.layers:
            raise ValueError("a bare substrate is already bulk; add a layer to ask how thick it must be")
        top = self.layers[-1]
        target = float(Mirror(top.material, roughness=self.roughness).reflectivity(graze_deg, wavelength_nm, polarization))
        depth = float(self.penetration_depth(graze_deg, wavelength_nm)[top.material])
        if not np.isfinite(depth):
            return float("inf")
        # Fine enough for the fringes (period lambda / 2 Re kz) and the decay.
        graze = np.radians(graze_deg)
        wl = np.float64(wavelength_nm)
        indices = self._indices(wl)
        fringe = wl / (2 * max(float(np.real(_kz(indices[0], graze))), 1e-12))
        step = min(fringe, depth) / 20
        thickness = np.arange(1, int(min(40 * depth / step, _MAX_SWEEP)) + 1) * step
        rest = [layer.thickness for layer in reversed(self.layers[:-1])]
        values = _intensity(*stack_amplitude(indices, [thickness, *rest], graze, wl, self.roughness), polarization)
        off = np.abs(values - target) > tolerance * max(target, 1e-300)
        if off[-1]:
            return float("inf")
        if not off.any():
            return float(thickness[0])
        return float(thickness[np.flatnonzero(off)[-1] + 1])


#: Most thicknesses :meth:`Mirror.bulk_thickness` evaluates in one sweep.
_MAX_SWEEP = 400_000


def _intensity(r_s, r_p, polarization: Polarization) -> NDArray[np.float64]:
    if polarization == "s":
        return np.abs(r_s) ** 2
    if polarization == "p":
        return np.abs(r_p) ** 2
    if polarization == "unpolarized":
        return 0.5 * (np.abs(r_s) ** 2 + np.abs(r_p) ** 2)
    raise ValueError(f"polarization must be 's', 'p' or 'unpolarized', got {polarization!r}")
