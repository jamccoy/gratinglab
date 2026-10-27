"""The layered structure RCWA solves, and how a profile becomes one.

A :class:`Stack` is the whole problem along the normal: a superstrate, any
homogeneous layers above the grating region, the grating's slabs, any
homogeneous layers below it, and a substrate. Every slab has the same shape,
a background index plus segments of other indices along ``x``, so a
homogeneous film, a lamellar ridge and a staircase step with a conformal
coating on it are one type. That is what keeps a layer stack from being a
retrofit (``docs/theory/rcwa.md`` section 6).

Positions along ``x`` are fractions of the period, ``xi = x / period``.
Indices are complex refractive indices at one wavelength; the backend builds
a stack per wavelength, since a dispersive material changes only the indices.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...profiles import Lamellar

__all__ = ["Slab", "Stack", "Groove", "groove_geometry", "fourier_matrix"]


@dataclass(frozen=True, slots=True)
class Slab:
    """One slab between two planes normal to the grating.

    ``segments`` are ``(xi_start, xi_end, index)``, non-overlapping, where
    the slab holds ``index`` instead of ``background``. An interval may
    start below 0 or end above 1; the structure is periodic, so it means the
    same thing wrapped. No segments is a homogeneous layer.
    """

    thickness: float  # nm
    background: complex
    segments: tuple[tuple[float, float, complex], ...] = ()

    @property
    def homogeneous(self) -> bool:
        return not self.segments or all(
            index == self.background for _, _, index in self.segments
        )


@dataclass(frozen=True, slots=True)
class Stack:
    """Superstrate, slabs top to bottom, substrate.

    ``above`` and ``below`` are the homogeneous layers on either side of the
    grating region. Nothing in the solve distinguishes them from the grating
    slabs; they are kept apart so a caller can say which is which and a
    reader of a stack can see the grating.
    """

    superstrate: complex
    grating: tuple[Slab, ...]
    substrate: complex
    above: tuple[Slab, ...] = ()
    below: tuple[Slab, ...] = ()

    @property
    def slabs(self) -> tuple[Slab, ...]:
        return self.above + self.grating + self.below


@dataclass(frozen=True, slots=True)
class Groove:
    """The grating region's geometry, independent of the materials.

    ``slabs`` are ``(thickness_nm, intervals)`` top to bottom, where
    ``intervals`` are the ``(xi_start, xi_end)`` spans holding the grating
    material. ``exact`` is False for a staircase.
    """

    slabs: tuple[tuple[float, tuple[tuple[float, float], ...]], ...]
    exact: bool

    def stack(self, ridge: complex, *, superstrate: complex = 1.0) -> Stack:
        """Vacuum above, ridges of ``ridge`` in vacuum, ``ridge`` below.

        The mapping a :class:`~gratinglab.problem.Problem` makes today: one
        material, cut into, with nothing on it (``rcwa.md`` section 1).
        """
        grating = tuple(
            Slab(
                thickness=thickness,
                background=superstrate,
                segments=tuple((a, b, ridge) for a, b in intervals),
            )
            for thickness, intervals in self.slabs
        )
        return Stack(superstrate=superstrate, grating=grating, substrate=ridge)


#: Bisection steps refining a staircase wall from its sampled position. The
#: bracket starts at a few samples wide (~1e-4 of the period), so 40 halvings
#: put the wall within ~1e-16, which is below anything a truncation can see.
_BISECTIONS = 40


def _refine(profile, guess: float, level: float, step: float) -> float:
    """Where ``height`` crosses ``level`` near ``guess``, by bisection.

    Falls back to the guess when the bracket holds no sign change, which
    happens only for a profile that touches the level without crossing it.
    """
    low, high = guess - 2.0 * step, guess + 2.0 * step
    f_low = float(profile.height(low)) - level
    f_high = float(profile.height(high)) - level
    if (f_low >= 0.0) == (f_high >= 0.0):
        return guess
    for _ in range(_BISECTIONS):
        middle = 0.5 * (low + high)
        if (float(profile.height(middle)) - level >= 0.0) == (f_low >= 0.0):
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def groove_geometry(profile, period: float, slices: int) -> Groove:
    """The grating region of ``profile``, as slabs top to bottom.

    A :class:`~gratinglab.profiles.Lamellar` groove is one exact slab, whatever
    ``slices`` says. Anything else is a staircase of ``slices`` equal slabs
    by the mid-height rule of ``Profile.slice_layers``, with each wall moved
    from its sampled position to the profile's own crossing.

    The profile's ``t`` runs along ``-d`` (``docs/conventions.md`` section 3),
    so ``t`` becomes ``xi = 1 - t``, the placement the integral solver uses.
    """
    if isinstance(profile, Lamellar):
        intervals = ((1.0 - profile.duty_cycle, 1.0),)
        return Groove(slabs=((profile.depth * period, intervals),), exact=True)

    layers = profile.slice_layers(slices)
    samples = max(2048, 64 * slices)  # the sampling slice_layers used
    step = 1.0 / samples
    slabs = []
    for layer in reversed(layers):  # slice_layers is bottom to top
        level = 0.5 * (layer.y_lower + layer.y_upper)
        intervals = []
        for start, end in layer.intervals:
            if end - start >= 1.0:
                intervals.append((0.0, 1.0))
                continue
            start = _refine(profile, start, level, step)
            end = _refine(profile, end, level, step)
            intervals.append((1.0 - end, 1.0 - start))
        slabs.append((layer.thickness * period, tuple(intervals)))
    return Groove(slabs=tuple(slabs), exact=False)


def _coefficients(slab: Slab, power: int, harmonics: int) -> NDArray[np.complex128]:
    """Fourier coefficients ``c_n``, ``n = -2M..2M``, of ``eps(xi)**power``.

    ``eps(xi) = sum_n c_n exp(2 pi i n xi)``, so ``c_n`` is the period average
    of ``eps**power * exp(-2 pi i n xi)``: exact for a piecewise-constant slab.
    """
    n = np.arange(-2 * harmonics, 2 * harmonics + 1)
    base = complex(slab.background) ** 2
    background = base**power
    coefficients = np.where(n == 0, background, 0.0).astype(np.complex128)
    nonzero = n != 0
    for start, end, index in slab.segments:
        value = (complex(index) ** 2) ** power - background
        if value == 0.0:
            continue
        integral = np.empty(len(n), dtype=np.complex128)
        integral[~nonzero] = end - start
        k = -2j * np.pi * n[nonzero]
        integral[nonzero] = (np.exp(k * end) - np.exp(k * start)) / k
        coefficients += value * integral
    return coefficients


def fourier_matrix(slab: Slab, harmonics: int, *, inverse: bool) -> NDArray[np.complex128]:
    """The ``(2M+1)``-square Toeplitz matrix for a product ``eps E``.

    ``inverse=False`` is Laurent's rule, ``[[eps]]``. ``inverse=True`` is the
    inverse rule, ``[[1/eps]]^-1``, which is the correct one for a field
    component that jumps where ``eps`` does (Li 1996; ``rcwa.md`` section 3).
    """
    size = 2 * harmonics + 1
    coefficients = _coefficients(slab, -1 if inverse else 1, harmonics)
    rows = np.arange(size)
    toeplitz = coefficients[(rows[:, None] - rows[None, :]) + 2 * harmonics]
    return np.linalg.inv(toeplitz) if inverse else toeplitz
