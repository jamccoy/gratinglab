r"""The profile as the C-method needs it: Fourier series, computed once.

The translation coordinate ``u = y - a(x)`` flattens the boundary, and the
profile then enters the method through two families of Fourier coefficients
over one period, both in the normalised position ``s = x / P``:

- ``slope(j)``, of ``a'(x)``, which fills the Toeplitz matrix of the
  eigenproblem;
- ``plane_wave(beta)``, of ``exp(i beta a(x))``, which is a vacuum plane wave
  written in the new coordinates, used for the propagating orders.

Three ways of computing them, chosen by the profile's type rather than guessed:

- **Polygon** (:class:`~gratinglab.profiles.Blazed`,
  :class:`~gratinglab.profiles.FromProfileData`, and a
  :class:`~gratinglab.profiles.Trapezoidal` with sharp corners). Each straight
  segment integrates in closed form, so a facet corner carries no sampling
  error at all. That matters: ``a'`` jumps there, and an FFT of a sampled jump
  converges at first order in the sample count.
- **Pieces** (a :class:`~gratinglab.profiles.Trapezoidal` with rounded
  corners). The straight runs integrate in closed form as for a polygon, and
  each fillet arc by Gauss-Legendre quadrature over its own interval, where
  the integrand is analytic, with enough nodes for the fastest phase it can
  carry. Sampled over the whole period instead, the curvature jump at every
  tangent point leaves a floor: 8.9e-5 in the slope coefficients of a
  0.005-period fillet at 4096 samples (measured 2026-10-29; ``findings.md``,
  "Trapezoidal through the four solvers").
- **Sampled** (:class:`~gratinglab.profiles.Sinusoidal`, and any other
  profile that offers ``height`` and ``slope``). The trapezoidal rule on a
  smooth periodic integrand converges exponentially, so an FFT of a few
  thousand samples is exact to rounding.

The profile parameter runs against the dispersion direction, ``t = 1 - s``
(``docs/conventions.md``, "Direction of the profile parameter"); this module
is the one place that reflection happens, as ``_boundary.py`` is for the
integral method.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ...problem import Problem
from ...profiles import (
    Blazed,
    FromProfileData,
    Lamellar,
    ProfileRepresentationError,
    Trapezoidal,
)
from ..base import UnsupportedConfiguration

__all__ = ["ProfileSeries", "profile_series"]

#: Fewest samples for the sampled path. The integrand is analytic for a
#: sinusoid, so this is far more than rounding needs; the floor exists so a
#: tiny truncation still sees the profile's own bandwidth.
_MIN_SAMPLES = 4096


def _phi(z: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """``(exp(z) - 1) / z``, accurate as ``z`` goes to zero."""
    z = np.asarray(z, dtype=np.complex128)
    small = np.abs(z) < 1e-4
    out = np.empty_like(z)
    zs = z[small]
    out[small] = 1.0 + zs / 2.0 + zs * zs / 6.0
    zl = z[~small]
    out[~small] = np.expm1(zl) / zl
    return out


@dataclass(frozen=True)
class ProfileSeries:
    """One period of ``a(x)`` in nm, reflected onto ``x = (1 - t) P``.

    ``kind`` is ``"polygon"`` (vertices ``s``, heights ``a``, the polygon
    closing on ``s[0] + 1``), ``"sampled"`` (``s`` equally spaced on
    ``[0, 1)``, with ``slopes`` the sampled ``a'``) or ``"pieces"`` (``s`` and
    ``a`` the vertices as for a polygon, where the segment starting at a
    vertex listed in ``arcs`` is that fillet rather than a straight run).

    ``arcs`` holds, per fillet, the index of the vertex it starts at and the
    profile's own arc piece in the profile's ``t``, so the quadrature reads
    the exact circle rather than a copy of it.
    """

    period: float
    kind: str
    s: NDArray[np.float64]
    a: NDArray[np.float64]
    slopes: NDArray[np.float64] | None = None
    arcs: tuple = ()

    @property
    def max_slope(self) -> float:
        """Largest ``|a'|``, dimensionless."""
        if self.kind == "polygon":
            return float(np.abs(self._segments()[3]).max())
        if self.kind == "pieces":
            # A fillet's slope is monotone, so its largest is at an end.
            ends = [
                np.abs(piece.slope(np.array([piece.x0, piece.x1]))).max()
                for _, piece in self.arcs
            ]
            return float(max(np.abs(self._straight()[3]).max(initial=0.0), *ends))
        return float(np.abs(self.slopes).max())

    @property
    def depth(self) -> float:
        return float(np.ptp(self.a))

    def _segments(self):
        s0 = self.s
        s1 = np.append(self.s[1:], self.s[0] + 1.0)
        a0 = self.a
        a1 = np.append(self.a[1:], self.a[0])
        width = s1 - s0
        slope = (a1 - a0) / (self.period * width)
        return s0, width, a0, slope

    def _straight(self):
        """The straight segments of a ``"pieces"`` series, as ``_segments``."""
        s0, width, a0, slope = self._segments()
        straight = np.ones(len(s0), dtype=bool)
        straight[[index for index, _ in self.arcs]] = False
        return s0[straight], width[straight], a0[straight], slope[straight]

    def _straight_coefficients(self, segments, beta, j, slope_weighted):
        s0, width, a0, slope = segments
        q = beta * self.period * slope[None, :] - 2.0 * np.pi * j[:, None]
        pieces = (
            np.exp(1j * beta * a0[None, :] - 2j * np.pi * j[:, None] * s0[None, :])
            * width[None, :]
            * _phi(1j * q * width[None, :])
        )
        if slope_weighted:
            pieces = pieces * slope[None, :]
        return pieces.sum(axis=1)

    def _arc_coefficients(self, piece, beta, j, slope_weighted):
        """One fillet, by Gauss-Legendre in the angle round its centre.

        ``s = 1 - t`` reverses the interval and ``ds = -dt``, so integrating
        in ``t`` from ``x0`` to ``x1`` carries no sign, and ``a'(x) dx`` is
        ``-dh``. On the circle ``t = cx + r cos(theta)``,
        ``h = cy + r sin(theta)``, so both integrands are entire functions of
        ``theta``, however steep the wall the fillet meets. In ``t`` they are
        not: ``h'`` has a square-root singularity where the circle turns
        vertical, just past the wall end of the arc, and on an 89 deg wall a
        rule in ``t`` was still moving by 7e-5 at four times its nodes. The
        phase changes by at most ``r (2 pi |j| + |beta| P)`` per radian, and
        one node per radian of it, plus a margin, converges to rounding.
        """
        ends = np.array([piece.x0, piece.x1])
        theta0, theta1 = np.arctan2(piece.height(ends) - piece.cy, ends - piece.cx)
        sweep = theta1 - theta0
        rate = piece.r * (
            2.0 * np.pi * float(np.abs(j).max(initial=0)) + abs(beta) * self.period
        )
        u, w = np.polynomial.legendre.leggauss(24 + int(np.ceil(abs(sweep) * rate)))
        theta = theta0 + 0.5 * sweep * (u + 1.0)
        weights = 0.5 * sweep * w
        t = piece.cx + piece.r * np.cos(theta)
        h = piece.cy + piece.r * np.sin(theta)
        if slope_weighted:
            values = weights * -piece.r * np.cos(theta)  # -dh/dtheta
        else:
            values = weights * -piece.r * np.sin(theta)  # dt/dtheta
        values = values * np.exp(1j * beta * self.period * h)
        return np.exp(-2j * np.pi * np.outer(j, 1.0 - t)) @ values

    def _coefficients(self, beta: complex, j: NDArray[np.int64], slope_weighted: bool):
        """``int_0^1 w(s) exp(i beta a - 2 pi i j s) ds`` per segment, summed."""
        if self.kind == "polygon":
            return self._straight_coefficients(self._segments(), beta, j, slope_weighted)
        if self.kind == "pieces":
            total = self._straight_coefficients(self._straight(), beta, j, slope_weighted)
            for _, piece in self.arcs:
                total = total + self._arc_coefficients(piece, beta, j, slope_weighted)
            return total
        samples = len(self.s)
        values = np.exp(1j * beta * self.a)
        if slope_weighted:
            values = values * self.slopes
        spectrum = np.fft.fft(values) / samples
        return spectrum[np.mod(j, samples)]

    def slope(self, j: NDArray[np.int64]) -> NDArray[np.complex128]:
        """Fourier coefficients of ``a'(x)`` at harmonic indices ``j``."""
        return self._coefficients(0.0, np.asarray(j), slope_weighted=True)

    def plane_wave(self, beta: float, j: NDArray[np.int64]) -> NDArray[np.complex128]:
        """Fourier coefficients of ``exp(i beta a(x))`` at indices ``j``."""
        return self._coefficients(beta, np.asarray(j), slope_weighted=False)


def _polygon(problem: Problem, t, y) -> ProfileSeries:
    """A single-valued polyline in ``(t, y)``, normalised, as a series."""
    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    # Repeated t: the same point twice is harmless, two heights at one
    # position is a vertical segment, which has no slope at all.
    steps = np.diff(t)
    vertical = (steps == 0.0) & (np.diff(y) != 0.0)
    if vertical.any():
        raise UnsupportedConfiguration(
            "the profile has a vertical segment (two heights at one position). "
            "The C-method describes the boundary by its slope, which is "
            "infinite there; use the integral method, which parametrises the "
            "curve, rather than a smoothed facet."
        )
    keep = np.append(True, steps > 0.0)
    t, y = t[keep], y[keep]
    # The polyline closes periodically, as np.interp(period=1) reads it;
    # a last point at t = 1 duplicates the first.
    if t[-1] - t[0] >= 1.0:
        if not np.isclose(y[-1], y[0]):
            raise UnsupportedConfiguration(
                "the profile's ends at t = 0 and t = 1 differ in height, so the "
                "periodic boundary has a vertical step there, which the C-method "
                "cannot represent."
            )
        t, y = t[:-1], y[:-1]
    s = np.mod(1.0 - t[::-1], 1.0)
    a = y[::-1] * problem.period
    order = np.argsort(s, kind="stable")
    return ProfileSeries(problem.period, "polygon", s[order], a[order])


def profile_series(problem: Problem, samples: int = _MIN_SAMPLES) -> ProfileSeries:
    """The profile of ``problem`` as a :class:`ProfileSeries`.

    Refuses, with :class:`UnsupportedConfiguration`, every profile the
    method cannot describe by a finite slope: vertical facets and sidewalls,
    and undercut boundaries (also caught earlier by ``Capabilities.check``).
    """
    profile = problem.profile
    if not profile.is_single_valued():
        raise UnsupportedConfiguration(
            "the C-method cannot represent an undercut profile: the boundary "
            "must be a single-valued height function of position"
        )
    if isinstance(profile, Lamellar):
        raise UnsupportedConfiguration(
            "Lamellar has vertical sidewalls; the C-method has no finite slope "
            "there and would have to smooth them. Use the integral method."
        )
    if isinstance(profile, Blazed):
        if profile.has_vertical_facet:
            raise UnsupportedConfiguration(
                "a Blazed profile with antiblaze_angle=90 has a vertical facet; "
                "the C-method has no finite slope there and would have to smooth "
                "it. Use the integral method, or a profile that is not vertical."
            )
        return _polygon(problem, [0.0, profile.apex], [0.0, profile.depth])
    if isinstance(profile, FromProfileData):
        return _polygon(problem, profile.t, profile.y)
    if isinstance(profile, Trapezoidal):
        return _trapezoidal(problem, profile)

    samples = max(int(samples), _MIN_SAMPLES)
    s = np.arange(samples) / samples
    t = np.mod(1.0 - s, 1.0)
    try:
        # d/dx of P h(1 - x / P) is -h'(t).
        slopes = -np.asarray(profile.slope(t), dtype=np.float64)
    except ProfileRepresentationError as exc:
        raise UnsupportedConfiguration(str(exc)) from exc
    a = np.asarray(profile.height(t), dtype=np.float64) * problem.period
    return ProfileSeries(problem.period, "sampled", s, a, slopes)


def _trapezoidal(problem: Problem, profile: Trapezoidal) -> ProfileSeries:
    """A trapezoid as a polygon, or as straight runs and fillet arcs.

    Until 2026-10-29 a Trapezoidal took the sampled path, whose FFT of a slope
    that jumps at four corners converges at first order in the sample count,
    so the answer sat on a floor the truncation knob could not move: 3.9e-4
    per order on a 72 deg, 0.425-period groove, with no corner warning
    (``findings.md``, "Trapezoidal through the four solvers").
    """
    if profile.has_vertical_sidewall:
        raise UnsupportedConfiguration(
            "Trapezoidal with sidewall_angle=90 has vertical sidewalls; the "
            "C-method has no finite slope there, and rounding the corners does "
            "not remove the vertical wall between them. Use the integral "
            "method or RCWA, or a slightly shallower sidewall_angle."
        )
    if profile.corner_radius == 0.0:
        corners = profile._corners()
        return _polygon(
            problem,
            [0.0, *(x for x, _ in corners), 1.0],
            [0.0, *(y for _, y in corners), 0.0],
        )

    # Pieces in t, left to right, each spanning [x0, x1]. A zero-width run
    # (a wall or flat the fillets consume entirely) drops out.
    pieces = [p for p in profile._pieces() if p.x1 > p.x0]
    # Reflected as _polygon reflects: the piece filling x0..x1 in t fills
    # 1 - x1 .. 1 - x0 in s, so its segment starts at the vertex s = 1 - x1,
    # at the height the profile has at x1, and ends where the next one starts.
    starts = np.mod(1.0 - np.array([p.x1 for p in pieces]), 1.0)
    heights = np.array([float(p.height(np.array([p.x1]))[0]) for p in pieces])
    order = np.argsort(starts, kind="stable")
    position = np.argsort(order)  # where each piece's vertex lands
    arcs = tuple(
        (int(position[k]), piece) for k, piece in enumerate(pieces) if piece.is_arc
    )
    return ProfileSeries(
        problem.period, "pieces", starts[order], heights[order] * problem.period,
        arcs=arcs,
    )
