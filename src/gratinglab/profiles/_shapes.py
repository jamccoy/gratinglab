""":class:`Blazed`, :class:`Lamellar`, :class:`Sinusoidal` and
:class:`FromProfileData`: every profile type except the trapezoid, whose corner
arithmetic has its own file.

Frozen pydantic models, built by callers and read by the solvers.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pydantic import Field, model_validator

from ._base import _BaseProfile
from ._polyline import _polyline_boundary
from ._types import BoundaryCurve, ProfileRepresentationError


class Blazed(_BaseProfile):
    """Sawtooth groove, the workhorse of blazed reflection gratings.

    ``depth = 1 / (cot(blaze) + cot(antiblaze))`` in units of the period, which
    reduces to ``tan(blaze)`` for the ideal sawtooth with a vertical anti-blaze
    facet (``docs/conventions.md`` 9).

    >>> round(Blazed(blaze_angle=30.0, antiblaze_angle=80.0).depth, 4)
    0.5241
    """

    blaze_angle: float = Field(gt=0.0, lt=90.0, description="Active facet angle, deg")
    antiblaze_angle: float = Field(
        default=90.0, gt=0.0, le=90.0, description="Opposite facet angle, deg"
    )

    @property
    def _runs(self) -> tuple[float, float]:
        """Horizontal extent of each facet, normalised. Sums to 1."""
        cot_b = 1.0 / np.tan(np.radians(self.blaze_angle))
        cot_a = 1.0 / np.tan(np.radians(self.antiblaze_angle))
        total = cot_b + cot_a
        return cot_b / total, cot_a / total

    @property
    def depth(self) -> float:
        cot_b = 1.0 / np.tan(np.radians(self.blaze_angle))
        cot_a = 1.0 / np.tan(np.radians(self.antiblaze_angle))
        return 1.0 / (cot_b + cot_a)

    @property
    def apex(self) -> float:
        """Normalised position of the facet crest."""
        return self._runs[0]

    @property
    def has_vertical_facet(self) -> bool:
        return bool(np.isclose(self.antiblaze_angle, 90.0))

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        apex, depth = self.apex, self.depth
        rising = depth * t / apex
        if self.has_vertical_facet:
            return np.where(t < 1.0, rising, 0.0)
        falling = depth * (1.0 - t) / (1.0 - apex)
        return np.where(t <= apex, rising, falling)

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        if self.has_vertical_facet:
            raise ProfileRepresentationError(
                f"{type(self).__name__} with antiblaze_angle=90 has a vertical "
                "facet and no finite slope there. The C-method cannot represent "
                "it; use a slightly shallower antiblaze_angle, or a method that "
                "handles discontinuities (RCWA, integral)."
            )
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        up = np.tan(np.radians(self.blaze_angle))
        down = -np.tan(np.radians(self.antiblaze_angle))
        return np.where(t <= self.apex, up, down)

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Exact sawtooth polygon.

        A vertical anti-blaze facet is the closing segment of the polygon, a
        piece of geometry that a height function cannot express.
        """
        if self.has_vertical_facet:
            x, y = np.array([0.0, 1.0]), np.array([0.0, self.depth])
        else:
            x, y = np.array([0.0, self.apex]), np.array([0.0, self.depth])
        return _polyline_boundary(x, y, n, grading=grading)


class Lamellar(_BaseProfile):
    """Rectangular (binary, laminar) groove."""

    depth_fraction: float = Field(gt=0.0, description="Depth / period")
    duty_cycle: float = Field(
        default=0.5, gt=0.0, lt=1.0, description="Fraction of the period at full height"
    )

    @property
    def depth(self) -> float:
        return self.depth_fraction

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        return np.where(t < self.duty_cycle, self.depth, 0.0)

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        raise ProfileRepresentationError(
            "Lamellar has vertical sidewalls and no finite slope. The C-method "
            "cannot represent it; use RCWA or the integral method."
        )

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Exact rectangular polygon; both sidewalls are genuine vertical segments."""
        d, w = self.depth, self.duty_cycle
        x = np.array([0.0, w, w, 1.0])
        y = np.array([d, d, 0.0, 0.0])
        return _polyline_boundary(x, y, n, grading=grading)


class Sinusoidal(_BaseProfile):
    """Holographic (sinusoidal) groove, smooth and therefore suited to the C-method."""

    depth_fraction: float = Field(gt=0.0, description="Peak-to-valley / period")

    @property
    def depth(self) -> float:
        return self.depth_fraction

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.asarray(t, dtype=np.float64)
        return 0.5 * self.depth * (1.0 - np.cos(2.0 * np.pi * t))

    def slope(self, t: ArrayLike) -> NDArray[np.float64]:
        t = np.asarray(t, dtype=np.float64)
        return np.pi * self.depth * np.sin(2.0 * np.pi * t)


class FromProfileData(_BaseProfile):
    """A measured profile, from AFM or from PCGrate's ``.ggp`` polygonal border.

    Points are normalised: ``t`` spans one period on ``[0, 1]`` and ``y`` is in
    units of the period. Use :meth:`from_measurements` to normalise raw nm data.

    The profile is stored as given and interpolated linearly, matching how
    PCGrate treats a polygonal border, so a scalar calculation and a PCGrate
    run can be driven from the identical geometry. One solver departs from
    this when asked: the integral method under ``conductivity="tabulated"``
    can smooth the polygon's slight vertices, which wreck its absorption in
    the X-ray regime (issue #6; its ``profile_smoothing`` option, off by
    default, and ``docs/theory/integral.md`` section 8).

    Undercut profiles are representable. If ``t`` doubles back, the boundary
    is not a function of position and :meth:`height` and
    :meth:`~_BaseProfile.slice_layers` refuse, but :meth:`boundary` still
    works, because an integral-equation solver parametrises the curve instead
    of assuming a height function. This asymmetry reflects the physics of the
    methods and is not a limitation of this class: RCWA and scalar theory
    cannot represent an undercut, and the integral method can.

    ``blaze_angle`` is optional and has a different meaning here from the one
    it has on :class:`Blazed`. There it is a definition: the facet has that
    angle because it was specified so. Here it is a measurement, obtained by
    fitting a line to a real facet, and it carries an uncertainty. Both feed
    ``solvers._scalar_groove._reflecting_graze`` identically, and the provenance string
    records which kind was used, because a reader deciding whether to trust an
    absolute efficiency needs to know.

    Leaving it ``None`` has consequences. The reflection then falls back to
    the mean surface, a zero-tilt facet, which for a sawtooth measured at
    grazing incidence is the wrong plane by the whole blaze angle. The default
    is ``None`` because most measured profiles arrive without a fit, and not
    because ``None`` is the better answer.
    """

    t: tuple[float, ...]
    y: tuple[float, ...]
    blaze_angle: float | None = None

    @model_validator(mode="after")
    def _check_points(self) -> "FromProfileData":
        if len(self.t) != len(self.y):
            raise ValueError(f"t has {len(self.t)} points, y has {len(self.y)}")
        if len(self.t) < 3:
            raise ValueError(f"need at least 3 points, got {len(self.t)}")
        t = np.asarray(self.t)
        if t.min() < 0.0 or t.max() > 1.0:
            raise ValueError(f"t must lie in [0, 1], got [{t.min()}, {t.max()}]")
        if self.blaze_angle is not None and not 0.0 < self.blaze_angle < 90.0:
            # The same open interval `Blazed` enforces. A fit that lands outside
            # it is not a shallow facet, it is a fit that found the wrong facet
            # or the wrong sign convention, and passing it on as geometry would
            # produce a confident, wrong reflectivity.
            raise ValueError(
                f"blaze_angle must lie in (0, 90) degrees, got {self.blaze_angle}"
            )
        return self

    @classmethod
    def from_measurements(
        cls,
        x_nm: ArrayLike,
        y_nm: ArrayLike,
        period_nm: float,
        *,
        blaze_angle: float | None = None,
    ) -> "FromProfileData":
        """Normalise a measured profile in nm, referencing its minimum to zero.

        Pass ``blaze_angle`` when the same measurement also produced a facet
        fit; see the class docstring for why it is worth carrying.
        """
        x = np.asarray(x_nm, dtype=np.float64)
        y = np.asarray(y_nm, dtype=np.float64)
        if period_nm <= 0:
            raise ValueError(f"period_nm must be positive, got {period_nm}")
        return cls(
            t=tuple((x - x.min()) / period_nm),
            y=tuple((y - y.min()) / period_nm),
            blaze_angle=blaze_angle,
        )

    @classmethod
    def along_dispersion(
        cls,
        x: ArrayLike,
        y: ArrayLike,
        *,
        blaze_angle: float | None = None,
    ) -> "FromProfileData":
        """A profile from points that advance along ``+d-hat``.

        PCGrate's ``.ggp`` and the metrology pipeline's ``x_norm`` both run
        along the dispersion direction, while ``t`` runs against it
        (``docs/conventions.md``, "Direction of the profile parameter"). This
        is the one place that difference is absorbed: ``t = 1 - x``, with the
        point order reversed so ``t`` still advances. :meth:`dispersion_axis`
        is the inverse.
        """
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return cls(
            t=tuple(1.0 - x[::-1]), y=tuple(y[::-1]), blaze_angle=blaze_angle
        )

    def dispersion_axis(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """The points as ``(x, y)`` advancing along ``+d-hat``, as a ``.ggp`` holds them."""
        t = np.asarray(self.t, dtype=np.float64)
        return 1.0 - t[::-1], np.asarray(self.y, dtype=np.float64)[::-1]

    @property
    def depth(self) -> float:
        return float(np.ptp(np.asarray(self.y)))

    def height(self, t: ArrayLike) -> NDArray[np.float64]:
        if not self.is_single_valued():
            raise ProfileRepresentationError(
                "this profile is undercut -- its boundary doubles back in t, so "
                "height is not a function of position. Use boundary() with an "
                "integral-equation solver; scalar theory and RCWA cannot "
                "represent an undercut."
            )
        query = np.mod(np.asarray(t, dtype=np.float64), 1.0)
        return np.interp(query, np.asarray(self.t), np.asarray(self.y), period=1.0)

    def is_single_valued(self) -> bool:
        """A measured boundary that doubles back in ``t`` is undercut."""
        return bool((np.diff(np.asarray(self.t)) >= 0).all())

    @property
    def apex(self) -> float:
        """Normalised position of the highest measured point."""
        return float(np.asarray(self.t)[int(np.argmax(np.asarray(self.y)))])

    def boundary(self, n: int, *, grading: int = 0) -> BoundaryCurve:
        """Resample the stored polygon at equal arc length, or graded.

        Works from the stored points directly instead of through ``height``,
        which allows it to handle an undercut boundary.
        """
        return _polyline_boundary(
            np.asarray(self.t, dtype=np.float64),
            np.asarray(self.y, dtype=np.float64),
            n,
            grading=grading,
        )
