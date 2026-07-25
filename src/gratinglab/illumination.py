"""Incident-wave geometry, stored once and expressible in every community's angles.

The canonical internal state is the pair ``(alpha_deg, gamma_deg)`` -- azimuthal
incidence angle and half-cone angle -- because every mount reduces to it and it
serialises readably. Direction cosines are available as a derived property.

Three constructors exist so that each community can pass the angles it actually
uses; they all resolve to the same internal state:

>>> Illumination.classical(alpha=10.0)                  # in-plane
>>> Illumination.conical(theta=10.0, phi=45.0)          # RCWA-style polar/azimuth
>>> Illumination.offplane(graze=1.5, azimuth=0.97)      # extreme off-plane, X-ray

See ``docs/conventions.md`` §3 for the frame and the sign of ``alpha``.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator

__all__ = ["Illumination", "Polarization"]

Polarization = Literal["TE", "TM", "unpolarized"]


class Illumination(BaseModel):
    """Direction and polarization of the incident plane wave.

    Wavelength is deliberately **not** a field: a single geometry is normally
    scanned over many wavelengths, so it is passed at solve time instead.

    Attributes
    ----------
    alpha_deg
        Azimuthal angle of incidence, measured from the grating normal in the
        plane containing the normal and the dispersion direction. Positive
        ``alpha`` tilts the incident wave toward :math:`-\\hat{d}`, so the
        specular order sits at :math:`\\beta_0 = -\\alpha`.
    gamma_deg
        Half-angle of the diffraction cone about the groove axis.
        ``gamma_deg = 90`` is the in-plane (classical) case; extreme off-plane
        mounts use small values (soft X-ray work runs ~1.5 deg).
    polarization
        ``TE`` has **E** parallel to the grooves, ``TM`` has **H** parallel to
        the grooves. In conical mounts these are referenced to the grooves, not
        to the plane of incidence (``docs/conventions.md`` §7).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    alpha_deg: float = Field(description="Azimuthal incidence angle, degrees")
    gamma_deg: float = Field(description="Half-cone angle, degrees; 90 is in-plane")
    polarization: Polarization = "TE"

    @field_validator("alpha_deg")
    @classmethod
    def _check_alpha(cls, v: float) -> float:
        if not -90.0 < v < 90.0:
            raise ValueError(
                f"alpha_deg must lie strictly within (-90, 90), got {v}. "
                "An incidence angle at or past grazing carries no power onto "
                "the grating."
            )
        return v

    @field_validator("gamma_deg")
    @classmethod
    def _check_gamma(cls, v: float) -> float:
        if not 0.0 < v <= 90.0:
            raise ValueError(
                f"gamma_deg must lie in (0, 90], got {v}. "
                "gamma = 90 is the in-plane case; gamma = 0 would send the wave "
                "along the grooves, which diffracts nothing."
            )
        return v

    # -- constructors ----------------------------------------------------

    @classmethod
    def classical(
        cls, alpha: float, polarization: Polarization = "TE"
    ) -> "Illumination":
        """In-plane mount: the plane of incidence is perpendicular to the grooves.

        Equivalent to ``gamma = 90``, where the diffraction cone opens into a
        full 2D arc.
        """
        return cls(alpha_deg=alpha, gamma_deg=90.0, polarization=polarization)

    @classmethod
    def conical(
        cls, theta: float, phi: float, polarization: Polarization = "TE"
    ) -> "Illumination":
        """General conical mount, in the polar/azimuth angles common to RCWA codes.

        Parameters
        ----------
        theta
            Polar angle from the grating normal, degrees.
        phi
            Azimuth of the plane of incidence, degrees, measured from the
            :math:`-\\hat{d}` (anti-dispersion) direction. ``phi = 0`` therefore
            reduces exactly to ``classical(alpha=theta)``, and ``phi = 90``
            puts the plane of incidence along the grooves.
        """
        th, ph = np.radians(theta), np.radians(phi)
        cos_gamma = np.sin(th) * np.sin(ph)
        gamma = np.arccos(np.clip(cos_gamma, -1.0, 1.0))
        alpha = np.arctan2(np.sin(th) * np.cos(ph), np.cos(th))
        return cls(
            alpha_deg=float(np.degrees(alpha)),
            gamma_deg=float(np.degrees(gamma)),
            polarization=polarization,
        )

    @classmethod
    def offplane(
        cls, graze: float, azimuth: float, polarization: Polarization = "TE"
    ) -> "Illumination":
        """Extreme off-plane mount, in the angles used by the X-ray community.

        Parameters
        ----------
        graze
            Half-cone opening angle :math:`\\gamma`, degrees -- the graze angle
            of the incident wave onto the grating surface.
        azimuth
            Azimuthal angle :math:`\\alpha`, degrees. Near-Littrow blazed
            operation has ``azimuth`` close to the facet angle.
        """
        return cls(alpha_deg=azimuth, gamma_deg=graze, polarization=polarization)

    @classmethod
    def from_direction_cosines(
        cls, ux: float, uy: float, uz: float, polarization: Polarization = "TE"
    ) -> "Illumination":
        """Build from an incident direction vector (need not be normalised).

        The vector must have ``uy < 0`` -- it travels *toward* the grating.
        """
        u = np.asarray([ux, uy, uz], dtype=float)
        norm = np.linalg.norm(u)
        if norm == 0:
            raise ValueError("direction vector must be non-zero")
        ux, uy, uz = u / norm
        if uy >= 0:
            raise ValueError(
                f"incident direction must travel toward the grating (uy < 0), got uy={uy}"
            )
        gamma = np.arccos(np.clip(uz, -1.0, 1.0))
        alpha = np.arctan2(-ux, -uy)
        return cls(
            alpha_deg=float(np.degrees(alpha)),
            gamma_deg=float(np.degrees(gamma)),
            polarization=polarization,
        )

    # -- derived quantities ----------------------------------------------

    @property
    def alpha(self) -> float:
        """Azimuthal incidence angle in radians."""
        return float(np.radians(self.alpha_deg))

    @property
    def gamma(self) -> float:
        """Half-cone angle in radians."""
        return float(np.radians(self.gamma_deg))

    @property
    def sin_alpha(self) -> float:
        return float(np.sin(self.alpha))

    @property
    def cos_alpha(self) -> float:
        return float(np.cos(self.alpha))

    @property
    def sin_gamma(self) -> float:
        return float(np.sin(self.gamma))

    @property
    def is_in_plane(self) -> bool:
        """True for a classical mount, to floating-point tolerance."""
        return bool(np.isclose(self.gamma_deg, 90.0))

    @property
    def direction_cosines(self) -> np.ndarray:
        r"""Unit vector along :math:`\mathbf{k}_i`, in the ``(x, y, z)`` frame.

        .. math::
            \hat{u}_i = \left[-\sin\alpha\sin\gamma,\;
                              -\cos\alpha\sin\gamma,\;
                              \cos\gamma\right]

        The ``y`` component is negative: the incident wave travels toward the
        grating (``docs/conventions.md`` §6).
        """
        sg = self.sin_gamma
        return np.array(
            [-self.sin_alpha * sg, -self.cos_alpha * sg, float(np.cos(self.gamma))]
        )

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        mount = "in-plane" if self.is_in_plane else f"gamma={self.gamma_deg:g}deg"
        return (
            f"Illumination(alpha={self.alpha_deg:g}deg, {mount}, "
            f"{self.polarization})"
        )
