"""The physical problem: what the grating *is*, independent of how it is solved.

``Problem`` carries **no solver-specific fields, ever**. That single constraint
is what makes the comparison harness possible -- every method receives the same
object, so a disagreement between methods is a disagreement about physics rather
than about how two codes were configured.

It is also serializable, so a benchmark case is a file rather than a script.
"""

from __future__ import annotations

from typing import Union

from pydantic import BaseModel, ConfigDict, Field

from .profiles import Blazed, FromProfileData, Lamellar, Sinusoidal

__all__ = ["Problem", "AnyProfile"]

#: Profiles that can round-trip through JSON. A ``Profile`` structurally typed
#: duck can be passed to a solver directly, but only these serialise.
AnyProfile = Union[Blazed, Lamellar, Sinusoidal, FromProfileData]


class Problem(BaseModel):
    """A grating, in physical terms.

    >>> from gratinglab.profiles import Blazed
    >>> p = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
    >>> round(p.depth, 2)
    92.38

    Attributes
    ----------
    period
        Groove spacing in nm (``docs/conventions.md`` §1).
    profile
        Groove shape, normalised to the period. Multiply by ``period`` for
        physical height, or use :attr:`depth`.
    coating, substrate
        Material identifiers. Optional for now -- scalar theory needs no optical
        constants, and the perfect-conductivity reference data it is first
        compared against needs none either. The materials layer lands next.
    roughness
        RMS surface roughness in nm, for the Névot–Croce / Debye–Waller factors.
    n_grooves
        Illuminated groove count, for the finite-N interference factor of
        ISSI eq. (8). ``None`` means the infinite-grating limit, where orders
        become delta functions.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    period: float = Field(gt=0.0, description="Groove spacing, nm")
    profile: AnyProfile = Field(discriminator=None)
    coating: str | None = None
    substrate: str | None = None
    roughness: float = Field(default=0.0, ge=0.0, description="RMS roughness, nm")
    n_grooves: int | None = Field(default=None, gt=0)

    @property
    def depth(self) -> float:
        """Peak-to-valley groove depth in nm."""
        return self.profile.depth * self.period

    def height_nm(self, t) -> "object":
        """Profile height in nm at normalised position ``t``.

        The one place the normalised profile is converted to physical units.
        Solvers should call this rather than multiplying by hand.
        """
        return self.profile.height(t) * self.period

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"Problem(period={self.period:g}nm, "
            f"{type(self.profile).__name__}, depth={self.depth:.3g}nm)"
        )
