# RCWA (Fourier modal method)

Implemented in [`solvers/rcwa/`](../../src/gratinglab/solvers/rcwa/__init__.py).
Symbols and sign conventions follow [`conventions.md`](../conventions.md). This
page is the design note the backend was written against (2026-10-27) and the
record of what has been validated so far; the section numbers are for citing.

The roadmap's requirements for this backend are not negotiable and each has a
section here: Li's factorization rules from the start (§3), S-matrix
propagation and never a plain transmittance matrix (§4), a fully vectorial
conical formulation (§2), and layer stacks above and below the grating region
in the data model from the start (§6).

---

## 1. What the model assumes

The grating is one-dimensional: the structure varies along `x` (the dispersion
direction `d̂`) and is invariant along the grooves. Along the normal it is a
stack of slabs, each of which is either homogeneous or piecewise constant in
`x`. Inside a slab the field is expanded in `2M + 1` Floquet harmonics,
`m = −M … M`, and Maxwell's equations become a matrix eigenproblem whose
eigenvectors are the slab's modes. The modes are matched across the slab
interfaces by a scattering matrix.

That is exact for a lamellar (binary) grating, whose single slab *is* the
geometry. Any other profile is replaced by a staircase of slabs (§5), which is
an approximation of the geometry that converges as the slabs are refined.

`truncation` is `M`, the accuracy knob, and `Capabilities.accuracy_knob` names
it, so `convergence.check_convergence` sweeps it. A truncation that cannot
represent every order the scan reports is refused (`ValueError`) rather than
truncating the orders.

### The problem it is given

`Problem` has no layer stack and no superstrate. The mapping today is therefore
the one the integral solver's tabulated mode uses, so that the same `Problem`
means the same grating to both:

| | |
|---|---|
| superstrate | vacuum, which is the medium the incidence angles are defined in |
| grating region | ridges of the named material, grooves filled with vacuum |
| substrate | the same material, semi-infinite, below the bottom of the grooves |
| material | `problem.coating`, falling back to `substrate` (as the integral solver does); naming neither is refused |

Efficiencies are **absolute** (conventions §5). The method needs a complex
index; there is no perfect-conductor mode, because a perfect conductor has no
permittivity to Fourier-transform. Comparisons with the integral solver use its
`conductivity="tabulated"` mode.

Where `coating` and `substrate` are both named and differ, the substrate is not
consulted and a provenance warning says so. A coated grating needs a coating
*thickness*, which `Problem` cannot carry; see §6 and the decisions at the end
of this page.

## 2. Formulation: vectorial and conical

Internally the frame is `(x, y, z) = (x, z_gl, −y_gl)`: `x` along the
dispersion direction as everywhere, `y` along the grooves, and `z` pointing
*down* into the grating, which is the usual RCWA orientation. It is
right-handed (`x̂ × ẑ_gl = −ŷ_gl`). The incident wave vector of
conventions §3, normalised by `k₀ = 2π/λ`, is

```
Kx0 = −sin α sin γ,   Ky = cos γ,   Kz0 = cos α sin γ  (> 0, downward)
```

and harmonic `m` carries `Kx_m = Kx0 + mλ/p`, which is the generalized grating
equation of conventions §4 with the same order labels. `Ky` is conserved.
Nothing in the formulation is specific to `γ = 90°`: the in-plane case is
`Ky = 0`, at which the equations below decouple into TE and TM by themselves,
and no branch of the code selects a different formulation for it.

With `U = η₀H` and `z̃ = k₀z`, Maxwell's curl equations for the transverse
components of one slab are, for `e^{−iωt}`,

```
d/dz̃ [Ex, Ey] = i P [Ux, Uy],      d/dz̃ [Ux, Uy] = i Q [Ex, Ey]

P = | Kx Ezz⁻¹ Ky       I − Kx Ezz⁻¹ Kx |      Q = | −Kx Ky        Kx² − Eyy |
    | Ky Ezz⁻¹ Ky − I   −Ky Ezz⁻¹ Kx    |          | Exx − Ky²      Ky Kx    |
```

where `Kx`, `Ky` are diagonal and `Exx`, `Eyy`, `Ezz` are the Fourier matrices
of the permittivity that §3 defines. Eliminating `U` gives
`d²E/dz̃² = −PQ E`, so a slab's modes are the eigenvectors `W` of `PQ` with
propagation constants `κ = √eig(PQ)`, the branch taken with `Im κ ≥ 0` (and
`Re κ > 0` on the real axis), so that a mode labelled forward decays or
propagates downward. The magnetic eigenvectors are `V = Q W κ⁻¹`. In a
homogeneous slab `PQ = (ε − Kx² − Ky²) I` and the modes are written down
directly rather than found by an eigensolver, which would otherwise have to
split a degenerate eigenspace.

The full `2N × 2N` problem (`N = 2M + 1`) is solved, not the two `N × N`
problems a lossless in-plane code would use. That is the price of conical
incidence and of TE–TM conversion, which off-plane is genuine physics
(conventions §7).

### Polarization, by the grooves

TE and TM are defined by the grooves (conventions §7), not by the plane of
incidence: TE is the incident wave whose magnetic field has no groove component
(`Uy = 0`), TM the one whose electric field has none (`Ey = 0`). In-plane these
are the usual s and p. Each order's efficiency is the total power it carries,
whatever its polarization, so a TE input that converts partly to TM off-plane
still reports one number per order. Unpolarized is `(TE + TM)/2`, from two
incident columns pushed through one S-matrix.

### Efficiency

Each order's efficiency is its Poynting flux normal to the grating,
`Re(Ex Uy* − Ey Ux*)`, divided by the incident flux. The flux is formed from the
transverse fields rather than from `|r|² Re κ`, so the same expression holds in
an absorbing substrate. Evanescent orders are set to exactly `0.0` with
`propagating = False` (conventions §4); their computed flux is zero to rounding
in a lossless superstrate, and the convention makes it exact.

## 3. Li's factorization rules

The permittivity of a slab, `ε(x)`, jumps at every ridge wall. The Fourier
coefficients of a product `εE` are correctly the Toeplitz matrix of `ε` times
those of `E` (Laurent's rule) only when `E` is continuous where `ε` jumps, and
the Toeplitz matrix of `1/ε`, *inverted*, when it is `εE` that is continuous and
`E` that jumps (the inverse rule). Li (1996) proved this, and it is what fixed
the notoriously slow convergence of TM gratings, metallic ones above all.

At a wall `x = const` the tangential components are `Ey` and `Ez`, and the
normal component of `D` is `Dx`. So:

| Fourier matrix | Component | Rule |
|---|---|---|
| `Exx` | `Dx = ε Ex`, `Dx` continuous, `Ex` jumps | inverse: `Exx = [[1/ε]]⁻¹` |
| `Eyy` | `Dy = ε Ey`, `Ey` continuous | Laurent: `Eyy = [[ε]]` |
| `Ezz` | `Dz = ε Ez`, `Ez` continuous | Laurent: `Ezz = [[ε]]`, inverted where `P` needs `Ezz⁻¹` |

That assignment is the conical one (Li 1997; Popov & Nevière 2000), and at
`Ky = 0` it reduces to the TM formulation of Lalanne & Morris and of Granet &
Guizal (1996), since TM then involves only `Exx` and `Ezz⁻¹`. In a homogeneous
slab all three are `ε I`.

It is used from the first commit and is not optional in `solve()`. The core
keeps the naive all-Laurent assignment behind a private argument for one
purpose: `tests/test_rcwa.py` shows that metallic TM converges with Li's rule
and visibly does not without it, which is the only evidence that the rule is
really in effect. A code that silently used Laurent's rule everywhere would
pass every energy-balance test.

## 4. Propagation: a scattering matrix

A plain transmittance matrix multiplies `e^{±iκL}` together and overflows as
soon as a slab carries an evanescent mode of any size, which a truncated
expansion always does. The scattering matrix keeps only the decaying
exponential, `X = e^{iκL}` with `|X| ≤ 1`, and composes slabs with the
Redheffer star product, so it is stable in slab thickness.

The formulation is Rumpf's symmetric S-matrix in a zero-thickness gap medium
(Rumpf, PIER B 35, 241, 2011). The gap is chosen mode by mode with `κ = 1`, so
its modes are `W₀ = I` and `V₀ = Q₀`, which is invertible at every `Kx` and so
has no Rayleigh-anomaly singularity of its own. The external regions (vacuum
above, substrate below) are homogeneous with analytic modes. The global matrix
is

```
S = S_ref ⊗ S_slab₁ ⊗ … ⊗ S_slabₙ ⊗ S_trn
```

and the reflected and transmitted amplitudes are `S₁₁` and `S₂₁` applied to the
incident amplitude. An external region has `V = Q κ⁻¹`, which is singular where
an order passes off (`κ = 0`, a Rayleigh anomaly). At that measure-zero point
`κ` is moved off zero by `10⁻¹²`, and any row with an order inside
`|cos β| < 10⁻³` is flagged in provenance, the guard the integral solver uses.

### Absorption, and what `R + A = 1` does and does not show

The scan's `absorption` is the power crossing into the structure through the
plane just below the vacuum, evaluated from the field inside the stack: the
amplitudes at that plane come from the partial S-matrices above and below it,
not from the reflected amplitudes. For a lossless substrate it is the
transmitted power; for an absorbing one, the power absorbed in the ridges plus
the power entering the substrate. `check_energy_balance(lossless=True)` then
checks `R + A = 1`.

That check is weaker here than for the integral solver, and the difference
matters. There, the absorption is an independent boundary integral and `R + A`
deviates from 1 by the discretisation error. In RCWA the truncated system
conserves flux by construction, because tangential `E` and `H` are continuous
across every interface the S-matrix joins, so `R + A = 1` holds to rounding at
*every* truncation, converged or not. It is a test of the implementation: of
signs, normalisation and the S-matrix composition. It is not an error estimate.
Convergence is shown by the truncation sweep alone, and the provenance notes say
so.

## 5. Profiles: exact slabs, and the staircase

| Profile | Slabs |
|---|---|
| `Lamellar` | one slab, exact: the ridge spans `t ∈ [0, duty_cycle)` |
| `Blazed`, `Sinusoidal`, single-valued `FromProfileData` | a staircase of `slices` equal-thickness slabs (below) |
| undercut `FromProfileData` | refused: `Capabilities.handles_undercut = False` |

The staircase uses the mid-height rule of `Profile.slice_layers`: a slab holds
material wherever the profile reaches its mid-height. The crossings are refined
by bisection on `Profile.height` to about `10⁻¹³` of the period rather than
taken from the sampling `slice_layers` uses to find them, because a sampling
error in a wall position is a floor no truncation can get under.

The profile's `t` runs along `−d̂` (conventions §3), so `t` maps to
`x = (1 − t) p`, as the integral solver places its boundary
(`solvers/integral/_boundary.py`). An asymmetric groove is therefore the same
groove to both solvers.

`slices` defaults to `max(8, 2 · truncation)`, so the one knob refines the
Fourier expansion and the staircase together and the convergence harness sees
both. That mirrors the integral solver's `spectral_terms`, which follows
`boundary_points` by default for the same reason. Passing `slices` explicitly
pins it, and a sweep then measures the Fourier convergence alone.

**The staircase is a known weakness of RCWA in TM, and it is not hidden.** Each
step has two right-angle corners, and for a metal in TM those corners carry
field singularities the smooth profile does not have (Popov, Nevière, Gralak &
Tayeb, JOSA A 19, 33, 2002). Convergence in `slices` is then slow and can be
non-monotone. Every staircased solve records the slice count in provenance, and
a staircased TM or unpolarized solve on an absorbing material carries a warning
naming the effect. **Staircase accuracy is not yet validated**; see §9.

## 6. Layer stacks

The core solves a `Stack`: a superstrate index, homogeneous layers *above* the
grating region, the grating slabs, homogeneous layers *below* it, and a
substrate index. Every slab, homogeneous or not, is a background index plus a
list of `x` segments each with its own index, so a conformal coating on a
staircased groove (ridge, coating and vacuum in one slab) is a slab with more
segments and not a new code path. The S-matrix composition is the same loop
whatever the stack holds, so a multilayer costs one more star product per layer
and no new physics.

What exists today, stated plainly:

- **In the core and tested**: homogeneous layers above and below the grating
  region, against the closed-form Airy formula for a single film, including a
  film beneath a lamellar slab whose energy balance is checked (§8).
- **Not reachable from `solve()`**, because `Problem` has no field for a layer:
  its `coating` is a material name without a thickness. A `Problem` therefore
  always becomes vacuum / grating / substrate with no layers. How to describe a
  stack on `Problem` is an interface decision for the whole package, since every
  solver must then either model the layers or refuse them, so it is listed in §9
  rather than settled inside the backend.

## 7. What is refused, rather than approximated

| Request | Response |
|---|---|
| an undercut profile | `UnsupportedConfiguration` (`Capabilities.check`): horizontal slabs cannot represent a boundary that doubles back, and a staircase would silently fill the overhang |
| no material named | `UnsupportedConfiguration`: there is no permittivity to expand |
| a material with gain (`Im n < 0`) | `ValueError`, as in the integral solver |
| `truncation` below the highest order the scan reports | `ValueError`, naming the order |
| a perfect conductor | not offered; use the integral solver |

Two things are *not* refused, and are recorded instead: a staircased profile
(§5, warned for an absorbing material in TM) and a Rayleigh anomaly (§4, flagged
per row).

## 8. Validation evidence

Everything here is in `tests/test_rcwa.py` and runs without the private corpus.

| Reference | What it pins | Independent of this code? |
|---|---|---|
| flat interface against `materials.fresnel`, in-plane, TE and TM, glass and a visible metal | signs and normalisation in the external regions | yes |
| flat interface, conical, unpolarized against `(Rs + Rp)/2` at the true angle of incidence | the conical kinematics, and that TE/TM by the grooves span the polarizations | yes |
| a single film on a substrate against the Airy formula, s and p, in-plane and conical | the layer stack | yes |
| lossless dielectric lamellar, `R + A = 1` with `A` the transmitted power, in-plane and conical | the S-matrix composition and flux bookkeeping (§4: not convergence) | partly |
| reciprocity, `checks.check_reciprocity`, metallic lamellar in TE and TM | the symmetry of the scattering problem | yes |
| truncation sweep, metallic lamellar in TM: Li's rule converges, Laurent's visibly does not | that §3 is in effect | yes |
| the integral solver (`conductivity="tabulated"`, corner-graded) on lamellar glass and metal, TE and TM, in-plane and conical | the whole method, against a different method | yes |

The measured numbers are in the test docstrings, beside the assertions they
justify.

No published benchmark table is committed yet. Moharam & Gaylord's and Li's
cases (`references.md`, Validation) are the obvious next ones, and their values
want transcribing from the papers rather than from memory.

## 9. Not yet validated, and decisions to make

Not yet validated:

- staircase accuracy for blazed and sinusoidal profiles, in TM above all (§5);
- multilayers beyond one film, and anything reachable through `solve()` with
  layers, since nothing is (§6);
- the soft-X-ray grazing regime. It is not refused, and it is expensive: a
  groove hundreds of wavelengths deep needs very many slices and harmonics,
  which is why the integral method is the method of record there;
- a published benchmark table (§8).

Decisions for the project rather than for this backend:

1. **How `Problem` describes a layer stack**, and what a coated grating means:
   a conformal coating of given thickness on the profile, a planar overcoat, or
   both. The core can represent all of these (§6); `Problem` can represent none.
2. **Whether a staircase is acceptable by default**, with a warning, or should
   need `slices` passed explicitly. The default here follows the integral
   solver's convention of one knob refining everything.
3. **Whether transmitted orders belong in `EfficiencyScan`.** A transmission
   grating's diffracted orders are the result a user wants; today their total
   appears only as `absorption` (§4).
