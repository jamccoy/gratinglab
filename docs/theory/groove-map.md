# Groove maps: a design draft

**Status: proposal, 2026-11-20. Nothing here is implemented.** It proposes the
first content for the design description in [`roadmap.md`](../roadmap.md)
("A design description, and where the fabrication boundary falls"): curved
and variable-line-space (VLS) grooves on curved substrates, from Fermat's
principle and holographic recording. The worked case is the author's slides,
*On Determining a Groove Map from Holographic Aberration Theory for a
Back-Diffracting Convex Cylinder* (2022), which follow Noda, Namioka & Seya,
*Geometric theory of the grating*, JOSA 64, 1031 (1974).

## 1. What it is for, and where it stops

The roadmap's design description is meant to be one shared description of
the grating someone intends to make. That description would feed design,
predicted efficiency, the written file (`io/txl.py`), metrology and measured
efficiency. A groove map is the part that exists before any of the rest: where
each groove goes on the substrate.

The boundary follows the roadmap's: **the grating, not the instrument.** This
package would describe the groove map, its local geometry, its aberrations
and its efficiency. Tracing a whole instrument stays with PyXFocus, and
exposure settings stay with `ebl_fracture`.

## 2. The object: a groove function on a surface

In Noda, Namioka & Seya's terms, a grating is a **surface** together with a
**groove function** n(P), where groove n passes through the points with
n(P) = n.

- **Surface.** The substrate as ξ(w, l), with w across the grooves at the
  centre O and l along them. A plane, a cylinder and a sphere are exact
  constructors. The general a_ij series of the slides' Eq. 1 is a fallback for
  shapes with no closed form.
- **Groove function.** It is supplied by a constructor:
  - a plane-wave or point-source **hologram**: recording sources C and D and
    the recording wavelength λ₀, with n λ₀ = (CP − DP) − (CO − DO) (slides,
    Eq. 7);
  - a **VLS polynomial** in w;
  - **radial** grooves (ISSI Eq. 17).
- **Local grating vector.** K(P) = 2π ∇_S n, the surface gradient. It gives the
  local period and the local groove direction at every point, which are what
  the solvers and the vector grating equation (ISSI Eq. 1) consume.

## 3. Two routes, and why both

- **Exact, numerical.** Groove positions are the roots of n(P) = integer on
  the exact surface, found by a bracketing root finder one groove at a time.
  Rays are traced through K(P) with the vector grating equation. Nothing is
  truncated, and no root has to be chosen by eye.
- **Series, analytic.** The light-path function F = AP + PB + nmλ expanded in
  F_ij (slides, Eq. 14), split into mechanical and holographic parts M_ij and
  H_ij (Eq. 16). This is what names the aberrations and shows how a recording
  geometry cancels them.

The exact route produces the map. The series explains it, and its closed
forms become tests of the exact code, as `sinc²` tests the scalar solver's
quadrature. The two routes share almost nothing, so agreement between them is
evidence. That is the lesson of the validation inventory: solvers that share
plumbing cannot check each other's shared parts.

## 4. The worked case, measured

The slides' case is plane-wave back-diffraction from a convex cylinder,
Rc = 155.1 mm, λ = 632.8 nm, m = −1, groove spacing d_O = 1000 nm at O, so
α = arcsin(mλ/2d_O) = −18.445°. The slides solve a quartic by Ferrari's method
and compare it with DeRoo's ray trace.

**The recording geometry has an exact solution, and it is simple.** The
recording angles of the slides' Eq. 23 make C and D counter-propagating plane
waves along the incident beam. Their fringes are planes perpendicular to the
beam, spaced |m|λ/2. The grooves are where the substrate crosses those planes.
For back-diffraction in Littrow, each groove then sends the light back
exactly: the light path F is the same through every point, so every F_i0
vanishes, not only F_20 and F_40. This holds **for any substrate shape**.

On the exact circle, the slides' groove function n d_O = w + cot α · ξ(w)
(Eq. 24, before truncation) becomes a **quadratic** in w once the square root
is isolated:

```
(1 + c²) w² − 2K w + K² − c² Rc² = 0,   c = cot α,   K = n d_O + c Rc
```

with the root through w₀ = 0. The slides' quartic appears only because the
circle was truncated at w⁴.

Measured 2026-11-20 with 40-digit arithmetic. The table gives the groove
spacing, in linear w, from the exact quadratic, and the errors of the slides'
two approximations against it:

| w (mm) | exact spacing (nm) | quartic − exact | quadratic (Eq. 29) − exact |
|---:|---:|---:|---:|
| −35 | 3272.98 | −10 nm | −270 nm |
| −22 | 1753.19 | −0.23 nm | −16 nm |
| −2.5 | 1050.80 | −1.4 fm | −7.0 pm |
| 0 | 1000.00 | 0 | 0 |
| +2.5 | 953.90 | +1.1 fm | +5.6 pm |
| +22 | 699.49 | +0.030 nm | +2.0 nm |
| +35 | 590.15 | +0.22 nm | +5.7 nm |

This reproduces the slides' own findings: a few picometres for the quadratic
at |w| ≤ 2.5 mm, a fraction of a nanometre for the quartic to 22 mm, and an
error that grows faster toward −22 mm than toward +22 mm. So DeRoo's ray
trace behaves like this exact solution.

As an independent check, every groove satisfies the local Littrow condition
2|sin α_local| · d_arc = |m|λ. Here α_local is α turned by the local surface
tilt and d_arc is the spacing along the surface. It holds to between 1e-9 and
1e-12 from −35 to +35 mm, which is the finite-difference error of the check
itself.

**What this means for the design.** The convex cylinder in Littrow is a test
case with a closed form, not a reason to need Ferrari. The general cases do
need the numerical route: finite source distances, a focusing spectrograph
where B differs from A, recording not along the beam, sources out of plane,
and a sphere. The slides' Eq. 13 series and quartic remain useful as tests of
the truncation regime.

The script behind this table is not committed yet. It uses mpmath, which is
not a dependency of the package, so it would need a guard or a float64 version
first.

## 5. Conventions to write down first

- **Linear position or arc length.** The slides found an ångström-scale
  systematic error from mixing them. EBL wants linear positions, and a ray
  trace works in arc length. The design should store one, convert in one
  function, and name which one at every interface. This belongs in
  `conventions.md`.
- **The w axis against d̂ and t̂.** `conventions.md` §3 records that the profile
  parameter runs along −d̂. The groove index n, the w axis and d̂ need the same
  explicit mapping.
- **Signs of m and α.** The slides' Eq. 25 is the package's
  `geometry.littrow_alpha` in plane, in the same sin α + sin β = mλ/d form as
  `conventions.md` §4. The code should call it, not restate it, so it inherits
  the conical check in #87.

## 6. How it would be checked

- **Closed forms.** The plane-wave Littrow solution above, exact for any
  surface. The Rowland circle and the Seya–Namioka mount. Holographic
  recordings that are stigmatic for chosen points, where named aberrations
  vanish exactly. Published F_ij values for a standard concave grating.
- **The slides' truncations.** The exact solver must reproduce the quadratic
  and quartic forms where their truncation holds, to the levels in §4.
- **DeRoo's ray trace, as data.** It is the only check here from outside the
  author's own derivation. Commit the traced positions with their provenance,
  if DeRoo agrees, so a test holds the code to them on every push.
- **The local grating equation**, at every groove, with the local tilt and
  arc-length spacing. A different identity from the groove function itself.
- **Validity measured, not estimated.** With the exact path function
  computable, the series residual F_exact − F_series across the aperture is a
  number, not an estimate. That answers the question the slides leave open
  (slide 11). With Rayleigh's λ/4 on the residual path, it gives an
  aberration-limited resolving power beside `gratinglab.resolution`.

## 7. Outputs

- **Groove positions to a `.txl`.** `io/txl.py`'s `groove_edges` takes a period
  callable and builds positions by stepping from each groove's leading edge.
  Over tens of thousands of grooves, any error in where the period is sampled
  accumulates. A groove map already has exact positions, so the writer should
  take positions directly. That is a small addition to `io/txl.py`. The size
  of the stepping error should be measured before deciding whether it matters
  for existing uses.
- **An efficiency map.** At each point, the local period, groove direction
  and surface tilt give a local α, β and γ, and `to_problem()` (the roadmap's
  design-description contract) gives a `Problem` for the existing solvers.
  This is the local-grating approximation. It needs the surface's radius of
  curvature to be much larger than the period, and it inherits every limit of
  the solver behind it, including the open question about coated scalar
  results in the validation inventory (#86).

## 8. The first slice

**In.**

- A plane and a cylinder.
- The plane-wave hologram constructor.
- Exact groove positions by root finding, with the closed form above as a
  test.
- The linear-position and arc-length conversion.
- Positions written to `.txl`.
- Tests:
  - the slides' quadratic and quartic in their regimes,
  - the local Littrow identity,
  - DeRoo's trace, if the data can be committed.

**Out**, for later slices:

- point sources at finite distance, and focusing (B different from A),
- the sphere and the general a_ij surface,
- out-of-plane sources and off-plane radial grooves,
- the F_ij analysis API,
- the efficiency map.

## 9. Open, for a person to decide

1. Whether this is the design description's first content, as proposed, or a
   module of its own that the design description later absorbs.
2. Whether the F_ij series is public API (aberration analysis for designers)
   or test-only.
3. Which position the design stores as primary: linear w or arc length.
4. Whether DeRoo's ray-trace data can be committed, and under what credit.
5. Whether off-plane radial grooves belong in the same groove-function object,
   or need their own because they are defined on a fan, not a plane.
