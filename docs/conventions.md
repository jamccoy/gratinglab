# Conventions

This document is normative: every solver backend, importer, and exporter in
`gratinglab` obeys it. Where an external code or a source paper disagrees with it, the
difference is absorbed in an adapter at the boundary and never inside a solver.

Projects that compare methods most often fail on conventions. Two codes that disagree
by a sign in the grating equation, or by an obliquity factor, produce plots that appear
to show a physical result but in fact reflect a bookkeeping error. These choices are
therefore fixed once, in this document, and are not revisited in code review.

---

## 1. Units

At the public API boundary, every length is in nanometres and every angle is in degrees.

| Quantity | Unit |
|---|---|
| Wavelength λ | nm |
| Period p | nm |
| Groove depth, coating thickness, RMS roughness σ | nm |
| All angles (α, β, γ, δ, ζ) | degrees |
| Wavenumber k = 2π/λ | nm⁻¹ |

Nanometres give human-readable magnitudes in both target regimes: soft X-ray work uses
λ ≈ 0.6–5 nm with p ≈ 150–350 nm, and visible work uses λ ≈ 400–700 nm with
p ≈ 1400 nm. In micrometres, soft-X-ray wavelengths would sit at 6 × 10⁻⁴, which is
error-prone both to read and to type.

Angles are given in degrees at the API and held in radians internally, and the
conversion happens once, at construction. No function takes an angle whose unit is
ambiguous from its name.

### Inside the metrology package

`gratinglab.metrology` reads instrument files, and instruments do not report in
nanometres. Four unit regimes therefore coexist in this package, each confined to its
own layer:

| Where | Unit | Set by |
|---|---|---|
| `core.processing.load_afm_data` returns heights | **metres** | the file formats (Nanoscope and Gwyddion both write SI) |
| lateral axis, scan width | **µm** | the same file formats, and the unit an AFM operator works in |
| `raw_data` onward, profile heights | **nm** | the project convention above |
| `BoundaryProfile.x_norm` / `y_norm`, and `.ggp` | **dimensionless**, fraction of the period | PCGrate's polygonal border |

Conversion occurs in exactly two places: `raw_data` (metres → nm) and
`normalize_profile` (nm → fraction of period). No function accepts a length whose
regime is ambiguous from its name, which is why the fields are named `y_avg_nm`,
`x_avg_um` and `y_norm` instead of all three being called `y`.

Because the normalised layer expresses heights as fractions of the period, the period
must be carried separately. `y_norm` is a fraction of the period, so physical depth is
`y_norm × period`. An incorrect period therefore produces an incorrect groove depth,
and groove depth enters the scalar phase term directly instead of serving only as a
label on the result. `BoundaryProfile.to_problem` carries the measured period for this
reason; a `.ggp` file cannot carry it.

---

## 2. Time convention and the sign of the extinction coefficient

The time dependence is `exp(-iωt)`, and a plane wave is written

$$\mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0 e^{i\mathbf{k}\cdot\mathbf{r}} e^{-i\omega t}$$

This follows McCoy, *Scalar Treatment of Diffraction Gratings* (2022), slide 7, and the
thesis throughout.

The consequence of this choice, and the only part of this section that most readers
need, is that a lossy material has

$$n = n' + i k, \qquad k > 0$$

The opposite time convention `exp(+iωt)` would require `n = n' - ik`. Optical constant
tables must be checked against this convention on import. CXRO/Henke tables give δ and
β with `n = 1 - δ + iβ`, which is already correct for `exp(-iωt)` with β > 0.
refractiveindex.info gives `n` and `k` with `k ≥ 0`, which is also already in this
convention; the few negative `k` values in its database (round-off, or Kramers-Kronig
artifacts) are clipped to zero on import and recorded in the table's `source`.

---

## 3. Coordinate frame

The coordinate frame is right-handed and follows Heilmann, Huenemoerder, McCoy &
McEntaffer (2024), §2.1:

| Vector | Direction | Meaning |
|---|---|---|
| `d̂` | `x̂` | dispersion direction; periodicity is along this axis |
| `n̂` | `ŷ` | grating normal, pointing into the vacuum half-space |
| `ĝ` | `-ẑ` | groove direction; grooves run along this axis |

with `d̂ × ĝ = n̂`.

The grating surface occupies `y ≤ g(x)`, where `g(x + p) = g(x)` is the groove profile.
Vacuum (or the incident medium) is `y > g(x)`.

### Wave vectors

$$\mathbf{k}_i = k\left[-\sin\alpha\,\sin\gamma\,\hat{x} - \cos\alpha\,\sin\gamma\,\hat{y} + \cos\gamma\,\hat{z}\right]$$

$$\mathbf{k}_m = k\left[\sin\beta_m\,\sin\gamma\,\hat{x} + \cos\beta_m\,\sin\gamma\,\hat{y} + \cos\gamma\,\hat{z}\right]$$

The incident wave travels toward the grating (negative ŷ component), and reflected
orders travel away from it (positive ŷ component). `k_z` is conserved for all orders.

Every wave vector carries `k_z = +k cos γ > 0` for γ < 90°, so all propagation is
toward `−ĝ`. The choice `ĝ = −ẑ` in the frame table is forced by `d̂ × ĝ = n̂` and is
not a free choice.

### Angles

- **α**: azimuthal angle of incidence, measured from the grating normal `n̂` in the
  plane containing `n̂` and `d̂`. Its range is α ∈ (−90°, 90°), enforced by
  `Illumination`, because an incidence angle at or past grazing carries no power onto
  the grating.
- **β_m**: azimuthal angle of diffracted order *m*, measured from the same reference.
- **γ**: half-angle of the diffraction cone about the groove axis. The cone's axis is
  the groove-axis line, and γ is the acute angle between a wave vector and that line.
  γ is the angle that sweeps out the cone, and α gives the azimuthal position around it.

  Stated as a directed angle, γ is the polar angle in the spherical system whose pole
  is `+ẑ = −ĝ`, since every wave vector above carries `k_z = +k cos γ`. At γ = 1.5°,
  for example, every ray lies 1.5000° from `+ẑ` and 178.5000° from `ĝ`, so the cone
  opens along `−ĝ` and not along `ĝ`. A drawing that opens it along `ĝ` places every
  ray 180° out. Code that consumes γ only through the scalars `sin γ` / `cos γ` cannot
  detect the difference, which is why this convention went unstated until something had
  to draw the cone (see `docs/findings.md`).

  In operational terms, `k_z = k cos γ` is conserved, and `sin γ = 1` (γ = 90°) is the
  in-plane case. Extreme off-plane mounts have small γ, and the soft-X-ray work uses
  γ ≈ 1.5°. The range is γ ∈ (0°, 90°], enforced by `Illumination`: γ = 90° is the
  in-plane case, and γ = 0° would send the wave along the grooves, which produces no
  diffraction.
- **δ**: groove-facet (blaze) angle of a sawtooth profile.
- **ζ**: graze angle onto the groove facet, `sin ζ = sin γ · cos(δ - α)`.

The geometry is stored as direction cosines instead of angles (see `Illumination`). The
three constructors (`classical`, `conical`, `offplane`) allow each community to pass
the angles it uses, and all three resolve to the same internal representation.

### Direction of the profile parameter

$$\hat{t} = -\hat{d}$$

`Profile.height(t)` describes a shape in a normalised parameter `t`. `profiles.py` does
not state which physical direction `+t` points in, because nothing required it until a
profile and a ray were drawn in the same frame. The wave conventions above determine
it: `+t` runs along `-d̂`, that is, against the dispersion direction.

Two independent derivations give this result.

**From the blaze direction.** `geometry.blaze_direction` returns `β_b = 2δ - α`.
Reflecting the projected incident direction `(-sin α, -cos α)` about a facet with
outward normal `n̂_f` gives `sin β = sin(2δ - α)` for `n̂_f = (+sin δ, cos δ)`, and
`sin β = -sin(2δ + α)` for `n̂_f = (-sin δ, cos δ)`. At the reference geometry
(δ = 29.5°, α = 25°) these give +34.0° and −84.0° respectively. The active facet's
outward normal therefore sits at azimuth `+δ`, which means that the facet descends with
increasing `x`. However, `Blazed.slope()` returns `+tan δ` for `t ≤ apex`, so the
active facet rises with `t`. A facet that rises in `t` and descends in `x` implies
`t̂ = -d̂`.

**From the Kirchhoff integral.** `solvers/scalar.py` forms
`G_m = ∫ e^{+iΦ_m(t)} e^{-2πimt} dt`, whereas the time convention of §2 and the `k`
vectors above give `∫ e^{-iΦ_m(t)} e^{-2πimt} dx/p` when `t ↔ +x`. Substituting
`t ↔ -x` maps one integral onto the complex conjugate of the other.

No numerical result depends on this choice. Conjugation leaves `|G_m|²`, the only
quantity any caller consumes, bit-identical, which is why the convention remained
consistent but unstated for so long. It matters only where a drawing places a facet
relative to a ray. If the direction is reversed, the incident beam strikes the
anti-blaze facet and the blaze arrow leaves through the back of the active facet.

The relation is recorded here instead of being enforced by flipping the signs in
`scalar.py`. That change would also preserve the modulus, but it would move the
module's headline formula away from the form transcribed from the two primary
references. The references themselves state no handedness and quote only `|G_m|²`; this
omission is itself the finding. The behaviour is pinned by
`tests/test_geometry.py::TestFacetHandedness`; see also `docs/findings.md`.

Not every consumer takes only `|G_m|²`. The integral solver places the curve
itself, and through 0.1.12 it placed it at `x = t P`, solving every asymmetric
groove as its mirror image. It now places it at `x = (1 − t) P`
(`solvers/integral/_boundary.py`), and
`tests/test_cross_method.py::TestAgreementOnAnAsymmetricGroove` holds both
solvers to the same blaze order. Files that run along `+d̂`, PCGrate's `.ggp`
and the metrology pipeline's `x_norm`, cross into `t` in one place,
`FromProfileData.along_dispersion` and its inverse `dispersion_axis`, as the
preamble of this document requires of an external convention.

#### Alternate names in the literature

Two naming conflicts in the literature are described here so that neither is resolved
silently.

In spherical coordinates with the groove axis as the pole, γ is strictly the polar
angle and α is azimuthal. The off-plane-grating literature, however, often calls α "the
polar" or "the incidence" angle. In the classical special case (in-plane, γ = 90°), α is
the ordinary incidence angle measured from the grating normal, and that name is carried
into the conical case by habit. Both usages exist; this document uses the strict
spherical reading.

Some sources parametrize the cone by its complement. Where this document uses γ
directly (small γ = grazing, γ = 90° = in-plane), others use `φ ≡ 90° − γ` (small φ =
in-plane, φ = 0° = grazing), and correspondingly may write `θ` where this document
writes `α`. The translation is as follows:

| Here | Elsewhere | Relation |
|---|---|---|
| α | θ | same angle |
| γ | φ | `φ ≡ 90° − γ`, so `sin γ = cos φ` |

The `θ` in this table is unrelated to `θ_c`, the total-external-reflection critical
angle used elsewhere in this project (`docs/theory/scalar.md` §7). The two never appear
in the same equation, but a reader of both documents should not conflate them.

### Conformal layers

`Problem.layers` lists the layers on a coated grating, innermost first. Layer
`k` occupies every point within its `thickness` of the surface beneath it and
outside that surface, so its outer surface is the true offset of the one below
along the normal: a convex corner becomes an arc of radius `thickness` and a
concave corner stays sharp (decided 2026-10-27, as the closest description of
an evaporated or sputtered film). Thickness is in nm, like every length here.
`coating` (or `substrate`) remains the material the grating is cut into. A
solver that cannot model layers refuses a problem that has any
(`Capabilities.handles_layers`).

---

## 4. The grating equation

$$\sin\alpha + \sin\beta_m = \frac{m\lambda}{p\,\sin\gamma}, \qquad m = 0, \pm1, \pm2, \dots$$

This is the *generalized (conical) grating equation*, ISSI eq. (7), which reduces to the
familiar in-plane form at `sin γ = 1`. It is the only form used internally.

The sign convention is the `+ sin β_m` form: reflected orders with `m > 0` diffract to
the side of specular reflection given by increasing `sin β_m`.

Transmission gratings are sometimes written `mλ/p = sin α − sin β_m` by redefining
`β_m → β_m − π`. This redefinition is not used internally. A transmitted order keeps
the same convention and is identified by `cos β_m < 0`, and the R/T flag is carried
explicitly on the order record.

### Propagating orders

Order *m* propagates if and only if

$$\left|\sin\alpha + \frac{m\lambda}{p\,\sin\gamma}\right| \leq 1$$

Orders that fail this condition are evanescent. They are retained in the record with
`efficiency = 0.0` and `propagating = False`; they are never dropped silently and never
represented as `NaN`. Dropping them would shift order indices between wavelengths,
which is the single most common source of off-by-one errors in scan post-processing.
(The existing prototype scripts index efficiency columns as `data[i+11]`, and this rule
prevents that kind of fragility.)

---

## 5. Efficiency normalization

Efficiencies are absolute, referenced to the power incident on the grating:

$$\mathscr{E}_m = \frac{I(m)}{I_0}$$

This matches PCGrate's `Eff.TE(m,R)` columns, so the on-disk reference corpus imports
without rescaling.

Relative (referenced-to-mirror) efficiency is available as a derived quantity but is
never the default. At γ = 1.5° the difference between absolute and relative efficiency
is a large factor, and mislabelling one as the other silently invalidates a throughput
budget.

For a lossless structure, `Σ_m over propagating E_m = 1`, and this identity is asserted
in the tests.

### Energy conservation as a check

Energy balance is treated asymmetrically, according to two rules:

- **A deficit is ordinary.** Power may go into absorption or into evanescent orders, or
  an approximate method may fail to capture it.
- **An excess is not ordinary.** No passive grating returns more power than it
  receives, whatever method is used to compute it. By default
  `checks.check_energy_balance` enforces only this weaker condition, because it has no
  legitimate exception. This check identified the unphysical finite-conductivity run in
  the reference corpus, in which the sum reached 3.6.

Scalar theory as normally written does not conserve energy, and this is accepted as a
deliberate trade-off. With the order-dependent phase of ISSI eq. (15) and thesis
Appendix-D.tex:651, the summed efficiency deviates from unity in both directions, with
measured values up to 1.61 across a broad parameter sweep.

An alternative exists and was rejected. Dropping `cos β_m` from the phase restores a
genuine Fourier pair and conserves energy to 1e-10, but it violates Lorentz reciprocity
(measured 0.44) and reproduces the blaze direction only at Littrow. Reciprocity requires
the phase to be symmetric under `α ↔ β_m`, whereas Parseval requires a single fixed
function that excludes `β_m`, so the two cannot hold together. Reciprocity is the
invariant that is kept.

The deviation does not indicate a fault in the implementation: in the shallow-groove
limit the sum is exactly 1.00000, and it degrades smoothly with depth (0.0001 → 1.00000;
0.10 → 0.92902). The deviation scales with the phase excursion across the groove, that
is, with depth relative to wavelength and hence with working order. It does not scale
with λ/p or with the number of propagating orders.

The full derivation and evidence are given in [`theory/scalar.md`](theory/scalar.md) §5.

The deviation is reported and never rescaled; it is recorded as a `Provenance` warning.
The extent to which an approximate theory departs from a conservation law is what a
scalar validity map should show, and normalising the deviation away would destroy that
information.

### The symmetric obliquity factor

Thesis Appendix D includes two factors that ISSI §2.1 does not. One of them has a
correct counterpart and the other does not:

| Form | Source | Status |
|---|---|---|
| `E_m = O_m ‖G_m‖²`, `O_m = 4 cos α cos β_m / (cos α + cos β_m)²` | first-order Rayleigh perturbation theory | **correct; what the code does** |
| `E_m = ‖G_m‖²` | ISSI §2.1, eq. groove_func | wrong by `1/O_m`; up to 64% off |
| `E_m ∝ [cos β_m / cos α] ‖G_m‖²` | thesis Appendix-D.tex:418 | removed: asymmetric, breaks reciprocity |
| renormalise so `Σ_m E_m = 1` | thesis Appendix-D.tex:420 | removed as an error |

The two factors differ in their symmetry. Appendix D's `cos β_m / cos α` is asymmetric
under `α ↔ β_m` and does violate Lorentz reciprocity, which is why it was rejected. The
conclusion that was then drawn, that the efficiency is the norm-squared Fourier
coefficient alone, went one step too far. `O_m` is symmetric, so reciprocity is
unaffected by it, and it is the factor that brings the shallow-groove limit into
agreement with a theory derived from the boundary condition instead of from a
transmittance function. See [`theory/scalar.md`](theory/scalar.md) §3 and
[`tests/test_perturbation.py`](../tests/test_perturbation.py).

By the AM–GM inequality `O_m ≤ 1`, and it equals 1 exactly at Littrow, at specular, and
for every `m = 0`. It can therefore only reduce an efficiency, and it leaves the
shallow-limit energy identity unchanged.

The renormalization has no comparable justification, since forcing the sum to unity
would erase the signal described above. It is not offered as an option, because an
error offered as a choice invites misuse.

---

## 6. Wave-vector signs

Wave-vector components are always stored with their signs:

- `k_y,i < 0`: the incident wave travels toward the grating.
- `k_y,m > 0`: the reflected order travels away from it.
- `k_y,m` is purely imaginary for an evanescent order.

Efficiency is formed from magnitudes:

$$\mathscr{E}_m = \frac{|k_{y,m}|}{|k_{y,i}|}\left\|\frac{A_m}{A_0}\right\|^2 = \frac{\cos\beta_m}{\cos\alpha}\left\|\frac{A_m}{A_0}\right\|^2$$

The convention is stated because the source material is ambiguous about it. In the
thesis, the Green's-theorem derivation at Chapter-2.tex:1006–1007 yields the stated
result only if `k_y < 0` for the incident wave, while the equation at line 905 writes `k''_{y,n}/k_y ≡
cos β / cos α`, which requires `k_y > 0`. Each can be recovered consistently with care,
but an ambiguity of this kind can flip a sign in code without causing any visible
failure.

The rule is therefore that the fields use signed components and the efficiency ratio
uses magnitudes. Any solver that needs the unsigned value writes `abs()` explicitly.

---

## 7. Polarization

- **TE**: electric field parallel to the grooves (`ĝ`).
- **TM**: magnetic field parallel to the grooves.
- **Unpolarized**: `(TE + TM) / 2`, computed and never assumed.

In conical mounts, TE and TM are defined with respect to the grooves and not with
respect to the plane of incidence; the two definitions coincide only in-plane. This is
stated explicitly because the other convention is common in the RCWA literature and is
a frequent source of cross-code disagreement in conical geometry.

In the `(E_z, B_z)` basis of the conical field decomposition (grooves along z), TE is a
pure incident `E_z` state (`(p_z, q_z) = (1, 0)`; Goray & Schmidt's "`B_z = 0`") and TM
is a pure incident `B_z` state (`(0, 1)`; their "`E_z = 0`"), and each reduces to the
in-plane definition above. With finite conductivity, a pure incident state diffracts
into orders that carry both components, since cross-polarization conversion does occur
off-plane. The per-order efficiency is then the flux-weighted sum
`(γ_m/γ₀)(|E_m|² + |B_m|²)`, which gives one number per order. The conversion is
identically zero only for a perfect conductor.

At grazing incidence in the soft X-ray, reflectivity is nearly polarization-independent
(thesis §C.1.2), which is why the existing PCGrate corpus carries only `Eff.TE`. That
approximation is justified in that regime but not in general, and it is not built into
the core.

---

## 8. Symbols

The two primary references use `d` for two different quantities:

| Symbol | McCoy thesis | ISSI chapter | **gratinglab** |
|---|---|---|---|
| period | `d` | `p` | **`period`** |
| groove depth | `h₀` | `d` | **`depth`** |
| blaze / facet angle | `δ` | `δ` | `blaze_angle` |
| facet graze angle | `ζ` | `ζ` | `facet_graze` |
| order index | `n` | `m` | `order` |
| index decrement | `δ` | `δ` | **`decrement`** |
| index absorption | `β` | `β` | **`absorption`** |

The last two rows present the same problem in a second place. The optics literature
writes `n = 1 - δ + iβ`, but both letters already have meanings in this project: `δ` is
the blaze angle and `β` is the diffracted-order angle throughout. `fresnel.reflectivity`
is called a few lines from `facet_graze(gamma, blaze_angle, alpha)`, and two meanings of
`δ` within one call stack would allow a sign error to be written without being noticed.
The materials layer therefore uses unambiguous names.

Bare `d` never appears in the code, whether as a variable, as a keyword argument, or in
a docstring formula without an accompanying definition. When an equation is transcribed
from either source, the first step is to rewrite `d` explicitly as `period` or `depth`.

---

## 9. Reference formulas, in gratinglab symbols

For an ideal sawtooth (blazed) reflection grating:

| Quantity | Formula |
|---|---|
| Facet depth | `depth = period · tan(blaze_angle)` |
| Facet graze angle | `sin ζ = sin γ · cos(blaze_angle − α)` |
| Sawtooth phase shift | `φ = k · depth · sin γ · [cos α + cos β_m]` |
| Flux obliquity | `O_m = 4 cos α cos β_m / (cos α + cos β_m)²` |
| Scalar DE | `E_m = O_m · sinc²(φ/2 − mπ)`, `sinc(x) ≡ sin(x)/x` |
| Blaze direction | `β_b = 2·blaze_angle − α` |
| Blaze wavelength | `m·λ_b = 2·period·sin ζ·sin(blaze_angle)` |
| Finite-N interference | `[sin(Ns)/(N sin s)]²`, `s ≡ [sin α + sin β_m]·sin γ·period·π/λ` |
| Resolving power (Rayleigh) | `R = |m|·N`; `Δλ = λ/R`. `N` is the illuminated groove count |

The phase shift includes the factor `sin γ`. The thesis expression at
Appendix-D.tex:651 omits it and is inconsistent with the thesis's own Φ_b at line 672,
which includes it. The thesis form is the in-plane special case and is valid only under
the substitution `λ → λ csc γ` introduced at line 446. This project uses the ISSI form,
which is general.

The finite-N factor is retained deliberately. Appendix D drops it in the `N → ∞` limit,
but the factor links efficiency to resolving power and therefore belongs in the API,
where it is now provided as `gratinglab.resolution`, with the factor itself in
`gratinglab.geometry`. The Rayleigh row above is derived from it. The first zeros sit at
`s = mπ ± π/N`, and placing the peak of one wavelength on the first zero of the other
gives `m·Δλ = λ/N` in any mount, because `sin γ` scales the line width and the
dispersion equally.

### Scalar validity guards

The following guards are reported alongside any scalar result:

- **Reduced wavelength ratio:** `λ/(period · sin γ) ≲ 0.1`. The smallness that scalar
  theory requires is judged at the reduced wavelength `λ/sin γ`, which is the wavelength
  of the transverse problem posed by the grooves (the same decoupling on which the
  integral method relies). It is not judged at `λ/period`, because that ratio passes
  extreme off-plane mounts that lie far outside the regime. In-plane, the two coincide.
- **Total external reflection** (ISSI §4): `θ_c ≈ √(2δ_n)`; the model is meaningful
  only for `ζ < θ_c`.
- **Fraunhofer smoothness** (ISSI §4): `λ > 32 · sin(ζ) · σ`.

A scalar result outside these bounds is returned with a warning flag on its provenance
record instead of being suppressed. Mapping where scalar theory breaks down is one of
the project's deliverables, so out-of-bounds values must remain available for
inspection.

---

## 10. Deviations from source references

The deviations below are recorded so that a deliberate choice is not reverted to match
a source that is in error.

| # | Source | Issue | Our choice |
|---|---|---|---|
| 1 | thesis Appendix-D.tex:651 | `Φ_g` missing `sin γ` | ISSI eq. (15) form |
| 2 | thesis Appendix-D.tex:418, 420 | obliquity factor `cos β/cos α`, plus renormalization so `Σ E = 1` | renormalization removed, and not available as an option. The obliquity factor is removed in that form (it is asymmetric, so it breaks reciprocity) and replaced by the symmetric `O_m` of §5, which perturbation theory requires |
| 2a | ISSI §2.1, groove function | `E_m = ‖G_m‖²` with no flux factor at all | wrong by `1/O_m`, by up to 64% away from Littrow. This is the one place where this project's normative reference is the source in error; see `theory/scalar.md` §3 for the evidence, which is independent of both sources |
| 3 | thesis Appendix-D.tex:699 | `∂Φ_b/∂α` missing factor `period` (dimensionally required) | corrected |
| 4 | thesis Appendix-D.tex:462, 622 | `for n = 0, ±1, …` on expressions singular at `n = 0` | `G₀` stated separately: `A₀·𝒲/period` (square), `A₀/2` (sawtooth) |
| 5 | thesis Chapter-2.tex:1016, 1045 | claims PCGrate internally multiplies by Fresnel reflectivity in perfect-conductivity mode, which contradicts the `ΣE = 1` result proved just above | Treat `ΣE ≈ R_F(ζ)` as an emergent result of an impedance-type (Leontovich) boundary condition. Follow Goray & Schmidt (2010) for the rigorous finite-conductivity case |
| 6 | thesis Chapter-2.tex:1006–1007 | `A_n'` should be `A_n''`; evanescent exponential missing its `y` | corrected |
| 7 | thesis Chapter-2.tex:905 vs 1006 | inconsistent `k_y` sign convention | §6 above |
| 8 | PyXFocus `grating.py` (`62454cd`), `blazeYaw`, `blazeAngle`, `blazeAngle2`, `litBetaAlpha` | a different frame, not an error: `inc` is the graze onto the grating surface, `sin inc = sin γ cos α`, and the yaw ψ is the beam's projection on the surface measured from the grooves, `tan ψ = sin γ sin α / cos γ`. Its dispersion cosine is signed opposite to ours, so `blazeAngle` returns −α where `blazeAngle2` returns α. `litBetaAlpha`'s "cone angle" `beta1` is `cos γ`, not γ | `geometry.littrow_alpha` in (α, γ); the mapping is pinned from direction cosines in `tests/test_littrow.py`, never by equating `cos inc` with `sin γ`, which is wrong by up to ~1 at grazing incidence |

Items 3, 4, 6 are transcription errors with no effect on published results. Items 1, 2,
5, 7 change what the code computes. Item 8 is no error at all, only a second frame, and
is recorded because the obvious reading of it is wrong.

---

## 11. Adapters

External formats are converted at the boundary, and an adapter is the only place in
which a foreign convention may appear.

| Format | Direction | Notes |
|---|---|---|
| PCGrate-SX `.txt` efficiency tables | import | 7 header lines; tab-separated; quoted `"Eff.TE(m,R)"` headers carry the order index and R/T flag; `--` marks missing/evanescent. Order mapping comes from the headers, never from a column offset |
| CXRO / Henke optical constants | import | `(energy_eV, δ, β)` → `n = 1 − δ + iβ` on a wavelength grid |
| refractiveindex.info YAML | import | Wavelength in µm. `n` and `k` tabulated, or `n` from one of nine dispersion formulas sampled adaptively to 1e-7; stored as `decrement = 1 − n`, `absorption = k`, so a glass has a negative decrement. `tabulated n2` is the nonlinear Kerr index and is refused |
| PCGrate `.ari` refractive index | export | `wavelength_nm  δ  β`, space-delimited |
| PCGrate `.ggp` polygonal border | import + export | Normalised: `t` on [0,1], `y` as a fraction of the period. Carries no period, no provenance and no version. Five control-string variants exist on disk; all are read. The points are written alone by default, and `header=True` writes the pair the PCGrate-S(X) 6.1 User's Guide specifies in 10.2.2 (p. 63): `3 0 - 20 harmonics - Polygonal type` then `Period: 1 PSC: 1`. Between 2 and 5000 data strings, x on [0,1] non-decreasing, y on [0,100], first point `0 0`, last `1 0` (10.2.3, p. 63). A measured profile should instead reach a solver through `BoundaryProfile.to_problem`, which keeps the period |
| Nanoscope `.spm` (AFM) | import | Binary. Heights in metres. The `Sens. Zsens` / `Sens. ZsensSens` distinction is load-bearing: the wrong one gives heights 5.13× too small |
| Gwyddion text export (AFM) | import | 4 header lines, tab-separated, heights in metres, scan width scraped from the header comments |
| GSolver, RETICOLO | planned | |

Round-tripping a `Problem` through an adapter and back must leave it unchanged. This
requirement is enforced by a test, which catches sign and normalization errors before
any physics is written.
