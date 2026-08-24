# Conventions

**This document is normative.** Every solver backend, importer, and exporter in
`gratinglab` obeys it. Where an external code or a source paper disagrees, the
difference is absorbed in an adapter at the boundary — never inside a solver.

Cross-method comparison projects usually fail on conventions, not physics. Two codes
that disagree by a sign in the grating equation, or by an obliquity factor, produce
plots that look like a physics result and are actually a bookkeeping error. So these
choices are fixed once, here, and are not revisited in code review.

---

## 1. Units

**Nanometres for every length, degrees for every angle, at the public API boundary.**

| Quantity | Unit |
|---|---|
| Wavelength λ | nm |
| Period p | nm |
| Groove depth, coating thickness, RMS roughness σ | nm |
| All angles (α, β, γ, δ, ζ) | degrees |
| Wavenumber k = 2π/λ | nm⁻¹ |

Nanometres serve both target regimes with human-readable magnitudes: soft X-ray runs
λ ≈ 0.6–5 nm against p ≈ 150–350 nm, visible runs λ ≈ 400–700 nm against p ≈ 1400 nm.
Micrometres would put soft-X-ray wavelengths at 6 × 10⁻⁴, which is error-prone to read
and to type.

Angles are degrees at the API and **radians internally**. Conversion happens once, at
construction. No function takes an angle whose unit is ambiguous from its name.

### Inside the metrology package

`gratinglab.metrology` reads instruments, and an instrument does not use nanometres.
Four unit regimes therefore coexist, each confined to its own layer:

| Where | Unit | Set by |
|---|---|---|
| `core.processing.load_afm_data` returns heights | **metres** | the file formats — Nanoscope and Gwyddion both write SI |
| lateral axis, scan width | **µm** | the same, and how an AFM operator thinks |
| `raw_data` onward, profile heights | **nm** | the project convention above |
| `BoundaryProfile.x_norm` / `y_norm`, and `.ggp` | **dimensionless**, fraction of the period | PCGrate's polygonal border |

Exactly two places convert: `raw_data` (metres → nm) and `normalize_profile`
(nm → fraction of period). No function accepts a length whose regime is ambiguous
from its name — hence `y_avg_nm`, `x_avg_um` and `y_norm` rather than three fields
called `y`.

**The normalised layer is why the period must travel separately.** `y_norm` is a
fraction of the period, so physical depth is `y_norm × period`: a wrong period is a
wrong groove depth, and groove depth sits inside the scalar phase term rather than
being a label on the result. `BoundaryProfile.to_problem` carries the measured
period for exactly this reason. A `.ggp` cannot.

---

## 2. Time convention and the sign of the extinction coefficient

**Time dependence is `exp(-iωt)`.**

A plane wave is written

$$\mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0 e^{i\mathbf{k}\cdot\mathbf{r}} e^{-i\omega t}$$

This follows McCoy, *Scalar Treatment of Diffraction Gratings* (2022), slide 7, and the
thesis throughout.

**Consequence, and the only thing most people need from this section:** a lossy material
has

$$n = n' + i k, \qquad k > 0$$

The opposite time convention `exp(+iωt)` would require `n = n' - ik`. Optical constant
tables must be checked against this on import. CXRO/Henke tables give δ and β with
`n = 1 - δ + iβ`, which is already correct for `exp(-iωt)` with β > 0.

---

## 3. Coordinate frame

Right-handed, following Heilmann, Huenemoerder, McCoy & McEntaffer (2024), §2.1:

| Vector | Direction | Meaning |
|---|---|---|
| `d̂` | `x̂` | dispersion direction — periodicity is along this axis |
| `n̂` | `ŷ` | grating normal, pointing into the vacuum half-space |
| `ĝ` | `-ẑ` | groove direction — grooves run along this axis |

with `d̂ × ĝ = n̂`.

The grating surface occupies `y ≤ g(x)`, where `g(x + p) = g(x)` is the groove profile.
Vacuum (or the incident medium) is `y > g(x)`.

### Wave vectors

$$\mathbf{k}_i = k\left[-\sin\alpha\,\sin\gamma\,\hat{x} - \cos\alpha\,\sin\gamma\,\hat{y} + \cos\gamma\,\hat{z}\right]$$

$$\mathbf{k}_m = k\left[\sin\beta_m\,\sin\gamma\,\hat{x} + \cos\beta_m\,\sin\gamma\,\hat{y} + \cos\gamma\,\hat{z}\right]$$

The incident wave travels toward the grating (**negative** ŷ component); reflected orders
travel away from it (**positive** ŷ component). `k_z` is conserved for all orders.

Every wave vector carries `k_z = +k cos γ > 0` for γ < 90°, so all propagation is toward
`−ĝ`. The frame table's `ĝ = −ẑ` is forced by `d̂ × ĝ = n̂` and is not a free choice.

### Angles

- **α** — azimuthal angle of incidence, measured from the grating normal `n̂` in the
  plane containing `n̂` and `d̂`. Range **α ∈ (−90°, 90°)**, enforced by
  `Illumination`: an incidence angle at or past grazing carries no power onto the
  grating.
- **β_m** — azimuthal angle of diffracted order *m*, same reference.
- **γ** — half-angle of the diffraction cone about the groove axis. The cone's axis is
  the groove-axis **line**, and γ is the acute angle between a wave vector and that
  line; γ is what sweeps out the cone, α is the azimuthal position around it.

  As a *directed* statement, γ is the polar angle in the spherical system whose pole is
  **`+ẑ = −ĝ`**, because every wave vector above carries `k_z = +k cos γ`. Concretely,
  at γ = 1.5° every ray lies **1.5000° from `+ẑ` and 178.5000° from `ĝ`: the cone opens
  along `−ĝ`, not along `ĝ`.** A drawing that opens it along `ĝ` points every ray 180°
  out. Nothing that consumes γ as the scalars `sin γ` / `cos γ` can tell the difference,
  which is why this went unstated until something had to draw it — see
  `docs/findings.md`.

  Operationally: `k_z = k cos γ` is conserved, and **`sin γ = 1` (γ = 90°) is the
  in-plane case**. Extreme off-plane mounts have small γ; the soft-X-ray work runs
  γ ≈ 1.5°. Range **γ ∈ (0°, 90°]**, enforced by `Illumination`: γ = 90° is in-plane;
  γ = 0° would send the wave along the grooves, which diffracts nothing.
- **δ** — groove-facet (blaze) angle of a sawtooth profile.
- **ζ** — graze angle onto the groove facet, `sin ζ = sin γ · cos(δ - α)`.

Storage is by **direction cosines**, not angles — see `Illumination`. The three
constructors (`classical`, `conical`, `offplane`) exist so that each community can pass
the angles it actually uses, and they all resolve to the same internal representation.

### Which way the profile parameter runs

$$\hat{t} = -\hat{d}$$

`Profile.height(t)` is a shape in a normalised parameter `t`; `profiles.py` never says
which physical direction `+t` points, because until something drew a profile and a ray in
the *same* frame, nothing needed to know. The wave conventions above fix it, and the
answer is that **`+t` runs along `-d̂`**, i.e. against the dispersion direction.

Two independent derivations agree:

**From the blaze direction.** `geometry.blaze_direction` returns `β_b = 2δ - α`.
Reflecting the projected incident direction `(-sin α, -cos α)` about a facet with outward
normal `n̂_f` gives `sin β = sin(2δ - α)` for `n̂_f = (+sin δ, cos δ)`, and
`sin β = -sin(2δ + α)` for `n̂_f = (-sin δ, cos δ)` — at the reference geometry
(δ = 29.5°, α = 25°), **+34.0°** versus **−84.0°**. So the active facet's outward normal
sits at azimuth `+δ`, meaning the facet *descends* with increasing `x`. But
`Blazed.slope()` returns `+tan δ` for `t ≤ apex`: the active facet *rises* with `t`.
Rising in `t`, descending in `x` ⟹ `t̂ = -d̂`.

**From the Kirchhoff integral.** `solvers/scalar.py` forms
`G_m = ∫ e^{+iΦ_m(t)} e^{-2πimt} dt`, whereas §2's time convention and the `k` vectors
above give `∫ e^{-iΦ_m(t)} e^{-2πimt} dx/p` when `t ↔ +x`. Substituting `t ↔ -x` maps one
integral onto the complex conjugate of the other.

**Nothing depends on the choice numerically.** Conjugation leaves `|G_m|²` — the only
quantity any caller consumes — bit-identical, which is why this was consistent but
unstated for so long. It matters only where a drawing places a facet relative to a ray:
get it backwards and the incident beam strikes the *anti-blaze* facet while the blaze
arrow leaves through the back of the active one.

Recorded rather than "fixed" by flipping the signs in `scalar.py`, which would also be
modulus-preserving but would move that module's headline formula away from the form
transcribed from the two primary references. The references themselves state no
handedness — they quote only `|G_m|²` — which is itself the finding. Pinned by
`tests/test_geometry.py::TestFacetHandedness`; see also `docs/findings.md`.

#### Alternate names in the literature

Two clashes worth flagging rather than silently picking a side on.

**Precisely, γ is the polar angle and α is azimuthal** — spherical coordinates with the
groove axis as the pole. Informally, though, off-plane-grating literature often calls α
"the polar" or "the incidence" angle instead, because in the classical (in-plane, γ = 90°)
special case α *is* the ordinary incidence angle read off the grating normal, and people
carry that name into the conical case out of habit. Both usages exist; this document uses
the strict spherical reading.

**Some sources parametrize the cone by its complement.** Where this document uses γ
directly (small γ = grazing, γ = 90° = in-plane), others use `φ ≡ 90° − γ` (small φ =
in-plane, φ = 0° = grazing) — and correspondingly may write `θ` where this document writes
`α`. The translation:

| Here | Elsewhere | Relation |
|---|---|---|
| α | θ | same angle |
| γ | φ | `φ ≡ 90° − γ`, so `sin γ = cos φ` |

Note `θ` here is unrelated to `θ_c`, the total-external-reflection critical angle used
elsewhere in this project (`docs/theory/scalar.md` §7) — the two never appear in the same
equation, but a reader skimming both documents should not conflate them.

---

## 4. The grating equation

$$\sin\alpha + \sin\beta_m = \frac{m\lambda}{p\,\sin\gamma}, \qquad m = 0, \pm1, \pm2, \dots$$

This is the **generalized (conical) grating equation**, ISSI eq. (7), and it reduces to
the familiar in-plane form at `sin γ = 1`. It is the only form used internally.

**Sign convention:** the `+ sin β_m` form. Reflected orders with `m > 0` diffract to the
side of specular reflection given by increasing `sin β_m`.

Transmission gratings are sometimes written `mλ/p = sin α − sin β_m` by redefining
`β_m → β_m − π`. **We do not do this internally.** A transmitted order keeps the same
convention and is distinguished by `cos β_m < 0`, with the R/T flag carried explicitly on
the order record.

### Propagating orders

Order *m* propagates iff

$$\left|\sin\alpha + \frac{m\lambda}{p\,\sin\gamma}\right| \leq 1$$

Orders failing this are evanescent. They are **retained in the record with
`efficiency = 0.0` and `propagating = False`**, never silently dropped and never
represented as `NaN`. Dropping them makes order indices shift between wavelengths, which
is the single most common source of off-by-one errors in scan post-processing. (The
existing prototype scripts index efficiency columns as `data[i+11]`; that fragility is
exactly what this rule prevents.)

---

## 5. Efficiency normalization

**Efficiencies are absolute**, referenced to the power incident on the grating:

$$\mathscr{E}_m = \frac{I(m)}{I_0}$$

This matches PCGrate's `Eff.TE(m,R)` columns, so the on-disk reference corpus imports
without rescaling.

Relative (referenced-to-mirror) efficiency is available as a derived quantity, never as
the default. At γ = 1.5° the absolute-vs-relative distinction is a large factor, not a
rounding detail, and mislabelling it silently invalidates a throughput budget.

For a lossless structure, `Σ_m over propagating E_m = 1`. This is asserted in tests, not
assumed.

### Energy conservation is a check, not a guarantee

Two rules, and they point in different directions:

- **A deficit is ordinary.** Power goes into absorption, into evanescent orders, or is
  simply missed by an approximate method.
- **An excess is not.** No passive grating returns more power than it receives, by any
  method. `checks.check_energy_balance` enforces only this weaker direction by default,
  because it never has a legitimate exception. It is what identified the unphysical
  finite-conductivity run in the reference corpus, where the sum reached 3.6.

**Scalar theory as normally written does not conserve energy**, and that is a deliberate
trade rather than a defect to fix. With the order-dependent phase of ISSI eq. (15) and
thesis Appendix-D.tex:651, summed efficiency deviates from unity in both directions —
measured up to 1.61 across a broad parameter sweep.

The alternative is available and rejected: dropping `cos β_m` from the phase restores a
genuine Fourier pair and conserves energy to 1e-10, but **violates Lorentz reciprocity**
(measured 0.44) and reproduces the blaze direction only at Littrow. Reciprocity requires
the phase to be symmetric under `α ↔ β_m`; Parseval requires a single fixed function that
forbids `β_m`. They cannot both hold. **Reciprocity is the invariant kept.**

The implementation is not at fault: in the shallow-groove limit the sum is 1.00000
exactly, degrading smoothly with depth (0.0001 → 1.00000; 0.10 → 0.92902). The deviation
scales with phase excursion across the groove — depth relative to wavelength, hence
working order — not with λ/p or propagating-order count.

Full derivation and evidence: [`theory/scalar.md`](theory/scalar.md) §5.

The deviation is **reported, never rescaled** — recorded as a `Provenance` warning. How
far an approximate theory strays from a conservation law is exactly what a scalar
validity map should show, and normalising it away would destroy that information.

### The obliquity factor — the symmetric one, and only that one

Thesis Appendix D carries two extra factors that ISSI §2.1 does not. One of them has a
correct counterpart; the other does not:

| Form | Source | Status |
|---|---|---|
| `E_m = O_m ‖G_m‖²`, `O_m = 4 cos α cos β_m / (cos α + cos β_m)²` | first-order Rayleigh perturbation theory | **correct — what the code does** |
| `E_m = ‖G_m‖²` | ISSI §2.1, eq. groove_func | wrong by `1/O_m`; up to 64% off |
| `E_m ∝ [cos β_m / cos α] ‖G_m‖²` | thesis Appendix-D.tex:418 | removed — asymmetric, breaks reciprocity |
| renormalise so `Σ_m E_m = 1` | thesis Appendix-D.tex:420 | removed — an error |

**The distinction is symmetry.** Appendix D's `cos β_m / cos α` is asymmetric under
`α ↔ β_m` and genuinely does violate Lorentz reciprocity, which is why it was rejected.
The conclusion drawn from that — that the efficiency is the norm-squared Fourier
coefficient and nothing else — went one step too far. `O_m` is symmetric, so reciprocity
survives it untouched, and it is the factor that makes the shallow-groove limit agree
with a theory derived from the boundary condition rather than from a transmittance
function. See [`theory/scalar.md`](theory/scalar.md) §3 and
[`tests/test_perturbation.py`](../tests/test_perturbation.py).

`O_m ≤ 1` by AM–GM and equals 1 exactly at Littrow, at specular, and for every `m = 0`,
so it can only reduce an efficiency and leaves the shallow-limit energy identity alone.

The renormalization has no such rehabilitation: forcing the sum to unity would erase
exactly the signal described above. It is not offered as an option, because offering an
error as a choice is a footgun, not a service.

---

## 6. Wave-vector signs — the rule that prevents silent errors

Wave-vector components are stored **signed**, always:

- `k_y,i < 0` — incident wave travels toward the grating.
- `k_y,m > 0` — reflected order travels away from it.
- `k_y,m` is **purely imaginary** for an evanescent order.

Efficiency is formed from **magnitudes**:

$$\mathscr{E}_m = \frac{|k_{y,m}|}{|k_{y,i}|}\left\|\frac{A_m}{A_0}\right\|^2 = \frac{\cos\beta_m}{\cos\alpha}\left\|\frac{A_m}{A_0}\right\|^2$$

This matters because the source material is ambiguous about it. In the thesis, the
Green's-theorem derivation at Chapter-2.tex:1006–1007 only yields the stated result if
`k_y < 0` for the incident wave, while eq. at line 905 writes `k''_{y,n}/k_y ≡
cos β / cos α`, which requires `k_y > 0`. Both are internally recoverable with care, and
that is precisely the problem: it is the kind of ambiguity that flips a sign in code
without anything failing loudly.

**The rule: signed components in the fields, magnitudes in the efficiency ratio.** Any
solver that needs the unsigned value writes `abs()` explicitly.

---

## 7. Polarization

- **TE** — electric field parallel to the grooves (`ĝ`).
- **TM** — magnetic field parallel to the grooves.
- **Unpolarized** — `(TE + TM) / 2`, computed, never assumed.

In conical mounts TE/TM are defined with respect to the **grooves**, not the plane of
incidence. The two coincide only in-plane. This is stated because the opposite
convention is common in the RCWA literature and is a frequent source of cross-code
disagreement in conical geometry.

At grazing incidence in the soft X-ray, reflectivity is nearly polarization-independent
(thesis §C.1.2), which is why the existing PCGrate corpus carries only `Eff.TE`. That is
a justified approximation *for that regime*, not a general one, and it is not baked into
the core.

---

## 8. Symbols — resolving a live collision

The two primary references use `d` to mean **opposite things**:

| Symbol | McCoy thesis | ISSI chapter | **gratinglab** |
|---|---|---|---|
| period | `d` | `p` | **`period`** |
| groove depth | `h₀` | `d` | **`depth`** |
| blaze / facet angle | `δ` | `δ` | `blaze_angle` |
| facet graze angle | `ζ` | `ζ` | `facet_graze` |
| order index | `n` | `m` | `order` |
| index decrement | `δ` | `δ` | **`decrement`** |
| index absorption | `β` | `β` | **`absorption`** |

The last two rows are the same trap in a second place. Optics literature writes
`n = 1 - δ + iβ`, and both letters are already spoken for here: `δ` is the blaze angle
and `β` is the diffracted-order angle, everywhere. `fresnel.reflectivity` is called a few
lines from `facet_graze(gamma, blaze_angle, alpha)`, and two meanings of `δ` in one call
stack is how a sign error gets written and not noticed. So the materials layer uses
unambiguous words.

**Bare `d` never appears in the code.** Not as a variable, not as a keyword argument, not
in a docstring formula without an accompanying definition. When transcribing an equation
from either source, the first step is to rewrite `d` into `period` or `depth` explicitly.

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

**Note the `sin γ` in the phase shift.** The thesis expression at Appendix-D.tex:651
omits it and is inconsistent with the thesis's own Φ_b at line 672, which carries it.
The thesis form is the in-plane special case, valid only under the substitution
`λ → λ csc γ` introduced at line 446. **We use the ISSI form, which is general.**

The finite-N factor is retained deliberately. Appendix D drops it in the `N → ∞` limit,
but it is the link between efficiency and resolving power and belongs in the API.

### Scalar validity guards

Reported alongside any scalar result:

- **Reduced wavelength ratio:** `λ/(period · sin γ) ≲ 0.1`. The smallness scalar
  theory needs is judged at the reduced wavelength `λ/sin γ` — the wavelength of
  the transverse problem the grooves pose (the same decoupling the integral
  method runs on) — not at `λ/period`, which passes extreme off-plane mounts
  that are far outside the regime. In-plane the two coincide.
- **Total external reflection** (ISSI §4): `θ_c ≈ √(2δ_n)`; the model is only
  meaningful for `ζ < θ_c`.
- **Fraunhofer smoothness** (ISSI §4): `λ > 32 · sin(ζ) · σ`.

A scalar result outside these bounds is returned with a warning flag on its provenance
record, not suppressed. Mapping where scalar theory breaks down is a deliverable, so the
out-of-bounds values must remain inspectable.

---

## 10. Deviations from source references

Recorded so nobody "fixes" a deliberate choice back to a source that is wrong.

| # | Source | Issue | Our choice |
|---|---|---|---|
| 1 | thesis Appendix-D.tex:651 | `Φ_g` missing `sin γ` | ISSI eq. (15) form |
| 2 | thesis Appendix-D.tex:418, 420 | obliquity factor `cos β/cos α`, plus renormalization so `Σ E = 1` | renormalization **removed**, not optional. The obliquity factor is removed *in that form* — asymmetric, so it breaks reciprocity — and replaced by the symmetric `O_m` of §5, which perturbation theory requires |
| 2a | ISSI §2.1, groove function | `E_m = ‖G_m‖²` with no flux factor at all | wrong by `1/O_m`, up to 64% away from Littrow. The one place this project's normative reference is the one in error — see `theory/scalar.md` §3 for the evidence, which is independent of both sources |
| 3 | thesis Appendix-D.tex:699 | `∂Φ_b/∂α` missing factor `period` (dimensionally required) | corrected |
| 4 | thesis Appendix-D.tex:462, 622 | `for n = 0, ±1, …` on expressions singular at `n = 0` | `G₀` stated separately: `A₀·𝒲/period` (square), `A₀/2` (sawtooth) |
| 5 | thesis Chapter-2.tex:1016, 1045 | claims PCGrate *internally multiplies* by Fresnel reflectivity in perfect-conductivity mode — contradicts the `ΣE = 1` result proved just above | Treat `ΣE ≈ R_F(ζ)` as an emergent result of an impedance-type (Leontovich) boundary condition. Follow Goray & Schmidt (2010) for the rigorous finite-conductivity case |
| 6 | thesis Chapter-2.tex:1006–1007 | `A_n'` should be `A_n''`; evanescent exponential missing its `y` | corrected |
| 7 | thesis Chapter-2.tex:905 vs 1006 | inconsistent `k_y` sign convention | §6 above |

Items 3, 4, 6 are transcription errors with no effect on published results. Items 1, 2,
5, 7 change what the code computes.

---

## 11. Adapters

External formats are converted at the boundary. An adapter is the *only* place a foreign
convention may appear.

| Format | Direction | Notes |
|---|---|---|
| PCGrate-SX `.txt` efficiency tables | import | 7 header lines; tab-separated; quoted `"Eff.TE(m,R)"` headers carry the order index and R/T flag; `--` marks missing/evanescent. Order mapping comes from the headers, never from a column offset |
| CXRO / Henke optical constants | import | `(energy_eV, δ, β)` → `n = 1 − δ + iβ` on a wavelength grid |
| PCGrate `.ari` refractive index | export | `wavelength_nm  δ  β`, space-delimited |
| PCGrate `.ggp` polygonal border | import + export | Normalised: `t` on [0,1], `y` as a fraction of the period. **Carries no period, no provenance, no version.** Three header variants exist on disk; all are read, only the canonical one is written. A measured profile should reach a solver through `BoundaryProfile.to_problem` instead, which keeps the period |
| Nanoscope `.spm` (AFM) | import | Binary. Heights in metres. The `Sens. Zsens` / `Sens. ZsensSens` distinction is load-bearing — the wrong one gives heights 5.13× too small |
| Gwyddion text export (AFM) | import | 4 header lines, tab-separated, heights in metres, scan width scraped from the header comments |
| GSolver, RETICOLO | planned | |

Round-tripping a `Problem` through an adapter and back must leave it unchanged. That is
a test, not an aspiration — it catches sign and normalization errors before any physics
is written.
