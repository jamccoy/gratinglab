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

`Problem` has no superstrate field. The mapping of a bare grating is the one the
integral solver's tabulated mode uses, so that the same `Problem` means the
same grating to both:

| | |
|---|---|
| superstrate | vacuum, which is the medium the incidence angles are defined in |
| grating region | ridges of the named material, grooves filled with vacuum |
| substrate | the same material, semi-infinite, below the bottom of the grooves |
| material | `problem.coating`, falling back to `substrate` (as the integral solver does); naming neither is refused |
| layers | `problem.layers`, conformal, on top of the grating (§6) |

Efficiencies are **absolute** (conventions §5). The method needs a complex
index; there is no perfect-conductor mode, because a perfect conductor has no
permittivity to Fourier-transform. Comparisons with the integral solver use its
`conductivity="tabulated"` mode.

Where `coating` and `substrate` are both named and differ, the substrate is not
consulted and a provenance warning says so. A coating of some thickness on the
grating is a conformal layer (§6), not a second material name.

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

### Transmission, absorption, and what `R + T + A = 1` does and does not show

Into a lossless substrate the transmitted orders are recorded per order as the
scan's `transmission` (`result.Transmission`), on an order axis of their own,
since the substrate admits more orders than the vacuum. They are formed like
the reflected ones, from `kz |E|²` of the transmitted amplitudes, and every
order the substrate admits must fit inside the harmonic window or the solve is
refused. An absorbing substrate has no transmitted orders: every wave in it
decays, so there is no far field to call an order.

The scan's `absorption` is what is *absorbed*. It starts from the power
crossing into the structure through the plane just below the vacuum,
evaluated from the field inside the stack (the amplitudes at that plane come
from the partial S-matrices above and below it, not from the reflected
amplitudes), and subtracts the transmitted orders. So it is the power absorbed
in the ridges, plus, for an absorbing substrate, the power entering the
substrate. `check_energy_balance(lossless=True)` then checks `R + T + A = 1`.

That check is weaker here than for the integral solver, and the difference
matters. There, the absorption is an independent boundary integral and the sum
deviates from 1 by the discretisation error. In RCWA the truncated system
conserves flux by construction, because tangential `E` and `H` are continuous
across every interface the S-matrix joins, so the sum is 1 to rounding at
*every* truncation, converged or not. It is a test of the implementation: of
signs, normalisation and the S-matrix composition. It is not an error estimate.
Convergence is shown by the truncation sweep alone, and the provenance notes say
so.

### Cost, and the two shortcuts

Per patterned slab the work is one eigenproblem of size `2N` and a handful
of `2N`-square products and inverses, `N = 2M + 1`. Two cases need less, and
the core takes both.

**In-plane, the system splits.** With `ky = 0`, `PQ` is block diagonal, so
every mode is pure `Ex` (TM) or pure `Ey` (TE), and `Q` and every `V` are
anti-diagonal, pairing `Ex` with `Uy` and `Ey` with `Ux`. Every product the
S-matrix algebra forms, `V⁻¹V0` and `V0⁻¹V`, is then block diagonal, so each
polarization runs the same code on its own `N`-square blocks: an eigenproblem
an eighth of the cost, products a quarter, and a polarization that nothing
excites is not solved at all. `decouple=False` on `solve_stack` keeps the
coupled solve, which is how the tests hold the two together (to `1e-12`).

**Only the top plane is needed.** The backend records `absorption`, which is
the flux through the first plane, and that needs the partial products from
the bottom and nothing from the top. `planes="top"` skips the top-down
products and the solve at every interior plane. The core's default stays
`planes="all"`, since a film's tests want every plane.

BLAS threading matters more than either. At these sizes OpenBLAS's threads
cost more than they compute, and on a busy machine a 122-square inverse took
a tenth of a second; the backend caps BLAS at one thread for the call
(`solvers/_execution.py`) and parallelises over wavelengths instead
(`workers`). The measurements are in the changelog and
`benchmarks/solver_speed.py` reproduces them.

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

### Conformal layers, through `Problem.layers`

A coated grating is `Problem.layers`, a list of `ConformalLayer(material,
thickness)` innermost first, each following the surface beneath it as its true
offset (conventions §3, decided 2026-10-27): convex corners rounded into arcs
of radius `thickness`, concave ones sharp. Every other solver refuses a problem
with layers (`Capabilities.handles_layers`), rather than solving the bare
grating beneath them.

Geometrically, surface `k` is the profile offset by the first `k` thicknesses,
and since offsets by disks add, that is one offset of the profile, not a chain.
Its height over `x` is computed exactly: for each segment of the profile's
polyline, the top of its stadium of radius `thickness` is concave along the
segment and so maximised in closed form, and the surface is the largest over
segments. Lamellar and blazed profiles use their exact vertices, so a vertical
wall stays a vertical segment; smooth profiles are sampled at 4096 points.

Each slab then takes, at its mid-height, the material of the innermost surface
it lies under, with every material boundary refined by bisection as the bare
walls are. Slab boundaries are placed at every height where some surface has a
flat (the profile's top and bottom, and their offsets), so on a lamellar groove
the flats are exact and only the rounded corners are staircased. A layered
grating is therefore always a staircase, `slices` applies, and the TM warning of
§5 applies to it.

### Homogeneous layers in the core

The core's `Stack` also takes homogeneous films above and below the grating
region, tested against the closed-form Airy formula for a single film (§8). No
`Problem` produces those yet: a planar film under a grating is not a conformal
layer, and nothing has needed one.

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

Everything here runs without the private corpus: `tests/test_rcwa.py`,
`tests/test_rcwa_benchmarks.py`, `tests/test_rayleigh_anomaly.py` and
`tests/test_conditioning.py`.

| Reference | What it pins | Independent of this code? |
|---|---|---|
| flat interface against `materials.fresnel`, in-plane, TE and TM, glass and a visible metal | signs and normalisation in the external regions | yes |
| flat interface, conical, unpolarized against `(Rs + Rp)/2` at the true angle of incidence | the conical kinematics, and that TE/TM by the grooves span the polarizations | yes |
| a single film on a substrate against the Airy formula, s and p, in-plane and conical | the core's homogeneous layers | yes |
| a conformal layer of the superstrate's own index against the bare grating, to 1e-12 | the layered geometry's bookkeeping | yes |
| a conformal layer of the grating's own material against the integral solver on the hand-written offset groove: 2e-5 dielectric, 2e-4 metal TE | the offset geometry and its staircase | yes |
| a conformal film on a vanishing groove against the Airy formula | layers in the flat limit | yes |
| lossless dielectric lamellar, `R + A = 1` with `A` the transmitted power, in-plane and conical | the S-matrix composition and flux bookkeeping (§4: not convergence) | partly |
| reciprocity, `checks.check_reciprocity`, metallic lamellar in TE and TM | the symmetry of the scattering problem | yes |
| truncation sweep, metallic lamellar in TM: Li's rule converges, Laurent's visibly does not | that §3 is in effect | yes |
| the integral solver (`conductivity="tabulated"`, corner-graded) on lamellar glass and metal, TE and TM, in-plane and conical | the whole method, against a different method | yes |
| **published**: Granet's metallic lamellar, TM, `R0 = 0.8484781` | metallic TM through `Problem` and the convergence harness | yes |
| **published**: Lalanne & Hugonin's dielectric lamellar on glass, TM, `T-1 = 51.062%` | transmitted orders, a substrate unlike the ridges | yes |
| **published**: Gralak's free-standing silicon rods, p, `R0 = 0.9487` | a high-index dielectric with vacuum on both sides | yes |
| **published**: Goray & Schmidt's dielectric lamellar, conical, circularly polarized, all eight orders | conical incidence, an arbitrary polarization state, transmitted orders | yes |
| **published**: Goray & Schmidt's metallic lamellar, conical, linear at 45° | a conical metal; agrees to 1.1e-3, with the gap traced to the table (§9) | yes |
| a sweep through a Rayleigh anomaly, to 1e-8 nm either side and on it | the external `kz -> 0`: finite, balanced, continuous | yes |

The measured numbers are in the test docstrings, beside the assertions they
justify.

### The published benchmarks

Each value was read in a primary text that is free to read, and is cited in
the test with its page. Measured 2026-10-27, error against the published value
by truncation:

| case | source | M = 20 | 40 | 80 | asserted |
|---|---|---|---|---|---|
| metal `n = 0.22 + 6.71i`, period = wavelength = depth, fill 0.5, 30°, TM `R0` | Granet 1999, as printed by Lalanne & Hugonin, JOSA A 17, 1033 (2000), Fig. 3, and by Popov (ed.) 2014, Tab. 5.1 | -4.2e-3 | -1.7e-3 | -6.3e-4 | 1e-3 at M = 80 |
| ridges of 2.3 on 1.5, period two wavelengths, fill 0.234, 30°, TM `T-1` | Lalanne & Hugonin 2000, Fig. 6 | +1.2e-3 | +2.0e-4 | +1.4e-5 | 1e-4 at M = 60 |
| silicon rods `eps = 12.96`, width 0.28 d, height `d/(2 sqrt 2)`, 45°, p `R0` | Gralak, ch. 10 of Popov (ed.) 2014, §10.6 | -2.2e-3 | -3.2e-4 | -6.6e-5 | 1e-4 at M = 80 |
| `eps = 2.25` lamellar, `2H = d/2`, fill 0.5, `lambda = d/2`, `theta = 35.264°`, `phi = 30°`, `delta = 45°`, `psi = 90°`, `R-2` (the slowest of eight orders) | Goray & Schmidt, ch. 12 of Popov (ed.) 2014, Table 12.2 | +3.3e-6 | +8.1e-7 | +4.1e-7 | every order to half its last printed digit at M = 80 |
| `eps = -24.99 + 1i` lamellar, `2H = d`, fill 0.5, `lambda = d/2`, `theta = 22.208°`, `phi = 20.705°`, `delta = 45°`, `psi = 0`, `R0` | the same, Table 12.1 | +5.0e-3 | -4.3e-4 | -1.1e-3 | 1.5e-3 at M = 80 (§9) |

The metallic TM case converges from below and slowly (-2.1e-4 still at
M = 200), as §3 leads one to expect. Lalanne & Hugonin's own RCWA column reads
0.84677 at 81 orders, the same error as this backend at M = 40, which is 81
orders. The dielectric case fixed two conventions its paper leaves to a figure
that could not be read: the fill is the ridge's, and the order is this code's
`-1`. The other three readings miss by more than 0.4, so the match identifies
the convention rather than being tuned to it.

Goray & Schmidt's two conical tables are the first published check of a
conical mount and of a polarization that is neither TE nor TM. Their chapter
defines everything needed, and nothing was guessed: the frame (§12.2.1,
Fig. 12.2) is conventions §3's, with `y` the normal into the cover and `z`
along the grooves; the incident wave vector is
`k (sin θ cos φ, −cos θ cos φ, sin φ)`, so `Kx0 = sin θ cos φ` and
`Ky = sin φ` here; the time factor is `e^{−iωt}`; orders are
`α_n = α + 2πn/d` (§12.2.2), this code's `m`; and the polarization basis
(§12.8.1) is `s = k × ŷ/|k × ŷ|`, `p = s × k/|k|`, with
`δ = arctan|E_s/E_p|` and `ψ = −arg(E_s/E_p)`. The incident field is built
from those and passed to `solve_stack` as it stands, since `Illumination`
carries only TE and TM. Each grating is one profile between the cover and the
substrate, so its ridges are the substrate's material, and with fill 0.5 which
of ridge or groove `c` measures moves only the orders' phases.

The dielectric table (12.2) is matched in all eight orders, three reflected and
five transmitted, to half a unit of the last printed digit from M = 80, with
the residuals no longer moving: what is left is the table's rounding. The
opposite handedness misses every order by more than a hundred such units, so
the match identifies the polarization convention. The metallic table (12.1) is
matched only to 1.1e-3, which §9 explains.

The canonical papers are **not** used. Moharam, Grann, Pommet & Gaylord, JOSA A
12, 1068 (1995); Lalanne & Morris, JOSA A 13, 779 (1996); Li, JOSA A 13, 1870
(1996); and Granet & Guizal, JOSA A 13, 1019 (1996) are all closed access, and
no number from them could be confirmed in the original. Granet 1999 is closed
too, so its values enter only as quoted by the two sources above, which agree
with each other. The two published values that do **not** match are in §9.

### Rayleigh anomalies

Where an order passes off, its `kz` goes to zero, and the external region's
`V = Q kz^-1` with it; `_KZ_FLOOR` keeps the exact point finite. The physics is
a square-root branch point, so the gap between the rows at `lambda* - d` and
`lambda* + d` should be `C sqrt(d / lambda*)` with `C` constant. Measured on a
metallic lamellar at M = 20, `C` is 0.96 in TE and 0.86 in TM, constant to 4%
from `d = 1e-2` to `1e-8` nm, and the tabulated integral solver finds 0.96 and
0.95 on the same grating. `R + A = 1` holds to rounding at every row, and the
exact point sits within 1.6e-6 of its neighbours. A real step would make `C`
grow a thousandfold across the sweep.

### Conditioning

Every dense solve and inverse in `_core` goes through
`solvers/_conditioning.py`: one LU factorisation serves both the solve and
LAPACK's `gecon` estimate of the 1-norm condition number, which costs `O(n^2)`
on top. The solves are numpy's to rounding. The worst estimate in a scan goes
to `Provenance.notes["condition"]` in two groups: `modes` (each slab's `W` and
`V`) and `smatrix` (each slab's `A` and `D`, the external connections' `A`, the
star products' `I - S S`, and the plane-flux solves). Measured 2026-10-27,
in-plane at 500 nm on period 600 unless stated:

| case | M | `modes` | `smatrix` |
|---|---|---|---|
| glass lamellar, depth 180 nm or 3 um | 20 | 3.7e2 | 7.4e3 |
| | 100 | 1.3e4 | 1.9e6 |
| metal `0.2 + 3i`, either depth | 20 | 9.2e2 | 2.2e4 |
| | 100 | 2.6e4 | 4.3e6 |
| | 200 | 1.3e5 | 4.4e7 |
| metal `0.03 + 20i` | 60 | 6.2e4 | 1.7e6 |
| a slab within 1e-11 of homogeneous, normal and conical incidence | 10 | 8.4e1 | 8.1e2 |
| metal, an order 1e-12 nm from grazing | 20 | 2.3e3 | 9.0e6 |
| metal, an order exactly at grazing (`kz` on its floor) | 20 | 2.3e3 | 5.0e11 |

Depth changes nothing, which is the S-matrix doing its job (§4), and neither
does the near-degenerate slab, where `numpy.linalg.eig` might have returned
nearly parallel eigenvectors. The estimate grows as `M^2` to `M^3`. The
warning threshold is `1e10`, more than two decades above every ordinary case,
and the one case past it is the exact anomaly. The answer there was still right
(above), because the near-singular direction is the grazing order's own
amplitude, so the warning means "unusual", not "wrong".

The table was measured on the coupled `2N` system. Since in-plane mounts are
solved as two decoupled systems, in-plane TE no longer meets that singularity
at all: its external region's `V` is `(kx^2 - eps) / kz = -kz`, which goes to
zero at grazing, so the same exact-grazing case reads 27. TM's `V` is
`eps / kz` and still crosses the threshold there, which is what
`tests/test_conditioning.py` now holds; conical mounts are still coupled.

## 9. Not yet validated, and decisions to make

Not yet validated:

- staircase accuracy for blazed and sinusoidal profiles, in TM above all (§5);
- conformal layers on a metal in TM: they converge, slowly, toward the integral
  solver (1.1e-2 at truncation 40 and 160 slices on a 40 nm layer), the
  staircase weakness of §5 met at the rounded corners;
- conformal layers of a *different* material from the grating: nothing
  independent checks them yet beyond energy balance, the vacuum-layer identity
  and the Airy limit, because no other solver here models two interfaces;
- the soft-X-ray grazing regime. It is not refused, and it is expensive: a
  groove hundreds of wavelengths deep needs very many slices and harmonics,
  which is why the integral method is the method of record there;
- the polarization angles `δ` and `ψ` of the diffracted orders, which Goray &
  Schmidt also tabulate; `solve_stack` returns powers, not order amplitudes.

Three published values this backend does **not** reproduce, recorded rather
than dropped:

- **Granet's metallic lamellar in TE** (`R0 = 0.7342789`, the same table as the
  TM benchmark in §8). This backend converges cleanly to `R0 = 0.13170` (M = 80
  to 200 agree to 8e-6), and the integral solver, a different method on the
  same `Problem`, finds 0.13070 at 400 graded nodes. The published number is
  matched instead, to 9e-6 at M = 160, by this code's **Littrow** order
  (`m = +1` here): 0.7342881. The groove is half a wavelength wide, at the TE
  slot-mode cutoff, so TE is very sensitive to the geometry; but no change of
  geometry moves TE alone, and TM's `R0` matches. Fig. 1 of Lalanne & Hugonin
  2000, which defines the parameters, could not be read, and neither could
  Granet 1999. Unresolved.
- **Gralak's silicon rods in s** (`R0 = 0.7323`, the p case's companion). This
  backend converges to 0.74230 (M = 40 to 120 agree to 6e-5), 0.0100 above it,
  and no other order is near the printed value. With vacuum on both sides the
  integral solver cannot check it. Possibly a one-digit misprint; not
  established.
- **Goray & Schmidt's metallic conical lamellar** (Table 12.1, ch. 12 of
  Popov (ed.) 2014; §8). This backend converges to `R-2, R-1, R0, R1 =
  7.5551, 13.2648, 44.1572, 31.1065%` at M = 240 against the printed
  `7.52, 13.25, 44.27, 31.05`: `R0` 1.1e-3 low, `R1` 5.7e-4 and `R-2` 3.5e-4
  high. The same state fed to the integral solver as `(E_z, B_z)`, corner
  grading 4, agrees with RCWA instead: `R0` is 44.249 at 400 points, 44.185
  at 600, 44.170 at 800 and 44.163 at 1200, and at 1200 every order is
  within 6e-5 of RCWA. At 400 points, the paper's count, it is within 2.2e-4 of the table in
  every order. So the printed values look like an integral solution not yet
  converged in its discretization, which a 1e-6 energy balance would not show
  (§4). That is an inference from this project's two solvers, not a statement
  of the chapter's; the dielectric companion table (12.2) is matched to its
  printed digits, which rules out a convention fault common to both.

Decisions for the project rather than for this backend:

1. ~~**How `Problem` describes a layer stack.**~~ Decided 2026-10-27: a coated
   grating means a layer that follows the groove, `Problem.layers` of
   `ConformalLayer`, its surface the true offset with rounded convex corners
   (§6).
2. ~~**Whether a staircase is acceptable by default.**~~ Decided 2026-10-27:
   yes, with the TM warning. RCWA can represent nothing but slabs, so the
   choice was only ever how many and when to warn.
3. ~~**Whether transmitted orders belong in `EfficiencyScan`.**~~ Decided
   2026-10-27: yes, as `EfficiencyScan.transmission`, filled by this backend
   and the integral solver (§4).
