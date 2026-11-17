# The coating as a flat mirror

Implemented in [`mirror.py`](../../src/gratinglab/mirror.py), and shown by the
Mirror tab. Symbols and sign conventions follow
[`conventions.md`](../conventions.md). This page says what the calculation is,
what its three outputs mean, and what has been checked; the section numbers
are for citing.

It is not a grating solver. It answers a question that comes before one: what
does the coating do on its own, and is its top layer thick enough that the
layers under it can be left out of the grating calculation?

---

## 1. What is modelled

A stack of flat, parallel, homogeneous layers under vacuum, on a substrate
that is taken as infinitely thick. The stack is the one a grating `Problem`
carries: `coating` (the dock's **Grating** material) is the substrate here,
and `layers` sit on it, innermost first, so the last layer is the surface the
light meets.

Angles are grazing, measured from the surface, as in
[`materials/fresnel.py`](../../src/gratinglab/materials/fresnel.py). The
polarizations are `s` and `p` with respect to the mirror's own plane of
incidence, which are not the groove-referenced TE and TM of `conventions.md`
§7; a flat mirror has no grooves to refer to.

What a flat stack leaves out is everything the grooves do. On a grating a
layer is conformal: it follows the groove and rounds its convex corners, and
the light meets each facet at its own angle. Only the
[RCWA](rcwa.md) models that (its §6). The mirror tells you about the material
stack at one angle at a time, which is why the Mirror tab marks the angles
your mount actually uses (§5).

## 2. Reflectivity: Parratt's recursion

In each medium `j` the wave has a normal wavevector component, divided by
`k = 2π/λ`,

```
q_j = sqrt(n_j² − cos²θ)
```

with the root that has a non-negative imaginary part, so the wave decays
going into the material (the `exp(−iωt)` convention of `conventions.md` §2).
At the interface between `j` and the medium below it the Fresnel coefficients
are

```
r_s = (q_j − q_j+1) / (q_j + q_j+1)
r_p = (n_j+1² q_j − n_j² q_j+1) / (n_j+1² q_j + n_j² q_j+1)
```

and the reflection seen from inside medium `j`, looking down at everything
beneath it, is built from the substrate upward:

```
X_j = (r + X_j+1 e^{2ik q_j+1 d_j+1}) / (1 + r X_j+1 e^{2ik q_j+1 d_j+1})
```

with `X = 0` in the substrate and `d` a layer's thickness. `X` in the vacuum
is the mirror's amplitude reflection coefficient, and the reflectivity is its
modulus squared; unpolarized is the mean of `s` and `p`.

The exponential only ever decays, since `Im q ≥ 0`, so the recursion is
stable for a layer of any thickness, an absorbing one many attenuation
lengths deep included. That is the same reason the RCWA propagates with a
scattering matrix (its §4).

Roughness is one RMS value `σ` applied at every interface as the
Névot–Croce factor on the amplitude,

```
r → r · exp(−2 k² q_j q_j+1 σ²)
```

which uses the wave on both sides, so it damps little below the critical
angle, where the transmitted wave is evanescent. One `σ` for every interface
is a simplification: a real stack's buried interfaces need not be as rough as
its surface.

## 3. How deep the light goes

`Mirror.penetration_depth` is, for each material in the stack taken as bulk,
the depth at which the transmitted wave's intensity has fallen to `1/e`:

```
depth = 1 / (2 k Im q)
```

Below the material's critical angle the wave is evanescent and the depth is a
few nanometres, set by the angle and nearly independent of absorption. Above
it the wave propagates and the depth is set by absorption, and is longer by
one or two orders of magnitude. The jump between the two is the critical
angle. A material with no absorption has an infinite depth above its critical
angle, and the function returns that.

This is a property of one material at one angle and wavelength. It is drawn
for every material in the stack, not only the top one, because a layer's own
depth says how far light gets into it, while the substrate's says whether
light that reaches it comes back.

## 4. When a layer may be taken as bulk

The practical question is whether a grating calculation needs the stack, or
may take the top layer as the grating material alone. The `1/e` depth
suggests an answer and does not give it: light reflected from the interface
under the layer interferes with light reflected from its surface, so the
reflectivity approaches the bulk value through fringes, not as a plain
exponential.

`Mirror.bulk_thickness` therefore sweeps the top layer's thickness, with the
real layers and substrate beneath it, and returns the smallest thickness
beyond which the reflectivity stays within a relative `tolerance` of that
layer's own bulk reflectivity. "Stays" is the point: a thickness where a
fringe happens to cross the bulk value does not count. The sweep steps finely
enough to follow both the fringes and the decay, runs to forty `1/e` depths,
and returns infinity if the reflectivity has not settled by then, as for a
lossless layer above its critical angle.

The Mirror tab compares that thickness with the layer as deposited, at every
point of the scan, and reports one of three things: the layer reads as bulk
across the scan, at some of it, or nowhere. For 10 nm of gold on silicon at
1 nm, with the typical-material model, the gold reads as bulk to one part in
a thousand only from 0.1° to 0.79° graze. That is a statement about a flat
mirror at that tolerance, and the tolerance is yours to set.

## 5. What it means for a grating

The mirror's angle is a graze angle onto a flat surface. A grating meets the
light at two such angles, and the Mirror tab draws both on an angle scan:

- onto the **grating plane**, `sin ζ = sin γ cos α`, and
- onto the **blaze facet**, `sin ζ = sin γ cos(δ − α)`, for a blazed profile
  ([`geometry.facet_graze`](../../src/gratinglab/geometry.py)).

If the top layer reads as bulk at those angles, across the wavelengths of
interest, then taking it as the grating material alone is supported for the
reflectivity. It is evidence, not proof, for the efficiency: a conformal
layer on a groove is thinner along some directions than its nominal
thickness, and multiple reflections between facets see other angles. Where
the answer matters, model the layers in the RCWA and compare.

## 6. Validation evidence

Checked in [`tests/test_mirror.py`](../../tests/test_mirror.py):

- A bare substrate reproduces the single-interface amplitudes of
  `materials/fresnel.py`, in `s` and `p`, to `10⁻¹²`, from 0.3° to 45°.
- A four-medium stack agrees with a characteristic-matrix calculation (Born
  and Wolf, §1.6) written out in the test, to `10⁻¹²` at four angles. The two
  share no code.
- A film of the substrate's own material changes nothing; a thick absorbing
  film reads as its own bulk; the reflectivity stays in `[0, 1]` over angle
  and wavelength; the layer order is innermost first.
- The evanescent depth below the critical angle matches the closed form for
  a lossless material, and the bulk thickness lies where the reflectivity
  just leaves the tolerance band.

## 7. Not yet validated

- Against a measured reflectivity curve, or against another code's output.
  The checks above are of the mathematics, and of consistency with this
  project's own Fresnel module.
- Roughness in a stack. A bare substrate agrees with `fresnel.py` in
  magnitude; nothing independent checks the factor at a buried interface.
- Whether "bulk on a flat mirror" predicts "bulk on a grating". §5 gives the
  reasoning; no comparison with layered RCWA runs has been made.
