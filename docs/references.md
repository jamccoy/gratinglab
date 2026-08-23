# References

The literature this implementation is built from, mapped to the part of the
project each one serves. Where two sources disagree, the resolution is recorded
in [`conventions.md` §10](conventions.md) and summarised in
[`findings.md`](findings.md#thesis-and-issi-chapter-reconciled).

---

## Normative

**Heilmann, R. K., Huenemoerder, D. P., McCoy, J. A. & McEntaffer, R. L.**
*Diffraction Gratings for X-ray Spectroscopy.* Springer ISSI Scientific Reports
(2024). [arXiv:2409.02297](https://arxiv.org/abs/2409.02297)

§2.1 gives the conical framework — the generalized grating equation, the
wave-vector construction, the Kirchhoff scalar form. §4 gives the blazed
reflection-grating scalar efficiency and the validity guards.

**This is normative where it disagrees with any other source here**: it is
published, peer-reviewed, and self-consistent under Parseval.

---

## Scalar theory

**McCoy, J. A.** *Scalar Treatment of Gratings*, PhD thesis Appendix D,
Pennsylvania State University.

The primary derivation, and the most detailed treatment of the closed forms:
square wave, sinusoid and sawtooth, with the off-plane substitution
λ → λ csc γ that lets one module serve both mounts. Four transcription errors
and one substantive difference from the ISSI chapter are recorded in
`conventions.md` §10.

**McCoy, J. A.** *On X-ray Reflection*, PhD thesis Appendix C.

Fresnel reflectivity, the soft-X-ray index of refraction, penetration depth, and
§C.1.4 surface roughness — the basis for the Névot–Croce and Debye–Waller
factors once the materials layer lands.

**Harvey, J. E. & Pfisterer, R. N.** *Understanding diffraction grating
behavior.* Opt. Eng. **58**(8), 087105 (2019); part II, **59**(1), 017103 (2020).

A modern parametric / linear-systems treatment that handles conical diffraction
natively. Not yet implemented; it is the natural basis for a second scalar
backend to compare against the Kirchhoff one.

**Born, M. & Wolf, E.** *Principles of Optics.*

Kirchhoff diffraction theory and the Fraunhofer limit.

---

## RCWA

**Pommet, D. A., Grann, E. B. & Moharam, M. G.** *Effects of process errors on
the diffraction characteristics of binary dielectric gratings.* Appl. Opt.
(1995).

Benchmark cases, and directly relevant to profile-error sensitivity — the
question of how much a measured deviation from an ideal profile changes
efficiency.

**Li, L.** — Fourier factorization rules.

Non-negotiable for the RCWA backend. Without the correct factorization, TM
metallic cases converge badly and cross-method plots would misrepresent RCWA
rather than reveal anything about physics.

---

## C-method

**Dusséaux, R., Faure, C. & Chandezon, J.** *New perturbation theory of
diffraction gratings and its application to the study of ghosts.* J. Opt. Soc.
Am. A **12**(6), 1271 (1995).

From the group that originated the coordinate-transformation method.

**Breidne, M. & Maystre, D.** *Variational theory of diffraction gratings and
its application to the study of ghosts.* J. Opt. Soc. Am. **72**(4), 499 (1982).

---

## Integral method

**Maystre, D. & Popov, E.** *Integral Method for Gratings*, ch. 4 in E. Popov
(ed.), *Gratings: Theory and Numeric Applications*, Institut Fresnel / CNRS /
AMU (2012). Freely available.

**The implementation guide — now the implemented guide.** 59 pages, including
§4.5 on conical mounting and §4.6 on the numerical tooling. The
perfect-conductivity solver follows it directly; the anchors actually used:
eqs. (4.31)–(4.34) TE, (4.37)–(4.42) TM, (4.65) the conical invariance
theorem, (4.70)–(4.74) the equal-arc rectangular rule, (4.82)–(4.87) the
Kummer kernel acceleration, (4.95)–(4.100) the log-singularity split, and
Tables 4.1/4.2 as literature benchmarks (Table 4.1 is reproduced in
`tests/test_integral_core.py`). Caveats worth knowing before extending:
§4.6.5 on edges (the measured first-order TM corner convergence), §4.6.6 on
graded meshes (the fix). One transcription warning: the printed sgn/sign
placement in eqs. (4.33)/(4.40)/(4.41) did not survive PDF text extraction —
the implementation derives the kernels from the spectral form and pins every
sign with the flat-mirror and Table 4.1 tests instead.

**Goray, L. I. & Schmidt, G.** *Solving conical diffraction grating problems
with integral equations.* J. Opt. Soc. Am. A (2010).

Conical plus integral equations is exactly the off-plane X-ray case. Goray is
the author of the standard commercial implementation; Schmidt supplied the
rigorous analysis.

**Maystre, D.** *Analytic Properties of Diffraction Gratings*, ch. 2 of the same
volume.

Anomalies and analytic structure — directly relevant to the Rayleigh-anomaly
instability documented in [`findings.md`](findings.md#a-reference-file-violates-energy-conservation).

**Kalhor, H. A. & Neureuther, A. R.** *Effects of conductivity, groove shape,
and physical phenomena on the design of diffraction gratings.* J. Opt. Soc. Am.
**62** (1972).

Early integral-equation treatment; useful for the perfectly-conducting limit.

**McCoy, J. A.** *Modeling Diffraction Efficiency*, PhD thesis Chapter 2 §2.2.

The same derivation in this project's notation: boundary value problem →
Helmholtz with boundary conditions → Green's function → Dirichlet enforcement →
the perfectly-conducting limit, with Σℰ = 1 proved via Green's theorem. §2.1
documents the ALS beamline 6.3.2 measurement method.

> §2.2.3 is **outdated** and should not be ported as written — it describes the
> commercial code as internally multiplying by Fresnel reflectivity in
> perfect-conductivity mode, which contradicts the Σℰ = 1 result proved
> immediately above it. See `conventions.md` §10, item 5.

---

## Validation

**Moharam, M. G. & Gaylord, T. K.** — canonical RCWA test cases.

**Li, L.** — crossed-grating benchmarks.

These live in the visible and near-IR, which is why that regime is the
project's *validation* regime even though soft X-ray is the application. See
[`roadmap.md`](roadmap.md#two-regimes-two-roles).

**Heuberger, G., Klepp, J., Guo, J., Tomita, Y. & Fally, M.** *Light diffraction
from a phase grating at oblique incidence in the intermediate diffraction
regime.* Appl. Phys. B **127**, 72 (2021).

Measured data on where scalar theory stops working — the empirical basis for a
scalar validity map.

---

## Application context

**Tutt, J. et al.** — diffraction efficiency testing of sinusoidal and blazed
off-plane reflection gratings. J. Astron. Instrum. (2016).

**Marlowe, H. et al.** — polarization dependence at grazing incidence. The
justification for treating soft-X-ray reflectivity as polarization-independent,
which is why the reference corpus carries only TE.

**McCoy, J. A. et al.** — J. Vac. Sci. Technol. B (2018); ApJ **891**, 114
(2020); OSA Continuum **3**(11), 3141 (2020).
