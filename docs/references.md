# References

This page lists the literature on which this implementation is based, with the
part of the project that each source serves. Where two sources disagree, the
resolution is recorded in [`conventions.md` §10](conventions.md) and summarised
in [`findings.md`](findings.md#thesis-and-issi-chapter-reconciled).

---

## Normative

**Heilmann, R. K., Huenemoerder, D. P., McCoy, J. A. & McEntaffer, R. L.**
*Diffraction Gratings for X-ray Spectroscopy.* Springer ISSI Scientific Reports
(2024). [arXiv:2409.02297](https://arxiv.org/abs/2409.02297)

§2.1 gives the conical framework: the generalized grating equation, the
wave-vector construction and the Kirchhoff scalar form. §4 gives the scalar
efficiency of the blazed reflection grating and the validity guards.

This source is normative wherever it disagrees with any other source listed
here, because it is published, peer-reviewed, and self-consistent under
Parseval.

---

## Scalar theory

**McCoy, J. A.** *Scalar Treatment of Gratings*, PhD thesis Appendix D,
Pennsylvania State University.

This appendix is the primary derivation and gives the most detailed treatment of
the closed forms for the square wave, sinusoid and sawtooth, together with the
off-plane substitution λ → λ csc γ that allows one module to serve both mounts.
Four transcription errors and one substantive difference from the ISSI chapter
are recorded in `conventions.md` §10.

**McCoy, J. A.** *On X-ray Reflection*, PhD thesis Appendix C.

This appendix covers Fresnel reflectivity, the soft-X-ray index of refraction,
penetration depth, and (in §C.1.4) surface roughness, which forms the basis for
the Névot–Croce and Debye–Waller factors in the materials layer.

**Harvey, J. E. & Pfisterer, R. N.** *Understanding diffraction grating
behavior.* Opt. Eng. **58**(8), 087105 (2019); part II, **59**(1), 017103 (2020).

These papers give a modern parametric, linear-systems treatment that handles
conical diffraction natively. It has not yet been implemented; it is the natural
basis for a second scalar backend against which the Kirchhoff one could be
compared.

**Born, M. & Wolf, E.** *Principles of Optics.*

Kirchhoff diffraction theory and the Fraunhofer limit.

---

## RCWA

**Pommet, D. A., Grann, E. B. & Moharam, M. G.** *Effects of process errors on
the diffraction characteristics of binary dielectric gratings.* Appl. Opt.
(1995).

This paper supplies benchmark cases and bears directly on profile-error
sensitivity, that is, on how much a measured deviation from an ideal profile
changes the efficiency.

**Li, L.** Fourier factorization rules.

These rules are required for the RCWA backend. Without the correct
factorization, TM metallic cases converge badly, and cross-method plots would
then misrepresent RCWA instead of revealing anything about the physics.

**Li, L.** *Use of Fourier series in the analysis of discontinuous periodic
structures.* J. Opt. Soc. Am. A **13**(9), 1870 (1996).

The paper behind the rule above, and the one `theory/rcwa.md` §3 follows: the
inverse rule for the field component normal to a discontinuity, Laurent's for
the tangential ones.

**Lalanne, P. & Morris, G. M.** *Highly improved convergence of the
coupled-wave method for TM polarization.* J. Opt. Soc. Am. A **13**(4), 779
(1996). **Granet, G. & Guizal, B.** *Efficient implementation of the
coupled-wave method for metallic lamellar gratings in TM polarization.* J. Opt.
Soc. Am. A **13**(5), 1019 (1996).

The in-plane TM formulation that Li's rules reduce to.

**Rumpf, R. C.** *Improved formulation of scattering matrices for
semi-analytical methods that is consistent with convention.* Prog.
Electromagn. Res. B **35**, 241 (2011).

The symmetric S-matrix in a zero-thickness gap medium that the backend uses for
propagation (`theory/rcwa.md` §4).

**Popov, E., Nevière, M., Gralak, B. & Tayeb, G.** *Staircase approximation
validity for arbitrary-shaped gratings.* J. Opt. Soc. Am. A **19**(1), 33
(2002).

Why a staircased metallic profile converges slowly in TM, and the reason the
backend warns when it staircases one (`theory/rcwa.md` §5).

---

## C-method

**Dusséaux, R., Faure, C. & Chandezon, J.** *New perturbation theory of
diffraction gratings and its application to the study of ghosts.* J. Opt. Soc.
Am. A **12**(6), 1271 (1995).

This paper comes from the group that originated the coordinate-transformation
method.

**Breidne, M. & Maystre, D.** *Variational theory of diffraction gratings and
its application to the study of ghosts.* J. Opt. Soc. Am. **72**(4), 499 (1982).

---

## Integral method

**Maystre, D. & Popov, E.** *Integral Method for Gratings*, ch. 4 in E. Popov
(ed.), *Gratings: Theory and Numeric Applications*, Institut Fresnel / CNRS /
AMU (2012). Freely available.

This chapter served as the implementation guide, and the solver now follows it.
It runs to 59 pages, including §4.5 on conical mounting and §4.6 on
the numerical tooling, and the perfect-conductivity solver follows it directly.
The anchors used are eqs. (4.31)–(4.34) for TE, (4.37)–(4.42) for TM, (4.65)
for the conical invariance theorem, (4.70)–(4.74) for the equal-arc rectangular
rule, (4.82)–(4.87) for the Kummer kernel acceleration, (4.95)–(4.100) for the
log-singularity split, and Tables 4.1/4.2 as literature benchmarks (Table 4.1
is reproduced in `tests/test_integral_core.py`). Two caveats apply before the
solver is extended: §4.6.5 on edges, which accounts for the measured
first-order TM corner convergence, and §4.6.6 on graded meshes, which is the
fix. The printed sgn/sign placement in eqs. (4.33)/(4.40)/(4.41) did not
survive PDF text extraction, so the implementation derives the kernels from the
spectral form and pins every sign with the flat-mirror and Table 4.1 tests.

**Goray, L. I. & Schmidt, G.** *Solving conical diffraction grating problems
with integral equations.* J. Opt. Soc. Am. A **27**, 585 (2010).

The combination of conical mounting and integral equations is exactly the
off-plane X-ray case. Goray is the author of the standard commercial
implementation, and Schmidt supplied the rigorous analysis. The method was
implemented at M20 (`conductivity="tabulated"`;
[`theory/integral.md`](theory/integral.md) §8) with the following anchors:
eq. (6) for the conical jump conditions, eqs. (15)–(21) for the potential
operators, eq. (23) for the coupled system, §2.C for the energy balance,
eq. (26) for the absorption integral, and the numerical-option note in §3 that
`V⁺∂ₜV⁻` may be discretised by numerical differentiation (the option taken here,
via the spectral `D_t` matrix). The coupled system was re-derived in project
conventions, because their normal points into the metal and their potentials
carry a factor 2, which makes transcription a sign trap; the flat-interface
Fresnel tests pin every sign instead. Their comparison Tables 3/5/6 are usable
as validation data. Table 4 is a publisher's duplication of Table 3 (see
[`findings.md`](findings.md)). The incident state of Table 2 requires Li's
polarization-angle convention, in which δ is measured so that tan δ = |E_s/E_p|
and ψ is the relative phase (decoded from the Table 3/4 footnotes). Its deep
lamellar profile is corner-limited on the current mesh, so that case is left to
the graded-mesh milestone.

**Li, L.** *Multilayer modal method for diffraction gratings of arbitrary
profile, depth, and permittivity* (and the conical-mount companion data).
J. Opt. Soc. Am. A **10**, 2581 (1993).

This is the independent cross-method source behind the comparison tables of
Goray & Schmidt, and it is what makes their Table 3 a genuine anchor between two
codes instead of a self-comparison.

**Maystre, D.** *Analytic Properties of Diffraction Gratings*, ch. 2 of the same
volume.

This chapter treats anomalies and analytic structure, and is directly relevant
to the Rayleigh-anomaly instability documented in
[`findings.md`](findings.md#a-reference-file-violates-energy-conservation).

**Kalhor, H. A. & Neureuther, A. R.** *Effects of conductivity, groove shape,
and physical phenomena on the design of diffraction gratings.* J. Opt. Soc. Am.
**62** (1972).

This is an early integral-equation treatment and is useful for the
perfectly-conducting limit.

**McCoy, J. A.** *Modeling Diffraction Efficiency*, PhD thesis Chapter 2 §2.2.

This section gives the same derivation in the notation of this project, proceeding
from the boundary value problem to the Helmholtz equation with boundary
conditions, then to the Green's function, Dirichlet enforcement and the
perfectly-conducting limit, with Σℰ = 1 proved via Green's theorem. §2.1
documents the measurement method used at ALS beamline 6.3.2.

> §2.2.3 is outdated and should not be ported as written. It describes the
> commercial code as internally multiplying by Fresnel reflectivity in
> perfect-conductivity mode, which contradicts the Σℰ = 1 result proved
> immediately before it. See `conventions.md` §10, item 5.

---

## Validation

**Moharam, M. G. & Gaylord, T. K.** Canonical RCWA test cases.

**Li, L.** Crossed-grating benchmarks.

Neither is yet used: the papers (Moharam, Grann, Pommet & Gaylord, JOSA A 12,
1068 (1995); Li, JOSA A 13, 1870 (1996); likewise Lalanne & Morris, JOSA A 13,
779 (1996) and Granet & Guizal, JOSA A 13, 1019 (1996)) are closed access, and
no value is transcribed from memory. The RCWA benchmarks in use come from:

**Lalanne, P. & Hugonin, J.-P.** *Numerical performance of finite-difference
modal methods for the electromagnetic analysis of one-dimensional lamellar
gratings.* J. Opt. Soc. Am. A **17**(6), 1033-1042 (2000). Open copy on HAL,
hal-00867708. Fig. 3 (Granet's metallic case, printed as 84.848% TM) and Fig. 6
(a dielectric case, 51.062%).

**Popov, E. (ed.)** *Gratings: Theory and Numeric Applications*, 2nd revisited
edition. Institut Fresnel (2014). Open copy on HAL, hal-01084458. Ch. 5 Tab.
5.1 prints Granet's metallic case to seven digits; ch. 10 §10.6 (Gralak) gives
the silicon-rod case.

**Granet, G.** *Reformulation of the lamellar grating problem through the
concept of adaptive spatial resolution.* J. Opt. Soc. Am. A **16**, 2510-2516
(1999). The origin of the metallic values; closed access, so used only as the
two sources above print them.

These cases lie in the visible and near-IR, which is why that regime serves as
the project's *validation* regime even though the application is in the soft
X-ray. See [`roadmap.md`](roadmap.md#regimes-and-the-methods-that-serve-them).

**Heuberger, G., Klepp, J., Guo, J., Tomita, Y. & Fally, M.** *Light diffraction
from a phase grating at oblique incidence in the intermediate diffraction
regime.* Appl. Phys. B **127**, 72 (2021).

This paper provides measured data on where scalar theory ceases to hold, and is
the empirical basis for a map of scalar validity.

---

## Optical constants

No optical constants ship with the package; these are the sources whose formats
its readers are written for.

**Henke, B. L., Gullikson, E. M. & Davis, J. C.** *X-ray interactions:
photoabsorption, scattering, transmission, and reflection at E = 50–30000 eV,
Z = 1–92.* At. Data Nucl. Data Tables **54**(2), 181–342 (1993).

This is the compilation behind the CXRO tables, whose format the X-ray reader
reads.

**Polyanskiy, M. N.** *Refractiveindex.info database of optical constants.*
Sci. Data **11**, 94 (2024).

This database defines the format of the ultraviolet-to-infrared reader and its
nine dispersion formulas. The database is CC0.

**Malitson, I. H.** *Interspecimen comparison of the refractive index of fused
silica.* J. Opt. Soc. Am. **55**(10), 1205–1209 (1965).

This is the source for the fused silica that the reader's tests hold to its
published value of n_d = 1.4585, a value that neither the reader nor the
database's own code had produced.

---

## Application context

**Tutt, J. et al.** Diffraction efficiency testing of sinusoidal and blazed
off-plane reflection gratings. J. Astron. Instrum. (2016).

**Marlowe, H. et al.** Polarization dependence at grazing incidence. This work
justifies treating soft-X-ray reflectivity as polarization-independent, which is
why the reference corpus carries only TE.

**McCoy, J. A. et al.** J. Vac. Sci. Technol. B (2018); ApJ **891**, 114
(2020); OSA Continuum **3**(11), 3141 (2020).

---

## Metrology

**Villarrubia, J. S.** *Algorithms for scanned probe microscope image
simulation, surface reconstruction, and tip estimation.* J. Res. Natl. Inst.
Stand. Technol. **102**(4), 425–454 (1997).

This paper underlies the tip-deconvolution milestone (roadmap §4): morphological
dilation and erosion, the certainty map that marks unrecoverable points, and
blind tip estimation. Erosion bounds the true surface but does not recover it.

---

## Sibling code

**endstation** ([github.com/jamccoy/endstation](https://github.com/jamccoy/endstation),
on PyPI). This package is shared scaffolding for the author's desktop science
tools and is gratinglab's one non-scientific dependency. It provides the
optional-dependency guards behind every missing-extra message, form-error
collection, the project and data-root lookups, the metrology wiki's page index,
progress throttling and cancellation, stale-result tokens, the guard that keeps
garbage collection on the GUI thread, Qt layout sizing, the toolkit-boundary
test and pytest plugin, and the icon, macOS app and wiki-publishing tools. Its core declares no
dependencies, so it adds nothing to a solver-only install. Fixes to any of these
components belong in endstation, not in this repository.

**PyXFocus** is a sequential X-ray raytracer (Allured; forks by kbuffo and by
this project's author). Its `grating.py` holds the analytic Littrow-mount
formulas (`blazeYaw`, `blazeAngle`) that are slated for a native port
(roadmap §3), and the order-centroid and spot-width estimator in its GUI serves
as the instrument-level cross-check of `gratinglab.resolution`. The raytracer
itself is deliberately not absorbed into this project.
