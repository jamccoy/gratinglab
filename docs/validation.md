# Validation inventory

**Draft, 2026-11-20. For the author to review, not yet a claim.** For each part
of gratinglab, this lists the evidence that it is correct and sorts it by how
independent that evidence is. Two uses: it shows where the author's own
reading is most needed, and it is the oversight record a JOSS AI-usage
disclosure can point at.

It was compiled by reading the tests, fixtures and docs. No test was run to
produce it. Five of its sharpest claims were checked by hand against the code
(marked ✓). The rest are as read from the tests, and §8 says what was not
looked at.

## 1. Tiers

| tier | meaning |
|---|---|
| **independent** | held by something this code did not produce: a published value or table, a closed form from a cited source, an analytic limit, a measurement, another program's output, or a solver here that uses a genuinely different method |
| **self-consistent** | checked only against itself: invariants (reciprocity, energy), convergence, round trips, synthetic data made under the code's own assumptions, or numbers the code produced and the test pinned |
| **unchecked** | nothing exercises the claim, or a test checks only that it runs or that a warning fires |

**Tautology risk** marks a test whose expected value was produced by the code
it tests.

## 2. Three caveats that apply throughout

- **The solvers share their plumbing.** All three rigorous solvers, and the
  scalar one, take `geometry.sin_beta`/`cos_beta`, `order_range`,
  `Illumination`, `Problem`/`profiles` and `materials.lookup` from the same
  modules. Every analytic oracle in the scalar tests also takes its order
  directions from `sin_beta`. So agreement between methods cannot catch an
  error in the angle convention, order bookkeeping, normalisation or material
  lookup. Only the published and analytic checks can. The in-plane grating
  equation has one textbook point (`test_geometry.py::TestGratingEquation::test_in_plane_reduces_correctly`).
  The conical λ csc γ form is checked only against itself.
- **The strongest external evidence never runs in CI.** Tests on the PCGrate
  corpus (`GRATINGLAB_REF_DIR`), the AFM archive (`GRATINGLAB_AFM_DIR`) and the
  ALS beamline tree skip without their data. CI only confirms they skip
  (`.github/workflows/test.yml`, "Confirm reference-data tests skipped").
  None of the three is on MC1. So the PCGrate comparisons, every real-scan
  metrology check and the ALS campaign tests have no recorded recent pass.
- **"Independent" is often the same author.** PyXFocus, the dissertation's
  tables, the in-test Born & Wolf transfer matrix and the SEM imaging model
  are independent of this code, but not of the person and sessions that wrote
  it. They catch a slip, not a shared misunderstanding.

## 3. Rigorous solvers

| feature | tier | evidence | notes |
|---|---|---|---|
| RCWA flat interface, homogeneous films | independent | `test_rcwa.py::TestFlatInterfaceIsFresnel`, `TestLayerStack::test_one_film_is_the_airy_formula` (1e-12) | Fresnel, Airy |
| RCWA metallic TM lamellar | independent | `test_rcwa_benchmarks.py::TestGranetMetallicTM`, 0.8484781 (Granet 1999 via L&H 2000, Popov Tab. 5.1), 1e-3 at M = 80 | |
| RCWA dielectric T−1; silicon rods p | independent | `::TestLalanneHugoninDielectricTM` 0.51062; `::TestGralakSiliconRodsTM` 0.9487; both 1e-4 | fill and order readings chosen by matching, with a negative control |
| RCWA conical dielectric | independent | `::TestGoraySchmidtDielectricConical`, Table 12.2, all eight orders, handedness | the strongest single check in the package |
| RCWA conical metal | independent, partial | `::TestGoraySchmidtMetallicConical`, Table 12.1, 1.5e-3 | tolerance covers a measured 1.1e-3 gap, attributed to the paper using this package's integral solver |
| RCWA Granet TE, Gralak s | **not reproduced** | `docs/theory/rcwa.md` §9: 0.13170 against a published 0.7342789; 0.74230 against 0.7323 | no test. The published Granet TE matches this code's Littrow order to 9e-6, which looks like a labelling or convention question |
| RCWA layers of a second material | self-consistent | vacuum identity, energy balance | rcwa.md §9 says so |
| RCWA staircased blazed or sinusoidal profiles, metal TM | unchecked for accuracy | geometry only (`::TestGrooveGeometry`) | admitted in rcwa.md §9 |
| RCWA Li factorization | self-consistent, tautology risk | `::TestLisRuleIsInEffect` against its own value at truncation 150 | |
| C-method profile series, flat interface | independent | `test_cmethod.py::TestProfileSeries`, `::TestFlatInterfaceIsFresnel` | Bessel, Fresnel |
| C-method conical finite conductivity | independent | `::TestGorayEtSchmidtTable3`, 5e-5 | Table 3 is itself Li's C-method, an external implementation |
| C-method perfect conductor TM | independent | `::test_marechal_stroke` (1e-12); Maystre & Popov Table 4.1 TM 0.9579 | |
| C-method perfect conductor TE | **tautology risk** ✓ | `test_cmethod.py:367` expects 0.466441, the code's value, not the printed 0.4659 | the defence is that the integral solver agrees, but they share plumbing |
| C-method steep-wall refusal, deep-profile digits | self-consistent, tautology risk | pinned integral-solver outputs (`STEEP_80_REFERENCE`, 0.844609) | |
| Integral flat mirror, finite-conductivity flat | independent | `test_integral_core.py::test_unit_specular`; `test_integral_finite.py::test_te_amplitude_is_r_s` (1e-4) | |
| Integral perfect conductor | independent, loose | `test_integral_core.py::TestMaystrePopovTable41`, 0.4659 at abs 1e-3 | passes TE only because 1e-3 exceeds the 5.4e-4 gap the C-method test records |
| Integral G&S Table 3 | independent, **looser than documented** ✓ | `test_integral_finite.py::test_reflected_and_transmitted_orders` at abs 2e-3 | `docs/theory/integral.md:315` says "within 6e-5" |
| Integral against PCGrate | independent, circular, corpus-gated | `test_integral_corpus.py::test_orderwise_agreement_on_the_identified_profile` (<2e-3) | profile, α and γ were recovered from the same table; `corpus.toml` has `confirmed = false` |
| Integral Green's function | independent | `test_integral_greens.py::test_accelerated_matches_brute_force_off_the_boundary` (1e-12) | |
| Integral graded mesh, undercut, kernels | self-consistent | convergence ladders, energy balance, GEMM against pointwise | undercut is energy balance only |
| Integral measured-polygon smoothing | independent (cross-method) | `test_integral_measured_polygon.py::test_reflectance_closes_on_the_c_method` | |
| Cross-method agreement | independent, shared plumbing | `TestAgreementWithTheIntegralSolver` in both RCWA and C-method tests; `test_four_way.py` | strong for the physics in each solver, blind to §2's shared parts |
| TASTE measured comparison | **unchecked by tests** ✓ | `test_measured.py:379`, "Agreement is not asserted" | the RMS and the three-solver 2e-4 agreement that "not experimental since 0.3" rests on are in `findings.md` and `roadmap.md` only |

## 4. Scalar solver and geometry

| feature | tier | evidence | notes |
|---|---|---|---|
| Absolute normalisation, obliquity factor | **independent** | `test_perturbation.py::TestAgainstFirstOrderPerturbationTheory`, Rayleigh first order, rtol 2 to 5e-3 | one depth (1e-5), TE-like only. The strongest model check here |
| Fourier quadrature: sawtooth, lamellar, sinusoid | independent (numerics only) | `test_scalar.py::TestSawtoothAgainstClosedForm` (1e-8), `::TestBinaryPhaseGratingAgainstClosedForm`, `::TestSinusoidAgainstBesselForm` | the closed forms share the scalar model, so they check the quadrature, not the physics |
| Agreement with the integral solver | independent (cross-method) | `test_cross_method.py::TestAgreementInTheOverlapRegime` (2e-3), `::TestAgreementOnAnAsymmetricGroove` | divergence tests assert a trend, not a size |
| Littrow angles | independent | `test_littrow.py::TestAgainstPyXFocus` | PyXFocus is the same author's |
| Facet graze, two-facet depth | independent, weak | single thesis values (Chapter-2.tex:1046, :1022) | |
| Reciprocity, energy behaviour | self-consistent | `test_checks.py::test_scalar_solver_is_reciprocal`, with a negative control | |
| **Coated scalar: groove-resolved reflectivity** (`"local"`, `"average"`) | self-consistent ✓ | `test_scalar.py::TestAbsoluteEfficiency` re-implements the model | the only rigorous comparison is a factor-of-two bound (`test_four_way.py:152`). This is what makes coated scalar numbers absolute |
| Measured or sampled profiles | self-consistent / unchecked | `::test_a_sampled_sawtooth_reproduces_the_analytic_sawtooth`; `::test_runs_on_a_measured_boundary` checks finite output only | |
| Validity guards (reduced ratio, Fraunhofer, critical angle) | unchecked as physics | `::TestProvenance` checks that they fire | thresholds are choices |
| Resolving power | independent / self-consistent | `test_resolution.py::test_r_is_order_times_n` (textbook); profile and Rayleigh checks are internal | |
| Convergence plateau rule | self-consistent, tautology risk | `test_convergence.py`, a pinned measured sequence | |

## 5. Materials, mirror and file formats

| feature | tier | evidence | notes |
|---|---|---|---|
| refractiveindex.info formulas and tables | independent | `test_refractiveindex.py::test_n_matches_the_database_tool` against `reference.json`; Malitson and N-BK7 d-line values | the script that made `reference.json` is not committed |
| CXRO reader | self-consistent | `test_materials.py::test_gold_at_1_nm` against `make_au.py`'s formula | **stale label** ✓: the module docstring and `class TestAgainstPublishedValues` still say "published Henke constants"; `Au.txt` is synthetic since 2026-11-11. `::test_and_is_a_real_table` asserts "CXRO" on it |
| Real CXRO file layout, `factors:` with real asf files | unchecked | hand-written `CXRO_TEXT`, invented `.nff` files | |
| `typical:` model | self-consistent, strongly | `test_typical.py::TestTheClosedFormIsTheKramersKronigTransform`, sum rule | its accuracy against CXRO (`docs/typical-materials.md`) was measured once on uncommitted tables and cannot be re-run |
| `factors:` number density | independent, weak | the author's dissertation, Appendix C, 1% | |
| Fresnel amplitudes | independent (limits) | normal incidence, total reflection, R ≤ 1 | no tabulated R(θ) |
| Roughness (Névot–Croce, Debye–Waller) | self-consistent | checked against each other | |
| Mirror stack, s polarization | independent | `test_mirror.py::test_a_stack_agrees_with_characteristic_matrices` (1e-12) | |
| Mirror stack, p polarization; buried roughness | unchecked | p only on a bare substrate; `mirror.md` admits the second | |
| `.ari`, `.nk`, boundary JSON/CSV, `.txl` | self-consistent | round trips; `.txl` checked by a parser written in the test | |
| `.ggp` writer | independent | `test_ggp.py::test_writes_the_exact_pcgrate_line_format`, bytes of a file PCGrate 6.7.1 loads | the source file is not committed |
| `.ggp` reader, efficiency-table import | independent, corpus-gated | `test_ggp.py::TestPcgrateOwnFiles`, `test_efficiency_table.py::TestCorpus` | skip without the corpus |
| ALS campaigns | independent, gated | `test_measured.py::test_round_3a_recovers_the_published_table_1` | needs `als_eff_meas` and the beamline tree |

## 6. Metrology and the GUI

| feature | tier | evidence | notes |
|---|---|---|---|
| Nanoscope `.spm` reader | independent, archive-gated | `test_spm_reader.py::test_against_gwyddions_own_numbers`, 70 pairs against Gwyddion | the strongest metrology check, and it skips everywhere it has been run lately |
| Asylum `.ibw`, Zeiss SEM readers | self-consistent | files written by the tests' own builders | no real file. The roadmap's "379 micrographs" has no record in the repo |
| Blaze angle and period | independent (idealised) | analytic sawtooth fixture; `test_orientation.py` with closed-form cos θ errors | no tip, drift or creep. Real-scan values are regression (`docs/metrology/BASELINES.md`) |
| Tip erosion (Villarrubia) | self-consistent, **tautology risk** | `synthetic_dilated_scan.txt` comes from `core.tip.dilate`, and erosion is tested against it | the roadmap's "validated pointwise" is an inverse checked against its own forward model |
| Tip presets | unchecked | constants restated from one reading of a Bruker datasheet | |
| **ICC estimator** | self-consistent, tautology risk ✓ | `test_icc.py::test_icc_matches_hand_computed_value` restates the formula | `stats/icc.py:86-90` uses var(group means)/(that + within). Without the −σ²_w/k correction it reads about 1/(k+1) on independent data, where ICC(1) reads about 0. Conservative (inflates the SEM) but not the textbook estimator |
| Design effect, corrected SEM | self-consistent | Kish's formula restated | real-data cases skip |
| SEM cross-section to trapezoid | independent geometry, same-author imaging model | `test_cross_section.py` on `xsec_synth.render` | no real micrograph run, and the roadmap says so |
| Copy as Python, examples, `docs/python.md` | self-consistent | `test_gui_script.py::test_it_reproduces_the_solve` (1e-12); `test_examples.py` string regressions | shows the GUI reports the solver's numbers. `examples/metrology/uncertainty_demo.py` is not run |

## 7. Where the author's reading is most needed

Ranked by what an error would cost and how little stands behind it now.

1. **The conical grating equation and the obliquity factor.** Every solver and
   every oracle shares `sin_beta`, and the obliquity factor rests on one
   perturbation file at one depth. A hand check against ISSI §2.1 or Werner,
   and of `4 k_z0 k_zm |ĝ|²`, covers the base the whole package stands on.
2. **The groove-resolved reflectivity model** (`scalar.md` §9,
   `_scalar_groove.py`). It is what makes coated scalar numbers absolute, and
   nothing independent supports it.
3. **The TASTE comparison is not in any test.** It is the evidence for "not
   experimental since 0.3". Commit it as a gated test, or record the run that
   produced it.
4. **The published values not matched.** Maystre & Popov TE (0.4659 against
   0.466441), Granet TE, Gralak s, and the 1.1e-3 G&S metallic conical gap.
   Read the tables. Granet TE matching the Littrow order to 9e-6 suggests a
   convention, not a misprint.
5. **The ICC estimator.** Decide whether the variance-ratio form is intended.
   Every corrected SEM depends on it.
6. **Run the gated tests where the data lives** (PCGrate corpus, AFM archive,
   ALS tree) and record the date and result. They are the only external
   evidence for several areas, and none has a recorded recent pass.
7. **Tip erosion against a published case**, one Villarrubia example or one
   Gwyddion certainty map.

## 8. Doc and test mismatches found

These are small and fixable, and they are the kind a reviewer would find. The
first three are fixed in pull request #85, not yet merged when this was
written. The rows in §3 and §5 that cite them describe `master` before it.

- `test_materials.py`: module docstring and `TestAgainstPublishedValues` claim
  published Henke constants; the table is synthetic. ✓
- `docs/theory/integral.md:315` claims G&S Table 3 within 6e-5; the test holds
  2e-3. ✓
- `test_cmethod.py:367` pins the code's TE value in a test named for a
  published table. ✓
- The roadmap's "validated pointwise" (tip erosion) and "379 micrographs"
  (SEM) describe evidence that is either self-referential or not in the repo.

## 9. Not covered by this draft

- No test was run. "Passes" here means "is asserted".
- Source modules were not read line by line against their theory documents.
- The cited papers were not opened, so transcribed values are taken as quoted.
- Not read in full: `test_rayleigh_anomaly.py`, the graded-mesh and operator
  internals, the analyzer, swath and boundary-pipeline tests, the Qt GUI
  tests, and the provenance panel's values.
- Whether the corpus-gated tests pass today is unknown.
