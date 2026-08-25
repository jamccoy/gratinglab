# Groove metrology: from an AFM scan to a boundary profile

Implemented in [`metrology/`](../../src/gratinglab/metrology/). Units and sign
conventions follow [`conventions.md`](../conventions.md), §1 in particular — the
metrology package spans four unit regimes and the table there is normative.

This page exists for one reason: **an efficiency computed from a measured
profile inherits every assumption made while producing that profile**, and those
assumptions are invisible in the resulting curve. Anyone reading such a curve
should be able to find them written down.

---

## 1. What the pipeline assumes

| Assumption | Consequence when it fails |
|---|---|
| The scan is of a periodic grating, and the period estimate is close enough for peak detection to find groove minima | Grooves are missed or spurious; the measured period is wrong, and with it the groove depth (§3) |
| Every groove in the scan is a sample of *the same* groove | Real variation across the scan is averaged into a shape that no individual groove has |
| Row-averaging the 2-D image to 1-D loses nothing that matters | True only if the grooves run parallel to the image rows. A rotated scan smears the facet |
| The averaged groove tiles seamlessly into a periodic boundary | Enforced rather than assumed — see §2 — but enforcement is itself an operation on the data |
| The AFM tip is small compared to the features | The steep facet is the tip shape, not the grating, and the apex and trough are rounded off — which removes depth while leaving the mid-facet slope a line fit lands on intact |
| The groove is faceted, so a fitted facet angle and the groove depth describe the same shape | False on a rounded groove, and then the two disagree. The solver uses the **depth** |

The last two are not mitigated anywhere in this package and have no diagnostic.
On the one measured grating in this repo they are both false: see `findings.md`,
"The measured groove is rounded, not faceted", where the fitted facet angle
(27.91° ± 2.13°) and the angle implied by the depth (20.33°) differ by 3.6σ.

**The cheapest available diagnostic does not exist yet**, and would be: compute
the blaze angle the measured depth implies, and compare it with the fitted one.
They agree on a sharp groove and diverge on a rounded one, and the divergence is
precisely the error that reaches the efficiency. Two numbers the pipeline
already has.

**Threshold it in degrees, not in sigma.** The fitted angle's scatter is small
on clean data, so a gap of a few tenths of a degree is already several sigma on
a groove that is perfectly sharp — the synthetic fixture shows 0.33°, which is
3.8σ and means nothing. The real scan's gap is 7.6°. Absolute size discriminates;
sigma does not.

Note what the diagnostic would *not* say. A rounded groove and a blunt tip
produce the same profile, so a divergence flags "this is not a sharp sawtooth"
rather than "your tip is bad". A **flat land** — a genuinely unfaceted flat
within the period, common on real blazed gratings — produces a depth deficit
too, and is the first alternative to exclude: it shows as a spike at zero in the
local-slope distribution, which rounding does not.

---

## 2. What `normalize_profile` does, in order

Four operations, in
[`boundary/average.py`](../../src/gratinglab/metrology/boundary/average.py),
kept in the original script's order because changing the order changes the
answer:

1. **Endpoint flattening.** Subtract the straight line through the first and
   last sample, so the trace starts and ends at zero. This is *not* the
   background flattening used for blaze-angle work, which removes a fitted
   background so a facet can be measured. Swapping one for the other leaves a
   step at the period boundary — the defect that makes a PCGrate efficiency
   curve wrong while the file still looks fine.
2. **Shift and normalise.** Reference to the mean of the two half-minima, then
   divide both axes by the period.
3. **Roll the trough to the boundary.** `x = 0` and `x = 1` are the deepest
   point; the apex sits somewhere in between. Only `y` is rolled — `x_norm`
   stays a uniform ramp.
4. **Optional smoothing**, `uniform_filter1d(size=5, mode='wrap')` by default,
   then the endpoints are forced back to exactly zero.

**Handedness is preserved and never mirrored.** The facet orientation matches
increasing scan-x in the source data. Whether that matches the sense the solver
expects is a question about the *mount*, not about this package — see
`conventions.md` §3, "Which way the profile parameter runs".

### The stretch

Step 2 divides by the *averaged window width*, and step 3's output spans exactly
`[0, 1]` regardless of what that width was. If a groove near the scan edge
narrows the window every groove shares, the profile is then **stretched to fill
a period it never occupied**. This is real and was found empirically:
`test_edge_rule_reduces_x_stretch` exists because of it, and
`ggp_min_half_width` exists to mitigate it by discarding clipped grooves.

A stretched profile is not obviously wrong in a plot. It is wrong in the
efficiency, because the facet angle it implies is wrong.

---

## 3. Why the period is load-bearing

The exported profile is dimensionless. Physical depth is `y_norm × period`, and
that depth enters the scalar phase term directly:

    phase = (2π / λ) · height · sin γ

so a fractional error in the period is a fractional error in the phase, at every
wavelength. The period is not metadata.

It is measured, not assumed: the mean spacing of detected groove centres,
`period_nm` on `BoundaryProfile`. The initial `period_est` from settings is only
used to seed peak detection, and is replaced as soon as there are two grooves to
measure between.

**A `.ggp` file cannot carry this number.**
[`BoundaryProfile.to_problem`](../../src/gratinglab/metrology/boundary/pipeline.py)
exists so that a profile reaching a solver in-process keeps it. The file-based
route needs the period supplied again from somewhere else, which is what
`benchmarks/corpus.toml` is.

---

## 4. What is measured, and what is not

| Quantity | Measured? | Where |
|---|---|---|
| Period | **yes**, from groove spacing | `BoundaryProfile.period_nm` |
| Groove depth, peak-to-valley | **yes**, as a fraction of the period | `metrics` |
| Blaze / anti-blaze facet angle | **yes**, by linear fit with a propagated uncertainty | `core.analysis.extract_blaze_angle` — the *other* branch of the pipeline |
| Groove-to-groove form spread | **yes** | `y_std_nm`, the ±1σ band |
| Max sidewall angle | derived, not fitted | `metrics['max_angle_deg']` — the steepest local slope of the normalised profile |
| Apex angle | **no** | not computed anywhere |
| **RMS surface roughness** | **no** | — |

### Roughness is absent, and the near-misses are not substitutes

`Problem.roughness` means high-spatial-frequency surface roughness in nm — the
quantity Névot–Croce damps a Fresnel coefficient with. Nothing here computes it.

Two quantities look like they might and do not:

- `y_std_nm` is the spread *between grooves*: form variation, at the spatial
  scale of a whole groove. Névot–Croce is about structure far below that.
- `metrics['rms_slope']` is a slope statistic of the normalised profile, in
  dimensionless units. It is not a length at all.

Passing either into `Problem.roughness` would produce a confident wrong
efficiency instead of an obviously missing one, so `to_problem` does not, and a
test pins the refusal. The honest quantity — the RMS residual of each groove
about the average — is computable from data the pipeline already holds and is
not yet computed. Until it is, an absolute efficiency from a measured profile
has **no roughness damping at all**, which overestimates it.

---

## 4b. Tip correction: what erosion can and cannot give back

An AFM image is the surface **dilated** by the tip: at every pixel the tip
descends until first contact, and the recorded height is where the apex
stopped — above the true surface wherever contact happened on the flank.
Villarrubia (1997) gives the exact morphological inverse this pipeline
implements in `core/tip.py`, opt-in through
`AnalysisSettings.tip_correction = 'erosion'` with a parametric tip: a cone of
`tip_half_angle_deg` (from the axis) capped by a spherical apex of
`tip_radius_nm`. It runs on the 2-D array right after image flattening, in
both branches, so the blaze fit and the boundary profile see the same surface.

Three statements, each load-bearing:

- **Erosion is a least upper bound, not a resurrection.** Everywhere the apex
  made contact, the bound *is* the surface and the recovery is exact — machine
  precision on clean data, the noise floor on the fixture
  (`tests/metrology/test_tip.py`). Everywhere else the tip physically could
  not reach, and no algorithm recovers what is there. The trough wedge between
  facets is the canonical case: a sphere of radius R stands off a corner by a
  distance that scales with R, and that depth is *gone from the data*.
- **The certainty map is the deliverable, not a by-product.** A pixel is
  certain when it is the unique contact point for some image position; the
  fraction lands in `metrics['tip_certain_fraction']` and the metrics sidecar,
  and the uncertain remainder must be read as "at most this high", not
  "this high". A corrected profile without that number would manufacture
  exactly the false confidence the correction exists to remove.
- **A facet steeper than the flank is unrecoverable outright.** The flank
  rises at `90 − half_angle` degrees; an 18° tip clears a 70.5° anti-blaze
  facet by 1.5°. That near-parallelism, not the nanometre apex, is what limits
  blazed-grating metrology at these geometries.

The correction is off by default because it changes measured numbers — depth
especially — and a corrected and an uncorrected depth are different
measurements. Where it is on, the settings and the certain fraction are
recorded with the output (metrics, sidecar, and the CLI's log line).

Validation is the dilated fixture: the same synthetic sawtooth as the sharp
fixture with the identical noise field, imaged through a known R = 20 nm tip
by `dilate` — committed truth next to committed image, so recovery is asserted
pointwise rather than hoped. What the milestone found on the real TASTE scan
is in `findings.md` ("A nominal tip does not explain the rounded groove"):
erosion with the nominal probe recovers nothing, and the rounding needs an
~80 nm apex to be a tip artefact at all.

---

## 5. Two branches, one front-end

Loading, image flattening, row-averaging and groove detection are shared. After
that the pipeline forks:

- **Blaze angle** — `core.analysis.extract_blaze_angle` fits a line to each
  facet independently, per groove, with row-group statistics and an ICC
  correction for the correlation between rows of the same scan.
- **Boundary profile** — `boundary.build_boundary_profile` averages the grooves
  into one and normalises it.

They were once separate scripts with their own copies of the front-end, and
drifted until a bug fixed in one persisted in the other. That is why the shared
front-end is shared.

**The two branches do not talk to each other.** `build_boundary_profile` does
not run the facet fit, so `to_problem(blaze_angle=...)` takes the angle as an
argument rather than finding it. Connecting them is real work and deliberately
not done: the fit has settings of its own (facet trimming in particular) and
quietly adopting defaults for them would bury a choice that changes the answer.
