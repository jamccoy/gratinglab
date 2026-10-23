# Groove metrology: from an AFM scan to a boundary profile

Implemented in [`metrology/`](../../src/gratinglab/metrology/). Units and sign
conventions follow [`conventions.md`](../conventions.md), §1 in particular. The
metrology package spans four unit regimes, and the table in that section is
normative.

An efficiency computed from a measured profile inherits every assumption made
while producing that profile, and those assumptions are not visible in the
resulting curve. This page records them so that anyone reading such a curve can
find them.

---

## 1. What the pipeline assumes

| Assumption | Consequence when it fails |
|---|---|
| The scan is of a periodic grating, and the period estimate is close enough for peak detection to find groove minima | Grooves are missed or spurious; the measured period is wrong, and with it the groove depth (§3) |
| Every groove in the scan is a sample of *the same* groove | Real variation across the scan is averaged into a shape that no individual groove has |
| Row-averaging the 2-D image to 1-D loses nothing that matters | True only if the grooves run parallel to the image rows. A rotated scan smears the facet |
| The averaged groove tiles seamlessly into a periodic boundary | Enforced rather than assumed (see §2), but enforcement is itself an operation on the data |
| The AFM tip is small compared to the features | The steep facet is the tip shape, not the grating, and the apex and trough are rounded off; this removes depth while leaving intact the mid-facet slope on which a line fit lands |
| The groove is faceted, so a fitted facet angle and the groove depth describe the same shape | False on a rounded groove, and then the two disagree. The solver uses the **depth** |

The last two assumptions are not mitigated anywhere in this package and have no
diagnostic. On the one measured grating in this repository both are false: see
`findings.md`, "The measured groove is rounded, not faceted", where the fitted
facet angle (27.91° ± 2.13°) and the angle implied by the depth (20.33°) differ
by 3.6σ.

The least expensive diagnostic available has not yet been implemented. It would
compute the blaze angle implied by the measured depth and compare it with the
fitted angle. The two agree on a sharp groove and diverge on a rounded one, and
the divergence is the same error that reaches the efficiency. Both numbers are
already available within the pipeline.

Such a diagnostic should be thresholded in degrees. The scatter of the fitted
angle is small on clean data, so a gap of a few tenths of a degree already
amounts to several sigma on a groove that is perfectly sharp: the synthetic
fixture shows 0.33°, which is 3.8σ and carries no meaning. The gap on the real
scan is 7.6°. The absolute size of the gap discriminates between the two cases,
and its size in sigma does not.

The diagnostic would also have limits. A rounded groove and a blunt tip produce
the same profile, so a divergence indicates that the groove is not a sharp
sawtooth without indicating that the tip is at fault. A *flat land*, a genuinely
unfaceted flat within the period that is common on real blazed gratings, also
produces a depth deficit and is the first alternative to exclude. It appears as
a spike at zero in the local-slope distribution, which rounding does not
produce.

---

## 2. What `normalize_profile` does, in order

The function performs four operations, in
[`boundary/average.py`](../../src/gratinglab/metrology/boundary/average.py).
They are kept in the order used by the original script, because changing the
order changes the answer.

1. **Endpoint flattening.** The straight line through the first and last sample
   is subtracted, so that the trace starts and ends at zero. This differs from
   the background flattening used for blaze-angle work, which removes a fitted
   background so that a facet can be measured. Substituting one for the other
   leaves a step at the period boundary, a defect that makes a PCGrate
   efficiency curve wrong while the file still appears correct.
2. **Shift and normalise.** The profile is referenced to the mean of the two
   half-minima, and both axes are then divided by the period.
3. **Roll the trough to the boundary.** The deepest point is placed at `x = 0`
   and `x = 1`, and the apex lies somewhere between them. Only `y` is rolled;
   `x_norm` remains a uniform ramp.
4. **Optional smoothing**, by default `uniform_filter1d(size=5, mode='wrap')`,
   after which the endpoints are forced back to exactly zero.

Handedness is preserved and never mirrored. The facet orientation matches
increasing scan-x in the source data. Whether that matches the sense the solver
expects is a question about the *mount* and lies outside this package; see
`conventions.md` §3, "Which way the profile parameter runs".

### The stretch

Step 2 divides by the *averaged window width*, and the output of step 3 spans
exactly `[0, 1]` whatever that width was. If a groove near the scan edge narrows
the window that every groove shares, the profile is consequently stretched to
fill a period it never occupied. This effect is real and was found empirically:
`test_edge_rule_reduces_x_stretch` was written because of it, and
`ggp_min_half_width` mitigates it by discarding clipped grooves.

A stretched profile does not look obviously wrong in a plot. Its error appears in
the efficiency, because the facet angle it implies is wrong.

---

## 3. Why the period is load-bearing

The exported profile is dimensionless. Physical depth is `y_norm × period`, and
that depth enters the scalar phase term directly:

    phase = (2π / λ) · height · sin γ

so a fractional error in the period is a fractional error in the phase, at every
wavelength. The period therefore enters the result directly and cannot be
treated as metadata.

The period is measured as the mean spacing of detected groove centres, stored as
`period_nm` on `BoundaryProfile`. The initial `period_est` from settings is used
only to seed peak detection, and it is replaced as soon as there are two grooves
between which to measure.

A `.ggp` file cannot carry this number.
[`BoundaryProfile.to_problem`](../../src/gratinglab/metrology/boundary/pipeline.py)
exists so that a profile reaching a solver in-process retains it. The file-based
route requires the period to be supplied again from another source, and
`benchmarks/corpus.toml` serves that purpose.

---

## 4. What is measured, and what is not

| Quantity | Measured? | Where |
|---|---|---|
| Period | **yes**, from groove spacing | `BoundaryProfile.period_nm` |
| Groove depth, peak-to-valley | **yes**, as a fraction of the period | `metrics` |
| Blaze / anti-blaze facet angle | **yes**, by linear fit with a propagated uncertainty | `core.analysis.extract_blaze_angle`, the *other* branch of the pipeline |
| Groove-to-groove form spread | **yes** | `y_std_nm`, the ±1σ band |
| Max sidewall angle | derived, not fitted | `metrics['max_angle_deg']`, the steepest local slope of the normalised profile |
| Apex angle | **no** | not computed anywhere |
| **RMS surface roughness** | **no** | (none) |

### Roughness is not computed, and related quantities are not substitutes

`Problem.roughness` denotes high-spatial-frequency surface roughness in nm, the
quantity with which Névot–Croce damps a Fresnel coefficient. Nothing in this
package computes it.

Two quantities appear related but do not measure it:

- `y_std_nm` is the spread *between grooves*: form variation at the spatial
  scale of a whole groove. Névot–Croce concerns structure far below that scale.
- `metrics['rms_slope']` is a slope statistic of the normalised profile, in
  dimensionless units. It is not a length.

Passing either into `Problem.roughness` would produce a confidently wrong
efficiency where an obviously missing one is expected, so `to_problem` refuses
to do so, and a test pins the refusal. The appropriate quantity, the RMS residual
of each groove about the average, is computable from data the pipeline already
holds but is not yet computed. Until it is, an absolute efficiency from a
measured profile has no roughness damping at all, and is therefore an
overestimate.

---

## 4b. Tip correction: what erosion can and cannot give back

An AFM image is the surface *dilated* by the tip. At every pixel the tip
descends until first contact, and the recorded height is the position at which
the apex stopped, which lies above the true surface wherever contact occurred on
the flank. Villarrubia (1997) gives the exact morphological inverse, which this
pipeline implements in `core/tip.py`. It is enabled through
`AnalysisSettings.tip_correction = 'erosion'` with a parametric tip: a spherical
apex of `tip_radius_nm` on two straight flanks, `tip_half_angle_deg` on the -x
side and `tip_half_angle_right_deg` on the +x side (both measured from the axis;
the second defaults to the first, giving a symmetric cone). The two flanks are
the cross-section of a pyramid in the scan plane, so they depend on how the
cantilever lay against the fast axis. A ScanAsyst-Air, for example, shows its
17.5° sides across the cantilever and its 15° front and 25° back along it. The
correction runs on the 2-D array immediately after image flattening, in both
branches, so that the blaze fit and the boundary profile see the same surface.

With unequal flanks the erosion uses the reflected tip,
`rec(y) = min_u [img(y - u) + tip(u)]`. A symmetric tip conceals the
reflection, which is why it is stated here.

Three properties of the correction determine how its output may be used:

- **Erosion gives a least upper bound.** Wherever the apex made contact, the
  bound coincides with the surface and the recovery is exact: to machine
  precision on clean data, and to the noise floor on the fixture
  (`tests/metrology/test_tip.py`). Elsewhere the tip physically could not reach
  the surface, and no algorithm can recover what lies there. The trough wedge
  between facets is the canonical case: a sphere of radius R stands off a corner
  by a distance that scales with R, and that depth is absent from the data.
- **The certainty map is a primary output.** A pixel is certain when it is the
  unique contact point for some image position. The certain fraction is stored
  in `metrics['tip_certain_fraction']` and in the metrics sidecar, and the
  uncertain remainder must be read as an upper bound on the height, not as the
  height itself. A corrected profile reported without that number would create
  the false confidence that the correction is intended to remove.
- **A facet steeper than the flank facing it cannot be recovered at all.** A
  flank rises at `90 − half_angle` degrees, and a wall rising towards +x is
  followed by the +x flank. An 18° tip clears a 70.5° anti-blaze facet by 1.5°.
  The 25° back face of a ScanAsyst-Air, a 65° flank, does not clear it at all,
  so along the cantilever the orientation decides whether that facet exists in
  the data. At these geometries it is this near-parallelism, and not the
  nanometre apex, that limits blazed-grating metrology.

The correction is off by default because it changes measured numbers, depth in
particular, and a corrected depth and an uncorrected depth are different
measurements. When it is enabled, the settings and the certain fraction are
recorded with the output (in the metrics, the sidecar, and the log line of the
CLI).

The correction is validated against the dilated fixture: the same synthetic
sawtooth as the sharp fixture, with the identical noise field, imaged through a
known R = 20 nm tip by `dilate`. The truth is committed next to the image, so
recovery is asserted pointwise. The result of applying it to the real TASTE scan
is recorded in `findings.md` ("A nominal tip does not explain the rounded
groove"): erosion with the nominal probe recovers nothing, and the rounding
would require an apex of ~80 nm to be a tip artefact at all.

---

## 5. Two branches, one front-end

Loading, image flattening, row-averaging and groove detection are shared. After
these steps the pipeline divides into two branches:

- **Blaze angle.** `core.analysis.extract_blaze_angle` fits a line to each facet
  independently, per groove, with row-group statistics and an ICC correction for
  the correlation between rows of the same scan.
- **Boundary profile.** `boundary.build_boundary_profile` averages the grooves
  into one and normalises it.

The two branches were once separate scripts with their own copies of the
front-end, and these copies drifted apart until a bug fixed in one persisted in
the other. The front-end is shared for that reason.

The two branches do not exchange results. `build_boundary_profile` does not run
the facet fit, so `to_problem(blaze_angle=...)` takes the angle as an argument
instead of obtaining it. Connecting them would require substantial work and has
deliberately not been done: the fit has settings of its own (facet trimming in
particular), and adopting their defaults silently would conceal a choice that
changes the answer.
