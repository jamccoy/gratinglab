# Changelog

Notable changes, newest first. The project is pre-1.0, so a minor version may
change an interface; any such change is stated here when it is made. The
versioning policy, including what 0.2 and 1.0 will mean, is in
[`docs/roadmap.md`](docs/roadmap.md#versioning).

The format loosely follows [Keep a Changelog](https://keepachangelog.com). The
version is read from `gratinglab.__version__`, and a release commit moves the
entries under "Unreleased" into a section named for the new version;
`tests/test_changelog.py` fails if it does not.

## 0.2.4 (2026-10-28)

### Fixed

- The README, which is also the PyPI page, still called RCWA and the
  C-method "planned". Both shipped in 0.2.0, with GUI tabs in 0.2.1; its
  Status and Scope sections and the method table now say so, with what
  each has and has not been checked against. No code changes.

## 0.2.3 (2026-10-28)

### Changed

- The AFM window's controls column is a splitter, not a fixed 272 px strip:
  it starts wide enough for its controls, drags wider or narrower down to a
  floor at which labels still fit, and gives the plot the rest. The same
  layout serves the Import, Analysis and Boundary tabs
  (`metrology.gui.qt.layout.controls_beside`).
- One title style for every AFM plot (`gui.theme.set_titles`): what the
  panel is, left-aligned and semibold, and the figures it reports (scan
  width, groove count, method, depth) on the right in the theme's muted
  colour. Titles no longer run to a second line or carry explanations.

## 0.2.2 (2026-10-28)

### Changed

- Requires endstation 0.2.6, and uses its copies of components this package
  had carried its own versions of. Nothing a user sees changes beyond the
  wording of one message:
  - The missing-PySide6 check in `gratinglab-gui`, `gratinglab-app` and
    `gratinglab.metrology.gui.require_qt` is endstation's `require_extra`,
    and its message is endstation's (`pip install -e '.[gui]'`, quoted for
    zsh). `QT_API` is set by importing `endstation.qt`. The console scripts
    still exit with the message rather than a traceback.
  - `metrology.config` finds the project root with
    `endstation.paths.find_project_root` (still `pyproject.toml` only) and
    reads `GRATINGLAB_AFM_DIR` with `data_root`, which also expands `~` in
    it. Both stay strings. The suite's corpus and scan helpers skip through
    `endstation.testing`.
  - The metrology wiki's page list, order and titles come from
    `endstation.docs.index.load_pages`. `WikiPage` and `is_empty` are
    unchanged.
  - `SolveCancelled` is now a subclass of `endstation.progress.Cancelled`.
    It keeps its name and module, so `except SolveCancelled` and
    `raise SolveCancelled` work as before; `except Cancelled` now catches it
    too. The efficiency window's progress callback is endstation's
    `cancelling(throttled(...))`, still throttled to whole percent with the
    cancel checked on every call, and it raises `Cancelled`.
  - Both windows tell a stale result from a current one with
    `endstation.qt.tokens.StaleGuard`.
  - The guard that keeps Python's cycle collector on the GUI thread is
    `endstation.qt.gc_guard`, which this package's copy moved into;
    `gratinglab.gui.qt.gc_guard` is removed.
  - `gratinglab.gui.state.FormErrors` is `endstation.forms.FormErrors`. This
    is the one interface change here: it is no longer a `ValueError` (nothing
    in this package caught it as one), and its `str()` puts each field error
    on its own line instead of joining them with `"; "`. `.errors` is
    unchanged.

## 0.2.1 (2026-10-28)

### Added

- A "C-method (Chandezon)" tab in the GUI: `truncation` (1 to 256, default
  32), the same `conductivity` selector as the integral tab, progress and a
  working Cancel, and the convergence sweep. The solver's warnings (corners,
  condition number, energy balance) print in the provenance panel. A too-small
  truncation is refused on its field before a solve starts, by
  `cmethod.truncation_floor`, the function the solver itself now applies.
  Solvers may declare `Capabilities.accuracy_ceiling`, and the default
  convergence ladder stops there: the C-method's ran on to 2048, so a sweep
  that had not settled by 256 failed on a refused rung instead of reporting.
- An RCWA tab in the efficiency window, the third solver tab, on the same
  skeleton as the integral one. It has the truncation (the knob the
  Check convergence button sweeps) and the staircase's slices (blank follows
  the truncation), and a note that says, once the geometry is known, whether
  the profile is one exact lamellar slab or a staircase, with the TM caveat
  for the latter. Its form refuses what the solver would: no coating, a scan
  outside the material's table, gain, and a truncation that cannot hold every
  reflected or transmitted order, the last through the solver's own
  `required_truncation`, now public, so the two cannot disagree.
- Coated gratings in the window: a Layers field in the geometry panel's
  Surface group, written `Au 20, Cr 5` (material and thickness in nm,
  innermost first), fills `Problem.layers`. Every bad entry is named in one
  error on the field, and a material must be installed, as the coating must.
  A tab whose solver cannot model layers (Scalar, Integral, C-method)
  refuses on the field before anything reaches the worker, rather than
  solving the bare grating; the RCWA tab solves them.

## 0.2.0 (2026-10-27)

A minor release, for two reasons the versioning policy names: the in-house
RCWA backend arrives, with the C-method beside it, and `absorption` changes
meaning (see Changed). The metrology milestone that 0.2 was planned around
has not been signed off yet; see `docs/roadmap.md`.

### Added

- A C-method (Chandezon) solver, registered as `cmethod`, with the integral
  solver's two boundary conditions: `conductivity="perfect"` in any mount and
  `"tabulated"`, one interface into a material read from the coating, coupled
  `E_z`/`H_z` off-plane, absolute efficiencies. Into a lossless material the
  transmitted orders are recorded as the scan's `transmission`, with
  `absorption` keeping the remainder. The knob is `truncation`. It reproduces
  Goray & Schmidt's Table 3, which is Li's C-method, to 1.4e-5 on all twelve
  orders, Fresnel exactly on a flat boundary, and the integral solver to that
  solver's own error on sinusoids in-plane, at 60° and off-plane at 1.5° on
  CXRO gold. Vertical facets, sidewalls and undercuts are refused. Two weak
  regions are measured and reported in provenance rather than hidden: on a
  sawtooth it converges at about first order and its energy balance closes
  long before its efficiencies do, and on a profile deeper than about half
  the period a high truncation exhausts double precision and makes the answer
  worse. The method note is `docs/theory/cmethod.md`; no GUI tab yet.
- Coated gratings: `Problem.layers`, a tuple of `ConformalLayer(material,
  thickness)`, innermost first, each following the surface beneath it as its
  true offset, so convex corners are rounded and concave ones sharp
  (conventions §3). The RCWA models them; every other solver refuses them,
  through a new `Capabilities.handles_layers`, rather than solving the bare
  grating. Checked: a vacuum layer reproduces the bare grating to 1e-12; a
  layer of the grating's own material matches the integral solver on the
  hand-written offset groove to 2e-5 (dielectric) and 2e-4 (metal, TE); a film
  on a vanishing groove is the Airy film. Metal in TM converges, slowly, as
  staircases do on rounded metal corners.

- Transmitted orders. `EfficiencyScan.transmission` is a new optional
  `Transmission` record (its own orders, efficiency and propagating flags,
  under the same zero-means-evanescent rules), filled by the RCWA and by the
  integral solver's tabulated mode whenever the material is lossless. A
  substrate admits more orders than the vacuum, so it has its own order axis.
  `check_energy_balance` adds its total, and the window's plot and provenance
  show `T` and `Σ + T + A`.
- A published conical benchmark for RCWA: Goray & Schmidt's lamellar Tables
  12.1 and 12.2 (ch. 12 of Popov (ed.), *Gratings*, 2014), with the chapter's
  own definitions of the conical angles and of the polarization angles `delta`
  and `psi` mapped onto `solve_stack`, which takes the elliptically polarized
  incident field directly. The dielectric table is matched in all eight
  orders, reflected and transmitted, to half a unit of its last printed digit
  at M = 80, and the opposite handedness misses every order by more than a
  hundred such units. The metallic table is matched only to 1.1e-3 (`R0`),
  and the integral solver, fed the same state, converges to RCWA's answer and
  away from the table (within 6e-5 of RCWA at 1200 points, within 2.2e-4 of
  the table at the paper's 400); recorded in `docs/theory/rcwa.md` section 9.
- RCWA against published lamellar benchmarks (`tests/test_rcwa_benchmarks.py`),
  each value read in a free primary text and cited with its page: Granet's
  metallic TM `R0 = 0.8484781` (as printed by Lalanne & Hugonin 2000 and in
  Popov's *Gratings*, 2014) to 6.3e-4 at M = 80, run through
  `check_convergence` so the one-call converged path is exercised on it;
  Lalanne & Hugonin's dielectric TM `T-1 = 51.062%` to 5.4e-5 at M = 60; and
  Gralak's silicon rods, p, `R0 = 0.9487` to 6.6e-5 at M = 80. Tolerances come
  from truncation sweeps recorded in the tests. Moharam et al. 1995, Lalanne &
  Morris 1996, Li 1996 and Granet & Guizal 1996 are closed access and no value
  from them is used. Two published values are not reproduced and are recorded
  in `docs/theory/rcwa.md` section 9: Granet's TE `R0`, which this code's
  Littrow order matches instead (and the integral solver agrees with RCWA),
  and Gralak's s `R0`, 0.0100 off.
- Condition estimates for the rigorous solvers' dense systems. Every solve
  and inverse in RCWA's core and the integral method's Nystrom systems goes
  through `solvers/_conditioning.py`, which takes LAPACK's `gecon` 1-norm
  estimate from the same LU factorisation for `O(n^2)` extra; the solutions
  are numpy's to rounding. The worst value per group lands in
  `Provenance.notes["condition"]` (`modes` and `smatrix` for RCWA, `nystrom`
  for the integral method), and each solver warns past its own threshold,
  placed by measurement above everything it reaches on ordinary problems:
  `1e10` for RCWA (ordinary at most 4.4e7) and `1e14` for the integral method
  (ordinary up to 4.5e12 on a corner-graded mesh at 3200 points). Measured at
  a Rayleigh anomaly, both stay accurate far past those values, so the warning
  says a problem is unusual, not that its answer is wrong.
- A Rayleigh-anomaly sweep for both rigorous solvers
  (`tests/test_rayleigh_anomaly.py`): RCWA and the integral method under both
  boundary conditions, approaching an order passing off to 1e-8 nm from each
  side and landing on it. Every efficiency stays finite, energy balance holds
  at the solver's usual accuracy, and the gap either side is a square-root
  branch point, `C sqrt(d / lambda*)` with `C` constant to 4% over six
  decades of `d`, which is what tells continuity from a step.
- The rigorous solvers, `integral` and `rcwa`, cap BLAS at one thread for
  the length of a solve and restore the process's setting on the way out,
  however the solve ends. OpenBLAS's thread-per-core default is slower at
  these matrix sizes, and on a busy machine much slower: RCWA paid 20 to 135x
  and the integral solver 1.1 to 3.7x. `blas_threads=None` opts out. The cap
  is process-wide, as a BLAS pool is, and reference-counted so overlapping
  solves on different threads share it. New dependency: `threadpoolctl`
  (pure Python, BSD-3), the maintained way to reach NumPy's and SciPy's
  separate OpenBLAS copies.
- `workers=N` on both rigorous solvers solves wavelengths on `N` threads.
  The default, 1, is the serial loop. The scan is identical to the serial
  one, and the progress contract holds: every call on the caller's thread,
  `(0, n)` first, `done` stepping by one, and `SolveCancelled` re-raised
  unchanged once the rows already running finish, with nothing left running.
  The integral solver gives each thread its own kernel memo over the scan's
  shared geometry (`OperatorAssembler.fork`).
- RCWA solves an in-plane mount as two decoupled `N`-square systems, TE and
  TM, instead of one coupled `2N` system, and skips a polarization nothing
  excites. The backend also forms the flux at the top plane only, which is
  all `absorption` needs; `solve_stack(planes="all")`, still the core's
  default, gives every plane. Both agree with the paths they skip to 1e-12
  (`tests/test_solver_execution.py`).
- `benchmarks/solver_speed.py` times each case in a fresh interpreter under
  inherited and single-thread BLAS. Before and after, on MC1 (16 cores) with
  another session's solver running alongside, best of 3, one wavelength
  unless noted:

  | case | before, default BLAS | before, 1 thread | after |
  |---|---|---|---|
  | integral, 512 points | 2.94 s | 1.43 s | 1.07 s |
  | integral, 1024 points | 10.93 s | 7.80 s | 9.63 s |
  | integral, tabulated conical, 512 | 3.12 s | 2.66 s | 2.65 s |
  | rcwa lamellar in-plane, M=30 | 4.93 s | 0.057 s | 0.024 s |
  | rcwa lamellar in-plane, M=80 | 15.01 s | 0.733 s | 0.232 s |
  | rcwa lamellar conical, M=30 | 3.16 s | 0.046 s | 0.059 s |
  | rcwa lamellar conical, M=80 | 14.18 s | 0.761 s | 0.947 s |
  | rcwa sinusoid in-plane, M=30, 16 slabs | 64.00 s | 0.470 s | 0.322 s |
  | integral, 8 wavelengths, 512 | 13.10 s | 3.55 s | 5.49 s; 1.87 s at `workers=8` |
  | rcwa conical sinusoid, 8 wavelengths | 324.3 s | 4.99 s | 5.18 s; 0.66 s at `workers=8` |

  "After" is with the inherited environment, so it is what a caller now
  gets without setting anything. The conical and integral rows run the same
  arithmetic as the single-thread column before and differ from it by the
  machine's load, which moved between runs (12% to 25% at the start of each,
  and higher while they ran); read differences under about 30% there as
  noise. The in-plane RCWA gain over the single-thread column, 2 to 3x, is
  the decoupling.

- An RCWA backend, `rcwa`, written in-house: the first working slice of
  roadmap item 2. It is the Fourier modal method with Li's factorization
  rules (the inverse rule for the field component normal to the ridge walls),
  Rumpf's symmetric S-matrix rather than a transmittance matrix, and the full
  vectorial conical formulation, so off-plane incidence and TE-TM conversion
  are exact. `truncation` is the accuracy knob. A lamellar groove is one exact
  slab; any other single-valued profile is a staircase whose `slices` follow
  the truncation by default, and an undercut is refused. A `Problem` means
  what it means to the integral solver's tabulated mode, so efficiencies are
  absolute and `absorption` is the power that enters the structure. The core
  solves a layer stack with films above and below the grating region;
  `Problem` cannot describe one yet, so that is tested at the core only.
  Validated without the corpus: Fresnel on a flat interface to 1e-13, the
  Airy formula for a film to 1e-12, reciprocity falling with truncation, Li's
  rule converging on metallic TM where Laurent's is off by a fifth, and the
  integral solver agreeing order by order (to 5e-7 on a dielectric and 8e-4 in
  metallic TM, the gap shrinking as both refine). It registers, so it appears
  in the conformance suite and the Help menu; it has no GUI tab yet. The
  design note and what is not yet validated are in `docs/theory/rcwa.md`.

### Changed

- `absorption` on a scan is now power absorbed, not power that entered the
  material. For the integral solver's tabulated mode on a lossless material it
  used to equal the transmitted power; that is now `transmission`, and the
  absorption is the small remainder, which is the discretisation error. The
  conserved sum, and so `check_energy_balance`, is unchanged in value. On an
  absorbing material nothing changes.
- Solver tabs carry a short label ("Integral" rather than "Integral
  (boundary integral equation)"), with the full name as the tab's tooltip and
  still in the Help menu. The full names overflowed the tab bar.
- The AFM window names the loaded file once, in its window title, instead of
  in the title of every plot on every tab.
- The Boundary tab's ±1σ band is drawn with edges, and more strongly on a
  dark theme, where it had been all but invisible.

### Fixed

- The theme choice is now in every window. The AFM window had no menu bar
  and so no View > Theme at all, and its plots ignored the theme: in Dark the
  controls went dark and every figure stayed white. Both windows now carry the
  same Theme menu (`theming.ThemeMenu`), checked in step whichever one it was
  used in, and the AFM window's five figures redraw in the new theme, with
  their black lines and grey captions taken from it.

## 0.1.13 (2026-10-27)

### Added

- The PCGrate header on a written `.ggp` is optional: `write_ggp(header=False)`,
  `GGP_WRITE_HEADER` in the metrology config, and a "Write the PCGrate header"
  checkbox on the Boundary tab, on by default. Whether PCGrate needs the header
  rests on the files found on disk rather than its documentation, and this is
  how to test it. Prompted by Aditi, whose script writes `.ggp` files with no
  header.
- Two-column line profiles (`x z`, as Gwyddion's profile export writes them)
  load as a one-line scan, with the width taken from the x column. The same
  surface exported as an image and as a profile averages to the same groove
  (`tests/metrology/test_line_profile.py`). From Aditi's changes to the
  groove-averaging script (`afm_scan_avg_a270.py`), which detect the two
  formats and take the scan width from the file.
- Groove orientation in the AFM pipeline: a sample mounted a degree or two
  off square is measured and, when it matters, corrected, with no rotating by
  eye in Gwyddion. `metrology.core.orientation` finds the grooves band by band,
  links them into tracks and fits one slope to all of them, giving the angle
  and its standard error; the correction shifts each row so the grooves line
  up and takes the pixel along the groove normal, which is exact for a
  grating. On a synthetic 30° grating rotated 3° the blaze angle comes back
  as 29.998° against 29.966° uncorrected, the period 315.00 nm against 315.43,
  and the whole-image groove depth to within 0.2% where uncorrected it reads
  32% shallow. The `groove_orientation` setting is `auto` by default: correct
  at 0.2° or more when that is three standard errors from zero; `always`,
  `detect` and `off` are the others. The angle is reported whether or not it
  was applied, in the analysis result, the results panel, the Import and
  Boundary tabs, a boundary profile's metrics (so the `.json` export) and the
  metrics sidecar. A new wiki page, *Groove orientation*, explains it.

- `load_scan_y_size` reads the slow-axis scan size from `.spm`, `.ibw` and
  Gwyddion text files, scaled to the lines actually recorded, and returns
  `None` rather than assuming a square scan when the file does not record it.
  A scan without one is reported as not measured.
- AFM scan-line scar repair, the step Gwyddion calls "Remove scars": a run of
  one to four lines standing off the lines either side of it is found after
  image flattening and either interpolated down its columns or dropped from row
  grouping (`scar_repair = 'none' | 'interpolate' | 'drop'`, default `none`,
  with `scar_threshold`, `scar_min_length`, `scar_max_width`). It is its own
  stage in `metrology.core.scars`, not an image-flattening method, because it
  is used together with one. The number of lines repaired is shown in the
  Import tab and written to the results panel, the summary file and the
  boundary metrics file. On the synthetic scan with seven lines scarred over a
  third of their length the mean blaze angle reads 29.731° against 29.896°
  clean and the averaged groove is 1.01 nm off; interpolated, 29.895° and
  0.006 nm. Lines are judged levelled by their median step rather than by
  `align_rows`, and a scar's shadow on a rotated grating is set aside.
- A region mask: drag on the image in the AFM window's Import tab to exclude a
  rectangle, whole scan lines or whole columns. It is stored as settings
  (`mask_rect`, `mask_rows`, `mask_cols`, in pixels; `MASK_*` in `config.py`) so
  a batch or CLI run reproduces it, and it is recorded wherever scar repair is.
  Masked lines join no row group and masked columns yield no grooves. Masked,
  a 150 nm particle that moved the mean blaze angle by 0.49 degrees has no
  influence left: the result equals the clean scan's under the same mask.
- Wiki page "Scars and masking".

- Each window carries a **Go** menu: it switches to the other mode and returns
  to the chooser, with `Ctrl+Shift+M` for the latter. The chooser was a
  one-way door — picking a mode closed it with no route onward — which is half
  a fix for the problem it was added to solve.
- Switching reuses a window that is already open instead of constructing a
  second one. The metrology window holds a loaded scan and its settings, so a
  switch that produced an empty copy would read as the work having been lost.
  A window that was closed comes back with its state, since closing does not
  destroy it while the launcher holds a reference.

### Changed

- `write_ggp` no longer ends the file with a line break after the last point.
  PCGrate rejects a `.ggp` that ends in one, as `np.savetxt` output does.
  Found by Aditi.
- The default analysis now corrects a groove rotation of 0.2° or more. A
  synthetic scan with none, including the committed fixture, is unaffected
  (it measures -0.002 ± 0.002° and is left alone). **The real scans have not
  been re-measured with it**, since none are on the machine this was written
  on, so the real-scan numbers quoted in the wiki (the 0.49° flattening span,
  the ICC range) predate it. Set `GROOVE_ORIENTATION = 'detect'` to see what it
  would do to a scan without changing anything.
- The mask is applied before image flattening, so `align_rows` and `plane`
  now skip NaN pixels; a scan with nothing masked levels exactly as before.
- The Import tab's y axis counts scan lines from the top, where line 0 was
  always drawn, so it reads the same line numbers the mask stores. It used to
  count from the bottom.
- With both of the above, a scan is prepared in this order: region mask,
  image flattening, scar repair, groove orientation, tip correction. Scars
  are repaired before the orientation correction shifts rows, and the mask
  is drawn and stored in pixels of the scan as loaded, which is why the
  Import tab draws the image before that correction and the profile after.

### Fixed

- The integral solver read the profile parameter backwards and solved every
  asymmetric groove as its mirror image. It placed the curve at `x = t P`
  where `docs/conventions.md` puts `t` along `-d̂`, so a `Blazed` grating
  sent its blaze to the wrong side of the normal while the scalar solver, on
  the same problem, sent it where `geometry.blaze_direction` says. At
  `Blazed(9)`, α = 40°, 550 nm, p = 1667 nm the integral put 0.56 in m = −1
  and now puts 0.80 in m = +1, beside scalar's 0.89. Symmetric profiles,
  including every sinusoid, are unchanged, which is why the cross-method
  tests never saw it; an asymmetric case now guards it.
- `.ggp` files and the metrology pipeline's `x_norm` advance along `+d̂`,
  which the PCGrate TASTE match fixes, so reading or writing one now reflects
  the points into `t` (`FromProfileData.along_dispersion` and
  `dispersion_axis`). A measured groove therefore reaches both solvers facing
  the way PCGrate sees it. The integral solver's results on a `.ggp` are
  unchanged, to 2e-14, because its two reflections cancel; the scalar
  solver's change, because it had been solving every measured asymmetric
  groove mirrored. `read_profile(path).t` is now `1 − x` reversed rather than
  the file's `x`, and the boundary `.json` keeps storing `x_norm`.

- **The slow-axis scan size ignored `\Aspect Ratio`, so the measured groove
  rotation was too small by that ratio on most scans.** Nanoscope writes the
  same number twice in `Scan Size` even for a rectangular scan and records the
  shape in `\Aspect Ratio`; `io/spm.py::_slow_size_um` took the second value
  at face value. 1013 of the 1245 scans in this project's archive are not
  1:1 -- 911 at 4:1, 68 at 8:1, 26 at 2:1 -- and since a rotation is an
  arctangent of a ratio of the two axes, each of those reported a rotation too
  small by its aspect ratio. Some then declined an orientation correction they
  needed, because the understated angle fell below the threshold.

  **This changes measured numbers.** Any slow-axis length, and anything
  derived from one, differs for a non-1:1 scan: a 1.5 um scan at 4:1 now
  reports 0.375 um of slow axis instead of 1.5. Blaze angles taken from the
  fast axis are unaffected, as are all 1:1 scans -- which is why nothing
  caught this, the TASTE fixture and the synthetic scan both being 1:1.

  Checked against Gwyddion, which reads the field and states the result as
  `# Height:` in its export: of the 70 archive scans with an export beside
  them, 3 agreed before the fix and 70 agree after. The 3 were the 1:1 ones.

  One trap is pinned by its own test. Two different fields are called
  `\Aspect Ratio`: the scan list states the shape as a ratio with a colon,
  while every plane section carries a bare number of the same name that is not
  a shape at all (0.352542 on the TASTE file). Reading the section's value in
  preference -- which is how every other field in this reader works -- turns
  that file's 1.57 um slow axis into 4.44 um.

## 0.1.12 (2026-10-24)

### Added

- Asylum Research `.ibw` files open directly, in the AFM window and through
  `load_afm_data`, with no Gwyddion export step: `metrology.io.ibw` reads Igor
  binary waves (version 5, either byte order) without a new dependency, takes
  the ZSensor channel (Height when there is none) in the chosen scan direction,
  and refuses a channel that is not in metres. **Not yet checked against a file
  from an instrument**: the tests build files to the specification, and a
  synthetic 25° blaze at 400 nm comes back as 24.81° and 400.29 nm through the
  full analysis.

- `gratinglab.measured.MeasuredEfficiency`, measured efficiency with an
  uncertainty on every point: sparse per order, as parallel `order`,
  `wavelength`, `value` and `sigma` float arrays, with a
  `MeasurementProvenance` naming the facility, campaign, polarization and
  mount. It is a separate type rather than an `EfficiencyScan` with error
  fields, because measured data has none of a model's dense-grid shape, and
  the `uncertainties` package does not become a dependency. A polarization
  that was not recorded stays `None` and is never defaulted, and a slightly
  negative value from a weak order is kept as measured.
  `MeasuredEfficiency.from_energies` converts photon energies with the
  constant the materials code uses.
- `compare.align` and `compare.records` take measurements beside scans, so
  measured, scalar and integral results share one `Comparison`.
  `Comparison.measured_order(m)` gives the points with their error bars for
  plotting, and `Comparison.residuals(measurement, method)` the model at each
  measured point with the residual, the pull in units of sigma and a
  `summary()` that reports chi-squared with its point count. The model must
  have been swept on the measured wavelengths (`wavelength_grid()`, or its
  union with a denser grid): nothing interpolates, and an off-grid point is
  refused by name. `sweep` refuses a measurement and says to pass it to
  `align`, since resampling sparse data onto a grid would invent points.
  Existing comparisons are unchanged.
- The PCGrate XML results dialect, ported from `als_eff_meas`:
  `io.efficiency_table.read_xml_table`, `xml_to_scan` and `read_xml_scan`.
  TE, TM and unpolarized all read; TM never could in the reader it came
  from. The file is read by position, as the original did, and its scan
  axis is not read, so the wavelengths come from the caller and must pair one
  to one with the scan points. An order a point does not list is recorded as
  non-propagating there.
- `tests/test_measured.py` holds a worked comparison against the ported ALS
  campaign (OGRE laminar grating I, October 2021). It needs `als_eff_meas`
  and the private beamline tree, and skips without them; neither MC1 nor
  MAC1 holds that tree, so it has not yet run.

## 0.1.11 (2026-10-24)

### Added

- `corner_grading=p` on the integral solver: a corner-graded boundary mesh
  (Kress's sigmoidal transformation, M&P §4.6.6) for perfectly conducting
  gratings with corners. Cornered TM converged at first order on the
  equal-arc mesh; with `corner_grading=4` its energy-balance deviation is about
  1e-5 at 256 points on a lamellar, an ideal sawtooth and a 29.5°/70.5° blaze,
  against 2e-3 to 3e-2 before, and reciprocity improves by the same factor. It
  is opt-in, the default results are unchanged bit for bit, and a profile with
  no corner keeps the equal-arc mesh whatever is passed. It serves
  `conductivity="tabulated"` too, where `R + A` reaches 2e-6 to 3e-5 at 512
  points on a visible-light metal and on conical X-ray lamellar and blazed
  gratings, against 5e-3 to 0.16 on the equal-arc mesh. That needed two
  corrections the perfectly conducting solve does not use, a singularity
  subtraction of the double layer at corners and the tangential derivative
  taken on the assembled trace. The measurements are in
  `docs/theory/integral.md` §5 and `docs/findings.md`.
- `BoundaryCurve.polygon`: a graded mesh carries the vertices it lies on.
- `Profile.boundary(n, grading=p)`, which every built-in profile accepts, and
  `BoundaryCurve.graded`. A profile written against the older `boundary(n)`
  signature still works with the default mesh.
- The Setup tab acts on the material library instead of only describing it.
  **Add material…** picks a CXRO, refractiveindex.info or PCGrate table, asks
  its coating name, reads it to prove it is usable and copies it into your
  library; **Rescan library** re-reads the directories after a change made
  outside the window; **Open library folder** opens your library, creating it
  first. Both changes reach the geometry dock's Coating dropdown at once, which
  before needed a restart, and a selected coating that disappears falls back to
  "(none)" and says so. The page's links now open in the system browser, and
  its provenance section cites the kinds of table actually installed rather
  than Henke for every one.
- `materials.install(path, name=None, overwrite=False)`, the same install from
  Python. An unreadable table is refused before it reaches the library.
- Themes: **View → Theme** offers System (follows the operating system's light
  or dark setting, live), Light, Dark, Solarized Light, Solarized Dark and High
  contrast, remembered between runs and shared by both windows. Light is the
  look the windows had before, colour for colour. Every theme's text meets
  WCAG AA (4.5:1) and every plotted line 3:1, which `tests/test_gui_theme.py`
  asserts. The palette, every plot, the provenance panel, the Setup and theory
  pages and the geometry views all restyle in place.

- A **Corner mesh** selector on the Integral tab: off (the equal-arc mesh,
  still the default) or a grading order from 2 to 6, with a note on what it
  does. The form checks the graded mesh's own resolution floor, whose
  coarsest spacing is about twice the mean, and refuses an order that would
  put nodes closer than double precision separates, each on the field that
  caused it, using the solver's own check.

### Changed

- Diffraction orders are coloured by sign, cool for negative, warm for
  positive and neutral for m = 0, instead of from matplotlib's ten-colour
  cycle, which drew a sixteen-order scan with repeats (m = -3 and m = +7 were
  the same orange). An order keeps its colour when others are hidden, and the
  order list shows each line's colour beside its entry. The legend moved
  outside the axes, where it no longer covers the curves it names.
- The geometry dock's four groups share one label column, so every input
  starts at the same offset (they sat up to 64 px apart), and switching mount
  no longer moves them.
- The Grating Geometry tab: the 3D box fills more of its panel (a uniform
  zoom, so angles stay exact) with five ticks a side instead of nine; the
  cross-section's order labels are pushed outward along their own rays until
  none overlaps another or the title, where they used to stack at the top of
  the arc; the 3D cone-axis label moved off the tip where the orders
  cluster; and the Wavelength box is two compact rows with the scan's ends
  shown, giving the captions beside it the width.

## 0.1.10 (2026-10-23)

### Added

- `gratinglab-app`, a chooser that asks which window you want, and the single
  application bundle that launches it. Both halves have always been this one
  package, but the desktop layer shipped two console scripts and so two icons,
  the AFM one still called "AFM Blaze Meas" with the bundle id
  `org.jamccoy.afmblazemeas`, all of it inherited from the repository the
  metrology was merged from. Spotlight indexed it as `AFMBlazeMeas`, which
  sorts nowhere near GratingLab and names a project that no longer exists, so
  the half of the package in regular use was the half that could not be found
  — including by its author, looking for the AFM analysis inside the
  efficiency window. `tools/apps.py bundle --app gratinglab` now builds one
  `GratingLab.app` (`org.gratinglab.app`) that opens the chooser, and the
  chooser opens either window in the same process. The two direct scripts are
  unchanged.
- Bundles that skip the chooser are still available as `--app efficiency` and
  `--app metrology`, renamed "GratingLab Efficiency" and "GratingLab AFM" with
  `org.gratinglab.*` identities. Neither is built by default, so the Dock gets
  one icon.

- `flatten_profile` takes `return_report=True` and returns a
  `FlatteningReport` beside the flattened profile: the method and polynomial
  order actually used, how many features the background was fitted through,
  whether edge exclusion was skipped, and the RMS scatter of those features
  about the background they produced. The keyword follows
  `find_groove_positions`' `return_n_edge_rejected` and is off by default, so
  two-value unpacking keeps working and the GUI needed no change.
- The blaze-angle analysis asks for a report per row group and the per-sample
  summary records the aggregate. `analysis_data_*.csv` gains `Flatten_method`,
  `Flatten_fallback_groups`, `Flatten_residual_nm` and
  `Flatten_residual_groups`, appended so earlier columns keep their positions,
  and the text summary says whether every row group got the method asked for.
  A row group levelled by a fallback can now be found in a file afterwards
  instead of only in the console of the run that produced it.
- `Flatten_residual_nm` is the new diagnostic: how far the lands lie from the
  background fitted through them, in nanometres, and so how much the plane
  every blaze angle is quoted against deserves to be trusted. It is the worst
  row group's rather than the mean, since the weakest reference is the one that
  matters. It is absent for `linear` and `polynomial`, whose residual would
  describe the groove depth rather than the reference, and groups whose fit had
  no spare features are excluded and counted separately — three lands and a
  quadratic give zero by arithmetic, not by flatness.

### Changed

- `flatten_profile`'s fallbacks raise a `UserWarning` instead of printing to
  standard output. A print scrolls past in a batch run with nothing able to
  catch it; a warning can be turned into an error, counted, or routed to a log
  by whatever drives the pipeline. The suite runs with `filterwarnings =
  error`, so no test can now pass on a background fitted by a method it did not
  ask for — which is also how this change was checked: the full suite still
  passes, so nothing in it was flattening by fallback.

## 0.1.9 (2026-10-23)

### Changed

- The PyPI keywords no longer lead with x-ray work or name PCGrate. Nothing in
  the solvers or the metrology assumes a wavelength band — the grating equation
  carries none, and the materials layer reads CXRO x-ray tables and
  refractiveindex.info ultraviolet-to-infrared tables alike — so the old list
  told a reader the in-plane half of the package was not meant for them.
  `x-ray` is still there, among the bands rather than ahead of them, and
  `optics`, `diffraction grating`, `diffraction efficiency`, `integral method`,
  `groove profile`, `ultraviolet` and `infrared` are new. `pcgrate` is gone:
  it is a program this project validates against, not a feature of this one.
  `CITATION.cff` gains the same balance.
- Three classifiers added: `Operating System :: OS Independent`, which CI
  demonstrates on three platforms, and `Topic :: Scientific/Engineering ::
  Information Analysis` and `:: Visualization` beside Physics. There is no
  trove classifier for optics, and no Qt environment classifier is claimed,
  because a classifier cannot be conditional on an extra and a default install
  pulls no toolkit.
- The README names [endstation](https://pypi.org/project/endstation/) as what
  gratinglab is built on and links its PyPI page rather than only its
  repository, in the dependency paragraph and in the licence table. A reader on
  gratinglab's project page who wants the shared scaffolding wants the package,
  not the source tree.
- The README's licence section, which is the PyPI project page, was the single
  line "BSD-3-Clause." It now states the copyright, what the licence permits
  (commercial and closed-source use, no copyleft reach), what it requires
  (notice retention, and the no-endorsement clause), and that there is no
  warranty — which matters more than usual while the solvers are experimental.
  It adds how to cite a version, a table of the licences each extra brings in,
  the LGPL position of the `gui` extra and why PySide6 rather than PyQt6, and
  the two things the licence does not cover: optical constant tables, whose
  terms belong to whoever published them, and the measurement data, which is
  not ours to publish.

- The metrology window's Flattening page explains the levelling it documents
  rather than only tabulating it: why a background has to be removed at all
  (the substrate plane is the reference a blaze angle is quoted against, and
  an error in it is an error in every angle), why fitting through the lands is
  not the same as fitting through the whole trace, and where the proof that
  image flattening is free stops holding. It also documents the third
  levelling, the endpoint flattening applied on PCGrate export, which was
  described only in `docs/theory/metrology.md`; the page previously said there
  were two stages.
- The same page now states the `level_grooves` feature-detection criteria
  (half a period of separation, 20% of the height range in prominence), the
  shipped defaults including `FLATTEN_FEATURE = 'peaks'`, the fallbacks taken
  when too few features are found, and the fact that `level_grooves` needs a
  period estimate, so the period used to level a profile is not independent of
  the period the run reports. `tests/metrology/test_wiki.py` pins each of
  these against the code.

### Fixed

- `flatten_profile`'s `level_grooves` method no longer fits its background to
  whatever survived edge exclusion while reporting that it had used every
  feature. Edge exclusion runs after feature detection, and on a trace whose
  features sit near the ends it could leave one point, or none; the arrays had
  already been overwritten, so the fit used the remnant and the message saying
  otherwise was wrong. The exclusion is now skipped, with an accurate warning,
  whenever it would leave fewer features than the requested order needs.
- The same method no longer asks numpy for a polynomial the feature set cannot
  determine. The guard admitted two features while `FLATTEN_POLY_ORDER`
  defaults to 2, and a degree-k fit needs k+1 points, so the background came
  back rank-deficient behind a `RankWarning` that nothing was reading. The
  order is reduced to what the features support and the reduction is reported.
  Both paths produced a background that looked like a fit and was not one,
  which matters because flattening sets the plane every blaze angle on the
  scan is quoted against. `tests/metrology/test_flatten_profile.py` covers
  them; the four tests that pin the fixes fail against the previous code.

- The Flattening page's method table marked `level_grooves` as the default in
  a row that also listed its three feature choices, which read as though
  `both` were the default feature. It is `peaks`.

- The README and the docstrings name `tools/apps.py`, the script that builds
  the icons and the macOS bundles. They still named `tools/make_app.py` and
  `tools/make_icon.py`, which were merged into it, so every command the
  documentation gave for building a bundle or regenerating an icon failed. The
  README also now says that a bundle execs the console script in the
  checkout's virtualenv rather than freezing one, so it does not survive being
  copied to another machine, and names `tools/apps.py check`, which reports
  that breakage.

## 0.1.8 (2026-10-22)

### Added

- This changelog, linked from the README and the PyPI project page, and a
  versioning policy in `docs/roadmap.md` saying what 0.2 and 1.0 will mean.
  `tests/test_changelog.py` fails if the current version has no section.

### Changed

- The wiki's generated Home page matches the README's opening and no longer
  says its pages come from `docs/` alone.

### Fixed

- On Python 3.11 the `gui` extra no longer installs PySide6 6.12.0, which
  aborts the interpreter at exit (`Fatal Python error: bool_dealloc`) after
  heavy use of the windows. Versions 6.9 to 6.11 are unaffected, and on Python
  3.12 and later the bug is harmless, so the exclusion applies to that one
  release on 3.11 only.

## 0.1.7 (2026-10-22)

### Changed

- The README, every page under `docs/`, the metrology window's help pages and
  the docstrings are rewritten in plain prose. No code changed: each module's
  syntax tree with its docstrings blanked is identical to 0.1.6's.
- The README describes the AFM metrology as the part in regular use and marks
  the scalar and integral solvers as experimental, since they pass their test
  suite but have not yet been used for real analysis.
- Figures derived from PCGrate's output are no longer quoted in the public
  documentation. The comparisons against it are still described, as against an
  independent integral-method code on the private reference corpus.

### Fixed

- `docs/metrology/README.md` said the SEM of the row-group analysis scales as
  1/√N and that the row groups are treated as independent. Comparisons use the
  ICC-corrected SEM and the effective sample size, and the page now says so.
- The disagreement between the fitted facet angle and the depth-implied blaze
  angle of the TASTE scan was quoted as 3.6σ in three places, although
  `findings.md` explains that sigma is the wrong measure for it. It is now
  given as the 7.6° gap.
- Four docstrings were wrong: an unclosed code span in `gui/integral_options`,
  a claim in `gui/qt/scalar_tab` that no second solver tab exists, "both
  actions" for four buttons in the Boundary tab, and a cross-reference to a
  renamed section.

## 0.1.6 (2026-10-22)

### Added

- Asymmetric tip shapes. Each flank of the tip has its own half angle
  (`tip_half_angle_right_deg`; `None` keeps the symmetric cone), so a pyramidal
  probe's cross-section in the scan plane can be modelled.
  `gratinglab.metrology.core.tip.TIP_PRESETS` holds Bruker's ScanAsyst-Air
  scanned across the cantilever and along it in both orientations, and the
  metrology window has a Tip shape picker.
- The GitHub wiki also carries the metrology window's help pages, published as
  `AFM-<slug>`.

### Fixed

- Tip erosion now uses the reflected tip. The difference is invisible for a
  symmetric tip; with an asymmetric one the previous form gave a looser bound
  whose re-imaging did not reproduce the image.

## 0.1.5 (2026-10-22)

### Added

- Swath statistics (`gratinglab.metrology.boundary.swath_statistics`): the
  averaged groove is built separately for each horizontal band of scan rows
  and compared, reporting the shape spread between swaths, the groove-to-groove
  spread within one, per-swath metrics, the lateral shift down the scan, and how
  far the whole-scan groove lies from the swaths' mean. The Boundary tab has a
  "Swath statistics…" dialog with a per-swath CSV export.

## 0.1.4 (2026-10-22)

### Added

- A versioned boundary format, `gratinglab.io.boundary`: `.json` carrying the
  period, groove count, metrics and provenance, plus a `.csv` table. The
  Boundary tab exports `.ggp`, `.json` or `.csv`, and "Send to efficiency
  window" opens the groove there with its period.
- Tip correction controls in the metrology window.

### Fixed

- The Qt windows run Python's cycle collector on the GUI thread only. Garbage
  collection on a worker thread could free Qt objects there, which segfaulted on
  macOS and aborted the interpreter at exit on Python 3.11.

## 0.1.3 (2026-10-22)

### Fixed

- Groove averaging skips grooves that the scan edge cuts into
  (`average_grooves(..., require_full_period=True)`). A partial groove at the
  edge of a scan had been included in the average and distorted the exported
  shape.

## 0.1.2 (2026-10-22)

### Changed

- The version is written once, in `gratinglab.__version__`; `pyproject.toml`
  reads it from there.

## 0.1.1 (2026-10-21)

### Changed

- Requires endstation 0.2.2, which fixed its own encoding warnings.
- The README shows how to install with uv and links to PyPI, the
  documentation, the releases and the issue tracker.

### Fixed

- The critical angle is reported only where the refractive index is near 1,
  and is computed exactly there.

## 0.1.0 (2026-10-21)

First release on PyPI.

- AFM groove metrology (`gratinglab.metrology`): blaze angles per groove from
  Nanoscope `.spm` or Gwyddion text scans, row-group statistics with the
  correlation (ICC) correction, an averaged groove exported as a PCGrate `.ggp`
  boundary, and tip erosion with a certainty map.
- Scalar and integral-method efficiency solvers, the latter with perfect and
  finite conductivity, a comparison harness that also accepts imported
  efficiency tables, and a convergence harness.
- A materials layer reading CXRO exports, refractiveindex.info files and
  PCGrate `.ari` tables from user libraries; no optical constant tables are
  distributed.
- Resolving power from the measured groove count.
- Qt windows for the efficiency solvers and for the metrology.
