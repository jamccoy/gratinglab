# Changelog

Notable changes, newest first. The project is pre-1.0, so a minor version may
change an interface; any such change is stated here when it is made. The
versioning policy, including what 0.2 and 1.0 will mean, is in
[`docs/roadmap.md`](docs/roadmap.md#versioning).

The format loosely follows [Keep a Changelog](https://keepachangelog.com). The
version is read from `gratinglab.__version__`, and a release commit moves the
entries under "Unreleased" into a section named for the new version;
`tests/test_changelog.py` fails if it does not.

## Unreleased

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

### Changed

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
