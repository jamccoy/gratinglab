# Vendored optical constants — provenance

Every file here is a **verbatim CXRO/Henke export**, not a re-derivation. They
are committed in their original format so that the reader shipped in
`optical.py` is the one exercised in anger, and so the header line carries the
density the tables were computed at.

| File | Material | Range | Points | Origin |
|---|---|---|---|---|
| `Au.txt` | Gold, density 19.32 g/cm³ | 200 eV to 2 keV (0.620 to 6.199 nm) | 500 | CXRO, <https://henke.lbl.gov/optical_constants/> |

That row said 6.171 nm and 498 points until 2026-10-21, and both were wrong:
1239.84 / 200 eV is 6.199 nm. The request that produced the file is recovered
rather than remembered. <https://henke.lbl.gov/optical_constants/getdb2.html>
with formula `Au`, density `-1` (tabulated), `Scan=Energy`, `Min=200`,
`Max=2000` and `Npts=500` returns it value for value.

## Citation

The underlying data is the Henke–Gullikson–Davis compilation. Work using these
constants should cite:

> B.L. Henke, E.M. Gullikson, and J.C. Davis, *X-ray interactions:
> photoabsorption, scattering, transmission, and reflection at E = 50–30000 eV,
> Z = 1–92*, Atomic Data and Nuclear Data Tables **54** (2), 181–342 (1993).

## Redistribution: none, for now

CXRO makes these tables freely available for scientific use and they are widely
redistributed, but the terms have not been confirmed for bundling inside a
BSD-3 package, and the site states no data licence. So **the wheel and sdist
ship no table**, decided 2026-10-21 before the first release. `Au.txt` stays in
the repository because the test suite needs a real table, and `release.yml`
fails a build whose wheel or sdist carries one.

An install therefore starts with no materials, and `lookup` refuses with a
message naming the way out. Nothing in `gratinglab` requires a table to exist;
`test_an_empty_data_directory_is_a_supported_state` holds it to that.

## Adding a material

For yourself, not here. Put the table in your own library,
`~/.gratinglab/materials` or a directory named in `GRATINGLAB_MATERIALS`, as
`<Name>.txt` for a CXRO export, `<Name>.yml` for a refractiveindex.info file or
`<Name>.ari` for PCGrate's format. That is
searched before this directory, survives an upgrade, and is never published.
See `materials/__init__.py`.

This directory holds only what the test suite needs and is never packaged, so
a table belongs here only if a test needs it, and only with its row above.

Do **not** hand-edit these files. If a table needs trimming or resampling, do it
at load time: an edited copy stops being the thing this file claims it is.
