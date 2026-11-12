# Optical constants: none here, and where to get them

This directory holds no tables, in the repository or in a release. gratinglab
commits only data that is its own or in the public domain, and no published
optical constants table is clearly either. CXRO's site states no data licence,
and refractiveindex.info's files, though dedicated CC0, report other people's
measurements.

Until 2026-11-11 the repository carried one, `Au.txt`: a verbatim CXRO/Henke
export of gold, 200 eV to 2 keV in 500 points, used by the test suite. It never
shipped in a wheel or sdist, and it remains in the git history and in the
source archives GitHub generates for releases up to 0.3.3. The suite now uses
an invented gold-sized table, `tests/data/materials/Au.txt`, written by
`make_au.py` beside it, and `coating="Au"` resolves to that file only inside
the tests.

An install therefore starts with no materials, and `lookup` refuses with a
message naming the way out. The one exception is not a table:
`coating="typical:Au"` and its siblings are an approximate X-ray model computed
from physical facts (atomic numbers and weights, core-level binding energies
from the X-ray Data Booklet, first ionization energies, densities), each stated
with its source in `materials/typical.py`. See `docs/typical-materials.md`. Nothing in `gratinglab` requires a table to exist;
`test_an_empty_data_directory_is_a_supported_state` holds it to that, and
`test_the_package_data_directory_has_no_table` holds this directory empty.

## Adding a material

Put the table in your own library, `~/.gratinglab/materials` or a directory
named in `GRATINGLAB_MATERIALS`, as `<Name>.txt` for a CXRO export, `<Name>.yml`
for a refractiveindex.info file or `<Name>.ari` for PCGrate's format. That is
searched before this directory, survives an upgrade, and is never published.
See `materials/__init__.py`.

| Kind | Where | Range |
|---|---|---|
| CXRO/Henke export | <https://henke.lbl.gov/optical_constants/>, `getdb2.html`, one material, energy scan | X-ray |
| refractiveindex.info | <https://refractiveindex.info>, the YAML file for a dataset | ultraviolet to infrared |

A table you add is governed by whoever published it, not by gratinglab's
licence.

## Citation

Work using CXRO tables should cite the compilation behind them:

> B.L. Henke, E.M. Gullikson, and J.C. Davis, *X-ray interactions:
> photoabsorption, scattering, transmission, and reflection at E = 50–30000 eV,
> Z = 1–92*, Atomic Data and Nuclear Data Tables **54** (2), 181–342 (1993).
