# refractiveindex.info fixtures

Verbatim copies of thirteen files from the refractiveindex.info database,
<https://github.com/polyanskiy/refractiveindex.info-database>, at commit
`c5c2f188e848453def5970e347399d653df2ffc2` (2026-09-04). The database is
dedicated to the public domain under CC0 1.0, so committing them needs no
permission; its `database/doc/permission.txt` says so.

One file per dispersion formula, the smallest the database has of each, plus
the three tabulated shapes and two well-known glasses with published indices to
check against. The name says what each exercises; `reference.json` gives each
one's path in the database.

`reference.json` holds the expected values, and they come from the database's
own code rather than from gratinglab's. For a formula entry they are the
formula block of the database's `database/tools/nkexplorer.py` run on five
log-spaced wavelengths inside the entry's range. For a tabulated entry they are
the middle row of the file. So the test compares the reader against an
implementation it does not share code with. The same comparison over every
formula entry in the database agreed to 1.4e-15 on 2026-10-21.

Do not edit these files. A fixture that has been tidied is no longer the
format the reader has to handle.
