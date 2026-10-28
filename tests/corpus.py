"""Where the reference corpus lives, for the tests that need it directly.

Separate from ``conftest.py`` on purpose. This is an ordinary helper that
three test modules import at module level, and ``conftest.py`` is not an
ordinary module: pytest imports it under a name of its own choosing, which
depends on the import mode and on whether ``tests/`` is a package. Importing
*from* it therefore works under a plain ``pytest`` run and breaks under tools
that drive pytest in-process from a copied tree -- mutmut is one, and that is
how this file came to exist.

``conftest.py`` keeps the fixtures and the Qt environment, which is what it is
for.
"""

from pathlib import Path

from endstation.paths import data_root

DATA = Path(__file__).parent / "data"

#: The corpus is reference data belonging to the research group and is
#: deliberately not committed. Point ``GRATINGLAB_REF_DIR`` elsewhere to
#: override.
DEFAULT_REF_DIR = Path.home() / "Documents" / "diffraction_efficiency"


def reference_dir() -> Path | None:
    """The PCGrate reference corpus, or ``None`` if unavailable.

    ``None`` rather than a skip, because the modules that import this decide a
    ``skipif`` at collection; a test that wants the directory itself takes the
    ``ref_dir`` fixture, which skips naming the variable.
    """
    candidate = data_root("GRATINGLAB_REF_DIR", DEFAULT_REF_DIR)
    return candidate if candidate.is_dir() else None
