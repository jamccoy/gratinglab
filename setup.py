"""Copy the Help pages into the built package.

Everything else about the build is declared in pyproject.toml; this file
exists only for the step setuptools cannot declare. The GUI's Help menu reads
``docs/conventions.md`` and ``docs/theory/*.md``, and those live at the
repository root, outside the package, because about two hundred comments,
docs and tests name them by that path. So ``docs/`` stays the one copy, and
``build_py`` copies the pages into ``gratinglab/_docs/`` in the wheel, in the
same layout, so their relative links still resolve. ``MANIFEST.in`` puts them
in the sdist, which is what a wheel is built from.

An editable install copies nothing: the GUI finds ``docs/`` by walking up
from its own source, as it always has in a checkout.
"""

import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

DOCS = Path(__file__).parent / "docs"

#: The pages the Help menu shows, relative to docs/.
HELP_PAGES = ("conventions.md", "theory/*.md")


class BuildWithHelpPages(build_py):
    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            return
        target = Path(self.build_lib) / "gratinglab" / "_docs"
        copied = 0
        for pattern in HELP_PAGES:
            for page in sorted(DOCS.glob(pattern)):
                destination = target / page.relative_to(DOCS)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(page, destination)
                copied += 1
        if not copied:
            # A wheel without them would show "not bundled" on every page.
            raise RuntimeError(f"no Help pages found under {DOCS}")


setup(cmdclass={"build_py": BuildWithHelpPages})
