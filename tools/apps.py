#!/usr/bin/env python
"""The two windows this project ships, and how to build a bundle for either.

    .venv/bin/python tools/apps.py icon --app metrology
    .venv/bin/python tools/apps.py bundle --app gratinglab [--dest ~/Applications]
    .venv/bin/python tools/apps.py check --app gratinglab

The bundling itself lives in :mod:`endstation.tools`, which is shared with the
other projects here. What stays is the part that is genuinely ours: which
windows exist, what they are called, and which artwork belongs to each.

It used to be two builders, then one builder with two tables that a comment had
to keep in step -- ``make_app.py`` holding names and bundle ids, ``make_icon.py``
holding icon modules and filenames, each listing the same two keys. Now one
table holds an entry per window, so adding a third window is one edit and there
is no second place to forget.

Bundles are **not committed**: the stub hard-codes an absolute path to this
checkout's virtualenv, which would be wrong on any other machine.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from endstation.tools.make_app import App, build_app, check_app
from endstation.tools.make_icon import build_icns

REPO = Path(__file__).resolve().parent.parent


@dataclass(frozen=True)
class Window:
    """One window's identity, across both the icon and the bundle."""

    key: str
    name: str
    stem: str
    icon_module: str
    icns: str
    script: str
    bundle_id: str

    @property
    def icon(self) -> Path:
        return REPO / "resources" / self.icns

    @property
    def launcher(self) -> Path:
        return REPO / ".venv" / "bin" / self.script

    def spec(self) -> App:
        return App(
            name=self.name,
            stem=self.stem,
            launcher=self.launcher,
            bundle_id=self.bundle_id,
            icon=self.icon if self.icon.exists() else None,
            version=_version(),
        )


WINDOWS = {
    window.key: window
    for window in (
        Window("gratinglab", "GratingLab", "GratingLab",
               "gratinglab.gui.icon", "GratingLab.icns",
               "gratinglab-gui", "org.gratinglab.gui"),
        Window("metrology", "AFM Blaze Meas", "AFMBlazeMeas",
               "gratinglab.metrology.gui.icon", "afm_blaze_meas.icns",
               "gratinglab-metrology-gui", "org.jamccoy.afmblazemeas"),
    )
}


def _version() -> str:
    sys.path.insert(0, str(REPO / "src"))
    try:
        from gratinglab import __version__

        return __version__
    except Exception:
        return "0.0.0"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("icon", "bundle", "check"))
    parser.add_argument("--app", choices=sorted(WINDOWS), default="gratinglab",
                        help="which window (default: gratinglab)")
    parser.add_argument("--dest", type=Path, default=REPO / "build",
                        help="where to write the .app (default: ./build)")
    args = parser.parse_args(argv)

    window = WINDOWS[args.app]

    if args.action == "icon":
        icns = build_icns(window.icon_module, window.icon)
        print(f"wrote {icns.relative_to(REPO)} ({icns.stat().st_size:,} bytes)")
        return 0

    if args.action == "check":
        problems = check_app(args.dest / f"{window.stem}.app")
        for problem in problems:
            print(f"  - {problem}")
        print("looks correct" if not problems else "needs rebuilding")
        return 1 if problems else 0

    if not window.icon.exists():
        raise SystemExit(
            f"{window.icon} not found -- run: "
            f"python tools/apps.py icon --app {window.key}"
        )

    args.dest.mkdir(parents=True, exist_ok=True)
    bundle = build_app(window.spec(), args.dest)
    print(f"wrote {bundle}")
    print(f"  launches: {window.launcher}")
    print(f"\nopen '{bundle}'   # or drag it to /Applications")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
