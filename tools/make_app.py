#!/usr/bin/env python
"""Build a macOS .app that launches the GUI.

    .venv/bin/python tools/make_app.py [--app metrology] [--dest ~/Applications]

Not a frozen bundle -- no PyInstaller, no embedded interpreter. The bundle is
three small files: a plist, a shell stub that execs the console script, and the
icon. That is enough for a Dock icon and a double-clickable launcher, and it
takes a second to rebuild.

The bundle is **not committed**: the stub hard-codes an absolute path to this
checkout's virtualenv, which would be wrong on any other machine. Committing
the generator instead keeps it reproducible.

There are two windows in this project and one builder. It used to be two
builders -- the metrology package brought its own copy, differing in four
constants and a filename -- which is exactly the arrangement that lets a fix
land in one and not the other.
"""

from __future__ import annotations

import argparse
import plistlib
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


class App:
    """One window's bundle identity.

    ``stem`` carries no spaces: it names the executable, the icon file and the
    iconset directory inside the bundle, where a space is a nuisance rather
    than an error.
    """

    def __init__(self, key, name, stem, icon, script, bundle_id):
        self.key, self.name, self.stem = key, name, stem
        self.icon = REPO / "resources" / icon
        self.launcher = REPO / ".venv" / "bin" / script
        self.bundle_id = bundle_id


APPS = {
    app.key: app
    for app in (
        App("gratinglab", "GratingLab", "GratingLab", "GratingLab.icns",
            "gratinglab-gui", "org.gratinglab.gui"),
        App("metrology", "AFM Blaze Meas", "AFMBlazeMeas", "afm_blaze_meas.icns",
            "gratinglab-metrology-gui", "org.jamccoy.afmblazemeas"),
    )
}


def build(destination: Path, spec: "App") -> Path:
    if sys.platform != "darwin":
        raise SystemExit(".app bundles are macOS-only")
    if not spec.launcher.exists():
        raise SystemExit(
            f"{spec.launcher} not found -- run: "
            ".venv/bin/pip install -e '.[dev,gui,metrology]'"
        )
    if not spec.icon.exists():
        raise SystemExit(
            f"{spec.icon} not found -- run: "
            f"python tools/make_icon.py --app {spec.key}")

    app = destination / f"{spec.stem}.app"
    contents = app / "Contents"
    if app.exists():
        shutil.rmtree(app)
    (contents / "MacOS").mkdir(parents=True)
    (contents / "Resources").mkdir(parents=True)

    version = _version()
    (contents / "Info.plist").write_bytes(
        plistlib.dumps(
            {
                "CFBundleName": spec.name,
                "CFBundleDisplayName": spec.name,
                "CFBundleExecutable": spec.stem,
                "CFBundleIconFile": spec.stem,
                "CFBundleIdentifier": spec.bundle_id,
                "CFBundlePackageType": "APPL",
                "CFBundleShortVersionString": version,
                "CFBundleVersion": version,
                # Not a guess and not a preference: PySide6-Essentials ships
                # as pyside6_essentials-6.11.1-cp310-abi3-macosx_13_0_universal2,
                # so the GUI extra cannot install below macOS 13. Claiming 11.0
                # would let Finder launch a bundle whose window can never open.
                "LSMinimumSystemVersion": "13.0",
                # Without this the window renders blurry on a Retina display.
                "NSHighResolutionCapable": True,
            }
        )
    )

    stub = contents / "MacOS" / spec.stem
    stub.write_text(f'#!/bin/sh\nexec "{spec.launcher}" "$@"\n')
    stub.chmod(0o755)

    shutil.copy2(spec.icon, contents / "Resources" / f"{spec.stem}.icns")

    # Finder caches icons aggressively; touching the bundle nudges it.
    subprocess.run(["touch", str(app)], check=False)
    return app


def _version() -> str:
    sys.path.insert(0, str(REPO / "src"))
    try:
        from gratinglab import __version__

        return __version__
    except Exception:
        return "0.0.0"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--app", choices=sorted(APPS), default="gratinglab",
        help="which window to bundle (default: gratinglab)",
    )
    parser.add_argument(
        "--dest", type=Path, default=REPO / "build",
        help="where to write the .app (default: ./build)",
    )
    args = parser.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)

    spec = APPS[args.app]
    app = build(args.dest, spec)
    print(f"wrote {app}")
    print(f"  launches: {spec.launcher}")
    print(f"\nopen '{app}'   # or drag it to /Applications")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
