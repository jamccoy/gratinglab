#!/usr/bin/env python
"""Build a macOS .app that launches the GUI.

    .venv/bin/python tools/make_app.py [--dest ~/Applications]

Not a frozen bundle -- no PyInstaller, no embedded interpreter. The bundle is
three small files: a plist, a shell stub that execs the console script, and the
icon. That is enough for a Dock icon and a double-clickable launcher, and it
takes a second to rebuild.

The bundle is **not committed**: the stub hard-codes an absolute path to this
checkout's virtualenv, which would be wrong on any other machine. Committing
the generator instead keeps it reproducible.
"""

from __future__ import annotations

import argparse
import plistlib
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ICON = REPO / "resources" / "GratingLab.icns"
LAUNCHER = REPO / ".venv" / "bin" / "gratinglab-gui"
BUNDLE_ID = "org.gratinglab.gui"


def build(destination: Path) -> Path:
    if sys.platform != "darwin":
        raise SystemExit(".app bundles are macOS-only")
    if not LAUNCHER.exists():
        raise SystemExit(
            f"{LAUNCHER} not found -- run: .venv/bin/pip install -e '.[dev,gui]'"
        )
    if not ICON.exists():
        raise SystemExit(f"{ICON} not found -- run: python tools/make_icon.py")

    app = destination / "GratingLab.app"
    contents = app / "Contents"
    if app.exists():
        shutil.rmtree(app)
    (contents / "MacOS").mkdir(parents=True)
    (contents / "Resources").mkdir(parents=True)

    version = _version()
    (contents / "Info.plist").write_bytes(
        plistlib.dumps(
            {
                "CFBundleName": "GratingLab",
                "CFBundleDisplayName": "GratingLab",
                "CFBundleExecutable": "GratingLab",
                "CFBundleIconFile": "GratingLab",
                "CFBundleIdentifier": BUNDLE_ID,
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

    stub = contents / "MacOS" / "GratingLab"
    stub.write_text(f'#!/bin/sh\nexec "{LAUNCHER}" "$@"\n')
    stub.chmod(0o755)

    shutil.copy2(ICON, contents / "Resources" / "GratingLab.icns")

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
        "--dest", type=Path, default=REPO / "build",
        help="where to write GratingLab.app (default: ./build)",
    )
    args = parser.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)

    app = build(args.dest)
    print(f"wrote {app}")
    print(f"  launches: {LAUNCHER}")
    print(f"\nopen '{app}'   # or drag it to /Applications")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
