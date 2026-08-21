#!/usr/bin/env python
"""Build the macOS .icns from the drawn artwork.

    .venv/bin/python tools/make_icon.py [--app metrology]

Renders a 1024 px master via each window's own ``gui.icon`` module, downsamples it with
``sips`` into the iconset macOS expects, then packs it with ``iconutil``. Both
tools ship with macOS, so there is nothing to install.

The resulting ``.icns`` is committed -- it is small, stable, and it *is* the
design. Re-run this only when the motif changes.

Two windows, one builder, matching ``make_app.py``. The artwork differs; the
sips/iconutil dance does not, and keeping two copies of it was how a fix could
reach one icon and not the other.
"""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

#: key -> (icon module, output filename). Kept in step with ``make_app.py``.
APPS = {
    "gratinglab": ("gratinglab.gui.icon", "GratingLab.icns"),
    "metrology": ("gratinglab.metrology.gui.icon", "afm_blaze_meas.icns"),
}

#: macOS expects exactly these, and iconutil is fussy about the names.
VARIANTS = [
    (16, "icon_16x16.png"),
    (32, "icon_16x16@2x.png"),
    (32, "icon_32x32.png"),
    (64, "icon_32x32@2x.png"),
    (128, "icon_128x128.png"),
    (256, "icon_128x128@2x.png"),
    (256, "icon_256x256.png"),
    (512, "icon_256x256@2x.png"),
    (512, "icon_512x512.png"),
    (1024, "icon_512x512@2x.png"),
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--app", choices=sorted(APPS), default="gratinglab",
        help="which window's artwork to render (default: gratinglab)",
    )
    args = parser.parse_args()
    module_name, icns_name = APPS[args.app]
    output = REPO / "resources" / icns_name

    if sys.platform != "darwin":
        print("iconutil is macOS-only; skipping .icns build", file=sys.stderr)
        return 0

    sys.path.insert(0, str(REPO / "src"))
    render = importlib.import_module(module_name).render

    with tempfile.TemporaryDirectory() as scratch:
        master = render(Path(scratch) / "master.png", size=1024)
        iconset = Path(scratch) / f"{output.stem}.iconset"
        iconset.mkdir()

        for pixels, name in VARIANTS:
            subprocess.run(
                ["sips", "-z", str(pixels), str(pixels), str(master),
                 "--out", str(iconset / name)],
                check=True, capture_output=True,
            )

        output.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["iconutil", "-c", "icns", str(iconset), "-o", str(output)],
            check=True, capture_output=True,
        )

        # Keep a PNG alongside for READMEs and non-macOS contexts.
        preview = output.with_suffix(".png")
        subprocess.run(
            ["sips", "-z", "512", "512", str(master), "--out", str(preview)],
            check=True, capture_output=True,
        )

    print(f"wrote {output.relative_to(REPO)} ({output.stat().st_size:,} bytes)")
    print(f"wrote {preview.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
