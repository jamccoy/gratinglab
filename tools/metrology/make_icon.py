#!/usr/bin/env python
"""Build the macOS .icns from the drawn artwork.

    .venv/bin/python tools/make_icon.py

Renders a 1024 px master via :mod:`afm_analysis.gui.icon`, downsamples it with
``sips`` into the iconset macOS expects, then packs it with ``iconutil``. Both
tools ship with macOS, so there is nothing to install.

The resulting ``.icns`` is committed -- it is small, stable, and it *is* the
design. Re-run this only when the motif changes.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUTPUT = REPO / "resources" / "afm_blaze_meas.icns"

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
    if sys.platform != "darwin":
        print("iconutil is macOS-only; skipping .icns build", file=sys.stderr)
        return 0

    sys.path.insert(0, str(REPO / "src"))
    from afm_analysis.gui.icon import render

    with tempfile.TemporaryDirectory() as scratch:
        master = render(Path(scratch) / "master.png", size=1024)
        iconset = Path(scratch) / "afm_blaze_meas.iconset"
        iconset.mkdir()

        for pixels, name in VARIANTS:
            subprocess.run(
                ["sips", "-z", str(pixels), str(pixels), str(master),
                 "--out", str(iconset / name)],
                check=True, capture_output=True,
            )

        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["iconutil", "-c", "icns", str(iconset), "-o", str(OUTPUT)],
            check=True, capture_output=True,
        )

        # Keep a PNG alongside for READMEs and non-macOS contexts.
        preview = OUTPUT.with_suffix(".png")
        subprocess.run(
            ["sips", "-z", "512", "512", str(master), "--out", str(preview)],
            check=True, capture_output=True,
        )

    print(f"wrote {OUTPUT.relative_to(REPO)} ({OUTPUT.stat().st_size:,} bytes)")
    print(f"wrote {preview.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
