#!/usr/bin/env python
"""Publish ``docs/`` into the GitHub wiki checkout.

    git clone https://github.com/jamccoy/gratinglab.wiki.git ../gratinglab.wiki
    .venv/bin/python tools/publish_wiki.py --wiki ../gratinglab.wiki

The copying and the page-link rewriting live in :mod:`endstation.tools`, shared
with the other projects here. What stays is the part that is genuinely ours, and
it is all one problem: **a wiki is flat and this project's docs are not.**

``docs/`` has subdirectories (``theory/``) and its pages link out of the docs
tree entirely (``../src/gratinglab/convergence.py``). The shared tool publishes
top-level pages and deliberately leaves every other link alone, because it
cannot know whether the target exists in the wiki. In a wiki those links resolve
against ``/wiki/`` and 404. So after the copy, this rewrites every
repository-relative target to an absolute blob URL on the default branch --
which is where that file actually lives, and where a reader following the link
wants to end up.

``Home.md`` and ``_Sidebar.md`` are generated rather than copied. They are
navigation, not prose: their whole content is a list of the pages that exist,
which is knowledge this script has and no document in ``docs/`` should have to
duplicate.

Nothing here commits or pushes -- same rule as the shared tool. It writes a
checkout and stops, so the diff can be read before anything is public.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from endstation.tools.publish_wiki import publish

REPO = Path(__file__).resolve().parent.parent
BLOB = "https://github.com/jamccoy/gratinglab/blob/master"

#: A markdown inline link with a local target: not an image, not a URL, not a
#: bare anchor. Anchors on the target are preserved.
_LOCAL_LINK = re.compile(
    r"(?<!\!)\[([^\]]*)\]\((?!https?://|#|mailto:)([^)\s#]+)(#[^)\s]*)?\)"
)

#: The wiki landing page. Deliberately short: it says what the project is and
#: hands the reader to the page that answers their question, because a landing
#: page that tries to be the documentation is the copy that goes stale.
HOME = """# gratinglab

An open platform for diffraction grating groove metrology and efficiency
analysis — from an AFM scan of a real grating to its diffraction efficiency,
without leaving the package.

These pages are published from [`docs/`]({blob}/docs) in the repository, which
is the source of truth. Edit them there; edits made here are overwritten by the
next publish.

## Start here

| If you want to | Read |
|---|---|
| Know where the project is and what is next | **[roadmap](roadmap)** |
| Write a solver, or file a bug about a sign | **[conventions](conventions)** — normative |
| See what has actually been established, with evidence | **[findings](findings)** |
| Find the literature behind a method | **[references](references)** |
| Run or read a mutation-testing sweep | **[mutation-testing](mutation-testing)** |

## The shape of it

```
AFM scan  ->  BoundaryProfile   (gratinglab.metrology)
                      |          measured period, averaged groove, facet fit
                      v
Problem       (period + profile + materials — no solver fields, ever)
Illumination  (wavelength + direction + polarization)
                      |
                      v
Solver (plugin)  ->  Result (+ provenance)  ->  resolving power
```

One `Problem` reaches every solver, so a disagreement between methods is a
disagreement about physics rather than about how two codes were configured.

Theory pages sit in the repository rather than here, one per method, each
stating its assumptions and the conditions under which it stops being
trustworthy: [scalar]({blob}/docs/theory/scalar.md),
[integral]({blob}/docs/theory/integral.md),
[metrology]({blob}/docs/theory/metrology.md).

Installation, the GUI, and where measurement data lives are in the
[README]({blob}/README.md).
"""


def _to_blob_urls(text: str, pages: set[str]) -> str:
    """Point repository-relative links at the repository, not at ``/wiki/``.

    A target naming a published page is left alone -- the shared tool has
    already made it a wiki link. Everything else is a path relative to
    ``docs/``, resolved against it and emitted as an absolute blob URL.
    """

    def replace(match: re.Match[str]) -> str:
        label, target, anchor = match.group(1), match.group(2), match.group(3) or ""
        if "/" not in target and target.removesuffix(".md") in pages:
            return match.group(0)
        # Relative to docs/, since that is where the page was written.
        resolved = Path("docs", target).as_posix()
        while "/../" in resolved:
            resolved = re.sub(r"[^/]+/\.\./", "", resolved, count=1)
        return f"[{label}]({BLOB}/{resolved}{anchor})"

    return _LOCAL_LINK.sub(replace, text)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wiki", type=Path, required=True,
                        help="the gratinglab.wiki checkout")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    docs = REPO / "docs"
    changed = publish(docs, args.wiki, dry_run=args.dry_run)
    pages = {path.stem for path in sorted(docs.glob("*.md"))}

    if not args.dry_run:
        for name in sorted(pages):
            page = args.wiki / f"{name}.md"
            page.write_text(_to_blob_urls(page.read_text(), pages))

        (args.wiki / "Home.md").write_text(HOME.format(blob=BLOB))
        (args.wiki / "_Sidebar.md").write_text(
            "### gratinglab\n\n"
            "- [Home](Home)\n"
            + "".join(f"- [{name}]({name})\n" for name in sorted(pages))
            + f"\n- [Repository]({BLOB.rsplit('/blob/', 1)[0]})\n"
        )

    print(f"{len(changed)} page(s) {'would change' if args.dry_run else 'written'}")
    if not args.dry_run:
        print(f"review and commit in {args.wiki}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
