#!/usr/bin/env python
"""Publish ``docs/`` and the metrology window's own pages into the GitHub wiki.

    git clone https://github.com/jamccoy/gratinglab.wiki.git ../gratinglab.wiki
    .venv/bin/python tools/publish_wiki.py --wiki ../gratinglab.wiki

The copying and the page-link rewriting live in :mod:`endstation.tools`, shared
with the other projects here. What remains here is specific to this project, and
it addresses one problem: a wiki is flat, and this project's docs are not.

``docs/`` has subdirectories (``theory/``) and its pages link out of the docs
tree entirely (``../src/gratinglab/convergence.py``). The shared tool publishes
top-level pages and deliberately leaves every other link alone, because it
cannot know whether the target exists in the wiki. In a wiki those links resolve
against ``/wiki/`` and 404. After the copy, this script therefore rewrites
every repository-relative target to an absolute blob URL on the default branch,
which is where that file lives and where a reader following the link expects
to arrive.

The metrology window carries a wiki of its own, in
``src/gratinglab/metrology/wiki/``, shipped as package data so an installed copy
has it. Those pages are the ones a person analysing a scan reads (tip
correction, swath statistics, reading the outputs), so they are published too,
as ``AFM-<slug>`` pages in the window's reading order. The prefix is needed
because a wiki names a page after its file, and a bare ``overview`` or
``flattening`` next to the project pages would say nothing about what it
covers. Their source stays inside the package, so an edit there is still one
edit.

``Home.md`` and ``_Sidebar.md`` are generated instead of copied. They are
navigation pages whose whole content is a list of the pages that exist, which
this script knows and which no document in ``docs/`` should have to duplicate.

Nothing here commits or pushes, following the same rule as the shared tool. It
writes a checkout and stops, so the diff can be read before anything is public.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from endstation.tools.publish_wiki import publish

REPO = Path(__file__).resolve().parent.parent
BLOB = "https://github.com/jamccoy/gratinglab/blob/master"
METROLOGY_WIKI = REPO / "src" / "gratinglab" / "metrology" / "wiki"

#: Prefix for the metrology window's pages, see the module docstring.
AFM_PREFIX = "AFM-"

#: A markdown inline link with a local target: not an image, not a URL, not a
#: bare anchor. Anchors on the target are preserved.
_LOCAL_LINK = re.compile(
    r"(?<!\!)\[([^\]]*)\]\((?!https?://|#|mailto:)([^)\s#]+)(#[^)\s]*)?\)"
)

#: The wiki landing page. Deliberately short: it says what the project is and
#: hands the reader to the page that answers their question, because a landing
#: page that tries to be the documentation is the copy that goes stale.
HOME = """# gratinglab

gratinglab is a Python package for measuring the groove profiles of
diffraction gratings from AFM scans and for computing the diffraction efficiency
of a grating from such a profile. The AFM metrology is the part in regular use;
the efficiency solvers are experimental.

These pages are published from [`docs/`]({blob}/docs) and from the metrology
window's help pages in the repository, which is the source of truth. Edit them
there; edits made here are overwritten by the next publish.

## Start here

| If you want to | Read |
|---|---|
| Know where the project is and what is next | **[roadmap](roadmap)** |
| Write a solver, or file a bug about a sign | **[conventions](conventions)** (normative) |
| See what has actually been established, with evidence | **[findings](findings)** |
| Find the literature behind a method | **[references](references)** |
| Run or read a mutation-testing sweep | **[mutation-testing](mutation-testing)** |
| Analyse an AFM scan: the metrology window's own guide | **[AFM overview](AFM-overview)**, then the AFM pages in the sidebar |
| Correct for the tip, or pick a ScanAsyst-Air preset | **[AFM tip correction](AFM-tip-correction)** |
| See how much the groove shape varies down a scan | **[AFM swath statistics](AFM-swath-statistics)** |

## The shape of it

```
AFM scan  ->  BoundaryProfile   (gratinglab.metrology)
                      |          measured period, averaged groove, facet fit
                      v
Problem       (period + profile + materials; no solver fields)
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
    """Point repository-relative links at the repository instead of ``/wiki/``.

    A target naming a published page is left alone, because the shared tool
    has already made it a wiki link. Everything else is a path relative to
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


def _metrology_pages() -> list[tuple[str, str, str]]:
    """``(wiki name, title, markdown)`` for each metrology page, in reading order.

    Read through ``gratinglab.metrology.wiki.pages``, so the order is the
    window's and a page the window would skip is skipped here too.
    """
    from gratinglab.metrology.wiki import pages

    return [(AFM_PREFIX + page.slug, page.title, page.markdown) for page in pages()]


def _publish_metrology(wiki: Path, dry_run: bool) -> list[Path]:
    """Copy the metrology pages in, comparing first so a dry run means "differs"."""
    changed = []
    for name, _title, markdown in _metrology_pages():
        target = wiki / f"{name}.md"
        text = markdown.replace("\r\n", "\n")
        if target.exists() and target.read_text(encoding="utf-8").replace("\r\n", "\n") == text:
            continue
        changed.append(target)
        if not dry_run:
            target.write_text(text, encoding="utf-8", newline="\n")
    return changed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wiki", type=Path, required=True,
                        help="the gratinglab.wiki checkout")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    docs = REPO / "docs"
    changed = publish(docs, args.wiki, dry_run=args.dry_run)
    pages = {path.stem for path in sorted(docs.glob("*.md"))}
    metrology = _metrology_pages()
    clash = pages & {name for name, _, _ in metrology}
    if clash:
        raise SystemExit(f"docs/ and the metrology wiki both publish {sorted(clash)}")
    changed += _publish_metrology(args.wiki, args.dry_run)

    if not args.dry_run:
        for name in sorted(pages):
            page = args.wiki / f"{name}.md"
            text = page.read_text(encoding="utf-8")
            page.write_text(_to_blob_urls(text, pages), encoding="utf-8")

        (args.wiki / "Home.md").write_text(HOME.format(blob=BLOB), encoding="utf-8")
        (args.wiki / "_Sidebar.md").write_text(
            "### gratinglab\n\n"
            "- [Home](Home)\n"
            + "".join(f"- [{name}]({name})\n" for name in sorted(pages))
            + "\n### AFM analysis\n\n"
            + "".join(f"- [{title}]({name})\n" for name, title, _ in metrology)
            + f"\n- [Repository]({BLOB.rsplit('/blob/', 1)[0]})\n",
            encoding="utf-8",
        )

    print(f"{len(changed)} page(s) {'would change' if args.dry_run else 'written'}")
    if not args.dry_run:
        print(f"review and commit in {args.wiki}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
