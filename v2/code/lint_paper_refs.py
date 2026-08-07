# -*- coding: utf-8 -*-
"""
Structural lint for the v2 manuscript: every \\ref resolves, every
\\cite has a \\bibitem, no label is orphaned, and the shared statement
counter is not written by hand anywhere.

The counter check matters because wall_v2.tex declares one counter that
every theorem, conjecture, corollary, proposition, lemma AND remark
advances. Writing "Proposition 15" by hand goes stale the moment a
remark is inserted; v1's PROVENANCE.md drifted by one from conj:L
onward for exactly this reason.

Usage: python lint_paper_refs.py [file.tex]
Exit status is the number of problems.
"""

import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "paper", "wall_v2.tex")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    src = io.open(path, encoding="utf-8", newline="").read()

    # \eqref counts as a reference; equation labels are reached that way.
    refs = set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", src))
    labs = set(re.findall(r"\\label\{([^}]*)\}", src))
    cites = set()
    for grp in re.findall(r"\\cite\{([^}]*)\}", src):
        cites.update(x.strip() for x in grp.split(","))
    bibs = set(re.findall(r"\\bibitem\{([^}]*)\}", src))

    # Statements written by hand: "Theorem 3", "Proposition 15", ... A
    # number attributed to a named source ("Huang--Li's Lemma 1") is not
    # this document's counter and is not flagged.
    kinds = "Theorem|Conjecture|Corollary|Proposition|Lemma|Remark"
    handwritten = []
    for m in re.finditer(r"(?<!\\)\b(" + kinds + r")~?\s(\d+)", src):
        back = src[max(0, m.start() - 26): m.start()]
        if re.search(r"'s\s*$|\\cite\{[^}]*\}[^.]{0,12}$", back):
            continue
        handwritten.append(f"{m.group(1)} {m.group(2)}")

    print("lint_paper_refs")
    print(f"target: {os.path.relpath(path, os.path.join(HERE, '..', '..'))}"
          f"   ({len(src.splitlines())} lines)")
    print("=" * 66)

    problems = []
    # counts_as_problem: an unused anchor or an unused background
    # reference breaks nothing; a dangling pointer or a hand-written
    # counter does.
    for name, items, fatal in (
        ("dangling \\ref (no matching \\label)", sorted(refs - labs), True),
        ("dangling \\cite (no \\bibitem)", sorted(cites - bibs), True),
        ("statement number written by hand", handwritten, True),
        ("orphaned \\label (never referenced)", sorted(labs - refs), False),
        ("uncited \\bibitem", sorted(bibs - cites), False),
    ):
        status = "ok" if not items else f"{len(items)} " + (
            "FOUND" if fatal else "(not fatal)")
        print(f"  {name:<42} {status}")
        for it in items:
            print(f"      - {it}")
        if fatal:
            problems.extend(items)

    print()
    print(f"  labels: {len(labs)}   refs: {len(refs)}   "
          f"cites: {len(cites)}   bibitems: {len(bibs)}")
    print(f"  problems: {len(problems)}")
    print()
    print("  note: orphaned labels and uncited bibitems are listed but not")
    print("  counted -- an unused anchor or a background reference breaks")
    print("  nothing. Dangling pointers and hand-written statement numbers")
    print("  do, and are counted.")
    return len(problems)


if __name__ == "__main__":
    sys.exit(main())
