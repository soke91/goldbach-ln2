# -*- coding: utf-8 -*-
"""Do the source papers' verify-pass script pointers resolve to files?

rem:maskstamp, rem:cellstamp and the c02 citation point at
verify/passN/code/*.py.  v2_log/verify/ is where the passes live in this
working copy.  This asks whether those paths hold anything.
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"z:\업무\goldbach-ln2-real"
PAT = re.compile(r"verify[/\\]pass\d[/\\]code[/\\][A-Za-z0-9_\\]+\.py")

print("STATISTIC: whether every verify-pass script pointer printed in "
      "v2/paper/*.md resolves to a file in this repository")
print("FIELD:     v2/paper/theorem_A.md and wall_v3.md; candidate roots "
      "are v2_log/, v2/ and the repository root")
print("CONSTANTS: none")
print("NULL:      a citable pointer resolves; a pointer to a path that "
      "holds no file is a reproducibility gap")
print("DENOM:     every such pointer in the two source papers")
print()

tot = 0
for f in ["theorem_A.md", "wall_v3.md"]:
    t = io.open(os.path.join(ROOT, "v2", "paper", f),
                encoding="utf-8").read()
    for i, line in enumerate(t.split("\n")):
        for m in PAT.finditer(line):
            tot += 1
            raw = m.group(0).replace("\\_", "_").replace("\\", "/")
            rel = raw.replace("/", os.sep)
            cands = [os.path.join(ROOT, "v2_log", rel),
                     os.path.join(ROOT, "v2", rel),
                     os.path.join(ROOT, rel)]
            ok = [c for c in cands if os.path.exists(c)]
            print("  %-13s:%-5d %-50s %s"
                  % (f, i + 1, raw, "RESOLVES" if ok else "NO SUCH FILE"))

print()
print("  pointers found: %d" % tot)
print()
print("  v2_log/verify/ contents:")
base = os.path.join(ROOT, "v2_log", "verify")
for d in sorted(os.listdir(base)):
    print("    %-8s %s" % (d, sorted(os.listdir(os.path.join(base, d)))))
