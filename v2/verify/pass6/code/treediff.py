# -*- coding: utf-8 -*-
"""Do the two evidence trees still agree?

deploy/code/ and deploy/results/ are copies of a subset of v2/code/ and
v2/results/.  If a file was regenerated on one side only, the deployed
paper and the working paper are quoting different runs of the same
script.  The commits alternate between the trees, so this is exactly
where a divergence would sit.
"""
import hashlib
import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"z:\업무\goldbach-ln2-real"
PAIRS = [("deploy/code", "v2/code"), ("deploy/results", "v2/results")]

# A deployed file with no counterpart in v2/code or v2/results is not
# automatically orphaned: the verification passes keep their own code and
# results under v2/verify/passN/, and the deploy tree draws on those too.
# Every "in deploy only" file is looked up there before being reported.
FALLBACK = ["v2/verify/pass%d/%s" % (n, sub)
            for n in range(1, 7) for sub in ("code", "results")]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def listing(rel):
    base = os.path.join(ROOT, rel.replace("/", os.sep))
    out = {}
    for d, _, fs in os.walk(base):
        if "__pycache__" in d:
            continue
        for f in sorted(fs):
            p = os.path.join(d, f)
            out[os.path.relpath(p, base).replace(os.sep, "/")] = sha(p)
    return out


print("STATISTIC: sha256 agreement, file by file, between the deployed "
      "evidence tree and the working evidence tree it was copied from")
print("FIELD:     deploy/code vs v2/code; deploy/results vs v2/results")
print("CONSTANTS: none")
print("NULL:      the deployed tree is a byte-identical subset of the "
      "working tree; any file present on both sides with different "
      "content means the two papers quote different runs")
print("DENOM:     every file in the deployed tree")
print()

bad = 0
for dep, src in PAIRS:
    D, S = listing(dep), listing(src)
    print("=" * 70)
    print("%s  (%d files)   vs   %s  (%d files)" % (dep, len(D), src,
                                                    len(S)))
    print("=" * 70)
    only = [k for k in sorted(D) if k not in S]
    diff = [k for k in sorted(D) if k in S and D[k] != S[k]]
    same = [k for k in sorted(D) if k in S and D[k] == S[k]]
    print("  identical            : %d" % len(same))
    print("  DIFFERENT content    : %d" % len(diff))
    print("  in deploy only       : %d" % len(only))
    for k in diff:
        print("     DIFF  %-40s %s -> %s"
              % (k, S[k][:12], D[k][:12]))
        bad += 1
    for k in only:
        found = []
        for fb in FALLBACK:
            p = os.path.join(ROOT, fb.replace("/", os.sep),
                             k.replace("/", os.sep))
            if os.path.exists(p):
                same = sha(p) == D[k]
                found.append("%s (%s)" % (fb, "identical" if same
                                          else "DIFFERENT"))
                if not same:
                    bad += 1
        print("     ONLY  %-34s %s"
              % (k, "; ".join(found) or "NOT ANYWHERE IN v2/"))
        if not found:
            bad += 1
    print()

print("total files differing or unaccounted for: %d" % bad)
