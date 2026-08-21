# -*- coding: utf-8 -*-
"""BRIEF sec.2.4 / sec.7: did pass4's confirmed fixes reach BOTH trees?

pass4 was blind and read only the projection.  Its findings were then
acted on.  This script asks, finding by finding, whether the corrected
text is present in deploy/papers/ and whether it is present in
v2/paper/ -- because the commits alternate between the two trees and a
fix applied to one is not a fix applied to the other.

Each probe is a pair (WANT, STALE): WANT is a string the corrected text
must contain, STALE a string only the uncorrected text contains.  The
probes are written from pass4/FINDINGS.md and both are printed, so the
screen is auditable and can be re-run.
"""
import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"z:\업무\goldbach-ln2-real"
DEP = ["deploy/papers/P1-mobius-fixed-class.tex",
       "deploy/papers/P2-no-go-divisor-switch.tex",
       "deploy/papers/P3-wall-second-moment.tex",
       "deploy/papers/P4-coherent-cell-floor.tex",
       "deploy/papers/P5-negative-map.tex"]
SRC = ["v2/paper/theorem_A.md", "v2/paper/wall_v3.md"]


def rd(rel):
    return io.open(os.path.join(ROOT, rel.replace("/", os.sep)),
                   encoding="utf-8", errors="replace").read()


D = "\n".join(rd(r) for r in DEP)
S = "\n".join(rd(r) for r in SRC)

# (id, subject, WANT = corrected marker, STALE = uncorrected marker)
P = [
 ("F1", "P2 eq.(5): the T_w term restored",
  "T_w", None),
 ("F2", "P2 Prop 21(ii): binary Goldbach called a classical sieve bound",
  None, "by classical sieve bounds"),
 ("F3", "P2 meas:direct: the two ratio lists share a numerator",
  "2.1519", "2.1591"),
 ("F4", "P1 eq.(9): B(N) missing the squarefree restriction",
  "\\mu^2(k)\\,(\\log k)\\,\\bigl|\\Emu(N;k)\\bigr|", None),
 ("F5", "P3 sec.5: Gumbel margin figures assume rho=1",
  "10^{22.842}", None),
 ("F6", "P5: route count (seventeen vs eighteen)",
  "eighteen", "seventeen"),
 ("F7", "P4 meas:scaleinv: two sample resolutions in one table",
  "+0.012633", "-0.000879"),
 ("F8", "P4: 'conservative by about two orders of magnitude'",
  None, "two orders of magnitude"),
 ("F9", "P4: the floor/sd ratio range (5.8-160 -> 7.3-158)",
  "158", "5.8$ to $160"),
 ("F10", "P4: the factor in N (140 -> 128)",
  "=128", "factor $140$"),
 ("F11", "P3: 'about 430 standard deviations'",
  None, "430"),
 ("F14", "P3: 0.6455 comes from a file not in the packet",
  "0.6455", None),
 ("F15", "P5 sec.6(2): the middle HB share (0.833180 -> 0.840039)",
  "0.840039", "0.833180"),
 ("F17", "P1: the comparison quoted as asymp 10^{-3}N",
  None, "10^{-3}N"),
 ("F18", "P3: A=0.787275 called the generic even-N value",
  "0.787275", "generic even $N"),
 ("F19a", "P4 meas:mc: 'agree to six digits'",
  None, "six digits"),
 ("F19b", "P4: 'the same 0.0013 throughout'",
  None, "same $0.0013$ throughout"),
 ("F19c", "P2 meas:direct: '2% of S' vs 2% of N",
  None, "$2\\%$ of $\\SS$"),
 ("F19d", "P2: sqrt(6/pi^2)=0.7797 reproduced 'exactly'",
  "0.7798", None),
]

print("STATISTIC: presence of each pass4 fix-marker, and of the "
      "corresponding stale marker, in the projection and in the source")
print("FIELD:     deploy/papers/*.tex (5 files) and v2/paper/*.md (2)")
print("CONSTANTS: %d probes, each a (want, stale) string pair taken "
      "from v2_log/verify/pass4/FINDINGS.md and printed below" % len(P))
print("NULL:      a fix that was acted on is present in BOTH trees and "
      "the stale marker is absent from both; anything else is a "
      "divergence between the two trees")
print("DENOM:     %d probes" % len(P))
print()
hdr = "  %-5s %-8s %-8s %-8s %-8s  %s" % (
    "id", "want/D", "want/S", "stale/D", "stale/S", "subject")
print(hdr)
print("  " + "-" * (len(hdr) - 2))
rows = []
for pid, subj, want, stale in P:
    wd = ("-" if want is None else ("yes" if want in D else "NO"))
    ws = ("-" if want is None else ("yes" if want in S else "NO"))
    sd = ("-" if stale is None else ("STALE" if stale in D else "gone"))
    ss = ("-" if stale is None else ("STALE" if stale in S else "gone"))
    rows.append((pid, wd, ws, sd, ss, subj, want, stale))
    print("  %-5s %-8s %-8s %-8s %-8s  %s" % (pid, wd, ws, sd, ss, subj))

print()
print("  probe strings:")
for pid, subj, want, stale in P:
    print("    %-5s want=%-46r stale=%r" % (pid, want, stale))

print()
print("=" * 74)
print("DIVERGENCES (fix in one tree only, or stale marker surviving)")
print("=" * 74)
n = 0
for pid, wd, ws, sd, ss, subj, want, stale in rows:
    bad = []
    if wd == "yes" and ws == "NO":
        bad.append("corrected in the PROJECTION only")
    if wd == "NO" and ws == "yes":
        bad.append("corrected in the SOURCE only")
    if sd == "STALE" and ss == "gone":
        bad.append("stale text survives in the PROJECTION only")
    if sd == "gone" and ss == "STALE":
        bad.append("stale text survives in the SOURCE only")
    if sd == "STALE" and ss == "STALE":
        bad.append("stale text survives in BOTH trees")
    if bad:
        n += 1
        print("  %-5s %s" % (pid, subj))
        for b in bad:
            print("         -> %s" % b)
print("  total probes showing a divergence or a surviving marker: %d" % n)
