# -*- coding: utf-8 -*-
"""
Which of this tree's findings does wall_v1_corrected.tex still carry?
(v1_verify2, Phase 2.)

`v1_verify/paper/wall_v1_corrected.tex` repairs the first pass's fifteen
findings. This tree's Phase 1 added fifteen more. The question this
script answers is the one a reader of the corrected paper needs
answered: is it clean?

Each entry below is a literal form from `FINDINGS.md`, in the LaTeX the
paper writes it in. A HIT means the corrected paper still carries the
defective form. Regexes are escaped for LaTeX -- the first pass's own
withdrawn-form auditor missed `0.39\\%` because it looked for `0.39%`,
which is finding #7's guard failure, so the same mistake is avoided here
by matching on the TeX source.

Exit status is the number of surviving findings.
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
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DEFAULT = os.path.join(ROOT, "v1_verify", "paper", "wall_v1_corrected.tex")
TARGET = sys.argv[1] if len(sys.argv) > 1 else DEFAULT

# (id, severity, one-line description, regex over the TeX source)
CHECKS = [
    ("A1", "high", "conj:wall item 2: gap +-0.0056 quoted as forty s.e.",
     r"0\.2238\\pm0\.0056"),
    ("M2", "high", "conj:wall item 1: -0.0005 on 6.3e6, cell means alone",
     r"Excess kurtosis \$-0\.0005\$"),
    ("M3", "high", "conj:wall item 3: the t=5 ratio and the extremes claim",
     r"extremes are attained at generic"),
    ("M14", "high", "sec:floor decay table: depth-0 exponent 0.6289+-0.0121",
     r"0\.6289.*\n?.*0\.0121|0\.6289"),
    ("M14b", "high", "sec:floor: chi^2/dof = 251 rejects a common exponent",
     r"chi\^2/\\mathrm\{dof\} = 251"),
    ("A2", "med", "sec:floor: 'steps at 5 to 30 standard errors'",
     r"\$5\$ to \$30\$ standard errors"),
    ("M7", "med", "sec:margin: the margin at N=1e8 is N^{0.454}",
     r"N\^\{0\.454\}"),
    ("M4", "med", "sec:coin: major-arc factors 8.40 and 15.16",
     r"8\.40.*\n?.*15\.16|15\.16"),
    ("M5", "med", "sec:coin: the 1.051-1.068 excess over sqrt(0.32264(X-h))",
     r"0\.32264"),
    ("A4", "low", "sec:closures: the undefined 'C4' and its 8.8%",
     r"C4 threshold"),
    # A6 is the ABSENCE of the completion step, not the presence of a
    # phrase: the defect is asserting m < N^{1-theta'} over k<K without
    # first completing the divisor sum. So it fires only when the
    # mechanism paragraph is there AND no completion language is.
    ("A6", "low", "thm:A: the one-line mechanism, completion step missing",
     ("squares itself away", r"complet(e|es|ing|ion)|complementary sum")),
    ("A7", "low", "prop:E: '>> N by Parseval'",
     r"\\gg N\$ by Parseval"),
    ("A8", "low", "rem:rho: the 1e8 conversion 0.810 -> 0.841",
     r"0\.841"),
    ("A9", "low", "prop:E: four values at seven abscissae, 'decaying'",
     r"below \$1\$ and decaying"),
    ("M10", "low", "lem:placebo: 'a fixed permutation of the label set'",
     r"fixed permutation of the label set"),
    ("A3", "low", "abstract: 'ten kill-tested technique designs'",
     r"ten\s+kill-tested"),
    ("M6", "low", "sec:coin: the shift shares, at an unstated X",
     r"48\.9\\%"),
]


def main():
    src = io.open(TARGET, encoding="utf-8", newline="").read()
    flat = re.sub(r"\s+", " ", src)

    print("lint_corrected_vs_findings   (v1_verify2 Phase 2)")
    print(f"target: {os.path.relpath(TARGET, ROOT)} "
          f"({len(src.splitlines())} lines)")
    print("=" * 74)
    print()
    print(f"  {'id':<6}{'sev':<7}{'still present?':<16}finding")
    print("  " + "-" * 86)
    alive = []
    for fid, sev, desc, pat in CHECKS:
        if isinstance(pat, tuple):
            # (form that must be present, language that cures it). The cure
            # must appear in the SAME passage, not anywhere in the file --
            # `wall_v1.tex` says "the complete divisor sum" three sections
            # later, about a different object, and a global search would
            # read that as a repair.
            need, cure = pat
            m = re.search(need, src)
            hit = False
            if m:
                window = src[max(0, m.start() - 1400): m.start()]
                hit = re.search(cure, window) is None
        else:
            hit = bool(re.search(pat, src)
                       or re.search(re.sub(r"\\n\?", "", pat), flat))
        if hit:
            alive.append((fid, sev, desc))
        print(f"  {fid:<6}{sev:<7}{'YES -- unfixed' if hit else 'no -- repaired':<16}{desc}")
    print()
    n_high = sum(1 for f, s, _ in alive if s == "high")
    print(f"  surviving: {len(alive)} of {len(CHECKS)}  "
          f"({n_high} of them high severity)")
    print(f"  repaired : {len(CHECKS) - len(alive)}")
    print()
    print("  For contrast, the first pass's own findings ARE repaired in")
    print("  this file; that is what it was written to do. What it does")
    print("  not carry is any repair for the findings it never made.")
    return len(alive)


if __name__ == "__main__":
    sys.exit(main())
