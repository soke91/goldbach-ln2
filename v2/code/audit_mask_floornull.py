# -*- coding: utf-8 -*-
r"""
"Nothing can be fitted there" -- is that true where the amplitude never clears?

WHAT IS AT STAKE

OPEN.md item 3 says the mask has no decay exponent because **at the
three shallowest depths the amplitude does not clear the exact floor
at any scale measured, so nothing can be fitted there.**  The premise
is exactly right: over the fifteen half-octave bands of
v1/results/wall/lab_mask_amplitude_law.txt no band reaches |z| = 2 at
depth 0 or depth 1, and only two do at depth 2, where the deeper
cells clear it in every band.

The conclusion does not follow from the premise, and the data says
why.  At depth 0 the amplitude runs 0.1005, 0.0847, 0.0632, 0.0512,
0.0411 across the first five bands -- a smooth fall, not a scatter --
while its error barely moves.  **A quantity can decay systematically
across fifteen bands without any single band distinguishing it from
zero**, and the fit sees the fifteen while the floor test sees them
one at a time.

So the question is whether the exponents quoted at those depths --
0.6289 at depth 0 and 0.3686 at depth 2, both at more than fifty
standard errors -- are the trend or the fitting of |noise|.  Under
pure noise dm ~ N(0, se^2) the fitted exponent has an expectation:
|dm| tracks se, so the exponent comes out at the exponent of se
itself, which {#rem:maskdmse} measured as 0.0420 at depth 0.  The
observed 0.6289 is fifteen times that.  The null makes the comparison
properly.

Nothing is measured.  Every value is read from v1's table.

BACKS: Remark {#rem:maskfloornull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  The control.  The exponents refitted here reproduce the six
      that v1/results/wall/lab_mask_exponent_se.txt prints, to the
      four decimals printed.
  F2  The premise holds: no band reaches |z| = 2 at depth 0 or at
      depth 1, and at most two do at depth 2, while every band does
      at depths 3, 4 and 5.
  F3  **And the conclusion does not.**  Under a null in which dm is
      drawn as N(0, se^2) at each band and fitted the same way, the
      observed exponent at depth 0 falls outside the null's 95 per
      cent range.
  F4  The same at depth 2.
  F5  And not at depth 1: there the observed exponent falls inside
      the null's range, so v1's "not measurable" stands where it was
      said.

REFUTATION RULE (fixed before the run)

  F1  REFUTED by any exponent outside the printed decimals.  THIS
      ONE GATES.
  F2  REFUTED if the counts differ.  Then the premise being argued
      against is not the one the data supports and F3 to F5 are
      about nothing.  THIS ONE GATES.
  F3  REFUTED if depth 0's exponent falls inside the null.  That is
      the outcome under which OPEN.md's sentence stands as written:
      the exponent would be what fitting the magnitude of noise
      produces, and quoting it at fifty standard errors would be
      quoting the error's own decay.  **Inside the null is not the
      same as consistent with zero** -- the null is centred near the
      se exponent, not at zero, and a result inside it says the
      amplitude adds nothing to what the error already explains.
  F4  REFUTED at depth 2 with the same reading.
  F5  REFUTED if depth 1's exponent falls outside the null.  Then
      the depth v1 called unmeasurable is measurable after all by
      this test, and the disagreement between the two tests would
      have to be resolved before either is used.

  F1 and F2 gate.  F3 to F5 are the measurement and do not gate.

  THE NULL IS RUN, and it is F3 to F5.  At each depth dm is redrawn
  as N(0, se^2) band by band with the recorded se, refitted by the
  same weighted least squares, and the exponent recorded; the range
  quoted is the 2.5 and 97.5 percentiles of those exponents.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "results")
V1 = os.path.abspath(os.path.join(ROOT, "..", "v1", "results", "wall"))
OUT = os.path.join(RES, "audit_mask_floornull.txt")

SEED = 20260823
DRAWS = 4000
DEC = 4
FLOOR = 2.0

ROW = (r"^  depth (\d)\s+\([^)]*\)\s+(\d+) points\s*\n[^\n]*N\s+n_c"
       r"[^\n]*\n((?:^\s+[\d.eE+-]+\s+\d+\s+[-+]?[\d.]+\s+[\d.]+\s+"
       r"[-+]?[\d.]+\s*\n)+)")


def read_bands():
    src = io.open(os.path.join(V1, "lab_mask_amplitude_law.txt"),
                  encoding="utf-8").read()
    out = {}
    for d, _n, body in re.findall(ROW, src, re.M):
        out[int(d)] = np.array(
            [[float(t) for t in ln.split()]
             for ln in body.strip().splitlines()])
    return out


def read_published():
    src = io.open(os.path.join(V1, "lab_mask_exponent_se.txt"),
                  encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^\s+(\d)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+"
                         r"([\d.]+)\s+([\d.]+)\s*$", src, re.M):
        out[int(m.group(1))] = float(m.group(3))
    return out


def read_published_se():
    src = io.open(os.path.join(V1, "lab_mask_exponent_se.txt"),
                  encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^\s+(\d)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+"
                         r"([\d.]+)\s+([\d.]+)\s*$", src, re.M):
        out[int(m.group(1))] = float(m.group(4))
    return out


def wls(x, y, w):
    A = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    inv = np.linalg.inv(A.T.dot(W).dot(A))
    c = inv.dot(A.T).dot(W).dot(y)
    return -float(c[1])


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    bands = read_bands()
    pub = read_published()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read %d depth tables from "
        "v1/results/wall/lab_mask_amplitude_law.txt; nothing is "
        "measured here" % len(bands))
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say("PRINTBOUND audit_mask_floornull %d %.8f" % (DEC, rnd))

    # -------------------------------------------------------------- F1
    say()
    say("F1  the control")
    f1 = True
    fits = {}
    for d in sorted(bands, reverse=True):
        R = bands[d]
        keep = R[:, 2] != 0.0
        x = np.log(R[keep, 0])
        dm, se = R[keep, 2], R[keep, 3]
        w = (dm / se) ** 2
        a = wls(x, np.log(np.abs(dm)), w)
        fits[d] = (x, dm, se, w, a)
        ok = abs(a - pub[d]) <= rnd
        f1 = f1 and ok
        say("  depth %d  a = %.4f against v1's %.4f  %s"
            % (d, a, pub[d], "yes" if ok else "NO"))
    say("  F1 %s   (cap: the four decimals v1 prints)"
        % ("hold" if f1 else "REFUTED"))

    # -------------------------------------------------------------- F2
    say()
    say("F2  the premise: where does the amplitude clear its floor?")
    cnt = {}
    for d in sorted(bands):
        z = np.abs(bands[d][:, 4])
        cnt[d] = int((z >= FLOOR).sum())
        say("  depth %d  |z| runs %.2f to %.2f, %d of %d bands reach "
            "%.1f" % (d, z.min(), z.max(), cnt[d], z.size, FLOOR))
    f2 = (cnt[0] == 0 and cnt[1] == 0 and cnt[2] <= 2
          and all(cnt[d] == bands[d].shape[0] for d in (3, 4, 5)))
    say("  F2 %s   (cap: none at 0 and 1, at most two at 2, all at "
        "3, 4, 5)" % ("hold" if f2 else "REFUTED"))
    d0 = bands[0]
    say("  and yet at depth 0 the amplitude falls smoothly: %s"
        % ", ".join("%.4f" % v for v in np.abs(d0[:5, 2])))
    say("  while its error barely moves: %s"
        % ", ".join("%.4f" % v for v in d0[:5, 3]))
    if not (f1 and f2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # --------------------------------------------------------- F3..F5
    say()
    say("F3/F4/F5  the null: what does fitting the magnitude of noise "
        "give?")
    say("  depth  observed  null 2.5%   null 50%   null 97.5%  "
        "outside?")
    rng = np.random.default_rng(SEED)
    verdict = {}
    spread = {}
    pubse = read_published_se()
    for d in (0, 1, 2):
        x, dm, se, w, a = fits[d]
        vals = []
        for _ in range(DRAWS):
            nd = rng.normal(0.0, se)
            nz = nd != 0.0
            if nz.sum() < 3:
                continue
            vals.append(wls(x[nz], np.log(np.abs(nd[nz])),
                            (nd[nz] / se[nz]) ** 2))
        v = np.array(vals)
        lo = float(np.percentile(v, 2.5))
        md = float(np.percentile(v, 50.0))
        hi = float(np.percentile(v, 97.5))
        out = a < lo or a > hi
        verdict[d] = out
        spread[d] = (lo, hi)
        say("  %-6d %-9.4f %-11.4f %-10.4f %-11.4f %s"
            % (d, a, lo, md, hi, "yes" if out else "no"))
        say("NULL mask_floornull_depth%d %.4f"
            % (d, float((v >= a).mean())))
    f3, f5, f4 = verdict[0], not verdict[1], verdict[2]
    say("  F3 %s   (cap: depth 0 outside the null)"
        % ("hold" if f3 else "REFUTED"))
    say("  F4 %s   (cap: depth 2 outside the null)"
        % ("hold" if f4 else "REFUTED"))
    say("  F5 %s   (cap: depth 1 inside the null)"
        % ("hold" if f5 else "REFUTED"))

    say()
    say("what the null implies about the published errors, derived "
        "and predicted by nothing")
    say("  the fit's own error assumes its weights are the true "
        "inverse variances of")
    say("  log|dm|, which breaks where |dm| is comparable to se -- "
        "the definition of")
    say("  these depths. The null's own spread is the honest one.")
    say("  depth  fit se    null-implied se   fit t     null t")
    for d in (0, 1, 2):
        x, dm, se, w, a = fits[d]
        lo, hi = spread[d]
        half = (hi - lo) / 2.0 / 1.96
        say("  %-6d %-9.4f %-17.4f %-9.1f %.1f"
            % (d, pubse[d], half, a / pubse[d], a / half))
        say("NULLSE mask_floornull_depth%d %.4f" % (d, half))

    say()
    say("what this settles")
    if f3 and f4:
        say("  the amplitude at depths 0 and 2 never clears its floor "
            "band by band and its")
        say("  decay is still outside what fitting the magnitude of "
            "noise produces, so")
        say("  \"nothing can be fitted there\" does not follow from "
            "\"no band clears\"")
    else:
        say("  at least one shallow depth's exponent is what noise "
            "would give, so the")
        say("  sentence stands where that is true")

    say()
    say("=" * 70)
    say("F1 %s  F2 %s  F3 %s  F4 %s  F5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (f1, f2, f3, f4, f5)))

    head = [
        "STATISTIC: per depth, how many of the fifteen half-octave",
        "           bands have |dm|/se at or above two; and the",
        "           fitted decay exponent of |dm| against the",
        "           distribution of that exponent when dm is redrawn",
        "           as N(0, se^2) band by band and fitted the same",
        "           way.",
        "NULL: RUN, and it is F3 to F5. At each shallow depth dm is",
        "      redrawn as N(0, se^2) with the recorded se, refitted",
        "      by the same weighted least squares, and the exponent",
        "      recorded; the range quoted is the 2.5 and 97.5",
        "      percentiles over 4000 draws with the fixed SEED. The",
        "      null is centred near the exponent of se itself, not",
        "      at zero, because the magnitude of noise tracks its",
        "      own scale.",
        "FIELD: the per-depth band tables of",
        "       v1/results/wall/lab_mask_amplitude_law.txt -- fifteen",
        "       half-octave bands from N ~ 1.189e5 to 1.431e7, six",
        "       depths by which of 3, 5, 7, 11, 13 divide N, with the",
        "       cell mean minus band mean and its exact error from",
        "       Lemma lem:cellmom. Nothing is measured here; the fit",
        "       is the weighted least squares identified in",
        "       {#rem:maskdmse} and checked again by F1 against",
        "       v1/results/wall/lab_mask_exponent_se.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (f1 and f2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
