# -*- coding: utf-8 -*-
r"""
The mask's decay, fitted on dm/se as v1 said to and nobody did.

WHAT IS AT STAKE

v1 handed one instruction forward with its open item on the mask's
decay exponent, and OPEN.md item 3 has carried it unexecuted:
**fit dm/se, not dm.**

The reason is in the data.  v1 fitted |dm| ~ N^-a per depth, where
dm = m_c - gm is the cell mean minus the band mean and se is its
exact error from Lemma lem:cellmom.  But se falls with N too, because
the cell population n_c grows with N -- at depth 5 it runs 2, 3, 4,
5, 8, 11, 16, 22 across the bands while se falls 0.7136 to 0.4727.
Over the same bands |dm| falls 7.0004 to 4.9979.  **So most of what
was fitted as the amplitude decaying is the error shrinking**, and
|z| = |dm|/se, the quantity that says whether the mask is there at
all, actually grows from 9.81 to 10.57.

The decomposition is exact: with the same weights, the exponent of
|dm/se| is the exponent of |dm| minus the exponent of se.  So this
run does not replace v1's number, it splits it, and the question is
how much of each depth's decay survives when the error's own decay is
taken out.

Nothing is measured.  Every value is read from the per-depth tables
of v1/results/wall/lab_mask_amplitude_law.txt, and the fit v1 used is
identified rather than assumed: weighted least squares of log|dm| on
log N with weights (dm/se)^2 and the covariance scaled by the
residual variance reproduces all six published exponents and five of
the six published standard errors, which E1 checks.

BACKS: Remark {#rem:maskdmse} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  The control.  Refitting |dm| that way returns the six
      exponents v1/results/wall/lab_mask_exponent_se.txt prints, to
      the four decimals printed, and their standard errors at the
      five depths where no dm vanishes.
  E2  The error decays too: at every depth the exponent of se
      against N is positive and resolved, |t| > 2.
  E3  So the amplitude exponent falls when the error's is taken out:
      at every depth the |dm/se| exponent is smaller than the |dm|
      one.
  E4  **And at least one depth loses its decay entirely**: a depth
      whose |dm| exponent v1 quotes as resolved has a |dm/se|
      exponent that is not, |t| < 2.
  E5  The monotone rise v1 reported -- the exponent rising as the
      cell gets shallower, excluding depth 1 -- does not survive in
      |dm/se|.

REFUTATION RULE (fixed before the run)

  E1  REFUTED by any exponent outside the printed decimals.  Then
      the fit identified here is not v1's and nothing below is about
      its numbers.  THIS ONE GATES.
  E2  REFUTED if any depth's se exponent is negative or unresolved.
      Unresolved would mean the error does not measurably decay
      there, so nothing is being subtracted and E3 is empty at that
      depth -- which is a fact about the design, not about the mask.
  E3  REFUTED at any depth where the |dm/se| exponent is not
      smaller.  That can only happen if se grows with N there, which
      E2 would already have caught.
  E4  REFUTED if every depth keeps a resolved decay.  Then the
      error's decay was a part of what v1 fitted but not the whole,
      each depth's mask really does decay, and the instruction
      changes the numbers without changing the reading.
  E5  REFUTED if the order survives.  The monotonicity would then be
      a property of the mask rather than of how the cell populations
      grow with depth.  **If instead it does not survive, the
      published reading "the mask decays faster where fewer small
      primes divide N" describes the cell counts and not the mask.**

  E1 gates.  E2 to E5 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  Nothing is sampled: two fits are
  compared on the same recorded table and the question is which
  exponent belongs to which quantity.  The coin arms for the mask
  were run in v1's lab_mask_coin_control.py.
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
OUT = os.path.join(RES, "audit_mask_dmse.txt")

DEC = 4

ROW = (r"^  depth (\d)\s+\([^)]*\)\s+(\d+) points\s*\n[^\n]*N\s+n_c"
       r"[^\n]*\n((?:^\s+[\d.eE+-]+\s+\d+\s+[-+]?[\d.]+\s+[\d.]+\s+"
       r"[-+]?[\d.]+\s*\n)+)")


def read_bands():
    src = io.open(os.path.join(V1, "lab_mask_amplitude_law.txt"),
                  encoding="utf-8").read()
    out = {}
    for d, npts, body in re.findall(ROW, src, re.M):
        rows = np.array([[float(t) for t in ln.split()]
                         for ln in body.strip().splitlines()])
        out[int(d)] = (int(npts), rows)
    return out


def read_published():
    src = io.open(os.path.join(V1, "lab_mask_exponent_se.txt"),
                  encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^\s+(\d)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+"
                         r"([\d.]+)\s+([\d.]+)\s*$", src, re.M):
        out[int(m.group(1))] = (float(m.group(3)), float(m.group(4)))
    return out


def wls(x, y, w):
    """v1's fit: weights (dm/se)^2, covariance scaled by the residual"""
    A = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    inv = np.linalg.inv(A.T.dot(W).dot(A))
    c = inv.dot(A.T).dot(W).dot(y)
    r = y - A.dot(c)
    s2 = float((w * r * r).sum()) / (x.size - 2)
    return -float(c[1]), math.sqrt(s2 * inv[1, 1])


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
        "v1/results/wall/lab_mask_amplitude_law.txt and the published"
        % len(bands))
    say("  exponents from lab_mask_exponent_se.txt; nothing is "
        "measured here")
    say("PRINTBOUND audit_mask_dmse %d %.8f" % (DEC, rnd))
    d5 = bands[5][1]
    say()
    say("the reason the instruction exists, echoed from the table it "
        "was left in:")
    say("  at depth 5 the cell population runs %s"
        % ", ".join("%d" % v for v in d5[:8, 1].astype(int)))
    say("  while se falls %.4f to %.4f and |dm| falls %.4f to %.4f"
        % (d5[0, 3], d5[7, 3], abs(d5[0, 2]), abs(d5[7, 2])))
    say("  so |z| = |dm|/se goes %.2f to %.2f -- it grows"
        % (abs(d5[0, 4]), abs(d5[7, 4])))

    # -------------------------------------------------------------- E1
    say()
    say("E1  the control: v1's fit, identified and reproduced")
    say("  depth  pts  a here    a pub     se here   se pub    ok")
    e1 = True
    fits = {}
    for d in sorted(bands, reverse=True):
        npts, R = bands[d]
        N, dm, se = R[:, 0], R[:, 2], R[:, 3]
        keep = dm != 0.0
        x = np.log(N[keep])
        w = (dm[keep] / se[keep]) ** 2
        a, sa = wls(x, np.log(np.abs(dm[keep])), w)
        pa, ps = pub[d]
        okA = abs(a - pa) <= rnd
        okS = abs(sa - ps) <= rnd
        e1 = e1 and okA
        fits[d] = (x, dm[keep], se[keep], w, a, sa)
        say("  %-6d %-4d %-9.4f %-9.4f %-9.4f %-9.4f %s"
            % (d, int(keep.sum()), a, pa, sa, ps,
               "yes" if okA and okS else
               ("a only" if okA else "NO")))
    say("  E1 %s   (cap: the four decimals v1 prints)"
        % ("hold" if e1 else "REFUTED"))
    if not e1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------- E2, E3, E4, E5
    say()
    say("E2/E3/E4  the error's own decay, and what is left of the "
        "mask's")
    say("  depth  a(|dm|)   s(se)     t(s)    a'(|dm/se|)  se(a')   "
        "t(a')  resolved?")
    e2 = e3 = True
    e4 = False
    rows = []
    for d in sorted(bands, reverse=True):
        x, dm, se, w, a, sa = fits[d]
        s, ss = wls(x, np.log(se), w)
        ts = s / ss if ss else float("inf")
        ap, sap = wls(x, np.log(np.abs(dm) / se), w)
        tap = ap / sap if sap else float("inf")
        e2 = e2 and (s > 0.0 and abs(ts) > 2.0)
        e3 = e3 and (ap < a)
        if abs(a / sa) > 2.0 and abs(tap) <= 2.0:
            e4 = True
        rows.append((d, a, s, ts, ap, sap, tap))
        say("  %-6d %+-9.4f %+-9.4f %-7.2f %+-12.4f %-8.4f %-6.2f %s"
            % (d, a, s, ts, ap, sap, tap,
               "yes" if abs(tap) > 2.0 else "no"))
        say("TSTAT mask_dmse_depth%d %.2f" % (d, tap))
        if abs(tap) < 2.0:
            say("UNRESOLVED SIGN mask_dmse_depth%d" % d)
        say("SPREAD mask_dmse_depth%d %.4f" % (d, x.max() - x.min()))
    say("  E2 %s   (cap: positive and |t| > 2 at every depth)"
        % ("hold" if e2 else "REFUTED"))
    say("  E3 %s   (cap: smaller at every depth)"
        % ("hold" if e3 else "REFUTED"))
    say("  E4 %s   (cap: at least one depth loses its resolution)"
        % ("hold" if e4 else "REFUTED"))

    # -------------------------------------------------------------- E5
    say()
    say("E5  does the monotone rise survive?")
    keep = [r for r in rows if r[0] != 1]
    # v1 reports the exponent rising as the cell gets shallower;
    # listed deepest first that is an increase, not a decrease
    mono_a = all(keep[i][1] < keep[i + 1][1]
                 for i in range(len(keep) - 1))
    mono_ap = all(keep[i][4] < keep[i + 1][4]
                  for i in range(len(keep) - 1))
    e5 = mono_a and not mono_ap
    say("  excluding depth 1, deepest first, the |dm| exponents run "
        "%s" % ", ".join("%+.4f" % r[1] for r in keep))
    say("  and they are monotone as v1 reports: %s"
        % ("yes" if mono_a else "NO"))
    say("  the |dm/se| exponents run %s"
        % ", ".join("%+.4f" % r[4] for r in keep))
    say("  monotone: %s" % ("yes" if mono_ap else "no"))
    say("SIGNRUN mask_dmse_monotone %d"
        % sum(1 for i in range(len(keep) - 1)
              if keep[i][4] < keep[i + 1][4]))
    say("  E5 %s   (cap: the order v1 reported)"
        % ("hold" if e5 else "REFUTED"))

    say()
    say("what this settles")
    if e4 and e5:
        say("  the error's decay carried a part of what v1 fitted, "
            "at least one depth")
        say("  keeps no resolved decay once it is out, and the "
            "monotone reading describes")
        say("  how the cell populations grow with depth rather than "
            "the mask")
    elif not e4:
        say("  every depth keeps a resolved decay, so the error's "
            "decay was a part of")
        say("  v1's number and not the whole of it; the instruction "
            "changes the numbers")
        say("  without changing the reading")
    else:
        say("  the order survives, so the monotonicity is the mask's "
            "and not the design's")

    say()
    say("=" * 70)
    say("E1 %s  E2 %s  E3 %s  E4 %s  E5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (e1, e2, e3, e4, e5)))

    head = [
        "STATISTIC: the decay exponent of the mask amplitude per",
        "           depth, fitted on |dm| as v1 fitted it and on",
        "           |dm/se| as v1's handed-forward instruction says",
        "           to, with the exponent of se itself so the",
        "           difference is accounted for; and whether the",
        "           monotone rise v1 reported across depths survives",
        "           the second fit.",
        "NULL: none is run and none applies. Nothing is sampled: two",
        "      fits are compared on the same recorded table and the",
        "      question is which exponent belongs to which quantity.",
        "      The coin arms for the mask were run in v1's",
        "      lab_mask_coin_control.py.",
        "FIELD: the per-depth band tables of",
        "       v1/results/wall/lab_mask_amplitude_law.txt -- fifteen",
        "       half-octave bands from N ~ 1.189e5 to 1.431e7, six",
        "       depths by which of 3, 5, 7, 11, 13 divide N, with the",
        "       cell mean minus band mean and its exact error from",
        "       Lemma lem:cellmom. Nothing is measured here; the fit",
        "       is weighted least squares of the log on log N with",
        "       weights (dm/se)^2 and the covariance scaled by the",
        "       residual variance, which E1 identifies against the",
        "       six exponents of",
        "       v1/results/wall/lab_mask_exponent_se.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not e1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
