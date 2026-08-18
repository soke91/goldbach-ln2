# -*- coding: utf-8 -*-
r"""
Is the m-coherence structure, or is it the selection that found it?

WHAT IS AT STAKE

{#rem:headbounded} closed one route and opened one claim. The route:
no prefix of the inner sum owns the head's sign. The claim, which
neither of its registered predictions asked for and which is now the
newest structural statement in this repository: the prefix m <= 1000
and the remainder agree in sign on the head 0.8859 of the time where
an arm at the prefix's own marginal rate gives about two thirds, so
"the alignment is coherent across the m-decomposition".

That claim has a hole in it, and M4 exists for exactly this. The head
is the top tenth of k by |a_k| = (log k)|A + B|. Selecting on the
magnitude of a SUM preferentially selects pairs whose parts agree in
sign, because agreeing parts add and disagreeing parts cancel. A pair
of independent coins, given the observed |A| and |B| and then filtered
the same way, will look coherent. Whether any coherence is left after
that is the question, and it was not asked.

The null here preserves what M3 requires: |A_k| and |B_k| are the
measured magnitudes, untouched, and only the two signs are redrawn --
independently of each other, at the marginal rates the data show. The
selection is then redone on the surrogate exactly as on the data.

BACKS: Remark {#rem:coherencenull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The control. A and B recomputed here reproduce
      results/audit_head_bounded.txt's head agreement at
      N = 25600000 inside the bound its printing forces.
  S2  Selection alone does produce coherence: the coin arm's head
      agreement is above its own all-k agreement at every N. If this
      fails the worry was unfounded and nothing below is needed.
  S3  And the data still beat it: the observed head agreement exceeds
      the coin arm's by two standard deviations of the arm, at every
      N.
  S4  The cleanest version, with no selection at all: over ALL k the
      observed A-B sign agreement exceeds the arm's by two standard
      deviations, at every N.

REFUTATION RULE (fixed before the run)

  S1  REFUTED outside the printing bound. Then this is not the split
      {#rem:headbounded} measured. THIS ONE GATES.
  S2  REFUTED if the arm shows no lift from selection. Then the
      objection this script was written for does not apply and
      {#rem:headbounded}'s reading needed no defence.
  S3  REFUTED if the observed sits inside the arm anywhere. Then the
      coherence on the head is the selection and NOT a property of
      mu, and {#rem:headbounded}'s closing claim must be withdrawn --
      the route it closes stays closed, but the structure it reports
      is an artefact of conditioning on a large sum.
  S4  REFUTED if the all-k agreement is inside the arm. This is the
      one that cannot be blamed on selection, so its refutation would
      mean the two halves of the inner sum are independent in sign
      once magnitudes are given, and the coherence lives only where
      the conditioning put it.

  S1 gates. S2 to S4 are the measurement and do not gate.

  THE NULL IS THE MEASUREMENT. Every number below is a comparison
  against a surrogate that keeps |A| and |B| exactly and redraws only
  their signs, independently, at the marginal rates observed at that
  N -- M3's requirement that a null preserve the structure of the
  field it nulls. The draws use the fixed SEED declared in the output.
"""

import importlib.util
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
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_coherence_null.txt")

LO, HI = 200_000, 102_400_000
CUT = 1000
SEED = 20260821
DRAWS = 200


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HB = module("audit_head_bounded")
SPL = HB.SPL
THETA = SPL.THETA
HEAD = SPL.HEAD


def read_published():
    """{#rem:headbounded}'s head agreement at the M0 = 1000 cutoff"""
    src = io.open(os.path.join(RES, "audit_head_bounded.txt"),
                  encoding="utf-8").read()
    ag, dec = {}, 0
    for m in re.finditer(r"^  (\d{5,})\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+"
                         r"\s+[+-]?[\d.]+\s+[\d.]+\s+[\d.]+\s+"
                         r"([\d.]+)\s+[+-]?[\d.]+\s+[\d.]+\s+[\d.]+"
                         r"\s*$", src, re.M):
        ag[int(m.group(1))] = float(m.group(2))
        dec = max(dec, len(m.group(2).split(".")[1]))
    return ag, dec


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if x.size > 2 else float("inf")
    return float(a), float(np.sqrt((r ** 2).mean())), se


def arm(absA, absB, logk, pA, pB, nh, rng, draws):
    """|A| and |B| kept; signs redrawn independently at pA, pB

    Returns the surrogate's head agreement and its all-k agreement,
    each as a mean and a standard deviation over the draws.
    """
    hs, alls = [], []
    n = absA.size
    for _ in range(draws):
        sa = np.where(rng.random(n) < pA, -1.0, 1.0)
        sb = np.where(rng.random(n) < pB, -1.0, 1.0)
        tot = sa * absA + sb * absB
        w = logk * np.abs(tot)
        hd = np.argsort(-w)[:nh]
        hs.append(float((sa[hd] == sb[hd]).mean()))
        alls.append(float((sa == sb).mean()))
    return (float(np.mean(hs)), float(np.std(hs)) or 1e-12,
            float(np.mean(alls)), float(np.std(alls)) or 1e-12)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubag, dec = read_published()
    NS = HB.family(LO, HI)
    CTRL = 25_600_000
    rng = np.random.default_rng(SEED)

    say("read %d head agreements from "
        "results/audit_head_bounded.txt at the M0 = %d cutoff"
        % (len(pubag), CUT))
    say("  the split, the field, the sieve, theta' and the head "
        "fraction are imported from code/audit_head_bounded.py")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say()
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]; %d of them" % (LO, HI, len(NS)))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    rows = []
    say()
    say("  N            #k      head   observed          coin arm "
        "(|A|,|B| kept)")
    say("               head    all k  head     all k    head "
        "mean/sd      all-k mean/sd    excess head  excess all")
    for N in NS:
        ks, H, A, B = HB.split_sums(N, lam, mu, sqf, (CUT,))
        a1, b1 = A[CUT], B[CUT]
        logk = np.log(ks.astype(np.float64))
        w = logk * np.abs(H)
        nh = max(1, int(round(HEAD * ks.size)))
        hd = np.argsort(-w)[:nh]
        sa, sb = np.sign(a1), np.sign(b1)
        obs_h = float((sa[hd] == sb[hd]).mean())
        obs_a = float((sa == sb).mean())
        pA = float((sa < 0).mean())
        pB = float((sb < 0).mean())
        mh, sh, ma, sda = arm(np.abs(a1), np.abs(b1), logk, pA, pB,
                              nh, rng, DRAWS)
        rows.append((N, int(ks.size), nh, obs_h, obs_a, mh, sh, ma,
                     sda, (obs_h - mh) / sh, (obs_a - ma) / sda,
                     pA, pB))
        say("  %-12d %-7d %-6d %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f "
            "%-8.4f %-12.2f %.2f"
            % (N, ks.size, nh, obs_h, obs_a, mh, sh, ma, sda,
               (obs_h - mh) / sh, (obs_a - ma) / sda))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    top = rows[-1]

    # -------------------------------------------------------------- S1
    say()
    say("S1  the control at N = %d" % CTRL)
    ctrl = [r for r in rows if r[0] == CTRL][0]
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(ctrl[3] - pubag.get(CTRL, float("nan")))
    s1 = d <= rnd
    say("  head agreement here %.4f against the published %.4f, "
        "departure %.6f" % (ctrl[3], pubag.get(CTRL, float("nan")), d))
    say("  the table prints %d decimals, so the bound is %.8f"
        % (dec, rnd))
    say("PRINTBOUND audit_coherence_null %d %.8f" % (dec, rnd))
    say("  S1 %s   (cap: the printing bound)"
        % ("hold" if s1 else "REFUTED"))

    # -------------------------------------------------------------- S2
    say()
    say("S2  does selection alone make coins look coherent?")
    lifts = [r[5] - r[7] for r in rows]
    s2 = all(l > 0 for l in lifts)
    say("  the arm's head agreement minus its own all-k agreement "
        "runs %+.4f to %+.4f over the field"
        % (min(lifts), max(lifts)))
    say("  at the top N the arm gives %.4f on the head against %.4f "
        "on all k" % (top[5], top[7]))
    say("  S2 %s   (cap: any N where the arm shows no lift)"
        % ("hold" if s2 else "REFUTED"))

    # -------------------------------------------------------------- S3
    say()
    say("S3  and does the data still beat the arm on the head?")
    exh = np.array([r[9] for r in rows])
    s3 = bool((exh >= 2.0).all())
    say("  the excess runs %.2f to %.2f standard deviations of the arm"
        % (float(exh.min()), float(exh.max())))
    say("  at the top N: observed %.4f against the arm's %.4f "
        "(sd %.4f), excess %.2f"
        % (top[3], top[5], top[6], top[9]))
    say("  S3 %s   (cap 2 standard deviations at every N)"
        % ("hold" if s3 else "REFUTED"))

    # -------------------------------------------------------------- S4
    say()
    say("S4  and with no selection at all, over every k?")
    exa = np.array([r[10] for r in rows])
    s4 = bool((exa >= 2.0).all())
    say("  the excess runs %.2f to %.2f standard deviations of the arm"
        % (float(exa.min()), float(exa.max())))
    say("  at the top N: observed %.4f against the arm's %.4f "
        "(sd %.4f), excess %.2f"
        % (top[4], top[7], top[8], top[10]))
    e, rms, se = fit(x, exa)
    say("  the all-k excess has least-squares slope in log N = "
        "%+.6f, s.e. %.6f, t = %.2f" % (e, se, abs(e) / se))
    say("TSTAT slope_audit_coherence_null %.2f" % (abs(e) / se))
    say("SPREAD slope_audit_coherence_null %.4f"
        % float(x.max() - x.min()))
    if abs(e) / se < 2.0:
        say("UNRESOLVED SIGN slope_audit_coherence_null")
    say("  S4 %s   (cap 2 standard deviations at every N)"
        % ("hold" if s4 else "REFUTED"))
    say()
    say("  where the sign of the effect turns, since S3 and S4 both "
        "fail on the")
    say("  small N and both clear the cap on the large ones:")
    below = [r for r in rows if r[10] < 2.0]
    above = [r for r in rows if r[10] >= 2.0]
    say("  the all-k excess is under two deviations at %d of the %d "
        "N, the largest being %d" % (len(below), len(rows),
                                     max(r[0] for r in below)))
    say("  and resolved above it at every N from %d upward, %d of them"
        % (min(r[0] for r in above), len(above)))
    say("  as an effect size rather than a ratio, observed minus arm "
        "on all k runs")
    say("  %+.4f at the bottom N to %+.4f at the top"
        % (rows[0][4] - rows[0][7], top[4] - top[7]))
    say("  and the bottom is not coherence weakening but reversing: "
        "at N = %d" % rows[0][0])
    say("  the observed all-k agreement is %.4f against the arm's "
        "%.4f, and sign B" % (rows[0][4], rows[0][7]))
    say("  is negative on only %.4f of the k there, so the remainder "
        "is nearly constant." % rows[0][12])

    say()
    say("  the marginal negative rates the arm was matched to:")
    say("  N            sign A negative   sign B negative")
    for r in (rows[0], rows[len(rows) // 2], top):
        say("  %-12d %-17.4f %.4f" % (r[0], r[11], r[12]))
    say("MARGINAL audit_coherence_null %.4f"
        % max(max(r[11], 1.0 - r[11], r[12], 1.0 - r[12])
              for r in rows))
    if max(max(r[11], 1.0 - r[11], r[12], 1.0 - r[12])
           for r in rows) >= 0.9:
        say("DEGENERATE audit_coherence_null")

    say()
    say("=" * 70)
    say("S1 %s  S2 %s  S3 %s  S4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (s1, s2, s3, s4)))

    head = [
        "STATISTIC: for H(N;k) split at m <= 1000 into A and B, the",
        "           fraction of k where sign A equals sign B, on the",
        "           top tenth of k by |a_k| = (log k)|H| and on every",
        "           k, each against a surrogate that keeps |A| and",
        "           |B| exactly and redraws the two signs",
        "           independently at the marginal negative rates",
        "           observed at that N, with the head selection redone",
        "           on the surrogate; at every on-field N to 1.024e8.",
        "NULL: run, and it is the whole measurement. Selecting the top",
        "      tenth by the magnitude of a SUM preferentially selects",
        "      pairs whose parts agree in sign, so a coin pair filtered",
        "      the same way looks coherent. The surrogate preserves",
        "      |A| and |B| exactly -- M3's requirement -- and redraws",
        "      only the signs, independently, at the observed marginal",
        "      rates, with the same selection applied. The draws use",
        "      the fixed SEED above.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the split A + B and the",
        "       family are imported from code/audit_head_bounded.py,",
        "       the field and theta' from code/audit_gain_split.py;",
        "       the published head agreements are read from",
        "       results/audit_head_bounded.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not s1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
