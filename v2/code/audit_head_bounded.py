# -*- coding: utf-8 -*-
r"""
Which m carry the head's sign -- a bounded set, or all of them?

WHAT IS AT STAKE

Item 4(b) is now a demand about the head: {#rem:splitreach} showed the
deficit sits there and nowhere else, {#rem:headaxis} showed the head is
selected and aligned by the imbalance I = H/T rather than by k, and
{#rem:destination} showed a coin is inside what the data allow for its
alignment but says nothing about the mechanism. Nothing yet asks the
one structural question the inner sum permits: H(N;k) is a sum over m,
and the m are ordered. WHICH m make it negative?

The question has a consequence rather than only an answer.
{#rem:provablehalf} calls a condition elementary when it is
multiplicative or of BOUNDED modulus, and {#rem:sievedepth} found that
nothing of bounded modulus carries the slope. But that was about which
k are negative. If the head's sign is carried by a bounded number of
LEADING terms -- m <= M0 for fixed M0, a signed sum of at most M0
values of Lambda -- then the obstruction is an object of exactly the
kind the elementary half already handles, and item 4(b) changes
character. If instead the sign survives deleting them, the alignment is
a property of the whole inner sum and no truncation reaches it.

BACKS: Remark {#rem:headbounded} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  The control. The split is exact -- A + B reproduces H to 1e-9
      relative at every k -- and at N = 25600000 the head's one-sign
      fraction reproduces results/audit_split_reach.txt.
  R2  The bounded part carries the sign: on the head, sign A agrees
      with sign H more often than an arm matched to sign A's own
      marginal rate, resolved at two standard deviations, at every N
      and for both cutoffs.
  R3  And it is what the alignment is made of: with the bounded part
      deleted, the one-sign fraction of B on the head falls to within
      two binomial deviations of one half at the top N.
  R4  The bounded part's share of |H| on the head decays: its
      least-squares exponent in log N is resolved negative, since the
      number of terms it omits grows like N/k.

REFUTATION RULE (fixed before the run)

  R1  REFUTED on any relative departure above 1e-9, or if the head's
      fraction misses the published by more than the bound that
      table's printing forces, computed here. THIS ONE GATES.
  R2  REFUTED if the excess over the matched arm fails to resolve at
      any N or either cutoff. Then the leading terms do not carry the
      head's sign and the alignment is not a truncation effect.
  R3  REFUTED if the remainder still clears one half. Then the head's
      alignment survives deleting a bounded number of leading terms,
      the obstruction is not a bounded object however M0 is chosen at
      this scale, and {#rem:sievedepth}'s verdict for the slope holds
      for the head's sign too.
  R4  REFUTED if the share's exponent is not resolved negative. Then
      the bounded part keeps a fixed proportion of the magnitude as N
      grows, which would be a stronger statement than R2 and would
      say the head is asymptotically a bounded object.

  R1 gates. R2 to R4 are the measurement and do not gate.

  THE NULL IS RUN AND IT IS THE POINT. A one-sign head makes any
  predictor look good, as G73 was added to enforce: the head is
  four-fifths one sign, so a constant predictor agrees four times in
  five. Every agreement below is therefore reported against an arm
  drawn at sign A's own marginal rate on the same set, and the
  MARGINAL line declares that rate. The arm uses the fixed SEED
  declared in the output.
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
OUT = os.path.join(RES, "audit_head_bounded.txt")

LO, HI = 200_000, 102_400_000
CUTS = (29, 1000)
SEED = 20260820
DRAWS = 400


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA
HEAD = SPL.HEAD


def read_published():
    """the head's one-sign fraction, and the decimals it is printed to"""
    src = io.open(os.path.join(RES, "audit_split_reach.txt"),
                  encoding="utf-8").read()
    ag, dec = {}, 0
    for m in re.finditer(r"^  (\d{5,})\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+"
                         r"\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s*$",
                         src, re.M):
        ag[int(m.group(1))] = float(m.group(2))
        dec = max(dec, len(m.group(2).split(".")[1]))
    return ag, dec


def family(lo, hi):
    out = []
    a = 1
    while 2 ** a <= hi:
        b = 1
        while 2 ** a * 5 ** b <= hi:
            v = 2 ** a * 5 ** b
            if v >= lo:
                out.append(v)
            b += 1
        a += 1
    return sorted(set(out))


def split_sums(N, lam, mu, sqf, cuts):
    """H(N;k) cut at each m <= M0, exactly audit_gain_split's inner sum"""
    PN = SPL.factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    H = np.zeros(ks.size)
    A = {c: np.zeros(ks.size) for c in cuts}
    B = {c: np.zeros(ks.size) for c in cuts}
    for i, k in enumerate(ks):
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in SPL.factor_set(k):
            ms = ms[ms % q != 0]
        vals = lam[N - ms * k] * mu[ms].astype(np.float64)
        H[i] = float(vals.sum())
        for c in cuts:
            # both halves summed from the mask, not one from the
            # other, so R1 checks the masking and not an identity
            A[c][i] = float(vals[ms <= c].sum())
            B[c][i] = float(vals[ms > c].sum())
    return ks, H, A, B


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if x.size > 2 else float("inf")
    return float(a), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubag, dec = read_published()
    NS = family(LO, HI)
    CTRL = 25_600_000
    rng = np.random.default_rng(SEED)

    say("read %d head one-sign fractions from "
        "results/audit_split_reach.txt" % len(pubag))
    say("  the field, the sieve, theta' and the head fraction are "
        "imported from code/audit_gain_split.py")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say()
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]; %d of them" % (LO, HI, len(NS)))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))
    say("  the cutoffs are m <= %s, fixed and not growing with N"
        % " and m <= ".join(str(c) for c in CUTS))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    rows = []
    worstrel = 0.0
    say()
    say("  N            #k      head   onesign  " +
        "  ".join("M0=%-5d agree  excess  share  rest" % c
                  for c in CUTS))
    for N in NS:
        ks, H, A, Bs = split_sums(N, lam, mu, sqf, CUTS)
        a = np.log(ks.astype(np.float64)) * H
        w = np.abs(a)
        nh = max(1, int(round(HEAD * ks.size)))
        hd = np.argsort(-w)[:nh]
        sg = np.sign(a[hd])
        onesign = max(float((sg > 0).mean()), float((sg < 0).mean()))
        got = [N, int(ks.size), nh, onesign]
        for c in CUTS:
            B = Bs[c]
            worstrel = max(worstrel,
                           float(np.abs((A[c] + B) - H).max()
                                 / max(float(np.abs(H).max()), 1e-12)))
            sa = np.sign(A[c][hd])
            sh = np.sign(H[hd])
            agree = float((sa == sh).mean())
            marg = max(float((sa > 0).mean()), float((sa < 0).mean()))
            # the matched arm: same marginal rate, drawn independently
            draws = []
            for _ in range(DRAWS):
                fake = np.where(rng.random(nh) < marg, 1.0, -1.0)
                if float((sa > 0).mean()) < 0.5:
                    fake = -fake
                draws.append(float((fake == sh).mean()))
            mu_a = float(np.mean(draws))
            sd_a = float(np.std(draws)) or 1e-12
            excess = (agree - mu_a) / sd_a
            share = float(np.abs(A[c][hd]).sum()
                          / np.abs(H[hd]).sum())
            sb = np.sign(B[hd])
            rest = max(float((sb > 0).mean()), float((sb < 0).mean()))
            got += [agree, excess, share, rest, marg]
        rows.append(tuple(got))
        say("  %-12d %-7d %-6d %-8.4f " % tuple(got[:4]) +
            "  ".join("%-8.4f %-7.2f %-6.4f %.4f"
                      % (got[4 + 5 * j], got[5 + 5 * j],
                         got[6 + 5 * j], got[7 + 5 * j])
                      for j in range(len(CUTS))))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    top = rows[-1]

    # -------------------------------------------------------------- R1
    say()
    say("R1  the control")
    ctrl = [r for r in rows if r[0] == CTRL][0]
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(ctrl[3] - pubag.get(CTRL, float("nan")))
    r1 = worstrel <= 1e-9 and d <= rnd
    say("  the split is exact to %.3e relative" % worstrel)
    say("  the head's one-sign fraction here %.4f against the "
        "published %.4f, departure %.6f"
        % (ctrl[3], pubag.get(CTRL, float("nan")), d))
    say("  the table prints %d decimals, so the bound is %.8f"
        % (dec, rnd))
    say("PRINTBOUND audit_head_bounded %d %.8f" % (dec, rnd))
    say("  R1 %s   (tol 1e-9 relative and the printing bound)"
        % ("hold" if r1 else "REFUTED"))

    # -------------------------------------------------------------- R2
    say()
    say("R2  does the bounded part carry the head's sign?")
    r2 = True
    marginals = []
    for j, c in enumerate(CUTS):
        ex = np.array([r[5 + 5 * j] for r in rows])
        ag = np.array([r[4 + 5 * j] for r in rows])
        mg = np.array([r[8 + 5 * j] for r in rows])
        okall = bool((ex >= 2.0).all())
        if not okall:
            r2 = False
        say("  M0 = %-6d agreement runs %.4f to %.4f, the matched "
            "arm's excess %.2f to %.2f, resolved everywhere: %s"
            % (c, float(ag.min()), float(ag.max()),
               float(ex.min()), float(ex.max()),
               "yes" if okall else "no"))
        say("  the predictor's own majority share runs %.4f to %.4f"
            % (float(mg.min()), float(mg.max())))
        marginals.append(float(mg.max()))
    say("MARGINAL audit_head_bounded %.4f" % max(marginals))
    if max(marginals) >= 0.9:
        say("DEGENERATE audit_head_bounded")
    say("  R2 %s   (cap 2 standard deviations of the matched arm)"
        % ("hold" if r2 else "REFUTED"))

    # -------------------------------------------------------------- R3
    say()
    say("R3  and does deleting it kill the alignment?")
    r3 = True
    for j, c in enumerate(CUTS):
        rest = top[7 + 5 * j]
        sd = math.sqrt(0.25 / top[2])
        cl = (rest - 0.5) <= 2.0 * sd
        if not cl:
            r3 = False
        say("  M0 = %-6d the remainder's one-sign fraction at the top "
            "N is %.4f, %.2f binomial deviations above one half"
            % (c, rest, (rest - 0.5) / sd))
    say("  the head there has %d members, so one deviation is %.6f"
        % (top[2], math.sqrt(0.25 / top[2])))
    say("  R3 %s   (cap 2 binomial deviations)"
        % ("hold" if r3 else "REFUTED"))

    # -------------------------------------------------------------- R4
    say()
    say("R4  does the bounded part's share of the magnitude decay?")
    r4 = True
    for j, c in enumerate(CUTS):
        sh = np.array([r[6 + 5 * j] for r in rows])
        e, rms, se = fit(x, np.log(sh))
        res = e < 0.0 and abs(e) / se >= 2.0
        if not res:
            r4 = False
        say("  M0 = %-6d share runs %.4f to %.4f, exponent %+.6f, "
            "s.e. %.6f, t = %.2f"
            % (c, float(sh.min()), float(sh.max()), e, se,
               abs(e) / se))
        say("TSTAT slope_headbounded_share%d %.2f" % (c, abs(e) / se))
        say("SPREAD slope_headbounded_share%d %.4f"
            % (c, float(x.max() - x.min())))
        if abs(e) / se < 2.0:
            say("UNRESOLVED SIGN slope_headbounded_share%d" % c)
    say("  R4 %s   (cap 2 standard errors)"
        % ("hold" if r4 else "REFUTED"))

    say()
    say("  and the head's own one-sign fraction over the field, for "
        "reference:")
    on = np.array([r[3] for r in rows])
    e, rms, se = fit(x, on)
    say("  it runs %.4f to %.4f with least-squares slope in log N = "
        "%+.6f, s.e. %.6f" % (float(on.min()), float(on.max()), e, se))
    say("TSTAT slope_audit_head_bounded %.2f" % (abs(e) / se))
    say("SPREAD slope_audit_head_bounded %.4f"
        % float(x.max() - x.min()))
    if abs(e) / se < 2.0:
        say("UNRESOLVED SIGN slope_audit_head_bounded")

    say()
    say("=" * 70)
    say("R1 %s  R2 %s  R3 %s  R4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (r1, r2, r3, r4)))

    head = [
        "STATISTIC: for H(N;k) = sum_m Lambda(N-mk) mu(m) over the m",
        "           coprime to k below N/k, the split H = A + B with A",
        "           the terms m <= M0 for the fixed cutoffs M0 = 29",
        "           and M0 = 1000; on the top tenth of k by",
        "           |a_k| = (log k)|H|, the fraction of k where sign A",
        "           equals sign H against an arm drawn at sign A's own",
        "           marginal rate on the same set, A's share of",
        "           sum|H|, and the majority-sign fraction of the",
        "           remainder B; at every on-field N to 1.024e8, with",
        "           the share's least-squares exponent in log N.",
        "NULL: run, and it is the point. A head that is four-fifths",
        "      one sign makes any predictor agree four times in five,",
        "      which is what G73 exists to stop being read as skill.",
        "      Every agreement is reported as an excess over an arm",
        "      drawn at the predictor's own marginal rate on the same",
        "      set, with the MARGINAL line declaring that rate, and",
        "      the arm uses the fixed SEED above.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the field, the sieve,",
        "       theta' and the head fraction are imported from",
        "       code/audit_gain_split.py; the published one-sign",
        "       fractions are read from results/audit_split_reach.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not r1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
