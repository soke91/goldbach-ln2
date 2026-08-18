# -*- coding: utf-8 -*-
r"""
What the surviving correlation is supported on.

WHAT IS AT STAKE

{#rem:denominator} showed sum a is sum_j Lambda(N-j) Lambda_K(j) where
Lambda_K is the mu * log convolution with the k-index restricted, and
{#rem:thetalaw} showed alpha rises towards 1 as theta' does, so the
truncation is what removes the main term and returns it. Neither says
what the surviving sum SITS ON.

The index tells us, and the arithmetic is exact rather than measured.
Write j = mk with k squarefree, 2 <= k < N^theta', k coprime to N and
m coprime to k. Then:

  * If j is prime, the only factorisations are (m, k) = (1, j) and
    (j, 1). The second has log 1 = 0 and k = 1 is outside the range,
    so the whole contribution is the m = 1 term with k = j prime,
    and it exists only for j < N^theta'.
  * If j = p^e with e >= 2, the only squarefree k > 1 dividing j is
    k = p, and then m = p^{e-1} shares p with k, so (m, k) = 1 fails.
    Prime powers contribute EXACTLY NOTHING.
  * Everything else -- j with at least two distinct prime factors --
    is the rest.

So the sum splits in two, one of them tiny and explicit: the primes
p < N^theta' with their Lambda(N - p), which is a Goldbach count with
one prime forced below the level, and the composites with at least two
distinct primes. Measuring the two says whether the negative bulk
{#rem:denominator} found is the composite part, and how the explicit
Goldbach piece compares in size to the whole.

BACKS: Remark {#rem:support} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control. The two parts sum to the repository's sum a to
      1e-12 relative at every N of the field.
  V2  Prime powers with e >= 2 contribute exactly zero -- the
      accumulated bucket is 0.0 and not merely small.
  V3  The signs separate: the prime part is positive at every N and
      the composite part is negative at every N.
  V4  And the explicit piece is the smaller one: the prime part's
      exponent in log N is below alpha at two standard errors, and is
      within two standard errors of theta' itself, which is what a
      Goldbach count with one prime below N^theta' should give.

REFUTATION RULE (fixed before the run)

  V1  REFUTED above 1e-12 relative anywhere. Then the split is not of
      the sum {#rem:denominator} measured. THIS ONE GATES.
  V2  REFUTED on any nonzero. Then the coprimality argument above is
      wrong and the two-way split is really three-way.
  V3  REFUTED if either sign fails anywhere. Then the negative total
      is not the composites' doing and the sign has to be attributed
      elsewhere.
  V4  REFUTED if the prime part is not resolved below alpha, or if
      its exponent is not theta' within two standard errors. The
      first would mean the explicit Goldbach piece is what the whole
      sum is, which would make item 4(b) a statement about Goldbach
      counting directly. The second would mean the piece is not of
      the size its own arithmetic predicts, and the prediction would
      have to be found wrong before anything else is read.

  V1 gates. V2 to V4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. The split is by the arithmetic
  class of the index and every part is an exactly computed sum; there
  is no background to detect against and no threshold a null would
  calibrate.
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
OUT = os.path.join(RES, "audit_support.txt")

LO, HI = 200_000, 102_400_000


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA


def read_published():
    """{#rem:denominator}'s alpha and its printing"""
    src = io.open(os.path.join(RES, "audit_denominator.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  \|sum a\|\s+([+-][\d.]+)\s+([\d.]+)\s+[\d.]+"
                  r"\s*$", src, re.M)
    dec = len(m.group(1).split(".")[1])
    return float(m.group(1)), float(m.group(2)), dec


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


def parts(N, lam, mu, sqf, isprime, ppow):
    """sum a split by the arithmetic class of j = mk"""
    PN = SPL.factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    pr = 0.0
    pp = 0.0
    co = 0.0
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in SPL.factor_set(k):
            ms = ms[ms % q != 0]
        vals = lam[N - ms * k] * mu[ms].astype(np.float64) * math.log(k)
        j = ms * k
        pmask = isprime[j]
        qmask = ppow[j] & (~pmask)
        pr += float(vals[pmask].sum())
        pp += float(vals[qmask].sum())
        co += float(vals[~(pmask | qmask)].sum())
    return int(ks.size), pr, pp, co


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

    puba, pubse, dec = read_published()
    NS = family(LO, HI)
    say("read alpha = %+.6f (s.e. %.6f) from "
        "results/audit_denominator.txt" % (puba, pubse))
    say("  the field, the sieve, theta' and the weighted sum are "
        "imported from code/audit_gain_split.py")
    say("  theta' = %.2f" % THETA)
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
    # Lambda is nonzero exactly on the prime powers, and equals log p
    # there; a prime is a prime power whose own Lambda is its log.
    ppow = lam > 0
    isprime = np.zeros(HI + 1, dtype=bool)
    pr_list = SPL.primes_upto(HI)
    isprime[pr_list] = True
    del pr_list

    rows = []
    say()
    say("  N            #k      sum a           prime part      "
        "prime powers  composite part")
    for N in NS:
        nk, pr, pp, co = parts(N, lam, mu, sqf, isprime, ppow)
        tot = pr + pp + co
        rows.append((N, nk, tot, pr, pp, co))
        say("  %-12d %-7d %+-15.1f %+-15.1f %-13.1f %+.1f"
            % (N, nk, tot, pr, pp, co))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))

    # -------------------------------------------------------------- V1
    say()
    say("V1  the control: does the split reproduce sum a?")
    worst = 0.0
    for N, nk, tot, pr, pp, co in rows:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        s = float(a.sum())
        worst = max(worst, abs(tot - s) / max(abs(s), 1e-12))
    v1 = worst <= 1e-12
    say("  worst relative departure over the %d N: %.3e"
        % (len(rows), worst))
    say("  V1 %s   (tol 1e-12 relative)"
        % ("hold" if v1 else "REFUTED"))

    # -------------------------------------------------------------- V2
    say()
    say("V2  do prime powers contribute anything?")
    mx = max(abs(r[4]) for r in rows)
    v2 = mx == 0.0
    say("  the largest absolute prime-power bucket over the field is "
        "%.1f" % mx)
    say("  V2 %s   (cap: exactly zero)"
        % ("hold" if v2 else "REFUTED"))
    if v2:
        say("  as the coprimality forces: the only squarefree k > 1 "
            "dividing p^e is p, and then m = p^{e-1} is not coprime "
            "to it.")

    # -------------------------------------------------------------- V3
    say()
    say("V3  do the signs separate?")
    npos = sum(1 for r in rows if r[3] > 0)
    nneg = sum(1 for r in rows if r[5] < 0)
    v3 = npos == len(rows) and nneg == len(rows)
    say("  the prime part is positive at %d of %d N and the composite "
        "part negative at %d of %d" % (npos, len(rows), nneg,
                                       len(rows)))
    say("  V3 %s" % ("hold" if v3 else "REFUTED"))

    # -------------------------------------------------------------- V4
    say()
    say("V4  how big is the explicit Goldbach piece?")
    ep, rp, sep = fit(x, np.log(np.array([r[3] for r in rows])))
    ec, rc, sec = fit(x, np.log(np.array([abs(r[5]) for r in rows])))
    et, rt, set_ = fit(x, np.log(np.array([abs(r[2]) for r in rows])))
    say("  part            exponent      s.e.        r.m.s.")
    for nm, e, se, rr in (("prime", ep, sep, rp),
                          ("composite", ec, sec, rc),
                          ("total |sum a|", et, set_, rt)):
        say("  %-15s %+-13.6f %-11.6f %.6f" % (nm, e, se, rr))
    say("TSTAT slope_audit_support %.2f" % (abs(ep) / sep))
    say("SPREAD slope_audit_support %.4f" % float(x.max() - x.min()))
    if abs(ep) / sep < 2.0:
        say("UNRESOLVED SIGN slope_audit_support")
    below = (et - ep) > 2.0 * math.sqrt(sep * sep + set_ * set_)
    attheta = abs(ep - THETA) <= 2.0 * sep
    v4 = below and attheta
    say("  the prime part sits below the total by %+.6f, which is "
        "%.2f standard errors of the difference"
        % (ep - et, abs(ep - et) / math.sqrt(sep * sep + set_ * set_)))
    say("  and against theta' = %.2f itself the departure is %+.6f, "
        "%.2f standard errors" % (THETA, ep - THETA,
                                  abs(ep - THETA) / sep))
    say("  V4 %s   (caps 2 standard errors on each)"
        % ("hold" if v4 else "REFUTED"))
    say("  the total here reproduces the published alpha %+.6f to "
        "%.6f" % (puba, abs(et - puba)))
    say("PRINTBOUND audit_support %d %.8f"
        % (dec, 0.5 * 10.0 ** (-dec)))

    say()
    say("=" * 70)
    say("V1 %s  V2 %s  V3 %s  V4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (v1, v2, v3, v4)))

    head = [
        "STATISTIC: sum a = sum_k (log k) H(N;k) split by the",
        "           arithmetic class of the index j = mk into the j",
        "           that are prime, the j that are higher prime",
        "           powers, and the j with at least two distinct",
        "           prime factors; each part's sign at every N and",
        "           each part's least-squares exponent in log N over",
        "           the on-field family to 1.024e8, with the prime",
        "           part's exponent against theta' and against the",
        "           exponent of the whole.",
        "NULL: none is run and none applies. The split is by the",
        "      arithmetic class of the index, every part is an exactly",
        "      computed sum, and there is no background to detect",
        "      against or threshold a null would calibrate.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; j = mk classified by",
        "       primality from the same integer sieve to 102400000",
        "       that gives Lambda and mu; the weighted sum, the sieve",
        "       and theta' are code/audit_gain_split.py's, imported;",
        "       alpha is read from results/audit_denominator.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not v1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
