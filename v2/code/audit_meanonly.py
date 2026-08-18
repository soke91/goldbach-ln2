# -*- coding: utf-8 -*-
r"""
Is the surviving correlation a correlation at all, or only a mean?

WHAT IS AT STAKE

{#rem:support} left item 4(b)'s denominator as one object: a
correlation of Lambda against Lambda_K, the mu * log convolution with
the k-index truncated, carried almost entirely by indices with two or
more distinct prime factors. For such j the untruncated convolution
is Lambda(j) = 0, so Lambda_K(j) there is exactly minus the part the
restrictions threw away -- a truncation defect of an identity that
sums to zero.

Reordering, sum a = sum_q (log q) Lambda_K(N - q) over prime powers q,
so it is Lambda_K sampled at the shifted primes. Nothing has asked
whether the sampling matters. If Lambda_K has a nonzero mean and the
shifted primes see only that mean, then

    sum a  ~  (psi(N)/N) * sum_j Lambda_K(j)

and the whole object is a divisor-sum average with no primes in it.
Item 4(b) would then be a statement about the mean of a truncated
convolution -- an elementary question -- rather than about a
correlation with the primes. If instead the shifted primes see
something the mean does not, the excess is the arithmetic content and
it is what has to be bounded.

The control is the one M4 asks for: compute the sum both ways, with
Lambda(N - j) as measured and with Lambda(N - j) replaced by its own
mean over the same range, and compare.

BACKS: Remark {#rem:meanonly} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The control. The sum recomputed here reproduces
      results/audit_support.txt's total inside the bound its printing
      forces, at every N of the field.
  W2  The mean-only value has the same sign as the sum at every N,
      so the mean of Lambda_K over the restricted index is negative
      as the sum is.
  W3  And it is most of it: |sum a - mean-only| is at most half of
      |sum a| at every N.
  W4  The residual is of smaller order: its least-squares exponent in
      log N is below alpha at two standard errors.

REFUTATION RULE (fixed before the run)

  W1  REFUTED outside the printing bound anywhere. Then this is not
      the sum {#rem:support} measured. THIS ONE GATES.
  W2  REFUTED if the signs differ anywhere. Then the mean does not
      even point the same way and the reading below is wrong from the
      start.
  W3  REFUTED if the residual exceeds half anywhere. Then the shifted
      primes see something the mean does not, the object is a genuine
      prime correlation, and item 4(b) keeps its arithmetic content.
      That is the outcome that leaves the problem where it was.
  W4  REFUTED if the residual is not resolved below alpha. Then the
      correlation part grows as fast as the whole and W3's accounting
      is a statement about this range only.

  W1 gates. W2 to W4 are the measurement and do not gate.

  THE NULL IS THE MEASUREMENT. The mean-only arm replaces
  Lambda(N - j) by its mean over the same range while keeping the
  index set, the weights mu(m) log k and the truncation exactly as
  they are -- M3's requirement that a null preserve the structure of
  the field it nulls. Nothing is randomised, so no seed is needed.
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
OUT = os.path.join(RES, "audit_meanonly.txt")

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
    """{#rem:support}'s totals and alpha, and the printing"""
    src = io.open(os.path.join(RES, "audit_support.txt"),
                  encoding="utf-8").read()
    tot, dec = {}, 0
    for m in re.finditer(r"^  (\d{5,})\s+\d+\s+([+-][\d.]+)\s+"
                         r"[+-][\d.]+\s+[\d.]+\s+[+-][\d.]+\s*$",
                         src, re.M):
        tot[int(m.group(1))] = float(m.group(2))
        dec = max(dec, len(m.group(2).split(".")[1]))
    m = re.search(r"^  total \|sum a\|\s+([+-][\d.]+)\s+([\d.]+)",
                  src, re.M)
    return tot, dec, float(m.group(1)), float(m.group(2))


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


def both(N, lam, mu, sqf, psibar):
    """the sum as measured, and with Lambda(N-j) at its own mean"""
    PN = SPL.factor_set(N)
    K = int(N ** THETA)
    obs = 0.0
    wsum = 0.0
    for k in range(2, K):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in SPL.factor_set(k):
            ms = ms[ms % q != 0]
        g = mu[ms].astype(np.float64)
        lk = math.log(k)
        obs += lk * float((lam[N - ms * k] * g).sum())
        wsum += lk * float(g.sum())
    return obs, psibar * wsum


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

    pubtot, dec, puba, pubse = read_published()
    NS = family(LO, HI)
    say("read %d totals and alpha = %+.6f (s.e. %.6f) from "
        "results/audit_support.txt" % (len(pubtot), puba, pubse))
    say("  the field, the sieve, theta' and the index set are "
        "imported from code/audit_gain_split.py")
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
    say("  N            sum a           mean-only       residual"
        "        residual/|sum a|  mean of Lambda")
    for N in NS:
        psibar = float(lam[1:N].sum()) / float(N - 1)
        obs, mo = both(N, lam, mu, sqf, psibar)
        res = obs - mo
        rows.append((N, obs, mo, res, abs(res) / abs(obs), psibar))
        say("  %-12d %+-15.1f %+-15.1f %+-15.1f %-17.4f %.6f"
            % (N, obs, mo, res, abs(res) / abs(obs), psibar))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))

    # -------------------------------------------------------------- W1
    say()
    say("W1  the control against results/audit_support.txt")
    rnd = 0.5 * 10.0 ** (-dec)
    worst = max(abs(r[1] - pubtot[r[0]]) for r in rows
                if r[0] in pubtot)
    w1 = worst <= rnd
    say("  worst departure over %d N: %.6f; the table prints %d "
        "decimals, so the bound is %.8f"
        % (len(rows), worst, dec, rnd))
    say("PRINTBOUND audit_meanonly %d %.8f" % (dec, rnd))
    say("  W1 %s   (cap: the printing bound)"
        % ("hold" if w1 else "REFUTED"))

    # -------------------------------------------------------------- W2
    say()
    say("W2  does the mean point the same way?")
    same = sum(1 for r in rows if (r[1] < 0) == (r[2] < 0))
    w2 = same == len(rows)
    say("  the two agree in sign at %d of %d N" % (same, len(rows)))
    say("  W2 %s" % ("hold" if w2 else "REFUTED"))

    # -------------------------------------------------------------- W3
    say()
    say("W3  how much of the sum is the mean?")
    fr = np.array([r[4] for r in rows])
    w3 = bool((fr <= 0.5).all())
    say("  the residual is %.4f to %.4f of the sum in absolute value"
        % (float(fr.min()), float(fr.max())))
    say("  at the top N: sum %+.1f, mean-only %+.1f, residual %+.1f"
        % (rows[-1][1], rows[-1][2], rows[-1][3]))
    say("  W3 %s   (cap: a half at any N)"
        % ("hold" if w3 else "REFUTED"))

    # -------------------------------------------------------------- W4
    say()
    say("W4  and how does the residual grow?")
    eo, ro, seo = fit(x, np.log(np.array([abs(r[1]) for r in rows])))
    em, rm, sem = fit(x, np.log(np.array([abs(r[2]) for r in rows])))
    er, rr, ser = fit(x, np.log(np.array([abs(r[3]) for r in rows])))
    say("  quantity      exponent      s.e.        r.m.s.")
    for nm, e, se, rrms in (("sum a", eo, seo, ro),
                            ("mean-only", em, sem, rm),
                            ("residual", er, ser, rr)):
        say("  %-13s %+-13.6f %-11.6f %.6f" % (nm, e, se, rrms))
    d = math.sqrt(ser * ser + seo * seo)
    w4 = (eo - er) > 2.0 * d
    say("  the residual sits below the sum by %+.6f, which is %.2f "
        "standard errors of the difference" % (er - eo, abs(er - eo) / d))
    say("TSTAT slope_audit_meanonly %.2f" % (abs(er) / ser))
    say("SPREAD slope_audit_meanonly %.4f" % float(x.max() - x.min()))
    if abs(er) / ser < 2.0:
        say("UNRESOLVED SIGN slope_audit_meanonly")
    say("  W4 %s   (cap 2 standard errors)"
        % ("hold" if w4 else "REFUTED"))
    say("  and the sum's exponent here against the published alpha "
        "%+.6f: %+.6f" % (puba, eo - puba))

    say()
    say("=" * 70)
    say("W1 %s  W2 %s  W3 %s  W4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (w1, w2, w3, w4)))

    head = [
        "STATISTIC: sum a = sum_k (log k) sum_m mu(m) Lambda(N - mk)",
        "           over the restricted index, computed as measured",
        "           and again with Lambda(N - mk) replaced by its own",
        "           mean over 1 <= n < N, the index set and the",
        "           weights mu(m) log k left exactly as they are; the",
        "           residual between the two, its share of the sum,",
        "           and each one's least-squares exponent in log N",
        "           over the on-field family to 1.024e8.",
        "NULL: the mean-only arm IS the null, and it is what M3 asks",
        "      for: the structure being nulled is the primality of",
        "      N - j, so that alone is replaced by its mean while the",
        "      index set, the truncation and the weights are kept",
        "      exactly. Nothing is randomised and no seed is needed.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the index set, the sieve",
        "       and theta' are code/audit_gain_split.py's, imported;",
        "       the totals and alpha are read from",
        "       results/audit_support.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not w1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
