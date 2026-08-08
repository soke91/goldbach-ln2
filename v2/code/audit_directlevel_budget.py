# -*- coding: utf-8 -*-
r"""
The one truncation exponent still standing, audited the same way the
other two were.

WHAT IS AT STAKE

audit_truncation_exponent.py withdrew two exponents -- M* for the
truncation in m and Q* for the combined modulus -- by sweeping the
tolerance that defined them and watching the fit fall apart.  Both were
first crossings of an ALTERNATING tail, which oscillates while it
decays, so the crossing read the last excursion rather than a trend.

Remark {#rem:directlevel} still quotes K*_H ~ N^{0.7361}, the crossing
of B_H(N;K) = sum_{k<K}(log k)|H(N;k)| above S(N)N, and the conclusion
that [eq:directcond] holds past the square-root barrier at a level
near theta' = 0.74.  That one is structurally different: B_H is a sum
of nonnegative terms, so it is monotone in K and the crossing is
unique, with no excursions to misread.  But it has a free parameter
too -- [eq:directcond] asks for B_H <= (1-eps) S(N) N, and K*_H was
taken at eps = 0.  Scaling the budget to c S(N) N and refitting at
several c is the same audit, and it is the one test that would
distinguish "monotone, so stable" from "stable-looking at the value
that was tried".

BACKS: Remark {#rem:budget} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  Sanity: K*_H(c) is strictly increasing in c at every N, which is
      what monotonicity of B_H forces.
  W2  The exponent is stable in the budget: fitting K*_H(c) ~ N^{e}
      at c = 0.3, 0.5, 1, 2, 3, the spread max(e) - min(e) is under
      0.10 -- against the 1.3033 the withdrawn Q* showed.
  W3  The fits stay good: the correlation exceeds 0.99 at every c --
      against correlations as low as 0.00931 for the withdrawn M*.
  W4  The conclusion survives the sweep: e > 0.5 at every c, so the
      level of [eq:directcond] exceeds one half however the budget is
      scaled.

REFUTATION RULE (fixed before the run)

  W1  REFUTED by a single non-increase. It is forced by nonnegativity,
      so a failure means the walk is wrong.
  W2  REFUTED if the spread reaches 0.10, in which case the third
      exponent goes the way of the other two and Remark
      [rem:directlevel] must be restated.
  W3  REFUTED if any correlation drops to 0.99 or below.
  W4  REFUTED if e drops to 0.5 or below at any c.

  All four gate.

  NO NULL IS RUN and none applies. This is not a detection against a
  background: B_H is a deterministic monotone function of K and the
  question is only whether a definition's free parameter moves the
  fitted exponent. The field's sign controls were run in
  lab_direct_level.py, whose no-cancellation reference established
  that the level is bought by cancellation at all.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_directlevel_budget.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
CLIM = 4_000_000
CS = [0.3, 0.5, 1.0, 2.0, 3.0]
THETAS = [0.56, 0.60, 0.70]


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return pr, lam, mu


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    rows = []
    for N in NS:
        PN = factor_set(N)
        S = twin
        for q in sorted(PN):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        ks = np.array([k for k in range(2, N // 2)
                       if sqf[k] and all(k % q for q in PN)],
                      dtype=np.int64)
        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)
        A = np.empty(ks.size, dtype=np.float64)
        for i, k in enumerate(ks):
            k = int(k)
            r = N % k
            A[i] = f0[r::k].sum() if r else f0[k::k].sum()
        del f0
        cum = np.cumsum(np.log(ks.astype(np.float64)) * np.abs(A))
        rows.append((N, S, ks, cum))
        say("  N = %-10d  #k = %-9d B_H at the cap = %.1f N"
            % (N, ks.size, cum[-1] / N))

    say()
    say("  K*_H at each budget c, the first K with B_H > c S(N) N")
    say("  c         " + "  ".join("N=%-9d" % N for N in NS))
    tab = {}
    mono = True
    for c in CS:
        row = []
        for N, S, ks, cum in rows:
            j = int(np.searchsorted(cum, c * S * N))
            row.append(int(ks[min(j, ks.size - 1)]))
        tab[c] = row
        say("  %-9.2f %s" % (c, "  ".join("%-11d" % v for v in row)))
    for i in range(len(CS) - 1):
        if not all(tab[CS[i]][j] < tab[CS[i + 1]][j]
                   for j in range(len(NS))):
            mono = False
    say("  W1  increasing in c at every N   %s"
        % ("hold" if mono else "REFUTED"))

    say()
    x = np.log(np.array(NS, dtype=float))
    say("  fitted exponents")
    say("  c         e(c)      correlation   K*/sqrt N at the top")
    es, rs = [], []
    for c in CS:
        y = np.log(np.array(tab[c], dtype=float))
        e = float(np.polyfit(x, y, 1)[0])
        r = float(np.corrcoef(x, y)[0, 1])
        es.append(e)
        rs.append(r)
        say("  %-9.2f %-9.4f %-13.5f %.4f"
            % (c, e, r, tab[c][-1] / math.sqrt(NS[-1])))

    say()
    sp = max(es) - min(es)
    w2 = sp < 0.10
    say("W2  spread of the exponent over the budgets: %.4f  (cap 0.10)"
        "   %s" % (sp, "hold" if w2 else "REFUTED"))
    w3 = min(rs) > 0.99
    say("W3  worst fit correlation: %.5f  (floor 0.99)   %s"
        % (min(rs), "hold" if w3 else "REFUTED"))
    w4 = all(e > 0.5 for e in es)
    say("W4  the exponent stays above 0.5: min %.4f   %s"
        % (min(es), "hold" if w4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The threshold-free version, which is")
    say("  what audit_truncation_exponent.py said should be reported")
    say("  alongside any crossing: B_H(N;K)/(S(N)N) at a FIXED cut.")
    say("  cut          " + "  ".join("N=%-9d" % N for N in NS))
    for th in THETAS:
        row = []
        for N, S, ks, cum in rows:
            j = int(np.searchsorted(ks, int(N ** th)))
            row.append(cum[min(j, cum.size - 1)] / (S * N))
        say("  K<N^%-8.2f %s" % (th, "  ".join("%-11.4f" % v
                                               for v in row)))
    say("  Read down: at a fixed exponent the ratio falls with N, so")
    say("  the crossing moves outward -- the same qualitative")
    say("  conclusion as the fit, reached with no free parameter. That")
    say("  is the difference from the withdrawn cases, where the")
    say("  parameter-free tables disagreed with the fitted exponents.")

    say()
    say("  Cross-check lines. lab_direct_level.py computes the same")
    say("  crossing at c = 1 by an independent walk.")
    for i, N in enumerate(NS):
        say("AGREE kstar_H N=%d %d 0.02" % (N, tab[1.0][i]))

    say()
    say("=" * 70)
    ok = mono and w2 and w3 and w4
    say("the crossing of the absolute demand is stable in its budget "
        "and its exponent survives" if ok else "REFUTED")

    head = [
        "STATISTIC: K*_H(N;c), the first K at which",
        "           B_H(N;K) = sum_{k<K}(log k)|H(N;k)| exceeds c S(N)N,",
        "           at c = 0.3, 0.5, 1, 2, 3; the exponent fitted from",
        "           K*_H against N at each c, with its correlation and",
        "           spread; and the parameter-free ratio B_H/(S(N)N) at",
        "           the fixed cuts K < N^0.56, N^0.60, N^0.70.",
        "NULL: none is run and none applies. B_H is a deterministic",
        "      monotone function of K, so there is no background to",
        "      detect against; the question is only whether a free",
        "      parameter in the definition moves the fitted exponent.",
        "      The sign controls for this field were run in",
        "      lab_direct_level.py, whose mu-squared reference",
        "      established that the level is bought by cancellation.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k over the squarefree",
        "       k < N/2 coprime to N, walked upward; S(N) from an Euler",
        "       product at the fixed bound 4e6.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
