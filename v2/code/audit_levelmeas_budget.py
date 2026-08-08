# -*- coding: utf-8 -*-
r"""
The withdrawn level measurement, audited on its budget and given the
control it should have had.

WHAT IS AT STAKE

Remark {#rem:levelmeas} defined

    K*(N) := max { K : sum_{k<K,(k,N)=1} (log k)|E_mu(N;k)|
                       <= S(N)(1 - A(N)) N },

the largest truncation at which [eq:nolog] still holds, measured
K* = 319, 537, 767, 1353, 2319 with K*/sqrt N rising through 1, and
then WITHDREW the measurement because a coin reached a higher K* than
mu at every N.

Two things have happened since.  Remark {#rem:whycoinwins} showed that
a coin is better than mu here BY CONSTRUCTION -- by [eq:dilate] the
coin's progression sum is a sum of independent signs and gets
square-root cancellation for free, while mu's is the dilated wall --
so the coin is a competitor and not a null.  And
audit_truncation_exponent.py withdrew two exponents by sweeping the
free parameter that defined them, while audit_directlevel_budget.py
saved a third the same way, the difference being that the surviving
crossing was of a MONOTONE sum.

B(N;K) = sum_{k<K}(log k)|E_mu(N;k)| is also a sum of nonnegative
terms, so this crossing is in the surviving class structurally.  This
script does both things at once: it sweeps the budget, and it runs the
control the withdrawal should have used -- mu replaced by mu^2, same
support, every sign +1, main term still subtracted, so no cancellation
is available anywhere.

BACKS: Remark {#rem:levelaudit} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  Sanity: K*(c) is strictly increasing in the budget factor c at
      every N, which nonnegativity of the summands forces.
  V2  The exponent is stable in the budget: fitting K*(c) ~ N^{e} at
      c = 0.3, 0.5, 1, 2, 3, the spread max(e) - min(e) is under 0.10.
  V3  The fits stay good: correlation above 0.99 at every c.
  V4  The control gets nowhere: with mu replaced by mu^2 the crossing
      has K*/sqrt N under 0.1 at every N, as it does for the other
      crossing in lab_direct_level.py. The level is bought by
      cancellation and not by the size of the terms.

REFUTATION RULE (fixed before the run)

  V1  REFUTED by a single non-increase; it is forced, so a failure
      means the walk is wrong.
  V2  REFUTED if the spread reaches 0.10, in which case N^{0.7057}
      joins the two withdrawn exponents.
  V3  REFUTED if any correlation drops to 0.99 or below.
  V4  REFUTED if the control reaches 0.1 at any N, in which case the
      measurement is not evidence about mu and the withdrawal stands
      on other grounds than the ones it gave.

  All four gate.

  THE CONTROL is V4 and it is deliberately NOT a coin. Remark
  {#rem:whycoinwins} established that a coin beats mu here by
  construction, so the coin comparison the withdrawal rested on could
  only ever come out the way it did. Replacing mu by mu^2 goes the
  other way: same support, no signs at all, so it says how much of the
  level is bought by cancellation rather than by size.
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
OUT = os.path.join(ROOT, "results", "audit_levelmeas_budget.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
CLIM = 4_000_000
CS = [0.3, 0.5, 1.0, 2.0, 3.0]
THETAS = [0.50, 0.56, 0.60]


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


def totients(n):
    phi = np.arange(n + 1, dtype=np.int64)
    for p in primes_upto(n):
        p = int(p)
        phi[p::p] -= phi[p::p] // p
    return phi


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


def cumB(N, ks, iph, lam, sign):
    """cumsum of (log k)|A(N;k) - C(N)/phi(k)| for the given sign array."""
    f = np.zeros(N, dtype=np.float64)
    idx = np.arange(1, N, dtype=np.int64)
    f[1:] = lam[1:N] * sign[N - idx]
    C = float(f.sum())
    A = np.empty(ks.size, dtype=np.float64)
    for i, k in enumerate(ks):
        k = int(k)
        r = N % k
        A[i] = f[r::k].sum() if r else f[k::k].sum()
    del f
    return np.cumsum(np.log(ks.astype(np.float64))
                     * np.abs(A - C / iph))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    phi = totients(NMAX)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    muf = mu.astype(np.float64)
    onef = np.zeros(NMAX + 1, dtype=np.float64)
    onef[sqf] = 1.0

    rows = []
    for N in NS:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))
        thr1 = S_ * (1.0 - A_) * N
        ks = np.array([k for k in range(2, N // 2)
                       if sqf[k] and all(k % q for q in PN)],
                      dtype=np.int64)
        iph = phi[ks].astype(np.float64)
        cm = cumB(N, ks, iph, lam, muf)
        co = cumB(N, ks, iph, lam, onef)
        rows.append((N, S_, A_, thr1, ks, cm, co))
        say("  N = %-10d  S(1-A) = %.6f   B at the cap = %.1f N"
            % (N, S_ * (1.0 - A_), cm[-1] / N))

    say()
    say("  K*(N) at each budget c, the first K with B > c S(1-A) N")
    say("  c         " + "  ".join("N=%-9d" % N for N in NS))
    tab = {}
    for c in CS:
        row = []
        for N, S_, A_, thr1, ks, cm, co in rows:
            j = int(np.searchsorted(cm, c * thr1))
            row.append(int(ks[min(j, ks.size - 1)]))
        tab[c] = row
        say("  %-9.2f %s" % (c, "  ".join("%-11d" % v for v in row)))
    v1 = all(all(tab[CS[i]][j] < tab[CS[i + 1]][j]
                 for j in range(len(NS)))
             for i in range(len(CS) - 1))
    say("  V1  increasing in c at every N   %s"
        % ("hold" if v1 else "REFUTED"))

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
    v2 = sp < 0.10
    say("V2  spread of the exponent over the budgets: %.4f  (cap 0.10)"
        "   %s" % (sp, "hold" if v2 else "REFUTED"))
    v3 = min(rs) > 0.99
    say("V3  worst fit correlation: %.5f  (floor 0.99)   %s"
        % (min(rs), "hold" if v3 else "REFUTED"))

    say()
    say("V4  the control: mu replaced by mu^2, same support, no signs")
    say("  N            K*_mu     K*_one    K*_one/sqrt N")
    v4 = True
    for i, (N, S_, A_, thr1, ks, cm, co) in enumerate(rows):
        j = int(np.searchsorted(co, thr1))
        ko = int(ks[min(j, ks.size - 1)])
        rr = ko / math.sqrt(N)
        if rr >= 0.1:
            v4 = False
        say("  %-12d %-9d %-9d %.4f" % (N, tab[1.0][i], ko, rr))
    say("  V4 %s" % ("hold" if v4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The parameter-free version:")
    say("  B(N;K)/(S(1-A)N) at a FIXED cut.")
    say("  cut          " + "  ".join("N=%-9d" % N for N in NS))
    for th in THETAS:
        row = []
        for N, S_, A_, thr1, ks, cm, co in rows:
            j = int(np.searchsorted(ks, int(N ** th)))
            row.append(cm[min(j, cm.size - 1)] / thr1)
        say("  K<N^%-8.2f %s" % (th, "  ".join("%-11.4f" % v
                                               for v in row)))
    say("  and where the exponent sits against one half:")
    say("  e(c) above 0.5 at every c: %s   min %.4f"
        % (all(e > 0.5 for e in es), min(es)))

    say()
    say("  Cross-check lines. lab_level_coin_null.py computes the same")
    say("  crossing at c = 1 by an independent walk.")
    for i, N in enumerate(NS):
        say("AGREE kstar_nolog N=%d %d 0.02" % (N, tab[1.0][i]))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the withdrawn level measurement is stable in its budget and "
        "its control is the one that was missing" if ok else "REFUTED")

    head = [
        "STATISTIC: K*(N;c), the first K at which",
        "           B(N;K) = sum_{k<K}(log k)|E_mu(N;k)| exceeds",
        "           c S(N)(1-A(N))N, at c = 0.3, 0.5, 1, 2, 3; the",
        "           exponent fitted from K* against N at each c with its",
        "           correlation and spread; the same crossing with mu",
        "           replaced by mu^2; and the parameter-free ratio",
        "           B/(S(1-A)N) at the fixed cuts K < N^0.50, N^0.56,",
        "           N^0.60.",
        "NULL: mu replaced by mu^2 -- same support, every sign +1, main",
        "      term C(N)/phi(k) still subtracted, so no cancellation is",
        "      available anywhere. It is deliberately not a coin:",
        "      [rem:whycoinwins] established that a coin beats mu here",
        "      by construction, so the coin comparison the original",
        "      withdrawal rested on could only come out one way.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k over the squarefree",
        "       k < N/2 coprime to N, walked upward; S(N) and A(N) from",
        "       Euler products at the fixed bound 4e6.",
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
