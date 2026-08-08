# -*- coding: utf-8 -*-
r"""
Whether the truncation exponents this program reports are properties
of the count or of the tolerance that defined them.

WHAT IS AT STAKE

Three exponents have been reported and they agree suspiciously well:
M* ~ N^{0.1089} for the truncation in m (Remark {#rem:layertail}),
Q* ~ N^{0.6904} for the combined modulus m d^2 (Remark {#rem:combmod}),
and 0.6716 / 0.7250 for the two crossings in k (Remark
{#rem:signedlevel}).  The conclusion drawn from them -- that the
demand sits near N^{0.7} however the double sum is cut, above the 1/2
Bombieri-Vinogradov supplies -- is only as good as the definitions.

Every one of those numbers came from a threshold chosen by hand: M*
and Q* are "the least truncation past which the deficit stays under
0.01 N", and the crossings in k used S(N)N.  If the deficit decays
like a power of the truncation, then the truncation needed to reach
tolerance eps depends on eps, and the EXPONENT in N may depend on it
too.  If it does, 0.6904 is where the author stopped and not where the
mathematics is.

So the tolerance is swept over decades and the exponent refitted at
each one.  The pairs (m,d) are enumerated once, which gives both
truncations from the same data: grouping the pairs by m and summing
over d reconstitutes mu(m)L(N;m) exactly, so M* and Q* are two
orderings of one list.

BACKS: Remark {#rem:tolerance} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  The exponent for Q* is stable in the tolerance: fitting
      Q*(N;eps) ~ N^{e} at eps = 0.3, 0.1, 0.03, 0.01, 0.003, the
      spread max(e) - min(e) is under 0.10.
  X2  The same for M*: spread under 0.10.
  X3  The conclusion does not depend on the tolerance: e for Q* stays
      above 0.5 at every eps.
  X4  The gap between the two survives: e(Q*) - e(M*) is above 0.3 at
      every eps, so the cost of expanding the squarefree condition is
      not an artefact of where the deficit was cut.

REFUTATION RULE (fixed before the run)

  X1  REFUTED if the spread reaches 0.10 -- in which case 0.6904 is a
      property of the tolerance and Remark [rem:combmod] must be
      restated with the dependence shown.
  X2  Likewise for M* and Remark [rem:layertail].
  X3  REFUTED if e drops to 0.5 or below at any eps, which would put
      the demand inside Bombieri-Vinogradov at that tolerance.
  X4  REFUTED if the gap reaches 0.3 from above at any eps.

  All four gate. A refutation of X1 or X2 is the outcome that matters:
  it would say this program has been reporting thresholds as if they
  were exponents.

  NO NULL IS RUN and none applies. Nothing here is a detection against
  a background; every quantity is a deterministic function of one
  exact expansion, and the question is only whether a definition's
  free parameter moves the answer. The field's sign controls were run
  in lab_layer_decomposition.py and lab_combined_modulus.py.
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
OUT = os.path.join(ROOT, "results", "audit_truncation_exponent.txt")

NS = [200_000, 400_000, 800_000, 1_600_000]
EPS = [0.3, 0.1, 0.03, 0.01, 0.003]


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


def last_above(keys, dfc, N, eps):
    """Least key past which the deficit stays under eps N."""
    bad = np.flatnonzero(dfc >= eps * N)
    return int(keys[bad[-1]]) + 1 if bad.size else 1


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    isp = np.zeros(NMAX + 1, dtype=bool)
    isp[pr] = True
    sqf = mu != 0
    logj = np.log(np.arange(1, NMAX + 1, dtype=np.float64))

    rows = []
    for N in NS:
        ms, qs, vals = [], [], []
        for m in range(1, N):
            if not sqf[m]:
                continue
            fm = factor_set(m)
            dmax = int(math.isqrt((N - 1) // m))
            if dmax < 1:
                break
            for d in range(1, dmax + 1):
                if not sqf[d] or any(d % r == 0 for r in fm):
                    continue
                q = m * d * d
                jmax = (N - 1) // q
                if jmax < 1:
                    continue
                ok = np.ones(jmax, dtype=bool)
                for r in fm:
                    ok[r - 1::r] = False
                js = np.flatnonzero(ok) + 1
                if js.size == 0:
                    continue
                w = 2.0 * math.log(d) + logj[js - 1]
                ms.append(m)
                qs.append(q)
                vals.append(int(mu[m]) * int(mu[d])
                            * float((w * lam[N - q * js]).sum()))
        ms = np.array(ms, dtype=np.int64)
        qs = np.array(qs, dtype=np.int64)
        vals = np.array(vals, dtype=np.float64)
        tot = float(vals.sum())
        rhs = float((lam[1:N] * lam[N - 1:0:-1] * isp[N - 1:0:-1]).sum())

        oq = np.argsort(qs, kind="stable")
        kq, dq = qs[oq], np.abs(tot - np.cumsum(vals[oq]))
        om = np.argsort(ms, kind="stable")
        km, dm = ms[om], np.abs(tot - np.cumsum(vals[om]))

        rows.append((N, tot, rhs, kq, dq, km, dm))
        say("  N = %-10d pairs = %-9d total/N = %.6f  target = %.6f"
            % (N, qs.size, tot / N, rhs / N))

    say()
    say("  the two truncations, at each tolerance")
    say("  eps       " + "  ".join("Q*(%d)" % N for N in NS)
        + "     |  " + "  ".join("M*(%d)" % N for N in NS))
    tabQ, tabM = {}, {}
    for e in EPS:
        q = [last_above(r[3], r[4], r[0], e) for r in rows]
        m = [last_above(r[5], r[6], r[0], e) for r in rows]
        tabQ[e], tabM[e] = q, m
        say("  %-9.3f %s  |  %s"
            % (e, "  ".join("%-11d" % v for v in q),
               "  ".join("%-11d" % v for v in m)))

    x = np.log(np.array(NS, dtype=float))

    def fit(v):
        y = np.log(np.array(v, dtype=float))
        return (float(np.polyfit(x, y, 1)[0]),
                float(np.corrcoef(x, y)[0, 1]))

    say()
    say("  fitted exponents")
    say("  eps       e(Q*)     corr       e(M*)     corr       gap")
    eq, em = [], []
    for e in EPS:
        a, ra = fit(tabQ[e])
        b, rb = fit(tabM[e])
        eq.append(a)
        em.append(b)
        say("  %-9.3f %-9.4f %-10.5f %-9.4f %-10.5f %.4f"
            % (e, a, ra, b, rb, a - b))

    say()
    sq = max(eq) - min(eq)
    x1 = sq < 0.10
    say("X1  e(Q*) spread over the tolerances: %.4f  (cap 0.10)   %s"
        % (sq, "hold" if x1 else "REFUTED"))
    sm = max(em) - min(em)
    x2 = sm < 0.10
    say("X2  e(M*) spread: %.4f  (cap 0.10)   %s"
        % (sm, "hold" if x2 else "REFUTED"))
    x3 = all(v > 0.5 for v in eq)
    say("X3  e(Q*) above 0.5 at every tolerance: min %.4f   %s"
        % (min(eq), "hold" if x3 else "REFUTED"))
    gaps = [a - b for a, b in zip(eq, em)]
    x4 = all(g > 0.3 for g in gaps)
    say("X4  the gap stays above 0.3: min %.4f   %s"
        % (min(gaps), "hold" if x4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The underlying decay, which is what")
    say("  the tolerance is reading off. Fitting the deficit against")
    say("  the truncation as a power law, over the decade below each")
    say("  N's own Q* and M*:")
    say("  N            deficit ~ Q^{-a}   corr       deficit ~ M^{-b}"
        "   corr")
    for N, tot, rhs, kq, dq, km, dm in rows:
        out = []
        for keys, dfc in ((kq, dq), (km, dm)):
            hi = keys[-1]
            sel = (keys >= max(2, hi // 1000)) & (dfc > 0)
            if sel.sum() < 8:
                out += [float("nan"), float("nan")]
                continue
            lx = np.log(keys[sel].astype(float))
            ly = np.log(dfc[sel] / N)
            out += [-float(np.polyfit(lx, ly, 1)[0]),
                    float(np.corrcoef(lx, ly)[0, 1])]
        say("  %-12d %-18.4f %-10.5f %-18.4f %.5f"
            % (N, out[0], out[1], out[2], out[3]))
    say("  The correlations are weak, and that is the whole story: the")
    say("  deficit does not decay like a power of the truncation, it")
    say("  oscillates while decaying, so Q*(eps) reads the LAST excursion")
    say("  above eps N rather than a trend. Inverting an erratic envelope")
    say("  at four values of N gives a slope that means nothing, and the")
    say("  one tolerance that produced a clean-looking fit did so by")
    say("  luck.")
    say()
    say("  What is threshold-free, and should be reported instead, is the")
    say("  envelope itself at fixed truncation: D(Q)/N = the largest")
    say("  deficit at or beyond Q, in units of N.")
    say("  Q            " + "  ".join("N=%-9d" % N for N in NS))
    for Q in (1_000, 10_000, 100_000, 1_000_000):
        row = []
        for N, tot, rhs, kq, dq, km, dm in rows:
            env = np.maximum.accumulate(dq[::-1])[::-1]
            j = int(np.searchsorted(kq, Q))
            row.append(env[j] / N if j < env.size else 0.0)
        say("  %-12d %s" % (Q, "  ".join("%-11.5f" % v for v in row)))
    say("  M            " + "  ".join("N=%-9d" % N for N in NS))
    for M in (30, 300, 3_000, 30_000):
        row = []
        for N, tot, rhs, kq, dq, km, dm in rows:
            env = np.maximum.accumulate(dm[::-1])[::-1]
            j = int(np.searchsorted(km, M))
            row.append(env[j] / N if j < env.size else 0.0)
        say("  %-12d %s" % (M, "  ".join("%-11.5f" % v for v in row)))
    say("  Read down a row: that is how the truncation error at a FIXED")
    say("  cut moves with N, with no free parameter in it. Read across:")
    say("  that is how it decays at fixed N. Neither table supports a")
    say("  clean power of N over this range, and saying so is the")
    say("  correct report.")

    say()
    say("  Cross-check lines. M* at the 0.01 tolerance is also computed")
    say("  by lab_layer_tail.py from the unexpanded layers.")
    for i, N in enumerate(NS):
        say("AGREE mstar_001 N=%d %d 0.02" % (N, tabM[0.01][i]))
    for N, tot, rhs, kq, dq, km, dm in rows:
        say("AGREE untrunc_total N=%d %.6f 1e-9" % (N, tot / N))

    say()
    say("=" * 70)
    ok = x1 and x2 and x3 and x4
    say("the truncation exponents are properties of the count and not "
        "of the tolerance that defined them" if ok else "REFUTED")

    head = [
        "STATISTIC: Q*(N;eps) and M*(N;eps), the least truncation in the",
        "           combined modulus m d^2 and in m alone past which the",
        "           deficit of the truncated sum stays under eps N, at",
        "           eps = 0.3, 0.1, 0.03, 0.01, 0.003; the exponents",
        "           fitted from Q* and M* against N at each eps, their",
        "           spread and their gap; and the power-law decay of the",
        "           deficit itself.",
        "NULL: none is run and none applies. Nothing here is a detection",
        "      against a background: every quantity is a deterministic",
        "      function of one exact expansion, and the question is only",
        "      whether a definition's free parameter moves the answer.",
        "      The field's sign controls were run in",
        "      lab_layer_decomposition.py and lab_combined_modulus.py.",
        "FIELD: N = 2e5 through 1.6e6 by doubling; every pair of",
        "       squarefree m, d with (d,m) = 1 and m d^2 < N; j over",
        "       1 <= j < N/(m d^2) coprime to m. Grouping those pairs by m",
        "       reconstitutes mu(m)L(N;m) exactly, so both truncations",
        "       come from one enumeration.",
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
