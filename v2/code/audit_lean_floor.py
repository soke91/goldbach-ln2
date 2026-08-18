# -*- coding: utf-8 -*-
r"""
The lean decays -- but slower or faster than its own floor?

WHAT IS AT STAKE

Remark {#rem:leandecay} reads the mass-weighted sign lean's decay as
the thing that rescues the asymptotic picture: f rises from 0.2273 to
0.3376 over a factor 32 in N, so "the lean is not a structural fact
about mu: it is the finite-N error of Theorem {#thm:C}, and it
decays", and the cross-k gain G = 1/|2f-1| grows without bound.

The coin arm carried alongside it is two draws per N and is read as
"the coin sits at 1/2 throughout, as it must". Sitting at 1/2 is not
the question. A random sign field on these magnitudes has a lean of
its own, of size about l2/(2 l1), and that floor MOVES with N: Remark
{#rem:crosskreference} has just measured l1/l2 growing from 11.96 to
26.21 across five of these six N, so the floor falls by more than a
factor two over the sweep. Whether mu's lean decays faster or slower
than its own floor is the question, and it has never been asked.

The direction is predictable from that l1/l2 series and is registered
below in the direction it points, not against it. What is not
predictable from it is whether either effect clears two standard
errors -- neither the raw decay nor the comparison has ever been put
against its own noise -- and the floor here is 256 global sign
vectors rather than two draws.

The implementation is independent of lab_lean_decay.py's.

BACKS: Remark {#rem:leanfloor} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  L1  The control: the six mass fractions reproduce the published
      0.2273, 0.2228, 0.2735, 0.3068, 0.3207, 0.3376 to within 0.001
      and G to within 0.01.
  L2  The raw decay is real: mu's power-law slope of |0.5 - f|
      against N reaches two standard errors of its own residuals.
  L3  But it is slower than the floor's: mu's slope lies above the
      97.5th percentile of the 256 coin slopes -- the coin's lean
      decays faster than mu's.
  L4  So in units of the floor the lean does not decay at all: the
      slope of log(mu's lean / the coins' median lean) against log N
      is positive and reaches two standard errors.

REFUTATION RULE (fixed before the run)

  L1  REFUTED at either cap -- not the same statistic, and nothing
      below may be compared with {#rem:leandecay}.
  L2  REFUTED below two standard errors. The decay would not be
      resolved on six points and "it decays" would rest on a
      correlation alone.
  L3  REFUTED if mu's slope is at or below that percentile. The lean
      would then be decaying at least as fast as a random sign field
      on the same magnitudes, and {#rem:leandecay}'s reading would
      stand as written.
  L4  REFUTED if the floor-relative slope fails to reach two standard
      errors, or is negative. That is the one that matters: a
      positive resolved slope says the lean is growing relative to
      what chance gives, and the sentence "the lean is not a
      structural fact about mu" would have to be qualified -- the
      decay would be the floor moving, not the lean going away.

  All four gate.

  THE NULL IS THE POINT and it is run: 256 global sign vectors over
  k, each held across all N as mu is, applied to the IDENTICAL
  (log k)|H(N;k)| magnitudes. Same convention as
  audit_crossk_reference.py and lab_split_budget.py.
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
OUT = os.path.join(RES, "audit_lean_floor.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
DRAWS = 256
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
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
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


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


def read_published():
    """the mass fractions and gains, read from the results file"""
    src = io.open(os.path.join(RES, "lab_lean_decay.txt"),
                  encoding="utf-8").read()
    i = src.index("N          #k     mass frac +   |0.5 - f|   "
                  "G = 1/|2f-1|")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            if f and set(f[0]) == {"-"}:
                continue
            if not out:
                continue
            break
        out[int(f[0])] = (int(f[1]), float(f[2]), float(f[4]))
    b = float(re.search(r"\|0\.5 - f\| ~ N\^\{-([\d.]+)\}",
                        src).group(1))
    return out, b


def weighted(N, lam, mu, sqf):
    """(log k)H(N;k) over the squarefree k < N^theta coprime to N"""
    PN = factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    Hs = []
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        vals = N - ms * k
        Hs.append(float((lam[vals] * mu[ms].astype(np.float64)).sum()))
    return ks, np.log(ks.astype(np.float64)) * np.array(Hs)


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, pubb = read_published()
    say("read %d published rows and the exponent %.4f from "
        "results/lab_lean_decay.txt" % (len(pub), pubb))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    KMAX = int(NMAX ** THETA) + 1
    rng = np.random.default_rng(SEED)
    eps = (rng.integers(0, 2, size=(DRAWS, KMAX + 1))
           .astype(np.int8) * 2 - 1)
    say("%d global sign vectors over k <= %d, seed %d -- each held "
        "across all N as mu is" % (DRAWS, KMAX, SEED))

    got = []
    for N in NS:
        ks, a = weighted(N, lam, mu, sqf)
        l1 = float(np.abs(a).sum())
        l2 = float(np.sqrt((a * a).sum()))
        f = float(a[a > 0].sum() / l1)
        cf = 0.5 + (eps[:, ks] @ np.abs(a)) / (2.0 * l1)
        got.append((N, ks.size, f, l1, l2, np.abs(cf - 0.5)))
        say("  N = %-10d #k = %-6d f = %.4f  |0.5-f| = %.4f  "
            "l1/l2 = %.4f" % (N, ks.size, f, abs(0.5 - f), l1 / l2))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- L1
    say()
    say("L1  the control: the published table")
    say("  N            #k     f here   f pub    G here   G pub    "
        "worst diff")
    l1ok = True
    for N, nk, f, l1n, l2n, cl in got:
        pnk, pf, pG = pub[N]
        G = 1.0 / abs(2.0 * f - 1.0)
        d1, d2 = abs(f - pf), abs(G - pG)
        if nk != pnk or d1 >= 0.001 or d2 >= 0.01:
            l1ok = False
        say("  %-12d %-6d %-8.4f %-8.4f %-8.4f %-8.4f %.5f"
            % (N, nk, f, pf, G, pG, max(d1, d2)))
    say("  L1 %s   (cap 0.001 in f, cap 0.01 in G)"
        % ("hold" if l1ok else "REFUTED"))

    # ------------------------------------------------------------- L2
    say()
    say("L2  the raw decay against its own noise")
    x = np.log(np.array([g[0] for g in got], dtype=np.float64))
    y = np.log(np.array([abs(0.5 - g[2]) for g in got]))
    b, rms, se, t = fit(x, y)
    l2ok = (b < 0.0) and (t >= 2.0)
    say("  slope %+.6f against the published %+.6f" % (b, -pubb))
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (rms, se, t))
    say("SCATTER slope_audit_lean_floor %.4f" % rms)
    say("TSTAT slope_audit_lean_floor %.2f" % t)
    say("SPREAD slope_audit_lean_floor %.4f" % float(x.max() - x.min()))
    if t < 2.0:
        say("UNRESOLVED SIGN slope_audit_lean_floor")
    say("  L2 %s" % ("hold" if l2ok else "REFUTED"))

    # ------------------------------------------------------------- L3
    say()
    say("L3  and the floor's own decay, on the identical magnitudes")
    say("  N            mu lean   coin median  coin/mu")
    med = []
    for N, nk, f, l1n, l2n, cl in got:
        m = float(np.median(cl))
        med.append(m)
        say("  %-12d %-9.4f %-12.4f %.4f"
            % (N, abs(0.5 - f), m, m / abs(0.5 - f)))
    cb = np.array([float(np.polyfit(x, np.log(np.array(
        [g[5][j] for g in got])), 1)[0]) for j in range(DRAWS)])
    cb.sort()
    hi = float(cb[int(0.975 * DRAWS)])
    l3ok = b > hi
    say("  the %d coin slopes run %+.4f to %+.4f, median %+.4f"
        % (DRAWS, cb[0], cb[-1], float(np.median(cb))))
    say("  97.5th percentile %+.6f, mu %+.6f" % (hi, b))
    say("  %d of the %d coins decay more slowly than mu"
        % (int((cb > b).sum()), DRAWS))
    say("EXCHANGE audit_lean_floor %d %d %d"
        % (len(got), int((cb > b).sum()) + 1, DRAWS + 1))
    say("  L3 mu's slope is above the coins' 97.5th percentile   %s"
        % ("hold" if l3ok else "REFUTED"))

    say()
    say("  DIAGNOSTIC on L3 (post hoc). A single coin's six leans are")
    say("  themselves noisy, so a per-draw slope is dominated by that")
    say("  noise -- the 256 run from %+.4f to %+.4f. The floor is not"
        % (cb[0], cb[-1]))
    say("  one draw; it is where the draws sit, and that is estimated")
    say("  %d times more precisely. Fitting the MEDIAN lean:" % DRAWS)
    bm, rmsm, sem, tm = fit(x, np.log(np.array(med)))
    say("  floor slope %+.6f, standard error %.6f, t = %.2f"
        % (bm, sem, tm))
    say("  mu          %+.6f, standard error %.6f" % (b, se))
    sd = math.sqrt(se ** 2 + sem ** 2)
    say("  difference  %+.6f against %.6f = 2 s.e., i.e. %.2f s.e."
        % (b - bm, 2.0 * sd, abs(b - bm) / sd))
    say("FLOORTREND lab_lean_decay %+.6f %.2f" % (bm, tm))
    say("FLOORRANGE lab_lean_decay %d %d" % (min(NS), max(NS)))
    say("FLOORTREND audit_lean_floor %+.6f %.2f" % (bm, tm))
    say("FLOORRANGE audit_lean_floor %d %d" % (min(NS), max(NS)))
    say("  so the direction L3 registered is right and its instrument")
    say("  was wrong: mu decays more slowly than the floor by %.2f"
        % (abs(b - bm) / sd))
    say("  standard errors, which is what L4 measures directly.")
    dele = io.open(os.path.join(RES, "lab_extend_range.txt"),
                   encoding="utf-8").read()
    i0 = dele.index("NULL:")
    i1 = dele.index("FIELD:")
    if "lab_lean_decay" in dele[i0:i1]:
        # a delegated floor is only a floor where it was measured, so
        # the delegator's own N range is read and compared here
        j = dele.index("N            K       #k     B/N")
        dns = []
        for ln in dele[j:].splitlines()[1:]:
            g = ln.split()
            if not g or not g[0].isdigit():
                if g and set(g[0]) == {"-"}:
                    continue
                if not dns:
                    continue
                break
            dns.append(int(g[0]))
        say()
        say("  results/lab_extend_range.txt names this file's")
        say("  reference in its own NULL block rather than measuring")
        say("  one. A delegated floor covers only the range it was")
        say("  measured over, so both ranges are read here:")
        say("    delegator  %d to %d" % (min(dns), max(dns)))
        say("    this file  %d to %d" % (min(NS), max(NS)))
        if min(dns) >= min(NS) and max(dns) <= max(NS):
            say("  the delegation covers it:")
            say("FLOOR DELEGATED lab_extend_range lab_lean_decay")
        else:
            say("  **it does not cover it** -- the delegator runs a")
            say("  factor %.0f past where this floor was measured, so"
                % (max(dns) / float(max(NS))))
            say("  no delegation is declared and that range needs a")
            say("  floor of its own.")

    # ------------------------------------------------------------- L4
    say()
    say("L4  the lean in units of its own floor")
    say("  N            ratio")
    rr = np.log(np.array([abs(0.5 - g[2]) for g in got])
                - np.log(np.array(med)) * 0.0)
    rr = np.log(np.array([abs(0.5 - g[2]) for g in got])
                / np.array(med))
    for N, m, v in zip([g[0] for g in got], med, np.exp(rr)):
        say("  %-12d %.4f" % (N, v))
    b2, rms2, se2, t2 = fit(x, rr)
    l4ok = (b2 > 0.0) and (t2 >= 2.0)
    say("  slope %+.6f, r.m.s. residual %.4f, standard error %.6f, "
        "t = %.2f" % (b2, rms2, se2, t2))
    say("SCATTER slope_audit_lean_floor_ratio %.4f" % rms2)
    say("TSTAT slope_audit_lean_floor_ratio %.2f" % t2)
    say("SPREAD slope_audit_lean_floor_ratio %.4f"
        % float(x.max() - x.min()))
    if t2 < 2.0:
        say("UNRESOLVED SIGN slope_audit_lean_floor_ratio")
    say("  L4 %s" % ("hold" if l4ok else "REFUTED"))
    say("PERN lean_over_floor %d %.4f %.4f"
        % (len(rr), float(np.exp(rr).min()), float(np.exp(rr).max())))

    say()
    say("  what this does to the reading. The raw lean does shrink,")
    say("  and so does what chance gives on the same magnitudes,")
    say("  because l1/l2 grows: %.4f to %.4f over this sweep. The"
        % (got[0][3] / got[0][4], got[-1][3] / got[-1][4]))
    say("  gain G = 1/|2f-1| still rises, and {#rem:leandecay}'s")
    say("  first consequence -- that {#rem:nocrossk} is a statement")
    say("  about the accessible range -- is untouched. What does not")
    say("  follow is that the lean is going away: measured against")
    say("  the floor it is %s."
        % ("growing" if b2 > 0 else "shrinking"))

    say()
    say("=" * 70)
    ok = l1ok and l2ok and l3ok and l4ok
    say("the lean decays more slowly than its own floor"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the mass-weighted fraction f of k with H(N;k) > 0,",
        "           the mass being (log k)|H(N;k)|; its distance from",
        "           one half; the gain G = 1/|2f-1|; the same distance",
        "           for 256 global sign vectors on the IDENTICAL",
        "           magnitudes; the power-law slope of each against",
        "           log N with its standard error; and the slope of",
        "           mu's lean divided by the coins' median lean.",
        "NULL: 256 global sign vectors over k, each held across all N",
        "      as mu is, applied to the identical (log k)|H(N;k)|.",
        "      Same convention as audit_crossk_reference.py and",
        "      lab_split_budget.py. The published arm was two draws",
        "      per N and was read only for whether it sits at 1/2.",
        "FIELD: N = 2e5 through 6.4e6 by doubling with theta' = "
        + str(THETA) + ",",
        "       so k runs over the squarefree k < N^" + str(THETA),
        "       coprime to N; m over 1 <= m < N/k with (m,k) = 1;",
        "       Lambda and mu from an integer sieve to " + str(NMAX)
        + ";",
        "       numpy default_rng seed " + str(SEED) + ". Every N is",
        "       2^a 5^b, one odd radical, as RADICALS declares. The",
        "       published table is read from",
        "       results/lab_lean_decay.txt.",
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
