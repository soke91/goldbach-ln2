# -*- coding: utf-8 -*-
r"""
Three exponents, one identity: what would make the lean go away?

WHAT IS AT STAKE

OPEN item 5 is that mu's sign lean grows relative to its own floor --
8.49 to 21.36 over a factor 128, at +0.159294 per unit log N
({#rem:leanextended}) -- and that nothing supports the lean vanishing
in the limit. That is stated as a measurement with no mechanism, and
it does not have to be. Writing a_k = (log k)H(N;k), all three
quantities in play are ratios of the same three norms:

    G          = l1 / |sum a|            ({#rem:nocrossk})
    ceiling    = l1 / l2                 ({#rem:crosskreference})
    lean/floor = |sum a| / (c * l2)

with c the constant a median sign sum sits on. So

    lean/floor = (l1/l2) / G / c

identically, and the three exponents must satisfy
e(lean/floor) = e(l1/l2) - e(G). The lean grows against its floor for
exactly one reason: the magnitude concentration l1/l2 grows faster
than the cancellation G does.

That turns item 5 into a question with a number attached. l1/l2 is
bounded by sqrt(#k) and #k is of order N^theta' on this k-range, so
e(l1/l2) cannot exceed theta'/2 = 0.28. If it sits AT 0.28 then the
concentration is asymptotically a fixed fraction of its own ceiling
and the only free exponent is G's -- and the lean goes away only if
e(G) rises from its measured value to 0.28.

BACKS: Remark {#rem:leanidentity} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The control: G reproduces {#rem:leandecay}'s 1.834, 1.804,
      2.207, 2.588, 2.789, 3.079 to within 0.01 at the six N it
      publishes.
  W2  The identity is exact: log(lean/floor) - log(l1/l2) + log G is
      constant across N to within 0.02.
  W3  So are the exponents: e(l1/l2) - e(G) - e(lean/floor) is within
      two standard errors of zero.
  W4  And the concentration runs at its ceiling rate: e(l1/l2) is
      within two standard errors of theta'/2 = 0.28.

REFUTATION RULE (fixed before the run)

  W1  REFUTED at 0.01 -- not the same statistic, and nothing below
      may be compared with {#rem:leandecay} or {#rem:nocrossk}.
  W2  REFUTED beyond 0.02. The three quantities would not be the
      ratios this remark claims and the identity would be wrong.
  W3  REFUTED beyond two standard errors, the same failure read
      through the fits rather than pointwise.
  W4  REFUTED beyond two standard errors of theta'/2. That is the one
      that matters: the concentration would then be losing ground
      against its own ceiling, and the lean's growth would have a
      second source besides G's shortfall -- item 5 would not reduce
      to one exponent.

  All four gate.

  NO NULL IS RUN for the identity, which is algebra checked
  numerically. The floor itself is the 256 global sign vectors of
  audit_lean_floor.py on the identical magnitudes, and the constant c
  is measured here rather than assumed.
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
OUT = os.path.join(RES, "audit_lean_identity.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
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
    """the published gain G at each N"""
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
        out[int(f[0])] = float(f[4])
    return out


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

    pub = read_published()
    say("read %d published gains from results/lab_lean_decay.txt"
        % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    KMAX = int(NMAX ** THETA) + 1
    rng = np.random.default_rng(SEED)
    eps = (rng.integers(0, 2, size=(DRAWS, KMAX + 1))
           .astype(np.int8) * 2 - 1)
    say("%d global sign vectors over k <= %d, seed %d"
        % (DRAWS, KMAX, SEED))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    rows = []
    kss, aas = {}, {}
    for N in NS:
        ks, a = weighted(N, lam, mu, sqf)
        kss[N], aas[N] = ks, a
        l1 = float(np.abs(a).sum())
        l2 = float(np.sqrt((a * a).sum()))
        sa = abs(float(a.sum()))
        w = np.abs(a)
        fl = float(np.median(np.abs(eps[:, ks] @ w)))
        rows.append((N, ks.size, l1, l2, sa, fl))
        say("  N = %-10d #k = %-6d G = %.4f  l1/l2 = %.4f  "
            "lean/floor = %.4f"
            % (N, ks.size, l1 / sa, l1 / l2, sa / fl))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))

    # ------------------------------------------------------------- W1
    say()
    say("W1  the control: the published gain")
    say("  N            G here     published  diff")
    w1 = True
    for N, nk, l1, l2, sa, fl in rows:
        if N not in pub:
            continue
        d = abs(l1 / sa - pub[N])
        if not (d < 0.01):
            w1 = False
        say("  %-12d %-10.4f %-10.4f %.5f" % (N, l1 / sa, pub[N], d))
    say("  W1 %s   (cap 0.01)" % ("hold" if w1 else "REFUTED"))

    # ------------------------------------------------------------- W2
    say()
    say("W2  the identity, pointwise")
    say("  lean/floor = (l1/l2) / G / c,  so log c is the residue")
    say("  N            log(lean/floor)  -log(l1/l2)  +log G   log c")
    cs = []
    for N, nk, l1, l2, sa, fl in rows:
        lf = math.log(sa / fl)
        lc = math.log(l1 / l2)
        lg = math.log(l1 / sa)
        cs.append(lc - lg - lf)
        say("  %-12d %-16.6f %-12.6f %-8.6f %.6f"
            % (N, lf, -lc, lg, cs[-1]))
    w2 = (max(cs) - min(cs)) < 0.02
    say("  log c runs %.6f to %.6f, a spread of %.6f   (cap 0.02)"
        % (min(cs), max(cs), max(cs) - min(cs)))
    say("  c itself is %.6f to %.6f; the median of |Z| for a standard"
        % (math.exp(min(cs)), math.exp(max(cs))))
    say("  normal would put it at %.6f"
        % float(np.median(np.abs(rng.standard_normal(400000)))))
    say("  W2 %s" % ("hold" if w2 else "REFUTED"))
    if not w2:
        say()
        say("  DIAGNOSTIC on W2 (post hoc). c is a MEDIAN over %d"
            % DRAWS)
        say("  draws, so it carries its own sampling error, and the")
        say("  cap was set without asking how large that is. Splitting")
        say("  the draws into 16 groups of %d and taking each group's"
            % (DRAWS // 16))
        say("  median gives the estimator's own scatter directly:")
        say("  N            group scatter  s.e. of the full median")
        ses = []
        for (N, nk, l1, l2, sa, fl), _ in zip(rows, rows):
            ks2, a2 = None, None
            gm = []
            for g in range(16):
                sl = eps[g * (DRAWS // 16):(g + 1) * (DRAWS // 16)]
                gm.append(float(np.median(np.abs(
                    sl[:, kss[N]] @ np.abs(aas[N])))))
            gm = np.array(gm) / l2
            sd = float(gm.std(ddof=1))
            ses.append(sd / 4.0)
            say("  %-12d %-14.6f %.6f" % (N, sd, sd / 4.0))
        obs = max(cs) - min(cs)
        samp = 4.0 * float(np.mean(ses)) / math.exp(np.mean(cs))
        say("  the observed log spread is %.6f; the sampling spread of"
            % obs)
        say("  eight independent medians at that standard error is")
        say("  about %.6f, so the two are %s."
            % (samp, "the same size" if obs <= 2.0 * samp
               else "not the same size"))
        say("CONSTSPREAD leanidentity_c %.6f %.6f" % (obs, samp))
        if obs > 2.0 * samp:
            say("CONST DRIFTS leanidentity_c")
        else:
            say("  So W2 fails as registered and the failure is the")
            say("  cap: c does not drift, it is estimated to about")
            say("  %.0f per cent by %d draws." % (100 * samp, DRAWS))

    # ---------------------------------------------------------- W3/W4
    say()
    say("W3/W4  the three exponents")
    eg, rg, seg, tg = fit(x, np.log(np.array(
        [r[2] / r[4] for r in rows])))
    ec, rc, sec, tc = fit(x, np.log(np.array(
        [r[2] / r[3] for r in rows])))
    el, rl, sel, tl = fit(x, np.log(np.array(
        [r[4] / r[5] for r in rows])))
    say("  quantity        exponent     s.e.       t")
    for nm, e, se, t in (("G", eg, seg, tg),
                         ("l1/l2", ec, sec, tc),
                         ("lean/floor", el, sel, tl)):
        say("  %-15s %+-12.6f %-10.6f %.2f" % (nm, e, se, t))
        say("SCATTER slope_leanidentity_%s %.4f"
            % (nm.replace("/", "_"),
               {"G": rg, "l1/l2": rc, "lean/floor": rl}[nm]))
        say("TSTAT slope_leanidentity_%s %.2f"
            % (nm.replace("/", "_"), t))
        say("SPREAD slope_leanidentity_%s %.4f"
            % (nm.replace("/", "_"), float(x.max() - x.min())))
        if t < 2.0:
            say("UNRESOLVED SIGN slope_leanidentity_%s"
                % nm.replace("/", "_"))
    sd3 = math.sqrt(sec ** 2 + seg ** 2 + sel ** 2)
    w3 = abs(ec - eg - el) <= 2.0 * sd3
    say("  e(l1/l2) - e(G) - e(lean/floor) = %+.6f against %.6f "
        "= 2 s.e." % (ec - eg - el, 2.0 * sd3))
    say("  W3 %s" % ("hold" if w3 else "REFUTED"))
    half = THETA / 2.0
    w4 = abs(ec - half) <= 2.0 * sec
    say("  theta'/2 = %.2f, and e(l1/l2) is %+.6f away, %.2f s.e."
        % (half, ec - half, abs(ec - half) / sec))
    say("  W4 %s" % ("hold" if w4 else "REFUTED"))

    say()
    say("  what would have to change. The lean stops growing against")
    say("  its floor exactly when e(G) reaches e(l1/l2). Measured,")
    say("  e(G) = %+.6f and it would have to reach %+.6f -- a factor"
        % (eg, ec))
    say("  %.2f in the exponent. And l1/l2 is bounded by sqrt(#k):"
        % (ec / eg))
    say("  N            l1/l2      sqrt(#k)   ratio")
    for N, nk, l1, l2, sa, fl in rows:
        say("  %-12d %-10.4f %-10.4f %.4f"
            % (N, l1 / l2, math.sqrt(nk), (l1 / l2) / math.sqrt(nk)))
    say("REFERENCE audit_lean_identity %d %.4f %.4f"
        % (len(rows),
           min((r[2] / r[3]) / math.sqrt(r[1]) for r in rows),
           max((r[2] / r[3]) / math.sqrt(r[1]) for r in rows)))

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("item 5 reduces to one exponent, G's" if ok else "REFUTED")

    head = [
        "STATISTIC: with a_k = (log k)H(N;k) over the squarefree",
        "           k < N^" + str(THETA) + " coprime to N, the three",
        "           ratios G = l1/|sum a|, l1/l2 and",
        "           |sum a| / median|sum eps a|; the constant relating",
        "           them; and each one's least-squares exponent",
        "           against log N with its standard error.",
        "NULL: none for the identity, which is algebra checked",
        "      numerically. The floor is " + str(DRAWS) + " global sign",
        "      vectors on the identical magnitudes, the convention of",
        "      audit_lean_floor.py, and the constant relating the",
        "      three is measured here rather than assumed.",
        "FIELD: N = 2e5 through 2.56e7 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + "; numpy default_rng",
        "       seed " + str(SEED) + ". Every N is 2^a 5^b, one odd",
        "       radical, as RADICALS declares. The published gains are",
        "       read from results/lab_lean_decay.txt.",
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
