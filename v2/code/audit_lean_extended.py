# -*- coding: utf-8 -*-
r"""
Does the lean keep growing against its floor over a factor 128?

WHAT IS AT STAKE

Remark {#rem:leanfloor} found that the sign lean's decay is the floor
moving: over N = 2e5 to 6.4e6 the median lean of 256 sign vectors on
the identical magnitudes falls at -0.315933 against mu's -0.167257,
and mu's lean divided by that floor RISES from 9.86 to 17.14. That is
six points, and this project's own record is that six points are
where a trend can be a short-sweep artefact -- Remark {#rem:slopes}
withdrew a reading on five and Remark {#rem:slopereach} settled it
only two octaves further out.

The longer lever already exists. lab_extend_range.py carries the same
f to N = 2.56e7, eight octaves, a factor 128, and publishes
f = 0.2273 ... 0.3608. Its own NULL block declines to run a control
on the grounds that "the coin reference level for f is 1/2 and was
measured there" -- in lab_lean_decay.py, which stops at 6.4e6. So the
floor it delegates to does not cover the range it uses, and the two
largest N in these papers have never had one.

The implementation is the one of audit_lean_floor.py extended, and
E1 recomputes lab_extend_range.py's table independently.

BACKS: Remark {#rem:leanextended} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  The control: the eight mass fractions reproduce the published
      0.2273 ... 0.3608 to within 0.001 with the same #k at every N.
  E2  And the published power exponent: the raw decay slope over the
      eight reproduces -0.1539 to within 0.001.
  E3  The floor still falls faster over the whole factor 128: the
      coin-median slope is below mu's by more than two standard
      errors of the difference.
  E4  And the floor-relative lean is still rising: the slope of
      log(mu's lean / the coins' median lean) against log N is
      positive and reaches two standard errors over all eight.

REFUTATION RULE (fixed before the run)

  E1  REFUTED at 0.001 or on any #k -- not the same statistic, and
      nothing below may be compared with {#rem:leandecay} or
      lab_extend_range.py.
  E2  REFUTED at 0.001, likewise.
  E3  REFUTED below two standard errors. {#rem:leanfloor}'s reversal
      would then be a six-point artefact that the longer lever does
      not support, and {#rem:leandecay}'s reading would be restored.
  E4  REFUTED if the floor-relative slope fails to reach two standard
      errors or turns negative -- the same conclusion by the direct
      route, and the one that matters, since it is the statement
      {#rem:leanfloor} rests on.

  All four gate.

  THE NULL IS THE POINT and it is run: 256 global sign vectors over k,
  each held across all N as mu is, applied to the IDENTICAL
  (log k)|H(N;k)| magnitudes. Same convention and seed as
  audit_lean_floor.py, so the two are comparable.
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
OUT = os.path.join(RES, "audit_lean_extended.txt")

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
    """the eight-N table and the published power exponent"""
    src = io.open(os.path.join(RES, "lab_extend_range.txt"),
                  encoding="utf-8").read()
    i = src.index("N            K       #k     B/N       |E_3|/N   "
                  "f         |0.5-f|   G")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 8 or not f[0].isdigit():
            if f and set(f[0]) == {"-"}:
                continue
            if not out:
                continue
            break
        out[int(f[0])] = (int(f[2]), float(f[5]))
    b = float(re.search(r"power N\^\{-([\d.]+)\}", src).group(1))
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
    say("read %d published rows and the power exponent %.4f from "
        "results/lab_extend_range.txt" % (len(pub), pubb))

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
        cf = np.abs((eps[:, ks] @ np.abs(a)) / (2.0 * l1))
        got.append((N, ks.size, f, l1, l2, cf))
        say("  N = %-10d #k = %-6d f = %.4f  |0.5-f| = %.4f  "
            "l1/l2 = %.4f" % (N, ks.size, f, abs(0.5 - f), l1 / l2))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- E1
    say()
    say("E1  the control: the published eight-N table")
    say("  N            #k here  #k pub   f here   f pub    diff")
    e1 = True
    for N, nk, f, l1n, l2n, cf in got:
        pnk, pf = pub[N]
        d = abs(f - pf)
        if nk != pnk or d >= 0.001:
            e1 = False
        say("  %-12d %-8d %-8d %-8.4f %-8.4f %.5f"
            % (N, nk, pnk, f, pf, d))
    say("  E1 %s   (cap 0.001)" % ("hold" if e1 else "REFUTED"))

    # ------------------------------------------------------------- E2
    say()
    say("E2  the published power exponent")
    x = np.log(np.array([g[0] for g in got], dtype=np.float64))
    y = np.log(np.array([abs(0.5 - g[2]) for g in got]))
    b, rms, se, t = fit(x, y)
    e2 = abs(b + pubb) < 0.001
    say("  slope %+.6f against the published %+.6f, diff %.6f"
        % (b, -pubb, abs(b + pubb)))
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (rms, se, t))
    say("SCATTER slope_audit_lean_extended %.4f" % rms)
    say("TSTAT slope_audit_lean_extended %.2f" % t)
    say("SPREAD slope_audit_lean_extended %.4f"
        % float(x.max() - x.min()))
    if t < 2.0:
        say("UNRESOLVED SIGN slope_audit_lean_extended")
    say("  E2 %s   (cap 0.001)" % ("hold" if e2 else "REFUTED"))

    # ------------------------------------------------------------- E3
    say()
    say("E3  the floor, over the whole factor 128")
    say("  N            mu lean   coin median  ratio")
    med, rat = [], []
    for N, nk, f, l1n, l2n, cf in got:
        m = float(np.median(cf))
        med.append(m)
        rat.append(abs(0.5 - f) / m)
        say("  %-12d %-9.4f %-12.4f %.4f" % (N, abs(0.5 - f), m,
                                             abs(0.5 - f) / m))
    bm, rmsm, sem, tm = fit(x, np.log(np.array(med)))
    sd = math.sqrt(se ** 2 + sem ** 2)
    e3 = (b - bm) > 2.0 * sd
    say("  floor slope %+.6f, standard error %.6f, t = %.2f"
        % (bm, sem, tm))
    say("  mu          %+.6f, standard error %.6f" % (b, se))
    say("  difference  %+.6f against %.6f = 2 s.e., i.e. %.2f s.e."
        % (b - bm, 2.0 * sd, abs(b - bm) / sd))
    say("FLOORTREND lab_extend_range %+.6f %.2f" % (bm, tm))
    say("FLOORRANGE lab_extend_range %d %d" % (min(NS), max(NS)))
    say("FLOORTREND audit_lean_extended %+.6f %.2f" % (bm, tm))
    say("FLOORRANGE audit_lean_extended %d %d" % (min(NS), max(NS)))
    say("  E3 %s" % ("hold" if e3 else "REFUTED"))

    # ------------------------------------------------------------- E4
    say()
    say("E4  the lean in units of its floor, over all eight")
    rr = np.log(np.array(rat))
    b2, rms2, se2, t2 = fit(x, rr)
    e4 = (b2 > 0.0) and (t2 >= 2.0)
    say("  ratio by N: " + ", ".join("%.2f" % v for v in rat))
    say("  slope %+.6f, r.m.s. residual %.4f, standard error %.6f, "
        "t = %.2f" % (b2, rms2, se2, t2))
    say("SCATTER slope_audit_lean_extended_ratio %.4f" % rms2)
    say("TSTAT slope_audit_lean_extended_ratio %.2f" % t2)
    say("SPREAD slope_audit_lean_extended_ratio %.4f"
        % float(x.max() - x.min()))
    if t2 < 2.0:
        say("UNRESOLVED SIGN slope_audit_lean_extended_ratio")
    say("PERN lean_over_floor_extended %d %.4f %.4f"
        % (len(rat), min(rat), max(rat)))
    say("  E4 %s" % ("hold" if e4 else "REFUTED"))

    say()
    say("  and against the six-point version, so the two are")
    say("  comparable rather than merely consistent:")
    six = [i for i, g in enumerate(got) if g[0] <= 6_400_000]
    b6, _r6, se6, _t6 = fit(x[six], rr[six])
    say("  six points  %+.6f (s.e. %.6f)" % (b6, se6))
    say("  eight       %+.6f (s.e. %.6f)" % (b2, se2))
    say("  the longer lever %s the slope by %+.6f, which is %.2f"
        % ("raises" if b2 > b6 else "lowers", b2 - b6,
           abs(b2 - b6) / math.sqrt(se6 ** 2 + se2 ** 2)))
    say("  standard errors of the difference.")

    say()
    say("=" * 70)
    ok = e1 and e2 and e3 and e4
    say("the lean grows against its floor over a factor 128"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the mass-weighted fraction f of k with H(N;k) > 0,",
        "           the mass being (log k)|H(N;k)|; its distance from",
        "           one half; the same distance for 256 global sign",
        "           vectors on the IDENTICAL magnitudes; the",
        "           power-law slope of each against log N with its",
        "           standard error; and the slope of mu's lean divided",
        "           by the coins' median lean, over eight N.",
        "NULL: 256 global sign vectors over k, each held across all N",
        "      as mu is, applied to the identical (log k)|H(N;k)|.",
        "      Same convention and seed as audit_lean_floor.py.",
        "FIELD: N = 2e5 through 2.56e7 by doubling with theta' = "
        + str(THETA) + ",",
        "       so k runs over the squarefree k < N^" + str(THETA),
        "       coprime to N; m over 1 <= m < N/k with (m,k) = 1;",
        "       Lambda and mu from an integer sieve to " + str(NMAX)
        + ";",
        "       numpy default_rng seed " + str(SEED) + ". Every N is",
        "       2^a 5^b, one odd radical, as RADICALS declares. The",
        "       published table is read from",
        "       results/lab_extend_range.txt.",
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
