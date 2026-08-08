# -*- coding: utf-8 -*-
r"""
The weight both predictors drop: log(N - mk).

WHAT IS AT STAKE

OPEN item 4 is that nothing elementary carries the sign lean where the
lean is measured. On k < N^0.56 the odd Mertens function overshoots by
1.73 to 2.96 with a flat decay ({#rem:oddmertensrange}) and the
sieve-weighted P = sum_m mu(m) w(m,k) undershoots by 0.66 to 0.87 and
decays too fast ({#rem:survivorrange}). They fail in opposite
directions, which says something is missing rather than that the idea
is wrong.

What is missing is elementary and is the same thing in both. H(N;k)
weights each surviving m by Lambda(N - mk) = log p, while P counts it.
Over a short inner sum that hardly matters: N - mk stays near N and the
weight is nearly constant, so counting and weighting agree up to a
factor. Over a long one N - mk runs from about N down to about k, a
factor of seven in the logarithm at these sizes, and the terms with the
largest weight are the ones with the smallest N - mk -- exactly the
ones {#rem:survivors} showed the prime density selects.

So the natural next predictor is the one that keeps the weight the
sieve already models the density of:

    P_log(N;k) = sum_m mu(m) w(m,k) log(N - mk).

It uses no primality, only the sieve weight and the logarithm, so it
is elementary in the same sense P is. Whether it carries the lean is
the question.

BACKS: Remark {#rem:logweightpredictor} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  K1  The control: on k < N^0.56 the sign agreement of P reproduces
      {#rem:survivorrange}'s 0.8129, 0.7632, 0.7446, 0.7367, 0.7759,
      0.7579 to within 0.001, and its predicted lean the published
      ratios to within 0.001.
  K2  The log weight helps: P_log's agreement exceeds P's at every N.
  K3  And it carries the lean: P_log's predicted |0.5 - f| is within a
      factor 1.2 of mu's at every N, where P ran 0.66 to 0.87.
  K4  And its trend: P_log's slope against log N is within two
      standard errors of mu's, where P's was 2.11 out.

REFUTATION RULE (fixed before the run)

  K1  REFUTED at 0.001 anywhere -- not the same statistic, and nothing
      below may be compared with {#rem:survivorrange}.
  K2  REFUTED if the agreement is not higher at every N. The missing
      weight would then not be what separates the predictors from H on
      long inner sums, and the opposite-direction failures of
      {#rem:oddmertensrange} and {#rem:survivorrange} would need
      another explanation.
  K3  REFUTED outside a factor 1.2 at any N. That is the one that
      matters: OPEN item 4 would stay open, with the log weight ruled
      out as the missing piece.
  K4  REFUTED beyond two standard errors, the same conclusion for the
      trend rather than the level.

  All four gate.

  THE NULL is the 256 global sign vectors of audit_lean_floor.py on
  the identical magnitudes, so P_log's lean is read on the same scale
  as mu's, as P's and as Modd's. A permutation control for the sign
  agreement is not repeated here: lab_survivor_selection.py ran it for
  this predictor family and it sits at 0.5372 to 0.5414.
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
OUT = os.path.join(RES, "audit_logweight_predictor.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
QSIEVE = 30
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


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


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
    """the agreements and lean ratios of {#rem:survivorrange}"""
    src = io.open(os.path.join(RES, "audit_survivor_range.txt"),
                  encoding="utf-8").read()
    i = src.index("N            #k     agreement  mu lean   P lean"
                  "    ratio")
    agr, rat = {}, {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 6 or not f[0].isdigit():
            break
        agr[int(f[0])] = float(f[2])
        rat[int(f[0])] = float(f[5])
    return agr, rat


def scan(N, lam, mu, oddsqf, sqf, vmask, qs):
    """H, the counting predictor P and the log-weighted P_log"""
    PN = factor_set(N)
    ks, Hs, Ps, Ls = [], [], [], []
    for k in range(2, int(N ** THETA)):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = N // k
        if M < 1:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[oddsqf[ms]]
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        kb = 0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(float(g[keep].sum()))
        Ls.append(float((g[keep]
                         * np.log(vals[keep].astype(np.float64)))
                        .sum()))
    return (np.array(ks, dtype=np.int64), np.array(Hs),
            np.array(Ps), np.array(Ls))


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

    pagr, prat = read_published()
    say("read %d agreements and lean ratios from "
        "results/audit_survivor_range.txt" % len(pagr))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    vmask = residue_mask(NMAX, qs)
    rng = np.random.default_rng(SEED)
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    # ---------------------------------------------------------- K1/K2
    say()
    say("K1/K2  the two predictors on k < N^%.2f" % THETA)
    say("  N            #k     P agree  published  P_log agree  "
        "better")
    k1 = k2 = True
    xs, mus, pl, pp, flo, als = [], [], [], [], [], []
    for N in NS:
        ks, H, P, L = scan(N, lam, mu, oddsqf, sqf, vmask, qs)
        sh, sp, sl = np.sign(H), np.sign(P), np.sign(L)
        ok = (sh != 0) & (sp != 0)
        ap = float((sh[ok] == sp[ok]).mean())
        ok2 = (sh != 0) & (sl != 0)
        al = float((sh[ok2] == sl[ok2]).mean())
        if abs(ap - pagr[N]) >= 0.001:
            k1 = False
        if al <= ap:
            k2 = False
        a = np.log(ks.astype(np.float64)) * H
        l1 = float(np.abs(a).sum())
        w = np.abs(a)
        f_mu = float(a[a > 0].sum() / l1)
        spz = np.where(sp == 0, 1.0, sp)
        slz = np.where(sl == 0, 1.0, sl)
        f_p = float(w[spz > 0].sum() / l1)
        f_l = float(w[slz > 0].sum() / l1)
        eps = (rng.integers(0, 2, size=(DRAWS, ks.size))
               .astype(np.int8) * 2 - 1)
        fl = float(np.median(np.abs((eps @ w) / (2.0 * l1))))
        als.append(al)
        xs.append(math.log(N))
        mus.append(abs(0.5 - f_mu))
        pp.append(abs(0.5 - f_p))
        pl.append(abs(0.5 - f_l))
        flo.append(fl)
        say("  %-12d %-6d %-8.4f %-10.4f %-12.4f %s"
            % (N, ks.size, ap, pagr[N], al,
               "yes" if al > ap else "NO"))
    say("  K1 P reproduces the published agreement   %s   (cap 0.001)"
        % ("hold" if k1 else "REFUTED"))
    say("  K2 the log weight raises it at every N   %s"
        % ("hold" if k2 else "REFUTED"))

    # ------------------------------------------------------------- K3
    say()
    say("K3  the lean each predictor gives on mu's own magnitudes")
    say("  N            mu lean   P lean    P ratio  P_log lean  "
        "P_log ratio")
    k3 = True
    for N, lm, vp, vl, r0 in zip(NS, mus, pp, pl,
                                 [prat[n] for n in NS]):
        rp, rl = vp / lm, vl / lm
        if abs(rp - r0) >= 0.001:
            k1 = False
        if not (1.0 / 1.2 <= rl <= 1.2):
            k3 = False
        say("  %-12d %-9.4f %-9.4f %-8.4f %-11.4f %.4f"
            % (N, lm, vp, rp, vl, rl))
    say("  K1 also covers the published P ratios   %s" % ("hold" if k1
                                                          else "REFUTED"))
    say("  K3 P_log within a factor 1.2   %s"
        % ("hold" if k3 else "REFUTED"))
    say("PERN logweight_over_mu %d %.4f %.4f"
        % (len(pl), min(v / m for v, m in zip(pl, mus)),
           max(v / m for v, m in zip(pl, mus))))

    # ------------------------------------------------------------- K4
    say()
    say("K4  and the trends")
    x = np.array(xs)
    bm, rm, sem, tm = fit(x, np.log(np.array(mus)))
    bp, rp2, sep, tp = fit(x, np.log(np.array(pp)))
    bl, rl2, sel, tl = fit(x, np.log(np.array(pl)))
    sd = math.sqrt(sem ** 2 + sel ** 2)
    k4 = abs(bm - bl) <= 2.0 * sd
    say("  mu        slope %+.6f, standard error %.6f" % (bm, sem))
    say("  P         slope %+.6f, standard error %.6f" % (bp, sep))
    say("  P_log     slope %+.6f, standard error %.6f" % (bl, sel))
    say("  mu - P_log %+.6f against %.6f = 2 s.e., i.e. %.2f s.e."
        % (bm - bl, 2.0 * sd, abs(bm - bl) / sd))
    say("SCATTER slope_audit_logweight_predictor %.4f" % rl2)
    say("TSTAT slope_audit_logweight_predictor %.2f" % tl)
    say("SPREAD slope_audit_logweight_predictor %.4f"
        % float(x.max() - x.min()))
    if tl < 2.0:
        say("UNRESOLVED SIGN slope_audit_logweight_predictor")
    say("  K4 %s" % ("hold" if k4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The two scores move in opposite")
    say("  directions, so the ledger is printed with both and the")
    say("  criterion fixed. Modd's row is read from")
    say("  results/audit_oddmertens_range.txt:")
    osrc = io.open(os.path.join(RES, "audit_oddmertens_range.txt"),
                   encoding="utf-8").read()
    oagr = [float(v) for v in re.findall(
        r"^  \d+\s+\d+\s+\d+\s+([\d.]+)\s+[\d.]+\s+[\d.]+"
        r"\s+[\d.]+\s*$", osrc, re.M)]
    orat = [float(v) for _n, v in re.findall(
        r"^  (\d+)\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+"
        r"\s+([\d.]+)\s*$", osrc, re.M)]
    oslope = float(re.search(r"Modd      slope ([-+][\d.]+)",
                             osrc).group(1))
    rows = [("modd", min(oagr), max(oagr), min(orat), max(orat),
             oslope),
            ("sieve_P", min(pagr.values()), max(pagr.values()),
             min(prat.values()), max(prat.values()), bp),
            ("sieve_P_log", min(als), max(als),
             min(v / m for v, m in zip(pl, mus)),
             max(v / m for v, m in zip(pl, mus)), bl)]
    say("  predictor      agreement        lean ratio       slope")
    for nm, a0, a1, r0, r1, sl in rows:
        say("  %-14s %.4f to %.4f  %.4f to %.4f  %+.6f"
            % (nm, a0, a1, r0, r1, sl))
        say("PREDICTOR ladder_lean %s %.4f %.4f %+.6f"
            % (nm, a0, r1, sl))
    say("PREDICTOR CRITERION ladder_lean agreement")
    best = max(rows, key=lambda t: t[1])[0]
    say("PREDICTOR BEST ladder_lean %s" % best)
    say("  mu's own slope is %+.6f, so on the lean ratio the log"
        % bm)
    say("  weight is the closest of the three and on the agreement it")
    say("  is not; the criterion registered here is the agreement,")
    say("  which is the score the predictors were built for.")

    say()
    say("  and all three against the floor of {#rem:leanfloor}:")
    say("  N            mu/floor   P/floor    P_log/floor")
    for N, lm, vp, vl, fl in zip(NS, mus, pp, pl, flo):
        say("  %-12d %-10.4f %-10.4f %.4f"
            % (N, lm / fl, vp / fl, vl / fl))
    b1, _r1, s1e, _t1 = fit(x, np.log(np.array(mus) / np.array(flo)))
    b3, _r3, s3e, _t3 = fit(x, np.log(np.array(pl) / np.array(flo)))
    say("  mu %+.6f (s.e. %.6f);  P_log %+.6f (s.e. %.6f); "
        "difference %.2f s.e."
        % (b1, s1e, b3, s3e,
           abs(b1 - b3) / math.sqrt(s1e ** 2 + s3e ** 2)))

    say()
    say("=" * 70)
    ok = k1 and k2 and k3 and k4
    say("the missing weight is log(N - mk)" if ok else "REFUTED")

    head = [
        "STATISTIC: on the squarefree k < N^" + str(THETA)
        + " coprime to N,",
        "           the sign agreement of H(N;k) with the counting",
        "           predictor P = sum_m mu(m) w(m,k) and with the",
        "           log-weighted P_log = sum_m mu(m) w(m,k)",
        "           log(N - mk); the mass-weighted lean each gives on",
        "           mu's own (log k)|H| magnitudes; their slopes",
        "           against log N with standard errors; and all three",
        "           against the median lean of " + str(DRAWS),
        "           sign vectors on those magnitudes.",
        "NULL: " + str(DRAWS) + " global sign vectors on the identical",
        "      magnitudes, the convention of audit_lean_floor.py, so",
        "      the three leans are read on one scale. The permutation",
        "      control for this predictor family was run in",
        "      lab_survivor_selection.py at 0.5372 to 0.5414 and is",
        "      not repeated.",
        "FIELD: N = 2e5 through 6.4e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m odd,",
        "       squarefree and coprime to k, m <= N//k, the convention",
        "       of lab_survivor_selection.py; the sieve weight over the",
        "       odd primes below " + str(QSIEVE) + "; Lambda and mu",
        "       from an integer sieve to " + str(NMAX) + "; numpy",
        "       default_rng seed " + str(SEED) + ". Every N is 2^a 5^b,",
        "       one odd radical, as RADICALS declares. The published",
        "       agreements and ratios are read from",
        "       results/audit_survivor_range.txt.",
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
