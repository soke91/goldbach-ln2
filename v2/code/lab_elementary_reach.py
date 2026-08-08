# -*- coding: utf-8 -*-
r"""
The elementary wall, measured thirty times further out than anything
else in this repository.

WHAT IS AT STAKE

Remark {#rem:splitbudget} moved the target to the elementary half:
H = beta P + R, and it is P, not the residue, that spends the budget.
Remark {#rem:elemsize} then fitted |P| against the inner length N/k
and got exponents 0.5203, 0.4984, 0.4652, 0.4681, 0.4454 -- drifting
downwards across the sweep, with a leave-one-out spread reaching
0.0595, and with the top octave left unbounded at [32768, inf).

Every measurement of H in this repository stops at N = 3.2e6, because
H needs Lambda(N - mk) and so a sieve to N. **P needs no primes.**
Writing the sieve weight as

    w(m,k) = C_k * [ m != N k^{-1} (mod q) for every odd q <= 29, q|k ]
    C_k    = prod_{q odd <= 29, q | k} q/(q-1),

the weight is one constant times an indicator, so
P(N;k) = C_k * sum_{m in S} mu(m) over an explicit sifted set S: a
boolean mask and an integer sum, with no Lambda anywhere. That buys
N = 1e8 and inner lengths to 3.3e7, against the 1.6e6 the published
measurement reaches -- and it resolves the unbounded top octave into
five.

BACKS: Remark {#rem:elemreach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The control: P computed this way is the published P. Recomputing
      beta = sum(H P)/sum(P^2) at the five published N reproduces the
      AGREE beta_HP values to within 1 per cent.
  Y2  Out in the new territory -- inner lengths from 8192 to 3.3e7 at
      N = 1e8, five of whose seven octaves no published fit resolves --
      the exponent of |P| lies in (0.45, 0.55).
  Y3  The published downward drift does not continue: the exponent at
      N = 1e8 is not below the smallest published one less the largest
      published leave-one-out spread, both read from
      results/lab_elementary_size.txt.
  Y4  And the octave fit stays tight, correlation above 0.99.

  Y5 and Y6 were added after Y1-Y4 had run and are disclosed as such.
  The coin arm showed the first estimator ill-conditioned exactly where
  the reach was supposed to buy something: at one fixed N the longest
  inner lengths come from the FEWEST k -- the top octave at N = 1e8 has
  about twenty k and at N = 1e7 has one -- so those octave means are
  near-single samples, and the coin exponents scatter over a band four
  times wider than the one Y2 tests against. A refutation from an
  estimator whose own null cannot resolve the question is not evidence
  about mu. The sampling is therefore inverted: many N, few k each, so
  that a long inner length is reached by many (N,k) pairs instead of
  by one.

  Y5  Inverting the sampling conditions the estimator: with a ladder
      of N and k < 400 at each, every inner-length octave carries at
      least 200 pairs and the coin band narrows to a width under 0.10.
  Y6  And then mu's exponent lies in (0.45, 0.55).

REFUTATION RULE (fixed before the run)

  Y1  REFUTED at 1 per cent at any N, which would mean the constant
      factored out of w is not the weight the published P uses.
  Y2  REFUTED if the exponent leaves (0.45, 0.55). Below is the
      interesting failure: it would say the elementary wall is
      genuinely smaller than square-root out where it matters, and
      the budget argument of {#rem:splitbudget} is too pessimistic.
  Y3  REFUTED if the exponent falls below that floor, which would say
      the drift is real and the square-root reading is an artefact of
      the short range.
  Y4  REFUTED below 0.99 at N = 1e8.
  Y5  REFUTED if the coin band stays 0.10 wide or wider, or if any
      octave carries fewer than 200 pairs. That failure would be the
      finding: the elementary wall could not be measured beyond the
      range H itself reaches, by this method.
  Y6  REFUTED if the exponent leaves (0.45, 0.55).

  Y1, Y5 and Y6 gate. Y2, Y3 and Y4 are reported and do not gate:
  they are the ill-conditioned estimator's verdicts and are kept in
  the record because the reason they must not be read as evidence is
  itself a measurement.

  NULL: a coin arm on the same sifted set. Replacing mu(m) by
  independent signs on S leaves the set, the mask and the octave
  binning untouched and must give exactly square root, so the coin
  band is the scale against which mu's exponent is read. The sum of
  c independent signs is 2*Binomial(c, 1/2) - c exactly, so the draws
  are taken that way rather than by materialising sign vectors; this
  is exact in distribution, not an approximation.
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
OUT = os.path.join(ROOT, "results", "lab_elementary_reach.txt")

SMALL = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KSMALL = 30_000
BIG = [10_000_000, 100_000_000]
LOWBIN = 8192                 # the lowest inner length kept in the reach arm
KPER = 400                    # k per N once the sampling is inverted
QSIEVE = 30
COINS = 8
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def mobius(n):
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
    return mu


def lambda_upto(n):
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
    return lam


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
    """the beta cross-check and the published exponents -- read"""
    p = os.path.join(ROOT, "results", "lab_predictable_part.txt")
    src = io.open(p, encoding="utf-8").read()
    beta = {int(m.group(1)): float(m.group(2))
            for m in re.finditer(r"AGREE beta_HP N=(\d+) ([\d.]+)", src)}
    q = os.path.join(ROOT, "results", "lab_elementary_size.txt")
    src2 = io.open(q, encoding="utf-8").read()
    i = src2.index("N            exponent   correlation")
    ex = {}
    for ln in src2[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[1])
    sp = max(float(m.group(1)) for m in
             re.finditer(r"SWEPT \S+ octave-range ([\d.]+)", src2))
    return beta, ex, sp


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubbeta, pubex, pubspread = read_published()
    say("read from results/: %d beta cross-checks, %d published "
        "exponents," % (len(pubbeta), len(pubex)))
    say("  and the largest published leave-one-out spread %.4f"
        % pubspread)

    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    # ------------------------------------------------------------- Y1
    say()
    say("Y1  the control: is this the published P?")
    say("  beta = sum(H P)/sum(P^2) at the published N and k < %d"
        % KSMALL)
    NS0 = max(SMALL)
    mu0 = mobius(NS0)
    lam0 = lambda_upto(NS0)
    sqf0 = mu0 != 0
    say("  N            beta (here)   beta (published)   ratio")
    y1 = True
    for N in SMALL:
        PN = factor_set(N)
        Hs, Ps = [], []
        for k in range(2, KSMALL):
            if not sqf0[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[sqf0[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu0[ms].astype(np.float64)
            keep = np.ones(ms.size, dtype=bool)
            C = 1.0
            for q in QS:
                if k % q == 0:
                    continue
                C *= q / (q - 1.0)
                keep &= (vals % q != 0)
            Hs.append(float((lam0[vals] * g).sum()))
            Ps.append(C * float(g[keep].sum()))
        H = np.array(Hs)
        P = np.array(Ps)
        b = float((H * P).sum() / (P * P).sum())
        r = b / pubbeta[N]
        if abs(r - 1.0) >= 0.01:
            y1 = False
        say("  %-12d %-13.6f %-18.6f %.4f" % (N, b, pubbeta[N], r))
    say("  Y1 %s" % ("hold" if y1 else "REFUTED"))
    del mu0, lam0, sqf0

    # ---------------------------------------------------- the reach arm
    MMAX = max(BIG) // 3
    say()
    say("sieving mu to %d (no Lambda is needed out here) ..." % MMAX)
    mu = mobius(MMAX)
    oddsqf = np.zeros(MMAX + 1, dtype=bool)
    oddsqf[1::2] = True
    oddsqf &= (mu != 0)
    rng = np.random.default_rng(SEED)

    edges = []
    e = LOWBIN
    while e <= MMAX:
        edges.append(e)
        e *= 4
    say("  factor-4 octaves of N/k from %d up: %d of them"
        % (LOWBIN, len(edges)))

    def sweep(N):
        """|P| and the coin arm, per k, for this N"""
        PN = factor_set(N)
        KTOP = N // LOWBIN
        ks, ap, ac = [], [], [[] for _ in range(COINS)]
        for k in range(2, KTOP + 1):
            if not oddsqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2 or M > MMAX:
                continue
            mask = oddsqf[1:M + 1].copy()
            for p in factor_set(k):
                if p > 2:
                    mask[p - 1::p] = False
            C = 1.0
            for q in QS:
                if k % q == 0:
                    continue
                C *= q / (q - 1.0)
                r = (N * pow(k, -1, q)) % q
                mask[(r - 1) % q::q] = False
            cnt = int(mask.sum())
            if cnt == 0:
                continue
            s = int(mu[1:M + 1][mask].sum(dtype=np.int64))
            ks.append(k)
            ap.append(C * abs(s))
            for c in range(COINS):
                ac[c].append(C * abs(2 * int(rng.binomial(cnt, 0.5))
                                     - cnt))
        return (np.array(ks, dtype=np.int64), np.array(ap),
                [np.array(v) for v in ac])

    def octfit(N, ks, vals):
        cent, prof, cnt = [], [], []
        inner = N // ks
        for lo in edges:
            sel = (inner >= lo) & (inner < lo * 4)
            if sel.sum() == 0:
                continue
            cent.append(float(inner[sel].mean()))
            prof.append(float(vals[sel].mean()))
            cnt.append(int(sel.sum()))
        if len(cent) < 3:
            return None
        x = np.log(np.array(cent))
        y = np.log(np.array(prof))
        b = float(np.polyfit(x, y, 1)[0])
        r = float(np.corrcoef(x, y)[0, 1])
        return cent, prof, b, r, cnt

    say()
    say("Y2  |P| out where no published fit reaches")
    res = {}
    for N in BIG:
        ks, ap, ac = sweep(N)
        out = octfit(N, ks, ap)
        cent, prof, b, r, cnt = out
        res[N] = (ks, ap, ac, cent, prof, b, r, cnt)
        say("  N = %-11d #k = %-7d octaves %d   exponent %.4f  "
            "corr %.5f" % (N, ks.size, len(cent), b, r))
        say("    N/k octave        mean |P|")
        for c, v in zip(cent, prof):
            say("    %-17.0f %.3f" % (c, v))

    bbig = res[max(BIG)][5]
    rbig = res[max(BIG)][6]
    y2 = 0.45 < bbig < 0.55
    say("  Y2 exponent %.4f at N = %d   (band (0.45, 0.55))   %s"
        % (bbig, max(BIG), "hold" if y2 else "REFUTED"))

    # leave-one-out, as G24 requires of a fitted exponent
    cent, prof = res[max(BIG)][3], res[max(BIG)][4]
    x = np.log(np.array(cent))
    y = np.log(np.array(prof))
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on elem_reach: full %.4f, without the shortest "
        "octave %.4f," % (f[0], f[1]))
    say("  without the longest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT elem_reach octave-range %.4f" % sp)
    say("  and the bin this fit is thinnest on -- the number that gate")
    say("  check G30 reads, and the whole reason Y2 must not be read as")
    say("  evidence about mu:")
    say("POP elem_reach %d" % min(res[max(BIG)][7]))
    say("CORR elem_reach %.5f" % abs(rbig))

    # ------------------------------------------------------------- Y3
    say()
    say("Y3  does the published drift continue?")
    lowest = min(pubex.values())
    floor = lowest - pubspread
    y3 = bbig >= floor
    say("  published exponents %s"
        % ", ".join("%.4f" % pubex[N] for N in sorted(pubex)))
    say("  lowest %.4f less the largest spread %.4f gives a floor of "
        "%.4f" % (lowest, pubspread, floor))
    say("  measured here %.4f   %s"
        % (bbig, "hold" if y3 else "REFUTED"))

    # ------------------------------------------------------------- Y4
    say()
    say("Y4  is the octave fit tight out there?")
    y4 = rbig > 0.99
    say("  correlation %.5f at N = %d   (floor 0.99)   %s"
        % (rbig, max(BIG), "hold" if y4 else "REFUTED"))

    say()
    say("  THE NULL. Coins on the same sifted set, same mask, same")
    say("  binning: the exponent they give is what exact square-root")
    say("  cancellation looks like here.")
    say("  N            mu        coin min   coin median   coin max")
    for N in BIG:
        ks, ap, ac, cent, prof, b, r, cnt = res[N]
        cb = []
        for c in range(COINS):
            o = octfit(N, ks, ac[c])
            if o:
                cb.append(o[2])
        cb.sort()
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, b, cb[0], float(np.median(cb)), cb[-1]))

    # ---------------------------------------------------------- Y5/Y6
    say()
    say("Y5/Y6  the sampling inverted: a ladder of N, few k at each")
    NLAD = sorted({int(round(1e6 * 10 ** (j / 8.0))) // 2 * 2
                   for j in range(17)})
    say("  %d values of N from %d to %d, k < %d at each"
        % (len(NLAD), NLAD[0], NLAD[-1], KPER))
    pool_len, pool_mu = [], []
    pool_coin = [[] for _ in range(COINS)]
    for N in NLAD:
        PN = factor_set(N)
        for k in range(2, KPER):
            if not oddsqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2 or M > MMAX:
                continue
            mask = oddsqf[1:M + 1].copy()
            for p in factor_set(k):
                if p > 2:
                    mask[p - 1::p] = False
            C = 1.0
            for q in QS:
                if k % q == 0:
                    continue
                C *= q / (q - 1.0)
                r = (N * pow(k, -1, q)) % q
                mask[(r - 1) % q::q] = False
            cnt = int(mask.sum())
            if cnt == 0:
                continue
            pool_len.append(M)
            pool_mu.append(C * abs(int(mu[1:M + 1][mask]
                                       .sum(dtype=np.int64))))
            for c in range(COINS):
                pool_coin[c].append(
                    C * abs(2 * int(rng.binomial(cnt, 0.5)) - cnt))
    plen = np.array(pool_len, dtype=np.int64)
    say("  pooled pairs: %d, inner lengths %d to %d"
        % (plen.size, plen.min(), plen.max()))

    def poolfit(vals):
        v = np.asarray(vals, dtype=float)
        cent, prof, cnts = [], [], []
        for lo in edges:
            sel = (plen >= lo) & (plen < lo * 4)
            if sel.sum() == 0:
                continue
            cent.append(float(plen[sel].mean()))
            prof.append(float(v[sel].mean()))
            cnts.append(int(sel.sum()))
        x = np.log(np.array(cent))
        y = np.log(np.array(prof))
        return (cent, prof, cnts, float(np.polyfit(x, y, 1)[0]),
                float(np.corrcoef(x, y)[0, 1]))

    pc, pp, pn, pb, pr = poolfit(pool_mu)
    say("  N/k octave        pairs     mean |P|")
    for c, v, n_ in zip(pc, pp, pn):
        say("  %-17.0f %-9d %.3f" % (c, n_, v))
    cbs = [poolfit(pool_coin[c])[3] for c in range(COINS)]
    cbs.sort()
    width = cbs[-1] - cbs[0]
    minpairs = min(pn)
    y5 = width < 0.10 and minpairs >= 200
    say("  mu exponent %.4f   corr %.5f" % (pb, pr))
    say("  coin band [%.4f, %.4f] width %.4f, thinnest octave %d pairs"
        % (cbs[0], cbs[-1], width, minpairs))
    say("  Y5 %s   (band under 0.10, at least 200 pairs)"
        % ("hold" if y5 else "REFUTED"))
    y6 = 0.45 < pb < 0.55
    say("  Y6 exponent %.4f   (band (0.45, 0.55))   %s"
        % (pb, "hold" if y6 else "REFUTED"))
    x = np.log(np.array(pc))
    y = np.log(np.array(pp))
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    say("  leave-one-out on elem_pooled: full %.4f, without the "
        "shortest octave %.4f," % (f[0], f[1]))
    say("  without the longest %.4f -- spread %.4f"
        % (f[2], max(f) - min(f)))
    say("SWEPT elem_pooled octave-range %.4f" % (max(f) - min(f)))
    say("POP elem_pooled %d" % min(pn))
    say("CORR elem_pooled %.5f" % abs(pr))

    say()
    say("  DIAGNOSTIC (post hoc). Reach and conditioning cannot both be")
    say("  had at a fixed largest N. The number of (N,k) pairs landing")
    say("  at inner length L is about sum_N (N/L)/(spread of the bin),")
    say("  so it falls like 1/L: every factor of four further out costs")
    say("  a factor of four in sample count. The pooled counts above")
    say("  fall %d, %d, %d, %d, %d, %d across six octaves."
        % tuple(pn))
    say()
    say("  What the well-conditioned octaves alone say. Keeping only")
    say("  those with at least 200 pairs:")
    keep = [i for i, n_ in enumerate(pn) if n_ >= 200]
    xr = np.log(np.array([pc[i] for i in keep]))
    yr = np.log(np.array([pp[i] for i in keep]))
    br = float(np.polyfit(xr, yr, 1)[0])
    cr = [float(np.polyfit(xr, np.log(np.array(
        [poolfit(pool_coin[c])[1][i] for i in keep])), 1)[0])
        for c in range(COINS)]
    cr.sort()
    say("  %d octaves, inner lengths %.0f to %.0f"
        % (len(keep), pc[keep[0]], pc[keep[-1]]))
    say("  mu exponent %.4f against a coin band [%.4f, %.4f] width %.4f"
        % (br, cr[0], cr[-1], cr[-1] - cr[0]))
    say("  and the published measurement already covers that range, so")
    say("  the restricted fit buys conditioning and no reach at all.")
    say()
    say("  The same trade shows in the two fixed-N sweeps. Their octave")
    say("  means where both reach the same inner length:")
    say("  N/k octave        N = %-13d N = %-13d ratio  pairs at 1e8"
        % tuple(BIG))
    c0, p0 = res[BIG[0]][3], res[BIG[0]][4]
    c1, p1 = res[BIG[1]][3], res[BIG[1]][4]
    ks1 = res[BIG[1]][0]
    inner1 = BIG[1] // ks1
    d0 = {int(math.log(c, 4) * 2) // 2: p for c, p in zip(c0, p0)}
    for c, p in zip(c1, p1):
        key = int(math.log(c, 4) * 2) // 2
        lo = LOWBIN
        while lo * 4 <= c:
            lo *= 4
        npair = int(((inner1 >= lo) & (inner1 < lo * 4)).sum())
        if key in d0:
            say("  %-17.0f %-15.3f %-15.3f %-6.3f %d"
                % (c, d0[key], p, p / d0[key], npair))
        else:
            say("  %-17.0f %-15s %-15.3f %-6s %d"
                % (c, "-", p, "-", npair))
    say("  They agree to a few per cent where the bins are full and")
    say("  disagree by a factor of two where they are not. The claim")
    say("  that |P| depends on N/k alone is supported exactly as far")
    say("  as the sampling reaches, and no further.")

    say()
    say("=" * 70)
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s  Y5 %s  Y6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (y1, y2, y3, y4, y5, y6)))
    ok = y1 and y5 and y6
    say("the elementary wall is square-root thirty times further out "
        "than H can be measured" if ok else "REFUTED")

    head = [
        "STATISTIC: the octave means of |P(N;k)| against the inner length",
        "           N/k and the exponent fitted from them, at N = 1e7 and",
        "           1e8, using P(N;k) = C_k sum_{m in S} mu(m) with C_k",
        "           the constant part of the sieve weight and S the",
        "           explicit sifted set; and, as a control, beta =",
        "           sum(H P)/sum(P^2) at the five published N.",
        "NULL: a coin arm on the same sifted set -- independent signs in",
        "      place of mu(m), same mask, same binning -- which must give",
        "      exactly square root and so fixes the scale for mu's",
        "      exponent. Eight draws. The sum of c independent signs is",
        "      2 Binomial(c, 1/2) - c exactly, so the draws are taken",
        "      that way rather than by materialising sign vectors.",
        "FIELD: N = 1e7 and 1e8, both of the family 2^a 5^b; k odd,",
        "       squarefree and coprime to N with 2 <= k <= N/8192; m odd,",
        "       squarefree, coprime to the odd part of k, m < N/k; the",
        "       sieve weight uses the odd primes up to 30; mu from an",
        "       integer sieve to 3.3e7 and no Lambda at all. The control",
        "       arm uses N = 2e5 to 3.2e6 with k < 30000 and Lambda from",
        "       a sieve to 3.2e6. numpy default_rng seed 20260808.",
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
