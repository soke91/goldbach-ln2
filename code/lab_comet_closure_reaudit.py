# -*- coding: utf-8 -*-
"""
Re-auditing a recorded CLOSURE: is the ln 2 comet really not measuring
C(N)? (increment 282)

WHAT WAS CLOSED. MEASUREMENTS section 9 refuted the hypothesis that
the ln 2 comet corpus of section 8 and the chain's final scalar C(N)
are the same object. Its pre-registered criterion was

    corr(C, D) >= 0.9 required;  0.60 pooled  ->  REFUTED

with D(N) := N - r(N)/S(N) the relative Goldbach deficit, and the
supporting narrative was that |R(N)| ~ N^0.599 grows faster than
|C(N)| ~ N^0.503, so D = C - R/S is dominated by R rather than by C.

WHY IT IS WORTH RE-OPENING. Three reasons, all visible in the record
without any new computation.

 (1) THE PER-GROUP CORRELATIONS RISE WITH N: 0.457, 0.746, 0.771,
     0.831, 0.726. A statistic compared against a fixed threshold
     while it is still climbing is hazard 5 -- a value below threshold
     at small N is not a refutation if the value is increasing.
 (2) THE SUPPORTING EXPONENTS COME FROM THE DESIGN INCREMENT 281 JUST
     REFUTED: five groups of 80 consecutive even N, replication sd
     0.043 on the exponent. The gap 0.599 - 0.503 = 0.096 needs to be
     judged against the sd of the DIFFERENCE, which nobody computed.
 (3) BOTH C AND D CONTAIN THE LOCATION MASK (found at increment 240,
     four increments after this measurement). A deterministic term
     common to both INFLATES their correlation. If the criterion was
     met partly by shared mask, the true correlation is lower and the
     closure is safer than it looks; if the mask sits in only one of
     them, it is not. Either way it was never checked.

Note that (1) argues the closure may be premature and (3) may argue it
is safe. This script is not built to reach a preferred answer.

PRE-REGISTRATION (fixed before the run).

  A. corr(C, D) per band, full census of every even N up to 1.6e7,
     RAW and with the mask removed by the same modular enumeration
     used at increments 280-281. Pooling across bands -- what section 9
     did -- mixes scales and is reported only for comparison, never as
     the estimate.
     DECISION RULE: the closure STANDS if the de-masked per-band
     correlation is flat or falling and stays below 0.9. It is
     PREMATURE if the correlation rises monotonically and either
     crosses 0.9 or extrapolates to cross it within the measured
     range's own extension.

  B. Exponents beta_C and beta_R from the full census, each with a
     standard error, and the DIFFERENCE with the standard error of
     the difference -- computed from the same bands so the two are
     paired, which is what makes the difference's error smaller than
     the sum of the individual ones.
     DECISION RULE: "R grows faster than C" stands only if
     beta_R - beta_C exceeds 3x its own standard error.

  C. The replication of section 9's design for the difference, the way
     increment 281 did it for beta_C alone: 500 replicates, and the
     spread of (beta_R - beta_C) under that design. This says what the
     original design could have resolved, which is the question its
     conclusion depended on.

IDENTITY USED. r = S*(N - C) + R and D = N - r/S give R = S*(C - D)
exactly, so D = C - R/S. This is checked numerically rather than
assumed, because an identity that is assumed is an identity that can
be wrong in the code.
"""
import numpy as np
import math
import random
import time

QS = [3, 5, 7, 11, 13, 17, 19, 23]


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam, primes


def conv(X, a, b):
    n = 1
    while n < 2 * (X + 1):
        n *= 2
    A = np.zeros(n); A[: X + 1] = a
    B = np.zeros(n); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n)[: X + 1]


def demask(v, key):
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq))
    tot = np.bincount(inv, weights=v, minlength=len(uniq))
    return v - (tot / cnt)[inv], len(uniq)


def fit_se(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    b, a = np.polyfit(x, y, 1)
    r = y - (a + b * x)
    n = len(x)
    s2 = float((r ** 2).sum()) / max(n - 2, 1)
    sxx = float(((x - x.mean()) ** 2).sum())
    return b, math.sqrt(s2 / sxx), r


def main():
    X = 16_000_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    C = conv(X, mu, lam)
    r = conv(X, lam, lam)
    print(f"sieve + 2 convolutions  t={time.time()-t0:.0f}s", flush=True)

    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    Sv = S[Ns]
    Cv = C[Ns]
    rv = r[Ns]
    Dv = Ns - rv / Sv
    Rv = rv - Sv * (Ns - Cv)

    # the identity, checked rather than assumed
    lhs = Rv
    rhs = Sv * (Cv - Dv)
    rel = float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(lhs)))
    print(f"\nidentity R = S*(C-D):  max relative deviation {rel:.3e}"
          f"   {'OK' if rel < 1e-9 else 'FAILS'}")

    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

    print("\n(A) corr(C, D) per band -- the closure's own criterion,")
    print("    at full census instead of 80 samples per group")
    print(f"{'band':>21} {'count':>9} {'corr raw':>10} "
          f"{'corr de-masked':>15} {'logNmid':>8}")
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        n = int(sel.sum())
        if n <= 1000:
            b = hi
            continue
        c_, d_, R_, k_ = Cv[sel], Dv[sel], Rv[sel], key[sel]
        cr = float(np.corrcoef(c_, d_)[0, 1])
        cd, _ = demask(c_, k_)
        dd, _ = demask(d_, k_)
        cm = float(np.corrcoef(cd, dd)[0, 1])
        Rd, _ = demask(R_, k_)
        mid = math.sqrt(b * hi)
        rows.append((n, cr, cm, math.log(mid),
                     float(np.abs(cd).mean()), float(np.abs(Rd).mean())))
        print(f"{b:>9}-{hi:>11} {n:>9} {cr:>10.4f} {cm:>15.4f} "
              f"{math.log(mid):>8.3f}")
        b = hi

    corr_raw = np.array([x[1] for x in rows])
    corr_dem = np.array([x[2] for x in rows])
    Ls = np.array([x[3] for x in rows])
    print(f"\n    section 9 pooled its five groups and got 0.60.")
    print(f"    pooling ALL bands here (same mistake, for comparison): "
          f"{float(np.corrcoef(Cv, Dv)[0,1]):.4f}")
    print(f"    per-band de-masked: first {corr_dem[0]:.4f}, "
          f"last {corr_dem[-1]:.4f}")
    sl, se, _ = fit_se(Ls, corr_dem)
    print(f"    trend in de-masked corr vs log N: "
          f"{sl:+.5f} +/- {se:.5f} per unit log N")
    rising = sl > 3 * se
    print(f"    {'RISING' if rising else 'not rising'} at 3 s.e.")
    if rising:
        need = (0.9 - corr_dem[-1]) / sl
        print(f"    at this rate it would reach 0.9 at log N = "
              f"{Ls[-1] + need:.1f}, i.e. N ~ 1e"
              f"{(Ls[-1]+need)/math.log(10):.0f}")
    print(f"    mask effect on the correlation: "
          f"{float(np.mean(corr_raw - corr_dem)):+.4f} "
          f"(raw minus de-masked)")

    print("\n(B) exponents from the full census, de-masked, paired")
    yC = np.log([x[4] for x in rows])
    yR = np.log([x[5] for x in rows])
    bC, seC, resC = fit_se(Ls, yC)
    bR, seR, resR = fit_se(Ls, yR)
    bD, seD, _ = fit_se(Ls, yR - yC)     # the paired difference
    print(f"    beta_C = {bC:.4f} +/- {seC:.4f}")
    print(f"    beta_R = {bR:.4f} +/- {seR:.4f}")
    print(f"    beta_R - beta_C = {bD:+.4f} +/- {seD:.4f}   "
          f"({bD/seD:.1f} s.e.)")
    print(f"    'R grows faster than C' requires > 3 s.e.:  "
          f"{'STANDS' if bD > 3*seD else 'NOT SUPPORTED'}")

    print("\n(C) what could section 9's design have resolved?")
    print("    500 replicates of five groups of 80 consecutive even N")
    base = [120_000, 240_000, 480_000, 960_000, 1_900_000]
    rng = random.Random(20260807)
    diffs, bCs, bRs = [], [], []
    for _ in range(500):
        xs, yc, yr = [], [], []
        for N0 in base:
            off = 2 * rng.randrange(0, 20_000)
            a = N0 + off
            idx = np.arange((a - lo) // 2, (a - lo) // 2 + 80)
            xs.append(math.log(float(Ns[idx].mean())))
            yc.append(math.log(float(np.abs(Cv[idx]).mean())))
            yr.append(math.log(float(np.abs(Rv[idx]).mean())))
        p1 = np.polyfit(xs, yc, 1)[0]
        p2 = np.polyfit(xs, yr, 1)[0]
        bCs.append(p1); bRs.append(p2); diffs.append(p2 - p1)
    diffs = np.array(diffs)
    print(f"    beta_C  {np.mean(bCs):.4f} +/- {np.std(bCs):.4f}")
    print(f"    beta_R  {np.mean(bRs):.4f} +/- {np.std(bRs):.4f}")
    print(f"    difference {diffs.mean():+.4f} +/- {diffs.std():.4f}"
          f"   -> the recorded gap 0.096 is "
          f"{0.096/diffs.std():.1f} s.e. of that design")
    print(f"    fraction of replicates with difference <= 0: "
          f"{float((diffs <= 0).mean()):.3f}")
    print("DONE")


if __name__ == "__main__":
    main()
