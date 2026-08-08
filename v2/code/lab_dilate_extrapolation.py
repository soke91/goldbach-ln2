# -*- coding: utf-8 -*-
r"""
Testing the extrapolation of Remark {#rem:dilateprofile} out of sample.

WHAT IS AT STAKE

Remark [rem:dilateprofile] measured rho(k) = |A(N;k)| k/N and found it
flat against sqrt(k), then extrapolated: writing rho(k) = c' sqrt(k/N),

    B(N)/N = sum_{k<K} (log k) rho(k)/k ~ 2 c' sqrt(K/N) log K,

so the Goldbach threshold is kept up to K ~ N/(log N)^2.  That step
assumed c' is independent of N, and NOTHING TESTED IT.  It is almost
certainly false: c' = |H(N;k)|/sqrt(M) with M = N/k, and for a coin
the summands Lambda(N-mk) eps(m) are nonzero at about M/log N places
with size log N each, so |H_eps| ~ sqrt(M log N) and c'_eps ~
sqrt(log N).  If c' grows like sqrt(log N) the conclusion survives but
the polylog power does not: K* ~ N/(log N)^3, not (log N)^2.

The remark also left a gap it named but did not close: the predicted
level exponent near 1 against the directly measured K* exponent of
0.7057.  A model that is right should predict K* out of sample, so
that is the test here -- fit the model on the three smallest N and
predict K* at the two largest, against the measured
1353 and 2319 from lab_level_of_distribution.

BACKS: Remark {#rem:dilateprofile} in paper/wall_v3.md, whose
extrapolation this either supports or corrects.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  c'(N) = median_k |A(N;k)| / sqrt(N/k) over the resolved range
      GROWS with N -- it is not the constant the remark assumed.
  R2  c'(N)/sqrt(log N) is constant across the five N to within 15%.
  R3  The same holds for a coin, to within 15%: this is the
      calibration, since a coin must give |H| ~ sqrt(M log N) exactly.
  R4  Out of sample: fitting the model |A| = gamma sqrt(N/k) sqrt(log N)
      on N = 2e5, 4e5, 8e5 only, and solving
      sum_{k<K}(log k)|A| = S(N)(1-A(N))N at N = 1.6e6 and 3.2e6,
      predicts K* within a factor 2 of the measured 1353 and 2319.

REFUTATION RULE (fixed before the run)

  R1  REFUTED if c'(N) fails to increase between consecutive N.
  R2  REFUTED if max/min of c'/sqrt(log N) exceeds 1.15.
  R3  REFUTED if the coin's max/min exceeds 1.15.
  R4  REFUTED if either predicted K* is off by a factor 2 or more.

  All four gate.  R4 is the one that matters: it is the only
  out-of-sample statement in this line of work, and the remark's
  conclusion rests on the model it tests.

CITED BY: {#rem:extrap} in paper/.
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
OUT = os.path.join(ROOT, "results", "lab_dilate_extrapolation.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
# The Euler products below must not be tied to the
# measurement range: audit_constants.py shows the
# truncation reaches the sixth printed decimal.
CLIM = 4_000_000
FIT_NS = NS[:3]
KSTAR_MEAS = {1_600_000: 1353, 3_200_000: 2319}
MMIN = 1000
DRAWS = 4
SEED = 20260808


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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    rng = np.random.default_rng(SEED)

    say()
    say("  N          #k     c'(N)      c'/sqrt(log N)   coin c'    "
        "coin c'/sqrt(log N)   threshold/N")
    say("  " + "-" * 92)
    cp, cpn, cpe, cpen, thr_n = {}, {}, {}, {}, {}
    kk_all = {}
    for N in NS:
        v, PN, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                PN.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            PN.add(v)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))
        thr_n[N] = S_ * (1.0 - A_)

        kmax = N // MMIN
        ks = np.array([k for k in range(2, kmax + 1)
                       if mu[k] != 0 and all(k % q for q in PN)])
        kk_all[N] = ks

        def prog(sig):
            f = np.zeros(N, dtype=np.float64)
            idx = np.arange(1, N, dtype=np.int64)
            f[1:] = lam[1:N] * sig[N - idx]
            out = np.empty(ks.size)
            for i, k in enumerate(ks):
                r = N % int(k)
                out[i] = f[r::int(k)].sum() if r else f[int(k)::int(k)].sum()
            return np.abs(out)

        Amu = prog(mu.astype(np.float64))
        norm = np.sqrt(N / ks.astype(float))
        cp[N] = float(np.median(Amu / norm))
        cpn[N] = cp[N] / math.sqrt(math.log(N))

        es = []
        for t in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            es.append(float(np.median(prog(sig) / norm)))
        cpe[N] = float(np.mean(es))
        cpen[N] = cpe[N] / math.sqrt(math.log(N))
        say("  %-10d %-6d %-10.4f %-16.4f %-10.4f %-21.4f %.4f"
            % (N, ks.size, cp[N], cpn[N], cpe[N], cpen[N], thr_n[N]))

    say()
    seq = [cp[N] for N in NS]
    r1 = all(seq[i] < seq[i + 1] for i in range(len(seq) - 1))
    say("R1  c'(N) increasing: %s   %s" % (r1, "hold" if r1 else "REFUTED"))
    vn = [cpn[N] for N in NS]
    r2 = max(vn) / min(vn) <= 1.15
    say("R2  c'/sqrt(log N) spread = %.4f   (cap 1.15)   %s"
        % (max(vn) / min(vn), "hold" if r2 else "REFUTED"))
    ve = [cpen[N] for N in NS]
    r3 = max(ve) / min(ve) <= 1.15
    say("R3  coin c'/sqrt(log N) spread = %.4f   (cap 1.15)   %s"
        % (max(ve) / min(ve), "hold" if r3 else "REFUTED"))

    gamma = float(np.mean([cpn[N] for N in FIT_NS]))
    say()
    say("R4  out of sample: gamma fitted on N = %s is %.4f"
        % (", ".join(str(N) for N in FIT_NS), gamma))
    say("  N          K* measured   K* predicted   ratio   status")
    r4 = True
    for N in KSTAR_MEAS:
        ks = kk_all[N]
        pred = gamma * math.sqrt(math.log(N)) * np.sqrt(N / ks.astype(float))
        run = np.cumsum(np.log(ks.astype(float)) * pred)
        target = thr_n[N] * N
        j = int(np.searchsorted(run, target))
        # A search that ends at the top of its own range has not found
        # anything; reporting the endpoint as a solution is how the
        # first version of this rule produced a false 'hold'.
        clipped = j >= ks.size
        Kp = int(ks[min(j, ks.size - 1)])
        rat = Kp / KSTAR_MEAS[N]
        if clipped or not (0.5 < rat < 2.0):
            r4 = False
        say("  %-10d %-13d %-14d %-7.3f %s"
            % (N, KSTAR_MEAS[N], Kp, rat,
               "CLIPPED at the top of k <= N/%d -- not a solution" % MMIN
               if clipped else "interior"))
    say("  R4 %s" % ("hold" if r4 else "REFUTED"))

    say()
    say("  what the corrected model says about the level. With")
    say("  |A(N;k)| = gamma sqrt(N/k) sqrt(log N),")
    say("  B(N)/N ~ 2 gamma sqrt(log N) sqrt(K/N) log K, so the threshold")
    say("  is kept up to K ~ N / (log N (log K)^2) -- an exponent of 1 in")
    say("  N with THREE powers of log, not the two the remark wrote.")
    for N in (10 ** 6, 10 ** 9, 10 ** 12):
        L = math.log(N)
        K = N / (L ** 3)
        say("    N = 10^%-3d  N/(log N)^3 = %.3e   as an exponent: %.4f"
            % (int(round(math.log10(N))), K, math.log(K) / L))

    say()
    say("=" * 70)
    ok = r1 and r2 and r3 and r4
    say("R1 %s  R2 %s  R3 %s  R4 %s"
        % tuple("hold" if v else "REFUTED" for v in (r1, r2, r3, r4)))
    say("the model predicts K* out of sample; the polylog power is 3, "
        "not 2" if ok else "REFUTED")

    head = [
        "STATISTIC: c'(N) = median over the resolved k of",
        "           |A(N;k)| / sqrt(N/k), for mu and for coin signs on the",
        "           same support; c'(N)/sqrt(log N); and an out-of-sample",
        "           prediction of K*(N) obtained by fitting",
        "           |A| = gamma sqrt(N/k) sqrt(log N) on the three",
        "           smallest N and solving sum_{k<K}(log k)|A| =",
        "           S(N)(1-A(N))N at the two largest.",
        "NULL: the coin is the control and the calibration -- a sum of",
        "      independent signs must give |H| ~ sqrt(M log N) exactly, so",
        "      its c'/sqrt(log N) must be flat; if it were not, the",
        "      normalisation and not the arithmetic would be at fault.",
        "      Four draws per N, same support and same k-range as mu.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6; k over the squarefree k",
        "       coprime to N with 2 <= k <= N/1000, so the dilate length",
        "       M = N/k never falls below 1000; the fit uses only",
        "       N = 2e5, 4e5, 8e5 and the prediction only N = 1.6e6 and",
        "       3.2e6, whose K* = 1353, 2319 come from",
        "       lab_level_of_distribution.py; numpy default_rng seed",
        "       20260808.",
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
