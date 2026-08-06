# -*- coding: utf-8 -*-
"""
Forge kill-test K2 (increment 173): does a MANUFACTURED congruence
buy anything?

Design K2: split the k-average into progressions mod d (cost: factor
d by Cauchy-Schwarz); within a class every pair satisfies k = k'
(mod d) and the pair determinant N(k'-k) = N j d carries the factor d
-- the divisibility that powers Lichtman's f-substitution. The design
is ALIVE only if congruent pairs exhibit d-dependent structure worth
more than the factor-d cost:
  (a) suppressed variance or coherent sign of C_{k,k+h} when d | h
      with d large (conductor-collapse shadow), or
  (b) autocorrelation of k -> C_{k,k+h} (secondary-averaging gain).
Conjecture L (pure) predicts: statistics of C_{k,k+h} are h-blind
(unit Gaussian on support, mean 0, autocorr 0) -> gain 0 -> DEAD.

Pre-registered thresholds: any per-h |mean|/SE >= 4, or m2_eff
deviation from its h=1 baseline by >= 5 sigma, or lag-1 autocorr
|rho| >= 0.15 across two h values => ALIVE coordinate; else DEAD.

Setup: N = 199_999_998, prime-indexed pairs, k in [2000, 4000],
h in {1, 2, 4, 6, 8, 16, 30, 32, 64, 210}; 600 consecutive-k samples
per h (consecutive so the autocorrelation is measurable).
"""
import numpy as np, time

def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool); pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    return mu

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def main():
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0 = 2000
    HS = [1, 2, 4, 6, 8, 16, 30, 32, 64, 210]
    P0 = N // (2 * 4400); P1 = 2 * P0
    ps = primes_in(P0, P1)

    base_m2 = None
    flags = []
    for h in HS:
        kvals = np.arange(K0, K0 + 600)
        C = np.zeros(len(kvals)); sup = np.zeros(len(kvals))
        for i, k in enumerate(kvals):
            k = int(k); kp = k + h
            pmax = (N - 2) // kp
            pp = ps[ps <= pmax]
            t = (mu[N - pp * k].astype(np.int64) *
                 mu[N - pp * kp].astype(np.int64))
            C[i] = t.sum(); sup[i] = np.count_nonzero(t)
        good = sup > 100
        Cg = C[good]; sg = sup[good]
        n = good.sum()
        mean_n = np.mean(Cg / np.sqrt(sg))
        se = np.std(Cg / np.sqrt(sg)) / np.sqrt(n)
        m2 = np.mean(Cg**2 / sg)
        m2_se = np.std(Cg**2 / sg) / np.sqrt(n)
        # lag-1 autocorr over consecutive viable k
        x = np.where(sup > 100, C / np.sqrt(np.maximum(sup, 1)), np.nan)
        pairs = [(x[i], x[i+1]) for i in range(len(x)-1)
                 if not (np.isnan(x[i]) or np.isnan(x[i+1]))]
        if len(pairs) > 50:
            a = np.array([p[0] for p in pairs])
            b = np.array([p[1] for p in pairs])
            rho = float(np.corrcoef(a, b)[0, 1])
        else:
            rho = float('nan')
        if base_m2 is None:
            base_m2 = (m2, m2_se)
        zm = mean_n / max(se, 1e-12)
        zdev = (m2 - base_m2[0]) / max(np.hypot(m2_se, base_m2[1]), 1e-12)
        flag = (abs(zm) >= 4) or (h > 1 and abs(zdev) >= 5) or \
               (abs(rho) >= 0.15)
        flags.append((h, flag))
        print(f"h={h:4d}  n={n:4d}  mean_n={mean_n:+.4f}({zm:+.1f}z)  "
              f"m2={m2:.3f}({zdev:+.1f}z vs h=1)  rho1={rho:+.3f}"
              f"{'  <-- FLAG' if flag else ''}", flush=True)
        print(f"  t={time.time()-t0:.0f}s", flush=True)

    nflag = sum(1 for _, f in flags if f)
    print(f"=== K2 KILL-TEST: {nflag} flagged h-values ===", flush=True)
    print("verdict:", "ALIVE coordinate found" if nflag >= 2 else
          ("MARGINAL (single flag -- repeat)" if nflag == 1 else
           "DEAD (no d-dependent structure; the manufactured "
           "congruence buys nothing against the factor-d cost)"),
          flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
