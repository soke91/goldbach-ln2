# -*- coding: utf-8 -*-
"""
Increment 149: does the factorization law close the thin-progression
stamp too?

Inc. 141 measured thin-progression Mobius sums (moduli L > sqrt(y)) at
0.738 of half-normal -- sub-random. Inc. 148 proved (empirically) the
dispersion field factorizes as local mask x exact Gaussian. Same test
here: for S(L, r) = Sum_{a<=y, a=r mod L} mu(a),
  - count the SUPPORT (nonzero mu terms) per class directly,
  - condition on gcd(r, L) and its square structure,
  - report r1_eff = |S|/sqrt(support) per cell.
Prediction if the law is universal: clean cells give exactly 0.798
and the 0.738 global figure is a mixture artifact.
"""
import numpy as np, time
from math import gcd

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

def squarefree(n):
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False
        p += 1
    return True

def main():
    rng = np.random.default_rng(20260810)
    y = 200_000_000
    t0 = time.time()
    mu = mobius_upto(y)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)

    target = 8000
    S = np.zeros(target); sup = np.zeros(target); nt = np.zeros(target)
    g = np.zeros(target, dtype=np.int64)
    gsf = np.zeros(target, dtype=bool)   # gcd squarefree?
    done = 0
    while done < target:
        L = int(rng.integers(20000, 60000))
        r = int(rng.integers(1, L))
        a = np.arange(r, y + 1, L, dtype=np.int64)
        if len(a) < 300:
            continue
        v = mu[a].astype(np.int64)
        S[done] = v.sum()
        sup[done] = np.count_nonzero(v)
        nt[done] = len(a)
        gg = gcd(r, L)
        g[done] = gg
        gsf[done] = squarefree(gg)
        done += 1
        if done % 1000 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    def report(mask, tag):
        n = mask.sum()
        if n < 60:
            return
        good = mask & (sup > 0)
        r1e = np.mean(np.abs(S[good]) / np.sqrt(sup[good]))
        m2e = np.mean(S[good]**2 / sup[good])
        kte = np.mean(S[good]**4 / sup[good]**2) / max(m2e, 1e-12)**2
        # raw benchmark: against full length (what inc-141 measured)
        print(f"  {tag:26s} n={n:6d}  r1_eff={r1e:.3f}  m2_eff={m2e:.3f}  "
              f"kurt_eff={kte:.2f}", flush=True)

    print("--- thin progressions, zero-accounted (target 0.798/1.0/3.0) ---",
          flush=True)
    report(np.ones(target, dtype=bool), "all")
    report(g == 1, "gcd(r,L)=1")
    report((g > 1) & gsf, "gcd>1 squarefree")
    report((g > 1) & ~gsf, "gcd non-squarefree")
    print("--- raw normalization (inc-141 style, |S|/sqrt(#terms)) ---",
          flush=True)
    def report_raw(mask, tag):
        n = mask.sum()
        if n < 60:
            return
        r1r = np.mean(np.abs(S[mask]) / np.sqrt(nt[mask]))
        print(f"  {tag:26s} n={n:6d}  r1_raw={r1r:.3f}  "
              f"support_frac={np.mean(sup[mask]/nt[mask]):.3f}", flush=True)
    report_raw(np.ones(target, dtype=bool), "all (vs inc-141: 0.738)")
    report_raw(g == 1, "gcd(r,L)=1")
    report_raw(g > 1, "gcd(r,L)>1")
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
