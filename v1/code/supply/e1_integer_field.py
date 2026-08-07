# -*- coding: utf-8 -*-
"""
Increment 155: does the factorization law extend to the INTEGER-indexed
field -- the actual E1 object D(k) = Sum_m mu(m) mu(N-mk)?

The prime-indexed field closed completely (inc. 144-154). But
verify_all's V3/V4 (integer-indexed) read 0.688/0.290 even with
support normalization: a residual deficit. Suspect: a computable
local SIGN correlation between mu(m) and mu(N-mk) (the mu-mu
singular-series analog; for every prime q the pair (m, N-mk) has a
joint local structure -- unlike the prime-indexed case where p was
restricted to units).

Test:
  (1) per-k per-term mean tbar_k = D(k)/support -- is it nonzero and
      a function of k's local data? (group by k mod 4, k mod 9,
      k mod 25 -- class consistency)
  (2) subtract the CLASS mean (not per-k mean): does
      |D - class_mean * support| / sqrt(support) return to 0.798?
  (3) enumeration check at q=2,3: local mean over m mod q^3 classes
      computed by brute force on a small window, compared to the
      class-mean map.
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

def main():
    N = 99_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    SQ = int(N ** 0.5)

    ks = np.arange(800, 1601)
    D = np.zeros(len(ks)); sup = np.zeros(len(ks))
    tbar = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        t = mu[ms].astype(np.int64) * mu[N - k * ms]
        D[i] = t.sum(); sup[i] = np.count_nonzero(t)
        tbar[i] = D[i] / max(sup[i], 1)
        if i % 200 == 199:
            print(f"k {i+1}/{len(ks)}  t={time.time()-t0:.0f}s",
                  flush=True)

    r1_raw = np.mean(np.abs(D) / np.sqrt(np.maximum(sup, 1)))
    print(f"(0) raw r1_eff = {r1_raw:.3f}  (V3-style; deficit vs 0.798)",
          flush=True)
    print(f"    global per-term mean = {tbar.mean():+.5f}  "
          f"(se {tbar.std()/np.sqrt(len(ks)):.5f})", flush=True)

    # (1) class structure of tbar
    print("(1) class means of tbar:", flush=True)
    for label, key in (("k mod 4", ks % 4), ("k mod 9", ks % 9),
                       ("k mod 25", ks % 25)):
        out = []
        for c in np.unique(key):
            m = key == c
            if m.sum() < 15:
                continue
            mm = tbar[m].mean(); se = tbar[m].std()/np.sqrt(m.sum())
            out.append(f"{c}:{mm:+.4f}({mm/max(se,1e-9):+.1f}z)")
        print(f"    {label}: " + "  ".join(out), flush=True)

    # (2) subtract class-mean model (k mod 36 = joint 4x9)
    key = (ks % 4).astype(int) * 9 + (ks % 9).astype(int)
    model = np.zeros(len(ks))
    for c in np.unique(key):
        m = key == c
        model[m] = tbar[m].mean()
    resid = D - model * sup
    r1_adj = np.mean(np.abs(resid) / np.sqrt(np.maximum(sup, 1)))
    print(f"(2) after class-mean (k mod 36) subtraction: r1_eff = "
          f"{r1_adj:.3f}  (target 0.798)", flush=True)

    # leave-one-out honesty check: model from other k's only
    resid2 = np.zeros(len(ks))
    for c in np.unique(key):
        m = np.nonzero(key == c)[0]
        if len(m) < 2:
            resid2[m] = D[m] - model[m] * sup[m]
            continue
        s_all = tbar[m].sum()
        for j in m:
            mj = (s_all - tbar[j]) / (len(m) - 1)
            resid2[j] = D[j] - mj * sup[j]
    r1_loo = np.mean(np.abs(resid2) / np.sqrt(np.maximum(sup, 1)))
    print(f"    leave-one-out version: r1_eff = {r1_loo:.3f}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
