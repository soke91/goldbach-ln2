# -*- coding: utf-8 -*-
"""
Forge R2 fix (increment 179b): part (a) rerun on odd k only
(even k with v2(N)=1 are annihilated cells -> 0/0 = nan in the
first run; the regression part was unaffected and read dead).
Original header follows.

Forge round-2 kill-test R2 (increment 179): does the determinant
phase see the pair field?

The pair constraint is the determinant equation; when DI/Kuznetsov
machinery controls such objects, the completed sums carry
Kloosterman-type phases e(N k'bar / k) (modular inverses across the
pair moduli). If the C-field has a component coherent with these
phases, a spectral construction has a handle; if it is phase-blind,
the determinant/automorphic channel closes at the visibility level.

Tests (600 pairs per fixed k, 40 k values, coprime pairs):
  (a) coherent twist gain: G(k) = |Sum_{k'} C~ e(-N inv(k',k)/k)|^2
      / Sum |C~|^2  -- under L this is ~1 (no coherence); a phase
      that matches hidden structure gives G >> 1. Also the mirror
      phase e(-N inv(k,k')/k').
  (b) per-pair regression of C~ on cos/sin of both phases vs a
      random-phase control.
Pre-registered: ALIVE iff mean G >= 2 or regression capture >= 2x
random control; else DEAD.
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

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def main():
    rng = np.random.default_rng(20260904)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 2000, 4000
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)

    odd = np.arange(K0 + 1, K1, 2)
    kvals = rng.choice(odd, size=40, replace=False)
    Gs, Gs_m, Gs_r = [], [], []
    allC, allph1, allph2 = [], [], []
    for ik, k in enumerate(kvals):
        k = int(k)
        kps = [int(x) for x in range(K0, K1)
               if x != k and gcd(x, k) == 1][:600]
        Cn = np.zeros(len(kps))
        ph1 = np.zeros(len(kps)); ph2 = np.zeros(len(kps))
        for i, kp in enumerate(kps):
            pmax = (N - 2) // max(k, kp)
            pp = ps[ps <= pmax]
            t = (mu[N - pp * k].astype(np.int64) *
                 mu[N - pp * kp].astype(np.int64))
            nz = int(np.count_nonzero(t))
            Cn[i] = t.sum() / np.sqrt(nz) if nz > 100 else np.nan
            ph1[i] = 2*np.pi*((N % k) * pow(kp, -1, k) % k) / k
            ph2[i] = 2*np.pi*((N % kp) * pow(k, -1, kp) % kp) / kp
        ok = ~np.isnan(Cn)
        v = Cn[ok]; p1 = ph1[ok]; p2 = ph2[ok]
        e1 = np.exp(-1j * p1); e2 = np.exp(-1j * p2)
        base = float(np.sum(v**2))
        G1 = abs(np.sum(v * e1))**2 / base
        G2 = abs(np.sum(v * e2))**2 / base
        rphase = rng.uniform(0, 2*np.pi, size=len(v))
        Gr = abs(np.sum(v * np.exp(-1j*rphase)))**2 / base
        Gs.append(G1); Gs_m.append(G2); Gs_r.append(Gr)
        allC.append(v); allph1.append(p1); allph2.append(p2)
        if ik % 10 == 9:
            print(f"k {ik+1}/40  t={time.time()-t0:.0f}s", flush=True)

    Gs = np.array(Gs); Gs_m = np.array(Gs_m); Gs_r = np.array(Gs_r)
    print(f"(a) coherent gain: phase1 mean G = {Gs.mean():.2f}  "
          f"phase2 = {Gs_m.mean():.2f}  random = {Gs_r.mean():.2f}  "
          f"(L prediction ~1)", flush=True)

    v = np.concatenate(allC)
    p1 = np.concatenate(allph1); p2 = np.concatenate(allph2)
    Xd = np.column_stack([np.cos(p1), np.sin(p1),
                          np.cos(p2), np.sin(p2)])
    coef, _, _, _ = np.linalg.lstsq(Xd, v, rcond=None)
    R2d = 1 - float(np.sum((v - Xd@coef)**2))/float(np.sum((v-v.mean())**2))
    r2rs = []
    for _ in range(8):
        rp = rng.uniform(0, 2*np.pi, size=(len(v), 2))
        Xr = np.column_stack([np.cos(rp[:,0]), np.sin(rp[:,0]),
                              np.cos(rp[:,1]), np.sin(rp[:,1])])
        c2, _, _, _ = np.linalg.lstsq(Xr, v, rcond=None)
        r2rs.append(1 - float(np.sum((v - Xr@c2)**2)) /
                    float(np.sum((v-v.mean())**2)))
    mr = float(np.mean(r2rs))
    print(f"(b) regression: R2_det = {R2d:.5f}  R2_random = {mr:.5f} "
          f"+- {np.std(r2rs):.5f}", flush=True)
    alive = (Gs.mean() >= 2) or (Gs_m.mean() >= 2) or \
            (mr > 0 and R2d >= 2*mr)
    print("=== R2 KILL-TEST ===", flush=True)
    print("verdict:",
          "ALIVE (determinant phase visible)" if alive else
          "DEAD (the pair field is phase-blind to its own determinant "
          "structure; the direct spectral-visibility channel closes)",
          flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
