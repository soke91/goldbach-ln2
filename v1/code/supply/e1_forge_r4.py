# -*- coding: utf-8 -*-
"""
Forge round-3 kill-test R4 (increment 194): does the dilate field
retain any residue of the EXACT divisor-switch cancellation?

New asset (increment 193, Theorem A). Applying the divisor switch to
the dilate field over its FULL ranges gives an exact identity:

    Sum_{k>=1} Sum_{m: mk <= N-1} mu(m) mu(N-mk)
       = Sum_{u<N} mu(N-u) * Sum_{m|u} mu(m)
       = Sum_{u<N} mu(N-u) * [u = 1]  =  mu(N-1).

That is PERFECT cancellation -- an O(1) bound on a double sum of
~N log N terms, far beyond square root. E1 imposes two restrictions
that make the inner divisor sum incomplete: m > sqrt(N) (the type-II
cut) and k in a dyadic band. The question this test settles is whether
any of the perfect cancellation survives those restrictions.

Under Conjecture L the field {D(k)} is mean-zero and INDEPENDENT
across k at unit-Gaussian scale, so summing D over a block of B
consecutive k must grow the variance by exactly B: the normalized
statistic

    ratio(B) = Sum_j |S_B(j)|^2 / Sum_k supp(k),
       S_B(j) = Sum_{k in block j} D(k),
       supp(k) = #{m : mu(m) != 0 and mu(N-mk) != 0},

is independent of B (and equals the E1 statistic ~1 at B = 1). A
surviving residue of the exact identity would instead show block sums
much SMALLER than independence predicts -- ratio(B) falling with B.

PRE-REGISTERED (fixed before the run):
  ALIVE  iff  ratio(B) <= 0.5 * ratio(1)  for at least two block sizes
              B >= 8, at BOTH values of N.
  DEAD   otherwise (the restrictions destroy the identity's coherence
              completely; the L1 route cannot fund an L2 bound).

Control decomposition. The same statistic is computed for the
FULL-m field D_full(k) = Sum_{m <= N/k} mu(m) mu(N-mk), which differs
from D(k) only by the short-m (type-I) part. If coherence appears in
D_full but not in D, the exact identity's cancellation lives entirely
in the classically controlled type-I slice -- informative, but not a
lever for E1.

Sanity check: the identity itself is verified by brute force at small N.
"""
import numpy as np
import time


def mobius_upto(X):
    """mu(0..X) as int8. int32 cofactor tracker to keep memory sane."""
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool)
    pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int32)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    mu[0] = 0
    return mu


def identity_check(Nsmall):
    """Brute-force check of Sum_k Sum_m mu(m) mu(N-mk) = mu(N-1)."""
    mu = mobius_upto(Nsmall)
    tot = 0
    for k in range(1, Nsmall):
        ms = np.arange(1, (Nsmall - 1) // k + 1, dtype=np.int64)
        if ms.size == 0:
            break
        tot += int((mu[ms].astype(np.int64)
                    * mu[Nsmall - k * ms].astype(np.int64)).sum())
    return tot, int(mu[Nsmall - 1])


def field(mu, N, ks, full):
    """D(k) and supp(k) for k in ks; full=True uses m from 1."""
    SQ = int(N ** 0.5)
    lo = 1 if full else SQ + 1
    D = np.zeros(len(ks))
    S = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(lo, N // k + 1, dtype=np.int64)
        a = mu[ms].astype(np.int64)
        b = mu[N - k * ms].astype(np.int64)
        D[i] = float((a * b).sum())
        S[i] = float(np.count_nonzero(a * b))
    return D, S


def ratios(D, S, blocks):
    out = {}
    tot_supp = float(S.sum())
    for B in blocks:
        nb = len(D) // B
        s = D[:nb * B].reshape(nb, B).sum(axis=1)
        # normalise by the support actually inside the used k-range
        used = float(S[:nb * B].sum())
        out[B] = float(np.sum(s ** 2)) / used
    out['_supp'] = tot_supp
    return out


def main():
    print("=== R4: sanity check of the exact identity ===", flush=True)
    for Ns in (5000, 20000):
        tot, pred = identity_check(Ns)
        print(f"  N={Ns}: Sum_k Sum_m mu(m)mu(N-mk) = {tot}   "
              f"mu(N-1) = {pred}   {'OK' if tot == pred else 'MISMATCH'}",
              flush=True)

    BLOCKS = [1, 8, 64, 512]
    K0, NK = 2000, 2048
    ks = np.arange(K0, K0 + NK)

    for N in (99_999_998, 99_960_002):
        t0 = time.time()
        mu = mobius_upto(N)
        print(f"\n=== N = {N}  (mu ready {time.time()-t0:.0f}s) ===",
              flush=True)
        print(f"k band [{K0}, {K0+NK}), sqrt(N) = {int(N**0.5)}",
              flush=True)
        for full in (False, True):
            tag = "D_full (all m)" if full else "D (m > sqrt N)"
            D, S = field(mu, N, ks, full)
            r = ratios(D, S, BLOCKS)
            line = "  ".join(f"B={B}: {r[B]:.4f}" for B in BLOCKS)
            print(f"  {tag:16s}  {line}   "
                  f"[mean D = {D.mean():.1f}, t={time.time()-t0:.0f}s]",
                  flush=True)
            if not full:
                base = r[1]
                hits = [B for B in BLOCKS if B >= 8 and r[B] <= 0.5 * base]
                print(f"      pre-registered hits (ratio <= 0.5*B1): "
                      f"{hits if hits else 'none'}", flush=True)
        del mu
    print("\nDONE", flush=True)


if __name__ == "__main__":
    main()
