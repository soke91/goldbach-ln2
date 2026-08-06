# -*- coding: utf-8 -*-
"""
Construction experiment C1 (increment 181): is the abelian spectrum
of the pair field mask-exact?

For fixed k, t_m = mu(m) mu(N - mk) on m in (sqrt N, N/k]. The
support mask is (quasi-)periodic in m modulo small prime powers, so
|FFT(t)|^2 necessarily shows structure at rationals a/q -- that much
the mask predicts. The construction question: does the REAL field
carry rational-spectrum energy EXCEEDING the mask x random-sign null?
Any excess = an abelian handle (seed of a C-I representation); none
=> the entire abelian spectrum is mask-exact and the representation
must be non-abelian.

Method per k: z-pad t to length 2^18, FFT, total rational-peak
energy E_Q = sum over a/q (q <= 32, gcd(a,q)=1) of |F|^2 in the bin
nearest a/q. Null: 8 draws of random signs on the real support.
Report z = (E_Q^real - mean)/std per k; pre-registered: alive iff
z >= 4 at >= 2 of 6 k-values.
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

def rational_bins(L, Q):
    idx = set()
    for q in range(1, Q + 1):
        for a in range(q):
            if gcd(a, q) == 1 or (a == 0 and q == 1):
                idx.add(int(round(a / q * L)) % L)
    return np.array(sorted(idx), dtype=np.int64)

def peak_energy(t, L, bins):
    F = np.fft.rfft(t, n=L)
    E = np.abs(F) ** 2
    tot = float(E.sum())
    b = bins[bins < len(E)]
    return float(E[b].sum()), tot

def main():
    rng = np.random.default_rng(20260905)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    SQ = int(N ** 0.5)
    L = 1 << 18
    bins = rational_bins(L, 32)
    print(f"rational bins (q<=32): {len(bins)}", flush=True)

    kvals = [2001, 2503, 3001, 3499, 3803, 2251]
    alive_hits = 0
    for k in kvals:
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        t = (mu[ms].astype(np.int64) *
             mu[N - k * ms]).astype(np.float64)
        Er, Tr = peak_energy(t, L, bins)
        sup = (t != 0)
        nulls = []
        for _ in range(8):
            s = rng.choice([-1.0, 1.0], size=t.shape)
            tn = np.where(sup, s, 0.0)
            En, _ = peak_energy(tn, L, bins)
            nulls.append(En)
        m = float(np.mean(nulls)); sd = float(np.std(nulls))
        z = (Er - m) / max(sd, 1e-9)
        frac = Er / Tr
        hit = z >= 4
        alive_hits += int(hit)
        print(f"k={k}: E_rat/E_tot={frac:.4f}  real={Er:.3e}  "
              f"null={m:.3e} +- {sd:.2e}  z={z:+.2f}"
              f"{'  <-- HIT' if hit else ''}  t={time.time()-t0:.0f}s",
              flush=True)

    print(f"=== C1: {alive_hits}/6 hits ===", flush=True)
    print("verdict:",
          "C-I ALIVE (abelian excess found)" if alive_hits >= 2 else
          ("MARGINAL (repeat)" if alive_hits == 1 else
           "C-I CLOSED (abelian spectrum is mask-exact; the "
           "representation must be non-abelian)"), flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
