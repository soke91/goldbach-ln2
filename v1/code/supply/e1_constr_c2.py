# -*- coding: utf-8 -*-
"""
Construction experiment C2 (increment 182): modular-inverse domain
visibility -- the threshold test for a non-abelian (C-II) basis.

Kloosterman/Kuznetsov territory begins exactly where the modular
INVERSE of the summation variable carries structure: S(a,b;c) =
Sum e((ar + b rbar)/c). Re-index the row field by the inverse:
    h_k(r) = Sum_{p : pbar = r (mod k)} mu(N - pk),   r in (Z/k)*
    g_k(a) = FFT_k(h)  -- the inverse-domain spectrum.
If the real field's inverse-domain spectrum matches the mask-null
(random signs on the same p-support), there is no row-level shadow
of any Kloosterman-type structure to build C-II on.

Per k (odd k, 300 values in [2001, 3999]):
  - E_special = |g(a)|^2 summed over a in {+-1, +-N, +-Nbar mod k}
    (the frequencies DI-type completions distinguish),
  - E_max = max_a |g(a)|^2 (largest inverse-domain line),
  - both compared to 8 mask-null draws (z-scores).
Pre-registered (family of 300, Bonferroni-aware): C-II row-shadow
ALIVE iff >= 3 k with z_special >= 4.5 or a systematic mean-shift
|mean z| >= 0.5 across the family; else absent.
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

def spec_energy(vals, invs, k, specials):
    h = np.zeros(k)
    np.add.at(h, invs, vals)
    g = np.fft.rfft(h)
    E = np.abs(g) ** 2
    es = 0.0
    for a in specials:
        aa = a % k
        aa = min(aa, k - aa)
        if aa < len(E):
            es += float(E[aa])
    return es, float(E[1:].max())  # exclude DC for max

def main():
    rng = np.random.default_rng(20260906)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    P0 = N // (2 * 4200); P1 = 2 * P0
    ps = primes_in(P0, P1)

    kvals = rng.choice(np.arange(2001, 4000, 2), size=300,
                       replace=False)
    zs_sp, zs_mx = [], []
    hits = 0
    for ik, k in enumerate(kvals):
        k = int(k)
        pmax = (N - 2) // k
        pp = ps[(ps <= pmax) & (ps % k != 0)]
        vals = mu[N - pp * k].astype(np.float64)
        # modular inverses of p mod k (vectorized via pow loop)
        invs = np.array([pow(int(p), -1, k) for p in pp],
                        dtype=np.int64)
        Nb = pow(N % k, -1, k) if np.gcd(N % k, k) == 1 else 1
        specials = [1, k - 1, N % k, (-N) % k, Nb, (-Nb) % k]
        Er, Mr = spec_energy(vals, invs, k, specials)
        sup = vals != 0
        eN, mN = [], []
        for _ in range(8):
            s = rng.choice([-1.0, 1.0], size=vals.shape)
            vn = np.where(sup, s, 0.0)
            a, b = spec_energy(vn, invs, k, specials)
            eN.append(a); mN.append(b)
        zsp = (Er - np.mean(eN)) / max(np.std(eN), 1e-9)
        zmx = (Mr - np.mean(mN)) / max(np.std(mN), 1e-9)
        zs_sp.append(zsp); zs_mx.append(zmx)
        if zsp >= 4.5:
            hits += 1
            print(f"  HIT k={k}: z_special={zsp:+.2f}", flush=True)
        if ik % 50 == 49:
            print(f"k {ik+1}/300  t={time.time()-t0:.0f}s", flush=True)

    zs_sp = np.array(zs_sp); zs_mx = np.array(zs_mx)
    msp = zs_sp.mean(); ssp = zs_sp.std() / np.sqrt(len(zs_sp))
    mmx = zs_mx.mean()
    print(f"family: mean z_special = {msp:+.3f} +- {ssp:.3f}  "
          f"mean z_max = {mmx:+.3f}  hits(z>=4.5): {hits}", flush=True)
    alive = (hits >= 3) or (abs(msp) >= 0.5)
    print("=== C2 ===", flush=True)
    print("verdict:",
          "C-II row-shadow ALIVE" if alive else
          "C-II row-shadow ABSENT (the inverse-domain spectrum is "
          "mask-exact; no Kloosterman-type structure is visible at "
          "row level)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
