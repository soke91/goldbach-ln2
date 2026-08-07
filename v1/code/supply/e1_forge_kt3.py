# -*- coding: utf-8 -*-
"""
Forge kill-test K3 (increment 175): the Wishart moment route.

Design K3: prove the E1 band ratio via trace moments of the Gram
matrix M = normalized C-matrix (lambda_max <= tr(M^{2j})^{1/2j}).
Structural risk (pre-registered): tr(M^{2j}) expands into 2j-cycles
of pair correlations -- 2j-fold mu-correlation sums, which are HARDER
binary objects unless they carry exploitable extra cancellation.
Therefore:
  ALIVE  <=> the real 4-cycle trace moment sits significantly BELOW
             the Wishart-null prediction (sub-Wishart cancellation =
             new structure the moment method could consume);
  DEAD   <=> real trace moments match the null (confirms Conjecture L
             at cycle level but makes the route circular).
Measurement: 400 x 400 band matrix (the inc-151/154 pipeline),
tr(M^2), tr(M^3), tr(M^4) real vs 12 null draws (random-sign rows on
the real support pattern).
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

def tr_moments(R):
    Rf = R.astype(np.float32)
    C = Rf @ Rf.T
    NZ = np.abs(Rf) @ np.abs(Rf.T)
    np.fill_diagonal(C, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        M = np.where(NZ > 50, C / np.sqrt(np.maximum(NZ, 1)), 0.0)
    np.fill_diagonal(M, 0.0)
    M = ((M + M.T) / 2).astype(np.float64)
    M2 = M @ M
    t2 = float(np.trace(M2))
    t3 = float(np.trace(M2 @ M))
    t4 = float(np.trace(M2 @ M2))
    return t2, t3, t4

def main():
    rng = np.random.default_rng(20260831)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 3000, 3400
    ks = np.arange(K0, K1); K = len(ks)
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)
    pp = ps[ps <= (N - 2) // (K1 - 1)]

    R = np.empty((K, len(pp)), dtype=np.int8)
    for i, k in enumerate(ks):
        R[i] = mu[N - pp * int(k)]
    t2r, t3r, t4r = tr_moments(R)
    print(f"real: tr(M^2)={t2r:.0f}  tr(M^3)={t3r:.0f}  "
          f"tr(M^4)={t4r:.0f}", flush=True)

    support = (R != 0)
    vals2, vals3, vals4 = [], [], []
    for d in range(12):
        signs = rng.choice(np.array([-1, 1], dtype=np.int8),
                           size=R.shape)
        Rn = np.where(support, signs, 0).astype(np.int8)
        a, b, c = tr_moments(Rn)
        vals2.append(a); vals3.append(b); vals4.append(c)
        if d % 4 == 3:
            print(f"  null {d+1}/12  t={time.time()-t0:.0f}s", flush=True)
    for tag, real, vals in (("tr(M^2)", t2r, vals2),
                            ("tr(M^3)", t3r, vals3),
                            ("tr(M^4)", t4r, vals4)):
        m = float(np.mean(vals)); s = float(np.std(vals))
        z = (real - m) / max(s, 1e-9)
        print(f"{tag}: real={real:.0f}  null={m:.0f} +- {s:.0f}  "
              f"z={z:+.2f}", flush=True)
    z4 = (t4r - np.mean(vals4)) / max(np.std(vals4), 1e-9)
    print("=== K3 KILL-TEST ===", flush=True)
    print("verdict:",
          "ALIVE (sub-Wishart 4-cycle cancellation)" if z4 <= -3 else
          "DEAD (cycle moments match the null -- the moment route is "
          "circular: tr(M^{2j}) consumes 2j-fold mu-correlations with "
          "no extra structure to exploit)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
