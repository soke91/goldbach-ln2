# -*- coding: utf-8 -*-
"""
Forge kill-test K1 (increment 171): can the ladder orbit REPRESENT
D(k)?

Design K1 (multiplicative Fejer kernel) requires that D(k) be
recoverable, up to log-power error, from its dilation orbit
{D(sk) : s squarefree over small primes}. The exact ladder gives the
p-divisible SUB-SUM exactly; the design question is whether the
orbit's span captures the whole of D(k) -- i.e., whether the
"p-coprime core" shrinks under orbit extension.

Test: least-squares R^2 of D(k) on the orbit vector
  {D(sk) : s | 2*3*5*7*11*13, s > 1, sk <= Kmax}
over 400 values of k, N = 2e8. Escalation thresholds (pre-registered):
  R^2 >= 0.95  -> K1 alive (write the reduction, adversarial review)
  R^2 in (0.9, 0.95) -> marginal (repeat at second N before deciding)
  R^2 <= 0.9   -> K1 dead (close; consistent with the cascade death)
Also reported: one-prime-layer R^2 (the old cascade number, sanity)
and the residual field's normalized size (what fraction of energy the
orbit can never see).
"""
import numpy as np, time
from itertools import combinations

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

PR = [2, 3, 5, 7, 11, 13]

def squarefree_divisors():
    out = []
    for r in range(1, len(PR) + 1):
        for c in combinations(PR, r):
            s = 1
            for q in c:
                s *= q
            out.append(s)
    return sorted(out)

def D_of(mu, N, k, SQ):
    ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
    t = mu[ms].astype(np.int64) * mu[N - k * ms]
    return float(t.sum()), int(np.count_nonzero(t))

def main():
    rng = np.random.default_rng(20260830)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    SQ = int(N ** 0.5)
    svals = squarefree_divisors()
    print(f"orbit divisors: {len(svals)} (s | 30030)", flush=True)

    # choose k so that s*k stays in a sane range: k in [500, 900],
    # max s*k = 30030*900 ~ 2.7e7 << N/  (m-range stays long)
    ks = rng.choice(np.arange(500, 900), size=400, replace=False)
    y = np.zeros(len(ks))
    X = np.zeros((len(ks), len(svals)))
    sup = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        d0, s0 = D_of(mu, N, k, SQ)
        y[i] = d0; sup[i] = s0
        for j, s in enumerate(svals):
            X[i, j], _ = D_of(mu, N, s * k, SQ)
        if i % 50 == 49:
            print(f"k {i+1}/400  t={time.time()-t0:.0f}s", flush=True)

    # full-orbit least squares
    coef, res, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    R2_full = 1 - ss_res / ss_tot
    # one-layer (single primes only) for comparison
    X1 = X[:, :len(PR)]
    # columns: first len(PR) svals sorted -- verify they are the primes
    prime_cols = [j for j, s in enumerate(svals) if s in PR]
    X1 = X[:, prime_cols]
    coef1, _, _, _ = np.linalg.lstsq(X1, y, rcond=None)
    R2_one = 1 - float(np.sum((y - X1 @ coef1) ** 2)) / ss_tot

    resid_energy = ss_res / float(np.sum(sup))
    print("=== K1 KILL-TEST RESULT ===", flush=True)
    print(f"R2 one-layer (primes only) = {R2_one:.3f}  "
          f"(cascade history ~0.68)", flush=True)
    print(f"R2 full orbit (63 divisors) = {R2_full:.3f}", flush=True)
    print(f"residual energy / support = {resid_energy:.3f}  "
          f"(fraction of unit-Gaussian budget the orbit cannot see)",
          flush=True)
    verdict = ("ALIVE" if R2_full >= 0.95 else
               "MARGINAL" if R2_full > 0.9 else "DEAD")
    print(f"pre-registered verdict: {verdict}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
