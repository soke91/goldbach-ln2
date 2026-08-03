"""정리 1의 오차항 구조 검증 — (share − 1/(1+ln2))·ln X 가 상수로 수렴하는가."""

import numpy as np

from goldbach.sieve import primes_upto, sieve

target = 1 / (1 + np.log(2))
print(f"{'X':>12} {'share':>8} {'(share-lim)*lnX':>16}")
for X in [10 ** 6, 10 ** 7, 10 ** 8]:
    is_p = sieve(X)
    y = int(round(X ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y)]
    W = 4_000_000 if X >= 10 ** 7 else 800_000
    m = np.arange(X - W + 1, X, 2)
    alive = np.ones(len(m), dtype=bool)
    for q in qs:
        if q > 2:
            alive &= (m % q != 0)
    surv = m[alive]
    share = float(is_p[surv].mean())
    print(f"{X:>12,} {share:>8.5f} {(share - target) * np.log(X):>16.4f}",
          flush=True)
print(f"limit 1/(1+ln2) = {target:.5f}")
