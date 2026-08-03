"""A5-3 상위 스케일 표본 증강 — 4e7(+1800), 1e8(+2400) 후 npz 병합·재판정.

기존 s_binom_{X}.npz(seed 7)와 겹치지 않게 seed 11로 새 n 표본 추출,
합산 R² 재계산. 각 스케일 완료 즉시 저장(세션-단절 내성).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

for X, W, M_add in [(40_000_000, 1_200_000, 1800), (100_000_000, 2_000_000, 2400)]:
    old = np.load(f"s_binom_{X}.npz")
    used = set(int(v) for v in old["cands"])

    is_p = sieve(X + W + 10)
    ps_all = primes_upto((X + W) // 2 + 10)
    y = int(round(X ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]

    pool = [int(n) for n in range(X, X + W, 2)
            if n % 6 == 2 and int(n) not in used]
    rng = np.random.default_rng(11)
    cands = sorted(rng.choice(pool, min(M_add, len(pool)),
                              replace=False).tolist())

    s_new = np.empty(len(cands))
    N_new = np.empty(len(cands))
    for i, n in enumerate(cands):
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        N_new[i] = len(surv)
        s_new[i] = float(is_p[surv].mean())
        if (i + 1) % 300 == 0:
            print(f"  X={X:,}: {i+1}/{len(cands)}", flush=True)

    s = np.concatenate([old["s"], s_new])
    N = np.concatenate([old["N"], N_new])
    call = np.concatenate([old["cands"], np.array(cands)])
    np.savez(f"s_binom_{X}.npz", cands=call, s=s, N=N)

    sig_obs = s.std()
    sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N)))
    R2 = (sig_obs / sig_bin) ** 2
    se = R2 * np.sqrt(2 / (len(s) - 1))
    print(f"X={X:>13,}  합산 {len(s):>5}  s {s.mean():.5f}  "
          f"R2 = {R2:.4f} +- {se:.4f}", flush=True)
