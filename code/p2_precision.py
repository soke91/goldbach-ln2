"""sub-binomial 비 r_sb 의 다중 스케일 정밀 판정 — √ln2 = 0.8326 인가?

스케일 4개 × 대형 표본. r_sb(n)이 스케일 불변 상수면 진짜 상수,
드리프트하면 유한크기 효과. 오차 막대 = r/√(2N_samp).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

rng = np.random.default_rng(832)

print(f"{'스케일':>13} {'표본':>5} {'s평균':>8} {'r_sb':>7} {'±':>7}")
results = []
for SCALE, NS in [(1_000_000, 4000), (3_000_000, 2500),
                  (10_000_000, 2500), (30_000_000, 800)]:
    is_p = sieve(SCALE + 9_000_000)
    ps_all = primes_upto(SCALE // 2 + 4_500_000)
    y = int(round(SCALE ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]
    cands = [int(v) for v in np.arange(SCALE, SCALE + 8_000_000, 2)
             if v % 6 != 0][:NS]
    s_vals, S_sizes = [], []
    for n in cands:
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        s_vals.append(float(is_p[surv].mean()))
        S_sizes.append(len(surv))
    s = np.array(s_vals)
    S = np.array(S_sizes, dtype=float)
    binom = np.sqrt(s.mean() * (1 - s.mean()) / S.mean())
    r = s.std() / binom
    se = r / np.sqrt(2 * len(s))
    results.append((SCALE, r, se))
    print(f"{SCALE:>13,} {len(s):>5} {s.mean():>8.5f} {r:>7.4f} {se:>7.4f}",
          flush=True)
print(f"\n대조: √ln2 = {np.sqrt(np.log(2)):.4f} | 1 = 이항(독립)")
