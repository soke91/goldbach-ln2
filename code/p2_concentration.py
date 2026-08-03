"""비율 집중도 — s(n) = 파트너-생존자 중 소수 비율의 n-분포.

스케일별 (10⁶/10⁷/10⁸): 평균, 표준편차, 최소/최대, 그리고
'전원 P2 음모'(s=0)까지의 거리 = mean/std [σ].
집중 속도: std(s) ~ ? (생존자 수 S ~ n/ln²n 의 이항 요동이면 std ~ 1/√S)
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

rng = np.random.default_rng(59)

for SCALE, NSAMP in [(1_000_000, 300), (10_000_000, 150), (100_000_000, 40)]:
    is_p = sieve(SCALE + 300_000)
    ps_all = primes_upto(SCALE // 2 + 150_000)
    y = int(round(SCALE ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]
    cands = [int(v) for v in np.arange(SCALE, SCALE + 250_000, 2)
             if v % 6 != 0][:NSAMP]
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
    binom = float(np.mean(np.sqrt(s.mean() * (1 - s.mean()) / S)))
    print(f"n ~ {SCALE:>11,} ({NSAMP}표본): s = {s.mean():.4f} ± {s.std():.4f} "
          f"[이항예측 σ = {binom:.4f}] 범위 [{s.min():.4f}, {s.max():.4f}]")
    print(f"   '전원 P2'(s=0)까지 거리 = {s.mean()/s.std():.0f} σ | "
          f"생존자 평균 {S.mean():,.0f}", flush=True)
print(f"\n이론 평균 1/(1+ln2) = {1/(1+np.log(2)):.5f}")
print("""독해: std가 이항예측(√(s(1-s)/S))과 같으면 비율 요동은 순수 유한표본
효과 — 진짜 s(n)은 사실상 결정론적 상수. std ≫ 이항이면 n-구조 존재.""")
