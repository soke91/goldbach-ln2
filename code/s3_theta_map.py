"""S3 — 거칠기 보간 장벽 지도.

y = X^theta 체질 후 생존자의 소수 비율 s_theta 실측 → price(theta) =
1/(1-s) (적: 비소수 생존자). 고전 선형 체: L = f(1/(2*theta)) (BV 수준
D = X^{1/2}), U = 2 (Selberg 관례). 격차 gamma = (U/L)/price — 1 미만이면
승리. 지도에서 최소 격차 지점 = 공략 최적 theta.

f(s) = 2 e^gamma ln(s-1)/s (2 <= s <= 4), s < 2 에서 0.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

EG = 0.5772156649015329
X = 2_000_000
W = 100_000
M = 800
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
rng = np.random.default_rng(3)
cands = sorted(rng.choice(cands, M, replace=False).tolist())

print(f"X={X:,}  표본 {M}  (U=2, D=X^0.5 관례)")
print(f"{'theta':>6} {'s평균':>8} {'s_min':>8} {'price':>7} {'L=f(s)':>8} "
      f"{'U/L':>7} {'gamma':>8}")

for th in [0.20, 0.22, 0.24, 0.25, 0.27, 0.30, 0.333]:
    y = X ** th
    qs = [int(q) for q in primes_upto(int(y) + 1) if q > 2]
    s_vals = np.empty(M)
    for i, n in enumerate(cands):
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        s_vals[i] = float(is_p[surv].mean())
    smean, smin = s_vals.mean(), s_vals.min()
    price = 1 / (1 - smean)
    s_sieve = 0.5 / th
    L = 2 * np.exp(EG) * np.log(s_sieve - 1) / s_sieve if s_sieve > 2 else 0.0
    UL = 2 / L if L > 0 else float("inf")
    gam = UL / price
    print(f"{th:>6.3f} {smean:>8.5f} {smin:>8.5f} {price:>7.4f} {L:>8.4f} "
          f"{UL:>7.3f} {gam:>8.3f}", flush=True)

print("\n비고: gamma < 1 이 승리 조건. 고전 상수로는 전 구간 >1 (패리티"
      " 장벽)이 예상 — 지도의 목적은 최소 격차 지점과 '필요 개선 배율' 특정.")
