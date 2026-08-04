"""추측 P의 직접 측정 — 체질된-소수 카운트의 s-프로파일.

n 표본 @ 10^8. z = (n/2)^beta 격자에서
  S_z(n) = #{p ≤ n/2: (n−p, P(z)) = 1, n−p > z}
모형: pi(n/2) · ∏_{2<q≤z, q∤n}(q−2)/(q−1) · (법칙 인자는 비율에서 상쇄)
— 정밀 모형 대신 **Buchstab-랜덤 기준**: 같은 밀도의 랜덤 홀수 집합을
같은 z로 체질한 기대치 대비 비율. 더 깨끗한 것: Mertens-곱 모형
  M_z(n) = pi(n/2) · ∏_{2<q≤z, q∤n} (q−2)/(q−1)
R_P(s) = S_z / M_z.  s = ln(n^{1/2})/ln z 관례로 s-축 표기.
P 예측: R_P → 1 (수정: 쌍상관 보정 후). s < 2에서의 거동이 관건.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 100_000_000
W = 2_000_000
M = 60
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)

rng = np.random.default_rng(97)
cands = sorted(int(X + (2 - X % 6) % 6 + 6 * k)
               for k in rng.integers(0, W // 6, M))

betas = [1/8, 1/6, 1/5, 1/4, 0.30, 1/3, 0.40, 0.45, 0.50]
print(f"n {len(cands)}개 @ 10^8;  s = 0.5/beta 관례", flush=True)
print(f"{'beta':>6} {'s':>6} {'R_P 평균':>9} {'±':>7} {'f(s)':>6} {'F(s)':>6}")

EG = 0.5772156649015329
def F_up(s):
    if s <= 0: return np.inf
    if s <= 3: return 2*np.exp(EG)/s
    return max(1.0, 2*np.exp(EG)/s*(1+0.06*(s-3)))
def f_low(s):
    if s < 2: return 0.0
    if s <= 4: return 2*np.exp(EG)*np.log(s-1)/s
    return min(1.0, 0.978+0.01*(s-4))

for beta in betas:
    ratios = []
    for n in cands:
        half = n // 2
        z = int(round(half ** beta))
        qs = [int(q) for q in primes_upto(z + 1) if q > 2]
        pl = ps_all[ps_all <= half]
        m = n - pl
        al = m > z
        prod = 1.0
        for q in qs:
            if n % q:
                al &= (m % q != 0)
                prod *= (q - 2) / (q - 1)
        S_z = int(al.sum())
        Mz = len(pl) * prod
        ratios.append(S_z / Mz)
    r = np.array(ratios)
    s = 0.5 / beta
    print(f"{beta:>6.3f} {s:>6.2f} {r.mean():>9.4f} {r.std():>7.4f} "
          f"{f_low(s):>6.3f} {F_up(s):>6.3f}", flush=True)
    np.savez(f"p_profile_{int(beta*1000)}.npz", r=r,
             cands=np.array(cands))
print("\nR_P = 실측/Mertens-모형. P 예측: 전 s에서 상수(쌍상관 보정치).")
print("[f, F] = 체 이론 허용 구간 — R_P가 구간 어디에 있는지가 P의 지도.")
