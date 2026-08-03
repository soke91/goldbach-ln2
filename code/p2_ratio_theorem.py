"""60% 정리 검증 — 소수:P2 = 1 : ln 2.

유도: 생존자(spf > x^{1/3}, m ~ x) 중
  소수 밀도 1/ln x, P2 밀도 = ∫_{1/3}^{1/2} du/(u(1-u))/ln x = ln2/ln x
  ⇒ 소수 비율 = 1/(1+ln2) = 0.59064

검증 2종:
  A. 비조건 m (무작위 홀수, m ~ X): 비율 → 0.5906 예측 (유도의 직접 검증)
  B. 골드바흐 파트너 m = n−p (기존 H4 측정 60.4~60.7%):
     초과분 = 쌍상관(𝔖) 조건화 보정 — A와의 차이로 정량화
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 100_000_000
is_p = sieve(X + 200_000)
y = int(round(X ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y)]
print(f"X = {X:,}, y = X^(1/3) ≈ {y}")

rng = np.random.default_rng(60)

# ── A. 비조건 무작위 홀수 ────────────────────────────────────
m = np.arange(X - 4_000_000 + 1, X, 2)  # 홀수 구간 (2백만 개)
alive = np.ones(len(m), dtype=bool)
for q in qs:
    if q == 2:
        continue
    alive &= (m % q != 0)
surv = m[alive]
prime_frac_A = float(is_p[surv].mean())
print(f"\nA. 비조건 생존자 {len(surv):,}개")
print(f"   소수 비율 = {prime_frac_A:.4f}   [예측 1/(1+ln2) = {1/(1+np.log(2)):.4f}]")

# ── B. 골드바흐 파트너 (얇은 띠 n 표본 20개) ─────────────────
ps_all = primes_upto(X // 2 + 100_000)
cands = [int(v) for v in np.arange(X, X + 40_000, 2) if v % 6 != 0][:20]
fr = []
for n in cands:
    mm = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(mm), dtype=bool)
    for q in qs:
        if q > 2 and n % q:
            al &= (mm % q != 0)
    al &= mm > 1
    S = mm[al]
    fr.append(float(is_p[S].mean()))
prime_frac_B = float(np.mean(fr))
print(f"\nB. 파트너 생존자 (n~10^8, 20표본)")
print(f"   소수 비율 = {prime_frac_B:.4f}")
print(f"\n조건화 보정(B−A) = {prime_frac_B - prime_frac_A:+.4f}")
print(f"이론 상수: 1/(1+ln2) = {1/(1+np.log(2)):.5f} | ln2/(1+ln2) = {np.log(2)/(1+np.log(2)):.5f}")
