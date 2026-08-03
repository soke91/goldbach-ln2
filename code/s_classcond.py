"""A5-2 — 계급-조건부 s 분석.

가설: s(n)의 산포는 𝔖₀(n)이 대부분 설명한다. 𝔖₀ 고정 계급 내 잔차 산포가
전체 산포보다 크게 작으면, 최악-케이스 통제가 얇은-띠 구조로 환원된다.

산출: corr(s, 𝔖₀), 회귀 후 잔차 σ 비교, mod-30 계급별 (평균, σ, min),
계급별 min의 계급-σ 기준 깊이.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 1_000_000
W = 180_000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]

cands = np.array([n for n in range(X, X + W, 2) if n % 6 == 2])
print(f"전수 {len(cands):,}개", flush=True)

def sing0(n):
    c, m_, p = 1.0, n, 3
    while m_ % 2 == 0: m_ //= 2
    while p * p <= m_:
        if m_ % p == 0:
            c *= (p - 1) / (p - 2)
            while m_ % p == 0: m_ //= p
        p += 2
    if m_ > 1: c *= (m_ - 1) / (m_ - 2)
    return c

s_vals = np.empty(len(cands))
for i, n in enumerate(cands):
    n = int(n)
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    s_vals[i] = float(is_p[surv].mean())
    if (i + 1) % 5000 == 0:
        print(f"  {i+1}/{len(cands)}", flush=True)

np.savez("s_classcond_data.npz", cands=cands, s=s_vals)
Sg = np.array([sing0(int(n)) for n in cands])
s = s_vals
mu, sd = s.mean(), s.std()

# 1) s ~ 𝔖₀ 상관/회귀
r = float(np.corrcoef(s, Sg)[0, 1])
A = np.vstack([np.ones_like(Sg), Sg]).T
coef, *_ = np.linalg.lstsq(A, s, rcond=None)
resid = s - A @ coef
print(f"\ncorr(s, 𝔖₀) = {r:+.4f}")
print(f"전체 σ = {sd:.5f} → 𝔖₀-회귀 잔차 σ = {resid.std():.5f} "
      f"(설명된 분산 {1 - resid.var()/s.var():.1%})")
print(f"회귀: s = {coef[0]:.4f} + {coef[1]:+.4f}·𝔖₀")

# 로그-회귀도 (곱구조 기대): s ~ a + b·ln𝔖₀
A2 = np.vstack([np.ones_like(Sg), np.log(Sg)]).T
c2, *_ = np.linalg.lstsq(A2, s, rcond=None)
r2 = s - A2 @ c2
print(f"ln𝔖₀-회귀 잔차 σ = {r2.std():.5f} (설명 {1 - r2.var()/s.var():.1%}), "
      f"b = {c2[1]:+.4f}")

# 2) mod-30 계급별
print(f"\n[mod 30 계급별]  전체: 평균 {mu:.5f} σ {sd:.5f} min {s.min():.5f}")
print(f"{'류':>4} {'개수':>6} {'평균':>8} {'σ':>8} {'min':>8} {'min깊이σ':>8}")
for cls in sorted(set(int(n) % 30 for n in cands)):
    mask = (cands % 30) == cls
    sc = s[mask]
    depth = (sc.mean() - sc.min()) / sc.std()
    print(f"{cls:>4} {mask.sum():>6} {sc.mean():>8.5f} {sc.std():>8.5f} "
          f"{sc.min():>8.5f} {depth:>8.2f}")

# 3) 잔차 기준 최악 재평가
worst_res = resid.min()
print(f"\n잔차 최악: {worst_res:+.5f} ({-worst_res/resid.std():.2f}σ_res)  "
      f"| 잔차 가우시안 기대 최심 {np.sqrt(2*np.log(len(s))):.2f}σ")
