"""A1+A5 — s(n) 최악-포락선 조밀 사냥 + 저-s 산술 프로파일.

10⁶ 스케일: 창 내 n≡2(6) 전수(3만 개). min-s 포락선, 하위 0.1% n들의
산술 특징(𝔖(n), 소인수 구조) vs 전체 — 저-s의 인상착의.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 1_000_000
W = 180_000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
print(f"전수 {len(cands):,}개 (n≡2 mod 6, [{X:,}, {X+W:,}])", flush=True)

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

s = s_vals
mu, sd = s.mean(), s.std()
print(f"\ns: 평균 {mu:.5f} ± {sd:.5f}")
print(f"전수 최소 {s.min():.5f} ({(mu-s.min())/sd:.2f}σ) / 최대 {s.max():.5f}")
print(f"s=0까지 {mu/sd:.0f}σ | 표본수 {len(s):,} — 가우시안 기대 최심 "
      f"{np.sqrt(2*np.log(len(s))):.2f}σ")

# 저-s 하위 0.2% 프로파일
k = max(20, len(s) // 500)
idx = np.argsort(s)[:k]
low_n = np.array(cands)[idx]
Sg_low = np.array([sing0(int(n)) for n in low_n])
Sg_all = np.array([sing0(int(n)) for n in np.random.default_rng(1).choice(cands, 500)])
print(f"\n[저-s 하위 {k}개 프로파일]")
print(f"  𝔖₀ 평균: 저-s {Sg_low.mean():.4f} vs 전체 {Sg_all.mean():.4f}")
mod30 = np.bincount(low_n % 30, minlength=30)
top = np.argsort(mod30)[::-1][:4]
print(f"  n mod 30 최빈: {[(int(t), int(mod30[t])) for t in top]}")
