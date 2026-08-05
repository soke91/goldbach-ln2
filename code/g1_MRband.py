"""G1 검증 엔진 — D2′ 완성-항등식의 실데이터 검산.

고정 (k,k′): C = Σ_{w≡c (mod L), w∈I} μ(w)μ(w′(w)),
  L-잉여 완성: C = (완전합의 잉여-투영) — 검산 ① 직접합 = 잉여-분해
  재합 (기계 정밀도). ② h-분해: 구간 지시자의 유한 푸리에 —
  주항(h=0) vs h≠0 위상항 크기 (√-상쇄 여부).
표본: (k,k′) 30쌍, k,k′ ~ 300..900, L = lcm ≤ 3×10⁵.
"""

import math

import numpy as np

X = 100_000_000

print("mu 계산...", flush=True)
mu = np.ones(X + 1, dtype=np.int8)
pm = np.ones(X + 1, dtype=bool)
pm[:2] = False
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        pm[p * p:: p] = False
        mu[p::p] *= -1
        mu[p * p:: p * p] = 0
val = np.arange(X + 1, dtype=np.int64)
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        val[p::p] //= p
        pp = p * p
        while pp <= X:
            val[pp::pp] //= p
            pp *= p
mu[val > 1] *= -1
mu[0] = 0
del val, pm
print("mu 완료됨", flush=True)

N = 99_999_998
rng = np.random.default_rng(191)
rows = []
t = 0
while t < 300:
    k = int(rng.integers(5000, 9999))
    kp = int(rng.integers(5000, 9999))
    if k == kp:
        continue
    g = math.gcd(k, kp)
    L = k * kp // g
    if L > 20_000_000:
        continue
    t += 1
    # w-합: w = N − pk, p 소수 ~ P — 여기선 단순화: p 정수 조건 대신
    # 정확 대상 재현: p ∈ 소수, w = N−pk, w′ = N−p·kp
    # (D2′의 w-표현과 동치: 직접 p-합으로 C 계산)
    P1 = min(110_000, (N - 2) // max(k, kp))
    P0 = P1 // 2
    if P1 - P0 < 5_000:
        t -= 1
        continue
    ps = np.arange(P0, P1, dtype=np.int64)
    # 소수 마스크 (에라토스: 상한 작아 즉석)
    sieve_p = np.ones(P1, dtype=bool)
    sieve_p[:2] = False
    for q in range(2, int(P1 ** 0.5) + 1):
        if sieve_p[q]:
            sieve_p[q * q:: q] = False
    ps = ps[sieve_p[P0:P1]]
    w = N - ps * k
    wp = N - ps * kp
    ok = (w > 1) & (wp > 1)
    C_direct = int((mu[w[ok]].astype(np.int16) * mu[wp[ok]]).sum())
    # 잉여-분해 재합 (검산 ①): p mod L 클래스별 부분합의 총합 = C
    resid = ps[ok] % L
    vals = (mu[w[ok]].astype(np.float64) * mu[wp[ok]])
    C_re = 0.0
    for a in np.unique(resid):
        C_re += vals[resid == a].sum()
    # h-분해 (검산 ②): p-구간 지시자 완성 — 주항 = (밀도)×전체,
    # h≠0 항 = 지수합. 여기선 h-항 크기의 직접 측정:
    # E_h = |Σ_p vals·e(h·p/L)| / √(#terms), h = 1..5
    n_terms = int(ok.sum())
    Eh = []
    for h in range(1, 6):
        ph = np.exp(2j * np.pi * h * ps[ok].astype(np.float64) / L)
        Eh.append(abs((vals * ph).sum()) / math.sqrt(n_terms))
    rows.append((k, kp, L, C_direct, abs(C_re - C_direct),
                 C_direct / math.sqrt(n_terms), max(Eh)))
    if t % 20 == 0:
        print(f"  {t}/30", flush=True)

d = np.array(rows)
print(f"\n[확장 검증, 쌍 300]")
print(f"검산① 잉여-재합 오차: 최대 {d[:,4].max():.1e} (0이어야 함)")
print(f"주항 규모 |C|/√n: 평균 {np.abs(d[:,5]).mean():.3f} 최악 "
      f"{np.abs(d[:,5]).max():.3f}")
print(f"h-위상항 |E_h|/√n (h=1..5 최대): 평균 {d[:,6].mean():.3f} 최악 "
      f"{d[:,6].max():.3f}")
print("→ h-항이 주항과 같은 √-급이면 완성 후 h-예산(x^o(1)개)로 통제 가능")
np.savez("g1_MRband.npz", rows=d)
print("전체완료", flush=True)
