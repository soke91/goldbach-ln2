"""V1 결정 실험 — Λ² 상계의 실현 손실 vs 오라클-λ 손실.

g(n) ≤ Σ_{p≤n/2} (Σ_{d|n−p, d≤D} λ_d)²  (λ_1=1, 임의 λ — 무조건 유효)

① 고전 λ_d = μ(d)ln(D/d)/lnD 실현값/g — 예상 ~4 (이론 손실).
② 오라클 λ* = argmin λᵀG(n)λ (실데이터 그람, λ_1=1) — Λ² 족의
   데이터-최적 한계. 비율이 ~4 그대로면 족 소진, 유의미하게 낮으면
   채널이 상수를 움직일 수 있음 (V1 생존).
G_{d1d2}(n) = #{p ≤ n/2: lcm(d1,d2) | n−p}.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 2_000_000
D = int(X ** 0.4)
is_p = sieve(X + 10)
ps = primes_upto(X // 2 + 10)
ps = ps[ps > 2]

# 홀수 무제곱 d ≤ D
def squarefree_odds(D):
    out = []
    for d in range(1, D + 1, 2):
        m, ok = d, True
        for q in range(3, int(d ** 0.5) + 1, 2):
            if m % (q * q) == 0:
                ok = False
                break
        if ok:
            out.append(d)
    return out

ds = np.array(squarefree_odds(D))
K = len(ds)
lam_cl = np.array([1.0 if d == 1 else
                   mu_ * np.log(D / d) / np.log(D)
                   for d, mu_ in [(int(d), 0) for d in ds]])

# 뫼비우스
def mobius(d):
    if d == 1:
        return 1
    m, cnt = d, 0
    q = 3
    while q * q <= m:
        if m % q == 0:
            m //= q
            if m % q == 0:
                return 0
            cnt += 1
        q += 2
    if m > 1:
        cnt += 1
    return (-1) ** cnt

mus = np.array([mobius(int(d)) for d in ds])
lam_cl = mus * np.log(D / ds) / np.log(D)
lam_cl[ds == 1] = 1.0

from math import gcd
lcms = np.zeros((K, K), dtype=np.int64)
for i in range(K):
    for j in range(i, K):
        l = int(ds[i]) * int(ds[j]) // gcd(int(ds[i]), int(ds[j]))
        lcms[i, j] = lcms[j, i] = l

rng = np.random.default_rng(83)
base = X + (2 - X % 6) % 6 - 200_000
cands = sorted(int(base + 6 * k) for k in rng.integers(0, 30000, 12))

print(f"D = {D}, 무제곱 홀수 d {K}개, n {len(cands)}개 @ ~2e6", flush=True)
print(f"{'n':>9} {'g':>6} {'고전Λ²/g':>9} {'오라클Λ²/g':>10}")
res = []
for n in cands:
    pl = ps[ps <= n // 2]
    m = n - pl
    # 그람: G_ij = #{m ≡ 0 mod lcm}
    G = np.zeros((K, K))
    # d별 residue 벡터 캐시로 lcm 카운트: 직접 lcm별 계산 (중복 lcm 캐시)
    cache = {}
    for i in range(K):
        for j in range(i, K):
            l = int(lcms[i, j])
            if l not in cache:
                if l > n:
                    cache[l] = 0
                else:
                    cache[l] = int(np.count_nonzero(m % l == 0))
            G[i, j] = G[j, i] = cache[l]
    g_true = int(is_p[m].sum())
    b_cl = float(lam_cl @ G @ lam_cl)
    # 오라클: min λᵀGλ, λ_1 = 1 → 블록 해: λ_rest = -G_rr^{-1} G_r1
    idx1 = int(np.where(ds == 1)[0][0])
    rest = [i for i in range(K) if i != idx1]
    Grr = G[np.ix_(rest, rest)] + 1e-9 * np.eye(K - 1)
    Gr1 = G[rest, idx1]
    lam_r = np.linalg.solve(Grr, -Gr1)
    lam_o = np.zeros(K)
    lam_o[idx1] = 1.0
    lam_o[rest] = lam_r
    b_or = float(lam_o @ G @ lam_o)
    res.append((n, g_true, b_cl / g_true, b_or / g_true))
    print(f"{n:>9} {g_true:>6} {b_cl/g_true:>9.3f} {b_or/g_true:>10.3f}",
          flush=True)

r = np.array(res)
print(f"\n고전 Λ²/g: 평균 {r[:,2].mean():.3f} ± {r[:,2].std():.3f}")
print(f"오라클 Λ²/g: 평균 {r[:,3].mean():.3f} ± {r[:,3].std():.3f}")
print("오라클이 ~1에 가까울수록 Λ² 족 내 개선 여지 큼 / 4 근처면 족 소진")
np.savez("v1_lambda2.npz", res=r)
