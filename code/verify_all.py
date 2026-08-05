"""재현성 종합검증 — 핵심 엔진 축약판 일괄 실행 (공개 리포용 CI-도장).

각 검증의 축약판 (표본 축소)을 순차 실행, 요약표 출력.
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
del val
print("mu 완료됨 (검증: Σμ(≤1e6) =",
      int(mu[:1_000_001].astype(np.int64).sum()), ")", flush=True)

N = 99_999_998
SQ = int(N ** 0.5)
rng = np.random.default_rng(211)
rows = []

# V1: 핵 총합 T(N) 반정규 (30 N)
ps_small = np.nonzero(pm[:X // 2])[0]
ps_small = ps_small[ps_small > 2].astype(np.int64)
logp = np.log(ps_small.astype(np.float64))
rs = []
for N2 in sorted(rng.integers(N - 3_000_000, N, 30) // 6 * 6 + 2):
    N2 = int(N2)
    muv = mu[N2 - ps_small[ps_small < N2 - 2]].astype(np.float64)
    lp = logp[: len(muv)]
    T = float((lp * muv).sum())
    V = float((lp ** 2 * (muv != 0)).sum())
    rs.append(abs(T) / math.sqrt(V))
rows.append(("핵 T(N) 반정규 (30N)", f"평균 r={np.mean(rs):.3f} (기준 0.80)"))
print(rows[-1], flush=True)

# V2: 사다리 항등식 (10쌍, p=3)
errs = []
for _ in range(10):
    k = int(rng.integers(500, 1200))
    ms = np.arange(SQ + 1, (N - 1) // k + 1, dtype=np.int64)
    sel = ms[ms % 3 == 0]
    A = int((mu[sel].astype(np.int16) * mu[N - k * sel]).sum(dtype=np.int64))
    m2 = np.arange(SQ // 3 + 1, (N - 1) // (3 * k) + 1, dtype=np.int64)
    m2 = m2[m2 % 3 != 0]
    B = int((mu[m2].astype(np.int16) * mu[N - 3 * k * m2]).sum(dtype=np.int64))
    errs.append(abs(A + B))
rows.append(("사다리 항등식 (10쌍)", f"최대 오차 {max(errs)} (0이어야)"))
print(rows[-1], flush=True)

# V3: 분산 엔진 (50쌍)
K0 = int(X ** 0.4)
rr = []
for _ in range(50):
    k1, k2 = int(rng.integers(K0 // 2, K0)), int(rng.integers(K0 // 2, K0))
    if k1 == k2:
        continue
    M = (N - 1) // max(k1, k2)
    msv = np.arange(1, M + 1, dtype=np.int64)
    prod = (mu[N - k1 * msv].astype(np.int16) * mu[N - k2 * msv]).astype(float)
    v = int(np.count_nonzero(prod))
    rr.append(abs(prod.sum()) / math.sqrt(max(v, 1)))
rows.append(("분산 비대각 (50쌍)", f"평균 r={np.mean(rr):.3f} (반정규 0.80)"))
print(rows[-1], flush=True)

# V4: 이음새 대역 (20쌍)
rr2 = []
for _ in range(20):
    k = int(rng.integers(252, 464))
    kp = int(rng.integers(252, 464))
    if k == kp:
        continue
    P1 = min(110_000, (N - 2) // max(k, kp))
    P0 = P1 // 2
    ps2 = np.arange(P0, P1, dtype=np.int64)
    ps2 = ps2[pm[P0:P1]]
    w = N - ps2 * k
    wp = N - ps2 * kp
    ok = (w > 1) & (wp > 1)
    vals = mu[w[ok]].astype(np.float64) * mu[wp[ok]]
    v = int(np.count_nonzero(vals))
    rr2.append(abs(vals.sum()) / math.sqrt(max(v, 1)))
rows.append(("이음새 대역 (20쌍)", f"평균 r={np.mean(rr2):.3f} (반정규 0.80)"))
print(rows[-1], flush=True)

# V5: 인수분해 법칙 (144~154차분) — 자유-클래스 가우시안 + 마스크 블라인드
from math import gcd
K0f, K1f = 2000, 4000
P0f = N // (2 * K1f)
ps_f = np.nonzero(pm[P0f:2 * P0f])[0] + P0f
m2s = []
n_free = 0
tries = 0
while n_free < 400 and tries < 20000:
    tries += 1
    k1 = int(rng.integers(K0f, K1f)); k2 = int(rng.integers(K0f, K1f))
    if k1 == k2 or gcd(k1 * k2, N) != 1:
        continue
    ppf = ps_f[ps_f <= (N - 2) // max(k1, k2)]
    if len(ppf) < 200:
        continue
    t = mu[N - ppf * k1].astype(np.int64) * mu[N - ppf * k2]
    nz = int(np.count_nonzero(t))
    if nz > 50:
        m2s.append(float(t.sum()) ** 2 / nz)
        n_free += 1
m2f = float(np.mean(m2s))
rows.append(("인수분해 법칙: 자유-클래스 분산비",
             f"m2_eff={m2f:.3f} (기준 1.0 ± 0.1, {n_free}쌍)"))
print(rows[-1], flush=True)

# V6: E1 비율 미니 도장 (100 k, 대역 [3000,4000))
ks6 = rng.choice(np.arange(3000, 4000), size=100, replace=False)
S2 = 0.0; SS = 0.0
for k6 in ks6:
    k6 = int(k6)
    ms6 = np.arange(SQ + 1, N // k6 + 1, dtype=np.int64)
    t6 = mu[ms6].astype(np.int64) * mu[N - k6 * ms6]
    S2 += float(t6.sum()) ** 2
    SS += int(np.count_nonzero(t6))
rows.append(("E1 비율 (100k, 참고용 — 실SE ~0.15)",
             f"{S2/SS:.3f} (기준 1.0, [0.7,1.4] 2σ-호환)"))
print(rows[-1], flush=True)

print("\n===== 종합검증 요약 =====")
for name, res in rows:
    print(f"  {name}: {res}")
print("전체완료", flush=True)
