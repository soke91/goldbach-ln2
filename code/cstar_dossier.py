"""(C*) 실증 서류 — 비가중/로그-가중 상쇄 비교 + 계수-크기 균일성.

쌍 (k, k') 로그-표본, k-크기 3다이애드 (N^0.2, N^0.3, N^0.45):
  r_unw = |Σ_m μμ| / √M_eff,   r_log = |Σ_m μμ/m| / √(Σ 1/m²... 규모)
비교: 두 상쇄가 동행하면 T1의 자연-측 간극 없음; r이 k-크기에 무관하면
T2의 자연-측 간극 없음. (정리의 간극만 남음을 확정하는 서류)
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
rng = np.random.default_rng(149)
bands = [(0.18, 0.22), (0.28, 0.32), (0.43, 0.47)]
print(f"{'대역':>12} {'쌍':>4} {'r_unw 평균':>10} {'r_log 평균':>10} "
      f"{'unw 최악':>8} {'log 최악':>8}")
all_rows = []
for lo, hi in bands:
    rows = []
    for t in range(120):
        k1 = int(rng.integers(int(X ** lo), int(X ** hi)))
        k2 = int(rng.integers(int(X ** lo), int(X ** hi)))
        if k1 == k2:
            continue
        M = (N - 1) // max(k1, k2)
        ms = np.arange(1, M + 1, dtype=np.int64)
        prod = (mu[N - k1 * ms].astype(np.int16) * mu[N - k2 * ms]
                ).astype(np.float64)
        M_eff = float(np.count_nonzero(prod))
        r_unw = abs(prod.sum()) / math.sqrt(max(M_eff, 1))
        w = 1.0 / ms
        # 로그-가중 규모: √(Σ w² · 밀도) — 반정규 비교 가능 규격화
        scale = math.sqrt(float((w * w * (prod != 0)).sum()))
        r_log = abs(float((prod * w).sum())) / max(scale, 1e-12)
        rows.append((k1, k2, r_unw, r_log))
    d = np.array(rows)
    all_rows.append(d)
    print(f"N^{lo:.2f}-{hi:.2f} {len(d):>4} {d[:,2].mean():>10.3f} "
          f"{d[:,3].mean():>10.3f} {d[:,2].max():>8.2f} {d[:,3].max():>8.2f}",
          flush=True)
np.savez("cstar_dossier.npz", b1=all_rows[0], b2=all_rows[1], b3=all_rows[2])
print("반정규 기준 0.798 — 두 열이 모두 ~0.8이고 대역 간 평탄하면")
print("T1(가중)·T2(균일성)의 자연-측 간극 부재 확정")
print("전체완료", flush=True)
