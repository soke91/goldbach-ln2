"""정면돌파 측정 — 타입 II의 모듈러스-대역 질량 프로파일.

B = Σ_{m>U} a_m Σ_{V<k<N/m} μ(N−mk).  내부합 inner(m) = Σ_k μ(N−mk).
쌍곡선 분할(m ≤ √N vs k ≤ √N 듀얼)에서 각 다이애딕 m-대역의
질량 Σ_m |inner(m)| 프로파일 — √N 경계 대역이 전체의 몇 %인가.
+ inner(m)의 상쇄 통계 (r = |inner|/√K_eff — 여기도 반정규인가).
m은 대역당 표본 추출로 추정 (전수 불필요).
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
U = int(N ** (1 / 3))
V = U
SQ = int(N ** 0.5)
rng = np.random.default_rng(151)

# k-면 (듀얼): dual(k) = Σ_{√N<m<N/k} μ(m)μ(N−mk) — 파리티 잔여의 실측
print(f"[k-면 듀얼] U=V={U:,}  √N={SQ:,}")
print(f"{'대역 (m~)':>16} {'표본평균|inner|':>14} {'평균 r':>8} "
      f"{'대역질량추정':>14} {'누적%':>7}")
bands = []
lo = SQ // 3
while lo < SQ:
    hi = min(lo * 2, SQ)
    ks_smp = sorted(set(int(v) for v in rng.integers(lo, hi, 60)))
    absr = []
    rr = []
    for k in ks_smp:
        M0 = SQ + 1
        M1 = (N - 1) // k
        if M1 <= M0:
            continue
        ms_arr = np.arange(M0, M1 + 1, dtype=np.int64)
        vals = mu[ms_arr].astype(np.int16) * mu[N - k * ms_arr]
        inner = int(vals.sum(dtype=np.int64))
        K_eff = int(np.count_nonzero(vals))
        absr.append(abs(inner))
        rr.append(abs(inner) / math.sqrt(max(K_eff, 1)))
    if absr:
        bands.append((lo, hi, float(np.mean(absr)), float(np.mean(rr)),
                      (hi - lo) * float(np.mean(absr))))
    lo = hi

tot = sum(b[4] for b in bands)
cum = 0.0
for lo, hi, ma, mr, mass in bands:
    cum += mass
    print(f"{lo:>7,}-{hi:>7,} {ma:>14.1f} {mr:>8.2f} {mass:>14.3e} "
          f"{100*cum/tot:>6.1f}%", flush=True)
top = bands[-1]
print(f"\n총질량(추정) Σ_m|inner| ≈ {tot:.3e}  vs N = {N:.1e}")
print(f"최상 대역(√N 경계) 점유: {100*top[4]/tot:.1f}%")
print(f"r(내부 상쇄): 전 대역 평균 "
      f"{np.mean([b[3] for b in bands]):.2f} (반정규 0.798 기준)")
print(f"절대합의 N-대비: {tot/N:.4f} — polylog 절약이면 승리 구도")
np.savez("typeII_kface_MR.npz", bands=np.array(bands))
print("전체완료", flush=True)
