"""S1-2 — 공통-모드(랭크-1) 제거 후 교차-섬유 잔차 상관.

모델: F_q(n) ~ a_q * b_n (곱 구조 = 공통 𝔖-모드). IPF 2회로 적합 후
z = F/(a_q b_n) 잔차의 교차상관. 음이면 예산 맞교환(전역 상계 슬롯),
0이면 섬유 독립(정보 없음), 양이면 미제거 구조 잔존.
"""

import numpy as np

d = np.load("s1_fiber_data.npz")
F, qs = d["F"], d["qs"]
M, K = F.shape

a = F.mean(axis=0)
for _ in range(30):
    b = (F / a).mean(axis=1)
    a = (F / b[:, None]).mean(axis=0)
model = np.outer(b, a)
z = F / model
print(f"랭크-1 적합: 설명 분산 {1 - (F-model).var()/F.var():.2%}")

C = np.corrcoef(z.T)
off = C[np.triu_indices_from(C, k=1)]
print(f"잔차 교차-섬유 상관: 평균 {off.mean():+.5f}  중앙값 "
      f"{np.median(off):+.5f}  범위 [{off.min():+.3f}, {off.max():+.3f}]")

# 잔차 총합 분산 vs 독립 기대 (음상관이면 <1)
tot_z = (z * a).sum(axis=1)
G_res = tot_z.var() / (z.var(axis=0) * a ** 2).sum()
print(f"잔차 집중 이득 G_res = {G_res:.3f}  (독립 1, 맞교환 <1)")

# 섬유 잔차의 포아송 비교: Var(F_q)/mean(F_q) 조건부(모델 제거 후)
fano = ((F - model) ** 2).mean(axis=0) / F.mean(axis=0)
print(f"섬유 Fano(잔차): 평균 {fano.mean():.3f}  중앙값 "
      f"{np.median(fano):.3f}  (포아송 1, sub-poisson <1)")
lo, hi = fano[: K // 3].mean(), fano[-K // 3:].mean()
print(f"  q-하위⅓ {lo:.3f} vs 상위⅓ {hi:.3f}")
