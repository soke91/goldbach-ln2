"""전선2 수치 — 보완 항등식 성분 실측 @ 10^6 (전체 배열, 빠름).

core(n) = S(n) − g(n) − mid(n).  (S = 생존자 총수, g = 소수 파트너,
mid = spf ≤ x^{1/6} 합성, core = spf > x^{1/6})
theta = 1/8 아날로그가 10^6에선 y = 5.6 (전-점근) — 대신 y = 소수 {3,5}
수준이라 왜곡. 그래서 여기선 u-변환: y = x^{1/6} = 10, 코어 = spf > x^{1/3}
구도로 동일 대수 구조를 낮은 u에서 검증 (성분 항등식과 교환율 산술 확인).

교환율 곡선: core ≤ S − g_lower − mid_lower 에서
필요 (이진 하한 γ, 준소수 하한 λ) 조합의 실측 격자.
"""

import numpy as np

import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(
        _os.path.join(_d, 'lib', 'goldbach')):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, 'lib'))
from goldbach.sieve import primes_upto, sieve

X = 1_000_000
W = 60_000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)

y = 10          # x^{1/6}-급
x13 = 100       # 코어 경계 spf > 100 (= x^{1/3})
small = [3, 5, 7]
med = [int(p) for p in primes_upto(x13 + 1) if p > 7]

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
rows = []
for n in cands:
    m = n - ps_all[ps_all <= n // 2]
    al = m > 1
    for q in small:
        if n % q:
            al &= (m % q != 0)
    surv = m[al]
    S = len(surv)
    pr = is_p[surv]
    g = int(pr.sum())
    comp = surv[~pr]
    mid_mask = np.zeros(len(comp), dtype=bool)
    for q in med:
        mid_mask |= (comp % q == 0)
    mid = int(mid_mask.sum())
    core = len(comp) - mid
    rows.append((n, S, g, mid, core))

d = np.array(rows, dtype=np.float64)
S, g, mid, core = d[:, 1], d[:, 2], d[:, 3], d[:, 4]
print(f"n {len(d)}개 @ 10^6, y=10(체질 {small}), 코어 spf>{x13}")
print(f"성분 평균 (S 정규화): g {np.mean(g/S):.4f}  mid {np.mean(mid/S):.4f}"
      f"  core {np.mean(core/S):.4f}")
print(f"항등식 검산 S = g+mid+core: 최대 오차 "
      f"{np.max(np.abs(S - g - mid - core)):.0f}")

# 교환율: core ≤ S − γ·g_model − λ·mid_model 이 1.74·core_model 이하가
# 되기 위한 (γ, λ) 격자 — 모델 = 실측 평균비 사용 (경험 프런티어)
gm, mm, cm = np.mean(g / S), np.mean(mid / S), np.mean(core / S)
print(f"\n[교환율 프런티어: 필요 조건 1 − γ·{gm:.3f} − λ·{mm:.3f} "
      f"≤ 1.74×{cm:.3f} = {1.74*cm:.3f}]")
print(f"{'λ (준소수 하한 질)':>18} {'필요 γ (이진 하한 질)':>20}")
for lam in [1.0, 0.9, 0.8, 0.7, 0.5, 0.3, 0.0]:
    gam = (1 - 1.74 * cm - lam * mm) / gm
    print(f"{lam:>18.2f} {gam:>20.3f}")
print("\nγ ≤ 0 이면 이진 하한 불필요 (준소수만으로 코어 상계 달성)")
np.savez("complement_identity.npz", data=d)
