"""유도 지원 수치 — 균형-예산 삼중분해 밀도.

k = qrs, 각 인자 ∈ [K^{1/3}/C, C·K^{1/3}] (C = 4, 16) 가능한 k의 밀도.
게이트 산술의 R~S~K^{1/3} 예산이 족의 몇 %에서 실제 가용한지 +
불가 잔여의 구조 (소수/준소수/거대인수) 분류.
K-다이애딕 3대역: 10⁴, 10⁵, 10⁶.
"""

import numpy as np

KMAX = 2_000_000
print("spf 체...", flush=True)
spf = np.zeros(KMAX + 1, dtype=np.int64)
for p in range(2, int(KMAX ** 0.5) + 1):
    if spf[p] == 0:
        spf[p * p:: p] = np.where(spf[p * p:: p] == 0, p, spf[p * p:: p])
        spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
for i in range(2, KMAX + 1):
    if spf[i] == 0:
        spf[i] = i
print("완료", flush=True)

def factors(k):
    fs = []
    while k > 1:
        p = int(spf[k])
        fs.append(p)
        k //= p
    return fs

def balanced3(k, C):
    t = round(k ** (1 / 3))
    lo, hi = t / C, t * C
    fs = factors(k)
    if len(fs) < 3:
        return False
    # 탐욕/전수: 인수들을 3그룹 곱으로 분할, 각 곱 ∈ [lo, hi]
    from itertools import product
    n = len(fs)
    if n > 12:
        return True  # 충분히 많은 인수는 균형 분할 사실상 항상 가능
    for assign in product(range(3), repeat=n):
        g = [1, 1, 1]
        for f, a in zip(fs, assign):
            g[a] *= f
        if all(lo <= x <= hi for x in g):
            return True
    return False

rng = np.random.default_rng(179)
for Kc in [10_000, 100_000, 1_000_000]:
    ks = rng.integers(Kc, 2 * Kc, 4000)
    for C in [4, 16]:
        ok = sum(1 for k in ks if balanced3(int(k), C))
        print(f"K~{Kc:,}  C={C}: 균형-삼중분해 밀도 {ok/len(ks):.3f}",
              flush=True)
    # 잔여 구조
    n_pr = sum(1 for k in ks if len(factors(int(k))) == 1)
    n_sp = sum(1 for k in ks if len(factors(int(k))) == 2)
    print(f"  잔여 참조: 소수 {n_pr/len(ks):.3f}, 2-인수 {n_sp/len(ks):.3f}",
          flush=True)
print("전체완료", flush=True)
