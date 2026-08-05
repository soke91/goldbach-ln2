"""전선 B — 매끄러운 모듈러스에서의 c(n) 불일치 (BFI/Zhang 족 적용성).

q = 무제곱 y-매끄러운 합성수 (y = q^{1/4} 급, 소인수 ≤ 50), q ∈ (10³, 10⁵·⁵).
R_c(q) = 고정-잉여류(a = N mod q) 불일치/√기준 — 소수-q 결과(0.6~1.0)와
비교. 동등 건강이면 잘-인수분해 확장의 실측 근거.
"""

import math

import numpy as np

from goldbach.sieve import primes_upto

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

ps = primes_upto(X - 100)
ps = ps[ps > 2].astype(np.int64)
logp = np.log(ps.astype(np.float64))

# 매끄러운 무제곱 q 생성: 소인수 {7..47}, 3~4개 곱
small_ps = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
qs = set()
import itertools
for k in (3, 4):
    for combo in itertools.combinations(small_ps, k):
        q = 1
        for p in combo:
            q *= p
        if 1000 < q < 10 ** 5.5:
            qs.add(q)
qs = sorted(qs)
print(f"매끄러운 무제곱 q {len(qs)}개 ({qs[0]:,} ~ {qs[-1]:,})", flush=True)

NS = [99_999_998, 99_990_002]
for N in NS:
    muv = mu[N - ps]
    w = logp * muv
    tot = float(w.sum())
    v_unit = float((logp ** 2 * (muv != 0)).mean())
    nz = float((muv != 0).sum())
    rows = []
    for q in qs:
        if N % q == 0 or math.gcd(N, q) > 1:
            continue
        from math import gcd
        a = N % q
        sel = (ps % q == a)
        S_a = float(w[sel].sum())
        # 유효 잉여류 수: (Z/q)* 중 소수가 사는 류 = φ(q); a=N mod q가
        # (a, q)=1 이면 기대 = tot/φ(q)
        phi = q
        for p0 in small_ps:
            if q % p0 == 0:
                phi = phi // p0 * (p0 - 1)
        exp_a = tot / phi
        cnt = nz / phi
        rw = math.sqrt(max(cnt, 1.0) * v_unit)
        rows.append((q, math.log(q) / math.log(N), abs(S_a - exp_a) / rw))
    d = np.array(rows)
    print(f"\nN = {N:,}  매끄러운-q 표본 {len(d)}", flush=True)
    for lo, hi in [(0.35, 0.45), (0.45, 0.55), (0.55, 0.62), (0.62, 0.69)]:
        mk = (d[:, 1] >= lo) & (d[:, 1] < hi)
        if mk.sum():
            print(f"  θ {lo:.2f}-{hi:.2f}: R_c 평균 {d[mk,2].mean():.2f}  "
                  f"최악 {d[mk,2].max():.2f}  (표본 {int(mk.sum())})",
                  flush=True)
    np.savez(f"smooth_moduli_{N}.npz", rows=d)
print("전체완료", flush=True)
