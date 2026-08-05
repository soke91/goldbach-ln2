"""최종 좌표 직접 측정 — Huang–Li EH_mu 대상 c(n) = Λ(n)μ(N−n)의
고정-잉여류(a = N mod q) 불일치 지형, √N 장벽 전후.

R_c(q) = |Σ_{p≡N(q)} log(p)μ(N−p) − (1/φ(q))Σ_p log(p)μ(N−p)|
        / √(V_q),  V_q = (1/φ(q))Σ_p log²p·μ²(N−p)  (랜덤워크 규모)
θ = lnq/lnN 프로파일 0.30 → 0.70 (장벽 1/2 관통), N 3개.
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

def is_prime_u64(n):
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

NS = [99_999_998, 99_990_002, 99_960_002]   # N ≡ 2 (mod 6)
targets = np.exp(np.linspace(math.log(300), math.log(10 ** 5.6),
                             200)).astype(np.int64)

for N in NS:
    assert N % 6 == 2
    w = logp * mu[N - ps]                    # c(p) 가중
    tot = float(w.sum())
    v_unit = float((logp ** 2 * (mu[N - ps] != 0)).mean())
    rows = []
    for t in targets:
        q = int(t) | 1
        while not is_prime_u64(q):
            q += 2
        if N % q == 0:
            continue
        a = N % q
        sel = (ps % q == a)
        S_a = float(w[sel].sum())
        exp_a = tot / (q - 1)
        cnt = float((mu[N - ps] != 0).sum()) / (q - 1)
        rw = math.sqrt(max(cnt, 1.0) * v_unit)
        rows.append((q, math.log(q) / math.log(N), abs(S_a - exp_a) / rw))
    d = np.array(rows)
    print(f"\nN = {N:,}  (Σc = {tot:,.0f} — 자체도 상쇄 대상)", flush=True)
    for lo, hi in [(0.30, 0.40), (0.40, 0.45), (0.45, 0.50),
                   (0.50, 0.55), (0.55, 0.60), (0.60, 0.65), (0.65, 0.70)]:
        mk = (d[:, 1] >= lo) & (d[:, 1] < hi)
        if mk.sum():
            print(f"  θ {lo:.2f}-{hi:.2f}: R_c 평균 {d[mk,2].mean():.2f}  "
                  f"최악 {d[mk,2].max():.2f}  (표본 {int(mk.sum())})",
                  flush=True)
    np.savez(f"ehmu_final_{N}.npz", rows=d)
print("전체완료", flush=True)
