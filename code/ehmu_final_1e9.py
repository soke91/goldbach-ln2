"""최종 좌표 스케일 사다리 — c(n) = Λ(n)μ(N−n) 지형 @ 10⁹.

μ는 세그먼트(10⁸ 단위) 계산으로 메모리 회피 (mu 전체 1GB int8만 유지).
N 2개, 소수 q 로그-표본 160개, θ = 0.35→0.70.
예상: μ-세그먼트 ~15분 + 소수생성 ~2분 + q-루프 ~5분 = ~25분. RAM ~2.5GB.
"""

import math

import numpy as np

from goldbach.sieve import primes_upto

X = 1_000_000_000
SEG = 100_000_000

print("기저 소수...", flush=True)
base = primes_upto(int(X ** 0.5) + 1)
mu = np.empty(X + 1, dtype=np.int8)

print("mu 세그먼트 계산...", flush=True)
for s0 in range(0, X + 1, SEG):
    s1 = min(s0 + SEG, X + 1)
    L = s1 - s0
    m_seg = np.ones(L, dtype=np.int8)
    val = np.arange(s0, s1, dtype=np.int64)
    for p in base:
        p = int(p)
        st = ((s0 + p - 1) // p) * p - s0
        m_seg[st::p] *= -1
        val[st::p] //= p
        pp = p * p
        while pp <= s1:
            st2 = ((s0 + pp - 1) // pp) * pp - s0
            m_seg[st2::pp] = 0
            sl = val[st2::pp]
            sl[sl % p == 0] //= p
            pp *= p
    big = val > 1
    m_seg[big] *= -1
    mu[s0:s1] = m_seg
    print(f"  세그 {s0:,}", flush=True)
mu[0] = 0
print("mu 완료됨. 검증 Σμ(≤1e6) =",
      int(mu[:1_000_001].astype(np.int64).sum()), "(참값 212)", flush=True)

print("소수 목록...", flush=True)
ps = primes_upto(X - 100)
ps = ps[ps > 2].astype(np.int64)
logp = np.log(ps.astype(np.float64))
print(f"소수 {len(ps):,}", flush=True)

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

NS = [999_999_998, 999_960_002]
targets = np.exp(np.linspace(math.log(2000), math.log(10 ** 6.3),
                             160)).astype(np.int64)

for N in NS:
    assert N % 6 == 2
    muv = mu[N - ps]
    w = logp * muv
    tot = float(w.sum())
    v_unit = float((logp ** 2 * (muv != 0)).mean())
    rows = []
    for i, t in enumerate(targets):
        q = int(t) | 1
        while not is_prime_u64(q):
            q += 2
        if N % q == 0:
            continue
        a = N % q
        sel = (ps % q == a)
        S_a = float(w[sel].sum())
        exp_a = tot / (q - 1)
        cnt = float((muv != 0).sum()) / (q - 1)
        rw = math.sqrt(max(cnt, 1.0) * v_unit)
        rows.append((q, math.log(q) / math.log(N), abs(S_a - exp_a) / rw))
        if (i + 1) % 40 == 0:
            print(f"  N={N}: {i+1}/160", flush=True)
            np.savez(f"ehmu_final9_{N}.npz", rows=np.array(rows))
    d = np.array(rows)
    print(f"\nN = {N:,}  (Σc = {tot:,.0f})", flush=True)
    for lo, hi in [(0.35, 0.45), (0.45, 0.50), (0.50, 0.55),
                   (0.55, 0.60), (0.60, 0.65), (0.65, 0.70)]:
        mk = (d[:, 1] >= lo) & (d[:, 1] < hi)
        if mk.sum():
            print(f"  θ {lo:.2f}-{hi:.2f}: R_c 평균 {d[mk,2].mean():.2f}  "
                  f"최악 {d[mk,2].max():.2f}  (표본 {int(mk.sum())})",
                  flush=True)
    np.savez(f"ehmu_final9_{N}.npz", rows=d)
print("전체완료", flush=True)
