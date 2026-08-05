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


rng = np.random.default_rng(131)
base = X - 50_000_000
NS = sorted(set(int(base + 6 * k + (2 - base % 6) % 6)
                for k in rng.integers(0, 8_000_000, 120)))
logp2 = logp ** 2
rows = []
for i, N in enumerate(NS):
    muv = mu[N - ps].astype(np.float64)
    T = float((logp * muv).sum())
    V = float((logp2 * (muv != 0)).sum())
    rows.append((N, T, abs(T) / math.sqrt(V)))
    if (i + 1) % 12 == 0:
        print(f"  {i+1}/{len(NS)}", flush=True)
        np.savez("core_sum_1e9.npz", rows=np.array(rows))
d = np.array(rows)
r = d[:, 2]
print(f"[T(N) @ 1e9, N {len(d)}개]")
print(f"r 평균 {r.mean():.3f}  중앙값 {np.median(r):.3f}  최악 {r.max():.3f}")
print(f"참조 10⁸: 평균 0.751 / 중앙값 0.599 / 최악 3.14; 반정규 0.798")
print(f"부호 + {int((d[:,1]>0).sum())} / − {int((d[:,1]<0).sum())}")
np.savez("core_sum_1e9.npz", rows=d)
print("전체완료", flush=True)
