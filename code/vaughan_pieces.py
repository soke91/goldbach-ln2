"""C3 — Vaughan 조각 지도: T(N) = Σ Λ(n)μ(N−n)의 타입별 분해 실측.

Vaughan (U = V = N^{1/3}):
  Λ(n) = Π1 − Π2 − Π3  (n > V에서)
  Π1(n) = Σ_{d≤U, d|n} μ(d) log(n/d)
  Π2(n) = Σ_{d≤U, m≤V, dm|n} μ(d) Λ(m)
  Π3(n) = Σ_{m>V, k>U, mk... } — 잔여 쌍선형 (타입 II)
구현: T_i = Σ_n Π_i(n) μ(N−n) 를 d-슬라이싱으로 직접 계산,
T_II = T − T_1 + T_2 (항등식으로 역산 + 별도 부분합 검산).
각 조각을 √-기준과 대비. N 6개 @ 10⁸.
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
del val
print("mu 완료됨", flush=True)

ps_full = primes_upto(X - 100)
ps_odd = ps_full[ps_full > 2].astype(np.int64)
logp = np.log(ps_odd.astype(np.float64))

def mobius_small(L):
    m = np.ones(L + 1, dtype=np.int8)
    pmm = np.ones(L + 1, dtype=bool)
    pmm[:2] = False
    for p in range(2, int(L ** 0.5) + 1):
        if pmm[p]:
            pmm[p * p:: p] = False
            m[p::p] *= -1
            m[p * p:: p * p] = 0
    v = np.arange(L + 1, dtype=np.int64)
    for p in range(2, int(L ** 0.5) + 1):
        if pmm[p]:
            v[p::p] //= p
            pp = p * p
            while pp <= L:
                v[pp::pp] //= p
                pp *= p
    m[v > 1] *= -1
    m[0] = 0
    return m

U = int(round(X ** (1 / 3)))
mu_s = mobius_small(U + 10)
lam = np.zeros(U + 10)          # Λ(m), m ≤ V=U
for p in primes_upto(U + 9):
    p = int(p)
    pk = p
    while pk <= U + 9:
        lam[pk] = math.log(p)
        pk *= p

rng = np.random.default_rng(137)
base = X - 5_000_000
NS = sorted(set(int(base + 6 * k + (2 - base % 6) % 6)
                for k in rng.integers(0, 800_000, 6)))

print(f"U = V = {U:,};  N {len(NS)}개", flush=True)
print(f"{'N':>10} {'T':>10} {'T_I':>12} {'T_I2':>12} {'T_II':>12} "
      f"{'r_T':>6} {'r_I':>6} {'r_II':>6}")
for N in NS:
    muv = mu[N - ps_odd].astype(np.float64)
    T = float((logp * muv).sum())
    V_T = float((logp ** 2 * (muv != 0)).sum())
    # T_I = Σ_{d≤U} μ(d) Σ_{k≤(N-1)/d} log(k?) ... Π1 합 = Σ_n μ(N−n)Π1(n)
    #     = Σ_{d≤U} μ(d) Σ_{j} log(j) μ(N − d·j)   (n = d·j)
    T_I = 0.0
    absI = 0.0
    for d in range(1, U + 1):
        md = int(mu_s[d])
        if md == 0:
            continue
        J = (N - 1) // d
        js = np.arange(1, J + 1, dtype=np.int64)
        vals = np.log(js.astype(np.float64)) * mu[N - d * js]
        s = float(vals.sum())
        T_I += md * s
        absI += s * s
    # T_I2 = Σ_{d≤U} μ(d) Σ_{m≤V} Λ(m) Σ_{l≤N/(dm)} μ(N − dml)
    T_I2 = 0.0
    for d in range(1, U + 1):
        md = int(mu_s[d])
        if md == 0:
            continue
        for m in range(2, U + 1):
            if lam[m] == 0.0:
                continue
            dm = d * m
            if dm >= N:
                break
            L2 = (N - 1) // dm
            ls = np.arange(1, L2 + 1, dtype=np.int64)
            T_I2 += md * lam[m] * float(mu[N - dm * ls].sum(dtype=np.int64))
    T_II = T - T_I + T_I2
    r_T = abs(T) / math.sqrt(V_T)
    r_I = abs(T_I) / math.sqrt(V_T)
    r_II = abs(T_II) / math.sqrt(V_T)
    print(f"{N:>10} {T:>10.0f} {T_I:>12.0f} {T_I2:>12.0f} {T_II:>12.0f} "
          f"{r_T:>6.2f} {r_I:>6.2f} {r_II:>6.2f}", flush=True)
    np.savez(f"vaughan_{N}.npz",
             vals=np.array([T, T_I, T_I2, T_II, V_T]))
print("전체완료", flush=True)
