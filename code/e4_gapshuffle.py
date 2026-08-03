"""간격 순열 실험 — 앨리어싱 없는 초근거리 프로브.

세계:
  G1 전역 간격 순열: 실제 소수 간격 멀티셋 보존, 순서 무작위
     (간격 분포 = 완벽 실제, 간격상관·잔차정렬 = 파괴)
  G2 국소 간격 순열 (K=8 간격 창 내 순열): 중거리 구조 보존, 초근거리만 파괴
  비교: 실제 0.257 / 베르누이 1.0
판정:
  G1 ≈ 0.26 → 운반자 = 간격 분포 (1점 통계!) — 대발견
  G1 ≈ 1, G2 ≈ 0.26 → 운반자 = 초근거리 간격 배열(상관)
  둘 다 ≈ 1 → 운반자 = 잔차 정렬 포함 정확 구조
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve
from goldbach.stats import TWIN_PRIME_C2

N = 2 ** 21
L = 2 * (N + 1)
lam = np.zeros(N + 1)
for p in primes_upto(N):
    p = int(p)
    pk = p
    while pk <= N:
        lam[pk] = np.log(p)
        pk *= p
e = lam - 1.0
e[:2] = 0.0
start = int(0.55 * N)
start += (2 - start % 6) % 6
band = np.arange(start, int(0.95 * N), 6)
M = len(band)

ps2000 = [int(p) for p in primes_upto(2000) if p > 2]
band_mod = {p: (band % p).astype(np.int16) for p in ps2000}
fac = {p: (p - 1) / (p - 2) for p in ps2000}
def sing_shift(s):
    v = np.full(M, 2 * TWIN_PRIME_C2)
    for p in ps2000:
        v[band_mod[p] == (s % p)] *= fac[p]
    return v
x = np.log(band.astype(float)); x = (x - x.mean()) / x.std()
cols = [np.ones(M), x, x**2, x**3]
nscale = band / band.mean()
for s in range(-16, 17, 2):
    cols.append((sing_shift(s) - 1.0) * nscale)
X = np.column_stack(cols)
pinv = np.linalg.pinv(X)

def resid_var(sig):
    E = np.fft.rfft(sig, L)
    Q = np.fft.irfft(E * E, L)[: N + 1][band]
    return float((Q - X @ (pinv @ Q)).var())

rng = np.random.default_rng(3232)
primes = np.flatnonzero(sieve(N))
primes = primes[primes >= 3]
gaps = np.diff(primes)

def world_from_gaps(g):
    pos = 3 + np.concatenate([[0], np.cumsum(g)])
    pos = pos[pos <= N]
    sig = np.full(N + 1, -1.0)
    sig[np.arange(0, N + 1, 2)] = e[np.arange(0, N + 1, 2)]
    sig[:2] = 0.0
    sig[pos] = np.log(pos) - 1.0
    return sig

v_real = resid_var(e)
ms = np.arange(3, N + 1, 2)
pv = np.minimum(2 / np.log(ms), 1.0)
v_bern = []
for _ in range(3):
    b = np.full(N + 1, -1.0)
    b[np.arange(0, N + 1, 2)] = e[np.arange(0, N + 1, 2)]
    b[:2] = 0.0
    pick = ms[rng.random(len(ms)) < pv]
    b[pick] = np.log(pick) - 1.0
    v_bern.append(resid_var(b))
v_bern = float(np.mean(v_bern))
print(f"실제 {v_real/v_bern:.4f} / 베르누이 1.0")

# G1 전역 간격 순열
v_g1 = []
for _ in range(3):
    g = gaps.copy()
    rng.shuffle(g)
    v_g1.append(resid_var(world_from_gaps(g)))
print(f"G1 전역 간격 순열      : {np.mean(v_g1)/v_bern:.4f} ± {np.std(v_g1)/v_bern:.4f}")

# G2 국소 간격 순열 (8-간격 창)
v_g2 = []
K = 8
for _ in range(3):
    g = gaps.copy().astype(float)
    key = (np.arange(len(g)) // K) * 10.0 + rng.random(len(g))
    g = gaps[np.argsort(key)]
    v_g2.append(resid_var(world_from_gaps(g)))
print(f"G2 국소 간격 순열 (K=8): {np.mean(v_g2)/v_bern:.4f} ± {np.std(v_g2)/v_bern:.4f}")
