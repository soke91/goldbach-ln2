"""Q3 결정전 — 스케일 불변 추정량으로 4개 세계 대칭 비교.

추정량 (모든 세계 동일):
  1) 형태 모델 ĝ(n)로 z = g/ĝ 계산
  2) z를 표본 평균으로 재조정 (스케일 불변화)  z ← z/mean(z)
  3) 옥타브 내 log n 1차 드리프트 제거
  4) 비율 = σ(z_resid)·√mean(g)   [포아송 = √(1-μ̄) ≈ 0.96]

세계:
  W1 실제 소수      (형태: HL2 적분 × 특이급수)
  W2 순수 크라메르   (형태: 적분)
  W3 바퀴 mod 30    (형태: 적분 × 쌍허용비)
  W4 강성 크라메르   (h=5000 개수 고정)
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve
from goldbach.stats import TWIN_PRIME_C2

N = 2_000_000
OCT_LO, OCT_HI = 2 ** 19, 2 ** 20
evens = np.arange(OCT_LO, OCT_HI, 2)
rng = np.random.default_rng(123)

samp = np.linspace(OCT_LO, OCT_HI, 300)
base_int = np.array([0.5 * np.trapezoid(
    1 / (np.log(ts) * np.log(n - ts)), ts)
    for n in samp for ts in [np.linspace(3.0, n - 3.0, 2000)]])
Eg_base = np.interp(evens, samp, base_int)


def ratio(g, shape):
    ok = shape > 0
    g2, s2 = g[ok].astype(float), shape[ok]
    z = g2 / s2
    z /= z.mean()
    x = np.log(evens[ok])
    z -= np.polyval(np.polyfit(x, z, 1), x)  # 드리프트 제거
    return float(z.std() * np.sqrt(g2.mean()))


# ── W1 실제 소수 ─────────────────────────────────────────────
is_p = sieve(N)
ind_real = np.zeros(N + 1)
ind_real[np.flatnonzero(is_p)] = 1.0
S = np.fft.rfft(ind_real, 2 * (N + 1))
g_real = np.fft.irfft(S * S, 2 * (N + 1))[: N + 1][evens] / 2
def singular(n):
    c, m, p = 1.0, n, 3
    while m % 2 == 0: m //= 2
    while p * p <= m:
        if m % p == 0:
            c *= (p - 1) / (p - 2)
            while m % p == 0: m //= p
        p += 2
    if m > 1: c *= (m - 1) / (m - 2)
    return c
sing = np.array([2 * TWIN_PRIME_C2 * singular(int(n)) for n in evens])
r1 = ratio(g_real, Eg_base * sing)

# ── W2/W3/W4 랜덤 세계들 ─────────────────────────────────────
def gconv(ind):
    Sf = np.fft.rfft(ind, 2 * (N + 1))
    return np.fft.irfft(Sf * Sf, 2 * (N + 1))[: N + 1][evens] / 2

ms_odd = np.arange(3, N + 1, 2)
p_odd = 2 / np.log(ms_odd)          # 홀수 중 밀도 2/ln = 실제 소수와 동일
p_odd = np.minimum(p_odd, 1.0)

WHEEL = 30
cop = {r for r in range(WHEEL) if np.gcd(r, WHEEL) == 1}
ms_w = ms_odd[np.isin(ms_odd % WHEEL, list(cop))]
p_w = np.minimum((WHEEL / len(cop)) / np.log(ms_w), 1.0)
pf = np.array([sum(1 for a in cop if (int(n) - a) % WHEEL in cop) / len(cop)
               for n in evens])

def world_pure():
    ind = np.zeros(N + 1)
    ind[ms_odd[rng.random(len(ms_odd)) < p_odd]] = 1.0
    return ind

def world_wheel():
    ind = np.zeros(N + 1)
    ind[ms_w[rng.random(len(ms_w)) < p_w]] = 1.0
    return ind

def rigidify(ind, h):
    for lo in range(3, N + 1, h):
        blk = np.arange(lo | 1, min(lo + h, N + 1), 2)
        cur = int(ind[blk].sum())
        target = round(float(np.minimum(2 / np.log(blk), 1).sum()))
        d = target - cur
        if d > 0:
            e = blk[ind[blk] == 0]
            ind[rng.choice(e, size=min(d, len(e)), replace=False)] = 1.0
        elif d < 0:
            f = blk[ind[blk] == 1]
            ind[rng.choice(f, size=min(-d, len(f)), replace=False)] = 0.0
    return ind

r2 = np.mean([ratio(gconv(world_pure()), Eg_base * 2) for _ in range(4)])
r3 = np.mean([ratio(gconv(world_wheel()), Eg_base * (WHEEL / len(cop)) * pf)
              for _ in range(4)])
r4 = np.mean([ratio(gconv(rigidify(world_pure(), 5000)), Eg_base * 2)
              for _ in range(4)])

print(f"{'세계':<28} {'비율':>7}   (포아송 이론치 ≈ 0.95~0.97)")
print(f"{'W2 순수 크라메르 (구조 0)':<26} {r2:>7.3f}")
print(f"{'W3 바퀴 mod30 (잉여류만)':<26} {r3:>7.3f}")
print(f"{'W4 강성 h=5000 (개수 고정)':<25} {r4:>7.3f}")
print(f"{'W1 실제 소수 (모든 구조)':<26} {r1:>7.3f}")
print("""
독해: 네 값의 상대 위치가 전부를 말한다.
  실제≈순수≈바퀴 → g-요동은 그냥 포아송 (H3' 상수는 추정량 인공물이었음)
  실제<순수, 강성<순수 → 강성 기원설(Q3) 부활
  실제>순수, 바퀴>순수 → 산술 상관 기원(H3'')""")
