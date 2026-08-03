"""최종 정식화의 심장 — φ(k) 프로파일 측정.

φ(k) = |û(k)|²/샷,  û = 창-소수 지시함수 − 평활 밀도 (마스킹 없음).
k-옥타브 빈별로 실제/베르누이 비 → 프로파일 φ̂(k).
판정: 프로파일이 평탄(상수)이면 f 잘 정의 → 그 값 = ln 2 인지.
      평탄 아니면 ⟨φ_α φ_π⟩_G 적분이 (ln2)²인지 (G-가중 결합).
N-스캔(2^22, 2^24)으로 프로파일의 스케일 안정성 동시 확인.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

rng = np.random.default_rng(48)

def phi_profile(lo, hi, is_p, n_bern=3):
    span = hi - lo
    L = 2 * span
    m = np.arange(lo, hi)
    dens = np.zeros(span)
    odd = m % 2 == 1
    dens[odd] = 2 / np.log(m[odd])
    shot = float((dens * (1 - dens)).sum())

    def spec(ind):
        u = ind - dens * (ind.sum() / dens.sum())
        return np.abs(np.fft.rfft(u, L)) ** 2

    P_real = spec(is_p[lo:hi].astype(float))
    P_bern = np.mean([spec((rng.random(span) < dens).astype(float))
                      for _ in range(n_bern)], axis=0)
    kmax = len(P_real)
    out = []
    k0 = 2048
    while k0 < kmax:
        k1 = min(2 * k0, kmax)
        r = float(P_real[k0:k1].mean()) / float(P_bern[k0:k1].mean())
        out.append((k0, k1, r, span / k0))
        k0 = k1
    return out

for e in [22, 24]:
    N = 2 ** e
    is_p = sieve(N)
    print(f"\n═══ N = 2^{e} ═══")
    for lo, hi, name in [(3, N // 2, "하반"), (N // 2, N, "상반")]:
        prof = phi_profile(lo, hi, is_p)
        print(f"[{name}] {'k-빈':>16} {'파장(사이트)':>12} {'φ(실/베)':>9}")
        for k0, k1, r, lam in prof:
            print(f"     [{k0:>7},{k1:>7}) {lam:>12.0f} {r:>9.4f}")
print(f"\nln2 = {np.log(2):.4f}, 1 = 평탄 랜덤")
