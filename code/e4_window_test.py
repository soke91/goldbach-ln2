"""결정적 반증 시험 — ln 2 상수의 창(window) 의존성 검사.

유도 스케치가 제기한 위험: 옥타브(비율 2) 창의 ln 2 폭이 상수의 기원일
가능성. 검사: 같은 중심 스케일에서 창 비율 r을 바꿔 순수요동 초과비 측정.
  r-불변 → ln 2는 소수 고유 상수 (추측 생존, 유도는 창-독립 경로로)
  r-의존(특히 값 ≈ ln r 계열) → 아티팩트 (추측 재구성 필요!)
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve
from goldbach.stats import TWIN_PRIME_C2

N = 12_000_000
is_p = sieve(N)
ind = np.zeros(N + 1)
ind[np.flatnonzero(is_p)] = 1.0
S = np.fft.rfft(ind, 2 * (N + 1))
conv = np.fft.irfft(S * S, 2 * (N + 1))[: N + 1]
evens = np.arange(4, N + 1, 2)
g = conv[evens] / 2

sing = np.full(len(evens), TWIN_PRIME_C2)  # HL1 = C2·∏·n/ln² (비순서)
for p in primes_upto(10_000):
    p = int(p)
    if p > 2:
        sing[(evens % p) == 0] *= (p - 1) / (p - 2)
hl = sing * evens / np.log(evens) ** 2
z = g / hl

CENTER = 4_000_000
print(f"중심 스케일 ≈ {CENTER:,} | 창 비율 r 스캔")
print(f"{'r':>5} {'창':>22} {'표본':>7} {'순수 초과비':>10}")
for r in [1.3, 1.5, 2.0, 3.0, 4.0, 6.0]:
    lo = int(CENTER / np.sqrt(r))
    hi = int(CENTER * np.sqrt(r))
    m = (evens >= lo) & (evens < hi) & (evens % 6 == 2)
    nn, zz, hh = evens[m].astype(float), z[m], hl[m]
    x = np.log(nn)
    resid = zz - np.polyval(np.polyfit(x, zz, 1), x)
    ratio = float(resid.std() / np.mean(1 / np.sqrt(hh)))
    print(f"{r:>5.1f} [{lo:>9,},{hi:>9,}) {int(m.sum()):>7,} {ratio:>10.4f}")
print(f"\n참고: 시리즈의 이 스케일 예상값 ≈ 0.58~0.60 (r=2 기준)")
print(f"ln r 계열이 나타나면: ln1.3={np.log(1.3):.3f} ln1.5={np.log(1.5):.3f} "
      f"ln2={np.log(2):.3f} ln3={np.log(3):.3f} ln4={np.log(4):.3f} ln6={np.log(6):.3f}")
