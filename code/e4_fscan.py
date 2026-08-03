"""f 상수(=minor-arc 전력비)의 N-스캔 — ln 2 수렴 판정."""

import numpy as np

from goldbach.sieve import primes_upto, sieve

rng = np.random.default_rng(693)

def frac(lo, hi, is_p, world):
    L = 2 * (hi - lo)
    m = np.arange(lo, hi)
    dens = np.zeros(hi - lo)
    odd = m % 2 == 1
    dens[odd] = 2 / np.log(m[odd])
    if world == "real":
        ind = is_p[lo:hi].astype(float)
    else:
        ind = (rng.random(hi - lo) < dens).astype(float)
    u = ind - dens * (ind.sum() / dens.sum())
    P = np.abs(np.fft.rfft(u, L)) ** 2
    kmax = len(P)
    maj = np.zeros(kmax, bool)
    for q in range(1, 61):
        stepf = L / q
        a = 0
        while a * stepf < kmax:
            c = int(round(a * stepf))
            maj[max(0, c - 64) : min(kmax, c + 65)] = True
            a += 1
    maj[:2048] = True
    shot = float((dens * (1 - dens)).sum())
    return float(P[~maj].mean()) / shot

print(f"{'N':>7} {'f하반':>8} {'f상반':>8}   ln2 = {np.log(2):.4f}")
for e in [21, 22, 23, 24]:
    N = 2 ** e
    is_p = sieve(N)
    fl = frac(3, N // 2, is_p, "real") / frac(3, N // 2, is_p, "bern")
    fu = frac(N // 2, N, is_p, "real") / frac(N // 2, N, is_p, "bern")
    print(f"2^{e:>2}   {fl:>8.4f} {fu:>8.4f}", flush=True)
