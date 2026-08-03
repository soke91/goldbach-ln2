"""A안 실현 — 잔차장 직접 스펙트럼 (파세발 정확 항등식).

F = 사영 잔차(밴드 위), 0-임베딩 후 FFT:
  ||F||² = (1/L)Σ_k |F̂(k)|²  (정확)
산출:
  측도 m(k) = |F̂_널(k)|² 분포 (λ-옥타브 빈)
  적분항 r(k) = |F̂_실|²/|F̂_널|² 프로파일
  Σ m·r = Fano (항등식 — 자동 성립; 관심은 m과 r의 형태)
  GM 대조: r(λ) ≈ [1 − lnλ/lnN]²-형인지, m의 집중 위치는 어디인지.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve
from goldbach.stats import TWIN_PRIME_C2

N = 2 ** 21
L = 2 * (N + 1)
is_p = sieve(N)
pi = np.zeros(N + 1)
pi[np.flatnonzero(is_p)] = 1.0
alpha = pi.copy()
alpha[2 ** 20 :] = 0.0

band = np.arange(2 ** 20 + 2 ** 18, 2 ** 21 - 2 ** 18, 2)
band = band[band % 6 == 2]
B = len(band)

ps2000 = [int(p) for p in primes_upto(2000) if p > 2]
fac = {p: (p - 1) / (p - 2) for p in ps2000}
T0 = np.full(B, 2 * TWIN_PRIME_C2)
for p in ps2000:
    T0[band % p == 0] *= fac[p]
x = np.log(band.astype(float)); x = (x - x.mean()) / x.std()
t = T0 - 1.0
cols = [np.ones(B), x, x**2, x**3, t, t*x, t*x**2, t*x**3, t*t, t*t*x]
for q in (5, 7, 11, 13):
    uq = (band % q == 0).astype(float)
    cols += [uq, uq * x]
X = np.column_stack(cols)
pinv = np.linalg.pinv(X)

def resid_field(full, small):
    R = np.fft.irfft(np.fft.rfft(small, L) * np.fft.rfft(full, L), L)[: N + 1][band]
    F = R - X @ (pinv @ R)
    emb = np.zeros(N + 1)
    emb[band] = F
    return emb, float(F @ F)

emb_r, v_r = resid_field(pi, alpha)
Fr = np.abs(np.fft.rfft(emb_r)) ** 2

rng = np.random.default_rng(2)
ms = np.arange(3, N + 1, 2)
pv = np.minimum(2 / np.log(ms), 1.0)
Fn = np.zeros(len(Fr))
v_n = 0.0
NW = 5
for _ in range(NW):
    b_ = np.zeros(N + 1)
    pick = ms[rng.random(len(ms)) < pv]
    b_[pick] = 1.0
    bs = b_.copy(); bs[2 ** 20 :] = 0.0
    emb, v = resid_field(b_, bs)
    Fn += np.abs(np.fft.rfft(emb)) ** 2 / NW
    v_n += v / NW

print(f"직접 Fano비 = {v_r/v_n:.4f}")
kmax = len(Fr)
kf = np.arange(kmax)
lam = np.maximum((N + 1) / np.maximum(kf, 1), 2.0)
lnN = np.log(N)
print(f"\n{'λ-옥타브':>18} {'측도 m':>8} {'적분항 r':>9} {'GM²모델':>8}")
tot_n = Fn.sum()
acc = 0.0
for j in range(1, 21):
    lo, hi = 2.0 ** j, 2.0 ** (j + 1)
    m_ = (lam >= lo) & (lam < hi)
    if not m_.any():
        continue
    w = float(Fn[m_].sum()) / tot_n
    if w < 1e-5:
        continue
    r = float(Fr[m_].sum()) / float(Fn[m_].sum())
    gm = max(1 - np.log(np.sqrt(lo * hi)) / lnN, 0.03) ** 2
    acc += w * r
    print(f"[{lo:>8.0f},{hi:>8.0f}) {w:>8.4f} {r:>9.4f} {gm:>8.4f}")
print(f"\nΣ m·r = {acc:.4f} (= Fano비 {v_r/v_n:.4f} — 파세발 검증)")
print(f"(ln2)² = {np.log(2)**2:.4f}")
