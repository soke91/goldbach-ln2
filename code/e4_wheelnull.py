"""돌파 2단계 — 휠-30 HL-널 세계와의 무부기 Fano 비교.

beyond-HL 상수 := VarF(실제) / VarF(휠-30 널), 동일 사영 추정량.
보정 괄호: 실제에만 있는 q≥7 쌍-구조 S₇₊ = Σ_d 2(κ₇₊−1)M_d 를
유도 공식으로 계산해 [비보정, 보정] 범위 제시.
대조 상수: (ln2)² = 0.4805, ln 2 = 0.6931.
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
XtX_inv = np.linalg.inv(X.T @ X)

def varF(ind_full, ind_small):
    R = np.fft.irfft(np.fft.rfft(ind_small, L) * np.fft.rfft(ind_full, L),
                     L)[: N + 1][band]
    F = R - X @ (XtX_inv @ (X.T @ R))
    return float(F @ F)

vr = varF(pi, alpha)
print(f"VarF(실제) = {vr:.4e}", flush=True)

# ── 휠-30 널 세계 ────────────────────────────────────────────
rng = np.random.default_rng(3030)
m_all = np.arange(N + 1)
allowed = (m_all % 2 == 1) & (m_all % 3 != 0) & (m_all % 5 != 0)
allowed[:7] = False
dens = np.zeros(N + 1)
dens[allowed] = np.minimum(3.75 / np.log(np.maximum(m_all[allowed], 3)), 1.0)
# 전체 소수 수를 실제와 정합시키는 상수 보정
dens *= pi.sum() / dens.sum()
vals = []
for _ in range(5):
    draw = (rng.random(N + 1) < dens).astype(float)
    d_small = draw.copy()
    d_small[2 ** 20 :] = 0.0
    vals.append(varF(draw, d_small))
vn = float(np.mean(vals))
print(f"VarF(휠-30 널) = {vn:.4e} ± {np.std(vals):.2e}", flush=True)

# ── q≥7 구조 보정 S₇₊ (유도 공식 + M_all) ────────────────────
Pf = np.fft.rfft(pi, L)
ind_B = np.zeros(N + 1); ind_B[band] = 1.0
mu_B = np.fft.irfft(np.fft.rfft(ind_B, L) * np.conj(Pf), L)[: N + 1]
y = alpha * (mu_B / B)
Yf = np.fft.rfft(y, L)
M_all = np.fft.irfft(np.conj(Yf) * Yf, L)[: N + 1] * B

C7 = 1.0
for q in primes_upto(100000):
    q = int(q)
    if q >= 7:
        C7 *= 1 - 1 / (q - 1) ** 2
spf = np.zeros(N + 1, dtype=np.int32)
for q in range(2, int(N ** 0.5) + 1):
    if spf[q] == 0:
        spf[q*q::q][spf[q*q::q] == 0] = q
def kappa7(d):
    v = C7
    dd = d
    for q0 in (2, 3, 5):
        while dd % q0 == 0:
            dd //= q0
    while dd > 1:
        q = spf[dd] if spf[dd] else dd
        v *= (q - 1) / (q - 2)
        while dd % q == 0:
            dd //= q
    return v
ds_all = np.arange(6, N, 6)
S7 = 0.0
for d in ds_all:
    S7 += 2 * (kappa7(int(d)) - 1.0) * M_all[d]
print(f"C₇ = {C7:.5f} | S₇₊(q≥7 쌍-구조 총량, 비사영) = {S7:+.4e}", flush=True)

r_unc = vr / vn
r_cor = (vr - S7) / vn
print(f"\nbeyond-HL 상수 (무부기 정의):")
print(f"  비보정: VarF_real/VarF_null           = {r_unc:.4f}")
print(f"  보정  : (VarF_real − S₇₊)/VarF_null   = {r_cor:.4f}")
print(f"  대조: (ln2)² = 0.4805 | ln2 = 0.6931 | ln2/2 = 0.3466")
