"""돌파 계산 — 쌍-조건부 HL 구조의 정확 소거와 beyond-HL 프로파일.

유도된 구조 모델 (밴드 n≡2(6) ⇒ 허용 p ≡ 1(6) ⇒ 6|d):
  κ(d) = C₅·∏_{q≥5, q|d}(q−1)/(q−2),  C₅ = ∏_{q≥5}(1−1/(q−1)²) = 0.88021
  C_struct(d) = (κ(d)−1)·M_d,  M_d = Σ_p α(p)μ_p·α(p+d)μ_{p+d}  [전 d FFT 일괄]
검증·산출:
  A. 측정 d에서 C_beyond(d) = C̃(d) − C_struct(d) — 원거리에서 ≈0이면 모델 완성
  B. 전 d 구조 합 Σ C_struct (정확) → beyond 총량 = Var(F) − T̃1 − ΣC_struct
  C. beyond의 d-프로파일 → ln 2 기전 판독
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
X = np.column_stack([np.ones(B), x, x**2, x**3, t, t*x, t*x**2, t*x**3,
                     t*t, t*t*x])
K = X.shape[1]
XtX = X.T @ X
XtX_inv = np.linalg.inv(XtX)

R = np.fft.irfft(np.fft.rfft(alpha, L) * np.fft.rfft(pi, L), L)[: N + 1][band]
F = R - X @ (XtX_inv @ (X.T @ R))
varF = float(F @ F)

Xb = np.zeros((K, N + 1))
for j in range(K):
    col = np.zeros(N + 1); col[band] = X[:, j]; Xb[j] = col
Pf = np.fft.rfft(pi, L)
Xu = np.empty((K, N + 1))
for j in range(K):
    Xu[j] = np.fft.irfft(np.fft.rfft(Xb[j], L) * np.conj(Pf), L)[: N + 1]
c = XtX_inv @ Xu
mu_B = Xu[0]
p_idx = np.flatnonzero(alpha > 0)
quad_diag = np.einsum("jp,jk,kp->p", c[:, p_idx], XtX, c[:, p_idx])
T1t = float(np.sum(mu_B[p_idx]) - quad_diag.sum())

# ── M_d 전 d 일괄 (자기상관 FFT) ─────────────────────────────
y = alpha * (mu_B / B)
Yf = np.fft.rfft(y, L)
M_all = np.fft.irfft(np.conj(Yf) * Yf, L)[: N + 1] * B  # Σ_p y(p)y(p+d)·B → ⟨⟩스케일 맞춤
# 주의 스케일: C̃(d)는 내적 스케일(⟨u,u⟩−보정), M_d 내적 스케일 = Σ_p (Bμ)(Bμ)/B = B·Σμμ
# → M_inner(d) = B·Σ_p μ_pμ_{p+d} = M_all (위 정의와 일치)

# ── κ(d) (spf 체) ────────────────────────────────────────────
C5 = 1.0
for q in primes_upto(100000):
    q = int(q)
    if q >= 5:
        C5 *= 1 - 1 / (q - 1) ** 2
print(f"C₅ = {C5:.5f}")
Dmax = N
spf = np.zeros(Dmax + 1, dtype=np.int32)
for q in range(2, int(Dmax ** 0.5) + 1):
    if spf[q] == 0:
        spf[q*q::q][spf[q*q::q] == 0] = q
# κ 계산 함수
def kappa(d):
    v = C5
    dd = d
    while dd % 2 == 0: dd //= 2
    while dd % 3 == 0: dd //= 3
    while dd > 1:
        q = spf[dd] if spf[dd] else dd
        v *= (q - 1) / (q - 2)
        while dd % q == 0:
            dd //= q
    return v

# ── B. 전 d 구조 합 ──────────────────────────────────────────
ds_all = np.arange(6, Dmax, 6)
kap = np.empty(len(ds_all))
for i, d in enumerate(ds_all):
    kap[i] = kappa(int(d))
Cstruct_all = (kap - 1.0) * M_all[ds_all]
S_struct = float(2 * Cstruct_all.sum())
print(f"B·Var(F) = {varF:.3e} | T̃1 = {T1t:.3e}")
print(f"Σ C_struct (전 d, ±) = {S_struct:+.3e}")
S_beyond = varF - T1t - S_struct
print(f"⇒ beyond-HL 총량 = Var − T̃1 − ΣC_struct = {S_beyond:+.3e}")
print(f"   비율로: T̃1 기준 beyond = {S_beyond/T1t:+.3f}, "
      f"Fano 분해: 1 {S_struct/T1t:+.3f} {S_beyond/T1t:+.3f} = {varF/T1t:.3f}")

# ── A/C. 측정 d에서 beyond 프로파일 ─────────────────────────
E = np.empty((K, L // 2 + 1), complex)
for j in range(K):
    E[j] = np.fft.rfft(alpha * c[j], L)
corrfull = np.zeros(N + 1)
for j in range(K):
    for k in range(K):
        corrfull += XtX[j, k] * np.fft.irfft(np.conj(E[j]) * E[k], L)[: N + 1]

ds = ([6, 12, 18, 24, 30, 36, 42, 60, 90, 120, 210, 330, 630, 1050, 2310,
      4620, 9240, 18480, 36960, 73920, 147840, 295680])
print(f"\n{'d':>7} {'C̃(d)':>11} {'C_struct':>11} {'beyond':>10}")
cum = 0.0
for d in ds:
    v = np.zeros(N + 1); v[d:] = pi[d:] * pi[: N + 1 - d]
    w = np.zeros(N + 1); w[: N + 1 - d] = alpha[: N + 1 - d] * alpha[d:]
    A = np.fft.irfft(np.fft.rfft(w, L) * np.fft.rfft(v, L), L)[: N + 1][band]
    Cd = float(A.sum()) - float(corrfull[d])
    Cs = float((kappa(d) - 1.0) * M_all[d])
    cum += 2 * (Cd - Cs)
    print(f"{d:>7} {Cd:>11.1f} {Cs:>11.1f} {Cd-Cs:>10.1f}  누적beyond {cum:+.2e}")
