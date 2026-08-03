"""z-스캔 v3 (상호작용 기저) — Fano 정규화 (세계별 T̃1 사영보정 포함 정확 계산).

Fano(세계) = VarF / T̃1,  T̃1 = Σ_p [B·μ_p − c_p'(X'X)c_p]  (사영 일관)
판정량: Fano_real / Fano_null(z) 의 z-플래토 → (ln2)² 대조.
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
uqs = {}
for q in (5, 7, 11, 13):
    uqs[q] = (band % q == 0).astype(float)
    cols += [uqs[q], uqs[q] * x, t * uqs[q]]
ql = [5, 7, 11, 13]
for i in range(len(ql)):
    for j in range(i + 1, len(ql)):
        cols.append(uqs[ql[i]] * uqs[ql[j]])
X = np.column_stack(cols)
K = X.shape[1]
XtX = X.T @ X
XtX_inv = np.linalg.inv(XtX)
Xb = np.zeros((K, N + 1))
for j in range(K):
    col = np.zeros(N + 1); col[band] = X[:, j]; Xb[j] = col
Xb_f = [np.fft.rfft(Xb[j], L) for j in range(K)]

def fano(full, small):
    R = np.fft.irfft(np.fft.rfft(small, L) * np.fft.rfft(full, L), L)[: N + 1][band]
    F = R - X @ (XtX_inv @ (X.T @ R))
    varF = float(F @ F)
    Pf = np.conj(np.fft.rfft(full, L))
    Xu = np.empty((K, N + 1))
    for j in range(K):
        Xu[j] = np.fft.irfft(Xb_f[j] * Pf, L)[: N + 1]
    c = XtX_inv @ Xu
    idx = np.flatnonzero(small > 0)
    quad = np.einsum("jp,jk,kp->p", c[:, idx], XtX, c[:, idx])
    T1t = float(np.sum(Xu[0][idx]) - quad.sum())
    return varF, T1t

vF_r, T1_r = fano(pi, alpha)
fano_r = vF_r / T1_r
print(f"실제: VarF {vF_r:.3e} T̃1 {T1_r:.3e} Fano {fano_r:.4f}", flush=True)

rng = np.random.default_rng(4242)
m_all = np.arange(N + 1)

print(f"\n{'z':>4} {'Fano(널)':>9} {'Fano비':>8} {'±':>7}   (ln2)²={np.log(2)**2:.4f}")
for z in [3, 5, 7, 11, 13]:
    qs = [q for q in [3, 5, 7, 11, 13] if q <= z]
    allowed = m_all % 2 == 1
    boost = 2.0
    for q in qs:
        allowed &= m_all % q != 0
        boost *= q / (q - 1)
    allowed[:15] = False
    dens = np.zeros(N + 1)
    dens[allowed] = np.minimum(boost / np.log(np.maximum(m_all[allowed], 3)), 1.0)
    dens *= pi.sum() / dens.sum()
    fs = []
    for _ in range(8):
        draw = (rng.random(N + 1) < dens).astype(float)
        ds = draw.copy(); ds[2 ** 20 :] = 0.0
        v, t1 = fano(draw, ds)
        fs.append(v / t1)
    fs = np.array(fs)
    ratio = fano_r / fs.mean()
    se = ratio * fs.std() / np.sqrt(len(fs)) / fs.mean()
    print(f"{z:>4} {fs.mean():>9.4f} {ratio:>8.4f} {se:>7.4f}", flush=True)
