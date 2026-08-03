"""① F_S 유도 — 휠 쌍상관 공식으로 Var(S) 예측 vs 실측.

공식: Var(S) = S̄(1−⟨w⟩) + Σ_{d≠0} N(d)·⟨w⟩²·(ρ_w(d)−1)
  ⟨w⟩ = S̄/π_allowed (실측 캘리브레이션)
  N(d) = 허용 소수(p≡1 mod 6, p ≤ X/2) 쌍 수 — FFT 자기상관
  ρ_w(d) = ∏_{5≤q≤y, q|d} [q/(q−1)] / [q(q−2)/(q−1)²]  (일반항은 공통상수 —
           C_y = ∏_{5≤q≤y}(q(q−2)/(q−1)²) 로 묶임: ρ_w(d) = C_y·∏_{q|d,5≤q≤y}(q−1)²/((q−2)q)·(q/(q−1))…
           정리하면 ρ_w(d) = C_y·∏_{q|d}( (q−1)/(q−2) )   [κ-형!]
예측 F_S vs 실측 0.137 (10⁷).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 10_000_000
y = int(round(X ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q >= 5]
is_p = sieve(X)

# 허용 소수: p ≡ 1 (mod 6), p ≤ X/2
alpha = np.flatnonzero(is_p[: X // 2])
alpha = alpha[alpha % 6 == 1]
A = len(alpha)
print(f"허용 소수 {A:,}개, y = {y}")

# N(d): FFT 자기상관
H = X // 2
ind = np.zeros(H)
ind[alpha] = 1.0
L = 2 * H
E = np.fft.rfft(ind, L)
N_d = np.fft.irfft(np.conj(E) * E, L)[:H]

# 상수들
C_y = 1.0
for q in qs:
    C_y *= q * (q - 2) / (q - 1) ** 2
print(f"C_y = ∏(1−1/(q−1)²) [5≤q≤{y}] = {C_y:.5f}")

# 실측 캘리브레이션 (S̄, ⟨w⟩) — 소표본
ps_all = primes_upto(X // 2 + 50_000)
n0s = [int(v + (2 - v % 6) % 6) for v in
       np.linspace(X, X * 1.05, 30).astype(np.int64)]
is_p_hi = sieve(int(X * 1.05) + 10)
Ss = []
for n in n0s:
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= (m % 3 != 0) & (m > 1)
    Ss.append(int(al.sum()))
Sbar = float(np.mean(Ss))
wbar = Sbar / A
print(f"S̄ = {Sbar:,.0f}, ⟨w⟩ = {wbar:.4f}")

# ρ_w(d) − 1 = C_y·∏_{q|d, 5≤q≤y}(q−1)/(q−2) − 1  : spf 인수분해
spf = np.zeros(H, dtype=np.int32)
for q in range(2, int(H ** 0.5) + 1):
    if spf[q] == 0:
        spf[q*q::q][spf[q*q::q] == 0] = q
def rho_w(d):
    v = C_y
    dd = d
    for q0 in (2, 3):
        while dd % q0 == 0:
            dd //= q0
    while dd > 1:
        q = spf[dd] if dd < H and spf[dd] else dd
        if 5 <= q <= y:
            v *= (q - 1) / (q - 2)
        while dd % q == 0:
            dd //= q
    return v

# 합산: d는 6의 배수만 (허용 소수 차)
offdiag = 0.0
ds = np.arange(6, H, 6)
for d in ds:
    if N_d[d] > 0:
        offdiag += 2 * N_d[d] * (rho_w(int(d)) - 1.0)
offdiag *= wbar ** 2

var_pred = Sbar * (1 - wbar) + offdiag
F_pred = var_pred / Sbar
print(f"\n대각항 (1−⟨w⟩)      = {1-wbar:.4f}")
print(f"오프대각 합/S̄       = {offdiag/Sbar:+.4f}")
print(f"**예측 F_S = {F_pred:.4f}**  vs 실측 0.137")
