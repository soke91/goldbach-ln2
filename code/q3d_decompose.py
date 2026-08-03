"""Q3' 최종 유도 — 구조 측: 분해 항등식 검증 + Q̃ 지배 확인.

정확한 항등식 (대수적 사실, 검증 대상은 구현):
  R(n) = Σ_{a+b=n} Λ(a)Λ(b),  e(m) = Λ(m) − 1 (m ≥ 1, Λ(1)=0 → e(1)=-1)
  R(n) = (n−1) + 2[ψ(n−1) − (n−1)] + Q(n),   Q(n) = Σ_{a+b=n} e(a)e(b)

분산 분해:
  F = R − n𝔖 = 2(ψ(n−1)−(n−1)) [장파장·ζ영점] + Q̃(n) [단파장·영점쌍] + 결정항
  Var(F) 에서 Q̃의 몫과 교차항을 실측.
"""

import numpy as np

from goldbach.sieve import primes_upto
from goldbach.stats import TWIN_PRIME_C2

N = 2_000_000
lam = np.zeros(N + 1)
for p in primes_upto(N):
    p = int(p)
    pk = p
    while pk <= N:
        lam[pk] = np.log(p)
        pk *= p

e = lam.copy() - 1.0
e[0] = 0.0  # a=0 은 합에 없음

L = 2 * (N + 1)
R = np.fft.irfft(np.fft.rfft(lam, L) ** 2, L)[: N + 1]
Q = np.fft.irfft(np.fft.rfft(e, L) ** 2, L)[: N + 1]
psi = np.cumsum(lam)

# 항등식 검증 (n 표본)
print("[항등식 검증] R(n) = (n−1) + 2(ψ(n−1)−(n−1)) + Q(n)")
errs = []
for n in [1000, 12346, 123456, 1234568, 1999998]:
    lhs = R[n]
    rhs = (n - 1) + 2 * (psi[n - 1] - (n - 1)) + Q[n]
    errs.append(abs(lhs - rhs))
    print(f"  n={n:>9,}: |좌−우| = {abs(lhs-rhs):.2e}")
print(f"  → 최대 오차 {max(errs):.2e} (FFT 반올림 수준이면 항등식 확인)\n")

# 분산 분해 (옥타브, 얇은 띠)
def singular(n):
    c, m, p = 1.0, n, 3
    while m % 2 == 0: m //= 2
    while p * p <= m:
        if m % p == 0:
            c *= (p - 1) / (p - 2)
            while m % p == 0: m //= p
        p += 2
    if m > 1: c *= (m - 1) / (m - 2)
    return 2 * TWIN_PRIME_C2 * c

evens = np.arange(2 ** 20, 2 ** 21 - 2 ** 18, 2)
sing = np.array([singular(int(n)) for n in evens])

F = R[evens] - sing * evens
P2 = 2 * (psi[evens - 1] - (evens - 1))        # 장파장 성분
# Q̃ = Q − E[Q]: 항등식에서 정확히 Q̃ = F − 2Δψ
QT = F - P2

def detrend(v):
    x = np.log(evens.astype(float))
    return v - np.polyval(np.polyfit(x, v, 3), x)

Fd, Pd, Qd = detrend(F), detrend(P2), detrend(QT)
vF, vP, vQ = Fd.var(), Pd.var(), Qd.var()
cov = float(np.mean(Pd * Qd))
print("[분산 분해] Var(F) = Var(2Δψ) + Var(Q̃) + 2Cov")
print(f"  Var(F)  = {vF:.4e}")
print(f"  Var(2Δψ)= {vP:.4e} ({vP/vF:.1%})  ← ζ-영점 단일합 (장파장)")
print(f"  Var(Q̃) = {vQ:.4e} ({vQ/vF:.1%})  ← 영점 쌍 간섭 (단파장)")
print(f"  2Cov    = {2*cov:.4e} ({2*cov/vF:+.1%})")
print("""
해석: Var(Q̃)가 지배하면 '억제 상수는 영점 쌍 간섭항 Q̃의 성질'로 확정 —
유도의 무대가 정확히 Q̃ = Σ_{γ,γ'} B(ρ,ρ')·n^{i(γ+γ')} 임을 실증.""")
