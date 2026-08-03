"""Q3' 해석적 유도 1부 — 골드바흐 분광학: 데이터에서 리만 영점 검출.

이론 (Fujii 1991):
  G(x) = Σ_{n≤x} R(n),  R(n) = Σ_{a+b=n} Λ(a)Λ(b)
  G(x) = x²/2 − 2 Σ_ρ x^{ρ+1}/(ρ(ρ+1)) + O(x ln³x)
RH 하에서 ρ = ½+iγ → 요동항 h(u) := (G(e^u) − e^{2u}/2)/e^{3u/2}
  = −2 Σ_γ Re[ e^{iγu} / (ρ(ρ+1)) ] + 잡음
→ h(u)의 주파수 스펙트럼은 γ = 14.13, 21.02, 25.01, ... 에서 피크.

실험:
  A. 스펙트럼 피크 위치 vs 실제 영점 (검출 시 '요동=영점' 수치 증명)
  B. 영점 30개로 h(u)를 이론 재구성 → 실측과 상관계수
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

N = 2_000_000
ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
         40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248,
         59.347044, 60.831779, 65.112544, 67.079811, 69.546402, 72.067158,
         75.704691, 77.144840, 79.337375, 82.910381, 84.735493, 87.425275,
         88.809111, 92.491899, 94.651344, 95.870634, 98.831194, 101.317851]

# ── Λ 배열과 R = Λ⋆Λ, 누적합 G ──────────────────────────────
lam = np.zeros(N + 1)
for p in primes_upto(N):
    p = int(p)
    pk = p
    while pk <= N:
        lam[pk] = np.log(p)
        pk *= p
S = np.fft.rfft(lam, 2 * (N + 1))
R = np.fft.irfft(S * S, 2 * (N + 1))[: N + 1]
G = np.cumsum(R)

# ── h(u) 구성 (u = ln x, x ∈ [10^4, 2×10^6]) ────────────────
M = 4096
us = np.linspace(np.log(10_000), np.log(N), M)
xs = np.exp(us)
xi = np.minimum(xs.astype(int), N)
h = (G[xi] - xs ** 2 / 2) / xs ** 1.5
# 부드러운 저차항 제거 (저주파만 죽임 — 영점 대역 γ>10 무손상)
h = h - np.polyval(np.polyfit(us, h, 5), us)

# ── A. 주기도(periodogram) → 피크 vs 영점 ───────────────────
win = np.hanning(M)
H = np.fft.rfft(h * win, 8 * M)
freqs = 2 * np.pi * np.fft.rfftfreq(8 * M, d=us[1] - us[0])
band = (freqs > 8) & (freqs < 105)
P = np.abs(H) ** 2
Pb, fb = P[band], freqs[band]
# 국소 최대 탐색
peaks = []
for i in range(2, len(Pb) - 2):
    if Pb[i] == Pb[i - 2 : i + 3].max() and Pb[i] > Pb.mean() * 3:
        peaks.append((float(Pb[i]), float(fb[i])))
peaks = sorted(peaks, reverse=True)[:12]
print("A. 골드바흐 스펙트럼 피크 vs 리만 영점")
print(f"{'피크 γ(검출)':>11} {'가장 가까운 영점':>12} {'오차':>7}")
hits = 0
for _, f in sorted(peaks, key=lambda t: t[1]):
    nearest = min(ZEROS, key=lambda z: abs(z - f))
    err = f - nearest
    hits += abs(err) < 0.5
    print(f"{f:>11.3f} {nearest:>12.4f} {err:>+7.3f}")
print(f"→ {len(peaks)}개 피크 중 {hits}개가 영점 ±0.5 이내\n")

# ── B. 영점 30개 이론 재구성과의 상관 ───────────────────────
recon = np.zeros(M)
for g0 in ZEROS:
    rho = 0.5 + 1j * g0
    coef = 1.0 / (rho * (rho + 1))
    # 켤레 영점 쌍 포함: Σ_ρ (전체) = 양의 γ마다 2·Re → 계수 -2 × 2
    recon += -4 * np.real(coef * np.exp(1j * g0 * us))
recon = recon - np.polyval(np.polyfit(us, recon, 5), us)
corr = float(np.corrcoef(h, recon)[0, 1])
amp = float(np.std(h) / np.std(recon))
print("B. 영점 30개 이론 재구성 (자유 매개변수 0개 — 진폭·위상 모두 이론값)")
print(f"   상관계수 corr(실측, 재구성) = {corr:.4f}")
print(f"   진폭비 실측/이론 = {amp:.3f} (1이면 완벽)")
print("""
판정: corr ≳ 0.8 + 피크 일치 → '골드바흐 총합 요동 = 리만 영점 파동'의
수치 증명 완성. Q3' 전제(요동의 영점 기원) 확립 → 2부(분산 회계)로.""")
