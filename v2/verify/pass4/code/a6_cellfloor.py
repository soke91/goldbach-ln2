# -*- coding: utf-8 -*-
"""
a6_cellfloor.py  --  pass4, blind mathematical re-verification.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.  P4 Lemma 1 and the numbers P4 Measurements 8 and 10
print, re-implemented here from the statements alone.

  Band: the even N in (2e6, 4e6] with V(N) > 0.
  Cells: depth d(N) = #{p in 3,5,7,11,13 : p | N}, six cells d=0..5.
  V(N)  = sum_{v<N} mu^2(v) Lambda(N-v)^2       (exact FFT convolution)
  C(N)  = sum_{n<N} Lambda(n) mu(N-n)           (exact FFT convolution)
  Z(N)  = C(N)/sqrt(V(N))
  u_c(v)= sum_{N in c} Lambda(N-v)/sqrt(V(N))   (exact FFT correlation)
  Q_cd  = sum_v mu^2(v) u_c(v) u_d(v)
  se_c  = sqrt( Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2 )
  z_c   = (m_c - mbar)/se_c

  and, for P4 Note 6's clause "at the top octave the exact floor
  exceeds sd(Z)/sqrt(n_c) by factors of 5.8 to 160, growing with the
  cell", the ratios se_c / (sd(Z)/sqrt(n_c)).

WHAT WOULD FALSIFY WHAT.
  (a) P4 Measurement 8's printed n_c, se_c and z_c for this octave are
      383617/421978/163568/28507/2263/67,
      1.2400e-1/4.6622e-2/1.4834e-1/2.4178e-1/3.3253e-1/4.3685e-1,
      +0.1069/+1.1133/-0.2314/-2.3990/-5.9997/-11.0258.
      FALSIFIED if any recomputed value differs in the fourth
      significant digit.
  (b) The "5.8 to 160" clause is CONFIRMED only if the recomputed
      ratio range brackets those two endpoints to within 10%; it is
      an evidence gap, not a defect, if it does not, since the
      packet prints no sd(Z) anywhere.

NULL.  The floor of Lemma 1 IS the independent-sign variance, computed
in closed form rather than simulated, so every z_c below is already
quoted against its own coin null; no further control is run here and
none of the thresholds above is set from an effect size.
"""
import numpy as np
import math

import sys
TOP = int(sys.argv[1]) if len(sys.argv)>1 else 4000000
LO = TOP//2
print(__doc__.strip())
print()

s = np.ones(TOP + 1, dtype=bool)
s[:2] = False
for i in range(2, int(TOP ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
primes = np.nonzero(s)[0]
mu = np.ones(TOP + 1, dtype=np.int64)
mu[0] = 0
lam = np.zeros(TOP + 1)
for p in primes:
    p = int(p)
    lp = math.log(p)
    mu[p::p] *= -1
    if p * p <= TOP:
        mu[p * p::p * p] = 0
    q = p
    while q <= TOP:
        lam[q] = lp
        q *= p
mu2 = (mu != 0).astype(np.float64)

n = 1
while n < 2 * (TOP + 1):
    n <<= 1
FL = np.fft.rfft(lam, n)
V = np.fft.irfft(np.fft.rfft(mu2 * 1.0, n) * np.fft.rfft(lam * lam, n), n)[:TOP + 1]
C = np.fft.irfft(np.fft.rfft(mu.astype(np.float64), n) * FL, n)[:TOP + 1]

Ne = np.arange(LO + 2, TOP + 1, 2)
Ne = Ne[V[Ne] > 0]
Z = C[Ne] / np.sqrt(V[Ne])
depth = np.zeros(TOP + 1, dtype=np.int64)
for p in (3, 5, 7, 11, 13):
    depth[p::p] += 1
dc = depth[Ne]
nband = len(Ne)
mbar = float(Z.mean())

# u_c(v) for each cell and for the whole band
def uvec(mask):
    G = np.zeros(TOP + 1)
    G[Ne[mask]] = 1.0 / np.sqrt(V[Ne[mask]])
    return np.fft.irfft(np.fft.rfft(G, n) * np.conj(FL), n)[:TOP + 1]

ua = uvec(np.ones(nband, dtype=bool))
Qaa = float(np.dot(mu2, ua * ua))
sdZ = float(Z.std(ddof=1))
print("  band (%d, %d]:  n = %d   mean Z = %+.6f   sd Z = %.6f"
      % (LO, TOP, nband, mbar, sdZ))
print()
print("  %-6s %-9s %-13s %-13s %-11s %-9s %-9s"
      % ("depth", "n_c", "se_c", "paper se_c", "z_c", "paper z_c", "se/(sdZ/sqrt n_c)"))
paper_se = ([1.2400e-1, 4.6622e-2, 1.4834e-1, 2.4178e-1, 3.3253e-1, 4.3685e-1]
              if TOP==4000000 else [1.1823e-1,4.4452e-2,1.4144e-1,2.3053e-1,3.1686e-1,4.1153e-1])
paper_z = ([+0.1069, +1.1133, -0.2314, -2.3990, -5.9997, -11.0258]
             if TOP==4000000 else [+0.0440,+0.6905,-0.0561,-1.5529,-4.4589,-9.0642])
ratios = []
for d in range(6):
    mask = (dc == d)
    nc = int(mask.sum())
    uc = uvec(mask)
    Qcc = float(np.dot(mu2, uc * uc))
    Qca = float(np.dot(mu2, uc * ua))
    var = Qcc / nc ** 2 - 2.0 * Qca / (nc * nband) + Qaa / nband ** 2
    se = math.sqrt(var)
    mc = float(Z[mask].mean())
    z = (mc - mbar) / se
    r = se / (sdZ / math.sqrt(nc))
    ratios.append(r)
    print("  %-6d %-9d %-13.5e %-13.5e %-11.4f %-9.4f %.2f"
          % (d, nc, se, paper_se[d], z, paper_z[d], r))
print()
print("  P4 Note 6's clause: floor / (sd(Z)/sqrt n_c) ranges %.2f to %.2f here"
      % (min(ratios), max(ratios)))
print("  (P4 Note 6 says 5.8 to 160 at the TOP octave (8e6,1.6e7], one octave above this one)")
