# -*- coding: utf-8 -*-
"""
a3_margin_extract.py  --  pass4, blind mathematical re-verification.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.

 (J1) P1 Measurement 17, recomputed from the definitions in P1 eq:(3)
      and eq:(4) rather than from the paper's script:
        A(N) = prod_{q not| N} (1 - 1/(q(q-1)))
        S(N) = 2 prod_{p>2}(1-1/(p-1)^2) prod_{p|N,p>2}(1+1/(p-2))
      over the even N <= 1.6e7, reporting the median of S(1-A), its
      0.001 quantile, its minimum and argmin, and the minimum of
      S(1-A) log N loglog N over N >= 1e5, N >= 1e3, N >= 16.
      The paper prints, in order:
        0.333459, 0.119639, 0.060890 at N=9699690,
        2.482019 at 510510, 1.737799 at 2310, 0.810202 at 30.
      FALSIFIED by any disagreement in the sixth decimal, or by any
      argmin different from the one printed.

 (J2) P2 Measurement 17's ratio-7 table.  At N = 99999998,
      K = floor(N^0.56) = 30199, for the single-divisor weights
      w_k = [d0 | k] (so b = delta_{d0}), the brute-force
        B_w = sum_{k<K, (k,N)=1} mu(k) w_k / phi(k)
      against the factorised mu(d0)/phi(d0) * rho_{d0 N}(K/d0) of
      P2 Lemma 12, and the set of values of |B_w| phi(d0) over the
      d0 with floor(K/d0) = 7, split by gcd(d0,15).
      The paper prints 0.250000, 0.750000, 0.500000, 1.000000 for
      gcd = 1, 3, 5, 15, and says there are 184 such d0.
      FALSIFIED if the factorisation misses the brute force by more
      than 1e-12, or if the four values or the count 184 differ.

 (J3) P3 Measurement 3's clause "the figure for S is unaffected [by
      the lower cutoff]".  The shape residual sd(z_c/mean(z_c)) with
      z_c = (V/W)/c(N), c in {A, S}, over even N in [cutoff, 1.6e7],
      at cutoffs 1e5, 5e4, 1e4, 1e3, 1e2.  The paper prints the A
      column (0.000323, 0.000346, 0.000398, 0.000470, 0.000529) and
      asserts the S column is flat.
      FALSIFIED if the S column moves by more than 10% of its value
      across those cutoffs.
      NOTE: V(N) here is recomputed as sum_{v<N} mu^2(v) Lambda(N-v)^2
      by exact FFT, W(N) = sum_{w<N} Lambda(w)^2 by prefix sum, i.e.
      independently of the paper's script.

NULL.  J1 and J2 are deterministic arithmetic with no sign input; a
sign control would randomise nothing.  J3 is a shape comparison
between two fixed candidate constants, so the control is the second
candidate itself, which is what the statistic already contrasts.
"""
import numpy as np
import math

TOP = 16000000
print(__doc__.strip())
print()

s = np.ones(TOP + 1, dtype=bool)
s[:2] = False
for i in range(2, int(TOP ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
primes = np.nonzero(s)[0]

# --------------------------------------------------------------- J1
# Artin constant and 2C_2 as Euler products over p <= TOP
Artin = 1.0
twoC2 = 2.0
for p in primes:
    p = int(p)
    Artin *= (1.0 - 1.0 / float(p * (p - 1)))
    if p > 2:
        twoC2 *= (1.0 - 1.0 / float((p - 1) ** 2))

Afac = np.ones(TOP + 1)      # prod_{p|N} (1-1/(p(p-1)))^{-1}
Sfac = np.ones(TOP + 1)      # prod_{p|N,p>2} (1+1/(p-2))
for p in primes:
    p = int(p)
    Afac[p::p] *= 1.0 / (1.0 - 1.0 / float(p * (p - 1)))
    if p > 2:
        Sfac[p::p] *= (1.0 + 1.0 / float(p - 2))

Ne = np.arange(2, TOP + 1, 2)
Aval = Artin * Afac[Ne]
Sval = twoC2 * Sfac[Ne]
marg = Sval * (1.0 - Aval)

print("J1  P1 Measurement 17, recomputed")
print("    Artin = %.12f   2C_2 = %.12f" % (Artin, twoC2))
print("    median            = %.6f   paper 0.333459" % np.median(marg))
print("    quantile 0.001    = %.6f   paper 0.119639"
      % np.quantile(marg, 0.001))
i0 = int(np.argmin(marg))
print("    minimum           = %.6f at N = %d   paper 0.060890 at 9699690"
      % (marg[i0], Ne[i0]))
lg = np.log(Ne.astype(np.float64))
llg = np.log(np.maximum(lg, 1e-9))
prod = marg * lg * llg
for cut, pub in ((100000, 2.482019), (1000, 1.737799), (16, 0.810202)):
    sel = Ne >= cut
    j = int(np.argmin(prod[sel]))
    print("    min over N >= %-7d = %.6f at N = %-8d   paper %.6f"
          % (cut, prod[sel][j], Ne[sel][j], pub))
print("    ten smallest margins:",
      ", ".join("%d:%.6f" % (Ne[k], marg[k])
                for k in np.argsort(marg)[:10]))

# --------------------------------------------------------------- J2
print()
print("J2  P2 Measurement 17, the ratio-7 table")
Nb = 99999998
K = int(math.floor(Nb ** 0.56))
sk = np.ones(K + 1, dtype=bool)
sk[:2] = False
for i in range(2, int(K ** 0.5) + 1):
    if sk[i]:
        sk[i * i::i] = False
pk = np.nonzero(sk)[0]
muk = np.ones(K + 1, dtype=np.int64)
muk[0] = 0
phik = np.arange(K + 1, dtype=np.int64)
for p in pk:
    p = int(p)
    muk[p::p] *= -1
    phik[p::p] -= phik[p::p] // p
    if p * p <= K:
        muk[p * p::p * p] = 0
pfN = []
mm = Nb
for p in pk:
    p = int(p)
    if p * p > mm:
        break
    if mm % p == 0:
        pfN.append(p)
        while mm % p == 0:
            mm //= p
if mm > 1:
    pfN.append(mm)
print("    N = %d   prime factors %s   K = %d" % (Nb, pfN, K))
kk = np.arange(1, K)
copN = np.ones(K - 1, dtype=bool)
for p in pfN:
    if p <= K:
        copN[p - 1::p] = False
term = np.where(copN, muk[1:K] / phik[1:K].astype(float), 0.0)

worst = 0.0
d0list = []
for d0 in range(1, K):
    if muk[d0] == 0 or math.gcd(d0, Nb) != 1:
        continue
    d0list.append(d0)
brute_all = {}
for d0 in d0list:
    if K // d0 != 7:
        continue
    brute = float(term[d0 - 1::d0].sum())
    # factorised: mu(d0)/phi(d0) * sum_{j < K/d0, (j, d0 N)=1} mu(j)/phi(j)
    x = K / d0
    tot = 0.0
    for j in range(1, int(math.ceil(x))):
        if j < x and muk[j] != 0 and math.gcd(j, d0) == 1 \
                and math.gcd(j, Nb) == 1:
            tot += muk[j] / float(phik[j])
    fact = muk[d0] / float(phik[d0]) * tot
    worst = max(worst, abs(brute - fact))
    brute_all[d0] = abs(brute) * phik[d0]
byg = {}
for d0, val in brute_all.items():
    byg.setdefault(math.gcd(d0, 15), set()).add(round(val, 6))
print("    d0 with floor(K/d0) = 7 : %d   paper says 184" % len(brute_all))
print("    worst |brute - factorised| = %.3e" % worst)
for g in sorted(byg):
    print("      gcd(d0,15) = %-3d -> %s   paper %s"
          % (g, sorted(byg[g]),
             {1: 0.25, 3: 0.75, 5: 0.5, 15: 1.0}.get(g, "-")))

# --------------------------------------------------------------- J3
print()
print("J3  P3 Measurement 3, is the S column flat in the cutoff?")
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
lam2 = lam * lam
n = 1
while n < 2 * (TOP + 1):
    n <<= 1
V = np.fft.irfft(np.fft.rfft(mu2, n) * np.fft.rfft(lam2, n), n)[:TOP + 1]
W = np.cumsum(lam2)
ratio = np.zeros(TOP + 1)
ok = W[:TOP + 1] > 0
ratio[ok] = V[ok] / W[ok]
print("    %-9s %-12s %-12s" % ("cutoff", "sd for A", "sd for S"))
for cut in (100000, 50000, 10000, 1000, 100):
    sel = (Ne >= cut)
    r = ratio[Ne[sel]]
    out = []
    for cand in (Aval[sel], Sval[sel]):
        z = r / cand
        out.append(float(np.std(z / z.mean())))
    print("    %-9d %-12.6f %-12.6f" % (cut, out[0], out[1]))
