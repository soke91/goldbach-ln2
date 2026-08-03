"""S1-4 — chi2/(K-1)의 스케일 안정성 (10^7, 표본 3000).

10^6 전수: 전형 0.515, 최악 0.796. 10^7에서 전형/최악(표본)이 유지되는가.
섬유 q in (X^{1/3}, 3000] — K ~ 320개.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 10_000_000
W = 600_000
M = 3000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs_sieve = [int(q) for q in primes_upto(y) if q > 2]
qs_fiber = np.array([int(q) for q in primes_upto(3000) if q > y])
K = len(qs_fiber)

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
rng = np.random.default_rng(13)
cands = sorted(rng.choice(cands, M, replace=False).tolist())
print(f"표본 {M}개, 섬유 {K}개 (q in ({y}, 3000])", flush=True)

F = np.zeros((M, K))
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs_sieve:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    comp = surv[~is_p[surv]]
    for j, q in enumerate(qs_fiber):
        F[i, j] = np.count_nonzero(comp % q == 0)
    if (i + 1) % 300 == 0:
        print(f"  {i+1}/{M}", flush=True)

mask = np.array([[n % q != 0 for q in qs_fiber] for n in cands])
a = np.where(mask, F, np.nan)
a_q = np.nanmean(a, axis=0)
chi2 = np.empty(M)
kk = np.empty(M)
for i in range(M):
    mk = mask[i]
    e_shape = a_q[mk]
    b = F[i, mk].sum() / e_shape.sum()
    E = e_shape * b
    chi2[i] = float(((F[i, mk] - E) ** 2 / E).sum())
    kk[i] = mk.sum() - 1

r = chi2 / kk
print(f"\nchi2/(K-1): 평균 {r.mean():.4f}  중앙값 {np.median(r):.4f}  "
      f"sigma {r.std():.4f}")
print(f"표본 최악 {r.max():.4f}  | 최소 {r.min():.4f}")
print(f"(10^6 전수: 전형 0.515+-0.062, 최악 0.796)")
np.savez("s1_chi2_scale_data.npz", cands=np.array(cands), chi2=chi2, kk=kk)
