"""chi2 10^7 표준화 재측정 — 섬유 (y, 10y] = (215, 2150] (3점 트렌드 정합).

기존 10^7 측정은 (215, 3000] = 13.9y — 범위 교란 제거판.
틈새 창 [211^3, 223^3] 내부 (함정 #21b 면역).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

q1 = 211
lo = 9_394_100
hi = 9_957_700
M = 3000
is_p = sieve(hi + 10)
ps_all = primes_upto(hi // 2 + 10)
qs_sieve = [int(q) for q in primes_upto(q1 + 1) if q > 2]
qs_fiber = np.array([int(q) for q in primes_upto(2150) if q > q1])
K = len(qs_fiber)

cands = [int(n) for n in range(lo, hi, 2) if n % 6 == 2]
rng = np.random.default_rng(41)
cands = sorted(rng.choice(cands, M, replace=False).tolist())
print(f"표본 {M}개, 섬유 {K}개 (q in ({q1}, 2150])", flush=True)

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
      f"sigma {r.std():.4f}  최악 {r.max():.4f}")
print(f"(기존 13.9y 범위: 0.5392 / 표준화 10y 3점: 10⁶ 0.515?, 10⁷ 이 값, "
      f"10⁸ 0.5463)")
print(f"(1-mean)*lnX = {(1 - r.mean()) * np.log(lo):.2f}")
np.savez("s1_chi2_1e7_std.npz", cands=np.array(cands), chi2=chi2, kk=kk)
