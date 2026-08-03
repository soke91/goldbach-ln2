"""chi2 스케일 3점째 — 10^8 틈새 창, 표준화 섬유 (y, 10y].

chi2/(K-1) 전형값 트렌드: 0.515(1e6) -> 0.539(1e7) -> ?(1e8).
K->1이면 chi2-슬롯 상수 우위 소멸, <1 정체면 증명 사슬 1고리 유지.
틈새 창 [463^3, 467^3] 내부라 체질 집합 상수 (함정 #21b 면역).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

q1 = 463
lo = 99_253_000
hi = 101_847_000
M = 2000
is_p = sieve(hi + 10)
ps_all = primes_upto(hi // 2 + 10)
qs_sieve = [int(q) for q in primes_upto(q1 + 1) if q > 2]
qs_fiber = np.array([int(q) for q in primes_upto(4640) if q > q1])
K = len(qs_fiber)

cands = [int(n) for n in range(lo, hi, 2) if n % 6 == 2]
rng = np.random.default_rng(37)
cands = sorted(rng.choice(cands, M, replace=False).tolist())
print(f"표본 {M}개, 섬유 {K}개 (q in ({q1}, 4640])", flush=True)

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
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{M}", flush=True)
    if (i + 1) % 500 == 0:
        np.save("s1_chi2_1e8_partial.npy", F[: i + 1])

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
print(f"표본 최악 {r.max():.4f}  | (1e6: 0.515 전수최악 0.796 / 1e7: 0.539"
      f" 최악 0.729)")
print(f"(1-mean)*lnX = {(1 - r.mean()) * np.log(lo):.2f}  "
      f"(1e6: 6.7, 1e7: 7.4 — R2 사가의 ~5.2와 비교)")
np.savez("s1_chi2_1e8_data.npz", cands=np.array(cands), chi2=chi2, kk=kk)
