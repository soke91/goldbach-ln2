"""S1-3 — chi2(n) 전수 포락선 (10^6 창, 3만 n).

chi2(n) = Sigma_q (F_q - a_q b_n)^2/(a_q b_n), b_n 자기정규화(n 내부 총합),
a_q = 보편 q-프로파일(전수 평균). 포아송이면 chi2/(K-1) ~ 1.
최악-n 포락선이 과제: max_n chi2/(K-1) — 작으면 chi2-슬롯 상수 K의 증거.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 1_000_000
W = 180_000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs_sieve = [int(q) for q in primes_upto(y) if q > 2]
qs_fiber = np.array([int(q) for q in primes_upto(1100) if q > y])
K = len(qs_fiber)

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
M = len(cands)
print(f"전수 {M:,}개, 섬유 {K}개", flush=True)

F = np.zeros((M, K), dtype=np.float64)
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
    if (i + 1) % 3000 == 0:
        print(f"  {i+1}/{M}", flush=True)

np.save("s1_chi2_F.npy", F)

# q|n 섬유는 구조적 영(별도 류) — 제외하고 각 n마다 유효 섬유만
mask = np.array([[n % q != 0 for q in qs_fiber] for n in cands])
a = np.where(mask, F, np.nan)
a_q = np.nanmean(a, axis=0)          # 보편 프로파일 (q nmid n 조건부)
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
print(f"전수 최악 {r.max():.4f} (n={cands[int(np.argmax(r))]})  "
      f"| 최소 {r.min():.4f}")
print(f"포아송 기준 1.0 — 최악조차 {r.max():.2f}면 sub-Poisson 균일성")
mu, sd = r.mean(), r.std()
print(f"최악 깊이 {(r.max()-mu)/sd:.2f}sigma | 가우시안 기대 최심 "
      f"{np.sqrt(2*np.log(M)):.2f}sigma")
np.savez("s1_chi2_data.npz", cands=np.array(cands), chi2=chi2, kk=kk)
