"""chi2 4점째 — 10^9 틈새 창 (997^3 뒤), 섬유 (y, 10y] = (997, 9970].

3점 감속 트렌드의 극한-미만-1 판정. 조잡 적합(하위 2점 쌍)의 외삽 예측:
chi2(1e9) ~ 0.55 부근 (a=0.597, b=0.936 모형이면 0.552).
포아송 수렴(1-c/lnX)이면 ~0.60+ 로 점프해야 함 — 갈림길 시험.
부분 저장으로 세션-단절 내성. 예상 3~4.5h.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

q1 = 997
lo = 991_100_000
hi = 1_016_000_000
M = 400

print("sieve 구축...", flush=True)
is_p = sieve(hi + 10)
ps_all = primes_upto(hi // 2 + 10)
qs_sieve = [int(q) for q in primes_upto(q1 + 1) if q > 2]
qs_fiber = np.array([int(q) for q in primes_upto(9970) if q > q1])
K = len(qs_fiber)
print(f"표본 {M}개, 섬유 {K}개 (q in ({q1}, 9970])", flush=True)

cands = [int(n) for n in range(lo, hi, 2) if n % 6 == 2]
rng = np.random.default_rng(43)
cands = sorted(rng.choice(cands, M, replace=False).tolist())

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
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{M}", flush=True)
    if (i + 1) % 50 == 0:
        np.savez("s1_chi2_1e9_partial.npz", F=F[: i + 1],
                 cands=np.array(cands[: i + 1]))

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
print(f"3점 참조: 0.5153/0.5392/0.5463 — 극한<1이면 ~0.55, 포아송행이면 0.60+")
print(f"(1-mean)*lnX = {(1 - r.mean()) * np.log(lo):.2f}  (6.7/7.4/8.35)")
np.savez("s1_chi2_1e9_data.npz", cands=np.array(cands), chi2=chi2, kk=kk)
