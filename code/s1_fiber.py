"""S1 — 반소수 생존자의 q-섬유 스펙트럼과 교차-섬유 상관.

theta=1/3 지점. 각 n의 생존자 중 합성수(=반소수, 정리 1)는 m = q*r,
y < q <= sqrt(m). 섬유 F_q(n) = #{생존자: q | m}. 측정:
  (1) 집중 이득 G = Var_n(총합)/Sigma_q Var_n(F_q) — 1이면 독립,
      <<1이면 강한 음의 교차상관 = 전역 상계의 여지.
  (2) 섬유 평균의 q-형상 (1/q 스케일링 확인).
  (3) 총합(적 총수)의 상대 산포 vs 개별 섬유의 상대 산포.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 1_000_000
W = 120_000
M = 2000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs_sieve = [int(q) for q in primes_upto(y) if q > 2]
qs_fiber = [int(q) for q in primes_upto(1100) if q > y]

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
rng = np.random.default_rng(5)
cands = sorted(rng.choice(cands, M, replace=False).tolist())

F = np.zeros((M, len(qs_fiber)))
comp_tot = np.zeros(M)
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs_sieve:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    comp = surv[~is_p[surv]]
    comp_tot[i] = len(comp)
    for j, q in enumerate(qs_fiber):
        F[i, j] = np.count_nonzero(comp % q == 0)
    if (i + 1) % 200 == 0:
        print(f"  {i+1}/{M}", flush=True)

tot = F.sum(axis=1)
cover = float(np.mean(tot / comp_tot))
print(f"\n표본 {M}  섬유 q in ({y}, 1100]  {len(qs_fiber)}개")
print(f"섬유 합/합성 총수 커버리지: {cover:.4f} (반소수 m=q*r, q<=sqrt(m)"
      f" 커버 시 ~1; q>1100 소인수 반소수는 미커버)")

var_tot = tot.var()
sum_var = F.var(axis=0).sum()
G = var_tot / sum_var
print(f"\n집중 이득 G = Var(총합)/SigmaVar(섬유) = {var_tot:.1f}/{sum_var:.1f}"
      f" = {G:.3f}")
print(f"  (독립이면 1, 음의 교차상관이면 <1)")

# 총합 vs 섬유 상대 산포
rel_tot = tot.std() / tot.mean()
rel_fib = float(np.mean(F.std(axis=0) / np.maximum(F.mean(axis=0), 1e-9)))
print(f"상대 산포: 총합 {rel_tot:.4f} vs 섬유 평균 {rel_fib:.4f} "
      f"(비 {rel_fib/rel_tot:.1f}배)")

# q-형상: 평균 F_q * q / 총합평균 — 1/q 스케일이면 평탄
shape = F.mean(axis=0) * np.array(qs_fiber) / tot.mean()
print(f"\nq-형상 (F_q*q/tot, 1/q면 평탄): 앞 {shape[:5].round(3)} "
      f"중간 {shape[len(shape)//2:len(shape)//2+5].round(3)} "
      f"끝 {shape[-5:].round(3)}")

# 이웃 섬유 상관 요약
C = np.corrcoef(F.T)
off = C[np.triu_indices_from(C, k=1)]
print(f"교차-섬유 상관: 평균 {off.mean():+.4f}  중앙값 "
      f"{np.median(off):+.4f}  범위 [{off.min():+.3f}, {off.max():+.3f}]")

np.savez("s1_fiber_data.npz", cands=np.array(cands), F=F,
         comp=comp_tot, qs=np.array(qs_fiber))
