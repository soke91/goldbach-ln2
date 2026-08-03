"""v3 연장 — 10^9 틈새 점 (997^3 = 991,026,973 뒤, 다음 1009^3까지 틈새).

R2 -> 1 외삽 시험: 예측 R2(1e9) = 1 - 5.2/ln(1e9) = 1 - 5.2/20.7 = 0.749?
아니, v3 곡선의 (1-R2)lnX ~ 5.2 => 10^9에서 R2 ~ 0.75. 상승 지속 여부 판정.
메모리 ~1.6GB, 예상 ~1시간.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

q1 = 997
lo = 991_100_000
hi = 1_016_000_000
M = 300

print("sieve 구축...", flush=True)
is_p = sieve(hi + 10)
ps_all = primes_upto(hi // 2 + 10)
qs = [int(q) for q in primes_upto(q1 + 1) if q > 2]
print(f"준비 완료: 체질 소수 {len(qs)}개, p-풀 {len(ps_all):,}", flush=True)

cands = [int(n) for n in range(lo, hi, 2) if n % 6 == 2]
rng = np.random.default_rng(31)
cands = sorted(rng.choice(cands, M, replace=False).tolist())

s_vals = np.empty(M)
N_vals = np.empty(M)
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    N_vals[i] = len(surv)
    s_vals[i] = float(is_p[surv].mean())
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{M}", flush=True)
    if (i + 1) % 100 == 0:
        np.savez("s_v3_1e9_partial.npz", cands=np.array(cands[: i + 1]),
                 s=s_vals[: i + 1], N=N_vals[: i + 1])

s = s_vals
ln = np.log(np.array(cands, dtype=float))
A = np.vstack([np.ones_like(ln), ln - ln.mean()]).T
coef, *_ = np.linalg.lstsq(A, s, rcond=None)
resid = s - A @ coef
sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N_vals)))
R2_det = (resid.std() / sig_bin) ** 2
se = R2_det * np.sqrt(2 / (M - 1))
print(f"\nX~1e9 (y=997)  표본 {M}  s {s.mean():.5f}")
print(f"R2_det = {R2_det:.4f} +- {se:.4f}  (외삽 예측 ~0.75)")
print(f"(1-R2)*lnX = {(1 - R2_det) * np.log(lo):.2f}  (v3 하위 5점 ~5.2)")
np.savez("s_v3_1e9.npz", cands=np.array(cands), s=s_vals, N=N_vals)
