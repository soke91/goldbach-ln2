"""A5-3 캠페인 v2 — 오염 제거 프로토콜: n별 y = n^{1/3}, 고정 W/X, 디트렌드.

함정 #21 교정: 창-고정 y가 만드는 P3 오염(n^{1/3} > y 구간)이 cross-n
분산을 부풀림 (W/X=0.3에서 R2 2.24 관측). v2:
  - 체질 소수를 n마다 q <= n^{1/3}로 정확히.
  - W/X = 0.06 전 스케일 통일.
  - sigma_obs는 ln n 선형 디트렌드 잔차로.
스케일별 R2_raw / R2_detrended 보고, npz 저장.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

for X, M in [(100_000, 1000), (1_000_000, 3000), (4_000_000, 3000),
             (10_000_000, 3000), (40_000_000, 1500), (100_000_000, 1000)]:
    W = int(X * 0.06)
    is_p = sieve(X + W + 10)
    ps_all = primes_upto((X + W) // 2 + 10)
    qs_full = np.array([int(q) for q in primes_upto(
        int(round((X + W) ** (1 / 3))) + 2) if q > 2])
    q3 = qs_full.astype(np.int64) ** 3

    cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
    rng = np.random.default_rng(23)
    if len(cands) > M:
        cands = sorted(rng.choice(cands, M, replace=False).tolist())

    s_vals = np.empty(len(cands))
    N_vals = np.empty(len(cands))
    for i, n in enumerate(cands):
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs_full[q3 <= n]:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        N_vals[i] = len(surv)
        s_vals[i] = float(is_p[surv].mean())
        if (i + 1) % 500 == 0:
            print(f"  X={X:,}: {i+1}/{len(cands)}", flush=True)

    s = s_vals
    ln = np.log(np.array(cands, dtype=float))
    A = np.vstack([np.ones_like(ln), ln - ln.mean()]).T
    coef, *_ = np.linalg.lstsq(A, s, rcond=None)
    resid = s - A @ coef
    sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N_vals)))
    R2_raw = (s.std() / sig_bin) ** 2
    R2_det = (resid.std() / sig_bin) ** 2
    se = R2_det * np.sqrt(2 / (len(s) - 1))
    print(f"X={X:>13,} 표본 {len(s):>5}  s {s.mean():.5f}  "
          f"R2_raw {R2_raw:.4f}  R2_det {R2_det:.4f} +- {se:.4f}  "
          f"drift/efold {coef[1]:+.5f}", flush=True)
    np.savez(f"s_v2_{X}.npz", cands=np.array(cands), s=s_vals, N=N_vals)

print("\n(v1 오염치 참고: 0.820/0.728/0.731/0.708/0.783/0.802)")
