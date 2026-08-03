"""A5-3 캠페인 v3 (최종 프로토콜) — 세제곱-틈새 창.

함정 #21b: n별 y도 창 내 q^3 경계에서 s가 ~1-2 sigma 계단 (q 활성화가
합성 생존자 제거) -> 분산 오염. v3: 창을 연속 소수 세제곱 사이 틈새
[q1^3, q2^3]에 배치 -> 체질 집합 상수, 계단 원천 제거.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

# (q1, 창시작, 창끝, 표본수) — 창 [시작,끝] ⊂ (q1^3, q2^3)
SCALES = [
    (47, 103_900, 110_130, 1000),
    (101, 1_030_400, 1_092_600, 3000),
    (157, 3_870_000, 4_102_200, 3000),
    (211, 9_394_100, 9_957_700, 3000),
    (337, 38_273_000, 40_569_000, 1500),
    (463, 99_253_000, 101_847_000, 1000),
]

for q1, lo, hi, M in SCALES:
    is_p = sieve(hi + 10)
    ps_all = primes_upto(hi // 2 + 10)
    qs = [int(q) for q in primes_upto(q1 + 1) if q > 2]
    assert qs[-1] == q1 or q1 == 2

    cands = [int(n) for n in range(lo, hi, 2) if n % 6 == 2]
    rng = np.random.default_rng(29)
    if len(cands) > M:
        cands = sorted(rng.choice(cands, M, replace=False).tolist())

    s_vals = np.empty(len(cands))
    N_vals = np.empty(len(cands))
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
        if (i + 1) % 500 == 0:
            print(f"  X~{lo:,}: {i+1}/{len(cands)}", flush=True)

    s = s_vals
    ln = np.log(np.array(cands, dtype=float))
    A = np.vstack([np.ones_like(ln), ln - ln.mean()]).T
    coef, *_ = np.linalg.lstsq(A, s, rcond=None)
    resid = s - A @ coef
    sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N_vals)))
    R2_raw = (s.std() / sig_bin) ** 2
    R2_det = (resid.std() / sig_bin) ** 2
    se = R2_det * np.sqrt(2 / (len(s) - 1))
    print(f"X~{lo:>12,} (y={q1})  표본 {len(s):>5}  s {s.mean():.5f}  "
          f"R2_raw {R2_raw:.4f}  R2_det {R2_det:.4f} +- {se:.4f}  "
          f"drift {coef[1]:+.5f}", flush=True)
    np.savez(f"s_v3_{lo}.npz", cands=np.array(cands), s=s_vals, N=N_vals)
