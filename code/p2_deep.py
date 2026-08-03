"""집중도 심층 — sub-binomial 비 정밀값 + 최소-s 포락선 (10⁷, 1500표본)."""

import numpy as np

from goldbach.sieve import primes_upto, sieve

SCALE = 10_000_000
NS = 1500
is_p = sieve(SCALE + 4_000_000)
ps_all = primes_upto(SCALE // 2 + 2_000_000)
y = int(round(SCALE ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]
cands = [int(v) for v in np.arange(SCALE, SCALE + 3_000_000, 2)
         if v % 6 != 0][:NS]
s_vals, S_sizes = [], []
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    s_vals.append(float(is_p[surv].mean()))
    S_sizes.append(len(surv))
    if (i + 1) % 300 == 0:
        print(f"  {i+1}/{NS}", flush=True)
s = np.array(s_vals)
S = np.array(S_sizes, dtype=float)
binom = np.sqrt(s.mean() * (1 - s.mean()) / S.mean())
ratio = s.std() / binom
se_ratio = ratio / np.sqrt(2 * NS)
print(f"\nn ~ 10^7, {NS}표본:")
print(f"s = {s.mean():.5f} ± {s.std():.5f}")
print(f"sub-binomial 비 = {ratio:.4f} ± {se_ratio:.4f}")
print(f"대조 상수: √ln2 = {np.sqrt(np.log(2)):.4f} | ln2 = 0.6931 | "
      f"1/(1+ln2) = 0.5906 | √(1/(1+ln2)) = {np.sqrt(1/(1+np.log(2))):.4f}")
print(f"\n최소-s 포락선: min = {s.min():.5f} (평균 대비 {(s.mean()-s.min())/s.std():.1f}σ)")
print(f"최대-s: {s.max():.5f} | s=0까지 {s.mean()/s.std():.0f}σ")
