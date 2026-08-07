# -*- coding: utf-8 -*-
"""
Forge round-3, R4 power settlement (increment 194b).

R4 came in DEAD on its pre-registered criterion (no block ratio at or
below 0.5 x the B=1 baseline at either N). But the D-field showed a
mild deficit in the SAME direction at both N (B=8: 0.87 / 0.94;
B=64: 0.77 / 0.80 against baselines 0.96 / 0.99) on only 2048 k, where
the block-ratio estimator's own SE is ~9% at B=8 and ~25% at B=64.
This program's history contains a +3.7 sigma band that regressed to
noise under power, so the deficit is settled rather than reported.

Power increase: the type-II field D(k) = Sum_{sqrt N < m <= N/k}
is non-empty only for k < sqrt(N), so the k-range is widened to the
largest available band, [1500, 9500) = 8000 k. At B = 8 that is 1000
blocks (SE ~ 4.5%), at B = 64 it is 125 blocks (SE ~ 12.6%).

Added diagnostic (better SE than any block ratio): the lag-l
autocorrelation of the normalised field d(k) = D(k)/sqrt(supp(k)),
l = 1..4. Block deficits and negative autocorrelation are the same
phenomenon seen two ways; with 8000 k the autocorrelation has
SE ~ 1/sqrt(8000) = 0.011, so a real coherence of the size the block
ratios hint at (a ~10% variance deficit at B=8 needs r1 ~ -0.06)
would show at ~5 sigma.

PRE-REGISTERED (fixed before the run):
  - ALIVE criterion is unchanged and remains ratio(B) <= 0.5*ratio(1).
  - ANOMALY CONFIRMED (not alive, but real and worth recording) iff
    ratio(B=8) <= 0.85 at BOTH N *and* lag-1 autocorrelation
    r1 <= -0.04 at BOTH N.
  - NOISE otherwise: the increment-194 deficit is an artifact of the
    block-count SE and R4 closes clean.
"""
import numpy as np
import time

from e1_forge_r4 import mobius_upto, field


def stats(D, S, blocks, maxlag=4):
    out = {}
    for B in blocks:
        nb = len(D) // B
        s = D[:nb * B].reshape(nb, B).sum(axis=1)
        used = float(S[:nb * B].sum())
        out[B] = (float(np.sum(s ** 2)) / used, nb)
    d = D / np.sqrt(np.maximum(S, 1.0))
    d = d - d.mean()
    v = float(np.dot(d, d))
    out['acf'] = [float(np.dot(d[:-l], d[l:]) / v) for l in range(1, maxlag + 1)]
    out['n'] = len(D)
    return out


def main():
    BLOCKS = [1, 8, 64, 512]
    K0, K1 = 1500, 9500
    ks = np.arange(K0, K1)
    verdict = {}

    for N in (99_999_998, 99_960_002):
        t0 = time.time()
        mu = mobius_upto(N)
        print(f"=== N = {N}  (mu ready {time.time()-t0:.0f}s) ===",
              flush=True)
        print(f"k band [{K0}, {K1}) = {len(ks)} k, sqrt(N) = "
              f"{int(N**0.5)}", flush=True)
        D, S = field(mu, N, ks, False)
        r = stats(D, S, BLOCKS)
        line = "  ".join(f"B={B}: {r[B][0]:.4f} ({r[B][1]} blk)"
                         for B in BLOCKS)
        print(f"  D (m > sqrt N)   {line}", flush=True)
        se = 1.0 / np.sqrt(len(ks))
        acf_s = "  ".join(f"r{l+1}={v:+.4f}" for l, v in
                          enumerate(r['acf']))
        print(f"  autocorrelation  {acf_s}   (SE ~ {se:.4f}, "
              f"z1 = {r['acf'][0]/se:+.2f})", flush=True)
        print(f"  [mean D = {D.mean():.2f}, t = {time.time()-t0:.0f}s]",
              flush=True)
        verdict[N] = (r[8][0], r['acf'][0])
        del mu, D, S

    print("\n=== R4b SETTLEMENT ===", flush=True)
    alive = all(v[0] <= 0.5 for v in verdict.values())
    anom = (all(v[0] <= 0.85 for v in verdict.values())
            and all(v[1] <= -0.04 for v in verdict.values()))
    for N, (r8, a1) in verdict.items():
        print(f"  N={N}: ratio(B=8) = {r8:.4f}, r1 = {a1:+.4f}",
              flush=True)
    print("verdict:",
          "ALIVE" if alive else
          ("ANOMALY CONFIRMED (real coherence, below the alive bar)"
           if anom else
           "NOISE -- the increment-194 deficit is block-count SE; "
           "R4 closes clean"), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
