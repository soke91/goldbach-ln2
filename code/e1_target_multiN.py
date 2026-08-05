# -*- coding: utf-8 -*-
"""
Increment 163: N-robustness of the E1 ratio (three more N, three bands each).

Target E1 asks Sum_{k ~ K} |D(k)|^2 << (log)^{-A} Sum_k M_k for
D(k) = Sum_{sqrt(N) < m <= N/k} mu(m) mu(N - mk). Under Conjecture L
the natural normalized ratio is
    E1_ratio(K) = Sum_k |D(k)|^2 / Sum_k support_k,
predicted 1.0 (unit Gaussian variance on support) -- E1's truth with
the (log)^{-A} margin corresponds to this ratio staying O(1) (it
needs only NOT to grow). Full segmented mu to N; three dyadic bands.
"""
import numpy as np, time

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.nonzero(s)[0]

def mobius_full(B, seg=20_000_000):
    bp = primes_upto(int(B ** 0.5) + 1)
    out = np.empty(B, dtype=np.int8)
    for lo in range(0, B, seg):
        hi = min(lo + seg, B)
        mu = np.ones(hi - lo, dtype=np.int8)
        val = np.arange(lo, hi, dtype=np.int64)
        if lo == 0:
            val[0] = 1  # avoid 0-division artifacts; mu[0] set later
        for p in bp:
            p = int(p)
            start = ((lo + p - 1) // p) * p - lo
            mu[start::p] *= -1
            val[start::p] //= p
            pp = p * p
            if pp < hi:
                start2 = ((lo + pp - 1) // pp) * pp - lo
                mu[start2::pp] = 0
        mu[val > 1] *= -1
        out[lo:hi] = mu
        if (lo // seg) % 10 == 9:
            print(f"  mu segment {hi/1e6:.0f}M", flush=True)
    out[0] = 0
    if B > 1:
        out[1] = 1
    return out

def main():
    rng = np.random.default_rng(20260823)
    NS = [899_999_998, 799_999_996, 699_999_994]
    t0 = time.time()
    for N in NS:
      print(f"=== N = {N} ===", flush=True)
      print("building full mu to N...", flush=True)
      mu = mobius_full(N + 1)
      print(f"mu ready {time.time()-t0:.0f}s", flush=True)
      SQ = int(N ** 0.5)

      for K0, K1, nk in ((1000, 2000, 150), (3000, 4000, 150),
                         (10000, 20000, 150)):
          ks = rng.choice(np.arange(K0, K1), size=nk, replace=False)
          S2 = 0.0; SS = 0.0
          r1s = []
          for i, k in enumerate(ks):
              k = int(k)
              ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
              t = mu[ms].astype(np.int64) * mu[N - k * ms]
              D = float(t.sum()); sup = int(np.count_nonzero(t))
              S2 += D * D; SS += sup
              r1s.append(abs(D) / np.sqrt(max(sup, 1)))
              if i % 50 == 49:
                  print(f"  band[{K0},{K1}) {i+1}/{nk}  "
                        f"t={time.time()-t0:.0f}s", flush=True)
          ratio = S2 / SS
          print(f"BAND K in [{K0},{K1}): E1_ratio = {ratio:.3f}  "
                f"(pred ~1.0)  mean r1 = {np.mean(r1s):.3f}  "
                f"({nk} k)", flush=True)
      del mu
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
