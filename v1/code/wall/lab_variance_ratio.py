# -*- coding: utf-8 -*-
"""
Does the wall cancel below square root? (increment 257)

Increment 256 measured the off-diagonal of the variance identity
Sum_N |C(N)|^2 = Sum_h r_W(h) S_W(h) and found it negative and GROWING
as a fraction of the diagonal: |off|/diag = 0.083, 0.126, 0.172 at
X = 1e5, 1e6, 4e6. Equivalently the ratio

    Q(N) := Sum_N C(N)^2 / Sum_N V(N),   V(N) = Sum_v mu^2(v) Lambda(N-v)^2

runs 1.006, 0.955, 0.923, 0.890, 0.874, 0.873 across dyadic bands
(increment 238's column). V is the exact random-sign second moment, so
Q = 1 is exactly square-root cancellation.

TWO POSSIBILITIES, AND THEY ARE NOT THE SAME KIND OF FACT.
  (a) Q -> c for some 0 < c < 1. The wall cancels at square-root scale
      with a constant slightly better than random. Unremarkable.
  (b) Q -> 0. The wall cancels BELOW square root. That would be
      extraordinary: nothing in this campaign, and nothing known,
      produces sub-square-root cancellation in a mu-Lambda correlation,
      and Conjecture L asserts the opposite (Gaussian fluctuation at
      exactly the second-moment scale).

The measured sequence is consistent with both, which is precisely the
situation hazard 5 was named for. So the criterion is written in terms
of POWER and not trend.

WHAT IS MEASURED. Q by dyadic band up to X = 8e6, one band more than
before, with:
 * the sampling error of Q from the band itself (Q is a ratio of sums
   over ~n independent-ish terms, so its relative SE is estimated by
   the jackknife over halves of the band -- computed, not assumed);
 * three fits, Q = c, Q = a (log N)^-b, Q = a N^-b, compared by
   weighted residual on the same bands;
 * and the decisive line: the PREDICTED SEPARATION between the
   constant fit and the decaying fit at the largest N a computation
   could reach, against the measurement error there. If they separate
   by less than the error, the question is not numerically decidable
   and is to be labelled so rather than answered.

NULL. Q = 1 is exact square-root cancellation; the deviation from 1 is
what is being tracked, and its own error bar is printed beside it.
"""
import numpy as np
import math
import time


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return mu, lam


def main():
    X = 8_000_000
    lo = 50_000
    t0 = time.time()
    mu, lam = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    A = np.zeros(n_fft); B = np.zeros(n_fft)
    A[: X + 1] = mu; B[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    A[: X + 1] = np.abs(mu); B[: X + 1] = lam ** 2
    V = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    del A, B
    print(f"convolutions t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    Cv, Vv = C[Ns], V[Ns]

    print(f"\n{'band':>22} {'n':>9} {'Q':>9} {'jackknife SE':>13} "
          f"{'(Q-1)/SE':>10}")
    mids, Qs, SEs = [], [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 5000:
            c_, v_ = Cv[sel], Vv[sel]
            Q = float((c_ ** 2).sum() / v_.sum())
            # jackknife over 20 blocks: an honest SE for a ratio of sums
            m = 20
            idx = np.array_split(np.arange(len(c_)), m)
            qs = []
            for j in range(m):
                keep = np.ones(len(c_), bool); keep[idx[j]] = False
                qs.append(float((c_[keep] ** 2).sum() / v_[keep].sum()))
            qs = np.array(qs)
            se = math.sqrt((m - 1) / m * ((qs - qs.mean()) ** 2).sum())
            mid = math.sqrt(b * hi)
            mids.append(mid); Qs.append(Q); SEs.append(se)
            print(f"{b:>10}-{hi:>11} {int(sel.sum()):>9} {Q:>9.5f} "
                  f"{se:>13.5f} {(Q-1)/se:>10.1f}")
        b *= 2

    mids = np.array(mids); Qs = np.array(Qs); SEs = np.array(SEs)
    # The first bands have Q > 1, which is the pre-asymptotic regime:
    # V is the random-sign second moment, so Q > 1 means the wall is
    # LESS cancelled than random, which cannot be the asymptotic
    # behaviour. Fitting a global power law across a curve that is
    # still turning over manufactures a decay rate. Both fits are
    # reported: all bands, and only those with Q < 1.
    asym = Qs < 1.0
    print(f"\n  bands with Q < 1 (asymptotic regime): "
          f"{int(asym.sum())} of {len(Qs)}")
    print(f"  last three Q: " + ", ".join(f"{q:.5f}" for q in Qs[-3:])
          + f"   changes: "
          + ", ".join(f"{(Qs[i+1]-Qs[i])/SEs[i+1]:+.1f} SE"
                      for i in range(len(Qs) - 3, len(Qs) - 1)))
    w = 1.0 / SEs ** 2

    def wres(pred):
        return float((w * (Qs - pred) ** 2).sum())

    c0 = float((w * Qs).sum() / w.sum())
    r_const = wres(np.full_like(Qs, c0))
    L = np.log(mids)
    co1 = np.polyfit(np.log(L), np.log(Qs), 1)
    p1 = np.exp(co1[1]) * L ** co1[0]
    co2 = np.polyfit(np.log(mids), np.log(Qs), 1)
    p2 = np.exp(co2[1]) * mids ** co2[0]

    print(f"\n  fits, weighted residual (lower is better)")
    print(f"    Q = c              c = {c0:.5f}"
          f"                 res = {r_const:.4g}")
    print(f"    Q = a (log N)^-b   b = {-co1[0]:.4f}"
          f"                res = {wres(p1):.4g}")
    print(f"    Q = a N^-b         b = {-co2[0]:.5f}"
          f"               res = {wres(p2):.4g}")

    Lbig = math.log(1e12)
    q_const = c0
    q_log = math.exp(co1[1]) * Lbig ** co1[0]
    q_pow = math.exp(co2[1]) * (1e12) ** co2[0]
    if asym.sum() >= 4:
        ma, qa, sa = mids[asym], Qs[asym], SEs[asym]
        wa = 1.0 / sa ** 2
        ca = float((wa * qa).sum() / wa.sum())
        ra = float((wa * (qa - ca) ** 2).sum())
        La = np.log(ma)
        c1 = np.polyfit(np.log(La), np.log(qa), 1)
        pa = np.exp(c1[1]) * La ** c1[0]
        r1 = float((wa * (qa - pa) ** 2).sum())
        print(f"\n  the same fits on the Q < 1 bands only")
        print(f"    Q = c              c = {ca:.5f}"
              f"                 res = {ra:.4g}")
        print(f"    Q = a (log N)^-b   b = {-c1[0]:.4f}"
              f"                res = {r1:.4g}")
        q_log_a = math.exp(c1[1]) * math.log(1e12) ** c1[0]
        print(f"    extrapolated to 1e12:  const {ca:.4f}"
              f"   decaying {q_log_a:.4f}")

    print(f"\n  extrapolated to N = 1e12 (about the practical ceiling)")
    print(f"    constant fit  Q = {q_const:.4f}")
    print(f"    (log N)^-b    Q = {q_log:.4f}")
    print(f"    N^-b          Q = {q_pow:.4f}")
    print(f"    separation const vs (log N)^-b = "
          f"{abs(q_const-q_log):.4f}")
    print(f"    measurement error at the last band = {SEs[-1]:.5f}")
    if abs(q_const - q_log) < 5 * SEs[-1]:
        print("  => NOT NUMERICALLY DECIDABLE: the fits separate by less")
        print("     than a few times the error even at 1e12")
    else:
        print("  => decidable in principle at 1e12; the fits separate by")
        print(f"     {abs(q_const-q_log)/SEs[-1]:.1f} error bars")
    print("\n  note: Q = 1 is exactly square-root cancellation, so the")
    print("  (Q-1)/SE column says how far from random the wall sits,")
    print("  and the fits say whether that distance is growing")
    print("DONE")


if __name__ == "__main__":
    main()
