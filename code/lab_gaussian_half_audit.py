# -*- coding: utf-8 -*-
"""
Is Conjecture L's Gaussian half actually true for the wall? (inc. 283)

WHAT IS CLAIMED. CONJECTURE_L.md: on the surviving support the
fluctuation is "**exactly Gaussian** at half-normal scale: mean 0.798,
variance ratio 1, kurtosis 3, no mean field, no class structure."
MEASUREMENTS section 12-13 extended that claim to the wall's own scalar
C(N), which is the object the entire chain reduces to.

WHY IT IS BEING TESTED NOW. Increment 281 measured, as a by-product,
E|C|/sd = 0.79145 against the Gaussian sqrt(2/pi) = 0.79788 -- a
consistent -0.81% across all eight bands, with 5*10^4 to 3.2*10^6
samples each. That is not noise. Increment 240 separately reported
excess kurtosis +0.014 after masking, which reads as "essentially 3"
but has standard error sqrt(24/n) ~ 0.005 at these n, i.e. ~3 sigma.
Two independent signs of non-Gaussianity, neither followed up.

THE QUESTION THAT DECIDES WHAT IT MEANS. A residual deviation can be

  (i) MASK LEAKAGE -- the enumeration stops at q <= 23 and the primes
      above it still carry a location mask, or
  (ii) A VARIANCE MASK -- the cells have different SCALES, not just
      different means. Pooling cells of unequal variance produces a
      mixture, and a mixture of Gaussians has POSITIVE excess kurtosis
      and E|X|/sd BELOW sqrt(2/pi). Both observed signs match this.
      Conjecture L's mask list has a location mask and a scale mask
      sqrt(S(N)), but S is already divided out in Z, so a residual
      scale mask would be a NEW mask the conjecture does not name.
  (iii) INTRINSIC non-Gaussianity, in which case the conjecture is
      wrong as stated rather than incomplete.

These are distinguishable and the script distinguishes them.

PRE-REGISTRATION (fixed before the run).

  A. Baseline, per octave band, de-masked on cells keyed by which of
     {3,...,23} divide N: skewness, excess kurtosis, E|X|/sd, each
     with its Gaussian standard error (sqrt(6/n), sqrt(24/n),
     sqrt((1-2/pi)/n)). Deviations are reported as z, not as "small".

  B. DEPTH LADDER. Repeat with the cell key built from the first d
     odd primes, d = 2, 4, 6, 8, 10. If the deviation falls toward
     zero with depth, it is (i) mask leakage. Cell counts and the
     minimum cell occupancy are printed, because removing more means
     from the same n shrinks variance by construction -- every
     estimate carries the exact (n-k) correction.

  C. SCALE DE-MASKING. Repeat with per-cell STANDARDISATION (subtract
     the cell mean AND divide by the cell sd) instead of centring
     only. If the deviation collapses here but not in B, it is (ii) a
     variance mask.
     DECISION RULE, fixed now: (ii) is declared if per-cell
     standardisation removes at least 70% of the excess kurtosis while
     the depth ladder at fixed centring removes less than 30%.

  D. N-TREND. Does the deviation shrink with N? "Exactly Gaussian" is
     an asymptotic claim, so a deviation decaying with N leaves the
     conjecture's spirit intact; one that is flat or growing does not.

  E. ROBUSTNESS. Octave bands mix scales within the band (sd varies
     by ~1.6% across a factor 2 in N), which is itself a variance
     mixture. Repeated on quarter-octave bands to show this is not the
     explanation -- the predicted contribution is ~0.001 in excess
     kurtosis, far below what is measured, but predicted is not
     measured.

WHAT THIS CANNOT DO. It tests the extension of Conjecture L to C(N),
not the original claim for the mu-families D(k), which was measured on
different objects and is not touched here.
"""
import numpy as np
import math
import time

ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam, primes


def wall(X, mu, lam):
    n = 1
    while n < 2 * (X + 1):
        n *= 2
    A = np.zeros(n); A[: X + 1] = mu
    B = np.zeros(n); B[: X + 1] = lam
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n)[: X + 1]


def cellkey(Ns, d):
    k = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(ODD_PRIMES[:d]):
        k |= ((Ns % q) == 0).astype(np.int64) << i
    return k


def stats(z, ncells):
    """moments with the exact dof correction for removed cell means."""
    n = len(z)
    dof = n - ncells
    v = float((z ** 2).sum()) / dof
    s = math.sqrt(v)
    zz = z / s
    sk = float((zz ** 3).mean())
    ku = float((zz ** 4).mean()) - 3.0
    ea = float(np.abs(zz).mean())
    return s, sk, ku, ea, n


def report(tag, z, ncells):
    s, sk, ku, ea, n = stats(z, ncells)
    se_sk = math.sqrt(6.0 / n)
    se_ku = math.sqrt(24.0 / n)
    se_ea = math.sqrt((1.0 - 2.0 / math.pi) / n)
    d_ea = ea - math.sqrt(2.0 / math.pi)
    return dict(tag=tag, n=n, cells=ncells, skew=sk, z_sk=sk / se_sk,
                kurt=ku, z_ku=ku / se_ku, ea=ea, d_ea=d_ea,
                z_ea=d_ea / se_ea)


def bands(Ns, lo, X, ratio):
    b = lo
    out = []
    while b < X:
        hi = min(int(b * ratio), X)
        if hi <= b:
            break
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 5000:
            out.append((b, hi, sel))
        b = hi
    return out


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    C = wall(X, mu, lam)
    print(f"sieve + convolution  t={time.time()-t0:.0f}s", flush=True)

    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    Ns = np.arange(lo, X + 1, 2)
    Z = C[Ns] / np.sqrt(S[Ns] * Ns)
    bs = bands(Ns, lo, X, 2.0)

    def demask(v, key, scale=False, minc=30):
        uniq, inv = np.unique(key, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
        tot = np.bincount(inv, weights=v, minlength=len(uniq))
        m = tot / cnt
        r = v - m[inv]
        if not scale:
            return r, len(uniq)
        ss = np.bincount(inv, weights=r * r, minlength=len(uniq))
        sd = np.sqrt(ss / np.maximum(cnt - 1.0, 1.0))
        ok = (cnt >= minc) & (sd > 0)
        keep = ok[inv]
        out = r[keep] / sd[inv][keep]
        return out, int(ok.sum())

    print("\n(A) baseline: de-masked on {3..23}, per octave band")
    print(f"{'band':>21} {'n':>9} {'skew':>8} {'z':>7} "
          f"{'exc.kurt':>9} {'z':>8} {'E|X|/sd':>9} {'z':>8}")
    baseA = []
    for b, hi, sel in bs:
        z = Z[sel]
        r, k = demask(z, cellkey(Ns[sel], 8))
        d = report(f"{b}", r, k)
        baseA.append(d)
        print(f"{b:>9}-{hi:>11} {d['n']:>9} {d['skew']:>8.4f} "
              f"{d['z_sk']:>7.1f} {d['kurt']:>9.4f} {d['z_ku']:>8.1f} "
              f"{d['ea']:>9.5f} {d['z_ea']:>8.1f}")
    print(f"    Gaussian: skew 0, excess kurtosis 0, "
          f"E|X|/sd = {math.sqrt(2/math.pi):.5f}")

    print("\n(B) depth ladder, centring only -- is it mask leakage?")
    print(f"{'depth':>6} {'primes up to':>13} {'cells':>7} "
          f"{'exc.kurt':>10} {'z':>8} {'E|X|/sd dev':>13} {'z':>8}")
    kur_by_depth = {}
    for d_ in (2, 4, 6, 8, 10):
        rs, ks = [], 0
        for b, hi, sel in bs:
            r, k = demask(Z[sel], cellkey(Ns[sel], d_))
            rs.append(r / r.std())
            ks += k
        allr = np.concatenate(rs)
        dd = report("d", allr, ks)
        kur_by_depth[d_] = dd['kurt']
        print(f"{d_:>6} {ODD_PRIMES[d_-1]:>13} {ks:>7} "
              f"{dd['kurt']:>10.4f} {dd['z_ku']:>8.1f} "
              f"{dd['d_ea']:>13.5f} {dd['z_ea']:>8.1f}")

    print("\n(C) per-cell STANDARDISATION -- is it a variance mask?")
    print(f"{'depth':>6} {'cells kept':>11} {'exc.kurt':>10} {'z':>8} "
          f"{'E|X|/sd dev':>13} {'z':>8}")
    kur_std = {}
    for d_ in (4, 6, 8):
        rs, ks = [], 0
        for b, hi, sel in bs:
            r, k = demask(Z[sel], cellkey(Ns[sel], d_), scale=True)
            rs.append(r)
            ks += k
        allr = np.concatenate(rs)
        dd = report("s", allr, ks)
        kur_std[d_] = dd['kurt']
        print(f"{d_:>6} {ks:>11} {dd['kurt']:>10.4f} {dd['z_ku']:>8.1f} "
              f"{dd['d_ea']:>13.5f} {dd['z_ea']:>8.1f}")

    k0 = kur_by_depth[4]
    drop_depth = 1.0 - kur_by_depth[10] / k0 if k0 else float('nan')
    drop_scale = 1.0 - kur_std[8] / kur_by_depth[8] if kur_by_depth[8] else float('nan')
    print(f"\n    excess kurtosis removed by DEPTH  (d=4 -> d=10): "
          f"{drop_depth:.1%}")
    print(f"    excess kurtosis removed by SCALING (d=8):        "
          f"{drop_scale:.1%}")
    print("    pre-registered rule: a variance mask is declared if")
    print("    scaling removes >= 70% while depth removes < 30%")
    verdict = ("(ii) VARIANCE MASK" if (drop_scale >= 0.70 and drop_depth < 0.30)
               else ("(i) MASK LEAKAGE" if drop_depth >= 0.70
                     else "NEITHER RULE MET -- see (D)"))
    print(f"    verdict: {verdict}")

    print("\n(D) does the deviation shrink with N?")
    print(f"{'band':>21} {'exc.kurt':>10} {'E|X|/sd dev':>13}")
    for d, (b, hi, sel) in zip(baseA, bs):
        print(f"{b:>9}-{hi:>11} {d['kurt']:>10.4f} {d['d_ea']:>13.5f}")
    ku = np.array([d['kurt'] for d in baseA])
    Ls = np.array([math.log(math.sqrt(b * hi)) for b, hi, _ in bs])
    sl = np.polyfit(Ls, ku, 1)[0]
    print(f"    trend in excess kurtosis vs log N: {sl:+.5f} per unit")
    print(f"    {'decaying' if sl < 0 else 'NOT decaying'}")

    # ---- (C2) which normaliser? ----------------------------------
    # Increment 240 used Z = C/sqrt(V) with V(N) = Sum_v mu^2(v)
    # Lambda(N-v)^2, the EXACT second moment, and reported excess
    # kurtosis +0.0143. Everything above used Z = C/sqrt(S(N)*N),
    # which is what CONJECTURE_L.md and MEASUREMENTS display. If the
    # "variance mask" of (C) is just the difference between those two
    # normalisers, then repeating (A)-(C) under V must make it vanish
    # WITHOUT any per-cell standardisation. That is a sharp test and
    # it decides whether a new mask has been found or an old one
    # rediscovered.
    print("\n(C2) which normaliser? S(N)*N (displayed in the docs)")
    print("     versus V(N) = Sum_v mu^2 Lambda^2 (exact, used by")
    print("     increment 240)")
    n2 = 1
    while n2 < 2 * (X + 1):
        n2 *= 2
    A2 = np.zeros(n2); A2[: X + 1] = np.abs(mu)
    B2 = np.zeros(n2); B2[: X + 1] = lam ** 2
    V = np.fft.irfft(np.fft.rfft(A2) * np.fft.rfft(B2), n2)[: X + 1]
    del A2, B2
    Zv = C[Ns] / np.sqrt(np.maximum(V[Ns], 1e-300))
    print(f"{'normaliser':>14} {'treatment':>22} {'exc.kurt':>10} "
          f"{'z':>8} {'E|X|/sd dev':>13} {'z':>8}")
    for nlab, ZZ in (("S(N)*N", Z), ("V(N) exact", Zv)):
        for tlab, sc in (("cell means only", False),
                         ("cell mean + scale", True)):
            rs, ks = [], 0
            for b, hi, sel in bs:
                r, k = demask(ZZ[sel], cellkey(Ns[sel], 8), scale=sc)
                rs.append(r if sc else r / r.std())
                ks += k
            dd = report("x", np.concatenate(rs), ks)
            print(f"{nlab:>14} {tlab:>22} {dd['kurt']:>10.4f} "
                  f"{dd['z_ku']:>8.1f} {dd['d_ea']:>13.5f} "
                  f"{dd['z_ea']:>8.1f}")
    rat = (S[Ns] * Ns) / np.maximum(V[Ns], 1e-300)
    print(f"    the two normalisers differ: (S*N)/V has mean "
          f"{float(rat.mean()):.4f}, sd {float(rat.std()):.4f}, "
          f"range [{float(rat.min()):.3f}, {float(rat.max()):.3f}]")
    kk, _ = demask(rat, cellkey(Ns, 8))
    print(f"    and {100*(1 - float(kk.var())/float(rat.var())):.1f}% of "
          f"its variance is explained by the SAME cells -- which is")
    print("    why per-cell standardisation imitated it")

    print("\n(E) robustness: quarter-octave bands (scale mixing within")
    print("    a band is ~0.001 in excess kurtosis, predicted; here it")
    print("    is measured)")
    bq = bands(Ns, lo, X, 2 ** 0.25)
    rs, ks = [], 0
    for b, hi, sel in bq:
        r, k = demask(Z[sel], cellkey(Ns[sel], 8))
        rs.append(r / r.std())
        ks += k
    dd = report("q", np.concatenate(rs), ks)
    print(f"    {len(bq)} bands, excess kurtosis {dd['kurt']:+.4f} "
          f"(z {dd['z_ku']:.1f}), octave value "
          f"{kur_by_depth[8]:+.4f}")
    print("DONE")


if __name__ == "__main__":
    main()
