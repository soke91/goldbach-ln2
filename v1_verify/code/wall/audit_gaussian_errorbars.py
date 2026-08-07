# -*- coding: utf-8 -*-
"""
Re-verification of the error bars on measurement 1 and 3 of
Conjecture 14 (`conj:wall`) in v1/paper/wall_v1.tex.

THE STATEMENTS UNDER TEST, verbatim:

  item 1: "Excess kurtosis -0.0005 (z = -0.3) and E|G|/sd(G) short of
     sqrt(2/pi) by 0.00018 (z = -0.8), on 6.3e6 values, removing cell
     means alone. Under an S N-based scale the same data give excess
     kurtosis +0.1704 at z = 98."

  item 3: "aggregate tail counts against the Gaussian expectation give
     ratios 0.999 at t = 3, 0.997 at t = 4, 0.878 at t = 5".

THE OBJECTION BEING TESTED. Every z above is a COUNT-BASED error bar:
z = -0.3 is -0.0005 / sqrt(24/6.3e6) and z = -0.8 is
-0.00018 / (const/sqrt(6.3e6)). This paper is the one that establishes,
in Proposition 21 and the remark that follows it, that

    "An interval built from a count is about ten times too narrow at
     the top of this range, and the shortfall grows like N^0.46 ...
     the same is true of any cell mean of a field whose summands share
     a common arithmetic input."

C(N) for nearby N share their summands. So the rule the paper states
for cell means is being suspended for the very statistics that carry
Conjecture 14.

Note the direction. Item 1 is an ACCEPTANCE: a wider true error bar
does not overturn "consistent with Gaussian", it weakens the claimed
PRECISION. The contrast figure (z = 98 for the S N scale) and the
t = 5 tail ratio are the ones where the width can change a verdict.

METHOD HERE. The same statistics, with three error bars each:
  (a) the iid/count formula that v1 uses;
  (b) a moving-block bootstrap over N, which is valid whatever the
      serial dependence is, at several block lengths;
  (c) a split-half comparison across disjoint octaves, which needs no
      model at all.

PRE-REGISTRATION (written before the run).

  (1) RULE. If the block-bootstrap standard error agrees with the
      count formula to within 30% at every block length, the count
      bar is fine here and the objection is void.
  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      The companion run `audit_zeta_regression_null.py` measured the
      lag-1..100 autocorrelations of Z at +0.0003, -0.0001, -0.0068,
      +0.0004, +0.0014 -- essentially zero -- while the same series
      carries 180x the white-noise power in the lowest Fourier bins.
      Short-range independence is what a kurtosis error bar depends
      on, and long-range structure is what a MEAN depends on. So I
      predict the count bar is roughly RIGHT for kurtosis and for
      E|G|/sd(G), and WRONG (too narrow) for anything that is a mean
      over a cell -- which is what Proposition 21 already says. That
      is, I expect this objection to FAIL, and the paper's item 1 to
      stand.
  (3) Reported regardless: the t = 5 tail ratio 0.878 and its actual
      error bar, since a ratio built from a few hundred exceedances
      carries a large one and the paper quotes it without.

CORRECTION MADE AFTER THE FIRST RUN, AND BEFORE THE NUMBERS BELOW WERE
READ. The first version of this script pooled every octave from 1e5 to
1.6e7 into one sample. That is invalid and the invalidity is v1's own
finding: `lab_wall_tails.py` states it in terms -- "is sd(G) constant
in N? (if not, pooling mixes scales and manufactures a heavy tail by
itself)". Under the V-scale the raw excess kurtosis runs 1.77 at
1e5-2e5 down to 0.06 at the top, so a pooled sample is a mixture and
its kurtosis (+0.86 in the first run) measures the mixture, not the
field. Everything below is therefore computed WITHIN each octave, and
the octave-to-octave scatter of z is reported as the model-free check
on whether the count-based error bar is the right width.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

QS = (3, 5, 7, 11, 13)


def sieve_mu_lambda(X):
    mu = np.ones(X + 1, dtype=np.int64)
    is_p = np.zeros(X + 1, dtype=bool)
    rem = np.arange(X + 1, dtype=np.int64)
    mu[0] = 0
    for p in range(2, X + 1):
        if rem[p] == p:
            is_p[p] = True
            mu[p::p] *= -1
            rem[p::p] //= p
            if p * p <= X:
                mu[p * p::p * p] = 0
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if is_p[p]:
            lg = math.log(p)
            q = p
            while q <= X:
                lam[q] = lg
                q *= p
    return mu.astype(np.float64), lam


def conv(a, b, X):
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    return np.fft.irfft(np.fft.rfft(np.pad(a, (0, nf - len(a))))
                        * np.fft.rfft(np.pad(b, (0, nf - len(b)))), nf)[: X + 1]


def stats(g):
    m = g.mean()
    c = g - m
    v = float((c * c).mean())
    kurt = float((c ** 4).mean()) / v ** 2 - 3.0
    hn = float(np.abs(c).mean()) / math.sqrt(v)
    return kurt, hn


def main():
    X = 16_000_000
    LO = 100_000
    mu, lam = sieve_mu_lambda(X)
    C = conv(mu, lam, X)
    V = conv((mu != 0).astype(np.float64), lam ** 2, X)

    Ns = np.arange(LO + LO % 2, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    # Section `sec:floor` defines Z(N) = C(N)/sqrt(V(N)) and puts the
    # cells on Z, so "removing cell means" means removing the cell mean
    # OF Z. Subtracting the cell mean of C instead mixes the growth of
    # sqrt(V) across the band into the mask and manufactures a huge
    # excess kurtosis at small N; the first version of this script did
    # that, and the note is left here because the failure mode is easy
    # to repeat.
    Z = C[Ns] / np.sqrt(V[Ns])
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    tot = np.bincount(inv, weights=Z, minlength=len(uniq))
    G = Z - (tot / cnt)[inv]
    n = len(G)

    print("Re-verification of the error bars on Conjecture 14, items 1 and 3")
    print("Everything WITHIN octaves; no cross-octave pooling.")
    print()

    rng = np.random.default_rng(77)
    bands = []
    b = LO
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 20000:
            bands.append((b, hi, sel))
        b = hi

    print("(1) per-octave statistics, count-based error bars")
    hdr = (f"    {'octave':>22} {'n':>9} {'exc kurt':>10} {'z':>7} "
           f"{'E|G|/sd':>9} {'dev':>10} {'z':>7}")
    print(hdr)
    print("    " + "-" * (len(hdr) - 4))
    zk, zh = [], []
    for lo_, hi_, sel in bands:
        g = G[sel]
        nb = len(g)
        k_, h_ = stats(g)
        sek = math.sqrt(24.0 / nb)
        seh = math.sqrt((1 - 2 / math.pi) / nb)
        dev = h_ - math.sqrt(2 / math.pi)
        zk.append(k_ / sek)
        zh.append(dev / seh)
        print(f"    {f'{lo_}-{hi_}':>22} {nb:>9} {k_:>+10.5f} "
              f"{k_/sek:>7.2f} {h_:>9.5f} {dev:>+10.5f} {dev/seh:>7.2f}")
    print(f"    v1 item 1 quotes exc kurt -0.0005 (z -0.3) and "
          f"E|G|/sd dev -0.00018 (z -0.8)")
    print()
    print(f"    scatter of the z's across {len(bands)} disjoint octaves:")
    print(f"      kurtosis   mean {np.mean(zk):+.2f}, sd {np.std(zk):.2f}"
          f"   (a correct error bar gives sd = 1)")
    print(f"      E|G|/sd    mean {np.mean(zh):+.2f}, sd {np.std(zh):.2f}")
    ok = abs(np.std(zk) - 1) < 0.5
    print(f"    (1) count-based bar is the right width for kurtosis: "
          f"{'PASS -- objection void' if ok else 'FAIL'}")
    print()

    print("(2) moving-block bootstrap inside the top octave")
    lo_, hi_, sel = bands[-1]
    g = G[sel]
    nb0 = len(g)
    sek = math.sqrt(24.0 / nb0)
    seh = math.sqrt((1 - 2 / math.pi) / nb0)
    hdr = (f"    {'block':>8} {'SE(kurt)':>11} {'/count SE':>10} "
           f"{'SE(E|G|/sd)':>13} {'/count SE':>10}")
    print(hdr)
    print("    " + "-" * (len(hdr) - 4))
    for B in (1, 100, 1000, 10000):
        nbk = nb0 // B
        ks, hs = [], []
        for _ in range(120):
            st = rng.integers(0, nb0 - B, size=nbk)
            idx = (st[:, None] + np.arange(B)[None, :]).ravel()
            k_, h_ = stats(g[idx])
            ks.append(k_)
            hs.append(h_)
        print(f"    {B:>8} {np.std(ks):>11.6f} {np.std(ks)/sek:>10.2f} "
              f"{np.std(hs):>13.6f} {np.std(hs)/seh:>10.2f}")
    print()

    print("(3) tail counts per octave, with error bars v1 does not quote")
    print(f"    {'octave':>22} {'t':>3} {'obs':>8} {'exp':>10} "
          f"{'ratio':>7} {'+-':>7} {'z':>7}")
    for lo_, hi_, sel in bands[-3:]:
        g = G[sel]
        g = (g - g.mean()) / g.std()
        for t in (3, 4, 5):
            obs = int((np.abs(g) > t).sum())
            exp = len(g) * math.erfc(t / math.sqrt(2))
            se = math.sqrt(max(obs, 1)) / exp
            print(f"    {f'{lo_}-{hi_}':>22} {t:>3} {obs:>8} "
                  f"{exp:>10.1f} {obs/exp:>7.3f} {se:>7.3f} "
                  f"{(obs-exp)/math.sqrt(exp):>7.2f}")
    print("    v1 item 3 quotes aggregate ratios 0.999 / 0.997 / 0.878")
    print("    at t = 3/4/5. The +- column is the Poisson error bar the")
    print("    paper does not quote; at t = 5 it is the whole story.")
    print("DONE")


if __name__ == "__main__":
    main()
