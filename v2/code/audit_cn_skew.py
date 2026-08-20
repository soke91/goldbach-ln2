# -*- coding: utf-8 -*-
r"""
A statistic whose control shrinks: the skewness of G(N)

WHAT IS AT STAKE

rem:cnkurtlimit closed the kurtosis route by finding its floor.  The
excess kurtosis of G(N) = C(N)/sqrt(V(N)) is real and outside its coin
ensemble from b = 17 to b = 23, follows a power law well enough to
forecast two octaves out of sample, and then meets the coin at
b = 25 -- because the coin's draw-to-draw spread in excess kurtosis is
0.018 at b = 23 and 0.018 at b = 25, set by the field's correlation
across N rather than by sample size, while the real signal decays like
N^-0.73.  **Pushing N cannot separate them, because only one of the
two is falling.**

What OPEN.md item 1 asks for next is therefore not more N.  It is a
statistic whose control shrinks with N, so that the separation widens
instead of closing.  The named candidate is the tail asymmetry:
rem:cnlaw printed, at the 0.1 and 99.9 per cent quantiles of the real
arm, -3.1965 and +2.8074 against a normal's -+3.0902 -- heavier on the
left and lighter on the right, which is a third moment and not a
fourth.

**And there is a structural reason to expect its control to behave
differently.**  Replacing eps by -eps sends C_coin to -C_coin and so
G_coin to -G_coin.  The skewness is odd, so the coin's skewness
distribution over draws is exactly symmetric about zero: its mean is
zero by construction, not by measurement.  Nothing of the kind holds
for the kurtosis, whose coin mean had to be measured and turned out
non-zero and drifting.  Whether the *spread* of that symmetric
distribution shrinks with N is the question this script exists to
answer, and it decides whether this branch has anywhere left to go.

BACKS: Remark {#rem:cnskew} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  THE GATE.  Band b = 23's excess kurtosis reproduces
      rem:cnkurt's POINT cnkurt_23 to three decimals, tying this run
      to that one.
  S2  The real skewness is negative at every one of the seven
      octaves and lies outside its coin ensemble at |z| > 3.
  S3  **The control shrinks.**  The coin's draw-to-draw standard
      deviation of the skewness falls with N: regressing its log on
      log2 N over the seven octaves gives a negative slope resolved
      at |t| > 2.  This is the property the kurtosis did not have.
  S4  And so the separation widens: the real point's z against the
      coin ensemble grows with N, a positive slope resolved at
      |t| > 2.  The kurtosis went the other way, 33.4 down to 5.0.
  S5  The control's centre is where the symmetry says: the coin's
      mean skewness is within two of its own standard errors of zero
      at every band.

REFUTATION RULE (fixed before the run)

  S1  REFUTED outside three decimals; nothing below is reported.
  S2  REFUTED if any octave has positive skewness or fails |z| > 3.
      A sign that flips across octaves would mean the tail asymmetry
      is not a property of C but of the band, and the quantiles
      rem:cnlaw printed were one band's.
  S3  **REFUTED if the slope is positive or unresolved.**  Unresolved
      is the outcome to expect if the coin's skewness spread behaves
      like its kurtosis spread did -- flat -- and then this branch
      has no route left that more computation can open, which is a
      real finding and the one that would close item 1's remaining
      line.
  S4  REFUTED if the z falls with N or the slope is unresolved.  S3
      and S4 can disagree: a shrinking control with a faster
      shrinking signal still closes.  If S3 holds and S4 does not,
      the honest reading is that the control shrinks and the signal
      shrinks faster, which is the same wall in a different place.
  S5  REFUTED if any band's coin mean is two standard errors from
      zero.  That would contradict an exact symmetry and would mean
      the ensemble is not what it is taken to be -- an instrument
      failure rather than a finding.

  WHAT THIS CANNOT DO.  Seven octaves and thirty-two draws.  If S3
  and S4 hold, this says the route is open, not how far; the reach
  is a separate measurement and none is attempted here.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_cn_skew.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

BLO, BHI = 17, 23
NMAX = 1 << (BHI + 1)
DRAWS = 32
SEED = 20260823
DEC = 3
ZKEEP = 3.0
ZCAP = 2.0


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def skew(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    return float((d ** 3).mean() / v ** 1.5)


def excess_kurt(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    return float((d ** 4).mean() / v ** 2 - 3.0)


def ols(x, y):
    A = np.column_stack([np.ones(len(y)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    s2 = float((r ** 2).sum()) / (len(y) - 2)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, np.sqrt(np.diag(cov))


def published_kurt():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT cnkurt_%d ([-\d.]+)\s*$" % BHI, src, re.M)
    if not m:
        raise SystemExit("no POINT cnkurt_%d in %s" % (BHI, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the skewness of G(N) = C(N)/sqrt(V(N)) over even N in",
    "           each of seven octaves, against an ensemble of %d sign"
    % DRAWS,
    "           patterns recomputed over the whole ladder; and how",
    "           that ensemble's own spread, and the real point's z",
    "           against it, move with N.",
    "FIELD: even N in (2^b, 2^(b+1)] for b = %d..%d, so N runs to"
    % (BLO, BHI),
    "       %d, with Lambda and mu sieved once to that top and both"
    % NMAX,
    "       convolutions taken by FFT over the whole range, then",
    "       sliced per band. The gate re-measures the excess kurtosis",
    "       of band b = %d against rem:cnkurt's own value." % BHI,
    "SEED: the coin draws from numpy default_rng at seed %d; without"
    % SEED,
    "      it the file does not reproduce its own control.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = published_kurt()
    say("READ audit_cn_kurt_drift.txt POINT cnkurt_%d %.5f" % (BHI, pub))
    say("  the band this run has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_skew %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d fresh sign patterns over the whole ladder" % DRAWS)

    say("sieving to %d" % NMAX)
    lam, mu = sieves(NMAX)
    mu2 = (mu.astype(np.float64)) ** 2
    L = 1 << (2 * NMAX).bit_length()
    say("  transform length %d" % L)
    FL = np.fft.rfft(lam, L)
    C = np.fft.irfft(FL * np.fft.rfft(mu.astype(np.float64), L),
                     L)[: NMAX + 1]
    V = np.fft.irfft(np.fft.rfft(lam ** 2, L) * np.fft.rfft(mu2, L),
                     L)[: NMAX + 1]
    del mu, lam
    rootV = np.sqrt(np.maximum(V, 1e-300))
    del V

    def band(b):
        return np.arange((1 << b) + 2, (1 << (b + 1)) + 1, 2,
                         dtype=np.int64)

    bs = list(range(BLO, BHI + 1))
    Gs = {b: C[band(b)] / rootV[band(b)] for b in bs}
    gate_k = excess_kurt(Gs[BHI])
    real = {b: skew(Gs[b]) for b in bs}
    del C, Gs

    # -------------------------------------------------------------- S1
    say()
    say("S1  does this run reproduce the band it shares?")
    say("  b = %d excess kurtosis: here %.5f against its %.5f"
        % (BHI, gate_k, pub))
    s1 = abs(round(gate_k, DEC) - round(pub, DEC)) < 10.0 ** (-DEC) / 2
    say("  S1 %s   (cap: %d decimals)"
        % ("hold" if s1 else "REFUTED", DEC))
    if not s1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    coin = {b: [] for b in bs}
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L), L)[: NMAX + 1]
        for b in bs:
            coin[b].append(skew(Cc[band(b)] / rootV[band(b)]))
        del Cc
    coin = {b: np.array(v) for b, v in coin.items()}

    say()
    say("       b      even N       real skew    coin mean    coin sd"
        "        z")
    zs, sds = [], []
    for b in bs:
        cm, cs = float(coin[b].mean()), float(coin[b].std(ddof=1))
        z = (real[b] - cm) / cs
        zs.append(z)
        sds.append(cs)
        say("      %2d  %10d  %+12.5f  %+10.5f  %9.5f  %7.2f"
            % (b, len(band(b)), real[b], cm, cs, z))
        say("POINT cnskew_%d %.5f" % (b, real[b]))
        say("POINT cnskewsd_%d %.5f" % (b, cs))
    say("SCALES %d" % len(bs))
    zs = np.array(zs)
    sds = np.array(sds)
    x = np.array([b + 0.5 for b in bs])

    # -------------------------------------------------------------- S2
    say()
    say("S2  is the asymmetry real and negative at every octave?")
    neg = all(real[b] < 0 for b in bs)
    out = all(abs(z) > ZKEEP for z in zs)
    say("  every octave negative: %s; smallest |z|: %.2f"
        % ("yes" if neg else "NO", float(np.abs(zs).min())))
    say("TSTAT cnskew_min_z %.2f" % zs[np.argmin(np.abs(zs))])
    say("SPREAD cnskew_min_z %.5f" % sds[int(np.argmin(np.abs(zs)))])
    s2 = neg and out
    say("  S2 %s   (cap: negative everywhere and |z| > %.0f)"
        % ("hold" if s2 else "REFUTED", ZKEEP))

    # -------------------------------------------------------------- S3
    say()
    say("S3  does the control shrink with N?")
    c, se = ols(x, np.log(sds))
    t = c[1] / se[1]
    say("  log(coin sd) on log2 N: slope %+.5f +- %.5f, t = %.2f"
        % (c[1], se[1], t))
    say("  so the coin's skew spread goes like N^%+.4f, against the "
        "kurtosis" % (c[1] / math.log(2)))
    say("  spread which was flat")
    say("TSTAT cnskew_coinsd_slope %.2f" % t)
    say("SPREAD cnskew_coinsd_slope %.5f" % (x.max() - x.min()))
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN cnskew_coinsd_slope")
    s3 = c[1] < 0 and abs(t) > ZCAP
    say("POINT cnskew_coin_power %.5f" % (c[1] / math.log(2)))
    say("  S3 %s   (cap: negative and |t| > %.0f)"
        % ("hold" if s3 else "REFUTED", ZCAP))

    # -------------------------------------------------------------- S4
    say()
    say("S4  does the separation widen?")
    cz, sez = ols(x, np.abs(zs))
    tz = cz[1] / sez[1]
    say("  |z| on log2 N: slope %+.5f +- %.5f, t = %.2f"
        % (cz[1], sez[1], tz))
    say("  |z| runs %.2f at b = %d to %.2f at b = %d"
        % (abs(zs[0]), bs[0], abs(zs[-1]), bs[-1]))
    say("TSTAT cnskew_z_slope %.2f" % tz)
    say("SPREAD cnskew_z_slope %.5f" % (x.max() - x.min()))
    if abs(tz) < 2.0:
        say("UNRESOLVED SIGN cnskew_z_slope")
    s4 = cz[1] > 0 and abs(tz) > ZCAP
    say("  S4 %s   (cap: positive and |t| > %.0f)"
        % ("hold" if s4 else "REFUTED", ZCAP))
    cr, ser = ols(x, np.log(np.abs([real[b] for b in bs])))
    say("  for comparison the signal itself goes like N^%+.4f"
        % (cr[1] / math.log(2)))
    say("POINT cnskew_signal_power %.5f" % (cr[1] / math.log(2)))

    # -------------------------------------------------------------- S5
    say()
    say("S5  is the control centred where the symmetry says?")
    worst, wb = 0.0, None
    for b in bs:
        r = abs(float(coin[b].mean())) / (float(coin[b].std(ddof=1))
                                          / math.sqrt(DRAWS))
        if r > worst:
            worst, wb = r, b
    say("  furthest from zero at b = %s, at %.2f of its own standard "
        "error" % (wb, worst))
    s5 = worst <= ZCAP
    say("  S5 %s   (cap: two standard errors)"
        % ("hold" if s5 else "REFUTED"))

    say()
    say("=" * 70)
    say("S1 %s  S2 %s  S3 %s  S4 %s  S5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (s1, s2, s3, s4, s5)))
    say()
    if s3 and s4:
        say("this branch is not closed. The skewness has what the "
            "kurtosis did")
        say("not: a control that shrinks with N, so the separation "
            "widens as the")
        say("field grows instead of closing at a fixed scale. How far "
            "it goes is")
        say("a different measurement and none is made here.")
    elif s3 and not s4:
        say("the control shrinks and the signal shrinks faster, so "
            "the separation")
        say("still closes -- the same wall as the kurtosis, in a "
            "different place.")
    else:
        say("the control does not shrink, so the skewness is the "
            "kurtosis again")
        say("and this branch has no route that more computation can "
            "open.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
