# -*- coding: utf-8 -*-
r"""
Where in the sum does C(N)'s non-Gaussianity live?

WHAT IS AT STAKE

Four routes are now closed -- rem:cnkurtlimit, rem:cnskew,
rem:cnshift, rem:cndilation -- and every one of them closed the same
way: independent signs on mu's support reproduce what the field does.
But one thing in that branch is not closed and is not small.  At
b = 20 the excess kurtosis of G(N) = C(N)/sqrt(V(N)) is +0.41549
against a coin ensemble at -0.01824 with spread 0.02628, a z of 16.5,
and at b = 17 the z is 33.4.  **The finite-N separation is real and
has no mechanism.**  Nothing measured so far says which part of the
sum produces it.

C(N) = sum over v < N of mu(v) Lambda(N-v), and with N-v = p this is
sum over primes p < N of log(p) mu(N-p).  Splitting the range of v
dyadically splits the sum by how far the prime sits below N: small v
means p near N, few terms, each of weight about log N; large v means
small p, many terms of small weight.  If the non-Gaussianity is
carried by one part, that is a mechanism; if it is spread evenly, the
absence of one is itself worth knowing.

The coin control is essential here and not decoration.  A piece with
few terms is non-Gaussian for any signs at all, so raw per-piece
kurtosis says nothing; what the coin cannot make is the *excess* over
its own ensemble at the same piece, with the same term count and the
same weights.

BACKS: Remark {#rem:cnkurtwhere} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  THE GATE.  Two checks.  The pieces sum to the whole: the
      largest absolute difference between C(N) and the sum of the
      piece fields, over the band, is below 1e-6.  And the whole
      band's excess kurtosis reproduces rem:cnkurt's POINT cnkurt_20
      to three decimals.
  W2  **It is concentrated.**  The largest per-piece z against that
      piece's own coin ensemble exceeds the whole field's z at this
      band, which rem:cnkurt's numbers put at 16.5.  A quantity
      spread evenly over pieces would give every piece a smaller z
      than the whole.
  W3  And it is carried by the small-v end: the piece with the
      largest z has its v-range in the lower half of the dyadic
      ladder used here.
  W4  The control is doing work: at least one piece has a coin
      excess kurtosis above 0.5, so that reading the per-piece
      kurtosis without a control would point somewhere the z does
      not.

REFUTATION RULE (fixed before the run)

  W1  REFUTED on either check; nothing below is reported.
  W2  REFUTED if every piece's z is below the whole field's.  Then
      the non-Gaussianity is spread across the sum rather than made
      somewhere, and there is no mechanism of this kind to find --
      which is a real answer to the question and closes it rather
      than leaving it open.
  W3  REFUTED if the largest z sits in the upper half.  W3 can fail
      while W2 holds; the concentration would then be at small primes
      rather than at primes near N, which is the opposite mechanism
      and points at a different part of the arithmetic.  W3 is not
      read at all if W2 fails.
  W4  REFUTED if no piece's coin kurtosis reaches 0.5.  That would
      mean the pieces are all large enough for the central limit
      theorem and the control was not needed -- harmless, and worth
      knowing because it would make the raw per-piece table
      interpretable on its own.

  Pieces where V_j(N) vanishes for some N are excluded at those N and
  the excluded fraction is printed per piece; a piece excluded at more
  than half its N is dropped and named.

  WHAT THIS CANNOT DO.  One band.  A mechanism located here is
  located at N ~ 2e6 and nothing here says it stays put as N grows;
  the drift of the location is a separate measurement.
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
OUT = os.path.join(ROOT, "results", "audit_cn_kurt_where.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
JLO, JHI = 8, BAND          # dyadic v-pieces (2^j, 2^(j+1)]
DRAWS = 32
SEED = 20260826
DEC = 3
TOLSUM = 1e-6
WHOLEZ = 16.5
COINBIG = 0.5
DROPFRAC = 0.5


assert (1 << (JHI + 1)) <= NMAX, "the top piece must fit under NMAX"


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


def excess_kurt(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    if v <= 0:
        return float("nan")
    return float((d ** 4).mean() / v ** 2 - 3.0)


def published_kurt():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT cnkurt_%d ([-\d.]+)\s*$" % BAND, src, re.M)
    if not m:
        raise SystemExit("no POINT cnkurt_%d in %s" % (BAND, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the excess kurtosis of G_j(N) = C_j(N)/sqrt(V_j(N))",
    "           over even N in one band, where C_j restricts the sum",
    "           over v to a dyadic piece, against an ensemble of %d"
    % DRAWS,
    "           sign patterns restricted to the same piece; and where",
    "           over the pieces the whole field's non-Gaussianity",
    "           sits.",
    "FIELD: even N in (2^%d, 2^%d], with v-pieces (2^j, 2^(j+1)] for"
    % (BAND, BAND + 1),
    "       j = %d..%d and one lumped piece v <= 2^%d; Lambda and mu"
    % (JLO, JHI, JLO),
    "       sieved once to %d and every piece's two convolutions" % NMAX,
    "       taken by FFT over the whole range. N with V_j(N) = 0 are",
    "       excluded at that piece and the fraction is printed. The",
    "       gate re-measures the whole band against rem:cnkurt.",
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
    say("READ audit_cn_kurt_drift.txt POINT cnkurt_%d %.5f"
        % (BAND, pub))
    say("  the band this run has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_kurt_where %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d fresh sign patterns per piece; the whole field's z at "
        "this band" % DRAWS)
    say("  is %.1f, from rem:cnkurt's own numbers" % WHOLEZ)

    say("sieving to %d, which is 2^%d" % (NMAX, BAND + 1))
    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    mu2 = muf ** 2
    del mu
    L = 1 << (2 * NMAX).bit_length()
    say("  transform length %d" % L)
    FL = np.fft.rfft(lam, L)
    FL2 = np.fft.rfft(lam ** 2, L)
    del lam

    Ns = np.arange((1 << BAND) + 2, (1 << (BAND + 1)) + 1, 2,
                   dtype=np.int64)

    pieces = [(0, 1 << JLO)] + [(1 << j, 1 << (j + 1))
                                for j in range(JLO, JHI + 1)]

    def restrict(arr, lo, hi):
        out = np.zeros_like(arr)
        out[lo + 1: hi + 1] = arr[lo + 1: hi + 1]
        return out

    def conv(f, g):
        return np.fft.irfft(f * np.fft.rfft(g, L), L)[: NMAX + 1]

    Ctot = conv(FL, muf)
    Vtot = conv(FL2, mu2)
    whole = Ctot[Ns] / np.sqrt(np.maximum(Vtot[Ns], 1e-300))
    gate_k = excess_kurt(whole)

    Csum = np.zeros(len(Ns))
    Cj, Vj = [], []
    for lo, hi in pieces:
        c = conv(FL, restrict(muf, lo, hi))[Ns]
        v = conv(FL2, restrict(mu2, lo, hi))[Ns]
        Csum += c
        Cj.append(c)
        Vj.append(v)
    worst = float(np.abs(Csum - Ctot[Ns]).max())

    # -------------------------------------------------------------- W1
    say()
    say("W1  do the pieces sum to the whole, and does the band "
        "reproduce?")
    say("  largest |sum of pieces - C| over the band %.3e" % worst)
    say("  whole-band excess kurtosis here %.5f against its %.5f"
        % (gate_k, pub))
    w1 = (worst <= TOLSUM
          and abs(round(gate_k, DEC) - round(pub, DEC))
          < 10.0 ** (-DEC) / 2)
    say("  W1 %s   (cap: %.0e and %d decimals)"
        % ("hold" if w1 else "REFUTED", TOLSUM, DEC))
    if not w1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)
    del Ctot, Csum

    keep = []
    for i, (lo, hi) in enumerate(pieces):
        bad = float((Vj[i] <= 0).mean())
        if bad > DROPFRAC:
            say("  piece v in (%d, %d] dropped: V_j = 0 at %.4f of N"
                % (lo, hi, bad))
            continue
        keep.append(i)

    rng = np.random.default_rng(SEED)
    coin = {i: [] for i in keep}
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        w = eps * mu2
        for i in keep:
            lo, hi = pieces[i]
            c = conv(FL, restrict(w, lo, hi))[Ns]
            m = Vj[i] > 0
            coin[i].append(excess_kurt(c[m] / np.sqrt(Vj[i][m])))
        del w

    say()
    say("      v-piece            terms      empty V     real kurt"
        "     coin mean    coin sd       z")
    zs, reals, coinmeans = {}, {}, {}
    for i in keep:
        lo, hi = pieces[i]
        m = Vj[i] > 0
        r = excess_kurt(Cj[i][m] / np.sqrt(Vj[i][m]))
        arr = np.array(coin[i])
        cm, cs = float(arr.mean()), float(arr.std(ddof=1))
        z = (r - cm) / cs
        zs[i], reals[i], coinmeans[i] = z, r, cm
        say("  (%9d, %9d]  %9d  %11.5f  %+11.5f  %+11.5f  %9.5f  %7.2f"
            % (lo, hi, int(mu2[lo + 1: hi + 1].sum()),
               1.0 - float(m.mean()), r, cm, cs, z))
        say("POINT kurtwhere_%d %.5f" % (hi, r))
        say("POINT kurtwherez_%d %.5f" % (hi, z))
    say("SCALES 1")

    top = max(keep, key=lambda i: abs(zs[i]))
    say()
    say("  largest |z| %.2f at v in (%d, %d]"
        % (abs(zs[top]), pieces[top][0], pieces[top][1]))
    say("POINT kurtwhere_maxz %.5f" % abs(zs[top]))

    # -------------------------------------------------------------- W2
    say()
    say("W2  is it concentrated?")
    w2 = abs(zs[top]) > WHOLEZ
    say("  largest per-piece |z| %.2f against the whole field's %.1f"
        % (abs(zs[top]), WHOLEZ))
    say("TSTAT kurtwhere_maxz %.2f" % abs(zs[top]))
    say("SPREAD kurtwhere_maxz %.5f"
        % float(np.array(coin[top]).std(ddof=1)))
    say("  W2 %s   (cap: above the whole field's z)"
        % ("hold" if w2 else "REFUTED"))

    # -------------------------------------------------------------- W3
    say()
    say("W3  is it at the small-v end?")
    half = keep[len(keep) // 2]
    w3 = keep.index(top) < len(keep) / 2.0
    say("  the carrying piece is %d of %d up the ladder, the midpoint "
        "being v = %d" % (keep.index(top) + 1, len(keep),
                          pieces[half][1]))
    say("  W3 %s   (cap: lower half)" % ("hold" if w3 else "REFUTED"))
    if not w2:
        say("  W2 failed, so W3 is recorded and not read")

    # -------------------------------------------------------------- W4
    say()
    say("W4  is the control doing work?")
    big = [i for i in keep if coinmeans[i] > COINBIG]
    say("  pieces whose coin excess kurtosis exceeds %.1f: %d of %d"
        % (COINBIG, len(big), len(keep)))
    if big:
        i = max(big, key=lambda i: coinmeans[i])
        say("  the largest is %+.5f at v in (%d, %d], where the real "
            "arm reads %+.5f"
            % (coinmeans[i], pieces[i][0], pieces[i][1], reals[i]))
    w4 = len(big) > 0
    say("  W4 %s   (cap: at least one piece)"
        % ("hold" if w4 else "REFUTED"))

    say()
    say("=" * 70)
    say("W1 %s  W2 %s  W3 %s  W4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (w1, w2, w3, w4)))
    say()
    if w2 and w3:
        say("the non-Gaussianity is made in one part of the sum, and "
            "it is the")
        say("part where the prime sits near N -- few terms, each of "
            "weight about")
        say("log N. That is a mechanism, and it is the first thing "
            "this branch")
        say("has said about how C(N) departs from a coin rather than "
            "whether it")
        say("does.")
    elif w2 and not w3:
        say("the non-Gaussianity is made in one part of the sum and "
            "it is the")
        say("small-prime end, not the end near N. That is the "
            "opposite mechanism")
        say("to the one expected and points at a different part of "
            "the arithmetic.")
    else:
        say("the non-Gaussianity is spread across the sum rather than "
            "made")
        say("somewhere. There is no mechanism of this kind to find, "
            "which answers")
        say("the question rather than leaving it open.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
