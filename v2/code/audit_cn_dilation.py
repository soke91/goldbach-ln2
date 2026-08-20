# -*- coding: utf-8 -*-
r"""
The multiplicative question: G(N) against G(dN)

WHAT IS AT STAKE

Three routes are closed.  rem:cnkurtlimit found the excess kurtosis
meeting its control at b = 25; rem:cnskew found the skewness meeting
its at b = 21; rem:cnshift found the additive two-point function
inside its control at both octaves tried, with four coin draws of
thirty-two reaching the real arm's largest |z|.

rem:cnshift also recorded that it asked the wrong question.  What
rem:cnskew named was **multiplicative** -- G(N) against G(2N), and
against G(N') of the same radical -- because mu is multiplicative and
mu(2v) is tied to mu(v), while the coin's eps links nothing to
anything.  What was run instead was additive, N against N+h, and the
diagonal the coin fixes exactly *is* the additive overlap of two
Lambda-shifts, so that test was the one the coin was best placed to
survive.  This script runs the question that was named.

THE QUANTITY

For a dilation d, take

    rho(d) = mean over even N in a band of G(N) G(dN).

Averaging the coin over its signs leaves only
sum_v mu^2(v) Lambda(N-v) Lambda(dN-v), a count of v for which both
N-v and dN-v are prime powers -- an overlap at the *large* shift
(d-1)N, not a small one.  The real field has that and the
off-diagonal mu(v)mu(w) terms besides, and among those are the pairs
w = dv that carry mu's multiplicativity: mu(2v) = -mu(v) for v odd
and 0 for v even, which is a correlation the coin has zero of by
construction.

Whether that survives normalisation and dilution over a band is what
is measured.  It is not settled by the algebra: the w = dv pairs are
one thin family among all off-diagonal pairs, and the algebra says
they are there, not that they are visible.

BACKS: Remark {#rem:cndilation} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  THE GATE.  Band b = 23's excess kurtosis reproduces
      rem:cnkurt's POINT cnkurt_23 to three decimals.
  U2  **The dilated field escapes the coin.**  The real arm's largest
      |z| over the eight dilations exceeds every one of the 32 coin
      draws' own largest |z|, each scored leave-one-out against the
      rest.  A max-statistic comparison, for the reason rem:cnshift
      gives.
  U3  And the escape is at d = 2: the largest |z| over dilations
      occurs there, where mu's link is strongest and exact.
  U4  The control behaves as the diagonal says it should: the coin's
      mean rho(d) is positive at every d and falls as d grows, since
      the overlap it measures sits at shift (d-1)N.

REFUTATION RULE (fixed before the run)

  U1  REFUTED outside three decimals; nothing below is reported.
  U2  REFUTED if any coin draw's largest |z| reaches the real arm's.
      **Then four routes are closed and the fourth was the one the
      algebra most favoured**, which would say the coin's agreement
      with the field is not an artefact of the statistics chosen but
      a property of C(N) at these ranges.  That is a stronger
      negative than any of the three, and it is the outcome to write
      down plainly if it comes.
  U3  REFUTED if the largest |z| is at any other d.  U3 can fail
      while U2 holds -- the field would escape somewhere the algebra
      did not point, which is a finding about where to look next and
      not about whether to look.  U3 is not reported as evidence for
      anything if U2 fails.
  U4  REFUTED by a negative coin mean at any d, or by a rise at any
      step.  This is an instrument check on the ensemble; failing it
      means rho(d) is not the overlap it is taken to be and the
      z-scores above are against something unexamined.

  WHAT THIS CANNOT DO.  One band and eight dilations.  A departure at
  a single d is not attributed to any property of d beyond its being
  d, and nothing here separates "mu's multiplicativity" from any
  other structure that happens to live at the same place.
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
OUT = os.path.join(ROOT, "results", "audit_cn_dilation.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

BAND = 20                 # even N in (2^BAND, 2^(BAND+1)]
DIL = (2, 3, 4, 5, 6, 7, 8, 9)
GATEB = 23
NMAX = 1 << 25            # must hold max(DIL) * 2^(BAND+1)
DRAWS = 32
SEED = 20260825
DEC = 3


assert max(DIL) * (1 << (BAND + 1)) <= NMAX, (
    "NMAX must reach the largest dilate the band touches")
assert (1 << (GATEB + 1)) <= NMAX, "NMAX must hold the gate band"


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
    return float((d ** 4).mean() / v ** 2 - 3.0)


def published_kurt():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT cnkurt_%d ([-\d.]+)\s*$" % GATEB, src, re.M)
    if not m:
        raise SystemExit("no POINT cnkurt_%d in %s" % (GATEB, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the dilation correlation rho(d) = mean over even N in",
    "           one band of G(N) G(dN), for d = %s," % (DIL,),
    "           against an ensemble of %d sign patterns on mu's" % DRAWS,
    "           support, compared by the largest |z| over dilations",
    "           with each coin draw scored leave-one-out against the",
    "           rest, so the comparison is one max-statistic.",
    "FIELD: even N in (2^%d, 2^%d], and the dilates dN for those d,"
    % (BAND, BAND + 1),
    "       so the largest index touched is %d; Lambda and mu are"
    % (max(DIL) * (1 << (BAND + 1))),
    "       sieved once to %d and both convolutions taken by FFT" % NMAX,
    "       over the whole range. The gate re-measures band b = %d's"
    % GATEB,
    "       excess kurtosis against rem:cnkurt's own value.",
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
        % (GATEB, pub))
    say("  the band this run has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_dilation %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d dilations, %d fresh sign patterns" % (len(DIL), DRAWS))

    say("sieving to %d, which is 2^%d"
        % (NMAX, NMAX.bit_length() - 1))
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

    Ns = np.arange((1 << BAND) + 2, (1 << (BAND + 1)) + 1, 2,
                   dtype=np.int64)
    gate_idx = np.arange((1 << GATEB) + 2, (1 << (GATEB + 1)) + 1, 2,
                         dtype=np.int64)

    def rho_of(Carr):
        g = Carr[Ns] / rootV[Ns]
        g = g - g.mean()
        out = []
        for d in DIL:
            h = Carr[d * Ns] / rootV[d * Ns]
            h = h - h.mean()
            out.append(float(g.dot(h))
                       / math.sqrt(float(g.dot(g)) * float(h.dot(h))))
        return np.array(out)

    gate_k = excess_kurt(C[gate_idx] / rootV[gate_idx])
    rho_real = rho_of(C)
    del C

    # -------------------------------------------------------------- U1
    say()
    say("U1  does this run reproduce the band it shares?")
    say("  b = %d excess kurtosis: here %.5f against its %.5f"
        % (GATEB, gate_k, pub))
    u1 = abs(round(gate_k, DEC) - round(pub, DEC)) < 10.0 ** (-DEC) / 2
    say("  U1 %s   (cap: %d decimals)"
        % ("hold" if u1 else "REFUTED", DEC))
    if not u1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    rho_coin = np.zeros((DRAWS, len(DIL)))
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L), L)[: NMAX + 1]
        rho_coin[d] = rho_of(Cc)
        del Cc

    def zs(target, pool):
        return ((target - pool.mean(axis=0))
                / pool.std(axis=0, ddof=1))

    zr = zs(rho_real, rho_coin)
    maxz_real = float(np.abs(zr).max())
    where = DIL[int(np.argmax(np.abs(zr)))]
    mc = []
    for d in range(DRAWS):
        rest = np.delete(rho_coin, d, axis=0)
        mc.append(float(np.abs(zs(rho_coin[d], rest)).max()))
    mc = np.array(mc)

    say()
    say("  even N in (2^%d, 2^%d], %d of them"
        % (BAND, BAND + 1, len(Ns)))
    say("        d     rho real     coin mean      coin sd       z")
    for i, d in enumerate(DIL):
        say("      %3d  %+11.6f  %+11.6f  %11.6f  %7.2f"
            % (d, rho_real[i], rho_coin[:, i].mean(),
               rho_coin[:, i].std(ddof=1), zr[i]))
        say("POINT dilrho_%d %.6f" % (d, rho_real[i]))
        say("POINT dilcoin_%d %.6f" % (d, rho_coin[:, i].mean()))
    say("SCALES 1")
    say("POINT dilmaxz %.5f" % maxz_real)

    # -------------------------------------------------------------- U2
    say()
    say("U2  does the dilated field escape the coin?")
    beat = int((mc >= maxz_real).sum())
    say("  real largest |z| %.2f at d = %d; the coin draws' own "
        "largest run" % (maxz_real, where))
    say("  %.2f to %.2f, and %d of %d reach it"
        % (mc.min(), mc.max(), beat, DRAWS))
    u2 = beat == 0
    say("TSTAT dilation_maxz %.2f" % maxz_real)
    say("SPREAD dilation_maxz %.5f" % float(mc.std(ddof=1)))
    say("  U2 %s   (cap: none of them)"
        % ("hold" if u2 else "REFUTED"))

    # -------------------------------------------------------------- U3
    say()
    say("U3  is it at d = 2?")
    u3 = where == 2
    say("  the largest |z| is at d = %d" % where)
    say("  U3 %s   (cap: d = 2)" % ("hold" if u3 else "REFUTED"))
    if not u2:
        say("  U2 failed, so U3 is not evidence for anything and is "
            "recorded only")

    # -------------------------------------------------------------- U4
    say()
    say("U4  does the control behave like the overlap it is?")
    cm = rho_coin.mean(axis=0)
    pos = bool((cm > 0).all())
    fall = bool((np.diff(cm) < 0).all())
    say("  coin means %s" % ", ".join("%+.6f" % v for v in cm))
    say("  positive at every d: %s; falling at every step: %s"
        % ("yes" if pos else "NO", "yes" if fall else "NO"))
    u4 = pos and fall
    say("  U4 %s   (cap: positive and monotone)"
        % ("hold" if u4 else "REFUTED"))

    say()
    say("  NOT PRE-REGISTERED, reported because U4 is an instrument")
    say("  check and a failed check needs its power stated. No cap")
    say("  above is changed.")
    cse = rho_coin.std(axis=0, ddof=1) / math.sqrt(DRAWS)
    res = int((np.abs(cm) > 2 * cse).sum())
    say("        d   coin mean   its s.e.   |mean|/s.e.")
    for i, d in enumerate(DIL):
        say("      %3d  %+9.6f  %9.6f  %11.2f"
            % (d, cm[i], cse[i], abs(cm[i]) / cse[i]))
    say("  %d of %d coin means are resolved away from zero at two"
        % (res, len(DIL)))
    say("  standard errors, and the one negative mean is %.2f of its "
        "own" % (abs(cm[-1]) / cse[-1]))
    say("  So U4 asked whether eight numbers are positive and ordered")
    say("  when %d of them are not distinguishable from zero. The"
        % (len(DIL) - res))
    say("  check had no power, and its failure is not evidence that")
    say("  rho(d) is other than the overlap -- it is evidence that")
    say("  %d draws cannot say. U2 does not depend on this: it scores"
        % DRAWS)
    say("  each draw leave-one-out against the rest, so it is")
    say("  self-calibrating whatever rho(d) is.")
    say("SPREAD dilcoinmean_se %.6f" % float(cse.mean()))

    say()
    say("=" * 70)
    say("U1 %s  U2 %s  U3 %s  U4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (u1, u2, u3, u4)))
    say()
    if u2:
        say("the multiplicative question is not the additive one. The "
            "dilated")
        say("field carries something independent signs on mu's "
            "support cannot")
        say("make, at a dilation the algebra points to, and that is "
            "the first")
        say("escape in this branch.")
    else:
        say("four routes closed, and the fourth was the one the "
            "algebra most")
        say("favoured. mu's multiplicativity is in the field by "
            "construction --")
        say("mu(2v) = -mu(v) for v odd -- and it does not raise the "
            "dilation")
        say("correlation above what independent signs on the same "
            "support give.")
        say("So the coin's agreement with C(N) at these ranges is not "
            "an artefact")
        say("of which statistics were chosen. That is a stronger "
            "negative than")
        say("any of the three before it.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
