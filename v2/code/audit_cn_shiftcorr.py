# -*- coding: utf-8 -*-
r"""
Off the margin: does the field's correlation across N escape the coin?

WHAT IS AT STAKE

Two statistics of G(N) = C(N)/sqrt(V(N)) have now been measured against
the coin and both closed.  rem:cnkurtlimit found the excess kurtosis
meeting its control at b = 25 because the control does not shrink;
rem:cnskew found the skewness meeting its control at b = 21 even
though the control does shrink, because the signal shrinks five times
faster.  Both decay like N^-0.6 to N^-0.7 against a control that decays
like N^-0.12 or not at all.

**Every one of those is a functional of the marginal collection
{G(N)}, and the coin has mu's support and independent signs, which is
all a marginal needs.**  What the coin does not have is mu's
multiplicative structure linking different N.  That is not a moment of
the margin; it is a correlation of the field, and this branch has
never measured one.

THE QUANTITY, AND WHY IT IS THE RIGHT ONE

Write C(N) = sum_v mu(v) Lambda(N-v).  Then

    C(N) C(N') = sum_{v,w} mu(v) mu(w) Lambda(N-v) Lambda(N'-w),

and averaging the coin over its signs kills every off-diagonal term:
E_eps[C_coin(N) C_coin(N')] is exactly the diagonal
sum_v mu^2(v) Lambda(N-v) Lambda(N'-v).  **The off-diagonal is where
mu's own correlations live, and it is precisely what the coin cannot
make.**  So the shift correlation

    rho(h) = mean over N in a band of G(N) G(N+h)

is a place the two arms can differ for a reason, rather than by
accident.

The design compares two octaves so the N-dependence is visible in the
same run: b = 21, where the skewness had already closed, and b = 23,
where the kurtosis was at z = 5.0 and the skewness at 1.20.  Shifts
run over the 32 even h from 2 to 64.

BACKS: Remark {#rem:cnshift} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  THE GATE.  Band b = 23's excess kurtosis reproduces
      rem:cnkurt's POINT cnkurt_23 to three decimals.
  T2  **The field escapes the coin.**  At b = 23 the real arm's
      largest |z| over the 32 shifts exceeds every one of the 32 coin
      draws' own largest |z|, each scored leave-one-out against the
      rest of the ensemble.  This is a max-statistic comparison and
      not 32 separate tests, because 32 shifts at |z| > 3 would fire
      on noise about one time in twelve.
  T3  And it does not close with N: the real arm's largest |z| at
      b = 23 is at least as large as at b = 21.  Every marginal
      statistic fell by a factor of about four over those two
      octaves.
  T4  **This statistic is sample-limited, not walled.**  The coin's
      spread in rho(h) falls from b = 21 to b = 23 by a factor
      between 1.4 and 2.8 -- the band holds four times as many N, so
      a sample-limited spread would fall by 2, and the interval is
      wide enough to allow the correlation between neighbouring N to
      cost something.  That is the property the kurtosis lacked
      outright and the skewness had only weakly.

REFUTATION RULE (fixed before the run)

  T1  REFUTED outside three decimals; nothing below is reported.
  T2  REFUTED if any coin draw's largest |z| reaches the real arm's.
      Then the field's shift correlation is what independent signs on
      mu's support already give, and the change of object bought
      nothing -- which would be the third closed route in this branch
      and would say the marginal was never the problem.
  T3  REFUTED if the largest |z| falls between the two octaves.  A
      fall would mean this statistic closes like the others and only
      the scale of the closing has changed.  It is the outcome to
      expect if the escape in T2 is a finite-N effect of the same
      family.
  T4  REFUTED outside 1.4 to 2.8.  Below 1.4 means the spread is set
      by correlation across N rather than by sample size, which is
      the kurtosis's wall again.  Above 2.8 would mean the shifts are
      more independent than the count of N allows, and the error bar
      is then not the one this design assumes.

  T2 and T3 can disagree.  If T2 holds and T3 does not, the reading
  is that the field has correlation the coin cannot make and that it
  is fading like everything else in this branch; the escape would be
  real and the route still finite.  Which of those it is matters more
  than whether any single |z| is large, and both are recorded.

  WHAT THIS CANNOT DO.  Two octaves and thirty-two shifts.  Nothing
  here identifies WHICH arithmetic structure produces a departure, and
  a departure at a shift h is not attributed to any property of h; the
  per-shift table is printed so a reader can look, not so this script
  can claim.
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
OUT = os.path.join(ROOT, "results", "audit_cn_shiftcorr.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

BLOW, BHIGH = 21, 23
NMAX = 1 << (BHIGH + 1)
SHIFTS = 32
DRAWS = 32
SEED = 20260824
DEC = 3
RLO, RHI = 1.4, 2.8


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


def shiftcorr(g, kmax):
    """rho at index shifts 1..kmax, i.e. even N-shifts 2..2*kmax"""
    d = g - g.mean()
    s = float((d ** 2).sum())
    return np.array([float(d[:-k].dot(d[k:])) / s
                     for k in range(1, kmax + 1)])


def published_kurt():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT cnkurt_%d ([-\d.]+)\s*$" % BHIGH, src, re.M)
    if not m:
        raise SystemExit("no POINT cnkurt_%d in %s" % (BHIGH, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the shift correlation rho(h) = mean over N in a band",
    "           of G(N) G(N+h), for the %d even shifts h = 2..%d,"
    % (SHIFTS, 2 * SHIFTS),
    "           against an ensemble of %d sign patterns on mu's" % DRAWS,
    "           support; compared by the largest |z| over shifts,",
    "           with each coin draw scored leave-one-out against the",
    "           rest so the comparison is one max-statistic and not",
    "           %d separate tests." % SHIFTS,
    "FIELD: even N in (2^b, 2^(b+1)] for b = %d and %d, with Lambda"
    % (BLOW, BHIGH),
    "       and mu sieved once to %d and both convolutions taken" % NMAX,
    "       by FFT over the whole range, then sliced per band. The",
    "       gate re-measures band b = %d's excess kurtosis against"
    % BHIGH,
    "       rem:cnkurt's own value.",
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
        % (BHIGH, pub))
    say("  the band this run has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_shiftcorr %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d shifts, %d fresh sign patterns" % (SHIFTS, DRAWS))

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

    bands = (BLOW, BHIGH)
    idx = {b: band(b) for b in bands}
    Greal = {b: C[idx[b]] / rootV[idx[b]] for b in bands}
    gate_k = excess_kurt(Greal[BHIGH])
    rho_real = {b: shiftcorr(Greal[b], SHIFTS) for b in bands}
    del C, Greal

    # -------------------------------------------------------------- T1
    say()
    say("T1  does this run reproduce the band it shares?")
    say("  b = %d excess kurtosis: here %.5f against its %.5f"
        % (BHIGH, gate_k, pub))
    t1 = abs(round(gate_k, DEC) - round(pub, DEC)) < 10.0 ** (-DEC) / 2
    say("  T1 %s   (cap: %d decimals)"
        % ("hold" if t1 else "REFUTED", DEC))
    if not t1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    rho_coin = {b: np.zeros((DRAWS, SHIFTS)) for b in bands}
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L), L)[: NMAX + 1]
        for b in bands:
            rho_coin[b][d] = shiftcorr(Cc[idx[b]] / rootV[idx[b]],
                                       SHIFTS)
        del Cc

    def zscores(target, pool):
        """target against pool's mean and sd, per shift"""
        m = pool.mean(axis=0)
        s = pool.std(axis=0, ddof=1)
        return (target - m) / s

    maxz_real, maxz_coin, spread = {}, {}, {}
    for b in bands:
        zr = zscores(rho_real[b], rho_coin[b])
        maxz_real[b] = float(np.abs(zr).max())
        mc = []
        for d in range(DRAWS):
            rest = np.delete(rho_coin[b], d, axis=0)
            mc.append(float(np.abs(zscores(rho_coin[b][d],
                                           rest)).max()))
        maxz_coin[b] = np.array(mc)
        spread[b] = float(rho_coin[b].std(axis=0, ddof=1).mean())
        say()
        say("  b = %d, %d values of N, mean coin spread in rho %.6f"
            % (b, len(idx[b]), spread[b]))
        say("POINT shiftspread_%d %.6f" % (b, spread[b]))
        say("        h     rho real     coin mean      coin sd      z")
        for k in range(SHIFTS):
            say("      %3d  %+11.6f  %+11.6f  %11.6f  %6.2f"
                % (2 * (k + 1), rho_real[b][k],
                   rho_coin[b][:, k].mean(),
                   rho_coin[b][:, k].std(ddof=1), zr[k]))
        say("POINT shiftmaxz_%d %.5f" % (b, maxz_real[b]))
        say("  largest |z| here %.2f; the coin draws' own largest run "
            "%.2f to %.2f" % (maxz_real[b], maxz_coin[b].min(),
                              maxz_coin[b].max()))
    say("SCALES 2")

    # -------------------------------------------------------------- T2
    say()
    say("T2  does the field escape the coin at the top octave?")
    beat = int((maxz_coin[BHIGH] >= maxz_real[BHIGH]).sum())
    say("  coin draws whose own largest |z| reaches the real arm's: "
        "%d of %d" % (beat, DRAWS))
    t2 = beat == 0
    say("  T2 %s   (cap: none of them)"
        % ("hold" if t2 else "REFUTED"))
    say("TSTAT shiftcorr_maxz %.2f" % maxz_real[BHIGH])
    say("SPREAD shiftcorr_maxz %.5f" % float(maxz_coin[BHIGH].std(ddof=1)))

    # -------------------------------------------------------------- T3
    say()
    say("T3  does it close with N?")
    say("  largest |z| at b = %d is %.2f, at b = %d is %.2f"
        % (BLOW, maxz_real[BLOW], BHIGH, maxz_real[BHIGH]))
    t3 = maxz_real[BHIGH] >= maxz_real[BLOW]
    say("  T3 %s   (cap: the top octave is at least the lower one)"
        % ("hold" if t3 else "REFUTED"))

    # -------------------------------------------------------------- T4
    say()
    say("T4  is the control sample-limited?")
    ratio = spread[BLOW] / spread[BHIGH]
    say("  coin spread falls by a factor %.4f between the octaves; "
        "the band" % ratio)
    say("  holds four times as many N, so sample-limited would be 2")
    t4 = RLO <= ratio <= RHI
    say("POINT shiftratio %.5f" % ratio)
    say("  T4 %s   (cap: %.1f to %.1f)"
        % ("hold" if t4 else "REFUTED", RLO, RHI))

    say()
    say("=" * 70)
    say("T1 %s  T2 %s  T3 %s  T4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (t1, t2, t3, t4)))
    say()
    if t2 and t3 and t4:
        say("the change of object worked. The field's correlation "
            "across N carries")
        say("something independent signs on mu's support cannot make, "
            "it does not")
        say("fade between the two octaves, and its control is limited "
            "by how many")
        say("N are in the band rather than by the field's own "
            "correlation. That is")
        say("the first quantity in this branch that more computation "
            "can sharpen.")
    elif t2 and not t3:
        say("the field's correlation carries something the coin "
            "cannot make, and it")
        say("is fading like everything else here. The escape is real "
            "and the route")
        say("is finite, which is a different sentence from either of "
            "the two")
        say("closed ones.")
    elif not t2:
        say("the shift correlation is what independent signs on mu's "
            "support")
        say("already give. Three routes closed, and the marginal was "
            "never the")
        say("problem -- what the coin reproduces is more than a "
            "margin.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
