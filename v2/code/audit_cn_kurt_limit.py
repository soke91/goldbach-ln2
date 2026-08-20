# -*- coding: utf-8 -*-
r"""
Does the decay have a floor?  Two octaves the fit has never seen

WHAT IS AT STAKE

rem:cnkurt measured the excess kurtosis of G(N) = C(N)/sqrt(V(N)) over
seven octaves and found it real at every one -- outside its coin
ensemble by at least 5.0 standard deviations -- and decaying, with the
seven points lying on a line in log-log at kurtosis ~ N^-0.7312.  That
made it a finite-N effect rather than a law, and OPEN.md item 1 says
what replaces the question: **seven octaves cannot tell decay to zero
from decay to a positive limit.**  A positive limit would make the
non-Gaussianity of C(N) a property of C and not of the range, which is
what a law needs.

Fitting kurtosis = A N^-a + L on seven points is a question about
functional form, and this repository has been shown four times over
(rem:ladderdegree, rem:deficitregion, rem:maskformreach,
rem:maskrivals) that it cannot settle those.  So the limit is not
approached by fitting.  It is approached the way this repository's
strongest results were: **the pure power law, with no limit, is made
to forecast two octaves it has never seen, and then they are
measured.**  If they land where a floorless law says, no floor is
needed at this range; if they land above, the floor is what put them
there.

The seven measured octaves are read from results/audit_cn_kurt_drift.txt
and not recomputed.  What is computed here is b = 24 and b = 25, so
N runs to 6.71e7, and the coin ensemble at those two bands.

BACKS: Remark {#rem:cnkurtlimit} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  THE GATE.  Re-measured at the new sieve top, band b = 23
      reproduces its POINT cnkurt_23 to three decimals.  The two runs
      differ in sieve top and transform length and must agree.
  E2  **The floorless law forecasts them.**  Fitted on b = 17..23
      alone, log(kurtosis) against log2 N predicts b = 24 and b = 25,
      and both measured values land inside the forecast's own
      prediction interval at two standard errors -- the interval
      built from the seven-point fit's residual scatter and the
      leverage of the extrapolation, not from the new points.
  E3  A limit does not resolve: fitting kurtosis = A N^-a + L on all
      nine octaves leaves L within two standard errors of zero.
      **This is the outcome this repository's history predicts** and
      it is the weaker one.
  E4  The control is a constant offset, not a drift: the coin's mean
      excess kurtosis at b = 24 and b = 25 stays inside the range
      rem:cnkurt measured at b = 17..23, which is -0.02103 to
      -0.01376.

REFUTATION RULE (fixed before the run)

  E1  REFUTED outside three decimals; nothing below is reported.
  E2  **REFUTED if either octave lands outside.**  Above the interval
      is the interesting direction and would say the decay is
      flattening, which is what a positive limit looks like from
      below.  Below is also possible and would say the decay is
      steepening, and then neither a power law nor a floor describes
      it.  The interval is two standard errors and is fixed here, so
      a near miss is a miss.
  E3  REFUTED if L resolves away from zero.  A resolved positive L is
      the strong outcome -- the non-Gaussianity would have a floor
      and C(N) would be non-Gaussian at every scale -- and a resolved
      negative L would mean the three-parameter family is fitting
      noise, since a negative floor is not a possible limit for a
      quantity that is positive at every octave measured.
  E4  REFUTED if either new band's coin mean leaves that range.  Then
      D4's failure in rem:cnkurt is a drift and not an offset, and
      every z in that remark is scored against a moving baseline.

  E2 and E3 can disagree, and if they do E2 is the one to believe.
  A forecast tested out of sample is evidence; a third parameter
  fitted inside the sample is a description.  If E2 holds and E3 is
  refuted, what that means is that a floor can be fitted but was not
  needed to predict the two new points, and the honest reading is
  that nine octaves still cannot see one.

  WHAT THIS CANNOT DO.  Two more octaves is a factor of four in N.
  A limit smaller than the fit's residual scatter at b = 25 --
  which is what a floor below about 0.02 would be -- cannot be seen
  by this design however it comes out, and no bound below that is
  claimed.
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
OUT = os.path.join(ROOT, "results", "audit_cn_kurt_limit.py".replace(".py", ".txt"))
SRC = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

OLDLO, OLDHI = 17, 23
NEWLO, NEWHI = 24, 25
NMAX = 1 << (NEWHI + 1)
DRAWS = 32
SEED = 20260822
DEC = 3
ZCAP = 2.0
COINLO, COINHI = -0.02103, -0.01376


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


def ols(x, y):
    A = np.column_stack([np.ones(len(y)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    dof = len(y) - 2
    s2 = float((r ** 2).sum()) / dof
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, np.sqrt(np.diag(cov)), math.sqrt(s2), cov


def read_points():
    src = io.open(SRC, encoding="utf-8").read()
    out = {}
    for b in range(OLDLO, OLDHI + 1):
        m = re.search(r"^POINT cnkurt_%d ([-\d.]+)\s*$" % b, src, re.M)
        if not m:
            raise SystemExit("no POINT cnkurt_%d in %s" % (b, SRC))
        out[b] = float(m.group(1))
    return out


HEAD = [
    "STATISTIC: the excess kurtosis of G(N) = C(N)/sqrt(V(N)) at two",
    "           octaves the seven-point power law has never seen,",
    "           against that law's out-of-sample forecast; and the",
    "           limit L in kurtosis = A N^-a + L fitted on all nine.",
    "FIELD: even N in (2^b, 2^(b+1)] for b = %d and %d, so N runs"
    % (NEWLO, NEWHI),
    "       to %d, with Lambda and mu sieved to that top and both"
    % NMAX,
    "       convolutions taken by FFT over the whole range. The seven",
    "       octaves b = %d..%d are READ from" % (OLDLO, OLDHI),
    "       results/audit_cn_kurt_drift.txt and not recomputed; band",
    "       b = %d is recomputed here as the gate." % OLDHI,
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

    old = read_points()
    for b in range(OLDLO, OLDHI + 1):
        say("READ audit_cn_kurt_drift.txt POINT cnkurt_%d %.5f"
            % (b, old[b]))
    say("  the seven octaves this forecast is fitted on, read from "
        "that file")
    say("PRINTBOUND audit_cn_kurt_limit %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))

    xo = np.array([b + 0.5 for b in range(OLDLO, OLDHI + 1)])
    yo = np.log(np.array([old[b] for b in range(OLDLO, OLDHI + 1)]))
    c, se, rms, cov = ols(xo, yo)
    say()
    say("  the floorless law, fitted on b = %d..%d only"
        % (OLDLO, OLDHI))
    say("  log(kurtosis) = %+.5f %+.5f * log2(N), residual r.m.s. "
        "%.5f" % (c[0], c[1], rms))
    say("  so kurtosis ~ N^%+.4f" % (c[1] / math.log(2)))
    say("SCATTER slope_audit_cn_kurt_limit %.5f" % rms)

    fc = {}
    for b in (NEWLO, NEWHI):
        xv = np.array([1.0, b + 0.5])
        mu_ = float(xv.dot(c))
        var = float(xv.dot(cov).dot(xv)) + rms ** 2
        fc[b] = (mu_, math.sqrt(var))
        say("  forecast at b = %d: log(kurtosis) = %+.5f +- %.5f, so "
            "%.5f to %.5f"
            % (b, mu_, math.sqrt(var),
               math.exp(mu_ - ZCAP * math.sqrt(var)),
               math.exp(mu_ + ZCAP * math.sqrt(var))))
        say("BRACKET cnkurtlimit_%d %.5f %.5f"
            % (b, math.exp(mu_ - ZCAP * math.sqrt(var)),
               math.exp(mu_ + ZCAP * math.sqrt(var))))

    say()
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

    gate_k = excess_kurt(C[band(OLDHI)] / rootV[band(OLDHI)])
    new = {b: excess_kurt(C[band(b)] / rootV[band(b)])
           for b in (NEWLO, NEWHI)}
    del C

    # -------------------------------------------------------------- E1
    say()
    say("E1  does this run reproduce the band it shares?")
    say("  b = %d: here %.5f against its %.5f"
        % (OLDHI, gate_k, old[OLDHI]))
    e1 = abs(round(gate_k, DEC) - round(old[OLDHI], DEC)) < 10.0 ** (-DEC) / 2
    say("  E1 %s   (cap: %d decimals)"
        % ("hold" if e1 else "REFUTED", DEC))
    if not e1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- E2
    say()
    say("E2  does the floorless law forecast the two new octaves?")
    say("       b      even N      forecast          measured      "
        "z")
    e2 = True
    for b in (NEWLO, NEWHI):
        m_, s_ = fc[b]
        z = (math.log(new[b]) - m_) / s_
        e2 &= abs(z) <= ZCAP
        say("      %2d  %10d  %.5f to %.5f  %.5f  %+6.2f"
            % (b, len(band(b)), math.exp(m_ - ZCAP * s_),
               math.exp(m_ + ZCAP * s_), new[b], z))
        say("POINT cnkurt_%d %.5f" % (b, new[b]))
        say("TSTAT cnkurtlimit_z%d %.2f" % (b, z))
        say("SPREAD cnkurtlimit_z%d %.5f" % (b, s_))
        if abs(z) < 2.0:
            say("UNRESOLVED SIGN cnkurtlimit_z%d" % b)
    say("SCALES 2")
    say("  E2 %s   (cap: both inside two standard errors)"
        % ("hold" if e2 else "REFUTED"))

    # -------------------------------------------------------------- E3
    say()
    say("E3  does a limit resolve on all nine octaves?")
    bs = list(range(OLDLO, OLDHI + 1)) + [NEWLO, NEWHI]
    xa = np.array([b + 0.5 for b in bs])
    ya = np.array([old.get(b, new.get(b)) for b in bs])
    best = None
    for a in np.linspace(0.05, 2.0, 3901):
        w = np.exp(-a * math.log(2.0) * xa)
        A = np.column_stack([w, np.ones(len(ya))])
        p, *_ = np.linalg.lstsq(A, ya, rcond=None)
        r = ya - A.dot(p)
        ss = float((r ** 2).sum())
        if best is None or ss < best[0]:
            best = (ss, a, p, A, r)
    ss, a_hat, p_hat, A_hat, r_hat = best
    dof = len(ya) - 3
    s2 = ss / dof
    covL = s2 * np.linalg.inv(A_hat.T.dot(A_hat))
    seL = math.sqrt(covL[1, 1])
    tL = p_hat[1] / seL
    say("  best exponent a = %.4f, A = %+.5f, L = %+.5f +- %.5f, "
        "t = %.2f" % (a_hat, p_hat[0], p_hat[1], seL, tL))
    say("TSTAT cnkurtlimit_L %.2f" % tL)
    say("SPREAD cnkurtlimit_L %.5f" % seL)
    if abs(tL) < 2.0:
        say("UNRESOLVED SIGN cnkurtlimit_L")
    e3 = abs(tL) < ZCAP
    say("  residual r.m.s. %.5f against the floorless fit's %.5f "
        "in log" % (math.sqrt(s2), rms))
    say("  E3 %s   (cap: |t| < %.0f)"
        % ("hold" if e3 else "REFUTED", ZCAP))

    # -------------------------------------------------------------- E4
    say()
    say("E4  is the control an offset or a drift?")
    say("  %d fresh sign patterns at each new band" % DRAWS)
    rng = np.random.default_rng(SEED)
    coin = {b: [] for b in (NEWLO, NEWHI)}
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L), L)[: NMAX + 1]
        for b in (NEWLO, NEWHI):
            coin[b].append(excess_kurt(Cc[band(b)] / rootV[band(b)]))
        del Cc
    e4 = True
    say("       b   coin mean    coin sd        z of the real point")
    for b in (NEWLO, NEWHI):
        cm = float(np.mean(coin[b]))
        cs = float(np.std(coin[b], ddof=1))
        e4 &= COINLO <= cm <= COINHI
        say("      %2d  %+9.5f  %9.5f  %18.2f"
            % (b, cm, cs, (new[b] - cm) / cs))
        say("POINT cnkurtcoin_%d %.5f" % (b, cm))
    say("  the range rem:cnkurt measured is %+.5f to %+.5f"
        % (COINLO, COINHI))
    say("  E4 %s   (cap: both inside it)"
        % ("hold" if e4 else "REFUTED"))

    say()
    say("=" * 70)
    say("E1 %s  E2 %s  E3 %s  E4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (e1, e2, e3, e4)))
    say()
    if e2 and e3:
        say("no floor is needed and none can be seen. The floorless "
            "law predicted")
        say("two octaves it had never seen and they landed, and a "
            "third parameter")
        say("fitted on all nine does not resolve. Nine octaves cannot "
            "tell decay")
        say("to zero from decay to a limit below the scatter, and "
            "this design")
        say("says so rather than choosing.")
    elif not e2:
        say("the forecast missed, which is the outcome that carries "
            "information:")
        say("the decay is not the power law the seven octaves "
            "described, and the")
        say("direction of the miss is printed above.")
    elif e2 and not e3:
        say("the forecast held and a floor still fits. A forecast "
            "tested out of")
        say("sample is evidence and a parameter fitted inside it is a "
            "description,")
        say("so the reading is that a floor can be fitted and was not "
            "needed.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
