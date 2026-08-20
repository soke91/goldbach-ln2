# -*- coding: utf-8 -*-
r"""
Does the 105 | N escape persist, and is it a shift or a spread?

WHAT IS AT STAKE

rem:cnclass found the first escape in this branch after five closed
routes.  Cutting even N by which of 3, 5, 7 divide them, the class
105 | N has E[G^2] = 4.37497 against a coin ensemble at 1.02257, a z
of 9.01 that none of 32 draws reaches, while every other class sits
inside its control.

**Most of that is not spread but shift.**  The same class has mean
G = -1.70728, so mean^2 = 2.9148 of the 4.37497, leaving a variance of
1.4602 and a standard deviation of 1.2084.  The class is not a few
enormous |G| dragging a second moment; it is 4993 values of N sitting
about 1.7 standard deviations low together.  C(N) has a systematic
component of the size of its own fluctuation scale at N divisible by
105.

Two things follow that rem:cnclass did not test, and they decide what
the escape is worth.  Does it survive at other N -- everything else
this branch found real turned out to be a finite-N effect decaying
like a power -- and does the shift grow, hold, or fade as N grows?  A
shift that holds is a structural statement about C(N).  A shift that
fades is one more finite-N separation on the pile.

BACKS: Remark {#rem:cnclassreach} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  A1  THE GATE.  At b = 20 the class 105 | N reproduces rem:cnclass's
      POINT classm2_7 to three decimals.
  A2  **The escape is not one band's.**  At every band b = 17..22 the
      105 | N class's E[G^2] lies outside the coin ensemble, with none
      of the 32 draws' own largest |z| over classes reaching the real
      arm's at that band.
  A3  The shift has a fixed sign: the class mean of G at 105 | N is
      negative at every band.  A sign that flips would make it a
      fluctuation of the class rather than a property of it.
  A4  **And it does not fade.**  Regressing log|mean G at 105 | N| on
      log2 N over the six bands gives a slope that is not resolved
      negative -- either positive, or |t| < 2.

REFUTATION RULE (fixed before the run)

  A1  REFUTED outside three decimals; nothing below is reported.
  A2  REFUTED if any band has a coin draw reaching the real arm.
      **The small bands are where to expect it**: at b = 17 the class
      holds about 624 values of N against 4993 at b = 20, so a
      failure there is about sample size and the class count is
      printed beside every z.
  A3  REFUTED by a positive class mean at any band.
  A4  REFUTED if the slope is negative and resolved at |t| > 2.  That
      is the outcome that would make this the sixth finite-N
      separation rather than a structural one, and it is the outcome
      this branch's history predicts -- rem:cnkurtlimit, rem:cnskew
      and rem:cnwherereach all found powers.  If A4 breaks, the
      exponent is the finding and is reported as such.
      **The unresolved case is named and it counts as A4 holding**,
      which is deliberate and is the weakest way this prediction can
      be right: six bands too noisy to resolve a slope would say only
      that this design cannot see a fade, not that there is none, and
      the remark must then say "not resolved" and never "does not
      fade".  (This sentence was added after the run, to satisfy the
      gate's M9 check; the rule it describes is unchanged -- the
      condition was already "negative AND resolved" -- and the
      measured |t| was 89.38, so no reading of the outcome turns on
      it.)

  A3 and A4 are different claims and can part: a shift that keeps its
  sign while shrinking is still a fading shift, and A4 is the one that
  decides whether the escape means anything asymptotically.

  WHAT THIS CANNOT DO.  Six bands over N to 8.4e6, and one class.
  Nothing here attributes the shift to a formula, tests whether it is
  a singular series, or asks what happens at moduli beyond 3, 5, 7.
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
OUT = os.path.join(ROOT, "results", "audit_cn_class_reach.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_class.txt")

BLO, BHI = 17, 22
GATEB = 20
NMAX = 1 << (BHI + 1)
QS = (3, 5, 7)
FULL = 7                   # the class index with all three dividing
DRAWS = 32
SEED = 20260829
DEC = 3
ZCAP = 2.0


assert (1 << (BHI + 1)) <= NMAX and BLO <= GATEB <= BHI


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


def ols(x, y):
    A = np.column_stack([np.ones(len(y)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    s2 = float((r ** 2).sum()) / (len(y) - 2)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, np.sqrt(np.diag(cov))


def published_m2():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT classm2_%d ([-\d.]+)\s*$" % FULL, src, re.M)
    if not m:
        raise SystemExit("no POINT classm2_%d in %s" % (FULL, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the first and second moments of",
    "           G(N) = C(N)/sqrt(V(N)) inside the class of even N",
    "           divisible by %d, at six bands, against an ensemble of"
    % (QS[0] * QS[1] * QS[2],),
    "           %d sign patterns on mu's support; and the slope of" % DRAWS,
    "           log|class mean| against log2 N.",
    "FIELD: even N in (2^b, 2^(b+1)] for b = %d..%d, cut into the"
    % (BLO, BHI),
    "       eight classes given by which of %s divide N; Lambda and" % (QS,),
    "       mu sieved once to %d and both convolutions taken by FFT" % NMAX,
    "       over the whole range, then sliced per band and class. The",
    "       gate re-measures band b = %d's full class against" % GATEB,
    "       rem:cnclass's own value.",
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

    pub = published_m2()
    say("READ audit_cn_class.txt POINT classm2_%d %.5f" % (FULL, pub))
    say("  the class this run has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_class_reach %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d fresh sign patterns over the whole ladder" % DRAWS)

    say("sieving to %d, which is 2^%d" % (NMAX, BHI + 1))
    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    mu2 = muf ** 2
    del mu
    L = 1 << (2 * NMAX).bit_length()
    say("  transform length %d" % L)
    FL = np.fft.rfft(lam, L)
    C = np.fft.irfft(FL * np.fft.rfft(muf, L), L)[: NMAX + 1]
    V = np.fft.irfft(np.fft.rfft(lam ** 2, L) * np.fft.rfft(mu2, L),
                     L)[: NMAX + 1]
    del lam

    bands = list(range(BLO, BHI + 1))
    Ns, rootV, masks = {}, {}, {}
    for b in bands:
        nn = np.arange((1 << b) + 2, (1 << (b + 1)) + 1, 2,
                       dtype=np.int64)
        Ns[b] = nn
        rootV[b] = np.sqrt(np.maximum(V[nn], 1e-300))
        lab = np.zeros(len(nn), dtype=np.int64)
        for i, q in enumerate(QS):
            lab |= ((nn % q) == 0).astype(np.int64) << i
        masks[b] = [lab == k for k in range(1 << len(QS))]
    del V

    def cls_moments(b, arr):
        g = arr[Ns[b]] / rootV[b]
        m1 = np.array([float(g[m].mean()) for m in masks[b]])
        m2 = np.array([float((g[m] ** 2).mean()) for m in masks[b]])
        return m1, m2

    real1, real2 = {}, {}
    for b in bands:
        real1[b], real2[b] = cls_moments(b, C)
    del C

    # -------------------------------------------------------------- A1
    say()
    say("A1  does this run reproduce the class it shares?")
    say("  b = %d, class 105|N: E[G^2] here %.5f against its %.5f"
        % (GATEB, real2[GATEB][FULL], pub))
    a1 = abs(round(real2[GATEB][FULL], DEC) - round(pub, DEC)) \
        < 10.0 ** (-DEC) / 2
    say("  A1 %s   (cap: %d decimals)"
        % ("hold" if a1 else "REFUTED", DEC))
    if not a1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("  shift or spread, at the gate band")
    mg = real1[GATEB][FULL]
    say("  mean %+.5f, so mean^2 %.5f of the E[G^2] %.5f, leaving a"
        % (mg, mg * mg, real2[GATEB][FULL]))
    say("  variance %.5f and a standard deviation %.5f"
        % (real2[GATEB][FULL] - mg * mg,
           math.sqrt(max(real2[GATEB][FULL] - mg * mg, 0.0))))
    say("POINT classreach_meansq %.5f" % (mg * mg))
    say("POINT classreach_var %.5f"
        % (real2[GATEB][FULL] - mg * mg))

    coin1 = {b: np.zeros((DRAWS, 1 << len(QS))) for b in bands}
    coin2 = {b: np.zeros((DRAWS, 1 << len(QS))) for b in bands}
    for d in range(DRAWS):
        eps = rng_draw(d)
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L),
                          L)[: NMAX + 1]
        for b in bands:
            coin1[b][d], coin2[b][d] = cls_moments(b, Cc)
        del Cc

    def zsc(t, pool):
        return (t - pool.mean(axis=0)) / pool.std(axis=0, ddof=1)

    say()
    say("       b   class N     mean G     E[G^2]    coin mean"
        "     z(m2)   coin max |z|   beat")
    beats, means, z2s = [], [], []
    for b in bands:
        zr = zsc(real2[b], coin2[b])
        mc = []
        for d in range(DRAWS):
            rest = np.delete(coin2[b], d, axis=0)
            mc.append(float(np.abs(zsc(coin2[b][d], rest)).max()))
        mc = np.array(mc)
        mz = float(np.abs(zr).max())
        beat = int((mc >= mz).sum())
        beats.append(beat)
        means.append(real1[b][FULL])
        z2s.append(zr[FULL])
        say("      %2d  %8d  %+9.5f  %+9.5f  %+9.5f  %8.2f  "
            "%5.2f-%5.2f  %5d"
            % (b, int(masks[b][FULL].sum()), real1[b][FULL],
               real2[b][FULL], coin2[b][:, FULL].mean(), zr[FULL],
               mc.min(), mc.max(), beat))
        say("POINT classreach_m1_%d %.5f" % (b, real1[b][FULL]))
        say("POINT classreach_m2_%d %.5f" % (b, real2[b][FULL]))
        say("POINT classreach_z_%d %.5f" % (b, zr[FULL]))
    say("SCALES %d" % len(bands))

    # -------------------------------------------------------------- A2
    say()
    say("A2  does the escape survive at every band?")
    a2 = all(x == 0 for x in beats)
    say("  coin draws reaching the real arm, band by band: %s"
        % ", ".join(str(x) for x in beats))
    say("TSTAT classreach_minz %.2f" % min(z2s, key=abs))
    say("SPREAD classreach_minz %.5f"
        % float(coin2[bands[int(np.argmin(np.abs(z2s)))]]
                [:, FULL].std(ddof=1)))
    say("  A2 %s   (cap: none at any band)"
        % ("hold" if a2 else "REFUTED"))

    # -------------------------------------------------------------- A3
    say()
    say("A3  does the shift keep its sign?")
    a3 = all(m < 0 for m in means)
    say("  class means %s" % ", ".join("%+.5f" % m for m in means))
    say("  A3 %s   (cap: negative at every band)"
        % ("hold" if a3 else "REFUTED"))

    # -------------------------------------------------------------- A4
    say()
    say("A4  does the shift fade with N?")
    x = np.array([b + 0.5 for b in bands])
    if a3:
        y = np.log(np.abs(np.array(means)))
        c, se = ols(x, y)
        t = c[1] / se[1]
        say("  log|class mean| on log2 N: slope %+.5f +- %.5f, "
            "t = %.2f" % (c[1], se[1], t))
        say("  so the shift goes like N^%+.4f" % (c[1] / math.log(2)))
        say("POINT classreach_power %.5f" % (c[1] / math.log(2)))
        say("TSTAT classreach_slope %.2f" % t)
        say("SPREAD classreach_slope %.5f" % (x.max() - x.min()))
        if abs(t) < 2.0:
            say("UNRESOLVED SIGN classreach_slope")
        a4 = not (c[1] < 0 and abs(t) > ZCAP)

        say()
        say("  NOT PRE-REGISTERED, reported because a t of %.0f with "
            "no residual" % abs(t))
        say("  printed is exactly what this repository flags "
            "elsewhere. No cap")
        say("  above is changed.")
        res = y - (c[0] + c[1] * x)
        say("  residuals in log|mean|, band by band: %s"
            % ", ".join("%+.5f" % r for r in res))
        say("  r.m.s. %.5f against a mean |log| of %.4f"
            % (float(np.sqrt((res ** 2).mean())),
               float(np.abs(y).mean())))
        say("SCATTER slope_audit_cn_class_reach %.5f"
            % float(np.sqrt((res ** 2).mean())))
        say("  Six bands cut from one sieve are not independent "
            "samples, so the")
        say("  standard error above is a lower bound and no "
            "significance is")
        say("  claimed from the t. What the fit supports is that the "
            "six points")
        say("  lie on a line to that r.m.s.")
        zz = np.array(z2s)
        cz, cse = ols(x, np.log(np.abs(zz)))
        say("  and the separation itself: |z| goes like N^%+.4f, so "
            "from %.2f at" % (cz[1] / math.log(2), abs(zz[-1])))
        say("  b = %d it would reach 3 near b = %.1f -- an "
            "extrapolation, not a"
            % (bands[-1],
               bands[-1] + math.log(abs(zz[-1]) / 3.0) / abs(cz[1])))
        say("  measurement, and printed so the next band to try is "
            "named.")
        say("POINT classreach_zpower %.5f" % (cz[1] / math.log(2)))
    else:
        say("  A3 failed, so no log is taken and A4 is not decided "
            "here")
        a4 = False
    say("  A4 %s   (cap: not a resolved negative slope)"
        % ("hold" if a4 else "REFUTED"))

    say()
    say("=" * 70)
    say("A1 %s  A2 %s  A3 %s  A4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (a1, a2, a3, a4)))
    say()
    if a2 and a3 and a4:
        say("the escape is structural as far as this reaches. C(N) "
            "sits low at")
        say("N divisible by 105 at every band measured, by an amount "
            "that does")
        say("not fade, and no coin draw makes it. That is the first "
            "thing in")
        say("this branch that growing N does not erase.")
    elif a2 and a3 and not a4:
        say("the escape is real at every band and fades with N, so it "
            "joins the")
        say("pile: a sixth finite-N separation. The exponent is above "
            "and is")
        say("what this run adds.")
    else:
        say("the escape does not survive at every band. Where it "
            "fails and by")
        say("how much is in the table, and the class counts beside "
            "it are the")
        say("first thing to read.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


_RNG = np.random.default_rng(SEED)


def rng_draw(_d):
    return _RNG.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1


if __name__ == "__main__":
    sys.exit(main())
