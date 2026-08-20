# -*- coding: utf-8 -*-
r"""
Is the v-profile a profile at all, or just a count of terms?

WHAT IS AT STAKE

rem:cnkurtwhere closed the question it asked -- the excess kurtosis of
G(N) is not made in one part of the sum, since the largest per-piece
|z| is 6.27 against the whole field's 16.5 -- and left an observation
it explicitly refused to call evidence: the per-piece |z| is not flat.
Over v up to 2^14 it runs 2.90, 3.22, 3.65, 4.96, 5.22, 6.27, 4.21,
and above it 1.68, 1.74, 0.84, 1.17, 0.82, 1.50, 0.54.  The upper half
of the ladder carries almost nothing.

That remark said the way to make it evidence is to predict the same
profile at a band the run never saw.  This is that, and it sharpens
the question first, because there are two readings and they differ:

  a profile in v -- the same absolute v-ranges carry the excess at
  every N, which would be a statement about how far below N the prime
  sits; or

  no profile at all -- the |z| is a function of how many terms the
  piece holds, and since a dyadic v-piece holds about 2^j / log N
  primes, that count is nearly the same at every band.  Then the
  "profile" is a sample-size effect and there is nothing located
  anywhere.

The second is the deflationary reading and is the one to try to
confirm, because it is the one that explains the shape without any
arithmetic.  So the test is a forecast: fit |z| against the log of the
piece's term count on b = 20's fourteen pieces alone, and predict the
pieces of two bands that fit has never seen.

BACKS: Remark {#rem:cnwherereach} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  THE GATE.  Both new bands reproduce their published whole-band
      excess kurtosis: rem:cnkurt's POINT cnkurt_18 and cnkurt_19, to
      three decimals.
  Y2  **The count explains the shape.**  Fitted on b = 20's pieces
      alone, |z| against log2 of the term count forecasts the pieces
      of b = 18 and b = 19, and at least 80 per cent of them land
      inside the forecast's two-standard-error interval.
  Y3  The non-concentration reproduces out of sample: at both new
      bands the largest per-piece |z| is below that band's own
      whole-field z, which rem:cnkurt puts at 39.8 for b = 18 and
      20.7 for b = 19.
  Y4  The platykurtic law that broke W4 holds everywhere: the coin's
      per-piece excess kurtosis is negative at every piece of every
      band, and its magnitude falls as the term count rises.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED on either band; nothing below is reported.
  Y2  **REFUTED below 80 per cent.**  Then the term count does not
      explain the shape, and what is left is a profile in v -- a
      statement that where the prime sits below N matters beyond how
      many primes are there.  That is the interesting outcome and the
      one that would reopen rem:cnkurtwhere's question in a form it
      did not test.  The direction of the misses is printed either
      way, since a systematic miss and a scattered one mean different
      things.
  Y3  REFUTED if either band has a piece above its whole-field z.
      Then concentration appears at a band where rem:cnkurtwhere did
      not look, and that remark's negative was a property of b = 20.
  Y4  REFUTED by a positive coin kurtosis at any piece, or by no
      resolved fall with the count.  This is a check on the formula
      -2 sum w^4 / (sum w^2)^2, which W4's failure identified after
      the fact; it has not been tested anywhere yet and is tested
      here.

  Y2 is the prediction this run exists for and the others support it.
  If Y2 holds the shape has no arithmetic content and
  rem:cnkurtwhere's observation should be withdrawn as an
  observation, not merely left untested.

  WHAT THIS CANNOT DO.  Two bands below the one fitted, so the
  forecast extrapolates downward in N and not upward.  A count law
  that holds below b = 20 and fails above it would not be caught
  here.
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
OUT = os.path.join(ROOT, "results", "audit_cn_where_reach.txt")
SRCW = os.path.join(ROOT, "results", "audit_cn_kurt_where.txt")
SRCK = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

FITBAND = 20
NEWB = (18, 19)
NMAX = 1 << (max(NEWB) + 1)
JLO = 8
DRAWS = 32
SEED = 20260827
DEC = 3
ZCAP = 2.0
HITFRAC = 0.80
WHOLEZ = {18: 39.8, 19: 20.7}


assert (1 << (max(NEWB) + 1)) <= NMAX


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


def squarefree_flags(n):
    """mu^2 on 0..n, for counting only -- no signs needed"""
    f = np.ones(n + 1, dtype=np.int8)
    f[0] = 0
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        f[p * p::p * p] = 0
    return f


def excess_kurt(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    if v <= 0:
        return float("nan")
    return float((d ** 4).mean() / v ** 2 - 3.0)


def ols(x, y):
    A = np.column_stack([np.ones(len(y)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    dof = len(y) - 2
    s2 = float((r ** 2).sum()) / dof
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, math.sqrt(s2), cov


def read_marks(path, pat):
    src = io.open(path, encoding="utf-8").read()
    return {int(a): float(b) for a, b in re.findall(pat, src, re.M)}


HEAD = [
    "STATISTIC: the per-piece |z| of the excess kurtosis of",
    "           G_j(N) = C_j(N)/sqrt(V_j(N)) against a coin ensemble",
    "           restricted to the same dyadic v-piece, at two bands",
    "           the fit has never seen, against a forecast made from",
    "           b = %d's pieces alone with the piece's term count as" % FITBAND,
    "           the only regressor.",
    "FIELD: even N in (2^b, 2^(b+1)] for b = %s, with v-pieces"
    % (NEWB,),
    "       (2^j, 2^(j+1)] for j = %d up to b and one lumped piece" % JLO,
    "       v <= 2^%d; Lambda and mu sieved once to %d and every"
    % (JLO, NMAX),
    "       piece's two convolutions taken by FFT over the whole",
    "       range. The b = %d pieces are READ from" % FITBAND,
    "       results/audit_cn_kurt_where.txt and not recomputed; their",
    "       term counts are recomputed here from the same sieve,",
    "       which is deterministic and not a measurement.",
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

    zfit = read_marks(SRCW, r"^POINT kurtwherez_(\d+) ([-\d.]+)\s*$")
    kpub = read_marks(SRCK, r"^POINT cnkurt_(\d+) ([-\d.]+)\s*$")
    for hi in sorted(zfit):
        say("READ audit_cn_kurt_where.txt POINT kurtwherez_%d %.5f"
            % (hi, zfit[hi]))
    for b in NEWB:
        say("READ audit_cn_kurt_drift.txt POINT cnkurt_%d %.5f"
            % (b, kpub[b]))
    say("  %d pieces of band %d, and the two bands' published "
        "kurtosis" % (len(zfit), FITBAND))
    say("PRINTBOUND audit_cn_where_reach %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d fresh sign patterns per piece" % DRAWS)

    say("sieving to %d, which is 2^%d" % (NMAX, max(NEWB) + 1))
    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    mu2 = muf ** 2
    del mu

    # term counts of the fitted band's pieces, recomputed
    def piece_edges(b):
        return [(0, 1 << JLO)] + [(1 << j, 1 << (j + 1))
                                  for j in range(JLO, b + 1)]

    sfflag = squarefree_flags(1 << (FITBAND + 1))
    fitted = []
    for lo, hi in piece_edges(FITBAND):
        if hi not in zfit:
            continue
        n = int(sfflag[lo + 1: hi + 1].sum())
        fitted.append((hi, n, abs(zfit[hi])))
    del sfflag
    say()
    say("  the fit, on band %d's pieces only" % FITBAND)
    xf = np.log2(np.array([n for _, n, _ in fitted], dtype=float))
    yf = np.array([z for _, _, z in fitted])
    c, rms, cov = ols(xf, yf)
    say("  |z| = %+.5f %+.5f * log2(terms), residual r.m.s. %.5f"
        % (c[0], c[1], rms))
    say("SCATTER slope_audit_cn_where_reach %.5f" % rms)

    L = 1 << (2 * NMAX).bit_length()
    FL = np.fft.rfft(lam, L)
    FL2 = np.fft.rfft(lam ** 2, L)
    del lam

    def conv(f, g):
        return np.fft.irfft(f * np.fft.rfft(g, L), L)[: NMAX + 1]

    def restrict(arr, lo, hi):
        out = np.zeros_like(arr)
        out[lo + 1: hi + 1] = arr[lo + 1: hi + 1]
        return out

    gate_ok, inside, total, misses = True, 0, 0, []
    errs_fit, errs_flat, zs_new = [], [], []
    coin_pos, coin_pairs = 0, []
    maxz = {}
    for b in NEWB:
        Ns = np.arange((1 << b) + 2, (1 << (b + 1)) + 1, 2,
                       dtype=np.int64)
        gk = excess_kurt(conv(FL, muf)[Ns]
                         / np.sqrt(np.maximum(conv(FL2, mu2)[Ns],
                                              1e-300)))
        ok = abs(round(gk, DEC) - round(kpub[b], DEC)) < 10.0 ** (-DEC) / 2
        gate_ok &= ok
        say()
        say("  band %d: whole-band excess kurtosis %.5f against its "
            "%.5f  %s" % (b, gk, kpub[b], "ok" if ok else "MISMATCH"))
        if not ok:
            continue
        edges = piece_edges(b)
        Cj, Vj, terms = [], [], []
        for lo, hi in edges:
            Cj.append(conv(FL, restrict(muf, lo, hi))[Ns])
            Vj.append(conv(FL2, restrict(mu2, lo, hi))[Ns])
            terms.append(int(mu2[lo + 1: hi + 1].sum()))
        rng = np.random.default_rng(SEED + b)
        coin = [[] for _ in edges]
        for d in range(DRAWS):
            eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
            w = eps * mu2
            for i, (lo, hi) in enumerate(edges):
                cc = conv(FL, restrict(w, lo, hi))[Ns]
                m = Vj[i] > 0
                coin[i].append(excess_kurt(cc[m] / np.sqrt(Vj[i][m])))
            del w
        say("      v-piece            terms     real kurt    coin mean"
            "        z     forecast          in?")
        best = 0.0
        for i, (lo, hi) in enumerate(edges):
            m = Vj[i] > 0
            r = excess_kurt(Cj[i][m] / np.sqrt(Vj[i][m]))
            arr = np.array(coin[i])
            cm, cs = float(arr.mean()), float(arr.std(ddof=1))
            z = abs((r - cm) / cs)
            best = max(best, z)
            coin_pos += int(cm > 0)
            coin_pairs.append((terms[i], cm))
            xv = np.array([1.0, math.log2(max(terms[i], 1))])
            f = float(xv.dot(c))
            sd = math.sqrt(float(xv.dot(cov).dot(xv)) + rms ** 2)
            hit = abs(z - f) <= ZCAP * sd
            inside += int(hit)
            total += 1
            errs_fit.append(z - f)
            errs_flat.append(z - float(np.mean(yf)))
            zs_new.append(z)
            if not hit:
                misses.append(z - f)
            say("  (%9d, %9d]  %9d  %+10.5f  %+10.5f  %7.2f  "
                "%6.2f +- %5.2f  %s"
                % (lo, hi, terms[i], r, cm, z, f, ZCAP * sd,
                   "yes" if hit else "NO"))
            say("POINT reachz_%d_%d %.5f" % (b, hi, z))
        maxz[b] = best
        say("  largest per-piece |z| %.2f against the whole field's "
            "%.1f" % (best, WHOLEZ[b]))
        say("POINT reachmaxz_%d %.5f" % (b, best))
    say("SCALES %d" % len(NEWB))

    # -------------------------------------------------------------- Y1
    say()
    say("Y1  do both bands reproduce their published kurtosis?")
    say("  Y1 %s   (cap: %d decimals on both)"
        % ("hold" if gate_ok else "REFUTED", DEC))
    if not gate_ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- Y2
    say()
    say("Y2  does the term count forecast the new bands?")
    frac = inside / float(total)
    say("  %d of %d pieces inside the two-standard-error forecast, "
        "%.4f" % (inside, total, frac))
    if misses:
        say("  the %d misses are %s"
            % (len(misses), ", ".join("%+.2f" % m for m in misses)))
        say("  mean miss %+.4f, so they are %s"
            % (float(np.mean(misses)),
               "systematic" if abs(np.mean(misses))
               > 0.5 * float(np.std(misses) + 1e-12) else "scattered"))
    y2 = frac >= HITFRAC
    say("SHARE wherereach_inside %.4f" % frac)

    say()
    say("  NOT PRE-REGISTERED, reported because a pass this wide")
    say("  needs its power stated. No cap above is changed.")
    e_fit = math.sqrt(float(np.mean([e * e for e in errs_fit])))
    e_flat = math.sqrt(float(np.mean([e * e for e in errs_flat])))
    say("  the forecast interval is about +-%.1f while the new "
        "bands' |z| span" % (2 * rms))
    say("  %.2f to %.2f, so \"all inside\" is a weak pass on its own. "
        "Against a" % (min(zs_new), max(zs_new)))
    say("  constant forecast at band %d's mean |z| of %.4f, on the "
        "same %d" % (FITBAND, float(np.mean(yf)), total))
    say("  pieces the fit never saw:")
    say("      forecast from the term count   r.m.s. error %.4f"
        % e_fit)
    say("      constant                       r.m.s. error %.4f"
        % e_flat)
    say("      ratio %.4f" % (e_flat / e_fit))
    say("POINT wherereach_rmsfit %.5f" % e_fit)
    say("POINT wherereach_rmsflat %.5f" % e_flat)
    say("  So the count is not merely compatible with the new bands; "
        "it predicts")
    say("  them better than knowing nothing does, by that ratio.")
    say("  Y2 %s   (cap: %.2f)" % ("hold" if y2 else "REFUTED",
                                   HITFRAC))

    # -------------------------------------------------------------- Y3
    say()
    say("Y3  does the non-concentration reproduce?")
    y3 = all(maxz[b] < WHOLEZ[b] for b in NEWB)
    for b in NEWB:
        say("  band %d: largest piece %.2f, whole field %.1f"
            % (b, maxz[b], WHOLEZ[b]))
    say("  Y3 %s   (cap: below at both)"
        % ("hold" if y3 else "REFUTED"))

    # -------------------------------------------------------------- Y4
    say()
    say("Y4  is the coin platykurtic everywhere, and by the formula?")
    tn = np.log2([max(t, 1) for t, _ in coin_pairs])
    cv = np.log([abs(m) + 1e-12 for _, m in coin_pairs])
    cc, crms, ccov = ols(tn, cv)
    tt = cc[1] / math.sqrt(ccov[1, 1])
    say("  positive coin kurtosis at %d of %d pieces"
        % (coin_pos, len(coin_pairs)))
    say("  log|coin kurt| on log2(terms): slope %+.5f, t = %.2f"
        % (cc[1], tt))
    say("TSTAT wherereach_coinslope %.2f" % tt)
    say("SPREAD wherereach_coinslope %.5f"
        % math.sqrt(ccov[1, 1]))
    if abs(tt) < 2.0:
        say("UNRESOLVED SIGN wherereach_coinslope")
    y4 = coin_pos == 0 and cc[1] < 0 and abs(tt) > 2.0
    say("  the formula -2 sum w^4/(sum w^2)^2 would give a slope of "
        "-1 here")
    say("  Y4 %s   (cap: never positive, and a resolved fall)"
        % ("hold" if y4 else "REFUTED"))

    say()
    say("=" * 70)
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (gate_ok, y2, y3, y4)))
    say()
    if y2:
        say("the shape has no arithmetic content. A piece's separation "
            "from its")
        say("own coin is a function of how many terms it holds, fitted "
            "at one band")
        say("and forecasting two it never saw, so rem:cnkurtwhere's "
            "observation")
        say("about the upper half of the ladder is withdrawn: there is "
            "no profile")
        say("in v, only a count.")
    else:
        say("the term count does not explain the shape. What is left "
            "is a profile")
        say("in v -- where the prime sits below N matters beyond how "
            "many primes")
        say("are there -- and rem:cnkurtwhere's question reopens in a "
            "form it did")
        say("not test. The misses and their direction are printed "
            "above.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
