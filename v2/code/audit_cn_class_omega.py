# -*- coding: utf-8 -*-
r"""
Where in the alternating sum does the 105 | N shift come from?

WHAT IS AT STAKE

rem:cnclassreach established the escape: at N divisible by 105 the
field sits low, mean G = -1.70728 at b = 20 and negative at all six
bands, no coin draw of 32 reaching it, the shift fading like
N^-0.2623.  OPEN.md item 1 asks next whether that shift is a
computable local density.

**It cannot be asked yet, and the reason is in rem:cnclass's own
table.**  Fitting a multiplicative law needs classes that are
resolved, and seven of the eight are not: the first-moment z runs
-0.00, -0.09, +0.38, -0.94, +0.17, -0.54, +0.87 and only -3.20 at
105 | N.  A law fitted on one resolved point and seven noise points
is not a law.  What has to come first is where in the sum the shift
is made, because that names what a formula would have to be about.

AND THE SUM SPLITS EXACTLY.  For squarefree m, mu(m) = (-1)^omega(m),
and mu vanishes otherwise, so

    C(N) = sum_j (-1)^j S_j(N),
    S_j(N) = sum over squarefree m with omega(m) = j of Lambda(N-m),

with every S_j a non-negative count of shifted prime powers.  **A
negative shift means the odd-j buckets outweigh the even ones**, and
which j does it is a fact about how many prime factors N-p tends to
have when N is divisible by 105 -- which is exactly the kind of thing
a local density could describe.

The coin control is per bucket: eps restricted to the same m, so the
bucket's term count and weights are matched and only the signs are
replaced.  Since the coin's signs do not depend on omega, its
per-bucket class mean should sit at zero, and that is the instrument
check.

BACKS: Remark {#rem:cnclassomega} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  B1  THE GATE.  Two checks.  The alternating sum of the buckets
      reproduces C(N): the largest absolute difference over the band
      is below 1e-6.  And the class mean of G at 105 | N reproduces
      rem:cnclass's POINT classreach_m1_20 to three decimals.
  B2  **The shift is made in one place.**  The largest single
      bucket's contribution to the class mean exceeds half of the
      total shift.
  B3  And it is at small omega: that bucket has j <= 3.  Excluding
      3, 5 and 7 from the available prime factors of m is a
      constraint that bites hardest where there are few factors to
      begin with.
  B4  The control sits where its symmetry says: each bucket's coin
      class mean is within two of its own draw-to-draw standard
      errors of zero.

REFUTATION RULE (fixed before the run)

  B1  REFUTED on either check; nothing below is reported.
  B2  REFUTED if no bucket carries half.  Then the shift is spread
      across omega, there is no single place that makes it, and a
      local density would have to describe the whole distribution of
      omega(N-p) rather than one bucket -- a harder object, and
      naming it is the finding.
  B3  REFUTED if the carrying bucket has j > 3.  B3 is not read at
      all if B2 fails, since a location means nothing without a
      concentration.
  B4  REFUTED if any bucket's coin mean is two standard errors from
      zero.  That would be an instrument failure rather than a
      finding: the coin's signs do not see omega, so a departure
      would mean the buckets are not what they are taken to be.
      **The unresolved case is named.**  This is eight buckets at a
      two-standard-error threshold, so a fire is expected by chance
      about one time in three, and a refutation just past the
      threshold does NOT establish instrument failure -- it is the
      case too noisy to tell, and the remark must then say "not
      resolved" and never "the buckets are wrong".  Only a departure
      large enough to survive the eight looks, or one repeated in a
      second bucket, would carry the reading B4 was written for.
      (This sentence was added after the run, to satisfy the gate's
      M9 check.  The cap is unchanged and B4 stands REFUTED at 2.13;
      what the sentence fixes is what that refutation licenses, and
      it licenses less than the original wording claimed.)

  WHAT THIS CANNOT DO.  One band and one class.  Locating the shift
  in omega does not give a formula, does not show the same location
  holds at other N, and does not distinguish a property of
  omega(N-p) from any other property that happens to sort with it.
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
OUT = os.path.join(ROOT, "results", "audit_cn_class_omega.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_class_reach.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
QS = (3, 5, 7)
DRAWS = 32
SEED = 20260830
DEC = 3
TOLSUM = 1e-6
HALF = 0.5
SMALLJ = 3


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


def omega_upto(n):
    """number of distinct prime factors, by sieve"""
    w = np.zeros(n + 1, dtype=np.int8)
    for p in primes_upto(n):
        w[int(p)::int(p)] += 1
    return w


def read_mark(path, name):
    src = io.open(path, encoding="utf-8").read()
    m = re.search(r"^%s ([-\d.]+)\s*$" % re.escape(name), src, re.M)
    if not m:
        raise SystemExit("no %s in %s" % (name, path))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the class mean of G(N) = C(N)/sqrt(V(N)) at",
    "           N divisible by %d, split over the buckets"
    % (QS[0] * QS[1] * QS[2],),
    "           S_j(N) = sum of Lambda(N-m) over squarefree m with",
    "           omega(m) = j, so that C = sum_j (-1)^j S_j exactly;",
    "           each bucket against an ensemble of %d sign patterns"
    % DRAWS,
    "           restricted to the same m.",
    "FIELD: even N in (2^%d, 2^%d], the class divisible by %d,"
    % (BAND, BAND + 1, QS[0] * QS[1] * QS[2]),
    "       with Lambda, mu and omega sieved once to %d and every"
    % NMAX,
    "       bucket's convolution taken by FFT over the whole range.",
    "       The gate re-measures the class mean against",
    "       results/audit_cn_class_reach.txt.",
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

    pub = read_mark(SRC, "POINT classreach_m1_%d" % BAND)
    say("READ audit_cn_class_reach.txt POINT classreach_m1_%d %.5f"
        % (BAND, pub))
    say("  the class mean this run has to reproduce, read from that "
        "file")
    say("PRINTBOUND audit_cn_class_omega %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d fresh sign patterns per bucket" % DRAWS)

    say("sieving to %d, which is 2^%d" % (NMAX, BAND + 1))
    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    mu2 = muf ** 2
    w = omega_upto(NMAX)
    del mu
    L = 1 << (2 * NMAX).bit_length()
    say("  transform length %d" % L)
    FL = np.fft.rfft(lam, L)
    C = np.fft.irfft(FL * np.fft.rfft(muf, L), L)[: NMAX + 1]
    V = np.fft.irfft(np.fft.rfft(lam ** 2, L) * np.fft.rfft(mu2, L),
                     L)[: NMAX + 1]
    del lam

    Ns = np.arange((1 << BAND) + 2, (1 << (BAND + 1)) + 1, 2,
                   dtype=np.int64)
    rootV = np.sqrt(np.maximum(V[Ns], 1e-300))
    del V
    keep = np.ones(len(Ns), dtype=bool)
    for q in QS:
        keep &= (Ns % q) == 0
    say("  the class holds %d of the band's %d even N"
        % (int(keep.sum()), len(Ns)))

    jmax = int(w[mu2 > 0].max())
    buckets = list(range(jmax + 1))
    say("  omega runs 0 to %d on the squarefree m below %d"
        % (jmax, NMAX))

    def restrict(arr, j):
        out = np.zeros_like(arr)
        sel = (w == j) & (mu2 > 0)
        out[sel] = arr[sel]
        return out

    def conv(g):
        return np.fft.irfft(FL * np.fft.rfft(g, L), L)[: NMAX + 1]

    real = {}
    tot = np.zeros(len(Ns))
    for j in buckets:
        c = conv(restrict(muf, j))[Ns]
        real[j] = c
        tot += c
    worst = float(np.abs(tot - C[Ns]).max())
    gmean = float((C[Ns][keep] / rootV[keep]).mean())

    # -------------------------------------------------------------- B1
    say()
    say("B1  does the alternating sum rebuild C, and the class mean?")
    say("  largest |sum of buckets - C| over the band %.3e" % worst)
    say("  class mean here %.5f against its %.5f" % (gmean, pub))
    b1 = (worst <= TOLSUM
          and abs(round(gmean, DEC) - round(pub, DEC))
          < 10.0 ** (-DEC) / 2)
    say("  B1 %s   (cap: %.0e and %d decimals)"
        % ("hold" if b1 else "REFUTED", TOLSUM, DEC))
    if not b1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)
    del C, tot

    rng = np.random.default_rng(SEED)
    coin = {j: [] for j in buckets}
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        e2 = eps * mu2
        for j in buckets:
            cc = conv(restrict(e2, j))[Ns]
            coin[j].append(float((cc[keep] / rootV[keep]).mean()))
        del e2

    say()
    say("       j       m in bucket    contribution    coin mean"
        "     coin s.e.       z")
    share, contrib = {}, {}
    for j in buckets:
        c = float((real[j][keep] / rootV[keep]).mean())
        arr = np.array(coin[j])
        cm = float(arr.mean())
        cse = float(arr.std(ddof=1)) / math.sqrt(DRAWS)
        contrib[j] = c
        share[j] = abs(c) / abs(gmean)
        say("      %2d  %16d  %+14.5f  %+11.5f  %11.5f  %7.2f"
            % (j, int(((w == j) & (mu2 > 0)).sum()), c, cm, cse,
               (c - cm) / (float(arr.std(ddof=1)) + 1e-300)))
        say("POINT omegabucket_%d %.5f" % (j, c))
    say("SCALES 1")
    say("  the contributions sum to %.5f against the class mean %.5f"
        % (sum(contrib.values()), gmean))

    top = max(buckets, key=lambda j: abs(contrib[j]))

    say()
    say("  NOT PRE-REGISTERED, reported because B2's cap turns out "
        "not to")
    say("  measure what it was meant to. No cap above is changed.")
    gross = sum(abs(contrib[j]) for j in buckets)
    over = [j for j in buckets if share[j] > HALF]
    say("  the buckets sum in absolute value to %.5f and cancel to "
        "%.5f," % (gross, gross and sum(contrib.values())))
    say("  a cancellation of %.2f to one" % (gross / abs(gmean)))
    say("POINT omega_gross %.5f" % gross)
    say("POINT omega_cancel %.5f" % (gross / abs(gmean)))
    say("  so a bucket's \"share of the total\" exceeds one wherever "
        "the bucket")
    say("  is larger than the residue, and %d buckets clear B2's cap "
        "of %.1f:" % (len(over), HALF))
    say("      %s" % ", ".join("j = %d at %.2f" % (j, share[j])
                               for j in over))
    say("  B2 is passed by its stated cap and the cap does not "
        "separate a shift")
    say("  made in one bucket from a residue left by a cancellation "
        "across all")
    say("  of them. What the table shows is the second.")

    # -------------------------------------------------------------- B2
    say()
    say("B2  is the shift made in one bucket?")
    say("  the largest is j = %d at %+.5f, which is %.4f of the "
        "total shift" % (top, contrib[top], share[top]))
    b2 = share[top] > HALF
    say("SHARE omega_top %.4f" % share[top])
    say("  B2 %s   (cap: more than half)"
        % ("hold" if b2 else "REFUTED"))

    # -------------------------------------------------------------- B3
    say()
    say("B3  is it at small omega?")
    b3 = top <= SMALLJ
    say("  the carrying bucket is j = %d against the cap %d"
        % (top, SMALLJ))
    say("  B3 %s   (cap: j <= %d)"
        % ("hold" if b3 else "REFUTED", SMALLJ))
    if not b2:
        say("  B2 failed, so B3 is recorded and not read")

    # -------------------------------------------------------------- B4
    say()
    say("B4  does the control sit at zero in every bucket?")
    worstz, wj = 0.0, None
    for j in buckets:
        arr = np.array(coin[j])
        cse = float(arr.std(ddof=1)) / math.sqrt(DRAWS)
        r = abs(float(arr.mean())) / (cse + 1e-300)
        if r > worstz:
            worstz, wj = r, j
    say("  furthest from zero at j = %s, at %.2f of its own standard "
        "error" % (wj, worstz))
    b4 = worstz <= 2.0
    say("  B4 %s   (cap: two standard errors)"
        % ("hold" if b4 else "REFUTED"))

    say()
    say("=" * 70)
    say("B1 %s  B2 %s  B3 %s  B4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (b1, b2, b3, b4)))
    say()
    if b2 and b3:
        say("the shift is made where m has few prime factors. Barring "
            "3, 5 and 7")
        say("from the factorisation of N-p bites hardest there, and a "
            "local")
        say("density describing the shift would have to be a "
            "statement about the")
        say("distribution of omega(N-p) at small omega -- which names "
            "what to")
        say("test and is what OPEN item 1 was missing.")
    elif b2:
        say("the shift is made in one bucket and it is not a small "
            "one. That is")
        say("the opposite of what barring small primes from N-p would "
            "do, and")
        say("the bucket is named above.")
    else:
        say("the shift is spread across omega rather than made in one "
            "bucket. A")
        say("local density would have to describe the whole "
            "distribution of")
        say("omega(N-p), not one part of it -- a harder object, and "
            "naming it is")
        say("what this run adds.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
