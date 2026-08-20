# -*- coding: utf-8 -*-
r"""
The escape itself, against the multiplicative null

WHAT IS AT STAKE

rem:cnmultnull built the null this branch should have been using --
f(p) = +-1 iid, f(m) = product over p | m for squarefree m, so that
|f| = mu^2 and f is multiplicative exactly as mu is -- and then
compared the wrong statistic.  rem:cnclass's escape was on the second
moment, E[G^2] = 4.37497 at 105 | N against an iid coin ensemble at
1.02257, a z of 9.01 that no draw of 32 reached.  Its first moment was
already marginal there, holding by 3.20 against a coin ceiling of
3.17, and rem:cnclass said to read that as untested.
rem:cnmultnull compared first moments only.  So it re-ran the marginal
comparison against a second null and left the 9.01 untouched, and said
so.  This is the comparison that was missing.

The null is exact in the same way the iid one is.  For a random
multiplicative f, E[f(m)f(m')] vanishes unless m = m', because
f(m)f(m') is the product of f(p) over the symmetric difference of
their prime sets, so

    E_f[C_f(N)^2] = sum_m mu^2(m) Lambda(N-m)^2 = V(N)

exactly, and E[G_f^2] = 1 in every class at every N.  The null is
therefore "no structure at all", not "no difference between classes",
for this ensemble as much as for the coin.

BACKS: Remark {#rem:cnmultm2} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  THE GATE.  The real second moment at 105 | N reproduces
      rem:cnclass's POINT classm2_7 to three decimals.
  D2  **The multiplicative ensemble covers the escape.**  At least
      one of the 32 draws has its own largest |z| over classes,
      leave-one-out, reaching the real arm's.  This follows the
      reading rem:cnmultnull's C2 supports and is the deflationary
      outcome.
  D3  The two nulls are of similar strength here, as
      rem:cnmultnull's C3 found on the first moment: the ensemble's
      spread in E[G^2] at 105 | N is within a factor of 1.5 either
      way of the iid coin's 0.3721.
  D4  The ensemble is centred where the algebra says: each class's
      ensemble mean E[G^2] is within two of its own standard errors
      of 1.

REFUTATION RULE (fixed before the run)

  D1  REFUTED outside three decimals; nothing below is reported.
  D2  **REFUTED if no draw reaches the real arm.**  Then mu is
      outside the strongest null this branch can build, on the
      statistic where the escape actually lives.  That is a large
      claim and it would need repeating at another band before
      anything is built on it; the remark must say so and must not
      treat one band as settling it.
  D3  REFUTED outside 1.5 either way.  A much wider ensemble would
      mean D2's outcome says less than it appears to; a much
      narrower one would mean the multiplicative structure makes the
      null easier to escape, which would be worth more than D2
      itself.
  D4  REFUTED if any class's ensemble mean is two standard errors
      from 1.  **The unresolved case is named**: eight classes at a
      two-standard-error threshold fire by chance about one time in
      three, so a single crossing just past the threshold is the
      too-noisy-to-tell case and must be read as "not resolved",
      never as "the ensemble is biased".  Only a crossing large
      enough to survive eight looks, or one repeated across classes,
      would carry the reading D4 was written for.

  WHAT THIS CANNOT DO.  One band, one class family, 32 draws of a
  random multiplicative function -- an object with heavy tails, so 32
  draws estimate its spread poorly.  A D2 that holds by one draw is
  not a D2 that holds by twenty and the count is reported either way.
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
OUT = os.path.join(ROOT, "results", "audit_cn_multnull_m2.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_class.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
QS = (3, 5, 7)
FULL = 7
DRAWS = 32
SEED = 20260901
DEC = 3
COINSD = 0.3721            # the iid coin's spread, from rem:cnclass
FACT = 1.5


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


def sign_bits(n, primes, rng, draws):
    """parity of the chosen prime-signs dividing each m, all draws at once"""
    t = np.zeros(n + 1, dtype=np.uint32)
    for p in primes:
        p = int(p)
        t[p::p] ^= np.uint32(rng.integers(0, 1 << draws,
                                          dtype=np.uint64))
    return t


def read_mark(path, name):
    src = io.open(path, encoding="utf-8").read()
    m = re.search(r"^%s ([-\d.]+)\s*$" % re.escape(name), src, re.M)
    if not m:
        raise SystemExit("no %s in %s" % (name, path))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the second moments E[G^2] of",
    "           G(N) = C(N)/sqrt(V(N)) over the eight classes cut by",
    "           which of %s divide N, scored against an ensemble of" % (QS,),
    "           %d random multiplicative sign functions by the" % DRAWS,
    "           largest |z| over classes, each draw taken",
    "           leave-one-out against the rest.",
    "FIELD: even N in (2^%d, 2^%d], with Lambda and mu sieved once to"
    % (BAND, BAND + 1),
    "       %d and every convolution taken by FFT over the whole" % NMAX,
    "       range. For a random multiplicative f the identity",
    "       E[C_f(N)^2] = V(N) is exact, so the ensemble's E[G^2] is",
    "       1 in every class at every N and the null is no structure",
    "       at all.",
    "SEED: the draws come from numpy default_rng at seed %d; without"
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

    pub = read_mark(SRC, "POINT classm2_%d" % FULL)
    say("READ audit_cn_class.txt POINT classm2_%d %.5f" % (FULL, pub))
    say("  the class second moment this run has to reproduce")
    say("PRINTBOUND audit_cn_multnull_m2 %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d random multiplicative draws; the iid coin's spread at "
        "this class" % DRAWS)
    say("  was %.4f" % COINSD)

    say("sieving to %d, which is 2^%d" % (NMAX, BAND + 1))
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

    Ns = np.arange((1 << BAND) + 2, (1 << (BAND + 1)) + 1, 2,
                   dtype=np.int64)
    rootV = np.sqrt(np.maximum(V[Ns], 1e-300))
    del V
    lab = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        lab |= ((Ns % q) == 0).astype(np.int64) << i
    masks = [lab == k for k in range(1 << len(QS))]
    sizes = [int(m.sum()) for m in masks]

    def cls_m2(arr):
        g = arr[Ns] / rootV
        return np.array([float((g[m] ** 2).mean()) for m in masks])

    def cls_m1(arr):
        g = arr[Ns] / rootV
        return np.array([float(g[m].mean()) for m in masks])

    real = cls_m2(C)
    real1 = cls_m1(C)
    del C

    # -------------------------------------------------------------- D1
    say()
    say("D1  does this run reproduce the class it shares?")
    say("  class %d E[G^2] here %.5f against its %.5f"
        % (FULL, real[FULL], pub))
    d1 = abs(round(real[FULL], DEC) - round(pub, DEC)) \
        < 10.0 ** (-DEC) / 2
    say("  D1 %s   (cap: %d decimals)"
        % ("hold" if d1 else "REFUTED", DEC))
    if not d1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    pr = primes_upto(NMAX)
    say()
    say("  building %d multiplicative sign functions in one pass over "
        "%d primes" % (DRAWS, len(pr)))
    t = sign_bits(NMAX, pr, rng, DRAWS)
    del pr

    pool = np.zeros((DRAWS, len(masks)))
    pool1 = np.zeros((DRAWS, len(masks)))
    for d in range(DRAWS):
        f = mu2 * (1.0 - 2.0 * ((t >> np.uint32(d)) & np.uint32(1)))
        Cf = np.fft.irfft(FL * np.fft.rfft(f, L), L)[: NMAX + 1]
        pool[d] = cls_m2(Cf)
        pool1[d] = cls_m1(Cf)
        del Cf, f

    def zsc(x, p):
        return (x - p.mean(axis=0)) / p.std(axis=0, ddof=1)

    zr = zsc(real, pool)
    maxr = float(np.abs(zr).max())
    where = int(np.argmax(np.abs(zr)))
    mc = []
    for d in range(DRAWS):
        rest = np.delete(pool, d, axis=0)
        mc.append(float(np.abs(zsc(pool[d], rest)).max()))
    mc = np.array(mc)

    say()
    say("      3|N 5|N 7|N        N      E[G^2]    ens. mean"
        "      ens. sd       z     mean G   ens. m1 sd")
    for k in range(len(masks)):
        say("       %d   %d   %d  %8d  %+10.5f  %+11.5f  %11.5f  "
            "%6.2f  %+9.5f  %10.5f"
            % (k & 1, (k >> 1) & 1, (k >> 2) & 1, sizes[k],
               real[k], pool[:, k].mean(), pool[:, k].std(ddof=1),
               zr[k], real1[k], pool1[:, k].std(ddof=1)))
        say("POINT multm2_z_%d %.5f" % (k, zr[k]))
    say("SCALES 1")

    # -------------------------------------------------------------- D2
    say()
    say("D2  does the ensemble cover the escape?")
    reach = int((mc >= maxr).sum())
    say("  real largest |z| %.2f at class %d; the draws' own largest "
        "run %.2f to %.2f," % (maxr, where, mc.min(), mc.max()))
    say("  and %d of %d reach it" % (reach, DRAWS))
    d2 = reach > 0
    say("TSTAT multm2_maxz %.2f" % maxr)
    say("SPREAD multm2_maxz %.5f" % float(mc.std(ddof=1)))
    say("  D2 %s   (cap: at least one draw)"
        % ("hold" if d2 else "REFUTED"))

    # -------------------------------------------------------------- D3
    say()
    say("D3  are the two nulls of similar strength here?")
    sd = float(pool[:, FULL].std(ddof=1))
    ratio = sd / COINSD
    say("  ensemble spread at class %d is %.5f against the iid coin's "
        "%.4f," % (FULL, sd, COINSD))
    say("  a factor of %.2f" % ratio)
    d3 = (1.0 / FACT) <= ratio <= FACT
    say("POINT multm2_spreadratio %.5f" % ratio)
    say("  D3 %s   (cap: within %.1f either way)"
        % ("hold" if d3 else "REFUTED", FACT))

    # -------------------------------------------------------------- D4
    say()
    say("D4  is the ensemble centred at one?")
    worst, wk = 0.0, None
    for k in range(len(masks)):
        se = float(pool[:, k].std(ddof=1)) / math.sqrt(DRAWS)
        r = abs(float(pool[:, k].mean()) - 1.0) / (se + 1e-300)
        if r > worst:
            worst, wk = r, k
    say("  furthest from one at class %s, at %.2f of its own standard "
        "error" % (wk, worst))
    d4 = worst <= 2.0
    say("  D4 %s   (cap: two standard errors)"
        % ("hold" if d4 else "REFUTED"))

    say()
    say("=" * 70)
    say("D1 %s  D2 %s  D3 %s  D4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (d1, d2, d3, d4)))
    say()
    if d2:
        say("the escape does not survive the stronger null. A random "
            "multiplicative")
        say("sign function reaches where mu goes, so what rem:cnclass "
            "measured is")
        say("what multiplicativity plus the constraint at 105 | N can "
            "produce, and")
        say("this branch has no escape left: seven routes, and the "
            "seventh needed")
        say("a control weaker than the object.")
    else:
        say("mu is outside the strongest null this branch can build, "
            "on the")
        say("statistic where the escape lives. No draw of a random "
            "multiplicative")
        say("sign function -- same support, same multiplicativity, "
            "same constraint")
        say("at 105 | N -- reaches it. That is a large claim on one "
            "band and it")
        say("needs repeating at another before anything is built on "
            "it.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
