# -*- coding: utf-8 -*-
r"""
A better null: is the escape mu's, or is it multiplicativity's?

WHAT IS AT STAKE

Every control in this branch has been the same one -- eps(m) drawn
iid on mu's support -- and it is a weak null for mu.  mu is
multiplicative and eps is not.  A field built from independent signs
at every m has no reason to show the class structure
rem:cnclass found, and it does not; that is most of why the escape
looked like an escape.

The natural null keeps the structure and changes only the values.
Let f(p) = +-1 iid over primes and f(m) = product of f(p) over p | m
for squarefree m, zero otherwise.  Then |f| = mu^2 exactly, f is
multiplicative exactly as mu is, and E[f(m)] = 0 for every m > 1.
**Anything the escape owes to multiplicativity plus the coprimality
constraint at 105 | N, this ensemble has too.**  What it does not have
is mu's particular choice of sign at each prime.

rem:cnclassomega makes the question urgent rather than academic.  The
class shift is the residue of a 197.70-to-one cancellation between
omega-buckets, and what survives such a cancellation is exactly the
kind of thing that depends on the multiplicative structure of the
signs rather than on their individual values.

The comparison is the max-statistic device used throughout: the real
arm's largest |z| over the eight classes against each draw's own
largest, taken leave-one-out.  The iid coin's numbers are read from
rem:cnclass so the two nulls can be set beside each other.

BACKS: Remark {#rem:cnmultnull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  THE GATE.  The real class mean at 105 | N reproduces
      rem:cnclassreach's POINT classreach_m1_20 to three decimals.
  C2  **The multiplicative ensemble reproduces the escape.**  At
      least one of the 32 draws has its own largest |z| over classes
      reaching the real arm's, so mu is not outside this null.  This
      is the deflationary prediction and it is the one this branch's
      history points at: six routes closed and the seventh escaped
      only a null that was weaker than the object.
  C3  The ensemble is visibly the stronger null: its draw-to-draw
      spread in the class mean at 105 | N is at least three times the
      iid coin's, which rem:cnclass's numbers put at 0.5646.
  C4  It is still centred: each class's ensemble mean is within two
      of its own standard errors of zero, as E[f(m)] = 0 requires.

REFUTATION RULE (fixed before the run)

  C1  REFUTED outside three decimals; nothing below is reported.
  C2  REFUTED if no draw reaches the real arm.  **Then the escape is
      mu's own and not multiplicativity's** -- it would survive the
      strongest null this branch can build, and that is a far larger
      claim than rem:cnclass made.  It is the outcome to state
      plainly and carefully if it comes, and the first thing to do
      with it would be to repeat it at another band before believing
      it.
  C3  REFUTED below three times.  Then the two nulls are of similar
      strength and C2's outcome, either way, says less than it
      appears to: a null that does not spread more has not been made
      harder to escape.
  C4  REFUTED if any class's ensemble mean is two standard errors
      from zero.  **The unresolved case is named**: eight classes at a
      two-standard-error threshold fire by chance about one time in
      three, so a single crossing just past the threshold is the
      too-noisy-to-tell case and must be read as "not resolved", not
      as "the ensemble is biased".  Only a crossing large enough to
      survive eight looks would carry that reading.

  WHAT THIS CANNOT DO.  One band, one class family, and 32 draws of a
  random multiplicative function -- an object whose own distribution
  is famously heavy-tailed, so 32 draws estimate its spread poorly.
  A C2 that holds by one draw out of 32 is not the same as one that
  holds by twenty, and the count is reported.
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
OUT = os.path.join(ROOT, "results", "audit_cn_multnull.txt")
SRCR = os.path.join(ROOT, "results", "audit_cn_class_reach.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
QS = (3, 5, 7)
FULL = 7
DRAWS = 32
SEED = 20260831
DEC = 3
COINSD = 0.5646            # the iid coin's spread, from rem:cnclass
TIMES = 3.0


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
    """parity of the chosen prime-signs dividing each m, all draws at once

    One pass over the primes carries every draw: prime p gets a
    draws-bit mask and the parity is accumulated by xor, so the cost
    is one slice per prime rather than one per prime per draw.
    """
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
    "STATISTIC: the class means of G(N) = C(N)/sqrt(V(N)) over the",
    "           eight classes cut by which of %s divide N, scored" % (QS,),
    "           against an ensemble of %d random multiplicative sign" % DRAWS,
    "           functions -- f(p) = +-1 iid, f(m) = product over",
    "           p | m for squarefree m, zero otherwise -- by the",
    "           largest |z| over classes with each draw taken",
    "           leave-one-out against the rest.",
    "FIELD: even N in (2^%d, 2^%d], with Lambda and mu sieved once to"
    % (BAND, BAND + 1),
    "       %d and every convolution taken by FFT over the whole" % NMAX,
    "       range. |f| = mu^2 exactly, so the null matches the",
    "       support and the multiplicativity and differs only in",
    "       which sign each prime carries.",
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

    pub = read_mark(SRCR, "POINT classreach_m1_%d" % BAND)
    say("READ audit_cn_class_reach.txt POINT classreach_m1_%d %.5f"
        % (BAND, pub))
    say("  the class mean this run has to reproduce, read from that "
        "file")
    say("PRINTBOUND audit_cn_multnull %d %.8f"
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

    def cls_mean(arr):
        g = arr[Ns] / rootV
        return np.array([float(g[m].mean()) for m in masks])

    real = cls_mean(C)
    del C

    # -------------------------------------------------------------- C1
    say()
    say("C1  does this run reproduce the class it shares?")
    say("  class %d mean here %.5f against its %.5f"
        % (FULL, real[FULL], pub))
    c1 = abs(round(real[FULL], DEC) - round(pub, DEC)) \
        < 10.0 ** (-DEC) / 2
    say("  C1 %s   (cap: %d decimals)"
        % ("hold" if c1 else "REFUTED", DEC))
    if not c1:
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
    for d in range(DRAWS):
        f = mu2 * (1.0 - 2.0 * ((t >> np.uint32(d)) & np.uint32(1)))
        Cf = np.fft.irfft(FL * np.fft.rfft(f, L), L)[: NMAX + 1]
        pool[d] = cls_mean(Cf)
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
    say("      3|N 5|N 7|N        N      mean G    ens. mean"
        "      ens. sd       z")
    for k in range(len(masks)):
        say("       %d   %d   %d  %8d  %+10.5f  %+11.5f  %11.5f  "
            "%6.2f" % (k & 1, (k >> 1) & 1, (k >> 2) & 1, sizes[k],
                       real[k], pool[:, k].mean(),
                       pool[:, k].std(ddof=1), zr[k]))
        say("POINT multnull_z_%d %.5f" % (k, zr[k]))
    say("SCALES 1")

    # -------------------------------------------------------------- C2
    say()
    say("C2  is mu inside the multiplicative ensemble?")
    reach = int((mc >= maxr).sum())
    say("  real largest |z| %.2f at class %d; the draws' own largest "
        "run %.2f to %.2f," % (maxr, where, mc.min(), mc.max()))
    say("  and %d of %d reach it" % (reach, DRAWS))
    c2 = reach > 0
    say("TSTAT multnull_maxz %.2f" % maxr)
    say("SPREAD multnull_maxz %.5f" % float(mc.std(ddof=1)))
    say("  C2 %s   (cap: at least one draw)"
        % ("hold" if c2 else "REFUTED"))

    # -------------------------------------------------------------- C3
    say()
    say("C3  is this null the stronger one?")
    sd = float(pool[:, FULL].std(ddof=1))
    say("  ensemble spread at class %d is %.5f against the iid coin's "
        "%.4f," % (FULL, sd, COINSD))
    say("  a factor of %.2f" % (sd / COINSD))
    c3 = sd >= TIMES * COINSD
    say("POINT multnull_spreadratio %.5f" % (sd / COINSD))
    say("  C3 %s   (cap: %.0f times)"
        % ("hold" if c3 else "REFUTED", TIMES))

    # -------------------------------------------------------------- C4
    say()
    say("C4  is the ensemble centred?")
    worst, wk = 0.0, None
    for k in range(len(masks)):
        se = float(pool[:, k].std(ddof=1)) / math.sqrt(DRAWS)
        r = abs(float(pool[:, k].mean())) / (se + 1e-300)
        if r > worst:
            worst, wk = r, k
    say("  furthest from zero at class %s, at %.2f of its own "
        "standard error" % (wk, worst))
    c4 = worst <= 2.0
    say("  C4 %s   (cap: two standard errors)"
        % ("hold" if c4 else "REFUTED"))

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  C3 %s  C4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, c3, c4)))
    say()
    if c2:
        say("the escape was the null's weakness, not mu's strength. "
            "A random")
        say("multiplicative sign function shows the same class "
            "structure, so what")
        say("rem:cnclass measured is what multiplicativity plus the "
            "coprimality")
        say("constraint at 105 | N produce, and not something "
            "particular to mu.")
        say("Six routes closed and the seventh escaped only a control "
            "weaker than")
        say("the object it was controlling.")
    else:
        say("mu is outside the strongest null this branch can build. "
            "A random")
        say("multiplicative sign function has mu's support, mu's "
            "multiplicativity")
        say("and mu's coprimality constraint at 105 | N, and none of "
            "the draws")
        say("reaches the field. That is a much larger claim than "
            "rem:cnclass made")
        say("and the first thing to do with it is repeat it at "
            "another band.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
