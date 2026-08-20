# -*- coding: utf-8 -*-
r"""
One coin draw is not an error bar: what audit_cn_law's controls can say

WHAT HAPPENED

audit_cn_law.py measured G(N) = C(N)/sqrt(V(N)) over even N in
(2e6, 4e6] and reported four things.  Two of its predictions broke:

  P2  sd(G) = 0.92953, outside the pre-registered [0.95, 1.05].
  P4  the coin arm failed to reproduce the real arm -- sd apart by
      8.81 blocked standard errors, excess kurtosis by 15.79.

Both readings rest on the same assumption, and it is wrong.  **The
coin arm was one draw.**  A single sign pattern eps is fixed and then
used at every N in the band, so the million values of G_coin share it;
the blocked standard error treats blocks as independent samples and
cannot see an offset common to all of them.  At a fixed N the identity
is exact -- averaging C_coin(N)^2 over 4000 fresh draws reproduces
V(N) to sampling error -- and the coin arm still came out at
sd = 0.91811, close to the real arm's 0.92953 rather than to 1.

So the question P2 and P4 were really asking cannot be answered by one
draw.  It can be answered by many: draw eps M times, recompute the
whole band each time, and use the draw-to-draw spread as the error bar
that the blocking could not supply.

This script does not revise audit_cn_law.py's predictions or its
recorded verdicts.  P2 and P4 are REFUTED and stay REFUTED.  What is
registered here is what those refutations mean.

BACKS: Remark {#rem:cnlaw} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  THE GATE.  At a fixed N, the mean of C_coin(N)^2 over M fresh
      draws reproduces V(N): the ratio is within 4/sqrt(M) of 1.  If
      the identity does not hold here the instrument is wrong and
      nothing else is reported.
  Q2  **A single draw is not an error bar.**  The draw-to-draw
      standard deviation of sd(G_coin) is at least ten times the
      blocked standard error audit_cn_law.py used for it (0.00083).
  Q3  And the deficit is not arithmetic: the real arm's
      sd(G) = 0.92953 lies inside the coin's draw-to-draw
      distribution, |z| < 2.  Then P2's refutation is a statement
      about the estimator -- V is the second moment in expectation
      over signs, not the mean square of any one sign pattern -- and
      not about mu.
  Q4  The kurtosis separation survives.  The real arm's excess
      kurtosis 0.26422 lies outside the coin's draw-to-draw
      distribution at |z| > 3, so P4's second half is a real
      separation and the first half is not.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED outside 4/sqrt(M).  Instrument failure; nothing below
      is reported.
  Q2  REFUTED below ten times.  Then the blocking was an adequate
      error bar after all, and P4's sd comparison stands as
      audit_cn_law.py read it.
  Q3  REFUTED at |z| >= 2.  Then the real field's mean square really
      does sit below what a sign pattern gives, and the deficit is
      arithmetic -- a stronger and more interesting outcome than the
      prediction, and the one that would make P2's refutation a
      finding about mu rather than about V.
  Q4  REFUTED at |z| <= 3.  Then the coin makes the kurtosis too and
      the bulk goes the way lem:coin sent the phase content, so
      **nothing in audit_cn_law.py separates mu from a coin** and the
      whole of P4's refutation was the single-draw artefact.

  Q3 and Q4 are opposite in what they hope for and that is
  deliberate: one of them says the sd result is an artefact, the other
  says the kurtosis result is not, and writing only the pair I expect
  would let either surprise pass unremarked.

  M is fixed at 64 before the run.  The draw-to-draw distribution is
  summarised by its mean and standard deviation and the z-scores are
  taken against those; with 64 draws a |z| of 3 is worth quoting and a
  |z| of 10 is not worth quoting precisely, so large z is reported as
  a bound rather than a number.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_cn_coin_spread.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_law.txt")

NMAX = 4_000_000
NLO = 2_000_000
DRAWS = 64
SEED = 20260820
FIXN = 3_141_592          # even, inside the band, for the gate
ZCAP = 2.0
ZKEEP = 3.0
RATIO = 10.0


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


def conv(a, b, nmax):
    L = 1 << (int(nmax + nmax).bit_length())
    return np.fft.irfft(np.fft.rfft(a, L) * np.fft.rfft(b, L),
                        L)[: nmax + 1]


def excess_kurt(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    return float((d ** 4).mean() / v ** 2 - 3.0)


def read_marker(name):
    """the value audit_cn_law.py printed, read not retyped"""
    src = io.open(SRC, encoding="utf-8").read()
    import re
    m = re.search(r"^%s\s+([-\d.]+)\s*$" % re.escape(name), src, re.M)
    if not m:
        raise SystemExit("no marker %s in %s" % (name, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the draw-to-draw distribution of sd(G_coin) and of",
    "           its excess kurtosis over %d independent sign" % DRAWS,
    "           patterns, each recomputed over the whole band, used",
    "           as the error bar that a single draw and a blocked",
    "           standard error could not supply; and the real arm's",
    "           two numbers scored against it.",
    "FIELD: even N in (%d, %d], the same band and the"
    % (NLO, NMAX),
    "       same Lambda and mu as audit_cn_law.py, whose sd and",
    "       excess kurtosis are read from",
    "       results/audit_cn_law.txt rather than retyped. The gate",
    "       is at N = %d." % FIXN,
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

    sd_real = read_marker("POINT cnlaw_sd")
    ku_real = read_marker("POINT cnlaw_kurt")
    sd_se = read_marker("SPREAD cnlaw_sd")
    say("READ audit_cn_law.txt POINT cnlaw_sd %.5f" % sd_real)
    say("READ audit_cn_law.txt POINT cnlaw_kurt %.5f" % ku_real)
    say("READ audit_cn_law.txt SPREAD cnlaw_sd %.5f" % sd_se)
    say("  the real arm's standard deviation, its fourth moment "
        "and the blocked")
    say("  standard error, read from that file rather than "
        "retyped")
    say("  real arm sd %.5f, excess kurtosis %.5f, blocked s.e. %.5f"
        % (sd_real, ku_real, sd_se))

    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    mu2 = muf ** 2
    V = conv(mu2, lam ** 2, NMAX)
    Ns = np.arange(NLO + 2, NMAX + 1, 2, dtype=np.int64)
    rootV = np.sqrt(V[Ns])
    rng = np.random.default_rng(SEED)

    sds, kus, fixed = [], [], []
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = conv(lam, eps * mu2, NMAX)
        g = Cc[Ns] / rootV
        sds.append(float(g.std(ddof=1)))
        kus.append(excess_kurt(g))
        fixed.append(float(Cc[FIXN]))
    sds = np.array(sds)
    kus = np.array(kus)
    fixed = np.array(fixed)

    # -------------------------------------------------------------- Q1
    say()
    say("Q1  does the identity hold at a fixed N over fresh draws?")
    ratio = float((fixed ** 2).mean() / V[FIXN])
    cap = 4.0 / math.sqrt(DRAWS)
    say("  N = %d: mean C_coin^2 over %d draws is %.4f of V(N)"
        % (FIXN, DRAWS, ratio))
    q1 = abs(ratio - 1.0) <= cap
    say("  Q1 %s   (cap: within %.4f of 1)"
        % ("hold" if q1 else "REFUTED", cap))
    if not q1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("  the coin's own spread over %d draws" % DRAWS)
    say("      statistic          mean        sd     min       max")
    say("      %-14s %9.5f %9.5f %8.5f %8.5f"
        % ("sd(G_coin)", sds.mean(), sds.std(ddof=1), sds.min(),
           sds.max()))
    say("      %-14s %+9.5f %9.5f %+8.5f %+8.5f"
        % ("excess kurt", kus.mean(), kus.std(ddof=1), kus.min(),
           kus.max()))
    say("SCALES 1")
    say("SPREAD coin_sd_draws %.5f" % sds.std(ddof=1))
    say("SPREAD coin_kurt_draws %.5f" % kus.std(ddof=1))

    # -------------------------------------------------------------- Q2
    say()
    say("Q2  is one draw an error bar?")
    r = sds.std(ddof=1) / sd_se
    say("  draw-to-draw sd %.5f against the blocked s.e. %.5f, a "
        "factor of %.1f" % (sds.std(ddof=1), sd_se, r))
    q2 = r >= RATIO
    say("  Q2 %s   (cap: at least %.0f times)"
        % ("hold" if q2 else "REFUTED", RATIO))

    # -------------------------------------------------------------- Q3
    say()
    say("Q3  is the deficit below 1 arithmetic?")
    z_sd = (sd_real - sds.mean()) / sds.std(ddof=1)
    say("  the coin's own sd averages %.5f, not 1: a fixed sign "
        "pattern does" % sds.mean())
    say("  not have mean square V over the band, only in expectation "
        "over patterns")
    say("  real %.5f scores z = %+.2f against the coin's draws"
        % (sd_real, z_sd))
    say("TSTAT cnspread_sd_z %.2f" % z_sd)
    say("SPREAD cnspread_sd_z %.5f" % sds.std(ddof=1))
    if abs(z_sd) < 2.0:
        say("UNRESOLVED SIGN cnspread_sd_z")
    q3 = abs(z_sd) < ZCAP
    say("  Q3 %s   (cap: |z| < %.0f)"
        % ("hold" if q3 else "REFUTED", ZCAP))

    # -------------------------------------------------------------- Q4
    say()
    say("Q4  does the kurtosis separation survive?")
    z_ku = (ku_real - kus.mean()) / kus.std(ddof=1)
    say("  the coin's excess kurtosis averages %+.5f over draws, "
        "spread %.5f" % (kus.mean(), kus.std(ddof=1)))
    say("  real %+.5f scores z = %+.1f" % (ku_real, z_ku))
    say("TSTAT cnspread_kurt_z %.2f" % z_ku)
    say("SPREAD cnspread_kurt_z %.5f" % kus.std(ddof=1))
    q4 = abs(z_ku) > ZKEEP
    say("  Q4 %s   (cap: |z| > %.0f)"
        % ("hold" if q4 else "REFUTED", ZKEEP))
    say("  with %d draws a z this size is a bound and not a "
        "measurement" % DRAWS)

    say()
    say("=" * 70)
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (q1, q2, q3, q4)))
    say()
    if q2 and q3:
        say("audit_cn_law.py's P2 and the sd half of its P4 are "
            "artefacts of a")
        say("single sign pattern. V is the second moment in "
            "expectation over")
        say("signs and no one pattern -- mu included -- has mean "
            "square exactly V")
        say("over a band. Both stay REFUTED as registered; what they "
            "refute is")
        say("the reading, not prop:V and not mu.")
    if q4:
        say("The kurtosis is different. No coin draw comes near the "
            "real arm's")
        say("value, so G's tails are not what a sign pattern on mu's "
            "support")
        say("produces. **That is one statistic of the withdrawn law "
            "that the")
        say("coin cannot make** -- the first thing in this branch to "
            "survive")
        say("lem:coin.")
    elif not q4:
        say("And the kurtosis goes the same way, so nothing measured "
            "here")
        say("separates mu from a coin and the withdrawn law is no "
            "better off")
        say("than when it was withdrawn.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
