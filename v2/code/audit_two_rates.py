# -*- coding: utf-8 -*-
r"""
The two rates side by side, and the one number between them

WHAT IS AT STAKE

rem:maintermremoval found the structure of item 5: A = I(1) is a main
term of order N log N, B = sum_{d>=2} mu(d) I(d) cancels it, and the
cancellation improves -- |B|/|A| runs 0.874261 to 0.976934 across the
field.  And yet |S|/l2 grows, 5.0350 to 9.5118.  Both are true because
two rates are racing: B closes on A at one rate, l2 shrinks relative
to A at another, and the demand is that the first beat the second.

The arithmetic fixes what to measure.  With r = |A+B|/|A| the relative
residue, |S| = |A| r exactly, so

    e(|S|) = e(|A|) + e(r)

and the demand |S| <~ l2 is

    e(r)  <=  e(l2) - e(|A|)          the required residue exponent

against the achieved e(r) = -0.342848.  **The difference of those two
is the single number rem:maintermremoval asked for**, and it is
identically e(|S|) - e(l2), which the field publishes as the deficit
+0.134019 -- so the run carries its own cross-check rather than
needing to be believed.

Nothing new is computed at large N: the same eight N as
rem:maintermremoval, over 2.1 decades.

BACKS: Remark {#rem:tworates} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  AA1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker, and l2 there reproduces rem:jbarrier's, both
      to a relative 1e-12.
  AA2 THE IDENTITY.  e(|S|) equals e(|A|) + e(r) to within 0.001 --
      it must, and a failure would mean the three fits are not of the
      quantities they are labelled with.
  AA3 **THE NUMBER.**  The gap between the required residue exponent
      e(l2) - e(|A|) and the achieved e(r) is positive and lands
      within 0.05 of the published deficit +0.134019.
  AA4 **And in plain terms at the top N.**  To meet the demand there,
      1 - |B|/|A| would have to be smaller than it is by the factor
      |S|/l2, which puts the required cancellation below 0.005 --
      that is, exact to under half a per cent where it is now exact
      to 2.31 per cent.

REFUTATION RULE (fixed before the run)

  AA1 REFUTED outside 1e-12 on either; nothing below is reported.
  AA2 REFUTED outside 0.001.  Then the labels are wrong and no
      reading follows.
  AA3 **REFUTED if the gap is negative, or outside 0.05 of
      +0.134019.**  A negative gap would say the residue already
      shrinks fast enough on this field, which would contradict
      |S|/l2 growing and put the disagreement inside this run.  A gap
      far from the published deficit would say this eight-point field
      and the published one are measuring different things, and
      neither could then be used without reconciling them.
  AA4 REFUTED at or above 0.005.

  A SECOND BLOCK, REGISTERED AFTER AA1 WAS REFUTED AND SAYING SO

  AA1 is refuted above and stays refuted.  Its docstring said the
  marker is read at full double precision "because ten digits was not
  enough last time" -- and that fixed the reader while leaving the
  **emitters** at seven significant figures.  POINT mainA carried
  %.6e and POINT jbarrier carried %.6e, so the relatives came out
  1.17e-07 and 5.33e-10 against a 1e-12 cap: TOL BELOW PRINT for the
  third tick running, and the first time it was anticipated in words
  and still missed in fact.

  audit_mainterm_removal.py is this repository's own and has exactly
  one consumer -- this script -- so its marker is widened to full
  double precision and it is re-run, with no timestamp cascade.
  audit_jbarrier.py cannot be treated that way: re-running it makes
  audit_jbarrier_reach.py and audit_which_floor.py older than what
  they read, and G22 would then require re-running those too.  **A
  print-width fix is not worth invalidating two measured runs**, so
  the l2 leg is judged at its own marker's print bound instead, which
  is what PRINTBOUND exists to express.

  So, pre-registered before the second run and after seeing the first:

  BB1 A at N = 200000 reproduces the widened POINT mainA marker to a
      relative 1e-12, and l2 there reproduces POINT jbarrier to that
      marker's own print bound of 5e-7 relative, which this run
      prints.  **The tolerance on the l2 leg is loosened to the
      source's precision and nowhere else**; the A leg keeps 1e-12.

  REFUTATION for the second block.  BB1 refuted outside either
  tolerance, and then the split computed here is not the split those
  runs published and nothing below may be read.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  Three exponents are
  fitted on eight points and their errors are printed.  The gap in
  AA3 is a difference of two of them fitted on the *same* points, so
  its error is the error of a fit to the ratio and not the sum of
  theirs -- **this run fits log(|S|/l2) directly and prints that
  error**, which is the mistake rem:jbarrier's K1 made and which is
  not repeated here.  If that error exceeds 0.05, AA3's verdict word
  stands without a reading.

  WHAT THIS CANNOT DO.  Eight N over 2.1 decades.  The exponents here
  are this field's, not the published field's, and AA3 is a check
  that the two agree rather than a new measurement of the deficit.
  Nothing here bounds anything, and the required exponent is a
  restatement of the demand, not a route to it.
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
OUT = os.path.join(ROOT, "results", "audit_two_rates.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")
SRCJ = os.path.join(ROOT, "results", "audit_jbarrier.txt")
SRCL = os.path.join(ROOT, "results", "audit_deficit_log.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000]
NGATE = 200_000
RELID = 1e-12
IDCAP = 0.001
GAPCAP = 0.05
NEEDCAP = 0.005
JPRINT = 5e-7


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
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
    rem = np.arange(n + 1, dtype=np.int64)
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


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def split(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    A = 0.0
    B = 0.0
    for d in range(1, (N - 1) // K + 1):
        md = int(mu[d])
        if md == 0 or any(d % q == 0 for q in PN):
            continue
        ms = np.arange(K, (N - 1) // d + 1, dtype=np.int64)
        if ms.size == 0:
            continue
        keep = sqf[ms]
        for q in factor_set(d) | PN:
            keep &= (ms % int(q)) != 0
        if d == 1:
            keep &= lam[ms] == 0.0
        ms = ms[keep]
        if ms.size == 0:
            continue
        c = float((lam[N - d * ms]
                   * np.log(ms.astype(np.float64))).sum())
        if d == 1:
            A = md * c
        else:
            B += md * c
        del ms, keep
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        l2sq += (math.log(k)
                 * float((lam[N - ms * k]
                          * mu[ms].astype(np.float64)).sum())) ** 2
        del ms
    return A, B, math.sqrt(l2sq)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    m = re.search(r"^POINT mainA_%d ([-+]?[\d.eE+-]+)\s*$" % NGATE,
                  io.open(SRCM, encoding="utf-8").read(), re.M)
    j = re.search(r"^POINT jbarrier_%d ([\d.eE+-]+) ([\d.eE+-]+)\s*$"
                  % NGATE, io.open(SRCJ, encoding="utf-8").read(),
                  re.M)
    d = re.search(r"whole-field deficit on the \d+ published: "
                  r"([\d.]+)",
                  io.open(SRCL, encoding="utf-8").read())
    if not m or not j or not d:
        raise SystemExit("a published value is missing")
    return float(m.group(1)), float(j.group(2)), float(d.group(1))


HEAD = [
    "STATISTIC: the exponents of |A|, of the relative residue",
    "           r = |A+B|/|A|, of |S| = |A| r and of l2 on one field,",
    "           the residue exponent the demand requires, and the gap",
    "           between required and achieved.",
    "FIELD: N = %s; the Type II split of rem:bilinear at K ="
    % NS,
    "       floor(N^%.2f). A at N = %d, l2 there, and the published"
    % (THETA, NGATE),
    "       deficit are READ from results/audit_mainterm_removal.txt,",
    "       results/audit_jbarrier.txt and",
    "       results/audit_deficit_log.txt.",
    "DERIVED: |S| = |A| r exactly, so e(|S|) = e(|A|) + e(r) and the",
    "         demand |S| <~ l2 is e(r) <= e(l2) - e(|A|).",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pubA, publ2, pubdef = read_pub()
    say("READ audit_mainterm_removal.txt %d %.17e" % (NGATE, pubA))
    say("READ audit_jbarrier.txt %d %.6f" % (NGATE, publ2))
    say("READ audit_deficit_log.txt deficit %.6f" % pubdef)
    say("PRINTBOUND audit_two_rates %d %.20f" % (17, 5e-18))
    say("  theta %.2f, identity cap %.3f, gap cap %.2f, need cap %.3f"
        % (THETA, IDCAP, GAPCAP, NEEDCAP))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rows = []
    for N in NS:
        A, B, l2 = split(N, lam, mu, sqf)
        S = abs(A + B)
        rows.append((N, A, B, S, l2))
        say("  N = %-9d 1-|B|/|A| %.6f   |S|/l2 %8.4f"
            % (N, abs(A + B) / abs(A), S / l2))
        say("POINT rate_%d %.6f" % (N, abs(A + B) / abs(A)))
        say("POINT sl2_%d %.6f" % (N, S / l2))
    say("SCALES %d" % len(rows))

    # ------------------------------------------------------------- AA1
    say()
    say("AA1  the gate")
    g = [r for r in rows if r[0] == NGATE][0]
    ra = abs(g[1] - pubA) / max(abs(pubA), 1.0)
    rl = abs(g[4] - publ2) / max(abs(publ2), 1.0)
    aa1 = ra <= RELID and rl <= RELID
    bb1 = ra <= RELID and rl <= JPRINT
    say("  A  here %.17e" % g[1])
    say("     its  %.17e   relative %.2e" % (pubA, ra))
    say("  l2 here %.6f against its %.6f   relative %.2e"
        % (g[4], publ2, rl))
    say("  AA1 %s   (cap: %.0e relative on both)"
        % ("hold" if aa1 else "REFUTED", RELID))
    say()
    say("BB1  the gate, registered after AA1 and judging the l2 leg "
        "at its")
    say("     source's own precision")
    say("PRINTBOUND audit_jbarrier 7 %.10f" % JPRINT)
    say("  POINT jbarrier carries seven significant figures, so its "
        "print bound")
    say("  is %.0e relative; the A leg keeps %.0e on a marker widened "
        "to full" % (JPRINT, RELID))
    say("  double precision. TOL BELOW PRINT is what AA1 hit and is "
        "not repeated.")
    say("  A  relative %.2e against %.0e   %s"
        % (ra, RELID, "ok" if ra <= RELID else "MISMATCH"))
    say("  l2 relative %.2e against %.0e   %s"
        % (rl, JPRINT, "ok" if rl <= JPRINT else "MISMATCH"))
    say("  BB1 %s" % ("hold" if bb1 else "REFUTED"))
    if not bb1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    x = np.array([math.log(r[0]) for r in rows])
    eA, seA = fit(x, [math.log(abs(r[1])) for r in rows])
    er, ser = fit(x, [math.log(r[3] / abs(r[1])) for r in rows])
    eS, seS = fit(x, [math.log(r[3]) for r in rows])
    el2, sel2 = fit(x, [math.log(r[4]) for r in rows])
    say()
    say("    quantity          exponent        s.e.")
    for nm, e, s in (("|A|", eA, seA), ("r = |A+B|/|A|", er, ser),
                     ("|S|", eS, seS), ("l2", el2, sel2)):
        say("  %-16s %+.6f     %.6f" % (nm, e, s))
        say("POINT exp_%s %.6f"
            % (nm.split()[0].replace("|", "").replace("=", ""), e))

    # ------------------------------------------------------------- AA2
    say()
    say("AA2  the identity")
    d2 = abs(eS - (eA + er))
    aa2 = d2 <= IDCAP
    say("  e(|S|) %+.6f against e(|A|) + e(r) %+.6f, difference %.2e"
        % (eS, eA + er, d2))
    say("  AA2 %s   (cap: %.3f)"
        % ("hold" if aa2 else "REFUTED", IDCAP))

    # ------------------------------------------------------------- AA3
    say()
    say("AA3  the number")
    need = el2 - eA
    gap = er - need
    egap, segap = fit(x, [math.log(r[3] / r[4]) for r in rows])
    aa3 = gap > 0 and abs(gap - pubdef) <= GAPCAP
    say("  required residue exponent e(l2) - e(|A|) = %+.6f" % need)
    say("  achieved                              e(r) = %+.6f" % er)
    say("  gap %+.6f, and e(|S|) - e(l2) fitted directly is "
        "%+.6f +- %.6f" % (gap, egap, segap))
    say("  against the published deficit %.6f" % pubdef)
    say("TSTAT tworates_gap %.2f" % (egap / segap))
    say("SPREAD tworates_gap %.6f" % segap)
    say("POINT tworates_need %.6f" % need)
    say("POINT tworates_gap %.6f" % gap)
    say("  AA3 %s   (cap: positive and within %.2f of the published)"
        % ("hold" if aa3 else "REFUTED", GAPCAP))
    if segap > GAPCAP:
        say("  UNRESOLVED: the direct fit's error exceeds the cap, so "
            "the verdict")
        say("  stands without a reading, as the rule says")

    # ------------------------------------------------------------- AA4
    say()
    say("AA4  what the demand asks at the top N, in plain terms")
    N, A, B, S, l2 = rows[-1]
    now = abs(A + B) / abs(A)
    needed = now / (S / l2)
    aa4 = needed < NEEDCAP
    say("  at N = %d the cancellation is exact to %.6f" % (N, now))
    say("  to meet the demand there it would have to be exact to "
        "%.6f" % needed)
    say("  a factor of %.4f closer" % (S / l2))
    say("POINT neededcancel %.6f" % needed)
    say("  AA4 %s   (cap: below %.3f)"
        % ("hold" if aa4 else "REFUTED", NEEDCAP))

    say()
    say("=" * 70)
    say("BB1 %s" % ("hold" if bb1 else "REFUTED"))
    say("AA1 %s  AA2 %s  AA3 %s  AA4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (aa1, aa2, aa3, aa4)))
    say()
    if aa3 and aa2:
        say("the two rates are measured and the number between them "
            "is the deficit,")
        say("reached here by a route that never fits |sum a| against "
            "anything: it")
        say("is the residue exponent the demand requires minus the "
            "one the")
        say("cancellation achieves. that agreement is a cross-check "
            "of the")
        say("published number and not a new one.")
    elif not aa3:
        say("the gap on this field is not the published deficit, so "
            "the eight N")
        say("here and the published field are not measuring one "
            "thing and neither")
        say("can be used until that is reconciled.")
    else:
        say("the identity fails, so the three fits are not of the "
            "quantities they")
        say("are labelled with and nothing here is readable.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
