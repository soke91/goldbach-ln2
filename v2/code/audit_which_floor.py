# -*- coding: utf-8 -*-
r"""
Which floor is item 5's demand actually measured against

WHAT IS AT STAKE

rem:jbarrierreach closed a candidate answer and left the question:
item 5 asks |sum a| <~ l2, but l2 is a norm over the dilations k while
sum a is a sum over the indices j, and the j-side barrier D is a
different quantity.  On family A the two barriers sit on top of each
other -- D/l2 flat near a half, drift +0.005736 at t = +1.73 -- so
that field cannot tell which of them |sum a| is related to.  **Family
B breaks the degeneracy**: there D/l2 falls, drift -0.068713 at
t = -7.90, so the two floors separate by a resolved amount and the
sum must track at most one of them.

That is the whole question, and it is one run away.  If |sum a|/l2
drifts at the same rate in both families, l2 is the floor the deficit
is measured against and item 5's +0.134019 is a radical-independent
statement about a k-side norm.  If |sum a|/D is the one that agrees
across families, then the published deficit is measured against the
wrong object and its number carries the radical.

rem:jbarrier computed |sum a| on family A and rem:jbarrierreach
dropped it to compute D and the shared diagonal; |sum a| on family B
has never been computed.  This run is blind in the half that decides.

BACKS: Remark {#rem:whichfloor} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  THE GATE.  On family A at N = 200000 this reproduces the
      published |sum a| -- the POINT marker of
      results/audit_deficit_direct.txt -- and the published D/l2 --
      the POINT marker of results/audit_jbarrier.txt -- to four
      decimals.

      *Disclosed.*  M1 as first written read |sum a| from the *table*
      of results/audit_jbarrier.txt, which prints it with two
      decimals, and compared it at four.  It came out 87895.3236
      against 87895.3200 and M1 was refuted -- by the print bound and
      not by a disagreement.  That is TOL BELOW PRINT in this
      repository's own vocabulary, and G75 exists for it.  The script
      exits at M1, so no verdict on M2 to M4 existed when this was
      found.  The gate now reads markers that carry the digits it
      asks for; the tolerance is unchanged.

      Worth recording beside it: this run's care went into naming M2's
      unresolved case before the fact, which is where the last three
      failures had been, and the defect appeared in the gate instead.
  M2  **l2 is the floor.**  The drift of log(|sum a|/l2) against
      log N agrees between the two families: the two slopes differ by
      less than 0.02.
  M3  **And D is not.**  The drift of log(|sum a|/D) differs between
      the families by more than 0.02.  This is not independent of M2
      -- the two differ by the D/l2 drift rem:jbarrierreach measured
      at -0.068713 -- and it is registered so that the pair is read
      together and not one without the other.
  M4  |sum a| stands above both barriers at every N of both families.

REFUTATION RULE (fixed before the run)

  M1  REFUTED outside four decimals on either; nothing below is
      reported.  The substitution of markers for the table row is
      disclosed above and changes which digits are read, not what is
      required.
  M2  **REFUTED outside 0.02, and then the deficit carries the
      radical.**  Item 5's +0.134019 would be a number of the
      primorial-free family it was measured on and not a property of
      the demand, and every remark quoting it would inherit that.
      That is the outcome that costs the most and it must be stated
      plainly if it comes.
  M3  REFUTED at or below 0.02.  Then both ratios agree across
      families, which cannot happen unless the D/l2 drift is not what
      rem:jbarrierreach measured, and the disagreement would be with
      that run rather than with either floor.
  M4  REFUTED by any N where either inequality fails.

  **THE UNRESOLVED CASE, NAMED.**  Eleven N and nine N give drifts
  with standard errors this run prints.  A difference of slopes
  smaller than the error of that difference is not a difference, and
  M2 holding inside such an error is not evidence that the drifts
  agree -- it is evidence that this field cannot tell.  The error of
  the difference is printed beside the caps, and if the two slopes
  differ by less than it **M2's verdict word stands but its reading
  is barred**, exactly as rem:jbarrier's J2 was barred.  This is the
  fourth time a cap in this branch has been written on a quantity
  whose error was not the one that matters, so it is named before the
  run and not after.

  WHAT THIS CANNOT DO.  Two radicals.  Nothing here supplies a bound
  on |sum a|, and D remains a heuristic floor rather than a theorem.
  No forecast, no closure N; rem:shapepower and rem:deficitlog stand.
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
OUT = os.path.join(ROOT, "results", "audit_which_floor.txt")
SRC = os.path.join(ROOT, "results", "audit_jbarrier.txt")
SRCD = os.path.join(ROOT, "results",
                    "audit_deficit_direct.txt")
SRCL = os.path.join(ROOT, "results",
                    "audit_deficit_log.txt")

THETA = 0.56
NSA = [25_000 * (1 << j) for j in range(11)]
NSB = [30_030 * (1 << j) for j in range(9)]
NGATE = 200_000
DEC = 4
SLOPECAP = 0.02


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


def triple(N, lam, mu, sqf):
    """|sum a|, l2 over k, D over j -- one pass over the k-range"""
    PN = factor_set(N)
    K = int(N ** THETA)
    lk = np.zeros(N, dtype=np.float64)
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        mm = mu[ms].astype(np.float64)
        l2sq += (lg * float((lam[N - ms * k] * mm).sum())) ** 2
        lk[ms * k] += lg * mm
        del ms, mm
    j = np.arange(1, N, dtype=np.int64)
    w = lam[N - j] * lk[1:]
    sa = abs(float(w.sum()))
    dsq = float((w * w).sum())
    del j, w, lk
    return sa, math.sqrt(l2sq), math.sqrt(dsq)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    m = re.search(r"^POINT jbarrier_%d ([\d.eE+-]+) ([\d.eE+-]+)\s*$"
                  % NGATE, io.open(SRC, encoding="utf-8").read(), re.M)
    q = re.search(r"^POINT deficitdirect_%d ([\d.eE+-]+) " % NGATE,
                  io.open(SRCD, encoding="utf-8").read(), re.M)
    if not m or not q:
        raise SystemExit("missing published values for N = %d" % NGATE)
    w = re.search(r"whole-field deficit on the \d+ published: "
                  r"([\d.]+)",
                  io.open(SRCL, encoding="utf-8").read())
    if not w:
        raise SystemExit("no whole-field deficit in audit_deficit_log")
    d, l2 = float(m.group(1)), float(m.group(2))
    return float(q.group(1)), d / l2, float(w.group(1))


HEAD = [
    "STATISTIC: |sum a| against each of the two square-root barriers,",
    "           l2 over the dilations k and D over the indices j, on",
    "           two radical families; the drift of each ratio within",
    "           each family and whether the drifts agree across them.",
    "FIELD: family A, N = 25000*2^j for j < %d, radical {2,5};"
    % len(NSA),
    "       family B, N = 30030*2^j for j < %d, radical" % len(NSB),
    "       {2,3,5,7,11,13}. k over the squarefree k < N^%.2f coprime"
    % THETA,
    "       to N; j over every index below N. |sum a| and D/l2 at",
    "       N = %d are READ from results/audit_jbarrier.txt as the"
    % NGATE,
    "       gate.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    psa, pr, pdef = read_pub()
    say("READ audit_deficit_direct.txt %d %.4f" % (NGATE, psa))
    say("READ audit_jbarrier.txt ratio %.4f" % pr)
    say("READ audit_deficit_log.txt deficit %.6f" % pdef)
    say("  |sum a| and D/l2 at the gate N, both from POINT markers")
    say("  NOTE, disclosed: M1 first read |sum a| from the table of")
    say("  audit_jbarrier.txt, which prints two decimals, and judged "
        "it at four.")
    say("  It read the two-decimal value there and M1 was refuted "
        "by the print")
    say("  bound, not by a disagreement -- TOL BELOW PRINT. The "
        "script exits at")
    say("  M1, so no")
    say("  verdict on M2 to M4 existed when this was found. The "
        "tolerance is")
    say("  unchanged; only the source of the digits is.")
    say("PRINTBOUND audit_which_floor %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, slope cap %.2f" % (THETA, SLOPECAP))
    say("TOL NOT FROM PRINT audit_which_floor")
    say("RADICALS 2")

    NMAX = max(NSA + NSB)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    fams = {}
    for tag, ns in (("A", NSA), ("B", NSB)):
        say()
        say("family %s, radical %s" % (tag, sorted(factor_set(ns[0]))))
        rows = []
        for N in ns:
            sa, l2, d = triple(N, lam, mu, sqf)
            rows.append((N, sa, l2, d))
            say("  N = %-10d |sum a| %13.2f  /l2 %8.4f  /D %8.4f"
                % (N, sa, sa / l2, sa / d))
            say("POINT floor%s_%d %.6f %.6f" % (tag, N, sa / l2, sa / d))
        fams[tag] = rows
    say("SCALES %d" % (len(NSA) + len(NSB)))

    # -------------------------------------------------------------- M1
    say()
    say("M1  the gate: does family A reproduce the published pair?")
    g = [r for r in fams["A"] if r[0] == NGATE][0]
    a = abs(round(g[1], DEC) - round(psa, DEC)) < 10.0 ** (-DEC)
    b = abs(round(g[3] / g[2], DEC) - round(pr, DEC)) < 10.0 ** (-DEC)
    m1 = a and b
    say("  |sum a| here %.4f against its %.4f  %s"
        % (g[1], psa, "ok" if a else "MISMATCH"))
    say("  D/l2    here %.4f against its %.4f  %s"
        % (g[3] / g[2], pr, "ok" if b else "MISMATCH"))
    say("  M1 %s   (cap: %d decimals on both)"
        % ("hold" if m1 else "REFUTED", DEC))
    if not m1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    sl = {}
    say()
    say("    family   ratio        drift          s.e.        t")
    for tag in ("A", "B"):
        rows = fams[tag]
        x = np.array([math.log(r[0]) for r in rows])
        for nm, idx in (("|sum a|/l2", 2), ("|sum a|/D", 3)):
            y = np.log(np.array([r[1] / r[idx] for r in rows]))
            bb, se = fit(x, y)
            sl[(tag, nm)] = (bb, se)
            say("       %s   %-11s %+.6f     %.6f   %+6.2f"
                % (tag, nm, bb, se, bb / se))
            say("TSTAT floor%s%s %.2f"
                % (tag, "L2" if idx == 2 else "D", bb / se))
            say("SPREAD floor%s%s %.6f"
                % (tag, "L2" if idx == 2 else "D", se))

    # -------------------------------------------------------------- M2
    say()
    say("M2  do the |sum a|/l2 drifts agree across the families?")
    ba, sa_ = sl[("A", "|sum a|/l2")]
    bb, sb = sl[("B", "|sum a|/l2")]
    dl = ba - bb
    sdl = math.sqrt(sa_ ** 2 + sb ** 2)
    m2 = abs(dl) < SLOPECAP
    say("  A %+.6f, B %+.6f, difference %+.6f" % (ba, bb, dl))
    say("  error of the difference %.6f" % sdl)
    say("POINT l2driftgap %.6f" % dl)
    say("  M2 %s   (cap: %.2f)"
        % ("hold" if m2 else "REFUTED", SLOPECAP))
    if abs(dl) < sdl:
        say("  UNRESOLVED: the gap is inside the error of the "
            "difference, so the")
        say("  verdict word stands and the reading is barred, as the "
            "rule says")

    # -------------------------------------------------------------- M3
    say()
    say("M3  and do the |sum a|/D drifts disagree?")
    ca, sca = sl[("A", "|sum a|/D")]
    cb, scb = sl[("B", "|sum a|/D")]
    dd = ca - cb
    sdd = math.sqrt(sca ** 2 + scb ** 2)
    m3 = abs(dd) > SLOPECAP
    say("  A %+.6f, B %+.6f, difference %+.6f" % (ca, cb, dd))
    say("  error of the difference %.6f" % sdd)
    say("POINT ddriftgap %.6f" % dd)
    say("  M3 %s   (cap: above %.2f)"
        % ("hold" if m3 else "REFUTED", SLOPECAP))

    # -------------------------------------------------------------- M4
    say()
    say("M4  does the sum stand above both barriers everywhere?")
    bad = [(t, r[0]) for t in ("A", "B") for r in fams[t]
           if not (r[1] > r[2] and r[1] > r[3])]
    m4 = not bad
    say("  failures: %s" % (bad if bad else "none"))
    say("  M4 %s   (cap: every N of both families)"
        % ("hold" if m4 else "REFUTED"))

    # a diagnostic, after the verdicts and predicted by nothing:
    # the two families span different x, and the deficit's own drift
    # is not constant, so the ranges could be the whole difference
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  the families span different log N and rem:deficitdirect "
        "measured the")
    say("  deficit's drift as itself drifting, so the ranges could "
        "be doing this.")
    say("  refitting both on the log N they share:")
    lo = max(math.log(NSA[0]), math.log(NSB[0]))
    hi = min(math.log(NSA[-1]), math.log(NSB[-1]))
    say("  shared window %.4f to %.4f" % (lo, hi))
    mt = {}
    for tag in ("A", "B"):
        rr = [r for r in fams[tag]
              if lo - 1e-9 <= math.log(r[0]) <= hi + 1e-9]
        xx = np.array([math.log(r[0]) for r in rr])
        yy = np.log(np.array([r[1] / r[2] for r in rr]))
        bb, se = fit(xx, yy)
        mt[tag] = (bb, se)
        say("       %s   %d N   |sum a|/l2 drift %+.6f +- %.6f"
            % (tag, len(rr), bb, se))
    gap2 = mt["A"][0] - mt["B"][0]
    se2 = math.sqrt(mt["A"][1] ** 2 + mt["B"][1] ** 2)
    say("  difference %+.6f, error of the difference %.6f, t %+.2f"
        % (gap2, se2, gap2 / se2))
    say("POINT matchedgap %.6f" % gap2)
    say("TSTAT matchedgap %.2f" % (gap2 / se2))
    say("SPREAD matchedgap %.6f" % se2)

    say()
    say("=" * 70)
    say("M1 %s  M2 %s  M3 %s  M4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (m1, m2, m3, m4)))
    say()
    if m2 and abs(dl) >= sdl:
        say("l2 is the floor. the deficit drifts at the same rate on "
            "two radicals")
        say("whose j-side barriers differ by a resolved amount, so "
            "item 5's number")
        say("is a statement about the k-side norm and does not carry "
            "the radical.")
    elif not m2:
        say("the deficit carries the radical. item 5's %.6f is a "
            "number of the" % pdef)
        say("family it was measured on and not a property of the "
            "demand, and every")
        say("remark quoting it inherits that. this is the outcome "
            "M2's rule named as")
        say("costing the most.")
    else:
        say("the gap is inside the error of the difference. this "
            "field does not tell")
        say("which floor the deficit is measured against, and no "
            "reading is drawn --")
        say("the degeneracy family B was meant to break is not "
            "broken by these N.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
