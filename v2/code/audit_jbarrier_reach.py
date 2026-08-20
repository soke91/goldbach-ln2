# -*- coding: utf-8 -*-
r"""
The half, on a second radical and two octaves further, and where it comes from

WHAT IS AT STAKE

rem:jbarrier measured the two square-root barriers of one sum -- l2
over the dilations k and D over the indices j -- and found D/l2 flat
at 0.498823 across 2.41 decades.  It then said plainly what that
measurement is worth: K2 and K3 were registered *after* the first run
printed those ratios, so they confirm a reading of data already seen,
and no mechanism was offered for the constant.  The remark named what
would establish it: **a second radical family and a longer field,
measured on the ratio rather than on the exponents.**  This run is
that test, and it is blind -- family B has never been computed.

And the mechanism is available by algebra rather than by fitting.
Expand each barrier and take its diagonal.  For D, the square of
LK(j) = sum_{k|j} mu(j/k) log k has diagonal sum_{k|j} mu^2(j/k)
(log k)^2, so

    D^2 diagonal = sum_k (log k)^2 sum_m mu^2(m) Lambda(N-mk)^2 .

For l2, the square of H(N;k) = sum_m mu(m) Lambda(N-mk) has diagonal
sum_m mu^2(m) Lambda(N-mk)^2, so

    l2^2 diagonal = sum_k (log k)^2 sum_m mu^2(m) Lambda(N-mk)^2 .

**They are the same quantity.**  Call it DIAG.  So a ratio of one half
means D^2 = DIAG and l2^2 = 4 DIAG, or some other split of the two
off-diagonals that lands there; measuring DIAG against both says which,
and a mechanism for the half is exactly the statement that l2 stands
at twice the common diagonal while D stands at it.

BACKS: Remark {#rem:jbarrierreach} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  L1  THE GATE.  On the nine N of rem:jbarrier this reproduces its
      published ratios -- 0.5090 first and 0.4998 last -- to four
      decimals.
  L2  **The extension keeps it.**  Family A carried two octaves
      further, eleven N to 2.56e7, has mean D/l2 within 0.02 of 0.5.
  L3  **The second radical gives the same constant.**  Family B,
      N = 30030*2^j with radical {2,3,5,7,11,13} against family A's
      {2,5}, has mean D/l2 within 0.02 of 0.5.
  L4  **And the constant does not depend on the radical**: the two
      families' means differ by less than 0.02.
  L5  D is the common diagonal: D^2/DIAG is within 10 per cent of 1.
  L6  And l2 stands at twice it: l2^2/DIAG is within 10 per cent of 4.

REFUTATION RULE (fixed before the run)

  L1  REFUTED outside four decimals; nothing below is reported.
  L2  REFUTED outside 0.02.  Then the flatness rem:jbarrier measured
      does not survive its own field being lengthened, and the
      constant was a feature of 2.41 decades.
  L3  REFUTED outside 0.02.
  L4  **REFUTED outside 0.02, and this is the one that would matter
      most.**  A radical-dependent constant would say the relation
      between the two barriers carries the same arithmetic factor
      prod_{p|N} p/(p+1) that rem:targetderived measured in #k -- not
      that the half is wrong, but that it is not a universal number
      and rem:jbarrier's reading of a single constant is.  **That
      case is named here so it cannot afterwards be presented as
      agreement**, and if it fires the two family means are to be
      reported beside their arithmetic factors.
  L5  REFUTED outside 10 per cent.  Then D is not the diagonal, its
      off-diagonal is not negligible, and the half has no mechanism
      from this route.
  L6  REFUTED outside 10 per cent.  L5 and L6 are independent: either
      can fail alone, and if L5 holds while L6 fails the half is not
      explained and the failure locates the unexplained part in l2.

  WHAT THIS CANNOT DO.  DIAG is an identity's diagonal, not a bound.
  Two radicals is two, not many, and a constant that survives both is
  not thereby universal.  Nothing here measures |sum a| against
  anything or moves item 5's demand; rem:shapepower and
  rem:deficitlog stand.
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
OUT = os.path.join(ROOT, "results", "audit_jbarrier_reach.txt")
SRC = os.path.join(ROOT, "results", "audit_jbarrier.txt")

THETA = 0.56
NSA = [25_000 * (1 << j) for j in range(11)]
NSB = [30_030 * (1 << j) for j in range(9)]
DEC = 4
HALF = 0.5
HALFCAP = 0.02
FAMCAP = 0.02
PCCAP = 10.0
FOUR = 4.0


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


def barriers(N, lam, mu, sqf):
    """l2 over k, D over j, and the diagonal both of them share"""
    PN = factor_set(N)
    K = int(N ** THETA)
    lk = np.zeros(N, dtype=np.float64)
    l2sq = diag = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        mm = mu[ms].astype(np.float64)
        lv = lam[N - ms * k]
        h = float((lv * mm).sum())
        l2sq += (lg * h) ** 2
        diag += lg * lg * float(((mm * lv) ** 2).sum())
        lk[ms * k] += lg * mm
        del ms, mm, lv
    j = np.arange(1, N, dtype=np.int64)
    w = lam[N - j] * lk[1:]
    dsq = float((w * w).sum())
    del j, w, lk
    return math.sqrt(l2sq), math.sqrt(dsq), diag


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^SERIES jbarrier_ratio (.+)$", src, re.M)
    if not m:
        raise SystemExit("no SERIES marker in audit_jbarrier.txt")
    return [float(v) for v in m.group(1).split()]


HEAD = [
    "STATISTIC: the ratio D/l2 of the two square-root barriers of one",
    "           sum, on two radical families and two octaves beyond",
    "           rem:jbarrier's field, and each barrier against the",
    "           diagonal DIAG that the two share by algebra.",
    "FIELD: family A, N = 25000*2^j for j < %d, radical {2,5};"
    % len(NSA),
    "       family B, N = 30030*2^j for j < %d, radical" % len(NSB),
    "       {2,3,5,7,11,13}. k over the squarefree k < N^%.2f coprime"
    % THETA,
    "       to N; j over every index below N. The nine ratios of",
    "       rem:jbarrier are READ from results/audit_jbarrier.txt as",
    "       the gate.",
    "DERIVED: the diagonals of D^2 and of l2^2 are the same quantity,",
    "         sum_k (log k)^2 sum_m mu^2(m) Lambda(N-mk)^2.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_jbarrier.txt first %.4f" % pub[0])
    say("READ audit_jbarrier.txt last %.4f" % pub[-1])
    say("  the published ratios, %d of them, this run's gate" % len(pub))
    say("PRINTBOUND audit_jbarrier_reach %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  half %.1f +- %.2f, family cap %.2f, per cent cap %.1f, "
        "four %.1f" % (HALF, HALFCAP, FAMCAP, PCCAP, FOUR))
    say("RADICALS 2")

    NMAX = max(NSA + NSB)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    fams = {}
    for tag, ns in (("A", NSA), ("B", NSB)):
        rad = sorted(factor_set(ns[0]))
        arith = 1.0
        for p in rad:
            arith *= p / (p + 1.0)
        say()
        say("family %s, radical %s, arithmetic factor %.6f"
            % (tag, rad, arith))
        rows = []
        for N in ns:
            l2, d, dg = barriers(N, lam, mu, sqf)
            rows.append((N, l2, d, dg))
            say("  N = %-10d l2 %12.2f  D %12.2f  D/l2 %.6f"
                % (N, l2, d, d / l2))
            say("POINT reachratio%s_%d %.6f" % (tag, N, d / l2))
        fams[tag] = (rows, arith)
    say("SCALES %d" % (len(NSA) + len(NSB)))

    # -------------------------------------------------------------- L1
    say()
    say("L1  the gate: do the shared N reproduce the published ratios?")
    shared = {r[0]: r for r in fams["A"][0]}
    l1 = True
    for pv, N in zip(pub, [25_000 * (1 << j) for j in range(len(pub))]):
        r = shared.get(N)
        if r is None:
            l1 = False
            say("  N = %-10d missing" % N)
            continue
        here = r[2] / r[1]
        g = abs(round(here, DEC) - round(pv, DEC)) < 10.0 ** (-DEC)
        l1 &= g
        say("  N = %-10d here %.4f against its %.4f  %s"
            % (N, here, pv, "ok" if g else "MISMATCH"))
    say("  L1 %s   (cap: %d decimals)"
        % ("hold" if l1 else "REFUTED", DEC))
    if not l1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    ma = float(np.mean([r[2] / r[1] for r in fams["A"][0]]))
    mb = float(np.mean([r[2] / r[1] for r in fams["B"][0]]))

    # -------------------------------------------------------------- L2
    say()
    say("L2  does the extension keep the half?")
    l2h = abs(ma - HALF) <= HALFCAP
    say("  family A over %d N to %d: mean %.6f, off %+.6f"
        % (len(NSA), NSA[-1], ma, ma - HALF))
    say("POINT meanA %.6f" % ma)
    say("  L2 %s   (cap: %.2f)"
        % ("hold" if l2h else "REFUTED", HALFCAP))

    # -------------------------------------------------------------- L3
    say()
    say("L3  does the second radical give the same constant?")
    l3 = abs(mb - HALF) <= HALFCAP
    say("  family B over %d N to %d: mean %.6f, off %+.6f"
        % (len(NSB), NSB[-1], mb, mb - HALF))
    say("POINT meanB %.6f" % mb)
    say("  L3 %s   (cap: %.2f)"
        % ("hold" if l3 else "REFUTED", HALFCAP))

    # -------------------------------------------------------------- L4
    say()
    say("L4  does the constant depend on the radical?")
    l4 = abs(ma - mb) <= FAMCAP
    say("  A %.6f (factor %.6f) against B %.6f (factor %.6f)"
        % (ma, fams["A"][1], mb, fams["B"][1]))
    say("  difference %+.6f" % (ma - mb))
    say("POINT famgap %.6f" % (ma - mb))
    say("ACROSS jbarrier_reach %.6f" % abs(ma - mb))
    say("  L4 %s   (cap: %.2f)"
        % ("hold" if l4 else "REFUTED", FAMCAP))

    # ---------------------------------------------------------- L5, L6
    say()
    say("L5, L6  where does the half come from?")
    say("    family  N            D^2/DIAG    l2^2/DIAG")
    rd, rl = [], []
    for tag in ("A", "B"):
        for N, l2, d, dg in fams[tag][0]:
            a = d * d / dg
            b = l2 * l2 / dg
            rd.append(a)
            rl.append(b)
            say("       %s    %-11d %9.5f  %11.5f" % (tag, N, a, b))
    md = float(np.mean(rd))
    ml = float(np.mean(rl))
    pd = 100.0 * abs(md - 1.0) / 1.0
    pl = 100.0 * abs(ml - FOUR) / FOUR
    l5 = pd <= PCCAP
    l6 = pl <= PCCAP
    say("  mean D^2/DIAG  %.5f, off %.2f per cent from 1" % (md, pd))
    say("  mean l2^2/DIAG %.5f, off %.2f per cent from %.1f"
        % (ml, pl, FOUR))
    say("POINT dsqdiag %.6f" % md)
    say("POINT l2sqdiag %.6f" % ml)
    say("  L5 %s   L6 %s   (cap: %.1f per cent each)"
        % ("hold" if l5 else "REFUTED",
           "hold" if l6 else "REFUTED", PCCAP))

    # a diagnostic, run after the verdicts above and predicted by
    # nothing: the statistic the caps should have been written on
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  L5's cap was on the mean of a series; the series is "
        "printed above and")
    say("  runs far wider than the cap. the drift is what decides "
        "whether the")
    say("  two barriers share a scale, and no prediction was written "
        "on it.")
    say("    family   ratio drift        s.e.        t     range")
    for tag in ("A", "B"):
        rr = fams[tag][0]
        xx = np.array([math.log(r[0]) for r in rr])
        yy = np.log(np.array([r[2] / r[1] for r in rr]))
        b, a0 = np.polyfit(xx, yy, 1)
        res = yy - (b * xx + a0)
        se = math.sqrt(float((res ** 2).sum() / (len(xx) - 2))
                       / float(((xx - xx.mean()) ** 2).sum()))
        rt = [r[2] / r[1] for r in rr]
        say("       %s     %+.6f     %.6f   %+6.2f   %.4f-%.4f"
            % (tag, b, se, b / se, min(rt), max(rt)))
        say("TSTAT reachdrift%s %.2f" % (tag, b / se))
        say("SPREAD reachdrift%s %.6f" % (tag, se))
        say("POINT reachdrift%s %.6f" % (tag, b))
    say("  a drift resolved in either family is the two barriers not "
        "sharing a")
    say("  scale, whatever the mean of their ratio is")

    say()
    say("=" * 70)
    say("L1 %s  L2 %s  L3 %s  L4 %s  L5 %s  L6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (l1, l2h, l3, l4, l5, l6)))
    say()
    if l2h and l3 and l4:
        say("the half survives a second radical and two more octaves, "
            "blind. it is")
        say("a constant of this construction and not of the field "
            "rem:jbarrier had.")
    elif not l4:
        say("the constant depends on the radical. the relation "
            "between the two")
        say("barriers carries arithmetic, as #k does, and "
            "rem:jbarrier's single")
        say("number is the wrong shape for it -- the two family means "
            "and their")
        say("arithmetic factors are printed above and that is the "
            "statement.")
    else:
        say("the half does not survive its own test. what "
            "rem:jbarrier measured")
        say("was a feature of its field and is withdrawn as a "
            "constant.")
    if l5 and l6:
        say("and it has a mechanism: D stands at the diagonal the two "
            "barriers")
        say("share and l2 stands at twice it, so the ratio is one "
            "half by algebra")
        say("and not by coincidence.")
    elif l5:
        say("D is the shared diagonal but l2 is not twice it, so the "
            "half is not")
        say("explained and the unexplained part sits in l2's "
            "off-diagonal.")
    else:
        say("D is not the shared diagonal, so this route offers no "
            "mechanism for")
        say("the ratio and none is claimed.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
