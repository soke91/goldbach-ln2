# -*- coding: utf-8 -*-
r"""
At what spacing the level's scatter appears

WHAT IS AT STAKE

rem:leveldense established the size: L(N) = log(|sum a| / l2) carries
scatter at fixed radical of about 0.03, seven of twelve close pairs
above a 0.01 cap with the largest at 0.064590 against a smooth
prediction of 0.003428, and that 0.03 is the precision floor for every
radical statement in this branch.  It could not establish what the
scatter is, because 2^a p^b cannot be packed tighter than a log gap of
0.0237 at these N -- the continued fraction of log(5)/log(2) does not
offer a closer approach until 5^28.

Radical {2,3,5} is a far finer ruler.  Enumerating N = 2^a 3^b 5^c
with all three exponents at least one gives, for instance,

    983040 = 2^16 3 5    and    984150 = 2 3^9 5^2

a log gap of 0.001128, one twenty-first of the previous run's finest.
Across a band the gaps run from about 0.001 to 0.05, so **whether
|dL| follows the gap is visible inside one run** rather than needing
two.

    if it follows the gap   L is smooth and the scatter belongs to a
                            scale, which can then be named
    if it does not          the scatter is already there between
                            adjacent N: it is arithmetic noise of L
                            itself, the 0.03 is irreducible, and no
                            refinement of any radical statement in
                            this branch can go below it

BACKS: Remark {#rem:levelfine} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  THE GATE.  N = 1800000 and 2025000, both radical {2,3,5},
      reproduce rem:levelmatched's L of +3.187272 and +3.232112 to
      six decimals.
  U2  **|dL| follows the gap.**  Over the adjacent pairs of the
      enumerated band, the correlation of |dL| with the log gap is
      above 0.5.
  U3  **THE ONE THE FINE RULER IS FOR.**  Pairs with a log gap under
      0.005 have |dL| under 0.01.
  U4  And the level's slope is still the drift: fitted over the band,
      it agrees with {2,3,5}'s published drift to within 0.05.

REFUTATION RULE (fixed before the run)

  U1  REFUTED outside six decimals on either; nothing below is
      reported.
  U2  **REFUTED at 0.5 or below.**  Then |dL| does not follow the
      spacing at all and the scatter is not a scale effect.
  U3  **REFUTED by any fine pair above 0.01.**  Then the scatter is
      present between adjacent N twenty times closer than the
      previous run could reach, which makes it arithmetic noise of L
      rather than structure at a scale.  **That is the outcome that
      settles what this branch can ever claim**: the 0.03 would be
      irreducible, no radical statement here could be sharpened below
      it, and the remark must say so in those words.
  U4  REFUTED outside 0.05.  Then this band's level disagrees with
      the drift measured on doubling families, and the two are not
      measuring one thing.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  U3 depends on how
  many pairs the enumeration puts under 0.005, and that count is not
  known before the run.  **This run prints it.**  With fewer than
  three such pairs U3 tests almost nothing and its verdict word
  stands without a reading -- the remark must then report the count
  and decline the conclusion, exactly as if the cap had not been met.
  The same for U2: the correlation is over however many adjacent
  pairs the band yields, and that count is printed beside it.

  WHAT THIS CANNOT DO.  One radical.  A scatter that is arithmetic
  noise here is not thereby noise at other radicals, and nothing in
  this run explains what the noise is -- only at what spacing it is
  already present.  Nothing here bounds |sum a| or moves item 5.
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
OUT = os.path.join(ROOT, "results", "audit_level_fine.txt")
SRCL = os.path.join(ROOT, "results", "audit_level_matched.txt")
SRCR = os.path.join(ROOT, "results", "audit_radical_law.txt")
SRCD = os.path.join(ROOT, "results", "audit_level_dense.txt")

THETA = 0.56
LO, HI = 700_000, 1_400_000
GATES = [1_800_000, 2_025_000]
DRIFTBASE = 27_000
DEC = 6
CORRCAP = 0.5
FINEGAP = 0.005
FINECAP = 0.01
SLOPECAP = 0.05
MINFINE = 3


def band():
    """every N = 2^a 3^b 5^c in [LO, HI] with a, b, c all >= 1"""
    out = []
    pa = 2
    while pa <= HI:
        pb = pa * 3
        while pb <= HI:
            pc = pb * 5
            while pc <= HI:
                if pc >= LO:
                    out.append(pc)
                pc *= 5
            pb *= 3
        pa *= 2
    return sorted(set(out))


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


def level(N, lam, mu, sqf):
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
    sa = abs(float((lam[N - j] * lk[1:]).sum()))
    del j, lk
    return math.log(sa / math.sqrt(l2sq))


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    return float(b), math.sqrt(float((r ** 2).mean()))


def read_pub():
    src = io.open(SRCL, encoding="utf-8").read()
    lv = {}
    for N in GATES:
        m = re.search(r"^POINT levelmatched_%d ([-+]?[\d.]+)\s*$" % N,
                      src, re.M)
        if not m:
            raise SystemExit("no levelmatched marker for %d" % N)
        lv[N] = float(m.group(1))
    m = re.search(r"^POINT raddrift_%d ([-+]?[\d.]+)\s*$" % DRIFTBASE,
                  io.open(SRCR, encoding="utf-8").read(), re.M)
    if not m:
        raise SystemExit("no raddrift marker for %d" % DRIFTBASE)
    return lv, float(m.group(1))


HEAD = [
    "STATISTIC: L(N) = log(|sum a| / l2) at every N = 2^a 3^b 5^c with",
    "           a, b, c >= 1 in a band, and the difference of",
    "           neighbouring L against the log gap between them.",
    "FIELD: N of radical {2,3,5} in [%d, %d]; k over the squarefree"
    % (LO, HI),
    "       k < N^%.2f coprime to N; j over every index below N."
    % THETA,
    "       Two L are READ from results/audit_level_matched.txt as the",
    "       gate and the radical's drift from",
    "       results/audit_radical_law.txt.",
    "NOTE: the enumeration is exhaustive in the band, so the gaps are",
    "      whatever the arithmetic gives and were not chosen.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    lv, drift = read_pub()
    for N in GATES:
        say("READ audit_level_matched.txt %d %.6f" % (N, lv[N]))
    say("READ audit_radical_law.txt drift %.6f" % drift)
    say("PRINTBOUND audit_level_fine %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, corr cap %.2f, fine gap %.3f, fine cap %.2f,"
        % (THETA, CORRCAP, FINEGAP, FINECAP))
    say("  slope cap %.2f, minimum fine pairs %d" % (SLOPECAP, MINFINE))
    say("RADICALS 1")

    ns = band()
    say("  the band holds %d N of radical {2,3,5}" % len(ns))
    for N in ns:
        if factor_set(N) != {2, 3, 5}:
            raise SystemExit("%d is not of radical {2,3,5}" % N)
    say("  every one checked to have that radical exactly")
    ex1, ex2 = 983_040, 984_150
    say("  the finest pair the enumeration reaches near a million is")
    say("  %d = 2^16*3*5 and %d = 2*3^9*5^2, log gap %.6f"
        % (ex1, ex2, math.log(ex2) - math.log(ex1)))
    say("POINT examplegap %.6f" % (math.log(ex2) - math.log(ex1)))
    say("  this band spans %.4f in log N; the two sets of "
        "rem:leveldense span"
        % (math.log(HI) - math.log(LO)))
    dn = [int(m) for m in re.findall(r"^POINT dense_(\d+) ",
                                     io.open(SRCD,
                                             encoding="utf-8").read(),
                                     re.M)]
    for pr_ in (5, 7):
        g = [n for n in dn if n % pr_ == 0
             and factor_set(n) == {2, pr_}]
        say("    {2,%d}  %.4f over %d N"
            % (pr_, math.log(max(g)) - math.log(min(g)), len(g)))
        say("READ audit_level_dense.txt span_%d %.4f"
            % (pr_, math.log(max(g)) - math.log(min(g))))
    say("COPRIME %d"
        % len(set(tuple(sorted(factor_set(n))) for n in ns)))
    say("  one coprimality class: a, b, c are all at least one, so no")
    say("  N here drops a prime and the k-range is the same for all "
        "of them")

    NMAX = max(ns + GATES)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    # -------------------------------------------------------------- U1
    say()
    say("U1  the gate")
    u1 = True
    for N in GATES:
        L = level(N, lam, mu, sqf)
        g = abs(L - lv[N]) < 10.0 ** (-DEC)
        u1 &= g
        say("  N = %-9d here %+.6f against its %+.6f  %s"
            % (N, L, lv[N], "ok" if g else "MISMATCH"))
    say("  U1 %s   (cap: %d decimals)"
        % ("hold" if u1 else "REFUTED", DEC))
    if not u1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rows = []
    for N in ns:
        L = level(N, lam, mu, sqf)
        rows.append((N, math.log(N), L))
        say("  N = %-9d log N %.6f   L %+.6f" % (N, math.log(N), L))
        say("POINT fine_%d %.6f" % (N, L))
    say("SCALES %d" % len(rows))

    gaps, dls = [], []
    say()
    say("    N1        N2         log gap    |dL|      smooth")
    for i in range(len(rows) - 1):
        n1, x1, l1 = rows[i]
        n2, x2, l2v = rows[i + 1]
        g = x2 - x1
        d = abs(l2v - l1)
        gaps.append(g)
        dls.append(d)
        say("  %-9d %-9d  %.6f  %.6f  %.6f  %s"
            % (n1, n2, g, d, g * drift,
               "fine" if g < FINEGAP else ""))
        say("POINT finepair_%d %.6f" % (n1, d))
    gaps = np.array(gaps)
    dls = np.array(dls)

    # -------------------------------------------------------------- U2
    say()
    say("U2  does |dL| follow the gap?")
    c = float(np.corrcoef(gaps, dls)[0, 1])
    u2 = c > CORRCAP
    say("  correlation %+.6f over %d adjacent pairs" % (c, len(gaps)))
    say("CORR levelfine_regressors %.6f" % abs(c))
    say("POINT finecorr %.6f" % c)
    say("  U2 %s   (cap: above %.2f)"
        % ("hold" if u2 else "REFUTED", CORRCAP))

    # -------------------------------------------------------------- U3
    say()
    say("U3  are the finest pairs quiet?")
    fine = [(g, d) for g, d in zip(gaps, dls) if g < FINEGAP]
    u3 = all(d <= FINECAP for _, d in fine)
    say("  %d pairs under a log gap of %.3f" % (len(fine), FINEGAP))
    for g, d in fine:
        say("    gap %.6f  |dL| %.6f  smooth %.6f  %s"
            % (g, d, g * drift, "ok" if d <= FINECAP else "ABOVE"))
    say("COUNT finepairs %d" % len(fine))
    if fine:
        say("POINT finemax %.6f" % max(d for _, d in fine))
    say("  U3 %s   (cap: %.2f on every fine pair)"
        % ("hold" if u3 else "REFUTED", FINECAP))
    if len(fine) < MINFINE:
        say("  UNDERPOWERED: fewer than %d fine pairs, so U3's verdict"
            % MINFINE)
        say("  stands without a reading, as the rule says")

    # -------------------------------------------------------------- U4
    say()
    say("U4  is the band's slope the drift?")
    sl, rms = fit([r[1] for r in rows], [r[2] for r in rows])
    u4 = abs(sl - drift) <= SLOPECAP
    say("  slope %+.6f against drift %+.6f, r.m.s. residual %.6f"
        % (sl, drift, rms))
    say("POINT fineslope %.6f" % sl)
    say("POINT finerms %.6f" % rms)
    say("  U4 %s   (cap: %.2f)"
        % ("hold" if u4 else "REFUTED", SLOPECAP))

    say()
    say("=" * 70)
    say("U1 %s  U2 %s  U3 %s  U4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (u1, u2, u3, u4)))
    say()
    if not u3 and len(fine) >= MINFINE:
        say("the scatter is already there between adjacent N twenty "
            "times closer")
        say("than the last run could reach. it is arithmetic noise of "
            "L itself and")
        say("not structure at a scale, the 0.03 is irreducible, and "
            "no radical")
        say("statement in this branch can be sharpened below it.")
    elif u3 and u2 and len(fine) >= MINFINE:
        say("|dL| follows the spacing and the finest pairs are quiet, "
            "so L is")
        say("smooth and the scatter belongs to a scale this band can "
            "name.")
    elif len(fine) < MINFINE:
        say("the band did not supply enough fine pairs to decide. "
            "what the")
        say("enumeration gives is what it gives, and this question "
            "needs a")
        say("radical whose multiples pack closer than {2,3,5} does "
            "here.")
    else:
        say("the finest pairs are quiet but |dL| does not follow the "
            "gap, so the")
        say("scatter is neither a clean scale effect nor present at "
            "the finest")
        say("spacing, and what it is remains unmeasured.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
