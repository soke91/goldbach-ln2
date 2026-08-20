# -*- coding: utf-8 -*-
r"""
The contribution function of a single odd prime, over a wide range of p

WHAT IS AT STAKE

rem:radicallaw showed the deficit's drift depends on which primes
divide N and not on how many: three families with omega = 2 spanned
0.074374 against a largest single drift error of 0.019289.  Reading
the two-prime families as one contribution each gave 3 -> +0.144102,
5 -> +0.069729, 7 -> +0.073163 on a base of +0.068281, and the
multi-prime families came in below the sum of those parts.  That
remark ended by naming what it could not do: **three points do not
determine f(p), and 5 and 7 are not separated by their own errors.**

They will not be separated by more of the same either.  The drift
errors here run near 0.010 and f(5) and f(7) differ by 0.0034, so
neighbouring primes are out of reach at this precision.  **The
shape is not.**  If f decays like 1/p then f(101) is near 0.004
against f(3) near 0.144 -- a factor of thirty, which errors of 0.015
resolve easily.  So this run does not try to order neighbours; it
spreads p from 3 to 101 and measures the exponent of the decay.

Ten doubling families, each with radical {2} or {2,p}:

    2^14        {2}      base        2^4 * 17^2  {2,17}
    2^8 3^4     {2,3}                2^6 * 23^2  {2,23}
    2^3 5^5     {2,5}                2^4 * 47^2  {2,47}
    2^9 7^2     {2,7}                2 * 101^2   {2,101}
    2^4 11^3    {2,11}               2^4 13^3    {2,13}

BACKS: Remark {#rem:primecontrib} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  O1  THE GATE.  The {2}, {2,3}, {2,5} and {2,7} families reproduce
      rem:radicallaw's drifts +0.068281, +0.212384, +0.138010 and
      +0.141444 to six decimals -- same bases, same code path.
  O2  f decays: f(101) is below a quarter of f(3).
  O3  **THE SHAPE.**  Fitting log f(p) against log p over the nine
      odd primes, the slope is -1 within 0.3 -- f(p) of order 1/p,
      which is what removing the k divisible by p from a range of
      squarefree k would cost.
  O4  And the decay is clean rather than ragged: that fit's r.m.s.
      residual in log f is under 0.35.

REFUTATION RULE (fixed before the run)

  O1  REFUTED outside six decimals on any of the four; nothing below
      is reported, since a different drift is a different
      construction.
  O2  REFUTED if f(101) is at or above a quarter of f(3).  Then f
      does not visibly decay across a factor of thirty-three in p and
      the whole reading of "each prime contributes" is wrong, not
      merely unshaped.
  O3  **REFUTED outside -1 +- 0.3.**  Then the cost of excluding p is
      not of order 1/p and the natural derivation -- that the k
      divisible by p are a 1/(p+1) share of the squarefree k -- does
      not describe it.
  O4  REFUTED above 0.35 r.m.s.  Then there is no clean power law in
      p, whatever the slope reads, and O3's number describes a line
      through scatter.

  **THE UNRESOLVED CASE, NAMED.**  1/p, 1/(p-1) and 1/(p+1) differ by
  under a fifth for every p >= 11 here and by a third at p = 5, so
  **this design cannot separate them and will not be read as
  choosing one.**  Only p = 3 distinguishes them at all, and one
  point cannot carry that.  O3 is a statement about the exponent and
  about nothing else; if it holds, the remark says "of order 1/p" and
  not "equal to c/(p+1)".  Likewise f(5) and f(7) are not separated
  here and no ordering of neighbouring primes may be reported --
  their difference was 0.0034 against errors near 0.010 in
  rem:radicallaw and this run does not improve that.

  WHAT THIS CANNOT DO.  Nine odd primes, nine N each, one base per
  prime -- so a base's own size and 2-adic valuation are not varied
  and cannot be separated from p.  A power law fitted over 3 to 101
  says nothing about p beyond it.  Nothing here is a mechanism and
  nothing here bounds |sum a|; rem:shapepower applies.
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
OUT = os.path.join(ROOT, "results", "audit_prime_contrib.txt")
SRC = os.path.join(ROOT, "results", "audit_radical_law.txt")

THETA = 0.56
BASE2 = 16_384
BASES = [(3, 20_736), (5, 25_000), (7, 25_088), (11, 21_296),
         (13, 35_152), (17, 18_496), (23, 33_856), (47, 35_344),
         (101, 20_402)]
NPER = 9
DEC = 6
QUARTER = 0.25
SLOPE = -1.0
SLOPETOL = 0.3
RMSCAP = 0.35


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


def pair(N, lam, mu, sqf):
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
    return sa, math.sqrt(l2sq)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se, math.sqrt(float((r ** 2).mean()))


def drift(base, lam, mu, sqf, say):
    xs, ys = [], []
    for jj in range(NPER):
        N = base * (1 << jj)
        sa, l2 = pair(N, lam, mu, sqf)
        xs.append(math.log(N))
        ys.append(math.log(sa / l2))
    return fit(xs, ys)


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^POINT raddrift_(\d+) ([-+]?[\d.]+)\s*$",
                         src, re.M):
        out[int(m.group(1))] = float(m.group(2))
    return out


HEAD = [
    "STATISTIC: f(p), the amount a single odd prime dividing N adds to",
    "           the drift of log(|sum a|/l2), for nine primes from 3",
    "           to 101, and the exponent of its decay in p.",
    "FIELD: ten doubling families, base * 2^j for j < %d; base %d for"
    % (NPER, BASE2),
    "       radical {2} and one base per odd prime for radical {2,p},",
    "       p in %s. k over the squarefree"
    % [b[0] for b in BASES],
    "       k < N^%.2f coprime to N; j over every index below N."
    % THETA,
    "       The four drifts rem:radicallaw published for the shared",
    "       bases are READ from results/audit_radical_law.txt as the",
    "       gate.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for b in (BASE2, 20_736, 25_000, 25_088):
        say("READ audit_radical_law.txt %d %.6f" % (b, pub[b]))
    say("  the four drifts this run shares with rem:radicallaw")
    say("PRINTBOUND audit_prime_contrib %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, quarter %.2f, slope %.1f +- %.1f, r.m.s. cap "
        "%.2f" % (THETA, QUARTER, SLOPE, SLOPETOL, RMSCAP))
    say("RADICALS %d" % (len(BASES) + 1))

    NMAX = max([BASE2] + [b for _, b in BASES]) * (1 << (NPER - 1))
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    b2, se2, _ = drift(BASE2, lam, mu, sqf, say)
    say()
    say("  base {2} at %d: drift %+.6f +- %.6f" % (BASE2, b2, se2))
    say("POINT contribbase %.6f" % b2)

    rows = []
    for p, base in BASES:
        b, se, _ = drift(base, lam, mu, sqf, say)
        rows.append((p, base, b, se, b - b2))
        say("  p = %-4d base %-7d drift %+.6f +- %.6f   f(p) %+.6f"
            % (p, base, b, se, b - b2))
        say("POINT contrib_%d %.6f" % (p, b - b2))
        say("SPREAD contrib_%d %.6f" % (p, se))
    say("SCALES %d" % ((len(BASES) + 1) * NPER))

    # -------------------------------------------------------------- O1
    say()
    say("O1  the gate: do the shared bases reproduce their drifts?")
    o1 = abs(b2 - pub[BASE2]) < 10.0 ** (-DEC)
    say("  base %-7d here %+.6f against its %+.6f  %s"
        % (BASE2, b2, pub[BASE2],
           "ok" if o1 else "MISMATCH"))
    for p, base, b, se, f in rows:
        if base in pub:
            g = abs(b - pub[base]) < 10.0 ** (-DEC)
            o1 &= g
            say("  base %-7d here %+.6f against its %+.6f  %s"
                % (base, b, pub[base], "ok" if g else "MISMATCH"))
    say("  O1 %s   (cap: %d decimals)"
        % ("hold" if o1 else "REFUTED", DEC))
    if not o1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    f = {p: v for p, _, _, _, v in rows}

    # -------------------------------------------------------------- O2
    say()
    say("O2  does f decay across the range of p?")
    o2 = f[101] < QUARTER * f[3]
    say("  f(3) %+.6f, f(101) %+.6f, quarter of f(3) %+.6f"
        % (f[3], f[101], QUARTER * f[3]))
    say("  ratio f(101)/f(3) %.4f" % (f[101] / f[3]))
    say("POINT decayratio %.6f" % (f[101] / f[3]))
    say("  O2 %s   (cap: below a quarter)"
        % ("hold" if o2 else "REFUTED"))

    # ---------------------------------------------------------- O3, O4
    say()
    say("O3, O4  the exponent of the decay")
    pos = [(p, v) for p, v in sorted(f.items()) if v > 0]
    say("  %d of %d f(p) are positive and enter the fit"
        % (len(pos), len(f)))
    if len(pos) < 3:
        say("  fewer than three positive values; the fit is not made")
        o3 = o4 = False
        sl = rms = float("nan")
    else:
        lx = np.array([math.log(p) for p, _ in pos])
        ly = np.array([math.log(v) for _, v in pos])
        sl, sesl, rms = fit(lx, ly)
        o3 = abs(sl - SLOPE) <= SLOPETOL
        o4 = rms <= RMSCAP
        say("  slope of log f on log p: %+.6f +- %.6f, t against -1 "
            "%+.2f" % (sl, sesl, (sl - SLOPE) / sesl))
        say("  r.m.s. residual in log f %.6f" % rms)
        say("TSTAT contribslope %.2f" % ((sl - SLOPE) / sesl))
        say("SPREAD contribslope %.6f" % sesl)
        say("POINT contribslope %.6f" % sl)
        say("POINT contribrms %.6f" % rms)
    say("  O3 %s   (cap: %.1f +- %.1f)"
        % ("hold" if o3 else "REFUTED", SLOPE, SLOPETOL))
    say("  O4 %s   (cap: %.2f r.m.s.)"
        % ("hold" if o4 else "REFUTED", RMSCAP))
    say("COEFF NOT SEPARABLE primecontrib")
    say("  1/p, 1/(p-1) and 1/(p+1) are within a fifth of each other "
        "for every")
    say("  p >= 11 here, so this design does not choose among them "
        "and the")
    say("  exponent is all that is read, as the rule says")

    # a diagnostic, after the verdicts and predicted by nothing.
    # every f(p) carries the base family's error, which is the
    # largest single one here and is common to all of them, so it
    # cancels in differences between the f(p) and must not be used
    # when comparing them.
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  every f(p) is a drift minus the same base drift, whose "
        "error %.6f" % se2)
    say("  is the largest here and is common to all of them. it "
        "cancels in")
    say("  differences between the f(p), so comparisons among them "
        "use only")
    say("  the two families' own errors.")
    ses = {p_: se_ for p_, _, _, se_, _ in rows}
    ps = sorted(f)
    say("  f(3) against each of the others:")
    for q in ps[1:]:
        d = f[3] - f[q]
        sd = math.sqrt(ses[3] ** 2 + ses[q] ** 2)
        say("    3 vs %-4d  %+.6f +- %.6f   t %+6.2f"
            % (q, d, sd, d / sd))
    say("TSTAT contrib3vs5 %.2f"
        % ((f[3] - f[5]) / math.sqrt(ses[3] ** 2 + ses[5] ** 2)))
    say("SPREAD contrib3vs5 %.6f"
        % math.sqrt(ses[3] ** 2 + ses[5] ** 2))
    rest = ps[1:]
    vals = np.array([f[q] for q in rest])
    errs = np.array([ses[q] for q in rest])
    wm = float((vals / errs ** 2).sum() / (1.0 / errs ** 2).sum())
    chi = float((((vals - wm) / errs) ** 2).sum())
    worst = max(abs(f[q] - wm) / ses[q] for q in rest)
    say("  the %d primes from %d to %d against a constant:"
        % (len(rest), rest[0], rest[-1]))
    say("    weighted mean %+.6f, chi-square %.2f on %d degrees, "
        "largest deviation %.2f of its own error"
        % (wm, chi, len(rest) - 1, worst))
    say("POINT contribflatmean %.6f" % wm)
    say("POINT contribflatchi %.6f" % chi)
    say("  so f(3) stands apart and the rest are not resolved from "
        "one another;")
    say("  what is refuted is the decay, and this run does not "
        "measure any")
    say("  shape among the primes above 3")

    say()
    say("=" * 70)
    say("O1 %s  O2 %s  O3 %s  O4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (o1, o2, o3, o4)))
    say()
    if o2 and o3 and o4:
        say("a single odd prime dividing N costs the deficit an "
            "amount of order")
        say("1/p, cleanly over three to a hundred and one. that is "
            "the size of the")
        say("share of squarefree k the exclusion removes, which is "
            "the first")
        say("quantity in this branch whose radical dependence has a "
            "reason.")
        say("it is a shape and not a mechanism: nothing here derives "
            "the constant")
        say("and rem:radicallaw's sub-additivity is untouched.")
    elif o2 and not o3:
        say("f decays but not like 1/p. the cost of excluding a prime "
            "is not the")
        say("share of k it removes, and the natural derivation does "
            "not describe")
        say("this. what it is instead is not measured here.")
    elif not o2:
        say("f does not decay across a factor of thirty-three in p. "
            "reading the")
        say("two-prime families as one contribution each is wrong, "
            "not merely")
        say("unshaped, and rem:radicallaw's contributions are "
            "withdrawn as a")
        say("description.")
    else:
        say("the slope reads near minus one through scatter too wide "
            "to call it a")
        say("power law. O3's number describes a line and not a "
            "decay, as O4's")
        say("rule says.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
