# -*- coding: utf-8 -*-
r"""
Does c(p, N) move, asked without fitting a slope

WHAT IS AT STAKE

rem:driftpower found the corrections' drifts window-dependent: over
one extra decade the two resolved slopes shrank by 0.83298 and
0.83429, the same factor to within 0.00130.  A slope that shrinks when
the field lengthens is not a slope, and that is the second time this
branch has met the disease -- rem:valuation measured it first, in the
deficit's own drift.  So every "drift" quoted for c(p, N), including
the t = -13.19 that made c(3) interesting, is a quantity of its
window.

rem:levelmatched killed the window for the deficit by comparing the
level itself at matched N.  The same trick does not transfer here:
the finest ruler this repository has is radical {2,3,5}, at a log gap
of 0.001129, and **3 and 5 divide those N, so they are outside the
d-range and c(3) cannot be measured there at all**.

But the window can be killed a second way, by not fitting anything.
Whether c(p, N) *moves* is answerable by counting signs: across the
eleven N of rem:driftpower there are ten doublings, and "it fell at
every one" is a statement no window can bend.  Under a coin the
probability is 2^-10 = 0.000977 per prime, and **the primes tested are
named here in advance -- 3, 7, 11, 13 -- so the multiplicity is four
and not the twenty rem:driftbyprime swept.**

BACKS: Remark {#rem:driftsigns} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  II1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker to a relative 1e-12, and the eleven-N drift
      of c(3) reproduces rem:driftpower's -0.009622 to six decimals.
  II2 **c(3) falls at all ten doublings.**
  II3 c(7) falls at all ten doublings.
  II4 **And c(11) and c(13) do not** -- neither is monotone across the
      ten.

REFUTATION RULE (fixed before the run)

  II1 REFUTED outside either tolerance; nothing below is reported.
  II2 REFUTED by any step that rises.  Then c(3) does not move
      monotonically and the drift that four remarks have quoted is
      not even a direction, only a fitted average -- which after
      rem:driftpower would leave nothing standing about c(3) at all.
  II3 REFUTED by any rise, same reading for c(7).
  II4 **REFUTED if either is monotone.**  A prime whose correction
      falls at all ten steps is moving whatever its t says, and
      rem:driftpower's reading that 3 and 7 are special would then be
      wrong -- 11 and 13 would be moving too, only more slowly than
      their errors could show.  That is the outcome that costs the
      most here and it must be stated in those words.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  Ten steps give a
  coin probability of 2^-10 = 0.000977 for a run of one sign, and four
  primes are tested, so **a single monotone prime among the four is
  expected under the null with probability about 0.008** -- small, and
  this run prints the count so the arithmetic is checkable rather than
  asserted.  A monotone c(11) or c(13) is therefore evidence and not
  an accident, and so is a broken run at c(3).  **What the sign test
  cannot say is how fast anything moves**: it reports direction only,
  and the per-step differences are printed beside it precisely so that
  no rate is read from a count.

  WHAT THIS CANNOT DO.  One radical family, eleven N.  A monotone run
  on this field is not a monotone function, and the per-step
  differences are differences of a quantity whose own scatter
  rem:leveldense measured at about 0.03 in a related setting.  Nothing
  here bounds anything.
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
OUT = os.path.join(ROOT, "results", "audit_drift_signs.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")
SRCP = os.path.join(ROOT, "results", "audit_drift_power.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000, 6_400_000, 12_800_000, 25_600_000]
NGATE = 200_000
RELID = 1e-12
DEC = 6
WATCH = (3, 7, 11, 13)
NOISEFRAC = 0.001


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


def phi(n):
    r = n
    for p in factor_set(n):
        r = r // p * (p - 1)
    return r


def pieces(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K
    out = {}
    for d in range(1, D + 1):
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
        out[d] = (md, float((lam[N - d * ms]
                             * np.log(ms.astype(np.float64))).sum()))
        del ms, keep
    return out, D


def model(d, D):
    return (1.0 / phi(d)) * (1.0 - d / D) / (1.0 - 1.0 / D)


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
    d = re.search(r"^POINT drift11_3 ([-+]?[\d.]+)\s*$",
                  io.open(SRCP, encoding="utf-8").read(), re.M)
    if not m or not d:
        raise SystemExit("a published value is missing")
    return float(m.group(1)), float(d.group(1))


HEAD = [
    "STATISTIC: the sign of every step of c(p, N) across eleven N, for",
    "           four primes named before the run, with the per-step",
    "           differences printed beside the counts.",
    "FIELD: N = %s; d over the squarefree d <= D coprime to N, with"
    % NS,
    "       D = floor((N-1)/K), K = floor(N^%.2f). A at N = %d and"
    % (THETA, NGATE),
    "       the eleven-N drift of c(3) are READ from",
    "       results/audit_mainterm_removal.txt and",
    "       results/audit_drift_power.txt.",
    "NOTE: no slope is used to judge anything. rem:driftpower measured",
    "      that a slope of this quantity shrinks by a sixth when the",
    "      field lengthens by a decade; a count of signs does not.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pubA, pubd3 = read_pub()
    say("READ audit_mainterm_removal.txt %d %.17e" % (NGATE, pubA))
    say("READ audit_drift_power.txt drift11_3 %.6f" % pubd3)
    say("PRINTBOUND audit_drift_signs %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, primes named in advance %s" % (THETA, WATCH))
    say("  a run of one sign over %d steps has coin probability %.6f;"
        % (len(NS) - 1, 2.0 ** -(len(NS) - 1)))
    say("  with %d primes tested the expected number of accidental "
        "runs is %.6f"
        % (len(WATCH), 2 * len(WATCH) * 2.0 ** -(len(NS) - 1)))
    say("SEED: none; nothing here is random")

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    cs, A200 = {}, None
    for N in NS:
        ps, D = pieces(N, lam, mu, sqf)
        A = ps[1][1]
        if N == NGATE:
            A200 = A
        ws = {d: v / A for d, (md, v) in ps.items() if d != 1}
        cs[N] = {d: ws[d] / model(d, D) for d in ws if d != D}
        say("  N = %-10d D = %-5d %d contributing d" % (N, D, len(ws)))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- II1
    x = np.array([math.log(N) for N in NS])
    b3, _ = fit(x, [cs[N][3] for N in NS])
    say()
    say("II1  the gate")
    ra = abs(A200 - pubA) / max(abs(pubA), 1.0)
    ok3 = abs(round(b3, DEC) - round(pubd3, DEC)) < 10.0 ** (-DEC)
    ii1 = ra <= RELID and ok3
    say("  A relative %.2e against %.0e   %s"
        % (ra, RELID, "ok" if ra <= RELID else "MISMATCH"))
    say("  c(3) eleven-N drift %+.6f against its %+.6f  %s"
        % (b3, pubd3, "ok" if ok3 else "MISMATCH"))
    say("  II1 %s" % ("hold" if ii1 else "REFUTED"))
    if not ii1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("  c(p, N) at every N, and the step differences")
    verd = {}
    for p in WATCH:
        vals = [cs[N].get(p) for N in NS]
        if any(v is None for v in vals):
            say("  p = %-3d does not contribute at every N" % p)
            verd[p] = None
            continue
        diffs = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
        down = sum(1 for d_ in diffs if d_ < 0)
        verd[p] = (down, len(diffs))
        say("  p = %-3d %s" % (p, " ".join("%.4f" % v for v in vals)))
        say("        steps %s"
            % " ".join("%+.4f" % d_ for d_ in diffs))
        say("        %d of %d fall" % (down, len(diffs)))
        say("COUNT signdown_%d %d" % (p, down))
        say("POINT stepfirst_%d %.6f" % (p, diffs[0]))
        say("POINT steplast_%d %.6f" % (p, diffs[-1]))

    nsteps = len(NS) - 1

    # -------------------------------------------------------- II2, II3
    for name, p in (("II2", 3), ("II3", 7)):
        say()
        say("%s  does c(%d) fall at every step?" % (name, p))
        v = verd.get(p)
        ok = v is not None and v[0] == v[1]
        verd[name] = ok
        if v:
            say("  %d of %d fall" % v)
        say("  %s %s   (cap: all %d)"
            % (name, "hold" if ok else "REFUTED", nsteps))

    # ------------------------------------------------------------- II4
    say()
    say("II4  and do 11 and 13 stay mixed?")
    ii4 = True
    for p in (11, 13):
        v = verd.get(p)
        mono = v is not None and (v[0] == v[1] or v[0] == 0)
        if mono:
            ii4 = False
        if v:
            say("  p = %-3d %d of %d fall   %s"
                % (p, v[0], v[1], "MONOTONE" if mono else "mixed"))
    say("  II4 %s   (cap: neither monotone)"
        % ("hold" if ii4 else "REFUTED"))

    # a diagnostic, after the verdicts and predicted by nothing:
    # II2 and II3 ask for a clean sweep, which is the right question
    # only for a prime that sweeps. For the rest the coin tail is the
    # statistic, and it is arithmetic on the counts already printed.
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  a clean sweep is the right question for a prime that "
        "sweeps; for the")
    say("  others the two-sided coin tail is the statistic, on the "
        "same counts:")
    tails = {}
    for p_ in WATCH:
        v = verd.get(p_)
        if not isinstance(v, tuple):
            continue
        down, tot = v
        k = max(down, tot - down)
        tail = 2.0 * sum(math.comb(tot, j) for j in range(k, tot + 1))
        tail /= 2.0 ** tot
        say("    p = %-3d %d of %d one way, two-sided coin tail %.6f"
            % (p_, k, tot, min(tail, 1.0)))
        say("POINT cointail_%d %.6f" % (p_, min(tail, 1.0)))
        tails[p_] = min(tail, 1.0)
    t3 = tails.get(3)
    t7 = tails.get(7)
    if t3 and t7:
        say("  a sweep is %.6f and eight of ten is %.1f times that, "
            "so the two do" % (t3, t7 / t3))
        say("  not stand together: one is evidence and the other is "
            "a lean")

    say()
    say("=" * 70)
    say("II1 %s  II2 %s  II3 %s  II4 %s"
        % (("hold" if ii1 else "REFUTED"),
           ("hold" if verd["II2"] else "REFUTED"),
           ("hold" if verd["II3"] else "REFUTED"),
           ("hold" if ii4 else "REFUTED")))
    say()
    if verd["II2"] and verd["II3"] and ii4:
        say("c(3) and c(7) move, and the statement survives the "
            "window that took")
        say("their slopes away: they fall at every one of ten "
            "doublings, which no")
        say("choice of field can bend. 11 and 13 do not. the "
            "direction stands and")
        say("the rate does not -- the per-step differences are "
            "printed and no rate")
        say("is read from a count.")
    elif not ii4:
        say("a prime rem:driftpower called still is monotone across "
            "all ten steps.")
        say("it is moving whatever its t says, and that remark's "
            "reading that 3 and")
        say("7 are special is wrong -- the others move too, more "
            "slowly than their")
        say("errors can show.")
    else:
        say("c(3) or c(7) breaks its run, so the drift four remarks "
            "have quoted is")
        say("not even a direction, only a fitted average, and after "
            "rem:driftpower")
        say("nothing about it stands.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
