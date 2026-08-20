# -*- coding: utf-8 -*-
r"""
Is the residue an elementary Mobius sum with no primes in it

WHAT IS AT STAKE

rem:tworates measured what item 5 asks: the cancellation between
A = I(1) and B = sum_{d>=2} mu(d) I(d) is exact to 0.023066 at
N = 3200000 and would have to be exact to 0.002425 -- a factor of
9.5118 -- and the residue exponent must improve from -0.342848 to
-0.478577.  It ended by asking where that factor lives among the d.

Writing this run produced a candidate answer that removes the primes
entirely.  I(d) = sum_{K <= m < N/d} mu^2(m) Lambda(N-dm) log m counts
prime powers N - dm, which lie in the residue class N mod d, so its
size is governed by 1/phi(d); and m runs over [K, N/d), a range that
shrinks to nothing as d approaches D = N/K.  Writing w(d) = I(d)/A,

    A + B = A * sum_d mu(d) w(d),
    model:  w(d) ~ (1/phi(d)) * (1 - d/D) / (1 - 1/D).

**If that model carries the residue, then the whole of item 5's
remaining difficulty is the rate at which the elementary sum
sum_{d <= D} mu(d) (1 - d/D) / phi(d) approaches zero** -- a Fejer-
weighted Mobius sum with no primes and no Lambda in it, and a
classical question rather than this repository's.

This is a model and not an identity.  The run tests whether it
reproduces the residue, and reports the failure if it does not.

BACKS: Remark {#rem:residuemodel} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  CC1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      widened POINT mainA marker to a relative 1e-12.
  CC2 **The shape.**  For the d carrying the largest ten weights,
      the measured w(d) agrees with the model to within 20 per cent.
  CC3 **The model carries the residue.**  |sum_d mu(d) w_model(d)|
      is within a factor of 2 of the measured residue
      |A + B| / |A| at every N.
  CC4 And it moves at the same rate: fitting both against log N, the
      two exponents differ by less than 0.05.

REFUTATION RULE (fixed before the run)

  CC1 REFUTED outside 1e-12; nothing below is reported.  The marker
      was widened to full double precision in the previous tick after
      three ticks of TOL BELOW PRINT, so a failure here is a real
      disagreement.
  CC2 REFUTED outside 20 per cent on any of the ten.  Then 1/phi(d)
      with a linear taper is the wrong shape for w and the reduction
      below does not follow, whatever CC3 says -- a sum can come out
      right for the wrong reasons.
  CC3 **REFUTED outside a factor of 2, and that is the one that
      decides.**  Then the residue is not the elementary sum, the
      primes do not drop out, and item 5's difficulty stays where it
      was: in a correlation of Lambda against a sieve weight.
  CC4 REFUTED outside 0.05.  A model that lands near the value but
      moves at a different rate is a coincidence at these N, and the
      remark must say so rather than claim the reduction.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  The model's own
  sum can be near the measured residue for the trivial reason that
  both are small; the comparison is a ratio and the ratio is printed
  at every N, not just fitted.  **If the model's sum and the measured
  residue differ in sign at any N, CC3's factor test is meaningless
  there** -- a magnitude within a factor of two of a quantity of the
  opposite sign is not agreement -- and this run prints both signs so
  that cannot be hidden.  The exponents in CC4 are fitted on eight
  points and their errors are printed; an error above 0.05 leaves
  CC4's verdict without a reading.

  WHAT THIS CANNOT DO.  Eight N over 2.1 decades.  A model that
  reproduces a residue on this field is not thereby the asymptotics,
  and the classical rate for the Fejer-weighted Mobius sum is not
  computed here -- naming the question is not answering it.  Nothing
  in this run bounds anything.
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
OUT = os.path.join(ROOT, "results", "audit_residue_model.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000]
NGATE = 200_000
RELID = 1e-12
SHAPEPC = 20.0
TOPW = 10
FACTOR = 2.0
EXPCAP = 0.05


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
    """(d, mu(d), I(d)) for every contributing d, and D"""
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K
    out = []
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
        out.append((d, md,
                    float((lam[N - d * ms]
                           * np.log(ms.astype(np.float64))).sum())))
        del ms, keep
    return out, D


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
    if not m:
        raise SystemExit("no mainA marker for N = %d" % NGATE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the weights w(d) = I(d)/I(1) of the Type II sum",
    "           against the model (1/phi(d))(1-d/D)/(1-1/D), and the",
    "           residue sum_d mu(d) w(d) against the model's own",
    "           elementary Mobius sum.",
    "FIELD: N = %s; d over the squarefree d <= D = floor((N-1)/K)"
    % NS,
    "       coprime to N and m over the squarefree m in [K, N/d)",
    "       coprime to dN, with K = floor(N^%.2f). A at N = %d is"
    % (THETA, NGATE),
    "       READ from results/audit_mainterm_removal.txt at full",
    "       double precision.",
    "MODEL: N - dm lies in the class N mod d so I(d) is governed by",
    "       1/phi(d), and m runs over a range that empties as d ->",
    "       D, which the linear taper stands for. This is a model,",
    "       not an identity.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pubA = read_pub()
    say("READ audit_mainterm_removal.txt %d %.17e" % (NGATE, pubA))
    say("PRINTBOUND audit_residue_model %d %.20f" % (17, 5e-18))
    say("  theta %.2f, shape %.0f per cent on the top %d, factor %.1f,"
        % (THETA, SHAPEPC, TOPW, FACTOR))
    say("  exponent cap %.2f" % EXPCAP)

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rows = []
    for N in NS:
        ps, D = pieces(N, lam, mu, sqf)
        A = [v for d, md, v in ps if d == 1][0]
        meas = sum(md * v for _, md, v in ps) / A
        mod = sum(md * (1.0 / phi(d)) * (1.0 - d / D)
                  for d, md, _ in ps) / (1.0 - 1.0 / D)
        rows.append((N, D, A, meas, mod, ps))
        say("  N = %-9d D = %-5d residue %+.6f   model %+.6f   "
            "ratio %8.4f"
            % (N, D, meas, mod, mod / meas if meas else float("nan")))
        say("POINT resmeas_%d %.6e" % (N, meas))
        say("POINT resmodel_%d %.6e" % (N, mod))
    say("SCALES %d" % len(rows))

    # ------------------------------------------------------------- CC1
    say()
    say("CC1  the gate")
    g = [r for r in rows if r[0] == NGATE][0]
    rel = abs(g[2] - pubA) / max(abs(pubA), 1.0)
    cc1 = rel <= RELID
    say("  A here %.17e" % g[2])
    say("    its %.17e   relative %.2e" % (pubA, rel))
    say("  CC1 %s   (cap: %.0e relative)"
        % ("hold" if cc1 else "REFUTED", RELID))
    if not cc1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------------- CC2
    say()
    say("CC2  the shape, on the top %d weights at N = %d"
        % (TOPW, NS[-1]))
    N, D, A, meas, mod, ps = rows[-1]
    ws = sorted(((abs(v / A), d, md, v) for d, md, v in ps),
                reverse=True)[:TOPW]
    cc2 = True
    say("      d     w measured    w model      per cent")
    for _, d, md, v in ws:
        wm = v / A
        wp = (1.0 / phi(d)) * (1.0 - d / D) / (1.0 - 1.0 / D)
        pc = 100.0 * (wm - wp) / wp
        cc2 &= abs(pc) <= SHAPEPC
        say("  %7d  %+.6f   %+.6f   %+8.2f" % (d, wm, wp, pc))
        say("POINT shape_%d %.4f" % (d, pc))
    say("  CC2 %s   (cap: %.0f per cent)"
        % ("hold" if cc2 else "REFUTED", SHAPEPC))

    # ------------------------------------------------------------- CC3
    say()
    say("CC3  does the model carry the residue?")
    cc3 = True
    signs = []
    for N, D, A, meas, mod, _ in rows:
        same = (meas > 0) == (mod > 0)
        signs.append(same)
        f = abs(mod / meas) if meas else float("inf")
        ok = same and (1.0 / FACTOR) <= f <= FACTOR
        cc3 &= ok
        say("  N = %-9d measured %+.6f  model %+.6f  signs %s  "
            "factor %.4f  %s"
            % (N, meas, mod, "same" if same else "OPPOSITE", f,
               "ok" if ok else "no"))
    say("  CC3 %s   (cap: same sign and within a factor %.1f)"
        % ("hold" if cc3 else "REFUTED", FACTOR))
    if not all(signs):
        say("  NOTE: a magnitude within a factor of two of a quantity "
            "of the")
        say("  opposite sign is not agreement, as the rule says")

    # ------------------------------------------------------------- CC4
    say()
    say("CC4  do they move at the same rate?")
    x = np.array([math.log(r[0]) for r in rows])
    e1, s1 = fit(x, [math.log(abs(r[3])) for r in rows])
    e2, s2 = fit(x, [math.log(abs(r[4])) for r in rows])
    cc4 = abs(e1 - e2) <= EXPCAP
    say("  measured %+.6f +- %.6f, model %+.6f +- %.6f, gap %+.6f"
        % (e1, s1, e2, s2, e1 - e2))
    say("TSTAT residuemodel_gap %.2f"
        % ((e1 - e2) / math.sqrt(s1 ** 2 + s2 ** 2)))
    say("SPREAD residuemodel_gap %.6f"
        % math.sqrt(s1 ** 2 + s2 ** 2))
    say("POINT resexp_meas %.6f" % e1)
    say("POINT resexp_model %.6f" % e2)
    say("  CC4 %s   (cap: %.2f)"
        % ("hold" if cc4 else "REFUTED", EXPCAP))
    if max(s1, s2) > EXPCAP:
        say("  UNRESOLVED: an exponent's own error exceeds the cap, "
            "so CC4's")
        say("  verdict stands without a reading, as the rule says")

    say()
    say("=" * 70)
    say("CC1 %s  CC2 %s  CC3 %s  CC4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (cc1, cc2, cc3, cc4)))
    say()
    if cc2 and cc3 and cc4:
        say("the residue is an elementary Mobius sum. the primes drop "
            "out: what")
        say("item 5 needs is the rate at which sum_{d<=D} mu(d)(1 - "
            "d/D)/phi(d)")
        say("approaches zero, a Fejer-weighted Mobius sum with no "
            "Lambda in it.")
        say("that is a classical question and not this repository's, "
            "which is the")
        say("whole of what this run buys -- naming it is not "
            "answering it.")
    elif not cc3:
        say("the model does not carry the residue. the primes do not "
            "drop out and")
        say("item 5's difficulty stays where it was, in a correlation "
            "of Lambda")
        say("against a sieve weight.")
    elif not cc2:
        say("the model's sum lands near the residue while its shape "
            "is wrong, so")
        say("it comes out right for the wrong reasons and the "
            "reduction does not")
        say("follow.")
    else:
        say("the model lands near the residue but moves at a "
            "different rate, so")
        say("the agreement is a coincidence at these N and no "
            "reduction is claimed.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
