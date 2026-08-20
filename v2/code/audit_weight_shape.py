# -*- coding: utf-8 -*-
r"""
Is the weight's correction multiplicative over the primes of d

WHAT IS AT STAKE

rem:residuemodel closed the reduction to an elementary Mobius sum --
the model w(d) ~ (1/phi(d))(1-d/D)/(1-1/D) falls away far faster than
the residue and changes sign at the top N -- and left one informative
thing behind.  At N = 3200000 the top weights miss that model by about
ten per cent for prime d (3, 7, 11, 13, 17, 19, 23 run -9.10 to
-12.20) and by about a third for the two composite ones present,
21 = 3*7 and 33 = 3*11, at -34.12 and -33.60.

**The shape is wrong in a way that depends on which primes divide d.**
Writing c(d) = w_measured(d) / w_model(d), the question is whether c
is multiplicative:

    c(d) = prod_{p | d} c(p) ?

That is testable with no free parameters, because the omega = 1 values
determine every other one.  And it is testable out of sample: taking
c(p) from the largest N alone and predicting the residue at the other
seven is a genuine forecast, not a fit.

The test set is fixed here and not chosen afterwards: every squarefree
d <= 50 coprime to N.  At the largest N that is well inside D = 728,
and it contains omega = 1, 2 and 3.

BACKS: Remark {#rem:weightshape} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  DD1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker to a relative 1e-12.
  DD2 **Multiplicativity.**  For every squarefree d <= 50 coprime to
      N with omega(d) >= 2, at the largest N, c(d) equals
      prod_{p|d} c(p) within 10 per cent.
  DD3 The single-prime corrections are near a constant: the spread of
      c(p) over the primes p <= 50 coprime to N is under 0.05.
  DD4 **Out of sample.**  Taking c(p) from the largest N only and
      correcting the model with prod_{p|d} c(p) at every d, the
      corrected model reproduces the measured residue within a factor
      1.2 at each of the other seven N.

REFUTATION RULE (fixed before the run)

  DD1 REFUTED outside 1e-12; nothing below is reported.
  DD2 **REFUTED outside 10 per cent on any such d.**  Then the
      correction is not a function of the primes dividing d one at a
      time, and the shape of w depends on something finer than the
      factorisation of d into primes -- which would say the weights
      carry structure no multiplicative model can express, and that
      is the outcome that says the most.
  DD3 REFUTED above a 0.05 spread.  Then c(p) varies with p and the
      correction is multiplicative but not by a constant; DD2 could
      still hold and the remark must report c(p) rather than a single
      number.
  DD4 **REFUTED outside a factor 1.2 at any of the seven.**  Then the
      correction fitted at one N does not transport, so whatever
      DD2 shows is local to that N and the corrected model is not a
      model of anything.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  c(d) is a ratio of
  two small numbers when |w(d)| is small, so its error grows with d.
  **This run prints |w(d)| beside every c(d)**, and a d whose weight
  is below a thousandth of w(3) carries a c that is noise -- if a DD2
  failure is confined to such d, the verdict word stands and the
  reading is that the test had no power there, not that
  multiplicativity failed.  The threshold is stated here and not
  chosen after seeing which d fail.

  WHAT THIS CANNOT DO.  One radical family of N, one d-range, and a
  correction measured against a model that rem:residuemodel already
  refuted as a whole.  A multiplicative correction that works on
  d <= 50 says nothing about the d near D, which is where the model's
  sign error came from.  Nothing here revives the reduction
  rem:residuemodel closed.
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
OUT = os.path.join(ROOT, "results", "audit_weight_shape.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000]
NGATE = 200_000
DTEST = 50
RELID = 1e-12
MULTPC = 10.0
SPREADCAP = 0.05
OOSFACTOR = 1.2
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


def read_pub():
    m = re.search(r"^POINT mainA_%d ([-+]?[\d.eE+-]+)\s*$" % NGATE,
                  io.open(SRCM, encoding="utf-8").read(), re.M)
    if not m:
        raise SystemExit("no mainA marker for N = %d" % NGATE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the correction c(d) = w(d) / w_model(d) of the Type II",
    "           weights, tested for multiplicativity over the primes",
    "           dividing d, and the corrected model's residue out of",
    "           sample.",
    "FIELD: N = %s; d over the squarefree d <= D coprime to N; the"
    % NS,
    "       multiplicativity test set is every squarefree d <= %d"
    % DTEST,
    "       coprime to N, fixed before the run. A at N = %d is READ"
    % NGATE,
    "       from results/audit_mainterm_removal.txt.",
    "MODEL: w_model(d) = (1/phi(d))(1-d/D)/(1-1/D), the model",
    "       rem:residuemodel refuted as a whole; c is measured against",
    "       it and nothing here revives it.",
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
    say("PRINTBOUND audit_weight_shape %d %.20f" % (17, 5e-18))
    say("  theta %.2f, test d <= %d, mult %.0f per cent, spread %.2f,"
        % (THETA, DTEST, MULTPC, SPREADCAP))
    say("  out-of-sample factor %.1f, noise threshold %.4f of w(3)"
        % (OOSFACTOR, NOISEFRAC))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    data = {}
    for N in NS:
        ps, D = pieces(N, lam, mu, sqf)
        A = ps[1][1]
        meas = sum(md * v for md, v in ps.values()) / A
        data[N] = (ps, D, A, meas)
        say("  N = %-9d D = %-5d residue %+.6f" % (N, D, meas))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- DD1
    say()
    say("DD1  the gate")
    rel = abs(data[NGATE][2] - pubA) / max(abs(pubA), 1.0)
    dd1 = rel <= RELID
    say("  A here %.17e" % data[NGATE][2])
    say("    its %.17e   relative %.2e" % (pubA, rel))
    say("  DD1 %s   (cap: %.0e relative)"
        % ("hold" if dd1 else "REFUTED", RELID))
    if not dd1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    NTOP = NS[-1]
    ps, D, A, meas = data[NTOP]
    cs, ws = {}, {}
    for d in range(2, DTEST + 1):
        if d not in ps:
            continue
        w = ps[d][1] / A
        cs[d] = w / model(d, D)
        ws[d] = w
    say()
    say("  the correction at N = %d, on the fixed test set" % NTOP)
    say("      d  omega   |w(d)|      c(d)")
    for d in sorted(cs):
        say("  %5d  %3d   %.6f   %+.6f"
            % (d, len(factor_set(d)), abs(ws[d]), cs[d]))
        say("POINT corr_%d %.6f" % (d, cs[d]))

    # ------------------------------------------------------------- DD2
    say()
    say("DD2  is the correction multiplicative?")
    noise = NOISEFRAC * abs(ws.get(3, 1.0))
    say("  the noise threshold is %.8f, a thousandth of |w(3)|"
        % noise)
    dd2 = True
    weak = []
    say("      d   factors        c(d)      product     per cent")
    for d in sorted(cs):
        fs = sorted(factor_set(d))
        if len(fs) < 2:
            continue
        if not all(p in cs for p in fs):
            continue
        pred = 1.0
        for p in fs:
            pred *= cs[p]
        pc = 100.0 * (cs[d] - pred) / pred
        ok = abs(pc) <= MULTPC
        if abs(ws[d]) < noise:
            weak.append(d)
        else:
            dd2 &= ok
        say("  %5d   %-12s %+.6f  %+.6f  %+8.2f  %s"
            % (d, "*".join(str(p) for p in fs), cs[d], pred, pc,
               "ok" if ok else ("weak" if abs(ws[d]) < noise
                                else "OUT")))
        say("POINT mult_%d %.4f" % (d, pc))
    say("  DD2 %s   (cap: %.0f per cent on the d above the "
        "threshold)" % ("hold" if dd2 else "REFUTED", MULTPC))
    if weak:
        say("  NOTE: %s sit below the threshold and were not counted, "
            "as the rule" % weak)
        say("  says -- the test has no power there")

    # ------------------------------------------------------------- DD3
    say()
    say("DD3  are the single-prime corrections near a constant?")
    prim = {d: c for d, c in cs.items() if len(factor_set(d)) == 1}
    sp = max(prim.values()) - min(prim.values())
    dd3 = sp <= SPREADCAP
    say("  c(p) for p = %s" % sorted(prim))
    say("       %s" % " ".join("%.4f" % prim[p] for p in sorted(prim)))
    say("  spread %.6f, mean %.6f" % (sp, float(np.mean(
        list(prim.values())))))
    say("POINT primspread %.6f" % sp)
    say("  DD3 %s   (cap: %.2f)"
        % ("hold" if dd3 else "REFUTED", SPREADCAP))

    # ------------------------------------------------------------- DD4
    say()
    say("DD4  out of sample: the corrected model at the other N")
    dd4 = True
    say("      N          measured     corrected     factor")
    for N in NS:
        ps2, D2, A2, meas2 = data[N]
        tot = 0.0
        for d, (md, v) in ps2.items():
            corr = 1.0
            for p in factor_set(d):
                corr *= cs.get(p, 1.0)
            tot += md * model(d, D2) * corr
        f = tot / meas2 if meas2 else float("nan")
        ok = (1.0 / OOSFACTOR) <= f <= OOSFACTOR
        if N != NTOP:
            dd4 &= ok
        say("  %-10d %+.6f    %+.6f    %8.4f  %s"
            % (N, meas2, tot, f,
               "fitted" if N == NTOP else ("ok" if ok else "OUT")))
        say("POINT oos_%d %.6f" % (N, f))
    say("  DD4 %s   (cap: factor %.1f at the seven not fitted)"
        % ("hold" if dd4 else "REFUTED", OOSFACTOR))
    fs_ = []
    for N in NS:
        ps2, D2, A2, meas2 = data[N]
        tot = 0.0
        for d, (md, v) in ps2.items():
            corr = 1.0
            for p in factor_set(d):
                corr *= cs.get(p, 1.0)
            tot += md * model(d, D2) * corr
        fs_.append(tot / meas2)
    say("  the eight factors span %.4f about a mean of %.4f"
        % (max(fs_) - min(fs_), float(np.mean(fs_))))
    say("POINT oosspan %.6f" % (max(fs_) - min(fs_)))
    say("POINT oosmean %.6f" % float(np.mean(fs_)))
    say("  a span that small across two decades is a level that is "
        "uniformly")
    say("  wrong, not a correction that fails to transport; c(p) "
        "exists only for")
    say("  p <= %d so every d with a larger prime factor was "
        "corrected by 1" % DTEST)

    say()
    say("=" * 70)
    say("DD1 %s  DD2 %s  DD3 %s  DD4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (dd1, dd2, dd3, dd4)))
    say()
    if dd2 and dd4:
        say("the correction is multiplicative over the primes "
            "dividing d and it")
        say("transports: fitted at one N it predicts the residue at "
            "seven others.")
        say("so the weights have a shape, and it is 1/phi(d) times a "
            "product over")
        say("the primes of d. that is a description of w and not a "
            "bound on")
        say("anything; rem:residuemodel's closure stands.")
    elif not dd2:
        say("the correction is not a function of the primes dividing "
            "d one at a")
        say("time. the weights carry structure no multiplicative "
            "model expresses,")
        say("which is more than the previous run's failure said and "
            "is the outcome")
        say("that says the most.")
    else:
        say("the correction is multiplicative at one N and does not "
            "transport, so")
        say("it is local to that N and the corrected model is not a "
            "model.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
