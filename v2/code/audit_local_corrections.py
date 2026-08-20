# -*- coding: utf-8 -*-
r"""
Are the corrections themselves functions of N

WHAT IS AT STAKE

rem:omegafactor found the weight model right in shape and wrong in
transport: with c(p) and kappa fitted at N = 3200000 the corrected
residue came in at 1.0448 there -- inside its cap, so the hole
rem:weightshape fell into was closed and the other seven were readable
-- and at 1.2747, 1.2853, 1.2743, 1.2096, 1.1226, 1.1175, 1.1085
elsewhere, falling monotonically toward one as N grows.  It ended by
asking what shape that convergence has.

The shape of a symptom is the wrong thing to fit when the cause is one
step away.  The constants were measured at a single N and used at
seven others; **if c(p) and kappa are themselves functions of N, the
transport failure is not a defect of the model but a fact about its
coefficients**, and it is measured by fitting them locally at every N
rather than by fitting the deviation.

That splits the question in two.  Either the local model is right
everywhere -- in which case the whole of what remains here is the
N-dependence of eight or so numbers -- or it is not, and the model is
incomplete in a way no coefficient can absorb.

BACKS: Remark {#rem:localcorrections} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  FF1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker to a relative 1e-12.
  FF2 **THE ONE THAT SPLITS IT.**  With c(p) and kappa taken at the
      same N they are used at, the corrected model reproduces the
      residue within a factor 1.05 at every one of the eight N.
  FF3 And the constants move: kappa's range across the eight N
      exceeds 0.02, so the transport failure of rem:omegafactor is
      accounted for by kappa alone moving.
  FF4 And it moves smoothly: fitting kappa on log N, the slope
      resolves at |t| above 3.

REFUTATION RULE (fixed before the run)

  FF1 REFUTED outside 1e-12; nothing below is reported.
  FF2 **REFUTED at any N outside 1.05.**  Then the model is
      incomplete in a way locally fitted coefficients cannot absorb,
      and the description rem:omegafactor gave of w is not a
      description of w at every N -- only near the one it was fitted
      at.  That is the outcome that costs the most and it must be
      stated in those words.
  FF3 REFUTED at or below a 0.02 range.  Then kappa does not move
      enough to explain rem:omegafactor's transport failure, and what
      does move is the c(p) or the model itself; the remark must then
      report the c(p) drifts rather than kappa's.
  FF4 REFUTED at |t| of 3 or below.  Then kappa moves without a
      resolved trend on these eight points, and no shape may be read
      into it -- the range in FF3 would stand as a fact and its
      direction would not.

  **A CRASH FIXED BEFORE ANY VERDICT EXISTED, DISCLOSED.**  The first
  execution died at N = 400000 with a division by zero: there
  D = 291 = 3*97 is squarefree and coprime to N, so d = D contributes,
  and the model's taper (1 - d/D) is exactly zero there.  At the
  earlier N it did not arise -- at N = 200000, D = 215 = 5*43 shares
  the factor 5 with N and is skipped.  **The model predicts w(D) = 0
  exactly and the measurement does not**, so c(D) is not a ratio that
  exists; d = D is excluded from the c fits and its measured weight is
  printed, so what the model drops is on the record rather than
  papered over.  No verdict had been produced when this was found.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  kappa at each N is
  a mean over the omega = 2 d above the noise threshold, so it carries
  its own error; **this run prints the standard error of each kappa**
  and FF3's range is to be read against the largest of them.  A range
  below that error is not a range, and FF3's verdict would then stand
  without a reading.  The threshold is the one rem:weightshape fixed
  and rem:omegafactor kept: a d whose |w(d)| is below a thousandth of
  |w(3)| is excluded.

  WHAT THIS CANNOT DO.  One radical family, eight N over 2.1 decades.
  A model that absorbs everything into locally fitted constants has
  explained nothing by itself -- it has moved the question to those
  constants, which is progress only if they are simpler than what
  they replace, and this run does not claim they are.  Nothing here
  bounds anything and rem:residuemodel's closure stands.
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
OUT = os.path.join(ROOT, "results", "audit_local_corrections.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000]
NGATE = 200_000
RELID = 1e-12
LOCALFACTOR = 1.05
KRANGE = 0.02
TCAP = 3.0
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
    if not m:
        raise SystemExit("no mainA marker for N = %d" % NGATE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the weight corrections c(p) and kappa fitted at each N",
    "           separately, the corrected model's residue with the",
    "           local constants, and how those constants move with N.",
    "FIELD: N = %s; d over the squarefree d <= D coprime to N, with"
    % NS,
    "       D = floor((N-1)/K), K = floor(N^%.2f). c(p) is measured"
    % THETA,
    "       for every prime in the range at every N. A at N = %d is"
    % NGATE,
    "       READ from results/audit_mainterm_removal.txt.",
    "MODEL: w_model(d) = (1/phi(d))(1-d/D)/(1-1/D) corrected by",
    "       prod_{p|d} c(p) times kappa^(omega-1), the description",
    "       rem:omegafactor gave; nothing here revives the reduction",
    "       rem:residuemodel closed.",
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
    say("PRINTBOUND audit_local_corrections %d %.20f" % (17, 5e-18))
    say("  theta %.2f, local factor %.2f, kappa range %.2f, |t| %.1f,"
        % (THETA, LOCALFACTOR, KRANGE, TCAP))
    say("  noise %.4f of |w(3)|" % NOISEFRAC)

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    loc = {}
    for N in NS:
        ps, D = pieces(N, lam, mu, sqf)
        A = ps[1][1]
        meas = sum(md * v for md, v in ps.values()) / A
        ws = {d: v / A for d, (md, v) in ps.items() if d != 1}
        cs = {d: ws[d] / model(d, D) for d in ws if d != D}
        atD = ws.get(D)
        noise = NOISEFRAC * abs(ws[3])
        prim = {d: c for d, c in cs.items()
                if len(factor_set(d)) == 1}
        ks = []
        for d in cs:
            fs = sorted(factor_set(d))
            if len(fs) != 2 or abs(ws[d]) < noise:
                continue
            if not all(p in prim for p in fs):
                continue
            ks.append(cs[d] / (prim[fs[0]] * prim[fs[1]]))
        ka = np.array(ks)
        kappa = float(ka.mean())
        kse = float(ka.std(ddof=1) / math.sqrt(len(ka)))
        tot = 0.0
        for d, (md, v) in ps.items():
            corr = kappa ** max(len(factor_set(d)) - 1, 0)
            for p in factor_set(d):
                corr *= prim.get(p, 1.0)
            tot += md * model(d, D) * corr
        loc[N] = (D, meas, kappa, kse, len(ka), tot / meas, prim, A)
        say("  N = %-9d D = %-5d kappa %.6f +- %.6f on %-3d d   "
            "local factor %.6f"
            % (N, D, kappa, kse, len(ka), tot / meas))
        if atD is not None:
            say("    d = D = %d contributes w = %+.8f where the model "
                "gives exactly 0;" % (D, atD))
            say("    it is excluded from the c fits and the corrected "
                "model drops it")
            say("POINT atD_%d %.8f" % (N, atD))
        say("POINT localkappa_%d %.6f" % (N, kappa))
        say("SPREAD localkappa_%d %.6f" % (N, kse))
        say("POINT localfactor_%d %.6f" % (N, tot / meas))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- FF1
    say()
    say("FF1  the gate")
    rel = abs(loc[NGATE][7] - pubA) / max(abs(pubA), 1.0)
    ff1 = rel <= RELID
    say("  A here %.17e" % loc[NGATE][7])
    say("    its %.17e   relative %.2e" % (pubA, rel))
    say("  FF1 %s   (cap: %.0e relative)"
        % ("hold" if ff1 else "REFUTED", RELID))
    if not ff1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------------- FF2
    say()
    say("FF2  is the local model right at every N?")
    ff2 = True
    for N in NS:
        f = loc[N][5]
        ok = (1.0 / LOCALFACTOR) <= f <= LOCALFACTOR
        ff2 &= ok
        say("  N = %-9d factor %.6f   %s"
            % (N, f, "ok" if ok else "OUT"))
    say("  FF2 %s   (cap: %.2f at every N)"
        % ("hold" if ff2 else "REFUTED", LOCALFACTOR))

    # ------------------------------------------------------------- FF3
    say()
    say("FF3  do the constants move?")
    kv = [loc[N][2] for N in NS]
    kerr = max(loc[N][3] for N in NS)
    rng = max(kv) - min(kv)
    ff3 = rng > KRANGE
    say("  kappa %s" % " ".join("%.4f" % k for k in kv))
    say("  range %.6f, largest standard error %.6f" % (rng, kerr))
    say("POINT kapparange %.6f" % rng)
    say("  FF3 %s   (cap: above %.2f)"
        % ("hold" if ff3 else "REFUTED", KRANGE))
    if rng < kerr:
        say("  UNRESOLVED: the range is below the largest standard "
            "error, so FF3's")
        say("  verdict stands without a reading, as the rule says")

    # ------------------------------------------------------------- FF4
    say()
    say("FF4  does kappa move smoothly?")
    x = np.array([math.log(N) for N in NS])
    b, se = fit(x, np.array(kv))
    ff4 = abs(b / se) > TCAP
    say("  slope %+.6f +- %.6f, t %+.2f" % (b, se, b / se))
    say("TSTAT localkappa_slope %.2f" % (b / se))
    say("SPREAD localkappa_slope %.6f" % se)
    say("  FF4 %s   (cap: |t| above %.1f)"
        % ("hold" if ff4 else "REFUTED", TCAP))

    say()
    say("  and the same for the three largest c(p), for the record")
    for p in (3, 7, 11):
        vals = [loc[N][6].get(p, float("nan")) for N in NS]
        if any(v != v for v in vals):
            continue
        bp, sep = fit(x, np.array(vals))
        say("    c(%d)  %s   slope %+.6f +- %.6f, t %+.2f"
            % (p, " ".join("%.4f" % v for v in vals), bp, sep,
               bp / sep))
        say("POINT cslope_%d %.6f" % (p, bp))

    say()
    say("=" * 70)
    say("FF1 %s  FF2 %s  FF3 %s  FF4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ff1, ff2, ff3, ff4)))
    say()
    if ff2 and ff3:
        say("the model is right at every N once its constants are "
            "read there, so")
        say("rem:omegafactor's transport failure was the constants "
            "moving and not")
        say("the shape being wrong. what is left of this branch's "
            "description is")
        say("the N-dependence of those constants -- which is progress "
            "only if they")
        say("are simpler than what they replace, and this run does "
            "not claim they")
        say("are.")
    elif not ff2:
        say("locally fitted constants do not absorb it. the model is "
            "incomplete in")
        say("a way no coefficient can carry, so what rem:omegafactor "
            "described is a")
        say("description of w near the N it was fitted at and not at "
            "every N.")
    else:
        say("the local model is right and kappa does not move enough "
            "to explain the")
        say("transport failure, so what moves is the c(p) or the "
            "model itself; the")
        say("c(p) slopes are printed above and that is the statement.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
