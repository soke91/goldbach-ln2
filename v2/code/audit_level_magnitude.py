# -*- coding: utf-8 -*-
r"""
The bounded level supplies the sign. Does it supply the magnitudes?

WHAT IS AT STAKE

Remark {#rem:tailpredictors} found that the sign of the largest
imbalances is a bounded-modulus object: replacing every sign by
sign P_29 on mu's own magnitudes reproduces the observed top-decile
negative share to 0.0147, and the excess over the predictor's own
marginals is resolved at 3.06 to 11.94 spreads at the larger N and
growing. It closed by naming what it does not say: the magnitude
|a_k| = (log k) T_k |I_k| has two factors, the sign is one of them, and
whether a bounded level supplies the other was not touched.

That question has an exact form, because {#rem:leanidentity} reduced
item 4(b) to one comparison. With a_k = (log k)H(N;k),

    slope/floor = (l1/l2) / G / c,   G = l1/|sum a|,

so the slope stops growing against its floor exactly when e(G) catches
e(l1/l2), measured at +0.153911 against +0.287798. Both sides of that
comparison are computable for ANY vector, not just for mu's. So the
question becomes: does a bounded-level surrogate have a smaller deficit
e(l1/l2) - e(G) than mu does? If it does, a bounded modulus cancels
better than mu relative to its own concentration, and the route is to
explain the difference. If it does not, then the bounded level supplies
the sign and nothing else, and item 4(b)'s demand is on the magnitudes
alone.

Three surrogates separate the two factors on the gain's own field:

  sign swap    |a_k| with sign P_29        -- level signs, mu sizes
  size swap    (log k)|P_29| with sign H   -- mu signs, level sizes
  full level   (log k) P_29                -- both from the level

Nothing in this repository has computed a gain for any of them.
{#rem:sievedepth} measured the mass-weighted LEAN each level gives and
{#rem:leveldemand} the residual demand, both of which are one-sided
sums; a gain is a ratio of two norms and is a different statistic.

BACKS: Remark {#rem:levelmagnitude} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control. The whole-range gain and its exponent reproduce
      results/audit_gain_split.txt to 0.01 and 0.001, so this is the
      field {#rem:gainsplit} and {#rem:leanidentity} measured.
  V2  The level's signs are enough to keep the gain. On mu's own
      magnitudes with the level's signs the exponent stays within two
      standard errors of mu's own.
  V3  The level's magnitudes are not. On mu's signs with the level's
      magnitudes the exponent departs from mu's by more than two
      standard errors.
  V4  And the level does not close the deficit. For the full level
      surrogate e(l1/l2) - e(G) is no smaller than mu's own, so a
      bounded modulus cancels no better than mu relative to its own
      concentration.

REFUTATION RULE (fixed before the run)

  V1  REFUTED at 0.01 on a gain or 0.001 on the exponent. Either would
      mean this is not that field and nothing below may be compared
      with the published exponents.
  V2  REFUTED if the sign-swap exponent departs by more than two
      standard errors. The sign would then not be the factor
      {#rem:tailpredictors} took it for, at least not for the gain.
  V3  REFUTED if the size-swap exponent stays within two standard
      errors. The level would then supply both factors, and item 4(b)
      would be a statement about a bounded modulus outright.
  V4  REFUTED if the full surrogate's deficit is smaller than mu's.
      That is the outcome worth having: a bounded-level object would be
      cancelling better than mu against its own concentration, and the
      programme would have something to explain rather than something
      to prove.

  V1 gates: without it this is not the field the identity was measured
  on.
  V2, V3 and V4 are the measurement and do not gate.

  NO NULL IS RUN and none applies to the swaps: each is a deterministic
  rearrangement of two measured vectors, and the comparison is between
  their exponents, not against a background. The coin arm for the gain
  itself was run in audit_crossk_reference.py, where random signs on
  mu's own magnitudes gave 9.94 to 12.98 times mu's gain, and that is
  the scale any statement about cancellation here is read against; the
  sign-swap surrogate is the same construction with the level's signs
  in place of random ones.
"""

import importlib.util
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
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_level_magnitude.txt")

QFIX = 29
DRAWS = 256
SEED = 20260811


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TP = module("audit_tail_predictors")
SPL = TP.SPL
NS = TP.NS
THETA = TP.THETA


def read_published():
    """the gain, its exponent, and the concentration exponent"""
    src = io.open(os.path.join(RES, "audit_gain_split.txt"),
                  encoding="utf-8").read()
    g = {}
    for m in re.finditer(r"^  N = (\d+)\s+#k = \d+\s+head \d+\s+"
                         r"G ([\d.]+)\s+head [\d.]+\s+tail [\d.]+\s+"
                         r"mass [\d.]+\s*$", src, re.M):
        g[int(m.group(1))] = float(m.group(2))
    m = re.search(r"^  whole\s+([+-][\d.]+)\s", src, re.M)
    eg = float(m.group(1))
    src2 = io.open(os.path.join(RES, "audit_lean_identity.txt"),
                   encoding="utf-8").read()
    m2 = re.search(r"l1/l2\s+([+-][\d.]+)", src2)
    ec = float(m2.group(1)) if m2 else float("nan")
    return g, eg, ec


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se


def gain(v):
    s = abs(float(v.sum()))
    return float(np.abs(v).sum()) / s if s > 0 else float("inf")


def conc(v):
    l2 = float(np.sqrt((v ** 2).sum()))
    return float(np.abs(v).sum()) / l2 if l2 > 0 else float("inf")


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubg, pubeg, pubec = read_published()
    say("read from results/audit_gain_split.txt: %d gains and the "
        "whole-range" % len(pubg))
    say("  exponent %+.6f; from results/audit_lean_identity.txt the "
        "concentration" % pubeg)
    say("  exponent %+.6f. The deficit to close is their difference, "
        "%.6f." % (pubec, pubec - pubeg))
    say("  the field and the sieve come through "
        "code/audit_tail_predictors.py,")
    say("  so the k-set is the one {#rem:headsign} and "
        "{#rem:tailpredictors} used.")
    ceil_ = THETA / 2.0
    say("  and the ceiling the concentration exponent cannot pass is "
        "theta'/2 = %.4f" % ceil_)

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    pr = SPL.primes_upto(QFIX + 1)
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    NAMES = ["mu", "sign swap", "size swap", "full level"]
    G = dict((nm, {}) for nm in NAMES)
    C = dict((nm, {}) for nm in NAMES)
    coinG, flipG, agr, nflip = {}, {}, {}, {}
    rng = np.random.default_rng(SEED)
    say()
    say("  N            #k      " + "".join("G %-10s" % nm
                                            for nm in NAMES))
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        surv = TP.survivor_masks(N, [QFIX], pr)
        ko, _Ho, Pq = TP.predictors(N, lam, mu, sqf, oddsqf, surv,
                                    [QFIX])
        assert ko.size == ks.size and bool((ko == ks).all())
        P = Pq[QFIX]
        lk = np.log(ks.astype(np.float64))
        sh = np.sign(a)
        sp = np.sign(P)
        sp = np.where(sp == 0, sh, sp)          # a zero predicts nothing
        vecs = {
            "mu": a,
            "sign swap": np.abs(a) * sp,
            "size swap": lk * np.abs(P) * sh,
            "full level": lk * P,
        }
        for nm in NAMES:
            G[nm][N] = gain(vecs[nm])
            C[nm][N] = conc(vecs[nm])

        # the two arms the sign swap has to be read against: a full
        # coin on the same magnitudes, and mu's own signs with exactly
        # as many of them flipped at random as the predictor gets wrong
        w = np.abs(a)
        nk = w.size
        ok = (sh != 0) & (np.sign(P) != 0)
        agr[N] = float((sh[ok] == np.sign(P)[ok]).mean())
        nf = int(round((1.0 - agr[N]) * nk))
        nflip[N] = nf
        gc, gf = [], []
        for _ in range(DRAWS):
            eps = rng.integers(0, 2, size=nk) * 2.0 - 1.0
            gc.append(gain(w * eps))
            f = sh.copy()
            idx = rng.choice(nk, size=nf, replace=False)
            f[idx] *= -1.0
            gf.append(gain(w * f))
        coinG[N] = float(np.median(gc))
        flipG[N] = float(np.median(gf))
        say("  %-12d %-7d " % (N, ks.size)
            + "".join("%-12.4f" % G[nm][N] for nm in NAMES))

    x = np.log(np.array(NS, dtype=np.float64))
    eG, sG, eC, sC = {}, {}, {}, {}
    for nm in NAMES:
        eG[nm], _r, sG[nm] = fit(x, np.log(
            np.array([G[nm][N] for N in NS])))
        eC[nm], _r2, sC[nm] = fit(x, np.log(
            np.array([C[nm][N] for N in NS])))

    # ------------------------------------------------------------- V1
    say()
    say("V1  the control: the gain and its exponent")
    say("  N            here      published   diff")
    worst = 0.0
    for N in NS:
        if N in pubg:
            d = abs(G["mu"][N] - pubg[N])
            worst = max(worst, d)
            say("  %-12d %-9.4f %-11.4f %.6f" % (N, G["mu"][N],
                                                 pubg[N], d))
    de = abs(eG["mu"] - pubeg)
    v1 = worst < 0.01 and de < 0.001
    say("  exponent %+.6f against the published %+.6f, diff %.6f"
        % (eG["mu"], pubeg, de))
    say("  V1 %s   (cap 0.01 on a gain, cap 0.001 on the exponent)"
        % ("hold" if v1 else "REFUTED"))

    # -------------------------------------------------------- the table
    say()
    say("  the four vectors, their gain exponent, their concentration")
    say("  exponent, and the deficit between them:")
    say("  vector       e(G)         s.e.       e(l1/l2)     s.e.       "
        "deficit")
    for nm in NAMES:
        say("  %-12s %+-12.6f %-10.6f %+-12.6f %-10.6f %+.6f"
            % (nm, eG[nm], sG[nm], eC[nm], sC[nm], eC[nm] - eG[nm]))

    # ------------------------------------------------------------- V2
    say()
    say("V2  are the level's signs enough to keep the gain?")
    d2 = abs(eG["sign swap"] - eG["mu"])
    s2 = math.hypot(sG["sign swap"], sG["mu"])
    v2 = d2 <= 2.0 * s2
    say("  sign swap %+.6f against mu's %+.6f, difference %+.6f at "
        "%.2f standard errors"
        % (eG["sign swap"], eG["mu"], eG["sign swap"] - eG["mu"],
           d2 / s2))
    say("  V2 %s   (cap 2 standard errors)"
        % ("hold" if v2 else "REFUTED"))

    # ------------------------------------------------------------- V3
    say()
    say("V3  are the level's magnitudes not?")
    d3 = abs(eG["size swap"] - eG["mu"])
    s3 = math.hypot(sG["size swap"], sG["mu"])
    v3 = d3 > 2.0 * s3
    say("  size swap %+.6f against mu's %+.6f, difference %+.6f at "
        "%.2f standard errors"
        % (eG["size swap"], eG["mu"], eG["size swap"] - eG["mu"],
           d3 / s3))
    say("  V3 %s   (cap 2 standard errors)"
        % ("hold" if v3 else "REFUTED"))

    # ------------------------------------------------------------- V4
    say()
    say("V4  does the level close the deficit?")
    dmu = eC["mu"] - eG["mu"]
    dfl = eC["full level"] - eG["full level"]
    v4 = dfl >= dmu
    say("  mu's deficit         %+.6f" % dmu)
    say("  the level's deficit  %+.6f" % dfl)
    say("  the level's is %s" % ("no smaller" if v4 else "SMALLER"))
    say("  V4 %s" % ("hold" if v4 else "REFUTED"))
    if v4:
        say("  so a bounded modulus cancels no better than mu against "
            "its own")
        say("  concentration. The level supplies the sign and not the")
        say("  cancellation, and item 4(b)'s demand is on the "
            "magnitudes.")
    else:
        say("  so a bounded-level object cancels BETTER than mu against "
            "its")
        say("  own concentration. That is a thing to explain rather "
            "than a")
        say("  thing to prove, and it is where the programme should "
            "look.")

    say()
    say("  and against the ceiling, which is what the deficit has to "
        "reach zero")
    say("  under: e(l1/l2) can be at most %.4f, so a vector whose "
        "concentration" % ceil_)
    say("  exponent is far below it has an easier deficit to close and "
        "a weaker")
    say("  statement to make. Each vector's distance from the ceiling:")
    for nm in NAMES:
        say("  %-12s %+.6f" % (nm, eC[nm] - ceil_))

    # ------------------------------------------- not pre-registered
    say()
    say("X1  is the sign swap a bounded level, or is it a coin?")
    say("  (written after V2 fell; not pre-registered). A random sign "
        "vector")
    say("  gives |sum| of order l2, so its gain is l1/l2 up to a "
        "constant and")
    say("  its deficit is zero BY CONSTRUCTION. A small deficit is "
        "therefore")
    say("  not evidence of a good predictor -- it is what noise looks "
        "like, and")
    say("  {#rem:whycoinwins} settled that the coin is a competitor "
        "rather than")
    say("  a null. Two arms decide it: a full coin on mu's magnitudes, "
        "and mu's")
    say("  own signs with exactly as many flipped at random as the "
        "predictor")
    say("  gets wrong. %d draws at seed %d, medians." % (DRAWS, SEED))
    say("  N            agree     flips  G sign swap  G coin    "
        "G matched flip")
    for N in NS:
        say("  %-12d %-9.4f %-6d %-12.4f %-9.4f %.4f"
            % (N, agr[N], nflip[N], G["sign swap"][N], coinG[N],
               flipG[N]))
    ec_, _rc, sc_ = fit(x, np.log(np.array([coinG[N] for N in NS])))
    ef_, _rf, sf_ = fit(x, np.log(np.array([flipG[N] for N in NS])))
    say("  vector           e(G)         s.e.       deficit against "
        "e(l1/l2)")
    say("  sign swap        %+-12.6f %-10.6f %+.6f"
        % (eG["sign swap"], sG["sign swap"],
           eC["sign swap"] - eG["sign swap"]))
    say("  coin             %+-12.6f %-10.6f %+.6f"
        % (ec_, sc_, eC["mu"] - ec_))
    say("  matched flip     %+-12.6f %-10.6f %+.6f"
        % (ef_, sf_, eC["mu"] - ef_))
    dsw = abs(eG["sign swap"] - ef_) / math.hypot(sG["sign swap"], sf_)
    say("  the sign swap and the matched flip differ by %+.6f, %.2f "
        "standard errors"
        % (eG["sign swap"] - ef_, dsw))
    say("TSTAT slope_levelmag_swapvsflip %.2f" % dsw)
    say("SPREAD slope_levelmag_swapvsflip %.4f"
        % float(x.max() - x.min()))
    if dsw < 2.0:
        say("UNRESOLVED SIGN slope_levelmag_swapvsflip")
        say("  so the sign swap is NOT distinguishable from mu's own "
            "signs with")
        say("  the predictor's error rate applied at random. Its small")
        say("  deficit is the error rate behaving like noise, not a "
            "bounded")
        say("  level cancelling. V2's refutation carries no hopeful "
            "reading.")
    else:
        say("  so the sign swap is separated from the matched flip, and "
            "the")
        say("  direction of that separation says whether the level's "
            "signs")
        say("  cancel better or worse than its error rate alone would.")

    say()
    say("=" * 70)
    say("V1 %s  V2 %s  V3 %s  V4 %s"
        % tuple("hold" if v else "REFUTED" for v in (v1, v2, v3, v4)))

    head = [
        "STATISTIC: on the squarefree k < N^theta' coprime to N, the gain",
        "           G = l1/|sum| and the concentration l1/l2 of four",
        "           vectors: a_k = (log k)H(N;k) with m over all m, the",
        "           gain's own convention; |a_k| carrying sign P_29;",
        "           (log k)|P_29| carrying sign H; and (log k)P_29, with",
        "           P_Q = sum_m mu(m) [N - mk has no odd prime factor at",
        "           or below Q] over the odd squarefree m. For each, the",
        "           least-squares exponent of each norm ratio against",
        "           log N with its standard error, and the deficit",
        "           e(l1/l2) - e(G) that {#rem:leanidentity} requires to",
        "           reach zero.",
        "NULL: none is run and none applies to the swaps: each is a",
        "      deterministic rearrangement of two measured vectors and",
        "      the comparison is between their exponents, not against a",
        "      background. The coin arm for the gain was run in",
        "      audit_crossk_reference.py, where random signs on mu's own",
        "      magnitudes gave 9.94 to 12.98 times mu's gain; the",
        "      sign-swap vector is that construction with the level's",
        "      signs in place of random ones.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta'; for",
        "       H the m run over ALL m < N/k coprime to k, the convention",
        "       of code/audit_gain_split.py, and for P_29 over the odd",
        "       squarefree m, the convention of code/audit_sieve_depth.py",
        "       -- {#rem:tailpredictors} measured that the two give the",
        "       same sign at every k; Lambda and mu from an integer sieve",
        "       to " + str(NMAX) + "; the field, the sieve and the",
        "       predictor come through code/audit_tail_predictors.py; the",
        "       published gains and exponent are read from",
        "       results/audit_gain_split.txt and the concentration",
        "       exponent from results/audit_lean_identity.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not v1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
