# -*- coding: utf-8 -*-
r"""
The predictors that are not degenerate, tested where the tail is.

WHAT IS AT STAKE

Remark {#rem:tailmertens} eliminated the Mertens predictor as an
account of the tail, and for a specific reason: on k < N^theta' it takes
one sign on 0.9829 to 0.9970 of the k and on every k of the top decile
by |I_k|, so its agreement there is sign H's own marginal rate read
back. Gate check G73 then made every agreement claim in the repository
publish that share, and the survey said the fault is narrow. Two files
are degenerate; the sieve predictors are not, declaring 0.5801, 0.5879,
0.6686, 0.7091 and 0.7726. Those have real predictive power on the
whole range -- audit_sieve_depth.py measures sign agreement 0.8129 down
to 0.7579 at the fixed level Q = 29 -- and **none of them has ever been
looked at on the tail.**

There is a second thing nobody has checked, and it has to come first.
Every sign predictor in this repository is built with m ODD:
audit_sieve_depth.py, audit_survivor_range.py and
audit_logweight_predictor.py all sum over m = 1, 3, 5, ... The gain that
item 4(b) is about is not: audit_gain_split.py, and so
{#rem:headsign}'s split and every |I_k| in it, sums over ALL m. Those
are different fields on their face. They should nearly coincide, because
N is even and k coprime to N is odd, so an even m makes N - mk even and
Lambda(N - mk) vanishes unless N - mk is a power of two -- but "should"
is not a measurement, and if the two signs part company anywhere that
matters, then the predictors were never predictors for the gain's field
and the last several remarks have been comparing across a boundary.

BACKS: Remark {#rem:tailpredictors} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The two conventions are one field, and this is the same
      predictor. sign H summed over all m agrees with sign H summed
      over the odd m at more than 0.999 of the k at every N, and the
      agreement of sign H with sign P_29 reproduces
      results/audit_sieve_depth.txt to 0.01 on the N the two sweeps
      share.
  S2  The predictor keeps its variance where it is needed. Its
      majority sign share on the top decile by |I_k| stays below 0.9 at
      every N -- unlike the Mertens predictor, which is at 1.0000
      there.
  S3  And it reaches the tail. Replacing every sign by sign P_29 on
      mu's own magnitudes reproduces the observed top-decile negative
      share to within 0.05 at every N.
  S4  And it is the predictor and not the marginals. The tail
      agreement exceeds what signs drawn independently at the
      predictor's own negative rate give, by more than three of that
      arm's spreads, at every N.

REFUTATION RULE (fixed before the run)

  S1  REFUTED at 0.999 on the field agreement or 0.01 on the published
      agreement. The first would mean the repository's predictors are
      not predictors for the gain's field and the comparison cannot be
      made at all; the second would mean this is not that predictor.
  S2  REFUTED if the tail majority share reaches 0.9 at any N. The
      predictor would then be degenerate exactly where it is being
      asked to work, like the Mertens one, and the test would be void
      before it started.
  S3  REFUTED if the predicted share misses the observed by more than
      0.05 at any N. That is the outcome worth having in the negative
      direction: the last non-degenerate candidate would have been
      tested on the tail and failed, and item 4(b) would have no
      computational candidate left.
  S4  REFUTED if the tail agreement fails to beat the matched-marginal
      arm by three spreads at any N. Whatever agreement there is would
      then be the marginals again, in a subtler form than
      {#rem:tailmertens} found.

  S1 gates: without it the fields are not the same field.
  S2, S3 and S4 are the measurement and do not gate.

  THE NULL is S4's matched-marginal arm. Signs are drawn independently
  at the predictor's OWN negative rate on the tail set, so the baseline
  already contains everything the predictor's marginal distribution can
  explain, and only the pairing between k and its prediction is
  destroyed. It is closed form on a fixed set: the tail is fixed once
  the magnitudes are given, so the arm is binomial and needs no draws.
  The permutation of {#rem:leanmertens} is unavailable here for the
  reason audit_sieve_depth.py gives -- P_Q is not a function of
  floor(N/k), so there is no coarse label to shuffle among.
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
OUT = os.path.join(RES, "audit_tail_predictors.txt")

BLOCKS = 10
QS = [29, 211]
QFIX = 29


def head_module():
    """the split of {#rem:headsign}, imported so it cannot drift"""
    p = os.path.join(CODE, "audit_head_sign.py")
    spec = importlib.util.spec_from_file_location("audit_head_sign", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HS = head_module()
SPL = HS.SPL
NS = HS.NS
THETA = HS.THETA


def read_published():
    """the sieve-depth agreements and the head split's tail share"""
    src = io.open(os.path.join(RES, "audit_sieve_depth.txt"),
                  encoding="utf-8").read()
    i = src.index("N            Q=29 ")
    ag = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 2 or not f[0].isdigit():
            break
        ag[int(f[0])] = float(f[1])
    src2 = io.open(os.path.join(RES, "audit_head_sign.txt"),
                   encoding="utf-8").read()
    m = re.search(r"top decile by \|I\| runs ([\d.]+) to ([\d.]+)", src2)
    return ag, (float(m.group(1)), float(m.group(2)))


def survivor_masks(N, qs, pr):
    """for each Q, the v <= N with no odd prime factor at or below Q"""
    out = {}
    for q in qs:
        s = np.ones(N + 1, dtype=bool)
        s[0] = False
        for p in pr:
            p = int(p)
            if p > q:
                break
            if p == 2:
                continue
            s[p::p] = False
        out[q] = s
    return out


def predictors(N, lam, mu, sqf, oddsqf, surv, qs):
    """H over the odd m, and P_Q, on exactly signed_parts' k-set

    The k-set and the inner bound are audit_head_sign.signed_parts'
    line for line, so the two runs are indexed the same and S1 can
    compare them term by term; only the m-parity and the weight
    differ.
    """
    PN = SPL.factor_set(N)
    ks = np.array([k for k in range(2, int(N ** THETA))
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    Ho = []
    Ps = dict((q, []) for q in qs)
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[oddsqf[ms]]
        for q in SPL.factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            Ho.append(0.0)
            for q in qs:
                Ps[q].append(0.0)
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        Ho.append(float((lam[vals] * g).sum()))
        for q in qs:
            Ps[q].append(float(g[surv[q][vals]].sum()))
    return (ks, np.array(Ho),
            dict((q, np.array(Ps[q])) for q in qs))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    MAJS = []
    pubag, pubtop = read_published()
    say("read from results/audit_sieve_depth.txt: %d agreements at "
        "Q = %d, %.4f to %.4f"
        % (len(pubag), QFIX, min(pubag.values()), max(pubag.values())))
    say("  from results/audit_head_sign.txt the top decile by |I| "
        "running %.4f to %.4f" % pubtop)
    say("  the all-m split is imported from code/audit_head_sign.py; "
        "the odd-m")
    say("  convention and the sieve predictor are rebuilt here on the "
        "same k-set.")

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    pr = SPL.primes_upto(max(QS) + 1)
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    obs_top, fieldag, ag_all, ag_tail = {}, {}, {}, {}
    maj_all, maj_tail, pred_top, arm_mu, arm_sd = {}, {}, {}, {}, {}
    say()
    say("  N            #k      field agree  obs top   agree all  "
        "agree tail  maj tail  predicted")
    for N in NS:
        ks, P, M, C = HS.signed_parts(N, lam, mu, sqf)
        T = P + M
        I = np.where(T > 0, (P - M) / np.where(T > 0, T, 1.0), 0.0)
        surv = survivor_masks(N, QS, pr)
        ko, Ho, Pq = predictors(N, lam, mu, sqf, oddsqf, surv, QS)
        assert ko.size == ks.size and bool((ko == ks).all())
        n = I.size

        sh = np.sign(I)
        sho = np.sign(Ho)
        both = (sh != 0) & (sho != 0)
        fieldag[N] = float((sh[both] == sho[both]).mean())

        order = np.argsort(-np.abs(I))
        nt = max(1, int(round(n / BLOCKS)))
        tail = order[:nt]
        obs_top[N] = float((I[tail] < 0).mean())

        sp = np.sign(Pq[QFIX])
        ok = (sh != 0) & (sp != 0)
        ag_all[N] = float((sh[ok] == sp[ok]).mean())
        tk = tail[ok[tail]]
        ag_tail[N] = float((sh[tk] == sp[tk]).mean())
        maj_all[N] = max(float((sp[ok] > 0).mean()),
                         float((sp[ok] < 0).mean()))
        maj_tail[N] = max(float((sp[tk] > 0).mean()),
                          float((sp[tk] < 0).mean()))
        MAJS += [maj_all[N], maj_tail[N]]
        for q in QS:
            s2 = np.sign(Pq[q])
            o2 = (sh != 0) & (s2 != 0)
            MAJS.append(max(float((s2[o2] > 0).mean()),
                            float((s2[o2] < 0).mean())))

        Ip = np.abs(I) * np.where(sp != 0, sp, 1.0)
        op = np.argsort(-np.abs(Ip))[:nt]
        pred_top[N] = float((Ip[op] < 0).mean())

        pneg = float((sp[tk] < 0).mean())
        pobs = float((sh[tk] < 0).mean())
        p = pneg * pobs + (1.0 - pneg) * (1.0 - pobs)
        arm_mu[N] = p
        arm_sd[N] = math.sqrt(max(p * (1.0 - p), 0.0) / max(tk.size, 1))

        say("  %-12d %-7d %-12.6f %-9.4f %-10.4f %-11.4f %-9.4f %.4f"
            % (N, n, fieldag[N], obs_top[N], ag_all[N], ag_tail[N],
               maj_tail[N], pred_top[N]))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- S1
    say()
    say("S1  are the two conventions one field, and is this that "
        "predictor?")
    fmin = min(fieldag.values())
    shared = [N for N in NS if N in pubag]
    dpub = max(abs(ag_all[N] - pubag[N]) for N in shared)
    s1 = fmin > 0.999 and dpub < 0.01
    say("  sign H over all m against sign H over the odd m: %.6f to "
        "%.6f agree" % (fmin, max(fieldag.values())))
    say("  N is even and k coprime to N is odd, so an even m makes "
        "N - mk even")
    say("  and Lambda vanishes unless N - mk is a power of two -- the "
        "measured")
    say("  agreement is what that argument is worth.")
    say("  the Q = %d agreement against the published, over %d shared "
        "N: worst %.6f" % (QFIX, len(shared), dpub))
    say("  S1 %s   (cap 0.999 on the field, cap 0.01 on the published)"
        % ("hold" if s1 else "REFUTED"))

    # ------------------------------------------------------------- S2
    say()
    say("S2  does the predictor keep its variance on the tail?")
    say("  N            maj all   maj tail  Mertens maj tail "
        "(published)")
    s2 = True
    for N in NS:
        if maj_tail[N] >= 0.9:
            s2 = False
        say("  %-12d %-9.4f %-9.4f %s" % (N, maj_all[N], maj_tail[N],
                                          "1.0000"))
    say("  S2 %s   (cap 0.9)" % ("hold" if s2 else "REFUTED"))

    # ------------------------------------------------------------- S3
    say()
    say("S3  does it reach the tail?")
    say("  N            observed  predicted  difference")
    s3 = True
    for N in NS:
        d = pred_top[N] - obs_top[N]
        if abs(d) > 0.05:
            s3 = False
        say("  %-12d %-9.4f %-10.4f %+.4f" % (N, obs_top[N],
                                              pred_top[N], d))
    say("  S3 %s   (cap 0.05)" % ("hold" if s3 else "REFUTED"))

    # ------------------------------------------------------------- S4
    say()
    say("S4  is it the predictor or the marginals?")
    say("  N            agree tail  matched arm  spread   in spreads")
    s4 = True
    for N in NS:
        g = ((ag_tail[N] - arm_mu[N]) / arm_sd[N]
             if arm_sd[N] > 0 else float("inf"))
        if g <= 3.0:
            s4 = False
        say("  %-12d %-11.4f %-12.4f %-8.4f %+.2f"
            % (N, ag_tail[N], arm_mu[N], arm_sd[N], g))
    say("  S4 %s   (cap 3 spreads)" % ("hold" if s4 else "REFUTED"))

    say()
    et, _rt, set_ = fit_(x, np.log(np.array([ag_tail[N] for N in NS])))
    ea, _ra, sea = fit_(x, np.log(np.array([ag_all[N] for N in NS])))
    say("  and where the two agreements are going:")
    say("  quantity      exponent     s.e.       t")
    say("  agree tail    %+-12.6f %-10.6f %.2f" % (et, set_,
                                                   abs(et) / set_))
    say("  agree all     %+-12.6f %-10.6f %.2f" % (ea, sea,
                                                   abs(ea) / sea))
    say("TSTAT slope_tailpredictors_tail %.2f" % (abs(et) / set_))
    say("SPREAD slope_tailpredictors_tail %.4f"
        % float(x.max() - x.min()))
    if abs(et) / set_ < 2.0:
        say("UNRESOLVED SIGN slope_tailpredictors_tail")
    say("TSTAT slope_tailpredictors_all %.2f" % (abs(ea) / sea))
    say("SPREAD slope_tailpredictors_all %.4f"
        % float(x.max() - x.min()))
    if abs(ea) / sea < 2.0:
        say("UNRESOLVED SIGN slope_tailpredictors_all")

    mj = max(MAJS)
    say()
    say("  the predictor's own majority sign share, at its worst over "
        "everything")
    say("  reported above: %.4f. An agreement is only a measurement "
        "where the" % mj)
    say("  predictor has variance.")
    say("MARGINAL audit_tail_predictors %.4f" % mj)
    if mj >= 0.9:
        say("DEGENERATE audit_tail_predictors")

    say()
    say("=" * 70)
    say("S1 %s  S2 %s  S3 %s  S4 %s"
        % tuple("hold" if v else "REFUTED" for v in (s1, s2, s3, s4)))

    head = [
        "STATISTIC: on the squarefree k < N^theta' coprime to N, the",
        "           sign agreement of H(N;k) -- summed over ALL m, the",
        "           gain's convention -- with sign P_Q for",
        "           P_Q = sum_m mu(m) [N - mk has no odd prime factor at",
        "           or below Q], Q = 29 and 211, on the whole range and",
        "           on the top decile of the k by the imbalance |I_k| of",
        "           {#rem:headsign}; the predictor's own majority sign",
        "           share on each; the negative share of that top decile",
        "           when every sign is replaced by the predictor's on",
        "           mu's own magnitudes; the same agreement expected from",
        "           signs drawn independently at the predictor's own",
        "           negative rate; and, as the control, the fraction of k",
        "           at which sign H over all m equals sign H over the odd",
        "           m.",
        "NULL: the matched-marginal arm of S4. Signs are drawn",
        "      independently at the predictor's OWN negative rate on the",
        "      tail set, so the baseline contains everything the",
        "      predictor's marginal distribution can explain and only the",
        "      pairing between k and its prediction is destroyed; the",
        "      tail is a fixed set once the magnitudes are given, so the",
        "      arm is binomial and needs no draws. The permutation of",
        "      {#rem:leanmertens} is unavailable because P_Q is not a",
        "      function of floor(N/k).",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta'; for",
        "       I_k and for sign H the m run over ALL m < N/k coprime to",
        "       k, the convention of code/audit_gain_split.py, and for",
        "       the odd-m control and for P_Q over the odd squarefree m,",
        "       the convention of code/audit_sieve_depth.py; Lambda and mu",
        "       from an integer sieve to " + str(NMAX) + "; the split and",
        "       the sieve are imported from code/audit_head_sign.py; the",
        "       published agreements are read from",
        "       results/audit_sieve_depth.txt and the published tail",
        "       share from results/audit_head_sign.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not s1:
        raise SystemExit(1)
    return 0


def fit_(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se


if __name__ == "__main__":
    sys.exit(main())
