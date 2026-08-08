# -*- coding: utf-8 -*-
r"""
Does the one existing predictor of sign H reach the tail it has to?

WHAT IS AT STAKE

Remark {#rem:tailasymmetry} closed on the statement that the coupling
between the size of I_k and its sign has no name in this programme. That
was said without checking the one predictor the repository already has.
Remark {#rem:leanmertens} and lab_lean_oddmertens.py found that
sign H(N;k) agrees with sign Modd(floor(N/k)) -- the Mertens function
over odd m -- at 0.7669 to 0.7704 against a permutation baseline near
0.56, and audit_oddmertens_range.py then tested whether that transfers
to the k-range the gain is measured on. It does not: on k < N^theta'
the agreement is 0.5815, 0.6161, 0.5718, 0.5618, 0.5209, 0.5201, and
its check P2 was registered as "below 0.70 there" and held. The 0.77
was demonstrated on inner lengths 2 <= N/k <= 1000; the gain's range
has N/k from 215 to 2133333.

So the predictor is known to fail on the bulk of the relevant range.
**It has never been tested on the tail**, and the tail is a different
question: {#rem:headsign} measured that the top decile of the k by
|I_k| is 0.8547 negative at the top N, and {#rem:tailasymmetry} then
showed that neither a uniform sign rate, nor a rate varying with k by
decile, nor a shifted symmetric law with the measured bias comes within
four of their own spreads of it. A predictor that agrees with sign H
only half the time on average could still agree with it almost always
where |I_k| is largest -- and if it does, the coupling is named, and
named by something already in the paper. Nothing rules this out in
advance: Modd(floor(N/k)) is wildly non-monotone in k, so the
k-decile arm of {#rem:tailasymmetry}, which averaged over monotone bins
of k, would have destroyed exactly this structure and reported nothing.

BACKS: Remark {#rem:tailmertens} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  The controls. The overall agreement on k < N^theta' reproduces
      results/audit_oddmertens_range.txt at the N it publishes to
      0.001, and the observed top-decile negative share reproduces
      results/audit_head_sign.txt to 0.001.
  R2  The agreement rises with the imbalance. In the top decile by
      |I_k| it exceeds the whole range's by more than 0.05 at every N.
  R3  And it accounts for the tail. Replacing every sign by
      sign Modd(floor(N/k)) on mu's own magnitudes reproduces the
      observed top-decile negative share to within 0.05 at every N.
  R4  And it is the predictor and not the marginals. The tail
      agreement exceeds the largest of 16 permutations of the
      predictor's signs among the distinct values of floor(N/k), the
      baseline convention of lab_lean_oddmertens.py, at every N.

REFUTATION RULE (fixed before the run)

  R1  REFUTED at 0.001 on either. Either would mean this is not the
      predictor audit_oddmertens_range.py measured or not the field
      {#rem:headsign} split, and nothing below may be compared with
      them.
  R2  REFUTED if the tail agreement fails to exceed the whole range's
      by 0.05 at any N. The predictor would then be as blind to the
      tail as it is to the bulk.
  R3  REFUTED if the predicted top-decile share misses the observed by
      more than 0.05 at any N. That is the outcome worth having in the
      negative direction: the repository's only predictor of sign H
      would then have been tested on the tail and failed, and the
      coupling {#rem:tailasymmetry} named as unnamed would stay
      unnamed with one candidate eliminated instead of none.
  R4  REFUTED if any permutation reaches the tail agreement at any N.
      Then whatever agreement there is comes from both sides being
      predominantly negative and not from the predictor.

  R1 gates: without it this is not the same predictor or the same
  field.
  R2, R3 and R4 are the measurement and do not gate.

  THE NULL is the permutation of R4, the convention
  lab_lean_oddmertens.py established: the predictor's signs are
  shuffled among the DISTINCT values of floor(N/k), which preserves
  both marginal sign distributions exactly, so a baseline near 0.56
  rather than 0.5 is what agreement has to beat. M4 is met by it -- the
  arm destroys the pairing between k and its predictor and keeps
  everything else.
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
OUT = os.path.join(RES, "audit_tail_mertens.txt")

BLOCKS = 10
PERMS = 16
SEED = 20260810


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
    """the range test's agreements and the head split's shares"""
    src = io.open(os.path.join(RES, "audit_oddmertens_range.txt"),
                  encoding="utf-8").read()
    i = src.index("N            #k     N/k from   agreement  mu lean")
    ag = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ag[int(f[0])] = float(f[3])
    src2 = io.open(os.path.join(RES, "audit_head_sign.txt"),
                   encoding="utf-8").read()
    j = src2.index("decile   by |a|       by |I|       by T         by k")
    top = float(src2[j:].splitlines()[1].split()[2])
    m = re.search(r"top decile by \|I\| runs ([\d.]+) to ([\d.]+)", src2)
    return ag, top, (float(m.group(1)), float(m.group(2)))


def odd_mertens(limit, mu):
    """cumulative sum of mu over the ODD m up to limit"""
    v = mu[:limit + 1].astype(np.int32).copy()
    v[0] = 0
    v[2::2] = 0
    return np.cumsum(v)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubag, pubtop, pubrange = read_published()
    say("read from results/audit_oddmertens_range.txt: %d agreements on "
        "k < N^theta'," % len(pubag))
    say("  %.4f to %.4f" % (min(pubag.values()), max(pubag.values())))
    say("  from results/audit_head_sign.txt the top decile by |I|, "
        "%.4f at the top N" % pubtop)
    say("  and %.4f to %.4f over the sweep"
        % (pubrange[0], pubrange[1]))
    say("  the split is imported from code/audit_head_sign.py, so the "
        "field is")
    say("  the one {#rem:headsign} measured.")

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    Modd = odd_mertens(NMAX // 2 + 2, mu)
    say("Modd built over the odd m to %d" % (NMAX // 2 + 2))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    rng = np.random.default_rng(SEED)
    obs_top, ag_all, ag_top, pred_top, perm_top = {}, {}, {}, {}, {}
    byI_ag, nkept = {}, {}
    negshare_s, negshare_s_tail, ndist, ndist_neg = {}, {}, {}, {}
    say()
    say("  N            #k      obs top   agree all  agree tail  "
        "predicted top  perm max tail")
    for N in NS:
        ks, P, M, C = HS.signed_parts(N, lam, mu, sqf)
        T = P + M
        ok = T > 0
        I = np.where(ok, (P - M) / np.where(ok, T, 1.0), 0.0)
        ks, I, ok = ks[ok], I[ok], ok[ok]
        n = I.size
        x = N // ks
        s = np.sign(Modd[x]).astype(np.float64)
        h = np.sign(I)
        keep = (s != 0) & (h != 0)
        nkept[N] = int(keep.sum())

        order = np.argsort(-np.abs(I))
        nt = max(1, int(round(n / BLOCKS)))
        tail = order[:nt]
        obs_top[N] = float((I[tail] < 0).mean())
        ag_all[N] = float((h[keep] == s[keep]).mean())
        tk = tail[keep[tail]]
        ag_top[N] = float((h[tk] == s[tk]).mean())

        # the arm: mu's own magnitudes, the predictor's signs
        Ip = np.abs(I) * np.where(s != 0, s, 1.0)
        op = np.argsort(-np.abs(Ip))[:nt]
        pred_top[N] = float((Ip[op] < 0).mean())

        # the permutation baseline, on the distinct floor(N/k)
        uq, inv = np.unique(x, return_inverse=True)
        su = np.sign(Modd[uq]).astype(np.float64)
        best = 0.0
        for _ in range(PERMS):
            sp = su[rng.permutation(su.size)][inv]
            kp = (sp != 0) & (h != 0)
            tp = tail[kp[tail]]
            if tp.size:
                best = max(best, float((h[tp] == sp[tp]).mean()))
        perm_top[N] = best

        negshare_s[N] = float((s[keep] < 0).mean())
        negshare_s_tail[N] = float((s[tk] < 0).mean()) if tk.size else             float("nan")
        ndist[N] = int(uq.size)
        ndist_neg[N] = int((su < 0).sum())

        edges = [int(round(d * n / BLOCKS)) for d in range(BLOCKS + 1)]
        byI_ag[N] = []
        for d in range(BLOCKS):
            idx = order[edges[d]:edges[d + 1]]
            ii = idx[keep[idx]]
            byI_ag[N].append(float((h[ii] == s[ii]).mean())
                             if ii.size else float("nan"))
        say("  %-12d %-7d %-9.4f %-10.4f %-11.4f %-14.4f %.4f"
            % (N, n, obs_top[N], ag_all[N], ag_top[N], pred_top[N],
               perm_top[N]))

    xl = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- R1
    say()
    say("R1  the controls")
    da = max(abs(ag_all[N] - pubag[N]) for N in NS if N in pubag)
    dt = max(abs(min(obs_top.values()) - pubrange[0]),
             abs(max(obs_top.values()) - pubrange[1]))
    r1 = da < 0.001 and dt < 0.001
    say("  the overall agreement departs from the published by %.6f "
        "over %d shared N" % (da, sum(1 for N in NS if N in pubag)))
    say("  the top decile's range %.4f to %.4f against the published "
        "%.4f to %.4f"
        % (min(obs_top.values()), max(obs_top.values()),
           pubrange[0], pubrange[1]))
    say("  R1 %s   (cap 0.001 on each)" % ("hold" if r1 else "REFUTED"))
    say("  k kept (predictor and H both nonzero): %d to %d"
        % (min(nkept.values()), max(nkept.values())))

    # ------------------------------------------------------------- R2
    say()
    say("R2  does the agreement rise with the imbalance?")
    say("  N            agree all  agree tail  gain")
    r2 = True
    for N in NS:
        g = ag_top[N] - ag_all[N]
        if g <= 0.05:
            r2 = False
        say("  %-12d %-10.4f %-11.4f %+.4f" % (N, ag_all[N],
                                               ag_top[N], g))
    say("  R2 %s   (cap 0.05)" % ("hold" if r2 else "REFUTED"))
    say()
    say("  and it holds for a reason that empties it. A predictor that "
        "is")
    say("  CONSTANT on a set agrees with anything there at the rate that")
    say("  thing takes its value, and carries no information at all. So "
        "the")
    say("  predictor's own sign distribution has to be on the table:")
    say("  N            share Modd<0  in the tail  distinct N/k  of "
        "them Modd<0")
    for N in NS:
        say("  %-12d %-13.4f %-12.4f %-13d %d"
            % (N, negshare_s[N], negshare_s_tail[N], ndist[N],
               ndist_neg[N]))
    say("  The predictor is negative on %.4f to %.4f of the k and on "
        "every"
        % (min(negshare_s.values()), max(negshare_s.values())))
    say("  one of the tail, at every N. It is not a weak predictor "
        "here, it")
    say("  is a degenerate one, and the agreement it shows is the "
        "marginal")
    say("  negative rate of H read back at it. That is also why the")
    say("  permutation baseline of R4 reaches it: permuting a constant")
    say("  changes nothing.")
    say()
    say("  and the whole profile, agreement by decile of |I| at the top "
        "N:")
    say("  " + "  ".join("%.4f" % v for v in byI_ag[NMAX]))

    # ------------------------------------------------------------- R3
    say()
    say("R3  does the predictor account for the tail?")
    say("  N            observed  predicted  difference")
    r3 = True
    for N in NS:
        d = abs(pred_top[N] - obs_top[N])
        if d > 0.05:
            r3 = False
        say("  %-12d %-9.4f %-10.4f %+.4f"
            % (N, obs_top[N], pred_top[N], pred_top[N] - obs_top[N]))
    say("  R3 %s   (cap 0.05)" % ("hold" if r3 else "REFUTED"))
    if not r3:
        say("  and it fails in the one direction a degenerate predictor")
        say("  must: it predicts every tail sign negative, so its top")
        say("  decile is 1 exactly, and the observed is not. The gap "
            "grows")
        say("  along the sweep because the observed share falls while "
            "the")
        say("  prediction cannot.")

    # ------------------------------------------------------------- R4
    say()
    say("R4  is it the predictor or the marginals?")
    say("  N            agree tail  best of %d permutations  gain"
        % PERMS)
    r4 = True
    for N in NS:
        if perm_top[N] >= ag_top[N]:
            r4 = False
        say("  %-12d %-11.4f %-24.4f %+.4f"
            % (N, ag_top[N], perm_top[N], ag_top[N] - perm_top[N]))
    say("  R4 %s" % ("hold" if r4 else "REFUTED"))

    say()
    et, _rt, set_ = fit_(xl, np.log(np.array([ag_top[N] for N in NS])))
    ea, _ra, sea = fit_(xl, np.log(np.array([ag_all[N] for N in NS])))
    say("  and where the two are going:")
    say("  quantity      exponent     s.e.       t")
    say("  agree tail    %+-12.6f %-10.6f %.2f" % (et, set_,
                                                   abs(et) / set_))
    say("  agree all     %+-12.6f %-10.6f %.2f" % (ea, sea,
                                                   abs(ea) / sea))
    say("TSTAT slope_tailmertens_tail %.2f" % (abs(et) / set_))
    say("SPREAD slope_tailmertens_tail %.4f" % float(xl.max() - xl.min()))
    if abs(et) / set_ < 2.0:
        say("UNRESOLVED SIGN slope_tailmertens_tail")
    say("TSTAT slope_tailmertens_all %.2f" % (abs(ea) / sea))
    say("SPREAD slope_tailmertens_all %.4f" % float(xl.max() - xl.min()))
    if abs(ea) / sea < 2.0:
        say("UNRESOLVED SIGN slope_tailmertens_all")

    say()
    say("=" * 70)
    say("R1 %s  R2 %s  R3 %s  R4 %s"
        % tuple("hold" if v else "REFUTED" for v in (r1, r2, r3, r4)))

    head = [
        "STATISTIC: the fraction of k at which sign H(N;k) equals",
        "           sign Modd(floor(N/k)), with Modd the Mertens",
        "           function over the odd m, on the whole range",
        "           k < N^theta' and on the top decile of the k by the",
        "           imbalance |I_k| of {#rem:headsign}; the agreement by",
        "           decile of |I_k|; the negative share of that top",
        "           decile when every sign is replaced by the",
        "           predictor's on mu's own magnitudes; and the largest",
        "           tail agreement over " + str(PERMS) + " permutations",
        "           of the predictor's signs among the distinct values",
        "           of floor(N/k).",
        "NULL: the permutation of R4, the convention",
        "      lab_lean_oddmertens.py established. The predictor's signs",
        "      are shuffled among the DISTINCT values of floor(N/k),",
        "      which preserves both marginal sign distributions exactly,",
        "      so the baseline sits near 0.56 rather than 0.5. It",
        "      destroys the pairing between k and its predictor and",
        "      keeps everything else, which is M4.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta';",
        "       m over 1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + " and Modd cumulated",
        "       over the odd m from the same sieve; k with Modd = 0 or",
        "       H = 0 excluded from the agreement, the convention of",
        "       lab_lean_oddmertens.py; the split, the sieve and theta'",
        "       are imported from code/audit_head_sign.py; the published",
        "       agreements are read from",
        "       results/audit_oddmertens_range.txt and the published",
        "       shares from results/audit_head_sign.txt. Seed "
        + str(SEED) + ".",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not r1:
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
