# -*- coding: utf-8 -*-
r"""
The gain is a weighted average of ten signed imbalances. Which one
has to move?

WHAT IS AT STAKE

Remark {#rem:leanidentity} reduced item 4(b) of the standing note to
one exponent. With a_k = (log k)H(N;k), the three quantities that
matter are three norms of the same vector -- the gain
G = l1/|sum a|, the concentration l1/l2, and slope over floor
= |sum a|/(c l2) -- and the identity slope/floor = (l1/l2)/G/c is
exact. So the slope stops growing relative to its floor exactly when
e(G) catches e(l1/l2). Measured: e(G) = +0.153911 against
e(l1/l2) = +0.287798, and l1/l2 is capped by sqrt(#k) with
#k ~ N^theta', so the ceiling is theta'/2 = 0.28 and e(G) has to rise
by a factor of nearly two to meet it.

Remark {#rem:gainsplit} then split the range by magnitude and found
the top tenth of the k at +0.077963 and the remaining nine tenths at
+0.340006 -- **above** the ceiling. Remark {#rem:signmasshead} added
that the head carries only 44 to 52 per cent of what has to be
explained, and the standing note's item 4(b) closes on the sentence
that raising e(G) needs cancellation "across the whole range and not
in one decile."

That sentence is a conclusion drawn from a two-way split. A two-way
split cannot say where on the magnitude axis the cancellation fails,
and it cannot say whether the whole-range gain is bad because each
part is bad or because the parts fail to oppose one another. Those are
different demands with different costs, and there is an exact
decomposition that separates them, which has not been written down.

Order the k by |a_k| and cut into ten blocks of equal count. Write
L_d for block d's l1, S_d for its signed sum, w_d = L_d/l1 for its
share of the mass and s_d = S_d/L_d in [-1,1] for its signed
imbalance, so that 1/G_d = |s_d|. Then

    1/G = |sum_d w_d s_d|                                    (exact)

because both sides are |sum a| / l1. The gain is therefore the
reciprocal of a mass-weighted average of ten signed imbalances, and
two things can make it small:

  * every |s_d| small -- each block cancels internally;
  * the s_d of opposite sign -- the blocks cancel against each other.

Suppressing the second gives sum_d w_d |s_d|, an upper bound on
1/G that uses no cross-block cancellation at all. Comparing the two
says which mechanism the measured decay is made of, and the per-block
exponents say where on the axis it lives.

BACKS: Remark {#rem:gainprofile} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The control and the identity. The whole-range gain reproduces
      results/audit_gain_split.txt at all eight N to 0.01 and its
      exponent to 0.001, and 1/G = |sum_d w_d s_d| holds to 1e-12 at
      every N.
  W2  Nothing is gained across the blocks. The no-cross-cancellation
      value sum_d w_d |s_d| stays within a factor two of
      |sum_d w_d s_d| at every N, and the two exponents agree to
      within two standard errors of the difference. The demand is then
      on the blocks' internal cancellation and not on their
      opposition.
  W3  The profile is graded. The block gains G_d increase from the top
      block to the bottom at every N, and so do their exponents, with
      the top block's near zero and the bottom block's at or above
      theta'/2.
  W4  And one block decides. At every N the top block supplies the
      largest single term w_d|s_d|, and its share of the total rises
      with N.

REFUTATION RULE (fixed before the run)

  W1  REFUTED at 0.01 on any gain, 0.001 on the exponent, or 1e-12 on
      the identity. Any of the three would mean this is not the field
      {#rem:gainsplit} measured, or that the decomposition below is
      not an identity, and nothing after it may be read.
  W2  REFUTED if the two differ by more than a factor two at any N, or
      if their exponents are separated at two standard errors. That is
      the outcome worth having: the measured decay would then be
      partly cross-block, which is a cheaper thing to ask a proof for
      than ten blocks each cancelling to the square-root rate.
  W3  REFUTED if the block gains are not increasing at any N, or if
      more than one block has an exponent within its own standard
      error of zero. Either would mean the obstruction is not a
      gradient along the magnitude axis and the head is not the end
      of it.
  W4  REFUTED if any block other than the top supplies the largest
      term at any N, or if the top block's share of the total falls
      across the eight N. Then the whole-range gain is not the head's
      to fix even in the weighted sense, and {#rem:signmasshead}'s
      "only half" becomes an understatement.

  W1 gates: without it this is not the same field and the
  decomposition is not exact.
  W2, W3 and W4 are the measurement and do not gate.

  NO NULL IS RUN for the decomposition, which is a deterministic
  partition of a measured sequence by its own magnitudes and an
  identity between two ways of writing |sum a| / l1. The coin arms for
  G were run in audit_crossk_reference.py, where random signs on mu's
  own magnitudes gave 9.94 to 12.98 times mu's gain, and that is the
  background any statement about cancellation here is read against.
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
OUT = os.path.join(RES, "audit_gain_profile.txt")

BLOCKS = 10
RESOLUTIONS = [2, 5, 10, 20, 50]


def split_module():
    """the field of {#rem:gainsplit}, imported so it cannot drift"""
    p = os.path.join(CODE, "audit_gain_split.py")
    spec = importlib.util.spec_from_file_location("audit_gain_split", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = split_module()
NS = SPL.NS
THETA = SPL.THETA


def read_split():
    """the published whole-range gains and exponents -- read, not typed"""
    src = io.open(os.path.join(RES, "audit_gain_split.txt"),
                  encoding="utf-8").read()
    g = {}
    for m in re.finditer(r"^  N = (\d+)\s+#k = \d+\s+head \d+\s+"
                         r"G ([\d.]+)\s+head ([\d.]+)\s+tail ([\d.]+)\s+"
                         r"mass ([\d.]+)\s*$", src, re.M):
        g[int(m.group(1))] = (float(m.group(2)), float(m.group(3)),
                              float(m.group(4)), float(m.group(5)))
    e = {}
    for m in re.finditer(r"^  (whole|head tenth|tail)\s+"
                         r"([+-][\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$",
                         src, re.M):
        e[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return g, e


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubg, pube = read_split()
    say("read from results/audit_gain_split.txt: %d whole-range gains "
        "and the" % len(pubg))
    say("  exponents whole %+.6f, head %+.6f, tail %+.6f"
        % (pube["whole"][0], pube["head tenth"][0], pube["tail"][0]))
    say("  the field is imported from code/audit_gain_split.py, not "
        "copied,")
    say("  so the blocks below are cut out of the same a_k.")
    ceil_ = THETA / 2.0
    say("  the ceiling e(l1/l2) can reach is theta'/2 = %.4f" % ceil_)

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    G, ident, W, S, dom, domshare = {}, {}, {}, {}, {}, {}
    nocross, actual, nk = {}, {}, {}
    byres = {B: {} for B in RESOLUTIONS}
    say()
    say("  N            #k      G whole   1/G        sum w|s|   "
        "ratio    top term  its share")
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        w = np.abs(a)
        order = np.argsort(-w)
        n = ks.size
        nk[N] = n
        edges = [int(round(d * n / BLOCKS)) for d in range(BLOCKS + 1)]
        L1 = float(w.sum())
        for B in RESOLUTIONS:
            ed = [int(round(d * n / B)) for d in range(B + 1)]
            byres[B][N] = sum(
                abs(float(a[order[ed[d]:ed[d + 1]]].sum()))
                for d in range(B)) / L1
        ws, ss = [], []
        for d in range(BLOCKS):
            idx = order[edges[d]:edges[d + 1]]
            Ld = float(w[idx].sum())
            Sd = float(a[idx].sum())
            ws.append(Ld / L1)
            ss.append(Sd / Ld if Ld > 0 else 0.0)
        ws, ss = np.array(ws), np.array(ss)
        W[N], S[N] = ws, ss
        inv = abs(float((ws * ss).sum()))
        G[N] = 1.0 / inv
        ident[N] = abs(inv - abs(float(a.sum())) / L1)
        terms = ws * np.abs(ss)
        nocross[N], actual[N] = float(terms.sum()), inv
        dom[N] = int(np.argmax(terms))
        domshare[N] = float(terms.max() / terms.sum())
        say("  %-12d %-7d %-9.4f %-10.6f %-10.6f %-8.4f %-9.6f %.4f"
            % (N, n, G[N], inv, nocross[N], nocross[N] / inv,
               float(terms.max()), domshare[N]))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- W1
    say()
    say("W1  the control, and is the decomposition an identity?")
    say("  N            G here     published   diff       identity gap")
    worstg, worsti = 0.0, 0.0
    for N in NS:
        if N in pubg:
            d = abs(G[N] - pubg[N][0])
            worstg = max(worstg, d)
            say("  %-12d %-10.4f %-11.4f %-10.6f %.3e"
                % (N, G[N], pubg[N][0], d, ident[N]))
        worsti = max(worsti, ident[N])
    ew, rw, sew = fit(x, np.log(np.array([G[N] for N in NS])))
    de = abs(ew - pube["whole"][0])
    w1 = worstg < 0.01 and de < 0.001 and worsti < 1e-12
    say("  exponent %+.6f against the published %+.6f, diff %.6f"
        % (ew, pube["whole"][0], de))
    say("  worst gain departure %.6f, worst identity gap %.3e"
        % (worstg, worsti))
    say("  W1 %s   (cap 0.01 on a gain, cap 0.001 on the exponent, "
        "cap 1e-12 on the identity)" % ("hold" if w1 else "REFUTED"))
    say("TSTAT slope_audit_gain_profile %.2f" % (abs(ew) / sew))
    say("SPREAD slope_audit_gain_profile %.4f"
        % float(x.max() - x.min()))
    if abs(ew) / sew < 2.0:
        say("UNRESOLVED SIGN slope_audit_gain_profile")

    # ------------------------------------------------------------- W2
    say()
    say("W2  is any of the decay cross-block?")
    ratios = np.array([nocross[N] / actual[N] for N in NS])
    ea, ra, sea = fit(x, np.log(np.array([actual[N] for N in NS])))
    en, rn, sen = fit(x, np.log(np.array([nocross[N] for N in NS])))
    sed = math.hypot(sea, sen)
    w2 = bool(ratios.max() <= 2.0) and abs(ea - en) <= 2.0 * sed
    say("  suppressing the cross-block signs raises 1/G by a factor "
        "%.4f to %.4f" % (float(ratios.min()), float(ratios.max())))
    say("  quantity              exponent     s.e.")
    say("  |sum w s|             %+-12.6f %.6f" % (ea, sea))
    say("  sum w |s|             %+-12.6f %.6f" % (en, sen))
    say("  they differ by %+.6f against %.6f on the difference, "
        "which is %.2f" % (ea - en, sed, abs(ea - en) / sed))
    say("  W2 %s   (cap 2 on the factor, cap 2 on the standard errors)"
        % ("hold" if w2 else "REFUTED"))
    if w2:
        say("  so the blocks do not oppose one another to any resolved")
        say("  degree: the whole-range gain is the reciprocal of a")
        say("  mass-weighted average of ten internal imbalances, and")
        say("  every one of them has to move on its own.")
    else:
        say("  so part of the decay is between the blocks and not "
            "inside")
        say("  them. Of the %.6f the whole gain grows at, %.6f is what"
            % (abs(ea), abs(en)))
        say("  survives when the blocks are forbidden to oppose one")
        say("  another, and the remaining %.6f is their opposition."
            % (abs(ea) - abs(en)))

    # -------------------------------------------- not pre-registered
    say()
    say("X1  and how much of that is the partition's doing")
    say("  (this section was written after W2 fell and is not "
        "pre-registered)")
    say("  the split into within and between is NOT canonical. With B")
    say("  blocks, forbidding cancellation coarser than a block gives")
    say("  sum_d |sum_d a| / l1, which is 1/G at B = 1 and exactly 1 "
        "at")
    say("  B = #k. So the share attributed to opposition has to be "
        "read")
    say("  as a function of B or it is not a measurement at all.")
    say("  blocks   exponent     s.e.       cross share of %.6f"
        % abs(ea))
    shares = {}
    for B in RESOLUTIONS:
        y = np.log(np.array([byres[B][N] for N in NS]))
        eB, _rB, seB = fit(x, y)
        shares[B] = (abs(ea) - abs(eB)) / abs(ea)
        say("  %-8d %+-12.6f %-10.6f %.4f" % (B, eB, seB, shares[B]))
        say("CROSSSHARE gain_opposition %d %.6f" % (B, shares[B]))
    vals = [shares[B] for B in RESOLUTIONS]
    spread = max(vals) / min(vals) if min(vals) > 0 else float("inf")
    say("  the share runs %.4f to %.4f across the resolutions, a "
        "factor %.4f" % (min(vals), max(vals), spread))
    if spread > 1.5:
        say("RESOLUTION DEPENDENT gain_opposition")
        say("  so the number is not one number. What is invariant is")
        say("  the direction: at every resolution the exponent falls")
        say("  when opposition is forbidden, so opposition is a real")
        say("  part of the decay at every scale of the partition, and")
        say("  the coarsest reading -- B = 2, which is the split")
        say("  {#rem:gainsplit} published -- is the one that claims")
        say("  least.")

    # ------------------------------------------------------------- W3
    say()
    say("W3  where on the magnitude axis does cancellation fail?")
    GD = {N: 1.0 / np.maximum(np.abs(S[N]), 1e-300) for N in NS}
    say("  block    " + "".join("%-11d" % N for N in NS) + " w (top N)")
    for d in range(BLOCKS):
        say("  %-8d %s%.4f"
            % (d + 1, "".join("%-11.4f" % GD[N][d] for N in NS),
               W[NS[-1]][d]))
    say()
    say("  block    exponent     s.e.       t        against %.2f"
        % ceil_)
    bex, bse = [], []
    for d in range(BLOCKS):
        y = np.log(np.array([GD[N][d] for N in NS]))
        e_, _r, s_ = fit(x, y)
        bex.append(e_)
        bse.append(s_)
        say("  %-8d %+-12.6f %-10.6f %-8.2f %+.6f"
            % (d + 1, e_, s_, abs(e_) / s_, e_ - ceil_))
    bex, bse = np.array(bex), np.array(bse)
    mono = all(all(GD[N][d] <= GD[N][d + 1]
                   for d in range(BLOCKS - 1)) for N in NS)
    nzero = int((np.abs(bex) <= bse).sum())
    w3 = mono and nzero <= 1 and bex[-1] >= ceil_
    say("  the block gains increase from top to bottom at every N: %s"
        % ("yes" if mono else "no"))
    say("  blocks whose exponent is within its own standard error of "
        "zero: %d" % nzero)
    say("  the bottom block's exponent %+.6f against the ceiling %.4f"
        % (bex[-1], ceil_))
    say("  W3 %s" % ("hold" if w3 else "REFUTED"))

    # ------------------------------------------------------------- W4
    say()
    say("W4  which block decides the whole?")
    say("  N            top term   block   its share of sum w|s|")
    for N in NS:
        terms = W[N] * np.abs(S[N])
        say("  %-12d %-10.6f %-7d %.4f"
            % (N, float(terms.max()), dom[N] + 1, domshare[N]))
    alltop = all(dom[N] == 0 for N in NS)
    shs = np.array([domshare[N] for N in NS])
    esh, _rsh, sesh = fit(x, np.log(shs))
    w4 = alltop and esh > 0.0
    say("  the largest term is the top block at every N: %s"
        % ("yes" if alltop else "no"))
    say("  its share runs %.4f to %.4f, exponent %+.6f at %.2f "
        "standard errors"
        % (float(shs.min()), float(shs.max()), esh, abs(esh) / sesh))
    say("  W4 %s" % ("hold" if w4 else "REFUTED"))
    if abs(esh) / sesh < 2.0:
        say("  and it holds only by the letter of the rule. The share")
        say("  does not fall and the top block dominates at every N,")
        say("  but the RISE is not resolved, so the honest reading is")
        say("  that the share is flat: the top block's grip on the")
        say("  whole neither tightens nor loosens over a factor 128 "
            "in N.")
        say("TSTAT slope_gainprofile_topshare %.2f" % (abs(esh) / sesh))
        say("SPREAD slope_gainprofile_topshare %.4f"
            % float(x.max() - x.min()))
        say("UNRESOLVED SIGN slope_gainprofile_topshare")

    say()
    say("  what the profile costs. The mass-weighted shortfall against")
    say("  the ceiling is what a proof has to supply:")
    short = float((W[NS[-1]] * np.maximum(ceil_ - bex, 0.0)).sum())
    below = [d + 1 for d in range(BLOCKS) if bex[d] < ceil_]
    massbelow = float(W[NS[-1]][[d - 1 for d in below]].sum()) \
        if below else 0.0
    say("  blocks under the ceiling: %s" % (", ".join(map(str, below))
                                            if below else "none"))
    say("  they hold %.4f of the mass at the top N, and the "
        "mass-weighted" % massbelow)
    say("  shortfall in the exponent, summed over them, is %.6f"
        % short)

    say()
    say("=" * 70)
    say("W1 %s  W2 %s  W3 %s  W4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (w1, w2, w3, w4)))

    head = [
        "STATISTIC: a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, ordered by |a_k| and cut",
        "           into " + str(BLOCKS) + " blocks of equal count. For",
        "           each block d: its share w_d of l1, its signed",
        "           imbalance s_d = (sum a)/(sum|a|) and its gain",
        "           G_d = 1/|s_d|. The identity 1/G = |sum_d w_d s_d|",
        "           against |sum a|/l1; the no-cross-cancellation value",
        "           sum_d w_d |s_d| against it; each block's",
        "           least-squares exponent of log G_d on log N with its",
        "           standard error, against the ceiling theta'/2; and",
        "           the largest single term w_d|s_d| with its share.",
        "NULL: none is run for the decomposition, which is a",
        "      deterministic partition of a measured sequence by its own",
        "      magnitudes together with an identity between two ways of",
        "      writing |sum a|/l1. The coin arms for G were run in",
        "      audit_crossk_reference.py, where random signs on mu's own",
        "      magnitudes gave 9.94 to 12.98 times mu's gain.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta';",
        "       m over 1 <= m < N/k with (m,k) = 1; Lambda and mu from",
        "       an integer sieve to " + str(NMAX) + "; the field, the",
        "       sieve and theta' are imported from",
        "       code/audit_gain_split.py and the published gains and",
        "       exponents are read from results/audit_gain_split.txt.",
        "       The cut is by a FIXED FRACTION of the k, so every block",
        "       has #S of order N^theta' and the same square-root",
        "       reference theta'/2.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not w1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
