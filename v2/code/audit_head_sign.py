# -*- coding: utf-8 -*-
r"""
Why are the largest dilations all one sign, and does that go away?

WHAT IS AT STAKE

Item 4(b) of the standing note is down to one demand:
{#rem:gainopposition} closed the opposition route as a bounded factor,
so the only thing left that can raise e(G) from the 0.098386 that
within-block cancellation supplies to the ceiling theta'/2 = 0.28 is
the internal cancellation of the top blocks by magnitude. Those blocks
do not cancel: {#rem:gainprofile} measured their gains at 1.4639 and
2.1866 at the top N, and {#rem:gainsplit} measured that 0.8274 to
1.0000 of the top decile carries ONE sign.

Nothing has asked why. The answer is not far, and the paper already
contains the beginning of it. With

    H(N;k) = sum_{m<N/k, (m,k)=1} Lambda(N-mk) mu(m),

the discussion around [eq:layers] notes that when N/k is small the m
coprime to rad(N) are almost all primes, so mu(m) = -1 dominates and
H(N;k) "acquires a sign". That is a statement about the inner sum, and
it is measurable exactly by splitting the sum where the sign is:

    P_k = sum_{mu(m)=+1} Lambda(N-mk),   M_k = sum_{mu(m)=-1} ...,
    H_k = P_k - M_k,   T_k = P_k + M_k,   I_k = H_k / T_k in [-1,1].

Then |a_k| = (log k) T_k |I_k| factors the magnitude into a mass T_k
and an imbalance I_k, and the question of why the head is one-signed
becomes two questions that have different consequences. If the head is
the k with the largest T_k, its one-signedness is an accident of which
k have many prime values of N-mk, and it should wash out. If the head
is the k with the largest |I_k| -- the k where the inner sum fails to
cancel -- then the head is selected for exactly the failure that
matters, and whether it washes out depends on whether |I_k| shrinks
with N. That last is the number that decides whether the top fifth's
internal cancellation can improve at all by computing further.

BACKS: Remark {#rem:headsign} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The control. The head's one-sign fraction reproduces
      results/audit_gain_split.txt at all eight N to 0.001, and the
      whole-range gain reproduces its published value to 0.01, so
      H = P - M is the same field.
  Q2  The head is selected by imbalance and not by mass. The top
      block's average |I_k| exceeds the whole range's by more than half
      again at every N, while its average T_k does not.
  Q3  And the imbalance is signed. The share of the top block with
      I_k < 0 matches the head's one-sign fraction to within 0.05 at
      every N, so the head's alignment is the inner sum's sign and
      nothing else.
  Q4  The imbalance shrinks with N. The whole range's average |I_k| has
      a resolved negative exponent against log N. If it is flat, the
      head's one-signedness does not wash out and the top fifth's
      internal cancellation cannot be improved by computing further --
      which would leave item 4(b) with no computational route at all.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED at 0.001 on any head fraction or 0.01 on the gain.
      Either would mean the split is not of the field
      {#rem:gainsplit} measured and nothing below may be compared
      with it.
  Q2  REFUTED if the top block's average |I_k| fails to exceed the
      whole range's by half again at any N, or if its average T_k
      exceeds the whole range's by half again. Then the head is a
      mass selection and its one-signedness is not about cancellation.
  Q3  REFUTED if the two shares differ by more than 0.05 at any N.
      Then something other than the inner sum's sign is aligning the
      head, and that something has not been named.
  Q4  REFUTED if the average |I_k| has no resolved negative exponent.
      That is the outcome worth having and the worse one: it would say
      the obstruction is a standing property of the inner sum, so the
      top fifth cannot be made to cancel by taking N larger, and the
      one remaining route in item 4(b) is closed on the computational
      side.

  Q1 gates: without it this is not the same field.
  Q2, Q3 and Q4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. P_k and M_k are a deterministic
  partition of a sum by the sign of mu(m), and H = P - M is an
  identity; there is no background to detect against. The coin arms
  for the gain this feeds were run in audit_crossk_reference.py and,
  on the block decomposition, in lab_gain_opposition.py, whose
  block-sign arm is the control any statement about the head's
  alignment is read against.
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
OUT = os.path.join(RES, "audit_head_sign.txt")

BLOCKS = 10
HEAD = 0.10


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
    """the published head fractions and gains -- read, not typed"""
    src = io.open(os.path.join(RES, "audit_gain_split.txt"),
                  encoding="utf-8").read()
    g = {}
    for m in re.finditer(r"^  N = (\d+)\s+#k = \d+\s+head \d+\s+"
                         r"G ([\d.]+)\s+head [\d.]+\s+tail [\d.]+\s+"
                         r"mass [\d.]+\s*$", src, re.M):
        g[int(m.group(1))] = float(m.group(2))
    i = src.index("N            top-decile share  same sign in the head")
    ag = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) != 3 or not f[0].isdigit():
            break
        ag[int(f[0])] = float(f[2])
    return g, ag


def signed_parts(N, lam, mu, sqf):
    """the inner sum of [eq:layers] split by the sign of mu(m)

    Every line is audit_gain_split.weighted's, with the one product
    replaced by two masked sums, so P - M is that function's H term by
    term and Q1 checks it.
    """
    PN = SPL.factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    Ps, Ms, Cs = [], [], []
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in SPL.factor_set(k):
            ms = ms[ms % q != 0]
        lm = lam[N - ms * k]
        g = mu[ms]
        Ps.append(float(lm[g > 0].sum()))
        Ms.append(float(lm[g < 0].sum()))
        Cs.append(int(np.count_nonzero(lm[g != 0])))
    return ks, np.array(Ps), np.array(Ms), np.array(Cs)


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

    pubg, pubag = read_split()
    say("read from results/audit_gain_split.txt: %d whole-range gains "
        "and the" % len(pubg))
    say("  head's one-sign fractions, %.4f to %.4f"
        % (min(pubag.values()), max(pubag.values())))
    say("  the sieve and theta' are imported from "
        "code/audit_gain_split.py;")
    say("  the inner sum is the same one, split by the sign of mu(m).")
    say("  the ceiling e(l1/l2) can reach is theta'/2 = %.4f"
        % (THETA / 2.0))

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    G, agree, hI, wI, hT, wT, hneg = {}, {}, {}, {}, {}, {}, {}
    wR, hR, wC = {}, {}, {}
    byI, byT, byK, topIneg = {}, {}, {}, {}
    mI, negshare, kmI = {}, {}, {}
    bI, bT, bneg, bagree = {}, {}, {}, {}
    say()
    say("  N            #k      G        head agree  |I| whole  "
        "|I| head  T whole    T head")
    for N in NS:
        ks, P, M, C = signed_parts(N, lam, mu, sqf)
        H = P - M
        T = P + M
        a = np.log(ks.astype(np.float64)) * H
        n = ks.size
        ok = T > 0
        I = np.zeros(n)
        I[ok] = H[ok] / T[ok]

        w = np.abs(a)
        order = np.argsort(-w)
        s = abs(float(a.sum()))
        G[N] = float(w.sum()) / s if s > 0 else float("inf")

        nh = max(1, int(round(HEAD * n)))
        hd = order[:nh]
        sh = np.sign(a[hd])
        agree[N] = max(float((sh > 0).mean()), float((sh < 0).mean()))
        hI[N] = float(np.abs(I[hd]).mean())
        wI[N] = float(np.abs(I[ok]).mean())
        hT[N] = float(T[hd].mean())
        wT[N] = float(T[ok].mean())
        hneg[N] = float((I[hd] < 0).mean())
        good = C > 0
        R = np.zeros(n)
        R[good] = 1.0 / np.sqrt(C[good].astype(np.float64))
        wR[N] = float(R[good & ok].mean())
        hR[N] = float(R[hd].mean())
        wC[N] = float(C[good].mean())

        # the same deciles cut on three axes on their own, so the
        # sign correlation can be attributed to one of them
        edges = [int(round(d * n / BLOCKS)) for d in range(BLOCKS + 1)]
        for nm, key, store in (("I", -np.abs(I), byI),
                               ("T", -T, byT),
                               ("k", -ks.astype(np.float64), byK)):
            o = np.argsort(key)
            store[N] = [float((I[o[edges[d]:edges[d + 1]]] < 0).mean())
                        for d in range(BLOCKS)]
        topIneg[N] = byI[N][0]
        mI[N] = float(I[ok].mean())          # the SIGNED mean
        negshare[N] = float((I[ok] < 0).mean())
        ok2 = np.argsort(ks)
        kmI[N] = [float(I[ok2[edges[d]:edges[d + 1]]].mean())
                  for d in range(BLOCKS)]
        bI[N], bT[N], bneg[N], bagree[N] = [], [], [], []
        for d in range(BLOCKS):
            idx = order[edges[d]:edges[d + 1]]
            bI[N].append(float(np.abs(I[idx]).mean()))
            bT[N].append(float(T[idx].mean()))
            bneg[N].append(float((I[idx] < 0).mean()))
            sg = np.sign(a[idx])
            bagree[N].append(max(float((sg > 0).mean()),
                                 float((sg < 0).mean())))
        say("  %-12d %-7d %-8.4f %-11.4f %-10.6f %-10.6f %-9.2f %.2f"
            % (N, n, G[N], agree[N], wI[N], hI[N], wT[N], hT[N]))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- Q1
    say()
    say("Q1  the control: the same field, split")
    say("  N            agree here  published   diff       G here    "
        "published")
    wa, wg = 0.0, 0.0
    for N in NS:
        da = abs(agree[N] - pubag[N])
        dg = abs(G[N] - pubg[N])
        wa, wg = max(wa, da), max(wg, dg)
        say("  %-12d %-11.4f %-11.4f %-10.6f %-9.4f %.4f"
            % (N, agree[N], pubag[N], da, G[N], pubg[N]))
    q1 = wa < 0.001 and wg < 0.01
    say("  worst departure in the fraction %.6f, in the gain %.6f"
        % (wa, wg))
    say("  Q1 %s   (cap 0.001 on a fraction, cap 0.01 on a gain)"
        % ("hold" if q1 else "REFUTED"))

    # ------------------------------------------------------------- Q2
    say()
    say("Q2  is the head selected by imbalance or by mass?")
    say("  N            |I| head/whole   T head/whole")
    q2 = True
    for N in NS:
        ri, rt = hI[N] / wI[N], hT[N] / wT[N]
        if not (ri > 1.5 and rt <= 1.5):
            q2 = False
        say("  %-12d %-16.4f %.4f" % (N, ri, rt))
    say("  the imbalance ratio has to exceed 1.5 (cap) and the mass "
        "ratio not to")
    say("  Q2 %s" % ("hold" if q2 else "REFUTED"))

    # ------------------------------------------------------------- Q3
    say()
    say("Q3  is the head's alignment the inner sum's sign?")
    say("  N            share I<0 in head  one-sign fraction  diff")
    q3 = True
    for N in NS:
        d = abs(hneg[N] - agree[N])
        if d > 0.05:
            q3 = False
        say("  %-12d %-18.4f %-18.4f %.4f" % (N, hneg[N], agree[N], d))
    say("  Q3 %s   (cap 0.05)" % ("hold" if q3 else "REFUTED"))
    if q3:
        say("  and the difference is exactly zero, which is a warning "
            "about")
        say("  the check and not a strength of it. T_k > 0, so")
        say("  sign(I_k) = sign(H_k) = sign(a_k) identically and the "
            "two")
        say("  columns are the same column whenever the head's majority "
            "is")
        say("  negative. Q3 therefore measures one thing only -- that "
            "the")
        say("  majority IS negative at every N, which is")
        say("  {#rem:headidentity}'s finding -- and the exact agreement "
            "is")
        say("  the identity showing, not a second confirmation.")

    # ------------------------------------------------------------- Q4
    say()
    say("Q4  does the imbalance shrink with N?")
    ew, _rw, sew = fit(x, np.log(np.array([wI[N] for N in NS])))
    eh, _rh, seh = fit(x, np.log(np.array([hI[N] for N in NS])))
    et, _rt, set_ = fit(x, np.log(np.array([wT[N] for N in NS])))
    q4 = ew < 0.0 and abs(ew) / sew >= 2.0
    say("  quantity          exponent     s.e.       t")
    say("  |I| whole range   %+-12.6f %-10.6f %.2f" % (ew, sew, abs(ew) / sew))
    say("  |I| head          %+-12.6f %-10.6f %.2f" % (eh, seh, abs(eh) / seh))
    say("  T whole range     %+-12.6f %-10.6f %.2f"
        % (et, set_, abs(et) / set_))
    say("TSTAT slope_headsign_imbalance %.2f" % (abs(ew) / sew))
    say("SPREAD slope_headsign_imbalance %.4f" % float(x.max() - x.min()))
    if abs(ew) / sew < 2.0:
        say("UNRESOLVED SIGN slope_headsign_imbalance")
    say("TSTAT slope_headsign_headimbalance %.2f" % (abs(eh) / seh))
    say("SPREAD slope_headsign_headimbalance %.4f"
        % float(x.max() - x.min()))
    if abs(eh) / seh < 2.0:
        say("UNRESOLVED SIGN slope_headsign_headimbalance")
    say("  the whole range's average |I| runs %.6f to %.6f"
        % (min(wI.values()), max(wI.values())))
    say("  the head's runs %.6f to %.6f"
        % (min(hI.values()), max(hI.values())))
    say("  Q4 %s" % ("hold" if q4 else "REFUTED"))

    say()
    say("  and the rate to compare that with is not zero. The inner "
        "sum has")
    say("  n_k terms that actually contribute, the m with Lambda(N-mk) "
        "nonzero")
    say("  and mu(m) nonzero, and signs that cancelled at random would "
        "leave")
    say("  an imbalance of order 1/sqrt(n_k). That reference is "
        "computable")
    say("  term by term and has never been put beside |I_k|:")
    er, _rr, ser = fit(x, np.log(np.array([wR[N] for N in NS])))
    ec, _rc, sec = fit(x, np.log(np.array([wC[N] for N in NS])))
    say("  N            n_k avg    1/sqrt(n) avg  |I| avg    ratio")
    for N in NS:
        say("  %-12d %-10.2f %-14.6f %-10.6f %.4f"
            % (N, wC[N], wR[N], wI[N], wI[N] / wR[N]))
    rat = np.array([wI[N] / wR[N] for N in NS])
    era, _rea, sera = fit(x, np.log(rat))
    say("  quantity          exponent     s.e.       t")
    say("  n_k average       %+-12.6f %-10.6f %.2f"
        % (ec, sec, abs(ec) / sec))
    say("  1/sqrt(n) average %+-12.6f %-10.6f %.2f"
        % (er, ser, abs(er) / ser))
    say("  |I| over it       %+-12.6f %-10.6f %.2f"
        % (era, sera, abs(era) / sera))
    say("TSTAT slope_headsign_sqrtratio %.2f" % (abs(era) / sera))
    say("SPREAD slope_headsign_sqrtratio %.4f"
        % float(x.max() - x.min()))
    if abs(era) / sera < 2.0:
        say("UNRESOLVED SIGN slope_headsign_sqrtratio")
    say("  the ratio runs %.4f to %.4f" % (float(rat.min()),
                                           float(rat.max())))
    if abs(era) / sera < 2.0:
        say("  so the inner sum cancels at the random-sign rate to "
            "within an")
        say("  unresolved drift: the within-k cancellation is already "
            "as good")
        say("  as independent signs would make it, and there is nothing "
            "to")
        say("  win there. The whole of the deficit against theta'/2 is "
            "across")
        say("  k, not inside k.")
    else:
        say("  so the inner sum does NOT cancel at the random-sign "
            "rate, and")
        say("  the drift says which way: a rising ratio means the "
            "inner sums")
        say("  are getting worse than independent signs, a falling one "
            "better.")

    say()
    say("X1  which axis carries the sign correlation")
    say("  (written after Q3 turned out to be an identity and Q2 fell; "
        "not")
    say("  pre-registered). Under independent signs every decile would "
        "sit at")
    say("  half. The magnitude blocks do not, but |a_k| mixes three "
        "factors,")
    say("  so the same deciles are cut on each factor alone. Share "
        "with")
    say("  I_k < 0, at the top N:")
    say("  decile   by |a|       by |I|       by T         by k")
    for d in range(BLOCKS):
        say("  %-8d %-12.4f %-12.4f %-12.4f %.4f"
            % (d + 1, bneg[NMAX][d], byI[NMAX][d], byT[NMAX][d],
               byK[NMAX][d]))
    say("  spread from half, top decile minus bottom decile:")
    for nm, store in (("|a|", bneg), ("|I|", byI), ("T", byT),
                      ("k", byK)):
        say("  %-5s %+.4f" % (nm, store[NMAX][0] - store[NMAX][-1]))
    say()
    say("  and a correlation between size and sign is what a nonzero "
        "MEAN")
    say("  looks like from the tail, so the signed mean has to be on "
        "the")
    say("  table beside it. In units of the fluctuation scale the same")
    say("  1/sqrt(n_k) gives:")
    say("  N            mean I      share I<0   |mean I| / 1/sqrt(n)")
    for N in NS:
        say("  %-12d %+-11.6f %-11.4f %.4f"
            % (N, mI[N], negshare[N], abs(mI[N]) / wR[N]))
    em, _rem, sem = fit(x, np.log(np.abs(np.array([mI[N] for N in NS]))))
    eb, _reb, seb = fit(x, np.log(np.array([abs(mI[N]) / wR[N]
                                            for N in NS])))
    say("  quantity            exponent     s.e.       t")
    say("  |mean I|            %+-12.6f %-10.6f %.2f"
        % (em, sem, abs(em) / sem))
    say("  the same over 1/sqrt(n) %+-8.6f %-10.6f %.2f"
        % (eb, seb, abs(eb) / seb))
    say("TSTAT slope_headsign_bias %.2f" % (abs(eb) / seb))
    say("SPREAD slope_headsign_bias %.4f" % float(x.max() - x.min()))
    if abs(eb) / seb < 2.0:
        say("UNRESOLVED SIGN slope_headsign_bias")
    say("  and the signed mean by k-decile at the top N, to see whether "
        "it is")
    say("  a large-k effect as the [eq:layers] discussion suggests:")
    say("  " + "  ".join("%+.4f" % v for v in kmI[NMAX]))
    say("  spread across the k-deciles %.6f against the mean itself "
        "%+.6f" % (max(kmI[NMAX]) - min(kmI[NMAX]), mI[NMAX]))
    ei, _rei, sei = fit(x, np.log(np.array([topIneg[N] for N in NS])))
    say("  and the top decile by |I| runs %.4f to %.4f over the eight "
        "N," % (min(topIneg.values()), max(topIneg.values())))
    say("  exponent %+.6f at %.2f standard errors"
        % (ei, abs(ei) / sei))
    say("TSTAT slope_headsign_topIneg %.2f" % (abs(ei) / sei))
    say("SPREAD slope_headsign_topIneg %.4f"
        % float(x.max() - x.min()))
    if abs(ei) / sei < 2.0:
        say("UNRESOLVED SIGN slope_headsign_topIneg")

    # the profile, for the record
    say()
    say("  and the whole profile at the top N, by magnitude block:")
    say("  block    |I|         T           share I<0   one-sign")
    for d in range(BLOCKS):
        say("  %-8d %-11.6f %-11.2f %-11.4f %.4f"
            % (d + 1, bI[NMAX][d], bT[NMAX][d], bneg[NMAX][d],
               bagree[NMAX][d]))

    say()
    say("=" * 70)
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s"
        % tuple("hold" if v else "REFUTED" for v in (q1, q2, q3, q4)))

    head = [
        "STATISTIC: for every squarefree k < N^theta' coprime to N, the",
        "           inner sum of [eq:layers] split by the sign of mu(m):",
        "           P_k over mu(m) = +1 and M_k over mu(m) = -1, both",
        "           weighted by Lambda(N-mk), with H_k = P_k - M_k,",
        "           T_k = P_k + M_k and the imbalance I_k = H_k/T_k. The",
        "           average |I_k| and the average T_k over the whole",
        "           range, over the top decile by |a_k| = (log k)|H_k|",
        "           and over each of " + str(BLOCKS) + " equal-count",
        "           magnitude blocks; the share of each with I_k < 0;",
        "           the one-sign fraction of a_k in each; and the",
        "           least-squares exponents of the averages against",
        "           log N.",
        "NULL: none is run and none applies. P_k and M_k are a",
        "      deterministic partition of a sum by the sign of mu(m) and",
        "      H = P - M is an identity, so there is no background to",
        "      detect against. The coin arms for the gain this feeds",
        "      were run in audit_crossk_reference.py and, on the block",
        "      decomposition, in lab_gain_opposition.py, whose",
        "      block-sign arm is the control the head's alignment is",
        "      read against.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta';",
        "       m over 1 <= m < N/k with (m,k) = 1, so m with mu(m) = 0",
        "       contribute to neither part; Lambda and mu from an integer",
        "       sieve to " + str(NMAX) + "; the sieve, the field's k-set",
        "       and theta' are imported from code/audit_gain_split.py and",
        "       the published gains and head fractions are read from",
        "       results/audit_gain_split.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not q1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
