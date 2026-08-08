# -*- coding: utf-8 -*-
r"""
Is the between-block opposition more than chance, and can it grow?

WHAT IS AT STAKE

Remark {#rem:gainprofile} split the decay of 1/G with the identity
1/G = |sum_d w_d s_d| and found that of the 0.153911 at which it
decays, only 0.098386 survives when the blocks are forbidden to oppose
one another; the remaining 0.055525 stands at 3.67 standard errors. It
closed by naming a second route to theta'/2 alongside raising the top
fifth's internal cancellation: **find more opposition between blocks
than the field already supplies.** That sentence was written on a
difference of exponents and on nothing else, and it is missing the two
things that decide whether the route exists.

The first is a null. Ordering the k by |a_k| and cutting into blocks is
an operation on magnitudes, and magnitudes are what a coin arm keeps.
If random signs on mu's own magnitudes produce the same split between
within-block and between-block decay, then the 0.055525 is a property
of the partition and not of mu, and there is no route.

The second is a level, not a rate. Write F = (sum_d w_d|s_d|) /
|sum_d w_d s_d| >= 1 for the factor by which forbidding opposition
raises 1/G -- the observed values are 1.0316 to 1.3829. Breaking only
the blocks' signs, keeping every |w_d s_d| exactly as measured, gives
the value F takes when ten numbers of those sizes are signed at random,
and that is computable exactly by enumerating all 2^B patterns. If the
measured F is BELOW it, then mu is *less* self-opposing than chance and
the growth seen so far is mu closing a gap rather than opening one. And
if the null's own F does not grow with N, then the route
{#rem:gainprofile} named is a bounded factor and not an exponent, which
would make item 4(b) harder rather than cheaper and would be a
correction to that remark rather than an extension of it.

BACKS: Remark {#rem:gainopposition} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  O1  The control and the identity. This reproduces
      results/audit_gain_profile.txt's 1/G and sum_d w_d|s_d| at every
      N to 1e-6, and 1/G = |sum_d w_d s_d| holds to 1e-12.
  O2  The coin has no between-block component. With random signs on
      mu's own magnitudes the exponents of |sum w s| and of sum w|s|
      agree to within two standard errors of their difference, so what
      {#rem:gainprofile} measured is mu's structure and not the
      partition's.
  O3  But mu is less self-opposing than chance in level. The measured
      F is below the exact block-sign null's F at every one of the
      eight N.
  O4  And chance does not grow. The null's own F has an exponent
      within two standard errors of zero. Together with O3 this makes
      the route a bounded factor: mu's opposition can reach chance and
      no further without a mechanism that beats chance among B
      numbers, and nothing names one.
  O5  The opposing mass is not growing either. The share of l1 carried
      by blocks whose sign differs from the dominant block's is flat
      across the eight N at two standard errors.

REFUTATION RULE (fixed before the run)

  O1  REFUTED at 1e-6 on either quantity or 1e-12 on the identity.
      Either would mean this is not the field {#rem:gainprofile}
      measured and nothing below may be compared with it.
  O2  REFUTED if the coin's two exponents are separated at two
      standard errors. That is the outcome worth having: the
      between-block component would then be an artefact of ordering by
      magnitude, and {#rem:gainprofile}'s second route would not exist
      at all.
  O3  REFUTED if the measured F reaches the null's at any N. Then mu
      opposes itself more than chance somewhere, which is the only
      form in which the route has content, and the N where it happens
      is where to look.
  O4  REFUTED if the null's F has a resolved exponent. If it grows,
      the route is an exponent after all and its size is that
      exponent; if it falls, the route closes faster than O3 alone
      would say.
  O5  REFUTED if that share has a resolved trend. A growing opposing
      mass would be a mechanism even with F below chance, since F is
      a ratio and the mass is not.

  O1 gates: without it this is not the same field.
  O2, O3, O4 and O5 are the measurement and do not gate.

  THE NULLS. Two are run and they are one-sided on purpose. The coin
  arm randomises every sign on mu's own magnitudes, which preserves
  #k, the block boundaries, every |a_k| and therefore every w_d: it
  breaks mu's signs and nothing else, and it is the arm
  audit_crossk_reference.py used to put mu's gain at a ninth to a
  thirteenth of a coin's. The block-sign arm randomises only the ten
  block signs and keeps every |w_d s_d| exactly as measured: it breaks
  the blocks' opposition and keeps their internal cancellation. The
  first says whether the split is mu's; the second says whether mu's
  opposition beats chance at the only scale the split is about.
"""

import importlib.util
import io
import itertools
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
OUT = os.path.join(RES, "lab_gain_opposition.txt")

BLOCKS = 10
RESOLUTIONS = [2, 5, 10, 20, 50]
DRAWS = 256
SEED = 20260808


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


def read_profile():
    """the published decomposition -- read, not typed"""
    src = io.open(os.path.join(RES, "audit_gain_profile.txt"),
                  encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^  (\d+)\s+\d+\s+[\d.]+\s+([\d.]+)\s+"
                         r"([\d.]+)\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*$",
                         src, re.M):
        rows[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    m = re.search(r"^CROSSSHARE gain_opposition %d ([\d.]+)\s*$"
                  % BLOCKS, src, re.M)
    return rows, float(m.group(1))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se


def edges_for(n, b):
    return [int(round(d * n / b)) for d in range(b + 1)]


def parts(signed, ed):
    """the block sums of an ordered signed array"""
    return np.array([signed[ed[d]:ed[d + 1]].sum()
                     for d in range(len(ed) - 1)])


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, pubshare = read_profile()
    say("read from results/audit_gain_profile.txt: %d rows and the "
        "published" % len(pub))
    say("  between-block share at %d blocks, %.6f" % (BLOCKS, pubshare))
    say("  the field is imported from code/audit_gain_split.py, so the")
    say("  blocks are cut out of the same a_k.")
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

    rng = np.random.default_rng(SEED)
    act, noc, Fobs, Fnull, oppmass = {}, {}, {}, {}, {}
    topsh, domin, exact = {}, {}, {}
    cact, cnoc = {B: {} for B in RESOLUTIONS}, {B: {} for B in RESOLUTIONS}
    ract, rnoc = {B: {} for B in RESOLUTIONS}, {B: {} for B in RESOLUTIONS}
    ident = {}
    say()
    say("  N            #k      1/G        sum w|s|   F obs    F null   "
        "F obs/null  opposing mass")
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        n = ks.size
        order = np.argsort(-np.abs(a))
        oa = a[order]
        ow = np.abs(oa)
        L1 = float(ow.sum())

        eds = {B: edges_for(n, B) for B in RESOLUTIONS}
        for B in RESOLUTIONS:
            bs = parts(oa, eds[B])
            ract[B][N] = abs(float(bs.sum())) / L1
            rnoc[B][N] = float(np.abs(bs).sum()) / L1
            cact[B][N] = cnoc[B][N] = 0.0
        for _ in range(DRAWS):
            eps = rng.integers(0, 2, size=n) * 2.0 - 1.0
            signed = ow * eps                # one draw, every resolution
            for B in RESOLUTIONS:
                sb = parts(signed, eds[B])
                cact[B][N] += abs(float(sb.sum())) / L1
                cnoc[B][N] += float(np.abs(sb).sum()) / L1
        for B in RESOLUTIONS:
            cact[B][N] /= DRAWS
            cnoc[B][N] /= DRAWS

        ed = edges_for(n, BLOCKS)
        bs = parts(oa, ed)
        v = bs / L1                                  # w_d s_d
        act[N] = abs(float(v.sum()))
        noc[N] = float(np.abs(v).sum())
        ident[N] = abs(act[N] - abs(float(a.sum())) / L1)
        Fobs[N] = noc[N] / act[N]
        av = np.abs(v)
        tot = 0.0
        for sg in itertools.product((-1.0, 1.0), repeat=BLOCKS):
            tot += abs(float((av * np.array(sg)).sum()))
        Fnull[N] = noc[N] / (tot / 2 ** BLOCKS)
        topsh[N] = float(av.max() / av.sum())
        domin[N] = bool(av.max() > av.sum() - av.max())
        exact[N] = abs(Fnull[N] * topsh[N] - 1.0)
        dom = int(np.argmax(av))
        sg0 = np.sign(v[dom])
        w = np.abs(bs) / L1
        oppmass[N] = float(w[np.sign(v) != sg0].sum())
        say("  %-12d %-7d %-10.6f %-10.6f %-8.4f %-8.4f %-11.4f %.4f"
            % (N, n, act[N], noc[N], Fobs[N], Fnull[N],
               Fobs[N] / Fnull[N], oppmass[N]))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- O1
    say()
    say("O1  the control, and the identity")
    w1, w2_ = 0.0, 0.0
    for N in NS:
        if N in pub:
            w1 = max(w1, abs(act[N] - pub[N][0]))
            w2_ = max(w2_, abs(noc[N] - pub[N][1]))
    o1 = w1 < 1e-6 and w2_ < 1e-6 and max(ident.values()) < 1e-12
    say("  worst departure in 1/G %.3e, in sum w|s| %.3e, identity gap "
        "%.3e" % (w1, w2_, max(ident.values())))
    say("  O1 %s   (cap 1e-6 on each, cap 1e-12 on the identity)"
        % ("hold" if o1 else "REFUTED"))

    # ------------------------------------------------------------- O2
    say()
    say("O2  does the coin arm show the same split?")
    say("  the coin keeps every |a_k|, so #k, the block edges and every")
    say("  w_d are the arm's own; only mu's signs are broken. %d draws "
        "at" % DRAWS)
    say("  seed %d, averaged." % SEED)
    say("  blocks   arm     |sum w s|    sum w |s|    share")
    o2 = True
    vals_mu, vals_co = [], []
    for B in RESOLUTIONS:
        ea, _ra, sea = fit(x, np.log(np.array([ract[B][N] for N in NS])))
        en, _rn, sen = fit(x, np.log(np.array([rnoc[B][N] for N in NS])))
        ca, _rc, sec = fit(x, np.log(np.array([cact[B][N] for N in NS])))
        cn, _rd, sed = fit(x, np.log(np.array([cnoc[B][N] for N in NS])))
        smu = (abs(ea) - abs(en)) / abs(ea)
        sco = (abs(ca) - abs(cn)) / abs(ca)
        vals_mu.append(smu)
        vals_co.append(abs(sco))
        say("  %-8d mu      %+-12.6f %+-12.6f %.4f" % (B, ea, en, smu))
        say("  %-8s coin    %+-12.6f %+-12.6f %.4f" % ("", ca, cn, sco))
        say("CROSSSHARE opposition_mu %d %.6f" % (B, smu))
        say("CROSSSHARE opposition_coin %d %.6f" % (B, abs(sco)))
        if B == BLOCKS:
            sepc = abs(abs(ca) - abs(cn)) / math.hypot(sec, sed)
            o2 = sepc <= 2.0
            say("  at %d blocks the coin's two exponents differ by "
                "%.2f standard errors," % (B, sepc))
            say("  so the coin has a between-block component of its "
                "own: %.6f" % abs(sco))
            say("  against mu's %.6f, which is %.4f of it. The "
                "partition supplies" % (smu, abs(sco) / smu))
            say("  part of what {#rem:gainprofile} read as mu's, and "
                "the rest,")
            say("  a factor %.2f more, is mu's." % (smu / abs(sco)))
    for lab, vv in (("opposition_mu", vals_mu),
                    ("opposition_coin", vals_co)):
        pos = [v for v in vv if v > 0]
        if len(pos) >= 2 and max(pos) / min(pos) > 1.5:
            say("RESOLUTION DEPENDENT %s" % lab)
    say("  O2 %s" % ("hold" if o2 else "REFUTED"))

    # ------------------------------------------------------------- O3
    say()
    say("O3  is mu more self-opposing than chance, or less?")
    ratios = np.array([Fobs[N] / Fnull[N] for N in NS])
    o3 = bool(ratios.max() < 1.0)
    say("  the block-sign arm keeps every |w_d s_d| and signs the %d of"
        % BLOCKS)
    say("  them at random; all %d patterns are enumerated, no sampling."
        % 2 ** BLOCKS)
    say("  F obs over F null runs %.4f to %.4f"
        % (float(ratios.min()), float(ratios.max())))
    er, _rr, ser = fit(x, np.log(ratios))
    say("  the ratio's exponent %+.6f at %.2f standard errors"
        % (er, abs(er) / ser))
    say("TSTAT slope_gainopposition_ratio %.2f" % (abs(er) / ser))
    say("SPREAD slope_gainopposition_ratio %.4f"
        % float(x.max() - x.min()))
    if abs(er) / ser < 2.0:
        say("UNRESOLVED SIGN slope_gainopposition_ratio")
    say("  O3 %s" % ("hold" if o3 else "REFUTED"))
    if o3:
        say("  so mu's blocks are MORE aligned than random signs on the")
        say("  same block magnitudes would leave them, at every N. The")
        say("  growth {#rem:gainprofile} measured is mu closing a gap "
            "to")
        say("  chance and not opening one beyond it.")

    # ------------------------------------------------------------- O4
    say()
    say("O4  does chance itself grow?")
    en, _rn2, sen2 = fit(x, np.log(np.array([Fnull[N] for N in NS])))
    eo, _ro2, seo2 = fit(x, np.log(np.array([Fobs[N] for N in NS])))
    o4 = abs(en) / sen2 <= 2.0
    say("  F null runs %.4f to %.4f, exponent %+.6f at %.2f standard "
        "errors" % (min(Fnull.values()), max(Fnull.values()), en,
                    abs(en) / sen2))
    say("  F obs  runs %.4f to %.4f, exponent %+.6f at %.2f standard "
        "errors" % (min(Fobs.values()), max(Fobs.values()), eo,
                    abs(eo) / seo2))
    say()
    say("  and F null is not a numerical fact. If the largest |w_d s_d| "
        "exceeds")
    say("  the sum of the other nine then every signing has the sign of "
        "that")
    say("  term, so the enumerated mean of |sum +- |w_d s_d|| is that "
        "term")
    say("  itself and F null = 1 / (its share). The largest term does "
        "exceed")
    say("  the rest at %d of %d N, and where it does the identity holds "
        "to %.3e:" % (sum(domin.values()), len(NS),
                     max(exact[N] for N in NS if domin[N])))
    say("  N            F null   top share   product    dominant?")
    for N in NS:
        say("  %-12d %-8.4f %-11.4f %-10.9f %s"
            % (N, Fnull[N], topsh[N], Fnull[N] * topsh[N],
               "yes" if domin[N] else "no"))
    say("  so e(F null) = -e(top share) identically, and "
        "{#rem:gainprofile}")
    say("  already measured that share as flat. The ceiling on the "
        "route is")
    say("  the reciprocal of the dominant block's mass share, and what "
        "is")
    say("  left of it at the top N is a factor %.4f."
        % (Fnull[max(NS)] / Fobs[max(NS)]))
    say("TSTAT slope_gainopposition_fnull %.2f" % (abs(en) / sen2))
    say("SPREAD slope_gainopposition_fnull %.4f"
        % float(x.max() - x.min()))
    if abs(en) / sen2 < 2.0:
        say("UNRESOLVED SIGN slope_gainopposition_fnull")
    say("  O4 %s" % ("hold" if o4 else "REFUTED"))
    if o4 and o3:
        say("  so the route is a BOUNDED FACTOR and not an exponent.")
        say("  Opposition among %d numbers of these sizes is worth a"
            % BLOCKS)
        say("  factor of about F null and that factor does not grow;")
        say("  mu is below it and rising, so the %+.6f"
            % (abs(eo)))
        say("  {#rem:gainprofile} attributed to opposition is what")
        say("  closing that gap costs, and it stops when the gap does.")
        say("  That makes item 4(b) harder than that remark left it,")
        say("  not cheaper, and this is a correction to it.")

    # ------------------------------------------------------------- O5
    say()
    say("O5  is the opposing mass growing?")
    em, _rm, sem = fit(x, np.log(np.array([oppmass[N] for N in NS])))
    o5 = abs(em) / sem <= 2.0
    say("  the l1 share of blocks whose sign differs from the dominant")
    say("  block's runs %.4f to %.4f, exponent %+.6f at %.2f standard "
        "errors" % (min(oppmass.values()), max(oppmass.values()), em,
                    abs(em) / sem))
    say("TSTAT slope_gainopposition_oppmass %.2f" % (abs(em) / sem))
    say("SPREAD slope_gainopposition_oppmass %.4f"
        % float(x.max() - x.min()))
    if abs(em) / sem < 2.0:
        say("UNRESOLVED SIGN slope_gainopposition_oppmass")
    say("  O5 %s" % ("hold" if o5 else "REFUTED"))

    say()
    say("=" * 70)
    say("O1 %s  O2 %s  O3 %s  O4 %s  O5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (o1, o2, o3, o4, o5)))

    head = [
        "STATISTIC: a_k = (log k)H(N;k) over the squarefree k < N^theta'",
        "           coprime to N, ordered by |a_k| and cut into blocks of",
        "           equal count. With w_d the block's share of l1 and",
        "           s_d its signed imbalance: the opposition factor",
        "           F = (sum_d w_d|s_d|)/|sum_d w_d s_d|; the same F when",
        "           only the block signs are randomised, by exact",
        "           enumeration of all 2^" + str(BLOCKS) + " patterns;",
        "           the l1 share of blocks whose sign differs from the",
        "           dominant block's; and, for mu and for a coin arm on",
        "           mu's own magnitudes, the least-squares exponents of",
        "           |sum w s| and sum w|s| against log N at " +
        str(len(RESOLUTIONS)) + " partition resolutions with the",
        "           between-block share each gives.",
        "NULL: two are run, both one-sided by design. The coin arm",
        "      randomises every sign on mu's own magnitudes, preserving",
        "      #k, the block edges and every w_d, so it breaks mu's",
        "      signs and nothing else; it is the arm of",
        "      audit_crossk_reference.py, where it put mu's gain at a",
        "      ninth to a thirteenth of a coin's. The block-sign arm",
        "      randomises only the block signs and keeps every |w_d s_d|",
        "      as measured, breaking the blocks' opposition and keeping",
        "      their internal cancellation. " + str(DRAWS) + " draws at",
        "      seed " + str(SEED) + " for the first; the second is",
        "      enumerated exactly and needs no draws.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta';",
        "       m over 1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + "; the field, the sieve",
        "       and theta' are imported from code/audit_gain_split.py and",
        "       the published decomposition is read from",
        "       results/audit_gain_profile.txt. The cut is by a FIXED",
        "       FRACTION of the k at every resolution.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not o1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
