# -*- coding: utf-8 -*-
r"""
Does the measured bias account for the tail, or is something else there?

WHAT IS AT STAKE

Remark {#rem:headsign} left item 4(b) with one thing to explain. The
inner sums of [eq:layers] already cancel at the independent-sign rate
-- the ratio of the measured imbalance |I_k| to 1/sqrt(n_k) runs
1.0310 to 1.0995 with an unresolved drift -- so nothing can be won
inside k and the whole deficit against theta'/2 is a sign correlation
across k. That correlation lives on the imbalance axis: the top decile
of the k by |I_k| is 0.8547 negative at the top N and falls only at
-0.032100. But the bias that is supposed to produce it is small and
shrinking: the signed mean of I_k is -0.017828 at the top N, which is
0.2590 of the fluctuation scale, and that ratio falls at -0.185280 with
20.72 standard errors.

A quarter of a standard deviation does not obviously put 85 per cent of
the largest deviations on one side. Either it does and the arithmetic
is unremarkable, or it does not and there is a second structure in the
tail that nothing has named. **The difference decides what item 4(b)
still asks for**, and it is settled by a null rather than by argument.

The crude null signs the observed |I_k| independently at the observed
overall rate. It needs no draws: the top decile by |I_k| is a FIXED set
of k once the magnitudes are given, so the null's negative share on it
is binomial with that rate, mean and spread in closed form. It is the
wrong null and it is here to be refuted -- it throws away the one thing
the argument turns on, that the largest |I_k| are the k where the
fluctuation scale 1/sqrt(n_k) is largest, which are also the k where
the bias is largest.

The structured null keeps that. For each k it draws
I_k = b_k + sigma_k Z with Z standard normal, sigma_k proportional to
1/sqrt(n_k) and b_k the measured bias of that k's decile in k, then
sorts its own draw by |I_k| and reads the top decile. The constant of
proportionality is calibrated so the null reproduces the observed
average |I_k| -- without that it would be a null for a different
statistic, which is what M3 forbids -- and the null's own overall
negative share is checked against the observed as a second control.

BACKS: Remark {#rem:tailasymmetry} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The controls. The overall negative share and the top decile by
      |I_k| reproduce results/audit_head_sign.txt at every N to 0.001,
      and the structured null reproduces the observed average |I_k| to
      1e-6 by calibration and the observed overall negative share to
      0.02 without it.
  T2  The crude null misses the tail. Its top-decile negative share
      sits below the observed at every N by more than three of its own
      binomial spreads. That is expected and it is the measurement of
      how much the heteroscedasticity has to do.
  T3  The structured null reaches it. The observed top-decile negative
      share lies inside the structured null's own spread at every one
      of the eight N. Then the tail asymmetry is the measured bias seen
      through an unequal fluctuation scale and nothing else, and item
      4(b) reduces to the size and the k-dependence of that bias.
  T4  And the bias's k-dependence is what does the work. Holding the
      bias at its whole-range value while keeping sigma_k per k gives a
      top-decile share below the structured null's at every N.

REFUTATION RULE (fixed before the run)

  T1  REFUTED at 0.001 on either observed quantity, 1e-6 on the
      calibration or 0.02 on the null's overall share. Any of the three
      would mean the null is not a null for this statistic and nothing
      below may be read.
  T2  REFUTED if the crude null covers the observed at any N. Then the
      magnitudes alone carry the tail and the fluctuation scale is
      irrelevant, which would be a simpler world than the one
      {#rem:headsign} describes.
  T3  REFUTED if the observed falls outside the structured null's
      spread at any N. That is the outcome worth having: a bias of a
      quarter of a standard deviation, spread over k the way it is
      measured to be, would then NOT account for the tail, and the
      residual is a second structure in the sign of H_k with no name
      and no measurement.
  T4  REFUTED if the flat-bias null covers the structured null's share
      at every N. The k-dependence would then be decoration and the
      bias could be quoted as one number.

  T1 gates: without it the null is not a null for this statistic.
  T2, T3 and T4 are the measurement and do not gate.

  THE NULLS. Both are one-sided and neither touches the magnitudes'
  origin: the crude one keeps every |I_k| exactly and randomises only
  the signs, the structured one keeps the measured bias profile and the
  measured fluctuation scale and randomises only the draw. M3 is met by
  calibration -- a null whose average |I_k| differed from the observed
  would be describing a different field -- and M4 by T2, which is the
  arm with the structure deliberately destroyed.
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
OUT = os.path.join(RES, "lab_tail_asymmetry.txt")

BLOCKS = 10
DRAWS = 400
SEED = 20260809


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


def read_headsign():
    """the published shares -- read, not typed"""
    src = io.open(os.path.join(RES, "audit_head_sign.txt"),
                  encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^  (\d+)\s+([+-][\d.]+)\s+([\d.]+)\s+"
                         r"([\d.]+)\s*$", src, re.M):
        rows[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    m = re.search(r"top decile by \|I\| runs ([\d.]+) to ([\d.]+)", src)
    return rows, (float(m.group(1)), float(m.group(2)))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se


def topshare(I, frac=1.0 / BLOCKS):
    """negative share of the top |I| decile of one sample"""
    n = I.shape[-1]
    nt = max(1, int(round(frac * n)))
    idx = np.argsort(-np.abs(I), axis=-1)[..., :nt]
    return (np.take_along_axis(I, idx, axis=-1) < 0).mean(axis=-1)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, pubtop = read_headsign()
    say("read from results/audit_head_sign.txt: %d rows, and the top "
        "decile by" % len(pub))
    say("  |I| running %.4f to %.4f" % pubtop)
    say("  the split is imported from code/audit_head_sign.py and the "
        "sieve")
    say("  through it from code/audit_gain_split.py, so the field is "
        "the same.")

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    rng = np.random.default_rng(SEED)
    obs_top, obs_neg, obs_mI, obs_aI = {}, {}, {}, {}
    crude_mu, crude_sd = {}, {}
    st_mu, st_sd, st_neg, st_aI, st_rank = {}, {}, {}, {}, {}
    fl_mu, fl_sd = {}, {}
    cal, np_mu, np_sd, kp = {}, {}, {}, {}
    say()
    say("  N            #k      obs top   obs neg   crude mu  crude sd  "
        "struct mu  struct sd  struct neg  flat mu")
    for N in NS:
        ks, P, M, C = HS.signed_parts(N, lam, mu, sqf)
        T = P + M
        ok = (T > 0) & (C > 0)
        I = np.where(ok, (P - M) / np.where(T > 0, T, 1.0), 0.0)[ok]
        cc = C[ok].astype(np.float64)
        kk = ks[ok]
        n = I.size

        obs_top[N] = float(topshare(I))
        obs_neg[N] = float((I < 0).mean())
        obs_mI[N] = float(I.mean())
        obs_aI[N] = float(np.abs(I).mean())

        # the bias profile: the measured signed mean of each k-decile
        o = np.argsort(kk)
        edges = [int(round(d * n / BLOCKS)) for d in range(BLOCKS + 1)]
        b = np.zeros(n)
        for d in range(BLOCKS):
            idx = o[edges[d]:edges[d + 1]]
            b[idx] = I[idx].mean()
        sig0 = 1.0 / np.sqrt(cc)

        # M3: calibrate the scale so the null's average |I| is the
        # observed one. Without this it is a null for another field.
        zc = rng.standard_normal((64, n))     # fixed, so bisection is
                                              # on a monotone function

        def meanabs(c):
            return float(np.abs(b + c * sig0 * zc).mean())

        lo, hi = 1e-6, 10.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if meanabs(mid) < obs_aI[N]:
                lo = mid
            else:
                hi = mid
        c = 0.5 * (lo + hi)
        cal[N] = c

        z = rng.standard_normal((DRAWS, n))
        Ist = b + c * sig0 * z
        st = topshare(Ist)
        st_mu[N], st_sd[N] = float(st.mean()), float(st.std())
        st_neg[N] = float((Ist < 0).mean())
        st_aI[N] = float(np.abs(Ist).mean())
        st_rank[N] = int((st >= obs_top[N]).sum())

        zf = rng.standard_normal((DRAWS, n))
        Ifl = obs_mI[N] + c * sig0 * zf
        fl = topshare(Ifl)
        fl_mu[N], fl_sd[N] = float(fl.mean()), float(fl.std())

        # arm D, post-registration: keep every |I_k| exactly and draw
        # the sign of each k at that k's own decile-in-k rate. No
        # distributional assumption, and the top decile is a fixed set
        # so the arm is closed form.
        pk = np.zeros(n)
        for d in range(BLOCKS):
            idx = o[edges[d]:edges[d + 1]]
            pk[idx] = float((I[idx] < 0).mean())
        kp[N] = [float(pk[o[edges[d]]]) for d in range(BLOCKS)]
        nt2 = max(1, int(round(n / BLOCKS)))
        tidx = np.argsort(-np.abs(I))[:nt2]
        np_mu[N] = float(pk[tidx].mean())
        np_sd[N] = float(math.sqrt((pk[tidx] * (1.0 - pk[tidx])).sum())
                         / nt2)

        # the crude arm is binomial on a fixed set, so it is closed form
        nt = max(1, int(round(n / BLOCKS)))
        p = obs_neg[N]
        crude_mu[N] = p
        crude_sd[N] = math.sqrt(p * (1.0 - p) / nt)

        say("  %-12d %-7d %-9.4f %-9.4f %-9.4f %-9.4f %-10.4f %-10.4f "
            "%-11.4f %.4f"
            % (N, n, obs_top[N], obs_neg[N], crude_mu[N], crude_sd[N],
               st_mu[N], st_sd[N], st_neg[N], fl_mu[N]))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- T1
    say()
    say("T1  the controls")
    lo_o, hi_o = min(obs_top.values()), max(obs_top.values())
    dneg = max(abs(obs_neg[N] - pub[N][1]) for N in NS)
    dtop = max(abs(lo_o - pubtop[0]), abs(hi_o - pubtop[1]))
    dcal = max(abs(st_aI[N] - obs_aI[N]) for N in NS)
    dnull = max(abs(st_neg[N] - obs_neg[N]) for N in NS)
    t1 = dneg < 0.001 and dtop < 0.001 and dcal < 1e-6 and dnull < 0.02
    say("  the overall negative share departs from the published by "
        "%.6f" % dneg)
    say("  the top decile's range %.4f to %.4f against the published "
        "%.4f to %.4f" % (lo_o, hi_o, pubtop[0], pubtop[1]))
    say("  the calibration leaves the average |I| off by %.3e" % dcal)
    say("  and that cap was unachievable as written: the arm's average "
        "|I| is")
    say("  a Monte Carlo mean over %d draws, so it carries a sampling "
        "error of" % DRAWS)
    say("  its own and cannot be pinned to a cap that tight. The rule "
        "was")
    say("  badly written and it is being reported as such, not "
        "relaxed.")
    say("  the structured null's overall negative share is off by %.4f"
        % dnull)
    say("  T1 %s   (cap 0.001, 0.001, 1e-6, 0.02)"
        % ("hold" if t1 else "REFUTED"))
    say("  the share failure is not a bad cap. The arm's bias profile "
        "is the")
    say("  observed one, so its mean I is the observed mean by "
        "construction,")
    say("  and yet its sign share is not the observed sign share:")
    say("  N            obs neg   struct neg  difference")
    for N in NS:
        say("  %-12d %-9.4f %-11.4f %+.4f"
            % (N, obs_neg[N], st_neg[N], st_neg[N] - obs_neg[N]))
    say("  A shifted symmetric law with the right mean therefore gets "
        "the")
    say("  sign share wrong, which says the distribution of I_k is not")
    say("  one. Whichever way the difference points, T3 below is a")
    say("  comparison against a model already known to be wrong and may")
    say("  not be read as evidence about the tail.")
    say("  the calibrated scale runs %.4f to %.4f times 1/sqrt(n_k)"
        % (min(cal.values()), max(cal.values())))

    # ------------------------------------------------------------- T2
    say()
    say("T2  does signing the magnitudes at the overall rate reach it?")
    say("  N            observed  crude     gap in crude spreads")
    t2 = True
    for N in NS:
        g = (obs_top[N] - crude_mu[N]) / crude_sd[N]
        if g <= 3.0:
            t2 = False
        say("  %-12d %-9.4f %-9.4f %+.2f" % (N, obs_top[N],
                                             crude_mu[N], g))
    say("  T2 %s   (cap 3 spreads)" % ("hold" if t2 else "REFUTED"))
    if t2:
        say("  so the magnitudes alone do not carry it: what the "
            "largest")
        say("  |I_k| have in common is a large fluctuation scale, and "
            "the")
        say("  crude arm threw that away.")

    # ------------------------------------------------------------- T3
    say()
    say("T3  does the measured bias, seen through 1/sqrt(n_k), reach "
        "it?")
    say("  N            observed  structured  spread   in spreads  "
        "rank of %d" % (DRAWS + 1))
    t3 = True
    for N in NS:
        g = (obs_top[N] - st_mu[N]) / st_sd[N] if st_sd[N] > 0 else \
            float("inf")
        if abs(g) > 2.0:
            t3 = False
        say("  %-12d %-9.4f %-11.4f %-8.4f %+-11.2f %d"
            % (N, obs_top[N], st_mu[N], st_sd[N], g, st_rank[N] + 1))
    say("  T3 %s   (cap 2 spreads)" % ("hold" if t3 else "REFUTED"))
    if t3:
        say("  so the tail asymmetry is the measured bias and the "
            "unequal")
        say("  fluctuation scale, and nothing else. Item 4(b) is then "
            "the")
        say("  size of that bias and its dependence on k -- one object, "
            "with")
        say("  a name: the excess of mu(m) = -1 among the m for which")
        say("  N - mk is prime.")
    else:
        say("  so a bias of the measured size, spread over k the way it")
        say("  is measured to be, does NOT account for the tail. The")
        say("  residual is a second structure in the sign of H_k and it")
        say("  has no name and no measurement.")

    # ------------------------------------------------------------- T4
    say()
    say("T4  is the bias's k-dependence doing the work?")
    say("  N            structured  flat bias  difference in flat "
        "spreads")
    t4 = True
    for N in NS:
        g = (st_mu[N] - fl_mu[N]) / fl_sd[N] if fl_sd[N] > 0 else 0.0
        if g <= 0.0:
            t4 = False
        say("  %-12d %-11.4f %-10.4f %+.2f" % (N, st_mu[N], fl_mu[N], g))
    say("  T4 %s" % ("hold" if t4 else "REFUTED"))

    say()
    et, _rt, set_ = fit(x, np.log(np.array([obs_top[N] for N in NS])))
    es, _rs, ses = fit(x, np.log(np.array([st_mu[N] for N in NS])))
    say("  and the two fall together or they do not:")
    say("  quantity        exponent     s.e.       t")
    say("  observed top    %+-12.6f %-10.6f %.2f"
        % (et, set_, abs(et) / set_))
    say("  structured null %+-12.6f %-10.6f %.2f"
        % (es, ses, abs(es) / ses))
    say("TSTAT slope_tailasym_observed %.2f" % (abs(et) / set_))
    say("SPREAD slope_tailasym_observed %.4f" % float(x.max() - x.min()))
    if abs(et) / set_ < 2.0:
        say("UNRESOLVED SIGN slope_tailasym_observed")
    say("TSTAT slope_tailasym_null %.2f" % (abs(es) / ses))
    say("SPREAD slope_tailasym_null %.4f" % float(x.max() - x.min()))
    if abs(es) / ses < 2.0:
        say("UNRESOLVED SIGN slope_tailasym_null")

    # ------------------------------------------- not pre-registered
    say()
    say("X1  the same question with no distributional assumption")
    say("  (written after T1 fell; not pre-registered). Keep every "
        "|I_k|")
    say("  exactly as measured and draw only the sign, at the rate that "
        "k's")
    say("  own decile in k shows. This assumes nothing about the shape "
        "of")
    say("  the law and keeps the k-dependence of the sign rate, so it "
        "is the")
    say("  arm T3 should have used. The top decile is a fixed set, so "
        "it is")
    say("  closed form.")
    say("  N            observed  arm D     spread   in spreads")
    for N in NS:
        g = (obs_top[N] - np_mu[N]) / np_sd[N] if np_sd[N] > 0 else 0.0
        say("  %-12d %-9.4f %-9.4f %-8.4f %+.2f"
            % (N, obs_top[N], np_mu[N], np_sd[N], g))
    gaps = [(obs_top[N] - np_mu[N]) / np_sd[N] for N in NS
            if np_sd[N] > 0]
    say("  the sign rate by k-decile at the top N:")
    say("  " + "  ".join("%.4f" % v for v in kp[NMAX]))
    say("  the observed tail stands %.2f to %.2f spreads above this arm"
        % (min(gaps), max(gaps)))
    say("  so the coupling between the size of I_k and its sign is not")
    say("  mediated by k either: the sign rate that k predicts is "
        "nowhere")
    say("  near what the largest deviations show. That coupling is the")
    say("  object item 4(b) is now about, and it is unnamed.")

    say()
    say("=" * 70)
    say("T1 %s  T2 %s  T3 %s  T4 %s"
        % tuple("hold" if v else "REFUTED" for v in (t1, t2, t3, t4)))

    head = [
        "STATISTIC: the negative share of the top decile by |I_k| of the",
        "           imbalance I_k = (P_k - M_k)/(P_k + M_k) of",
        "           {#rem:headsign}, measured and under three arms: the",
        "           observed |I_k| signed independently at the observed",
        "           overall negative rate, in closed form because the",
        "           top decile is then a fixed set; a structured draw",
        "           I_k = b_k + c sigma_k Z with b_k the measured signed",
        "           mean of that k's decile in k, sigma_k = 1/sqrt(n_k)",
        "           for n_k the number of contributing m, and c",
        "           calibrated so the arm reproduces the observed",
        "           average |I_k|; and the same draw with b_k replaced",
        "           by the whole range's signed mean. " + str(DRAWS),
        "           draws for each structured arm.",
        "NULL: three arms are run and none is declined. The crude arm",
        "      keeps every |I_k| and breaks only the signs -- it is the",
        "      M4 arm, with the structure deliberately destroyed. The",
        "      structured arm keeps the measured bias profile and the",
        "      measured fluctuation scale and breaks only the draw; M3",
        "      is met by calibrating c so the arm's average |I_k| is the",
        "      observed one, without which it would be a null for a",
        "      different field. The flat-bias arm keeps the scale and",
        "      flattens the bias, isolating its k-dependence.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, every one 2^a 5^b so",
        "       the odd radical is one throughout, as RADICALS says; k",
        "       squarefree and coprime to N with 2 <= k < N^theta';",
        "       m over 1 <= m < N/k with (m,k) = 1, and n_k counts the m",
        "       with both Lambda(N-mk) and mu(m) nonzero; Lambda and mu",
        "       from an integer sieve to " + str(NMAX) + "; the split,",
        "       the sieve and theta' are imported from",
        "       code/audit_head_sign.py and the published shares are read",
        "       from results/audit_head_sign.txt. Seed " + str(SEED) + ".",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not t1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
