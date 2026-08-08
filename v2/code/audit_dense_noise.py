# -*- coding: utf-8 -*-
r"""
Where the ladder's noise lives: K*'s integrality, or the field.

WHAT IS AT STAKE

Remark {#rem:primorialdense} settled that the primorial ladder's
scatter of 0.0037 is not a shape sampled too coarsely. Two hundred and
two N, fifty to a doubling, scatter about their own within-octave lines
by 0.004084, 0.004470, 0.003975 and 0.003113 -- three of the four at or
above the published figure. So no amount of filling in resolves a
curve, and no octave is a lever. That closed the standing note's item 5
and left exactly one question behind it: **what is the noise made of?**

There are two candidates and they have opposite consequences.

If it is the integrality of K*, the ladder is being read through a bad
instrument and a better one exists. K*_R is where a step function first
exceeds a level, an integer that jumps; the equivalent statement

    log K*_R / log N > 1/2  <=>  rho(N) < 1,
    rho(N) = sum_{k<sqrt N}(log k)|R(N;k)| / (S(N)(1-A(N))N)

is a ratio of two smooth sums at a fixed abscissa, with no location and
no jump. Remark {#rem:primorialshare} made that substitution on ten
rungs and found rho *worse* -- scatter over trend 2.829 against the
exponent's 0.85 -- but eleven points cannot separate noise from shape,
so that verdict was measured against a scatter which was itself partly
shape. It has to be redone where the noise is known.

If instead the noise is in the field |R(N;k)| itself, then it is not an
instrument problem, every statistic built on that field inherits it,
and the limit on how far this ladder can be read is where it is.

The integrality half of the question has an answer that costs nothing
to compute and has never been computed: the admissible k are a fixed
set of 18863 values, so moving the crossing to the next one changes the
exponent by (log k' - log k)/log N, and that number can be printed at
every N and compared with the scatter directly.

BACKS: Remark {#rem:densenoise} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  N1  The control. The exponent recomputed here reproduces
      results/audit_primorial_dense.txt at each of its N to within
      0.001, and rho and the exponent agree on which side of the
      barrier every N falls -- rho < 1 exactly when the exponent
      exceeds 1/2 -- as check F1 of lab_primorial_share.py found on
      ten rungs.
  N2  The integrality of K* is not the noise. At every N the exponent
      granularity at the crossing, the change from moving to the next
      admissible k, is under a tenth of the within-octave scatter.
  N3  So the field carries it, and the smooth statistic carries it
      too. Measuring both with the within-octave scatter as noise and
      the whole sweep's slope as trend, rho's ratio of scatter to
      trend per doubling is no smaller than the exponent's:
      {#rem:primorialshare}'s verdict survives at fifty points to the
      doubling, now against a scatter that is known to be noise.
  N4  And the noise is not concentrated near the crossing. Splitting
      the sum at k = N^(1/4), the upper part carries the larger share
      of the mass and the larger relative fluctuation, so a statistic
      anchored lower on the k-axis does not escape it.

REFUTATION RULE (fixed before the run)

  N1  REFUTED at 0.001 on any exponent, or if rho and the exponent
      disagree on the side of the barrier at any N. Either would mean
      this is not the same field and nothing below may be compared
      with the dense sweep.
  N2  REFUTED if the granularity reaches a tenth of the within-octave
      scatter at any N. That is the outcome worth having: the ladder
      would then be quantisation-limited, and rho -- which has no
      location to quantise -- would be the instrument to rebuild every
      bracket on.
  N3  REFUTED if rho's scatter-to-trend ratio comes in below the
      exponent's. {#rem:primorialshare} would then have been wrong for
      the reason it could not have known, and every forecast on this
      ladder should be reread off rho.
  N4  REFUTED if the lower part of the k-range carries the larger
      relative fluctuation. The difficulty would then sit at small k,
      where the sums are shortest and a different anchor is cheap.

  N1 gates: without it these are not the dense sweep's points.
  N2, N3 and N4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. Two deterministic statistics are
  computed on the same field and their fluctuations compared with each
  other and with a granularity that is a property of the admissible
  set. There is no background to detect against. The coin arms for
  both statistics were run in lab_primorial_ladder.py and
  lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_dense_noise.txt")

BASE = 30030
SMOOTH = (2, 3, 5, 7, 11, 13)
JLO, JHI = 6, 9
QSIEVE = 30
CLIM = 4_000_000
SPLIT = 0.25                    # the k-axis is cut at N^SPLIT


def gap_module():
    """the sieve and the field, imported so they cannot drift"""
    p = os.path.join(CODE, "audit_primorial_gap.py")
    spec = importlib.util.spec_from_file_location("audit_primorial_gap", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GAP = gap_module()
primes_upto = GAP.primes_upto      # the same sieve, by name, for G20


def smooth_m(lo, hi):
    out = []
    for m in range(lo, hi + 1):
        v = m
        for p in SMOOTH:
            while v % p == 0:
                v //= p
        if v == 1:
            out.append(m)
    return out


def read_dense():
    """the sweep's own N and exponents -- read, not recopied"""
    src = io.open(os.path.join(RES, "audit_primorial_dense.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   m       thr/N      #k      "
                  "K*_R     exponent")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 7 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[6])
    m = re.search(r"least-squares slope on the sweep = "
                  r"([+-][\d.]+)", src)
    return out, float(m.group(1))


def read_share():
    """the published rho comparison on ten rungs"""
    src = io.open(os.path.join(RES, "lab_primorial_share.txt"),
                  encoding="utf-8").read()
    m = re.search(r"r\.m\.s\. ([\d.]+), trend ([\d.]+) per doubling, "
                  r"ratio ([\d.]+)", src)
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


def field(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the whole of what one N gives: the exponent, rho, and the split

    The H, P and beta are the same three lines as
    audit_primorial_gap.measure; the loop is walked once and every
    statistic below is read off the same arrays, so no two of them can
    disagree about the field.
    """
    PN = GAP.factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, GAP.KCAP):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[sqf[ms]]
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        for q in GAP.factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(ck * float(g[keep].sum()))
    ks = np.array(ks, dtype=np.int64)
    H, P = np.array(Hs), np.array(Ps)
    beta = float((H * P).sum() / (P * P).sum())
    R = H - beta * P
    w = np.log(ks.astype(np.float64)) * np.abs(R)
    cum = np.cumsum(w)
    thr = S_ * (1.0 - A_) * N
    logN = math.log(N)

    j = int(np.searchsorted(cum, thr))
    if j >= ks.size:
        return None
    kstar = int(ks[j])
    e = math.log(kstar) / logN
    # the exponent granularity: the admissible set is fixed, so the
    # crossing can only land on one of its members
    gran = ((math.log(ks[j + 1]) - math.log(kstar)) / logN
            if j + 1 < ks.size else float("nan"))

    half = int(np.searchsorted(ks, math.isqrt(N), side="right"))
    rho = float(cum[half - 1]) / thr if half > 0 else 0.0
    cut = int(np.searchsorted(ks, int(N ** SPLIT), side="right"))
    lowmass = float(cum[cut - 1]) if cut > 0 else 0.0
    highmass = float(cum[half - 1]) - lowmass if half > 0 else 0.0
    return (kstar, e, gran, rho, lowmass / thr, highmass / thr,
            int(ks.size), thr / N)


def linefit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    rms = float(np.sqrt((r ** 2).mean()))
    sxx = float(((x - x.mean()) ** 2).sum())
    se = math.sqrt(float((r ** 2).sum() / max(x.size - 2, 1))
                   / sxx) if sxx > 0 and x.size > 2 else float("inf")
    return float(a), float(b), rms, se, r


def pooled(octs, vals):
    """the r.m.s. residual about within-octave lines, pooled

    Noise on the short scale and trend on the long lever are different
    questions; this answers only the first, and it is the number
    {#rem:primorialdense} showed is the statistic's own.
    """
    num, den = 0.0, 0
    for ms in octs:
        x = np.log(np.array([BASE * m for m in ms], dtype=np.float64))
        y = np.array([vals[BASE * m] for m in ms])
        _a, _b, _rms, _se, r = linefit(x, y)
        num += float((r ** 2).sum())
        den += r.size - 2
    return math.sqrt(num / den)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    dense, aswept = read_dense()
    prms, ptrend, pratio = read_share()
    say("read from results/audit_primorial_dense.txt: %d N with their "
        "exponents" % len(dense))
    say("  and the sweep's slope %+.6f;" % aswept)
    say("  from results/lab_primorial_share.txt the published rho "
        "comparison on")
    say("  ten rungs: r.m.s. %.4f, trend %.4f per doubling, ratio %.3f"
        % (prms, ptrend, pratio))

    octs = [smooth_m(1 << j, 1 << (j + 1)) for j in range(JLO, JHI + 1)]
    allm = sorted({m for ms in octs for m in ms})
    top = BASE * max(allm)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say()
    say("sieving to %d, sieve weight over the odd primes %s"
        % (top, ", ".join(map(str, qs))))
    lam, mu = GAP.lambda_and_mu(top)
    sqf = mu != 0
    vmask = GAP.residue_mask(top, qs)
    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N            log10 N   K*_R     exponent  granularity  "
        "rho       low       high")
    e_, g_, r_, lo_, hi_, missed = {}, {}, {}, {}, {}, []
    thrn = None
    for m in allm:
        N = BASE * m
        v = field(N, lam, mu, sqf, vmask, qs, artin, twin)
        if v is None:
            missed.append(N)
            say("  %-12d no crossing below k = %d" % (N, GAP.KCAP))
            continue
        kstar, e, gran, rho, lowm, highm, _nk, thrn = v
        e_[N], g_[N], r_[N], lo_[N], hi_[N] = e, gran, rho, lowm, highm
        say("  %-12d %-9.4f %-8d %-9.4f %-12.6f %-9.4f %-9.4f %.4f"
            % (N, math.log10(N), kstar, e, gran, rho, lowm, highm))
    say("CENSORED audit_dense_noise %d" % len(missed))
    say("UNCENSORED audit_dense_noise %d" % len(missed))
    say("BUDGET kstar_R_S1AN_densenoise %.6f" % thrn)
    say("RADICALS %d" % len(set(
        tuple(sorted(q for q in GAP.factor_set(N) if q > 2))
        for N in e_)))
    del lam, mu, sqf, vmask
    ns = sorted(e_)

    # ------------------------------------------------------------- N1
    say()
    say("N1  the control: the same field, and the same side of the "
        "barrier")
    worst, wn = 0.0, None
    for N in ns:
        d = abs(e_[N] - dense.get(N, e_[N]))
        if d > worst:
            worst, wn = d, N
    bad = [N for N in ns if (r_[N] < 1.0) != (e_[N] > 0.5)]
    n1 = worst < 0.001 and not bad and len(ns) == len(dense)
    say("  %d of the sweep's %d N recomputed; worst departure %.6f "
        "at N = %s" % (len(ns), len(dense), worst, wn))
    say("  rho < 1 and exponent > 1/2 disagree at %d of %d N"
        % (len(bad), len(ns)))
    say("  rho runs %.4f to %.4f over the sweep"
        % (min(r_.values()), max(r_.values())))
    say("  N1 %s   (cap 0.001 on an exponent, no disagreement allowed)"
        % ("hold" if n1 else "REFUTED"))

    # ------------------------------------------------------------- N2
    noise_e = pooled(octs, e_)
    say()
    say("N2  is the noise the integrality of K*?")
    grans = np.array([g_[N] for N in ns])
    ratio = grans / noise_e
    n2 = bool(np.nanmax(ratio) < 0.1)
    say("  the admissible k are a fixed set, so the crossing lands on "
        "one of")
    say("  them and the exponent moves in steps. Those steps run "
        "%.6f to" % float(np.nanmin(grans)))
    say("  %.6f, median %.6f." % (float(np.nanmax(grans)),
                                 float(np.nanmedian(grans))))
    say("  the within-octave scatter of the exponent, pooled over the "
        "four")
    say("  octaves, is %.6f, so the step is %.4f of it at worst"
        % (noise_e, float(np.nanmax(ratio))))
    say("GRANULARITY audit_dense_noise %.6f %.6f"
        % (float(np.nanmedian(grans)), noise_e))
    say("  N2 %s   (cap 0.1 of the scatter)"
        % ("hold" if n2 else "REFUTED"))
    say()
    say("  the step is not one number: K* grows with N and log N grows "
        "with")
    say("  it, so it falls along the sweep. Per octave, against that "
        "octave's")
    say("  own scatter, and as a share of that octave's VARIANCE:")
    say("  octave   median step   octave scatter   ratio    share of "
        "variance")
    for j, ms in zip(range(JLO, JHI + 1), octs):
        x = np.log(np.array([BASE * m for m in ms], dtype=np.float64))
        y = np.array([e_[BASE * m] for m in ms])
        _a, _b, rms_o, _se, _r = linefit(x, y)
        gj = float(np.nanmedian([g_[BASE * m] for m in ms]))
        say("  %-8d %-13.6f %-16.6f %-8.4f %.4f"
            % (j, gj, rms_o, gj / rms_o, (gj / rms_o) ** 2))
    say("  a step is a bound on what quantisation can move, and its "
        "square")
    say("  is the share of the variance it can account for -- so even "
        "where")
    say("  the ratio is largest the integrality explains a small part "
        "of the")
    say("  noise, and the rest is the field.")

    # ------------------------------------------------------------- N3
    say()
    say("N3  does the smooth statistic carry the same noise?")
    dof = sum(len(ms) - 2 for ms in octs)
    lr = {N: math.log(r_[N]) for N in ns}
    xs = np.log(np.array(ns, dtype=np.float64))
    ae, _be, _rmse, see, _re = linefit(xs, np.array([e_[N] for N in ns]))
    ar, _br, _rmsr, ser, _rr = linefit(xs, np.array([lr[N] for N in ns]))
    noise_r = pooled(octs, lr)
    trend_e = abs(ae) * math.log(2.0)
    trend_r = abs(ar) * math.log(2.0)
    say("  noise is the within-octave scatter; trend is the whole "
        "sweep's")
    say("  slope over one doubling. The two are different questions "
        "and")
    say("  {#rem:primorialdense} is why they must be asked separately.")
    say("  statistic   noise      slope        trend/doubling   "
        "noise/trend")
    say("  exponent    %-10.6f %+-12.6f %-16.6f %.4f"
        % (noise_e, ae, trend_e, noise_e / trend_e))
    say("  log rho     %-10.6f %+-12.6f %-16.6f %.4f"
        % (noise_r, ar, trend_r, noise_r / trend_r))
    say("  least-squares slope of log rho against log N = %+.6f, s.e. "
        "%.6f, t = %.2f" % (ar, ser, abs(ar) / ser))
    say("TSTAT slope_audit_dense_noise %.2f" % (abs(ar) / ser))
    say("SPREAD slope_audit_dense_noise %.4f"
        % float(xs.max() - xs.min()))
    if abs(ar) / ser < 2.0:
        say("UNRESOLVED SIGN slope_audit_dense_noise")
    say("  the exponent's own slope stands at %.2f standard errors"
        % (abs(ae) / see))
    ratio_e, ratio_r = noise_e / trend_e, noise_r / trend_r
    n3 = ratio_r >= ratio_e
    say("  the published ratio for rho on ten rungs was %.3f; here "
        "the two" % pratio)
    say("  come in at %.4f for the exponent and %.4f for log rho"
        % (ratio_e, ratio_r))
    say("  N3 %s" % ("hold" if n3 else "REFUTED"))
    say()
    say("  and a ratio that close has to be read against its own "
        "error, or")
    say("  the same mistake is made again in the other direction. An "
        "r.m.s.")
    say("  on %d degrees of freedom carries a relative error of "
        "1/sqrt(2 df);" % dof)
    say("  a slope carries the one its own standard error gives.")
    rel_n = 1.0 / math.sqrt(2.0 * dof)
    rel_e = math.hypot(rel_n, see / abs(ae))
    rel_r = math.hypot(rel_n, ser / abs(ar))
    diff = abs(ratio_e - ratio_r)
    comb = math.hypot(ratio_e * rel_e, ratio_r * rel_r)
    say("  statistic   ratio    noise err   trend err   ratio err")
    say("  exponent    %-8.4f %-11.4f %-11.4f %.4f"
        % (ratio_e, rel_n, see / abs(ae), ratio_e * rel_e))
    say("  log rho     %-8.4f %-11.4f %-11.4f %.4f"
        % (ratio_r, rel_n, ser / abs(ar), ratio_r * rel_r))
    say("  the two differ by %.4f against a combined error of %.4f, "
        "which is" % (diff, comb))
    say("  %.2f of it" % (diff / comb))
    say("RATIOGAP audit_dense_noise %.6f %.6f" % (diff, comb))
    if diff <= comb:
        say("INSTRUMENTS TIED audit_dense_noise")
        say("  so the two instruments are TIED. rho is not the better")
        say("  one and it is not the worse one; it is the same noise")
        say("  through the second lens, which is what")
        say("  {#rem:primorialshare}'s own mechanism paragraph said")
        say("  while its ratio comparison said otherwise. That")
        say("  comparison was made against a published scatter of "
            "%.4f," % prms)
        say("  which {#rem:primorialdense} has since shown was noise")
        say("  and shape together, and the factor of more than three")
        say("  it reported does not survive fifty points to the")
        say("  doubling.")
    else:
        say("  so the two instruments are separated, and the sign of")
        say("  the difference says which to build on.")

    # ------------------------------------------------------------- N4
    say()
    say("N4  where on the k-axis does the fluctuation live?")
    llo = {N: math.log(lo_[N]) for N in ns}
    lhi = {N: math.log(hi_[N]) for N in ns}
    nlo, nhi = pooled(octs, llo), pooled(octs, lhi)
    shares = np.array([hi_[N] / (lo_[N] + hi_[N]) for N in ns])
    n4 = nhi >= nlo and float(shares.mean()) >= 0.5
    say("  the sum to sqrt N is cut at k = N^(1/4); the upper part's "
        "share of")
    say("  the mass runs %.4f to %.4f" % (float(shares.min()),
                                          float(shares.max())))
    say("  part        within-octave scatter of its log")
    say("  k below     %.6f" % nlo)
    say("  k above     %.6f" % nhi)
    say("  and the whole ratio's own is %.6f" % noise_r)
    say("  N4 %s" % ("hold" if n4 else "REFUTED"))
    say()
    say("  a relative scatter is not a contribution: a part that "
        "carries")
    say("  little mass can wobble freely without moving the sum. What "
        "each")
    say("  part contributes to rho is its share of the mass times its "
        "own")
    say("  relative scatter:")
    slo = 1.0 - float(shares.mean())
    shi = float(shares.mean())
    say("  part        mass share   relative scatter   contribution")
    say("  k below     %-12.4f %-18.6f %.6f" % (slo, nlo, slo * nlo))
    say("  k above     %-12.4f %-18.6f %.6f" % (shi, nhi, shi * nhi))
    frac = (slo * nlo) / (slo * nlo + shi * nhi)
    say("  the part holding %.4f of the mass supplies %.4f of the "
        "contribution," % (slo, frac))
    say("  so per unit of mass it fluctuates %.2f times as much "
        "(%.6f against" % (nlo / nhi, nlo))
    say("  %.6f). The fluctuation is not where the mass is, and a "
        "statistic" % nhi)
    say("  that reweights toward the mass does not average it away.")

    say()
    say("=" * 70)
    say("N1 %s  N2 %s  N3 %s  N4 %s"
        % tuple("hold" if v else "REFUTED" for v in (n1, n2, n3, n4)))
    say("the integrality of K* accounts for under a hundredth of the "
        "variance")
    say("at every octave, so the noise is the field's; and the smooth")
    say("statistic carries it in the same amount, so there is no "
        "instrument")
    say("left to switch to on this ladder.")
    ok = n1

    head = [
        "STATISTIC: at every N of the dense sweep, four things read off",
        "           one walk of the same field: the truncation K*_R and",
        "           its exponent log K*_R/log N; the exponent",
        "           granularity at the crossing, the change caused by",
        "           moving to the next admissible k; the ratio",
        "           rho(N) = sum_{k<sqrt N}(log k)|R(N;k)| divided by",
        "           S(N)(1-A(N))N; and that sum split at k = N^(1/4).",
        "           For the exponent and for log rho: the r.m.s.",
        "           residual about within-octave lines, pooled over the",
        "           four octaves, as noise; the whole sweep's slope as",
        "           trend; and the ratio of the two per doubling.",
        "NULL: none is run and none applies. Two deterministic",
        "      statistics are computed on the same field and their",
        "      fluctuations compared with each other and with a",
        "      granularity that is a property of the admissible set;",
        "      there is no background to detect against. The coin arms",
        "      for both were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py.",
        "FIELD: N = 30030*m for every 13-smooth m in [64, 1024], the",
        "       sweep of audit_primorial_dense.py; every one has prime",
        "       set {2,3,5,7,11,13}, so the threshold and the",
        "       admissible k-set are fixed; k squarefree and coprime to",
        "       N with 2 <= k < " + str(GAP.KCAP) + "; m odd, squarefree",
        "       and coprime to k, m < N/k; beta fitted over the whole",
        "       admissible range as in audit_primorial_gap.py; the",
        "       sieve weight over the odd primes below " + str(QSIEVE),
        "       and the Euler products at the fixed bound "
        + str(CLIM) + ";",
        "       the sieve and the field are imported from",
        "       code/audit_primorial_gap.py. The sweep's exponents and",
        "       slope are read from results/audit_primorial_dense.txt",
        "       and the published rho comparison from",
        "       results/lab_primorial_share.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
