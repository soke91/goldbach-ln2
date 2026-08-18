# -*- coding: utf-8 -*-
r"""
The largest single loss in the chain: is it a level, or a decline?

WHAT IS AT STAKE

Remark {#rem:residuesigned} sizes what [eq:directcond] throws away
when it drops the signs across k in the residue: gains of +0.2932,
+0.2526, +0.2346, +0.2121 in theta', against H's 0.053. That is the
biggest number in the whole chain -- "of everything the programme
discards, the signs across k in the residue are worth more than the
split itself" -- and it rests on one script.

Two things about it have never been asked. First, those four numbers
FALL, monotonically, and nothing has tested whether that decline is
above its own noise. If the gain is shrinking, then the largest loss
in the chain is a finite-N quantity and the sentence above is about
the accessible range, not about the method. Second, the fifth N is
censored: the signed walk does not cross below the k = 100000 at which
beta is fitted, so it is printed as "none" and the trend rests on four
points.

The censoring is an artefact of one convention. beta is fitted on
k < 100000 and the walk is truncated at the same k, but nothing needs
those to be the same number: beta is a property of the split and the
walk is a sum. Keeping beta on the published range and continuing the
walk to k < 400000 leaves the statistic comparable and uncensors the
fifth N.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of lab_residue_signed.py's, so W1 and W2 are cross-checks
of the published crossings and not reruns.

BACKS: Remark {#rem:signedgain} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The control: the absolute crossings reproduce the published
      993, 1447, 2019, 3319, 5923 exactly.
  W2  The control: the four uncensored signed crossings reproduce the
      published 35597, 37623, 48957, 68669 exactly.
  W3  Continuing the walk to k < 400000 with beta unchanged locates
      the fifth crossing, so the trend has five points and not four.
  W4  The gain is a decline and not a level: the least-squares slope
      of the gain against log N is negative and reaches two standard
      errors.

REFUTATION RULE (fixed before the run)

  W1  REFUTED on any mismatch -- not the same statistic, and nothing
      below may be compared with {#rem:residuesigned}.
  W2  REFUTED on any mismatch, likewise.
  W3  REFUTED if the signed walk still does not cross below 400000.
      The fifth N stays censored and W4 is decided on four points,
      which by {#rem:slopes}'s standard is not enough range to
      resolve a slope.
  W4  REFUTED if the slope fails to reach two standard errors. Then
      the decline is not measurable and the published 0.21-0.29 must
      be read as a level: the largest loss in the chain would be a
      standing feature of the method rather than something the range
      is working off.

  All four gate.

  NO NULL IS RUN here and none applies to W1-W4: crossings of a
  deterministic walk against a computed threshold are located, and
  there is no background to detect against. The sign control for this
  field is lab_residue_signed.py's own sixteen-draw sign
  randomisation, which refuted its rule S3 and established that mu's
  signs are worse than random -- that control is not repeated, it is
  the reason the gain is worth sizing at all.
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
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_signed_gain.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KFIT = 100_000                      # beta is fitted here, as published
KWALK = 400_000                     # the walk is continued to here
QSIEVE = 30
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
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
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


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


def read_H_gain():
    """H's own gain in theta' -- read from the results file"""
    src = io.open(os.path.join(RES, "lab_signed_level.txt"),
                  encoding="utf-8").read()
    return float(re.search(r"i\.e\. about ([\d.]+) in the",
                           src).group(1))


def read_published():
    """the absolute and signed crossings, read from the results file"""
    src = io.open(os.path.join(RES, "lab_residue_signed.txt"),
                  encoding="utf-8").read()
    i = src.index("N            K*_R abs   K*_R signed  factor   "
                  "exponent  clears")
    ab, sg = {}, {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ab[int(f[0])] = int(f[1])
        if f[2].isdigit():
            sg[int(f[0])] = int(f[2])
    return ab, sg


def walk(N, lam, mu, sqf, vmask, qs, artin, twin):
    """H, P and R over the extended k-range, beta fitted on KFIT"""
    PN = factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, KWALK):
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
        for q in factor_set(k):
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
    H = np.array(Hs)
    P = np.array(Ps)
    fit = ks < KFIT
    beta = float((H[fit] * P[fit]).sum() / (P[fit] * P[fit]).sum())
    R = H - beta * P
    w = np.log(ks.astype(np.float64))
    thr = S_ * (1.0 - A_) * N
    absc = np.cumsum(w * np.abs(R))
    sgn = np.cumsum(w * R)
    ja = int(np.searchsorted(absc, thr))
    out = np.flatnonzero(np.abs(sgn) > thr)
    kabs = int(ks[ja]) if ja < ks.size else None
    ksig = int(ks[out[0]]) if out.size else None
    return kabs, ksig, beta, ks.size, int(fit.sum()), thr / N


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubab, pubsg = read_published()
    say("read %d published absolute and %d uncensored signed "
        "crossings from results/lab_residue_signed.txt"
        % (len(pubab), len(pubsg)))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    vmask = residue_mask(NMAX, qs)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    got = []
    for N in NS:
        kabs, ksig, beta, nk, nfit, bpn = walk(
            N, lam, mu, sqf, vmask, qs, artin, twin)
        got.append((N, kabs, ksig, beta, nk, nfit))
        say("  N = %-10d #k %-7d (%d below the fit cap)  beta %.6f"
            % (N, nk, nfit, beta))
        say("BUDGET kstar_R_signedgain_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- W1
    say()
    say("W1  the control: the absolute crossing")
    say("  N            here       published")
    w1 = True
    for N, kabs, ksig, beta, nk, nfit in got:
        if kabs != pubab.get(N):
            w1 = False
        say("  %-12d %-10s %s" % (N, kabs, pubab.get(N)))
    say("  W1 %s" % ("hold" if w1 else "REFUTED"))

    # ------------------------------------------------------------- W2
    say()
    say("W2  the control: the four uncensored signed crossings")
    say("  N            here       published")
    w2 = True
    for N, kabs, ksig, beta, nk, nfit in got:
        if N not in pubsg:
            continue
        if ksig != pubsg[N]:
            w2 = False
        say("  %-12d %-10s %d" % (N, ksig, pubsg[N]))
    say("  W2 %s" % ("hold" if w2 else "REFUTED"))

    # ------------------------------------------------------------- W3
    say()
    say("W3  the censored fifth N, with the walk continued to "
        "k < %d" % KWALK)
    cens = [g for g in got if g[0] not in pubsg]
    w3 = all(g[2] is not None for g in cens)
    for N, kabs, ksig, beta, nk, nfit in cens:
        say("  %-12d signed crossing %s   (published: none below %d)"
            % (N, ksig, KFIT))
    say("  W3 %s" % ("hold" if w3 else "REFUTED"))
    say("  the censoring is not at random: the k-cap hides exactly")
    say("  the LARGEST crossings, so a trend fitted on the survivors")
    say("  is biased towards whatever the small ones do. Counted and")
    say("  resolved:")
    say("CENSORED lab_residue_signed %d" % len(cens))
    say("UNCENSORED lab_residue_signed %d"
        % sum(1 for g in cens if g[2] is not None))

    # ------------------------------------------------------------- W4
    say()
    say("W4  the gain in theta', and whether it is falling")
    say("  N            log K*/log N abs   signed    gain")
    xs, gains = [], []
    for N, kabs, ksig, beta, nk, nfit in got:
        if kabs is None or ksig is None:
            continue
        ea = math.log(kabs) / math.log(N)
        es = math.log(ksig) / math.log(N)
        xs.append(math.log(N))
        gains.append(es - ea)
        say("  %-12d %-18.4f %-9.4f %+.4f" % (N, ea, es, es - ea))
    x = np.array(xs)
    y = np.array(gains)
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    rms = float(np.sqrt((r ** 2).mean()))
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    t = abs(float(a)) / se
    say("  least-squares slope against log N = %+.6f" % a)
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (rms, se, t))
    say("SCATTER slope_audit_signed_gain %.4f" % rms)
    say("TSTAT slope_audit_signed_gain %.2f" % t)
    say("SPREAD slope_audit_signed_gain %.4f" % float(x.max() - x.min()))
    if t < 2.0:
        say("UNRESOLVED SIGN slope_audit_signed_gain")
    w4 = (a < 0.0) and (t >= 2.0)
    say("  two-sigma interval [%+.6f, %+.6f]" % (a - 2 * se, a + 2 * se))
    say("  W4 the slope is negative and above two standard errors   %s"
        % ("hold" if w4 else "REFUTED"))

    say()
    hgain = read_H_gain()
    say("  What this rate is and is not. Per octave of N the gain")
    say("  falls by %.4f, against the %.4f in theta' that Remark"
        % (abs(a * math.log(2.0)), hgain))
    say("  {#rem:signedlevel} measured for H and this file reads")
    say("  from results/lab_signed_level.txt. No N at which")
    say("  the gain would vanish is quoted here: that is an")
    say("  extrapolation over a factor 16 in N and Remark")
    say("  {#rem:forecastbracket} is the reason not to make it. What")
    say("  is measured is that the largest single discard in the")
    say("  chain is not a constant of the method.")

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("the largest loss in the chain is a decline, not a level"
        if ok else "REFUTED")

    head = [
        "STATISTIC: at each N, the first k at which",
        "           sum_{k<K}(log k)|R(N;k)| exceeds S(N)(1-A(N))N and",
        "           the first at which the SIGNED walk",
        "           sum_{k<K}(log k)R(N;k) leaves the same interval;",
        "           their exponents log K*/log N; the difference, which",
        "           is what discarding the signs costs in theta'; and",
        "           the least-squares slope of that difference against",
        "           log N with its standard error.",
        "NULL: none is run and none applies to W1-W4. Crossings of a",
        "      deterministic walk against a computed threshold are",
        "      located; there is no background to detect against. The",
        "      sign control for this field is lab_residue_signed.py's",
        "      sixteen-draw randomisation, which refuted its rule S3.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < " + str(KWALK) + ", while beta is",
        "       fitted as sum(H P)/sum(P^2) on k < " + str(KFIT) + ",",
        "       the published range -- the two caps are separated here",
        "       because the fit is a property of the split and the walk",
        "       is a sum; m odd, squarefree and coprime to k, m < N/k;",
        "       the sieve weight over the odd primes below "
        + str(QSIEVE) + ";",
        "       the Euler products at the fixed bound " + str(CLIM)
        + ". Every",
        "       N is 2^a 5^b, one odd radical, as RADICALS declares.",
        "       The published crossings are read from",
        "       results/lab_residue_signed.txt.",
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
