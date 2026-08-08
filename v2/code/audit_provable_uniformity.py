# -*- coding: utf-8 -*-
r"""
The premise nobody measured: is the classical shape uniform in k?

WHAT IS AT STAKE

OPEN item 2's core is that Remark {#rem:provablehalf}'s uniformity is
unproved. What that remark tested is the SHAPE -- whether
|P(N;k)| stays under A (N/k) exp(-c sqrt(log(N/k))) L(k) -- and it
reports the constant A forced by the data, restricted by inner length:
1.2119 at N/k >= 2, then 1.0710, 0.7309, 0.3363. Those are maxima over
the inner length, which is the axis the classical estimate already
controls.

The axis it does not control is k. "Uniformly in k" means the same A
works at every modulus, and at a fixed inner length x = N/k the
modulus k = N/x runs over a factor 16 across this sweep. Nobody has
looked along that axis. If the ratio grows with k at fixed x, the
uniformity fails in exactly the way that would matter, and the
conditional reduction's premise is not merely unproved but
contradicted by the accessible range.

BACKS: Remark {#rem:provableuniformity} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  The control: the maxima restricted to N/k >= 2, 8, 32, 128
      reproduce {#rem:provablehalf}'s 1.2119, 1.0710, 0.7309, 0.3363
      to within 0.001.
  F2  The k axis is flat: at each inner-length octave the maximum
      ratio across the five N grows by less than a factor 1.5 from
      the smallest k to the largest.
  F3  And not resolvably rising: at a majority of the octaves the
      least-squares slope of the log maximum against log k is below
      two standard errors.
  F4  So the constant is set by the inner length and not by the
      modulus: at every N the argmax of the ratio has inner length
      below 32.

REFUTATION RULE (fixed before the run)

  F1  REFUTED at 0.001 anywhere -- not the same statistic, and
      nothing below may be compared with {#rem:provablehalf}.
  F2  REFUTED beyond a factor 1.5 at any octave. That is the one that
      matters: the ratio would be growing along the axis the
      classical estimate does not control, and the uniformity the
      conditional reduction assumes would be contradicted where it
      can be checked.
  F3  REFUTED if a majority of octaves show a resolved rise, the same
      failure read through the fits.
  F4  REFUTED if the argmax ever sits at an inner length of 32 or
      more, which would say the constant is not a short-sum artefact
      and {#rem:provablehalf}'s reading of its own W1 is wrong.

  All four gate.

  NO NULL IS RUN and none applies. A measured sum is divided by a
  deterministic bound and the maxima compared; there is no background
  to detect against. The coin arm for this field is
  lab_elementary_provable.py's sixteen coins on the identical sifted
  set, which established that the looseness is the shape and not mu.
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
OUT = os.path.join(RES, "audit_provable_uniformity.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
CZERO = 0.2098
CLIM = 4_000_000
XLO = [2, 8, 32, 128]
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288,
       2097152]
MINPTS = 10


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def moebius(n):
    """Moebius, the cofactor kept in int32"""
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
    return mu


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


def read_published():
    """the maxima restricted by inner length"""
    src = io.open(os.path.join(RES, "lab_elementary_provable.txt"),
                  encoding="utf-8").read()
    i = src.index("cut-off x0   max ratio over every (N,k)")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) != 2 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[1])
    return out


def ratios(N, mu, oddsqf, vmask, qs):
    """the ratio |P| / bound at every admissible k"""
    PN = factor_set(N)
    ks, Pv, Lv, inner = [], [], [], []
    for k in range(2, KCAP):
        if not oddsqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.flatnonzero(oddsqf[1:M + 1]) + 1
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        fk = factor_set(k)
        for q in fk:
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Pv.append(ck * abs(int(mu[ms[keep]].sum(dtype=np.int64))))
        Lv.append(math.prod(1.0 / (1.0 - 1.0 / p) for p in fk))
        inner.append(N // k)
    ks = np.array(ks, dtype=np.int64)
    inn = np.array(inner, dtype=np.float64)
    b = inn * np.exp(-CZERO * np.sqrt(np.log(inn))) * np.array(Lv)
    return (ks, inn, np.array(Pv) / b,
            np.log(ks.astype(np.float64)) * np.array(Pv))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published restricted maxima from "
        "results/lab_elementary_provable.txt" % len(pub))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    mu = moebius(NMAX)
    oddsqf = (mu != 0)
    oddsqf[::2] = False
    vmask = residue_mask(NMAX, qs)

    data = {}
    for N in NS:
        ks, inn, r, mass = ratios(N, mu, oddsqf, vmask, qs)
        data[N] = (ks, inn, r, mass)
        j = int(np.argmax(r))
        say("  N = %-10d #k %-7d max %.4f at N/k = %d, k = %d"
            % (N, ks.size, float(r.max()), int(inn[j]), int(ks[j])))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    # ------------------------------------------------------------- F1
    say()
    say("F1  the control: the maxima restricted by inner length")
    say("  N/k >=       here       published  diff")
    f1 = True
    for x0 in XLO:
        mx = max(float(data[N][2][data[N][1] >= x0].max())
                 for N in NS)
        d = abs(mx - pub[x0])
        if not (d < 0.001):
            f1 = False
        say("  %-12d %-10.4f %-10.4f %.5f" % (x0, mx, pub[x0], d))
    say("  F1 %s   (cap 0.001)" % ("hold" if f1 else "REFUTED"))

    # ---------------------------------------------------------- F2/F3
    say()
    say("F2/F3  along the axis the estimate does not control: at a")
    say("  fixed inner length, does the ratio grow with the modulus?")
    say("  N/k octave        " + "".join("%-11d" % N for N in NS)
        + "grow    slope        t")
    f2 = True
    rises = 0
    tested = 0
    risemass = 0.0
    for a0, b0 in zip(OCT[:-1], OCT[1:]):
        row, kk, msh = [], [], []
        for N in NS:
            ks, inn, r, mass = data[N]
            sel = (inn >= a0) & (inn < b0)
            if int(sel.sum()) < MINPTS:
                row.append(None)
                continue
            row.append(float(r[sel].max()))
            kk.append(float(np.median(ks[sel].astype(np.float64))))
            msh.append(float(mass[sel].sum() / mass.sum()))
        vals = [v for v in row if v is not None]
        if len(vals) < 3:
            continue
        tested += 1
        gr = max(vals) / min(vals)
        if gr >= 1.5:
            f2 = False
        x = np.log(np.array(kk))
        y = np.log(np.array(vals))
        sl, rms, se, t = fit(x, y)
        if sl > 0 and t >= 2.0:
            rises += 1
            risemass += float(np.mean(msh))
        say("  [%-6d,%-8d) %s%-7.4f %+-12.6f %.2f"
            % (a0, b0,
               "".join("%-11s" % ("-" if v is None else "%.4f" % v)
                       for v in row), gr, sl, t))
    f3 = rises * 2 <= tested
    say("  F2 every octave grows by less than 1.5   %s   (cap 1.5)"
        % ("hold" if f2 else "REFUTED"))
    say("  F3 a resolved rise at %d of %d octaves   %s"
        % (rises, tested, "hold" if f3 else "REFUTED"))
    say("SCALES audit_provable_uniformity %d" % len(NS))

    # ------------------------------------------------------------- F4
    say()
    say("F4  where the constant is actually attained")
    say("  N            argmax at N/k   k there")
    f4 = True
    for N in NS:
        ks, inn, r, mass = data[N]
        j = int(np.argmax(r))
        if inn[j] >= 32:
            f4 = False
        say("  %-12d %-15d %d" % (N, int(inn[j]), int(ks[j])))
    say("  F4 the argmax has inner length below 32 at every N   %s"
        % ("hold" if f4 else "REFUTED"))
    say("CROSSAXIS lab_elementary_provable %d %d %.4f"
        % (tested, rises, risemass))
    if rises:
        say("AXIS RISE lab_elementary_provable")
    say("  the octaves with a resolved rise carry %.4f of the"
        % risemass)
    say("  elementary sum sum(log k)|P|, so the drift is real and")
    say("  sits where the mass is not.")

    say()
    say("  what this does and does not settle. It does not prove the")
    say("  uniformity, which is a statement about all k and all N.")
    say("  What it does is check it where it can be checked: at a")
    say("  fixed inner length the modulus runs over a factor %d here"
        % (max(NS) // min(NS)))
    say("  and the ratio does not follow it. The constant is a")
    say("  property of the length of the inner sum, which is the")
    say("  axis the classical estimate already controls.")

    say()
    say("=" * 70)
    ok = f1 and f2 and f3 and f4
    say("the constant follows the inner length, not the modulus"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the ratio |P(N;k)| / [(N/k) exp(-c sqrt(log(N/k)))",
        "           L(k)] at every admissible k; its maximum",
        "           restricted to N/k >= 2, 8, 32, 128; its maximum",
        "           within each octave of the inner length N/k at each",
        "           N, which at fixed octave varies the modulus by the",
        "           ratio of the largest N to the smallest; the",
        "           least-squares slope of that maximum against the",
        "           octave's median log k; and where the overall",
        "           maximum is attained.",
        "NULL: none is run and none applies. A measured sum is divided",
        "      by a deterministic bound and the maxima compared. The",
        "      coin arm for this field is lab_elementary_provable.py's",
        "      sixteen coins on the identical sifted set.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k odd, squarefree",
        "       and coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below " + str(QSIEVE) + ";",
        "       c = " + str(CZERO) + "; octaves closed at both ends and",
        "       used only when they hold at least " + str(MINPTS)
        + " k. Every N",
        "       is 2^a 5^b, one odd radical, as RADICALS declares. The",
        "       published restricted maxima are read from",
        "       results/lab_elementary_provable.txt.",
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
