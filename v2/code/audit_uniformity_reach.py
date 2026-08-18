# -*- coding: utf-8 -*-
r"""
The uniformity axis, with five more doublings of the modulus.

WHAT IS AT STAKE

OPEN item 2 rests on {#rem:provablehalf}'s uniformity, which is
unproved. {#rem:provableuniformity} checked it where it can be
checked -- at a fixed inner length N/k, vary the modulus k -- and
found the axis is not flat: two octaves of six show a resolved rise,
+0.167112 at t = 3.37 and +0.227768 at t = 2.18. Its own caution was
that the modulus only ran over a factor 16 there, and that those two
octaves carry 0.0673 of the elementary sum.

A factor of 16 is a short lever for a claim about all k, and the
computation turns out to be cheap: the same statistic at
N = 1.024e8 costs seconds, so the modulus can run over a factor 512
instead. Three things follow or fail. The rise may survive the longer
lever or wash out. The mass on the rising side may grow, which is what
would make the drift matter. And the overall maximum -- the number the
bound is really about -- may keep falling or stop.

One caveat is structural and is measured rather than argued: the
k-cap of 30000 is fixed, so at large N the long-inner-length octaves
are cut by it, and the share cut is reported for each.

BACKS: Remark {#rem:uniformityreach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  A1  The control. At the five published N the overall maximum ratio
      reproduces results/audit_provable_uniformity.txt inside the
      bound its printing forces.
  A2  The rise survives the lever: the two octaves that were resolved
      rising there are still resolved rising over the extended set.
  A3  And it starts to matter: the share of sum(log k)|P| carried by
      the resolved-rising octaves exceeds the published 0.0673.
  A4  But the bound gets safer, not tighter: the overall maximum
      ratio's least-squares slope in log N is resolved negative.

REFUTATION RULE (fixed before the run)

  A1  REFUTED outside the printing bound. Then this is not the
      statistic that remark measured. THIS ONE GATES.
  A2  REFUTED if either octave loses its resolved rise. Then the
      drift was a short-lever effect and the uniformity axis is flat
      as far as this can see -- which is what {#rem:provablehalf}
      needs.
  A3  REFUTED if the share does not grow. Then the drift stays where
      the mass is not, and it remains a curiosity rather than a
      threat to the reduction.
  A4  REFUTED if the maximum stops falling. That is the one that
      would matter: the classical bound's constant is the maximum,
      and a maximum that stops falling is a constant that stops
      improving. Note that "not resolved negative" includes "too
      noisy to tell", which is not the same as "stops falling" -- M9.

  A1 gates. A2 to A4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. A measured sum is divided by a
  deterministic bound and its maxima compared across N; there is no
  background to detect against. The coin arm for this field is
  lab_elementary_provable.py's sixteen coins on the identical sifted
  set.
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
OUT = os.path.join(RES, "audit_uniformity_reach.txt")

NS = [200_000 * 2 ** j for j in range(10)]
OCTS = [(8, 32), (32, 128), (128, 512), (512, 2048), (2048, 8192),
        (8192, 32768), (32768, 131072), (131072, 524288)]
MINK = 10


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PU = module("audit_provable_uniformity")
KCAP = PU.KCAP


def read_published():
    """the overall maxima and the two rising octaves' slopes"""
    src = io.open(os.path.join(RES, "audit_provable_uniformity.txt"),
                  encoding="utf-8").read()
    mx, dec = {}, 0
    for m in re.finditer(r"^  N = (\d+)\s+#k \d+\s+max ([\d.]+) at",
                         src, re.M):
        mx[int(m.group(1))] = float(m.group(2))
        dec = max(dec, len(m.group(2).split(".")[1]))
    slopes = {}
    for m in re.finditer(r"^  \[(\d+)\s*,(\d+)\s*\).*?([+-][\d.]+)\s+"
                         r"([\d.]+)\s*$", src, re.M):
        slopes[(int(m.group(1)), int(m.group(2)))] = (
            float(m.group(3)), float(m.group(4)))
    m = re.search(r"CROSSAXIS \S+ \d+ \d+ ([\d.]+)", src)
    return mx, dec, slopes, float(m.group(1))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if n > 2 else float("inf")
    return float(a), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubmx, dec, pubsl, pubshare = read_published()
    say("read %d published maxima and %d published octave slopes from"
        % (len(pubmx), len(pubsl)))
    say("  results/audit_provable_uniformity.txt, whose rising "
        "octaves carry %.4f of the elementary sum" % pubshare)
    say("  the statistic, the sieve, the k-cap %d and c are imported "
        "from code/audit_provable_uniformity.py" % KCAP)
    say()
    say("the doublings: N = %d to %d, %d of them, so the modulus at a"
        % (NS[0], NS[-1], len(NS)))
    say("  fixed inner length runs over a factor %d against the "
        "published 16" % (NS[-1] // NS[0]))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in PU.factor_set(N) if q > 2))
                  for N in NS)))

    NMAX = max(NS)
    qs = [int(q) for q in PU.primes_upto(PU.QSIEVE) if q > 2]
    say()
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    mu = PU.moebius(NMAX)
    oddsqf = (mu != 0)
    oddsqf[::2] = False
    vmask = PU.residue_mask(NMAX, qs)

    data, mx = {}, []
    say()
    say("  N            #k      max ratio   at N/k    k there   "
        "k-range cut by the cap")
    for N in NS:
        ks, inn, r, mass = PU.ratios(N, mu, oddsqf, vmask, qs)
        data[N] = (ks, inn, r, mass)
        j = int(np.argmax(r))
        mx.append(float(r.max()))
        # what the cap removes: the k an uncapped sweep would add
        # below the shortest inner length this N still reaches
        want = N // int(inn.min()) if inn.size else 0
        cut = max(0, want - KCAP)
        say("  %-12d %-7d %-11.4f %-9d %-9d %d"
            % (N, ks.size, float(r.max()), int(inn[j]), int(ks[j]),
               cut))

    x = np.log(np.array(NS, dtype=np.float64))

    # -------------------------------------------------------------- A1
    say()
    say("A1  the control at the published N")
    rnd = 0.5 * 10.0 ** (-dec)
    worst = max(abs(m - pubmx[N]) for N, m in zip(NS, mx)
                if N in pubmx)
    a1 = worst <= rnd
    say("  worst departure %.6f over %d shared N; the table prints "
        "%d decimals, so the bound is %.8f"
        % (worst, len(pubmx), dec, rnd))
    say("PRINTBOUND audit_uniformity_reach %d %.8f" % (dec, rnd))
    say("  A1 %s   (cap: the printing bound)"
        % ("hold" if a1 else "REFUTED"))

    # ------------------------------------------------- the octaves
    say()
    say("A2  the octave table, at a fixed inner length across the "
        "doublings")
    say("  N/k octave          points  ratio at the ends        "
        "slope        s.e.        t")
    rows = {}
    for lo, hi in OCTS:
        xs, ys = [], []
        for N in NS:
            ks, inn, r, mass = data[N]
            sel = (inn >= lo) & (inn < hi)
            if int(sel.sum()) < MINK:
                continue
            xs.append(math.log(N))
            ys.append(float(r[sel].max()))
        if len(xs) < 3:
            continue
        a, rr, se = fit(np.array(xs), np.array(ys))
        rows[(lo, hi)] = (len(xs), ys[0], ys[-1], a, se)
        say("  [%-7d,%-7d)  %-7d %.4f to %-16.4f %+-12.6f %-11.6f "
            "%.2f" % (lo, hi, len(xs), ys[0], ys[-1], a, se,
                      abs(a) / se))
    rising = [oc for oc, v in rows.items()
              if v[3] > 0 and abs(v[3]) / v[4] >= 2.0]
    a2 = all(oc in rising for oc in ((2048, 8192), (8192, 32768))
             if oc in rows)
    say("  resolved rising: %s"
        % (", ".join("[%d,%d)" % oc for oc in rising)
           if rising else "none"))
    say("  the two that were rising on the published set: %s"
        % ", ".join("[%d,%d) %s" % (lo, hi,
                                    "still" if (lo, hi) in rising
                                    else "no longer")
                    for lo, hi in ((2048, 8192), (8192, 32768))))
    say("  A2 %s   (cap 2 standard errors on each)"
        % ("hold" if a2 else "REFUTED"))
    say("SCALES audit_uniformity_reach %d" % len(NS))

    # -------------------------------------------------------------- A3
    say()
    say("A3  does the drift sit where the mass is?")
    shares = []
    for N in NS:
        ks, inn, r, mass = data[N]
        tot = float(mass.sum())
        s = 0.0
        for lo, hi in rising:
            sel = (inn >= lo) & (inn < hi)
            s += float(mass[sel].sum())
        shares.append(s / tot if tot > 0 else 0.0)
    a3 = max(shares) > pubshare
    say("  the resolved-rising octaves carry %.4f to %.4f of "
        "sum(log k)|P| across the doublings"
        % (min(shares), max(shares)))
    say("  against the published %.4f" % pubshare)
    say("CROSSAXIS audit_uniformity_reach %d %d %.4f"
        % (len(rows), len(rising), max(shares)))
    if rising:
        say("AXIS RISE audit_uniformity_reach")
    say("  A3 %s   (cap: the published share)"
        % ("hold" if a3 else "REFUTED"))

    # -------------------------------------------------------------- A4
    say()
    say("A4  and what the bound's constant is doing")
    em, rm, sem = fit(x, np.log(np.array(mx)))
    a4 = em < 0.0 and abs(em) / sem >= 2.0
    say("  the overall maximum runs %.4f down to %.4f over the ten "
        "doublings" % (mx[0], mx[-1]))
    say("  its least-squares slope in log N is %+.6f, s.e. %.6f, "
        "t = %.2f" % (em, sem, abs(em) / sem))
    say("TSTAT slope_audit_uniformity_reach %.2f" % (abs(em) / sem))
    say("SPREAD slope_audit_uniformity_reach %.4f"
        % float(x.max() - x.min()))
    if abs(em) / sem < 2.0:
        say("UNRESOLVED SIGN slope_audit_uniformity_reach")
    say("  A4 %s   (cap 2 standard errors)"
        % ("hold" if a4 else "REFUTED"))

    say()
    say("  what this settles and what it does not. It does not prove "
        "the uniformity,")
    say("  which is a statement about all k and all N. It checks it "
        "on the axis the")
    say("  classical estimate does not control, over a modulus range "
        "%d times longer" % (NS[-1] // NS[0] // 16))
    say("  than the published one, and reports where the drift sits "
        "against where the")
    say("  mass sits. The k-cap of %d is held fixed, so the shortest "
        "inner lengths" % KCAP)
    say("  drop out as N grows; the column above says how many k an "
        "uncapped sweep")
    say("  would have added at each N.")

    say()
    say("=" * 70)
    say("A1 %s  A2 %s  A3 %s  A4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (a1, a2, a3, a4)))

    head = [
        "STATISTIC: the ratio |P(N;k)| / [(N/k) exp(-c sqrt(log(N/k)))",
        "           L(k)] of code/audit_provable_uniformity.py,",
        "           imported unchanged, at every admissible k; its",
        "           overall maximum at each of ten doublings from",
        "           2e5 to 1.024e8; its maximum within each octave of",
        "           the inner length N/k, and the least-squares slope",
        "           of that maximum in log N, which at a fixed octave",
        "           varies the modulus by the ratio of the largest N",
        "           to the smallest; the share of sum(log k)|P|",
        "           carried by the octaves whose slope is resolved",
        "           positive; and the overall maximum's own slope.",
        "NULL: none is run and none applies. A measured sum is divided",
        "      by a deterministic bound and its maxima compared across",
        "      N; there is no background to detect against. The coin",
        "      arm for this field is lab_elementary_provable.py's",
        "      sixteen coins on the identical sifted set.",
        "FIELD: N = 2e5 through 1.024e8 by doubling, every one",
        "       2^a 5^b with one odd radical as RADICALS declares;",
        "       k odd, squarefree and coprime to N with",
        "       2 <= k < 30000, the cap held fixed at",
        "       code/audit_provable_uniformity.py's value; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below 30; c = 0.2098;",
        "       octaves closed at both ends and used only when they",
        "       hold at least 10 k. The published maxima, slopes and",
        "       mass share are read from",
        "       results/audit_provable_uniformity.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not a1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
