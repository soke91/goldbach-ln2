# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Proposition {#prop:scaleinv} and the paragraph above
it -- what predicts the floor's size.

WHY THIS HAS TO BE RUN

The paragraph says the floor's magnitude tracks

    D_c := E_same,c[S_2] - E_all[S_2],

with S_2(h) the Hardy-Littlewood singular series of the shift h,
E_same,c the mean over pairs N, N' both in cell c, and E_all the mean
over all pairs in the band; and Proposition [prop:scaleinv] says D_c
depends only on the distribution of h in the residue classes mod
3,5,7,11,13 that the cell fixes, hence is scale-invariant and predicts
an exponent of zero at every depth.

Neither S_2 nor D_c is computed anywhere in this repository.  The
proposition's evidence marker pointed at lab_cell_floor.py, which
computes the floor and the z-scores and never forms S_2 at all -- an
evidence marker can name an existing script that does not compute the
statement, and G1 and G4 only check that the file and its result exist.
So this is the claim in either paper with the weakest backing, and it
is the one that explains the floor.

    S_2(h) = 2 C_2 prod_{p | h, p > 2} (p-1)/(p-2)   for even h,

with C_2 the twin-prime constant.  Cells are indexed by depth, the
number of 3,5,7,11,13 dividing N, so a deeper cell makes p | h more
likely for those p and should raise E_same,c[S_2].

BACKS: Proposition {#prop:scaleinv} and Remark {#rem:copiedfloor}
in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  D_c > 0 at every depth d >= 1, and D_c increases with depth.
  M2  D_c is scale-invariant: over the three octaves (1e6,2e6],
      (2e6,4e6], (4e6,8e6] it varies by less than 2% at every depth,
      beyond sampling error.
  M3  The floor tracks it: across depths, the correlation between the
      exact floor se_c of Lemma [lem:cellmom] and D_c exceeds 0.8.
  M4  Fitting D_c ~ N^{-e} across the three octaves gives |e| < 0.01 at
      every depth -- the exponent zero the proposition predicts.

REFUTATION RULE (fixed before the run)

  M1  REFUTED if D_c <= 0 at some depth >= 1, or if the sequence is not
      increasing in depth.
  M2  REFUTED if the spread (max-min)/mean across octaves exceeds 0.02
      at any depth by more than three sampling standard errors.
  M3  REFUTED if the correlation is 0.8 or below.
  M4  REFUTED if |e| >= 0.01 at any depth.

  M1, M2 and M4 gate -- they are the proposition.  M3 is the
  paragraph's "tracks" claim and is reported.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_cell_singular.txt")

OCTS = [(1_000_000, 2_000_000), (2_000_000, 4_000_000),
        (4_000_000, 8_000_000)]
HMAX = 8_000_000
CELLP = (3, 5, 7, 11, 13)
SAMPLES = 400_000
SEED = 20260808
# The exact floor at (2e6,4e6] is computed by lab_mask_placebo.py and
# is READ from its result file, not copied. A hand-copied table is a
# dependency no check can see: G18 compares a script with its own
# result and G22 compares a script with what it reads, and neither
# sees a number that was typed in. Reading it makes the dependency
# visible and puts it under G22.
def read_floor():
    p = os.path.join(ROOT, "results", "lab_mask_placebo.txt")
    src = io.open(p, encoding="utf-8").read()
    blk = src[src.index("true labelling"):]
    out = {}
    for ln in blk.splitlines()[2:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[2])
    return [out[d] for d in sorted(out)]


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("building S_2(h) to %d ..." % HMAX)
    pr = primes_upto(HMAX)
    twin = 2.0
    for p in pr:
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    S2 = np.zeros(HMAX + 1, dtype=np.float64)
    S2[2::2] = twin
    for p in pr:
        p = int(p)
        if p == 2:
            continue
        S2[p::p] *= (p - 1.0) / (p - 2.0)
    S2[1::2] = 0.0
    S2[0] = 0.0
    say("  2*C_2 = %.6f;  S_2(2) = %.6f;  S_2(6) = %.6f  (should be "
        "2*C_2 * 2)" % (twin, S2[2], S2[6]))

    depth = np.zeros(HMAX + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    rng = np.random.default_rng(SEED)
    say()
    say("  octave              depth  n_c        E_same,c[S_2]   D_c"
        "          se(D_c)")
    say("  " + "-" * 76)
    D = {}
    for lo, hi in OCTS:
        Ns = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
        dep = depth[Ns]
        i1 = rng.integers(0, Ns.size, SAMPLES)
        i2 = rng.integers(0, Ns.size, SAMPLES)
        ok = i1 != i2
        hall = np.abs(Ns[i1[ok]] - Ns[i2[ok]])
        Eall = float(S2[hall].mean())
        for d in range(6):
            idx = np.flatnonzero(dep == d)
            if idx.size < 2:
                continue
            j1 = rng.integers(0, idx.size, SAMPLES)
            j2 = rng.integers(0, idx.size, SAMPLES)
            k = j1 != j2
            hs = np.abs(Ns[idx[j1[k]]] - Ns[idx[j2[k]]])
            v = S2[hs]
            Es = float(v.mean())
            se = float(v.std(ddof=1)) / math.sqrt(v.size)
            D.setdefault(d, []).append(Es - Eall)
            say("  (%9d,%9d] %-6d %-10d %-15.6f %-12.6f %.6f"
                % (lo, hi, d, idx.size, Es, Es - Eall, se))
        say("      E_all[S_2] = %.6f" % Eall)

    say()
    ds = sorted(D)
    m1 = all(D[d][1] > 0 for d in ds if d >= 1) and all(
        D[ds[i]][1] < D[ds[i + 1]][1] for i in range(len(ds) - 1))
    say("M1  D_c at (2e6,4e6] by depth: %s"
        % ", ".join("%.6f" % D[d][1] for d in ds))
    say("    positive for d>=1 and increasing in depth: %s   %s"
        % (m1, "hold" if m1 else "REFUTED"))

    say()
    say("M2/M4  scale invariance")
    say("  depth  D_c per octave                       spread   "
        "fitted exponent")
    m2 = m4 = True
    absp = {}
    for d in ds:
        vals = np.array(D[d])
        absp[d] = (vals.max() - vals.min(), abs(float(vals.mean())))
        spread = (vals.max() - vals.min()) / abs(vals.mean())
        xs = np.log(np.array([0.5 * (lo + hi) for lo, hi in OCTS]))
        e = -float(np.polyfit(xs, np.log(np.abs(vals)), 1)[0])
        if spread > 0.02:
            m2 = False
        if abs(e) >= 0.01:
            m4 = False
        say("  %-6d %-36s %-8.4f %+.6f"
            % (d, ", ".join("%.6f" % v for v in vals), spread, e))
    say("  M2 %s   M4 %s" % ("hold" if m2 else "REFUTED",
                             "hold" if m4 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). M2 and M4 fail only where D_c is small:")
    say("  the spread across octaves is nearly the same size at every")
    say("  depth while D_c itself grows, so the relative spread the")
    say("  rule looks at is large only where D_c is small:")
    say("  depth  spread across octaves   |D_c|      spread/|D_c|")
    for d in ds:
        a, m = absp[d]
        say("  %-6d %-23.6f %-10.6f %.4f" % (d, a, m, a / m))
    say("  Note also")
    say("  that the pre-registration allowed three sampling standard")
    say("  errors and the code applied the 2% band with no allowance --")
    say("  the code is the stricter of the two and its verdict stands.")
    say("  Re-sampling depths 0 and 1 ten times harder:")
    BIG = SAMPLES * 10
    for d in (0, 1):
        vals = []
        for lo, hi in OCTS:
            Ns = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
            dep = depth[Ns]
            i1 = rng.integers(0, Ns.size, BIG)
            i2 = rng.integers(0, Ns.size, BIG)
            ok = i1 != i2
            Eall = float(S2[np.abs(Ns[i1[ok]] - Ns[i2[ok]])].mean())
            idx = np.flatnonzero(dep == d)
            j1 = rng.integers(0, idx.size, BIG)
            j2 = rng.integers(0, idx.size, BIG)
            k = j1 != j2
            v = S2[np.abs(Ns[idx[j1[k]]] - Ns[idx[j2[k]]])]
            vals.append(float(v.mean()) - Eall)
        vals = np.array(vals)
        xs = np.log(np.array([0.5 * (lo + hi) for lo, hi in OCTS]))
        e = -float(np.polyfit(xs, np.log(np.abs(vals)), 1)[0])
        say("    depth %d: %s   spread %.4f   exponent %+.6f"
            % (d, ", ".join("%.6f" % v for v in vals),
               (vals.max() - vals.min()) / abs(vals.mean()), e))

    say()
    dv = np.array([D[d][1] for d in ds])
    SE_TRUE = read_floor()
    sv = np.array([SE_TRUE[d] for d in ds])
    rho = float(np.corrcoef(dv, sv)[0, 1])
    # M3's control. lab_mask_placebo.py permutes depth labels to test
    # the DETECTION z_c; nothing tested this correlation. With six
    # depths the permutation null is exact -- all 720 orderings.
    import itertools
    perm = sorted(abs(float(np.corrcoef(dv, np.array(q))[0, 1]))
                  for q in itertools.permutations(sv.tolist()))
    above = sum(1 for v in perm if v >= abs(rho))
    m3 = rho > 0.8
    say("M3  corr(se_c, D_c) across depths = %.4f   (floor 0.8)   %s"
        % (rho, "hold" if m3 else "REFUTED"))
    say("    se_c  = %s" % ", ".join("%.5f" % v for v in sv))
    say("    control: over all %d permutations of the depth"
        % len(perm))
    say("    labels, %d reach |r| = %.4f; the null median is"
        % (above, abs(rho)))
    say("    %.4f and its 95th percentile %.4f"
        % (perm[len(perm) // 2], perm[int(0.95 * len(perm))]))
    say("    -- the control lab_mask_placebo.py does not"
        " supply, since it permutes labels for z_c and not")
    say("    for this pair.")
    say("    D_c   = %s" % ", ".join("%.5f" % v for v in dv))

    say()
    say("=" * 70)
    ok = m1 and m2 and m4
    say("M1 %s  M2 %s  M3 %s  M4 %s"
        % tuple("hold" if v else "REFUTED" for v in (m1, m2, m3, m4)))
    say("Proposition {#prop:scaleinv} stands" if ok else "REFUTED")

    head = [
        "STATISTIC: D_c = E_same,c[S_2] - E_all[S_2] with",
        "           S_2(h) = 2 C_2 prod_{p|h,p>2}(p-1)/(p-2) the",
        "           Hardy-Littlewood singular series of the shift, E_same,c",
        "           the mean over sampled pairs N != N' both in cell c and",
        "           E_all the mean over sampled pairs in the band; its",
        "           spread across three octaves and its fitted exponent in",
        "           N; and its correlation across depths with the exact",
        "           floor se_c of Lemma {#lem:cellmom}.",
        "NULL: an exact permutation null on M3. D_c itself is a",
        "      deterministic arithmetic functional of the cell with no",
        "      sign input, so Lemma {#lem:coin} does not bite on it; but",
        "      M3 is a correlation across depths between D_c and the",
        "      exact floor, and that needed a control which nothing",
        "      supplied. With six depths the null is enumerable: all 720",
        "      permutations of the depth labels are run and the observed",
        "      |r| is placed in that distribution. The sampling error of",
        "      the pair average is reported alongside.",
        "FIELD: octaves (1e6,2e6], (2e6,4e6], (4e6,8e6]; even N; cells",
        "       indexed by depth = #{p in 3,5,7,11,13 dividing N};",
        "       400000 sampled ordered pairs per cell and per band, numpy",
        "       default_rng seed 20260808; S_2 sieved to 8e6; se_c taken",
        "       from lab_mask_placebo.py at (2e6,4e6].",
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
