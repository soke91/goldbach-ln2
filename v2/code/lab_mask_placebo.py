# -*- coding: utf-8 -*-
r"""
The placebo control for "the mask exists" -- wall_v3.md, Section
{#sec:floor}, run against Lemma {#lem:placebo}.

WHY THIS HAS TO BE RUN

Section [sec:floor] reports that the cell means clear the exact floor
by a wide margin -- max_c |z_c| between 9.1 and 13.0 over every octave,
carried by the deep cells -- and adds that "under the permutation of
Lemma [lem:placebo] the floor collapses to the independent-sign value,
so what is being detected is the correspondence between cells and
divisibility and not the cell sizes."

That permutation control is CLAIMED and was never run: lab_cell_floor
computes the floor and the z-scores but states in its own NULL: line
that the placebo is not run there.  After the level measurement of
Remark [rem:levelmeas] turned out to be a support statement once its
null was run, the remaining un-nulled detection claim in either paper
is this one, and it is the strongest positive claim in the program.

The control: replace the labelling l(N) by l(pi(N)) for a random
permutation pi of the even N in the band.  Cell sizes are preserved
exactly; the correspondence between a cell and the arithmetic of its
members is destroyed.  Everything else -- Z(N), V(N), the band, the
exact floor formula of Lemma [lem:cellmom] -- is held identical, and
the floor is RECOMPUTED for each permutation rather than reused, since
u_c(v) depends on which N sit in c.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  L1  With the true labelling on the octave (2e6, 4e6], max_c |z_c| is
      at least 5.
  L2  Over 10 label permutations, max_c |z_c| stays below 5 every time,
      and its mean over permutations is below 3.
  L3  The floor itself is a property of the cell sizes: the standard
      error se_c changes by less than 10% between the true labelling
      and the permutation mean, at every depth.
  L4  The true labelling's z_c is monotone decreasing in depth; under
      permutation the rank correlation between depth and z_c has mean
      near zero, |mean| < 0.3.

REFUTATION RULE (fixed before the run)

  L1  REFUTED if max_c |z_c| < 5 for the true labelling.
  L2  REFUTED if any permutation reaches 5, or if the mean is 3 or
      more.  This is the one that matters: if a permutation reproduces
      the detection, the mask is a property of the cell sizes and
      Section [sec:floor]'s central claim must be withdrawn.
  L3  REFUTED if any depth moves by 10% or more.
  L4  REFUTED if the permuted mean rank correlation has |mean| >= 0.3.

  L1, L2 and L3 gate.  L4 is reported.

BACKS: Proposition {#prop:placebo} and Remark {#rem:floorsignal} in
paper/wall_v3.md, and Lemma {#lem:placebo}, whose control it runs.
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
OUT = os.path.join(ROOT, "results", "lab_mask_placebo.txt")

LO, HI = 2_000_000, 4_000_000
CELLP = (3, 5, 7, 11, 13)
PERMS = 10
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def pow2(n):
    L = 1
    while L < n:
        L <<= 1
    return L


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    X = HI
    say("sieving to %d ..." % X)
    pr = primes_upto(X)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(X + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > X:
            break
        q = p * p
        while q <= X:
            lam[q] = lgp[i]
            if q > X // p:
                break
            q *= p

    mu = np.ones(X + 1, dtype=np.int8)
    rem = np.arange(X + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(X))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
        q = p
        while q <= X:
            rem[q::q] //= p
            if q > X // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    sqf = (mu != 0)

    say("convolving V and C ...")
    n = pow2(2 * (X + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:X + 1] = lam ** 2
    FL2 = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = sqf
    V = np.fft.irfft(FL2 * np.fft.rfft(a), n)[:X + 1]
    del FL2
    a[:] = 0.0
    a[:X + 1] = lam
    FL = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = mu
    C = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
    del a, FL

    depth = np.zeros(X + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    Ns = np.arange(LO + 2, HI + 1, 2, dtype=np.int64)
    Ns = Ns[V[Ns] > 0]
    g = 1.0 / np.sqrt(V[Ns])
    Z = C[Ns] / np.sqrt(V[Ns])
    nb = Ns.size
    zbar = float(Z.mean())
    say("band (%d, %d]: %d even N" % (LO, HI, nb))

    m = pow2(2 * (HI + 1))
    b = np.zeros(m, dtype=np.float64)
    b[:HI + 1] = lam[:HI + 1]
    FLam = np.conj(np.fft.rfft(b))
    w = sqf[:HI + 1].astype(np.float64)

    def ucorr(sel):
        b[:] = 0.0
        b[Ns[sel]] = g[sel]
        return np.fft.irfft(FLam * np.fft.rfft(b), m)[:HI + 1]

    ua = ucorr(np.ones(nb, dtype=bool))
    Qaa = float((w * ua * ua).sum())

    def run(lab):
        """z_c and se_c for a labelling, floor recomputed from scratch."""
        zs, ses, ds, ns = [], [], [], []
        for d in range(6):
            sel = lab == d
            nc = int(sel.sum())
            if nc == 0:
                continue
            uc = ucorr(sel)
            Qcc = float((w * uc * uc).sum())
            Qca = float((w * uc * ua).sum())
            var = Qcc / nc ** 2 - 2.0 * Qca / (nc * nb) + Qaa / nb ** 2
            se = math.sqrt(max(var, 0.0))
            mc = float(Z[sel].mean())
            zs.append((mc - zbar) / se if se > 0 else 0.0)
            ses.append(se)
            ds.append(d)
            ns.append(nc)
        return np.array(ds), np.array(ns), np.array(zs), np.array(ses)

    true_lab = depth[Ns]
    ds, ncs, z_true, se_true = run(true_lab)
    say()
    say("true labelling")
    say("  depth  n_c        se_c          z_c")
    for i in range(len(ds)):
        say("  %-6d %-10d %-13.4e %+.4f" % (ds[i], ncs[i], se_true[i],
                                            z_true[i]))
    mx_true = float(np.abs(z_true).max())
    say("  max |z_c| = %.4f" % mx_true)
    l1 = mx_true >= 5.0
    say("  L1 %s" % ("hold" if l1 else "REFUTED"))

    rng = np.random.default_rng(SEED)
    say()
    say("placebo: %d label permutations (cell sizes preserved exactly)"
        % PERMS)
    say("  draw   max |z_c|   z by depth")
    mxs, rhos, ses_p = [], [], []
    for t in range(PERMS):
        lab = rng.permutation(true_lab)
        d2, n2, z2, s2 = run(lab)
        mxs.append(float(np.abs(z2).max()))
        ses_p.append(s2)
        rhos.append(float(np.corrcoef(d2.astype(float), z2)[0, 1]))
        say("  %-6d %-11.4f %s"
            % (t + 1, mxs[-1], " ".join("%+.2f" % v for v in z2)))
    mxs = np.array(mxs)
    l2 = bool((mxs < 5.0).all() and mxs.mean() < 3.0)
    say("  max over draws = %.4f;  mean = %.4f   (caps 5 and 3)   %s"
        % (mxs.max(), mxs.mean(), "hold" if l2 else "REFUTED"))

    se_p = np.mean(np.array(ses_p), axis=0)
    dev = np.abs(se_p - se_true) / se_true
    l3 = bool((dev < 0.10).all())
    say()
    say("L3   the floor under permutation, by depth")
    say("  depth  se_true       se_perm       |rel dev|")
    for i in range(len(ds)):
        say("  %-6d %-13.4e %-13.4e %.4f"
            % (ds[i], se_true[i], se_p[i], dev[i]))
    say("  L3 %s   (cap 0.10)" % ("hold" if l3 else "REFUTED"))

    rho_true = float(np.corrcoef(ds.astype(float), z_true)[0, 1])
    rm = float(np.mean(rhos))
    l4 = abs(rm) < 0.3
    say()
    say("L4   depth-vs-z correlation: true %.4f, permuted mean %.4f   %s"
        % (rho_true, rm, "hold" if l4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = l1 and l2 and l3
    say("L1 %s  L2 %s  L3 %s  L4 %s"
        % tuple("hold" if v else "REFUTED" for v in (l1, l2, l3, l4)))
    say("the mask survives its placebo: the detection is the "
        "cell-arithmetic correspondence" if ok else
        "REFUTED -- Section {#sec:floor}'s central claim does not "
        "survive the permutation it cites")

    head = [
        "STATISTIC: z_c = (m_c - mbar)/se_c with se_c the exact floor of",
        "           Lemma {#lem:cellmom}, computed for the true depth",
        "           labelling and for 10 random permutations of that",
        "           labelling over the band; max_c |z_c| in each case; the",
        "           floor se_c itself under permutation against the true",
        "           one; and the correlation between depth and z_c.",
        "NULL: this file is the null -- the placebo of Lemma",
        "      {#lem:placebo}, which Section {#sec:floor} cites but which",
        "      no script had run. Cell sizes are preserved exactly and the",
        "      floor is recomputed from scratch for every permutation, so",
        "      the cell-to-arithmetic correspondence is the only thing",
        "      destroyed.",
        "FIELD: the octave (2e6, 4e6], even N with V(N) > 0; cells indexed",
        "       by depth = #{p in 3,5,7,11,13 dividing N}; Lambda, mu and",
        "       the squarefree indicator from an integer sieve to 4e6;",
        "       V = mu^2 * Lambda^2 and C = mu * Lambda by exact FFT",
        "       convolution; u_c by FFT cross-correlation; numpy",
        "       default_rng seed 20260808.",
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
