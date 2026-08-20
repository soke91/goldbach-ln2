# -*- coding: utf-8 -*-
r"""
pass2, third step: the branch of the mask no stamp on that grid can see

WHY THIS RUNS SEPARATELY

verify_conjL_mask_zeros.py registered W3 and W4 on a branch of the
annihilation rule, and its W1 gate refused to report them: W1 was
REFUTED, and W1's rule said that a refuted W1 means the rule under
test may not be the repository's, so nothing downstream verifies
anything.  verify_conjL_gap_resolve.py then settled what the
disagreement was.  It is not a rule disagreement at all --
**R1 held**: the field has empty support on exactly the 116000 pairs
the derived rule names and on no others, both directions zero.  The
gap of 50 is a grid difference, and R3 found 100 grids consistent
with the two facts the stamp publishes that give the stamp's own
115950.

So W3 and W4 can be asked, but not as verifications of the
repository's stamp -- that is what W1's gate correctly forbade, and
this script does not undo it.  They are asked here as a **new
measurement about the mask**, and reported as such.

THE BRANCH

The field is D(N,k) = sum over odd primes p of mu(N - pk).  Write
a = v_2(N), b = v_2(k).  If a = b >= 1 then N = 2^a N1 and k = 2^a k1
with N1, k1 odd, so N1 - p k1 is even for every odd p and
v_2(N - pk) >= a + 1 >= 2: the support is empty.  When a >= 2 this is
already inside "min(a,b) >= 2".  When **a = b = 1** it is not, and it
is not of the form "q^2 divides gcd(N,k)" for any q.

Every N on the repository's grid is divisible by 4.  So on that grid
a >= 2 always, the branch never fires, and the two rules agree on all
401000 pairs -- a stamp there cannot distinguish them however exact
it is.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  B1  On the repository's grid the derived rule and the weaker rule
      "some q^2 divides gcd(N,k)" agree on every pair.
  B2  On a grid with 2 || N they disagree, and the field sides with
      the derived rule: zero disagreements against it, and a positive
      number against the weaker one.
  B3  The branch is not a rounding error: on the 2 || N grid it
      annihilates at least a tenth of all pairs by itself.
  B4  It is the same phenomenon that makes the repository's own note
      true -- "4 | gcd(k,N) alone annihilates a quarter of the k" --
      one step down in the 2-adic valuation, so on the 2 || N grid the
      branch's share should be about a quarter as well.

REFUTATION RULE (fixed before the run)

  B1  REFUTED by one disagreeing pair; then the grid does see the
      branch and there is nothing here to report.
  B2  REFUTED if the two rules agree on the 2 || N grid (the branch is
      empty), or if the field disagrees with the derived rule on any
      pair (the derivation is wrong and this pass reports its own
      error).
  B3  REFUTED below a tenth.
  B4  REFUTED outside 0.20 to 0.30.  This is the one prediction here
      that could plausibly miss while B1--B3 hold, and missing it
      would mean the branch and the note's quarter are not the same
      mechanism seen twice.

  None of these can come out "too noisy to tell": every quantity is
  an exact integer count of a deterministic set, which is why this
  stamp was worth re-verifying at all.

  WHAT THIS DOES NOT CLAIM.  It does not claim the repository's mask
  omits the branch.  This pass is blind and has not read v1/code/;
  what it establishes is that **no stamp on a grid with 4 | N can
  tell whether the branch is there**, so the single witness OPEN.md
  item 4 names is a witness to less than the whole rule.

  NO NULL IS RUN and none applies.
"""

import io
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.abspath(os.path.dirname(__file__))
OUT = os.path.join(HERE, "..", "results", "verify_conjL_mask_branch.txt")

NN, NK, PMAX = 401, 1000, 2000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def mobius_upto(n):
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes_upto(int(n ** 0.5) + 1):
        p = int(p)
        mu[p::p] *= -1
        mu[p * p::p * p] = 0
    rest = np.ones(n + 1, dtype=np.int64)
    for p in primes_upto(int(n ** 0.5) + 1):
        p = int(p)
        rest[p::p] *= p
        q = p * p
        while q <= n:
            rest[q::q] *= p
            q *= p
    big = np.arange(n + 1, dtype=np.int64) // np.maximum(rest, 1)
    mu[(big > 1) & (mu != 0)] *= -1
    return mu


def vq(n, q):
    n = np.array(n, dtype=np.int64, copy=True)
    v = np.zeros(n.shape, dtype=np.int64)
    while True:
        d = (n % q == 0) & (n != 0)
        if not d.any():
            break
        v[d] += 1
        n[d] //= q
    return v


def qlist(kmax):
    return [int(q) for q in primes_upto(int(kmax ** 0.5) + 2)]


def bits(xs, qs):
    b = np.zeros(len(xs), dtype=np.int32)
    v2 = None
    for i, q in enumerate(qs):
        v = vq(xs, q)
        b |= (v >= 2).astype(np.int32) << i
        if q == 2:
            v2 = v
    return b, v2


def rules(Ns, ks, qs):
    bN, v2N = bits(Ns, qs)
    bk, v2k = bits(ks, qs)
    weak = (bN[:, None] & bk[None, :]) != 0
    branch = (v2N[:, None] == v2k[None, :]) & (v2N[:, None] >= 1)
    return weak | branch, weak, branch


def empty_support(Ns, ks, ps, mu):
    sup = np.zeros((len(Ns), len(ks)), dtype=np.int32)
    amu = np.abs(mu).astype(np.int32)
    for i, N in enumerate(Ns):
        idx = int(N) - np.outer(ps, ks)
        ok = idx >= 1
        sup[i] = (amu[np.where(ok, idx, 0)] * ok).sum(axis=0)
    return sup == 0


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    ks = np.arange(1, NK + 1, dtype=np.int64)
    qs = qlist(NK)
    ps = primes_upto(PMAX)
    ps = ps[ps > 2]

    # -------------------------------------------------------------- B1
    say("B1  the repository's grid cannot separate the two rules")
    Np = np.arange(0, NN, dtype=np.int64) * 2500 + 3_000_000
    full_p, weak_p, _ = rules(Np, ks, qs)
    dis_p = int((full_p != weak_p).sum())
    say("  N = 3e6 + 2500 j, j = 0..%d; v_2(N) >= 2 for %d of %d"
        % (NN - 1, int((vq(Np, 2) >= 2).sum()), NN))
    say("  the derived rule and \"some q^2 divides gcd(N,k)\" "
        "disagree on")
    say("  %d of %d pairs" % (dis_p, full_p.size))
    b1 = dis_p == 0
    say("  B1 %s   (cap: one pair)" % ("hold" if b1 else "REFUTED"))
    say("COUNT branch_visible_published %d" % dis_p)

    # -------------------------------------------------------------- B2
    say()
    say("B2  a grid with 2 || N, and which rule the field obeys")
    N2 = np.arange(0, 101, dtype=np.int64) * 10_000 + 3_000_002
    say("  N = 3000002 + 10000 j, j = 0..100; v_2(N) = 1 for %d of %d"
        % (int((vq(N2, 2) == 1).sum()), len(N2)))
    full_2, weak_2, br_2 = rules(N2, ks, qs)
    n2 = full_2.size
    mu = mobius_upto(int(N2.max()))
    obs = empty_support(N2, ks, ps, mu)
    d_full = int((full_2 != obs).sum())
    d_weak = int((weak_2 != obs).sum())
    say("  derived rule annihilates %d of %d (%.4f)"
        % (int(full_2.sum()), n2, full_2.sum() / n2))
    say("  weaker rule annihilates %d of %d (%.4f)"
        % (int(weak_2.sum()), n2, weak_2.sum() / n2))
    say("  field disagrees with the derived rule on %d pairs, with "
        "the weaker" % d_full)
    say("  rule on %d, over odd primes p <= %d" % (d_weak, PMAX))
    b2 = d_full == 0 and d_weak > 0
    say("  B2 %s   (cap: zero against the derived rule, positive "
        "against the" % ("hold" if b2 else "REFUTED"))
    say("        weaker one)")
    say("COUNT branch_field_vs_derived %d" % d_full)
    say("COUNT branch_field_vs_weaker %d" % d_weak)

    # -------------------------------------------------------------- B3
    say()
    say("B3  how much the branch carries on its own")
    extra = int((full_2 & ~weak_2).sum())
    frac = extra / n2
    b3 = frac >= 0.10
    say("  %d of %d pairs, %.4f of the grid" % (extra, n2, frac))
    say("  B3 %s   (cap: a tenth)" % ("hold" if b3 else "REFUTED"))
    say("SHARE branch_alone %.4f" % frac)

    # -------------------------------------------------------------- B4
    say()
    say("B4  is it the repository's own quarter, one valuation down?")
    b4 = 0.20 <= frac <= 0.30
    say("  the repository's note reads \"4 | gcd(k,N) alone "
        "annihilates a")
    say("  quarter of the k\"; on its grid that share is %.4f"
        % (int((vq(ks, 2) >= 2).sum()) / len(ks)))
    say("  the branch's share on the 2 || N grid is %.4f" % frac)
    say("  B4 %s   (cap: 0.20 to 0.30)"
        % ("hold" if b4 else "REFUTED"))

    say()
    say("=" * 70)
    say("B1 %s  B2 %s  B3 %s  B4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (b1, b2, b3, b4)))
    say()
    if b1 and b2:
        say("the mask has a branch that is not \"q^2 divides "
            "gcd(N,k)\", the field")
        say("obeys it exactly, and no stamp on a grid with 4 | N can "
            "tell whether")
        say("a rule contains it. The single witness OPEN.md item 4 "
            "names is a")
        say("witness to less than the whole rule -- which is a "
            "statement about")
        say("what the stamp can show, not a claim that the "
            "repository's mask")
        say("omits the branch. This pass is blind and has not read "
            "v1/code/.")

    head = [
        "STATISTIC: the annihilation rule with and without its",
        "           2-adic branch, each as a set equality against",
        "           the field, and the branch's share of the grid.",
        "FIELD: two grids at 1000 values of k each, with the field",
        "       D(N,k) = sum of mu(N-pk) over the 302 odd primes",
        "       p <= 2000: N = 3*10^6 + 2500 j for j = 0..400,",
        "       where 4 | N throughout, and N = 3000002 + 10000 j",
        "       for j = 0..100, where 2 || N throughout.",
        "STATEMENT: Conjecture conj:L, paper/wall_v3.md -- M is",
        "           computed exactly by finite modular enumeration",
        "           from the v_q-data. This step asks what a stamp on",
        "           a grid with 4 | N can and cannot show about that",
        "           enumeration.",
        "METHOD HERE: the annihilation rule derived from the odd-p",
        "           field, against the same rule with its 2-adic",
        "           branch removed, both checked against the field",
        "           rebuilt over odd primes p <= %d." % PMAX,
        "REPOSITORY'S NUMBER: none. The repository publishes no stamp",
        "           on a grid with 2 || N, which is the finding.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
