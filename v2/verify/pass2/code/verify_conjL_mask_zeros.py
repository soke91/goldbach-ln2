# -*- coding: utf-8 -*-
r"""
pass2, blind: does conj:L's mask predict its own annihilations?

TARGET

Conjecture [conj:L] in paper/wall_v3.md says every Moebius family in
the program factorizes as field = M x G, where M is "the deterministic
local mask, computed exactly by finite modular enumeration from the
$v_q$-data of $(N,k,k')$".  The paper then separates its evidence:
E1 has been reproduced independently, and

    "The remaining stamps -- pair statistics, the exact $(v_2,v_3)$
     cells, the blind mask prediction, and the Wishart pair spectrum
     -- are recorded in the repository and have *not* been
     independently re-verified."

OPEN.md's wall item 4 is exactly that sentence.  This pass takes the
sharpest of the four: **the blind mask prediction**, whose repository
stamp reads

    zero-support k: 115950 of 401000 (28.915%)
    M.3 check: predicted 115950, observed 115950,
               predicted-but-nonzero 0, unpredicted zeros 0
    M.3 EXACT

on the prime-indexed field, over 401 values of N and 1000 values of k.
It is the sharpest because it is not a statistic: it is a
deterministic set equality, so a second implementation either lands on
the same set or does not.

WHAT IS RE-DERIVED HERE, AND FROM WHAT

Nothing is read from v1/code/ and nothing is imported from v2/code/.
The annihilation rule is derived from the statement -- M is a function
of the $v_q$-data -- and from the field's own definition, here

    D(N,k) = sum over primes p in the declared range of mu(N - p k),

whose support is empty exactly when q^2 | N - pk for every p in range,
for some prime q.  Fix q and write a = v_q(N), b = v_q(k).  Because p
runs over ODD primes:

  q odd.  If min(a,b) >= 2 then q^2 | N - pk for every p.  If
          min(a,b) <= 1 then v_q(N - pk) = min(a,b) except when
          a = b, and a = b = 1 forces p into one residue class mod q,
          which the primes do not all satisfy.  So: annihilate iff
          min(a,b) >= 2.

  q = 2.  Same, with one more branch.  If a = b >= 1, write
          N = 2^a N1 and k = 2^a k1 with N1, k1 odd; then N1 - p k1 is
          even for EVERY odd p, so v_2(N - pk) >= a + 1 >= 2 and the
          field is annihilated.  So: annihilate iff min(a,b) >= 2
          **or** a = b >= 1.

The q = 2 branch a = b = 1 is not of the form "q^2 divides gcd(N,k)".
**It is invisible on the repository's grid**, whose own note records
that every N there is divisible by 4 -- with a >= 2 the extra branch
collapses into min(a,b) >= 2 and the two rules agree exactly.  That is
the reason W3 below exists.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  THE GATE.  On the repository's grid -- N = 3*10^6 + 2500 j for
      j = 0..400 and k = 1..1000 -- the rule derived above predicts
      **115950** annihilated pairs of 401000, the published count, to
      the integer.
  W2  And the rule is exact in both directions against the field
      actually built here: zero predicted-but-nonzero and zero
      unpredicted zeros, reproducing "M.3 EXACT".
  W3  **The repository's grid cannot see one branch of the mask.**
      On that grid the derived rule and the weaker rule "some q^2
      divides gcd(N,k)" agree on every one of the 401000 pairs; on a
      grid with 2 || N they disagree, the derived rule predicts more
      annihilations, and the field agrees with the derived rule and
      not with the weaker one.
  W4  On the 2 || N grid the extra branch is not a rounding error: it
      annihilates at least a tenth of all pairs by itself.

REFUTATION RULE (fixed before the run)

  W1  REFUTED at any count other than 115950.  Then the rule derived
      here is not the rule the repository used, and the disagreement
      is the finding rather than either side being right.  GATES: if
      W1 fails, W2--W4 are not reported as verifications of anything,
      because they would be testing a different rule.
  W2  REFUTED by a single pair in either direction.  A predicted-but-
      nonzero pair refutes the derivation; an unpredicted zero refutes
      only the converse, which is the empirical half and depends on
      the prime range -- if that half fails it is reported as a
      property of the range chosen here and NOT as a fault in the
      repository's stamp, whose range is not published.
  W3  REFUTED if the two rules disagree on the repository's grid
      (then the grid does see the branch and there is nothing to
      report), or if they agree on the 2 || N grid (then the branch
      is empty and the derivation has a redundant clause), or if the
      field sides with the weaker rule (then the derivation is
      wrong).
  W4  REFUTED below a tenth.  This one cannot come out "too noisy to
      tell": every quantity in W1--W4 is an exact integer count of a
      deterministic set, so the failure mode of an underpowered
      measurement does not exist here.  That is the whole reason this
      stamp was chosen over the other three.

  NO NULL IS RUN and none applies.  There is no signal to detect
  against a background: the claim is a set equality and the check is
  whether two sets are equal.
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
OUT = os.path.join(HERE, "..", "results", "verify_conjL_mask_zeros.txt")

PUBLISHED_ZEROS = 115950          # published in the repository stamp
PUBLISHED_PAIRS = 401000          # published in the repository stamp

KMAX = 1000
PMAX = 2000                       # this pass's own choice; see W2


def mobius_upto(n):
    """mu on 0..n by a sieve written here"""
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes_upto(int(n ** 0.5) + 1):
        mu[int(p)::int(p)] *= -1
        mu[int(p) * int(p)::int(p) * int(p)] = 0
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


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def vq(n, q):
    """q-adic valuation, elementwise on an int array."""
    n = np.array(n, dtype=np.int64, copy=True)
    v = np.zeros(n.shape, dtype=np.int64)
    live = n != 0
    while live.any():
        d = (n % q == 0) & live
        if not d.any():
            break
        v[d] += 1
        n[d] //= q
        live = d
    return v


def qlist(kmax):
    """q can only matter if q^2 divides some k in range"""
    return [int(q) for q in primes_upto(int(kmax ** 0.5) + 2)]


def _bits(xs, qs):
    """bit q set when v_q(x) >= 2; also return v_2"""
    b = np.zeros(len(xs), dtype=np.int32)
    v2 = None
    for i, q in enumerate(qs):
        v = vq(xs, q)
        b |= (v >= 2).astype(np.int32) << i
        if q == 2:
            v2 = v
    return b, v2


def derived_rule(Ns, ks, qs):
    """annihilate iff min(a,b)>=2 for some q, or a=b>=1 for q=2"""
    bN, v2N = _bits(Ns, qs)
    bk, v2k = _bits(ks, qs)
    kill = (bN[:, None] & bk[None, :]) != 0
    kill |= (v2N[:, None] == v2k[None, :]) & (v2N[:, None] >= 1)
    return kill


def weaker_rule(Ns, ks, qs):
    """the rule without the odd-p branch: some q^2 divides gcd(N,k)"""
    bN, _ = _bits(Ns, qs)
    bk, _ = _bits(ks, qs)
    return (bN[:, None] & bk[None, :]) != 0


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    ks = np.arange(1, KMAX + 1, dtype=np.int64)
    qs = qlist(KMAX)

    # ------------------------------------------------------------- W1
    say("W1  the derived rule against the published count")
    Ns_pub = np.arange(0, 401, dtype=np.int64) * 2500 + 3_000_000
    kill = derived_rule(Ns_pub, ks, qs)
    n_pairs = kill.size
    n_kill = int(kill.sum())
    say("  grid N = 3e6 + 2500 j, j = 0..400, and k = 1..%d" % KMAX)
    say("  pairs here %d against the published %d"
        % (n_pairs, PUBLISHED_PAIRS))
    say("  annihilated here %d against the published %d"
        % (n_kill, PUBLISHED_ZEROS))
    say("  rate %.5f%%" % (100.0 * n_kill / n_pairs))
    w1 = (n_pairs == PUBLISHED_PAIRS and n_kill == PUBLISHED_ZEROS)
    say("  W1 %s   (cap: the published integer, exactly)"
        % ("hold" if w1 else "REFUTED"))
    say("COUNT conjL_zeros_derived %d" % n_kill)
    say("COUNT conjL_zeros_published %d" % PUBLISHED_ZEROS)
    if not w1:
        say()
        say("  W1 gates. The rule derived here is not the rule the")
        say("  repository used, so W2 to W4 are not reported.")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------------- W2
    say()
    say("W2  the rule against the field, both directions")
    ps = primes_upto(PMAX)
    ps = ps[ps > 2]
    say("  odd primes p <= %d, %d of them; this pass's own choice, "
        "the" % (PMAX, len(ps)))
    say("  repository's range is not published")
    mu = mobius_upto(int(Ns_pub.max()))
    sup = support(Ns_pub, ks, ps, mu)
    obs = sup == 0
    bad_a = int((kill & ~obs).sum())      # predicted, yet nonzero
    bad_b = int((~kill & obs).sum())      # zero, yet unpredicted
    say("  observed empty support %d, predicted %d" % (int(obs.sum()),
                                                       n_kill))
    say("  predicted-but-nonzero %d, unpredicted zeros %d"
        % (bad_a, bad_b))
    w2 = bad_a == 0 and bad_b == 0
    say("  W2 %s   (cap: a single pair either way)"
        % ("hold" if w2 else "REFUTED"))
    say("COUNT conjL_predicted_nonzero %d" % bad_a)
    say("COUNT conjL_unpredicted_zero %d" % bad_b)

    # ------------------------------------------------------------- W3
    say()
    say("W3  the branch the published grid cannot see")
    weak_pub = weaker_rule(Ns_pub, ks, qs)
    same_pub = int((weak_pub != kill).sum())
    say("  on the published grid the derived rule and the weaker "
        "rule")
    say("  \"some q^2 divides gcd(N,k)\" disagree on %d of %d pairs"
        % (same_pub, n_pairs))
    say("  because every N there is divisible by 4: v_2(N) >= 2 "
        "holds for")
    say("  %d of %d, so the branch a = b = 1 never fires"
        % (int((vq(Ns_pub, 2) >= 2).sum()), len(Ns_pub)))

    Ns_2 = np.arange(0, 101, dtype=np.int64) * 10_000 + 3_000_002
    say("  now a grid with 2 || N: N = 3000002 + 10000 j, j = 0..100")
    say("  v_2(N) = 1 for %d of %d"
        % (int((vq(Ns_2, 2) == 1).sum()), len(Ns_2)))
    kill2 = derived_rule(Ns_2, ks, qs)
    weak2 = weaker_rule(Ns_2, ks, qs)
    n2 = kill2.size
    say("  derived rule annihilates %d of %d (%.3f%%)"
        % (int(kill2.sum()), n2, 100.0 * kill2.sum() / n2))
    say("  weaker rule annihilates %d of %d (%.3f%%)"
        % (int(weak2.sum()), n2, 100.0 * weak2.sum() / n2))
    sup2 = support(Ns_2, ks, ps, mu)
    obs2 = sup2 == 0
    d_bad = int((kill2 != obs2).sum())
    w_bad = int((weak2 != obs2).sum())
    say("  the field disagrees with the derived rule on %d pairs "
        "and with" % d_bad)
    say("  the weaker rule on %d" % w_bad)
    w3 = (same_pub == 0 and int((weak2 != kill2).sum()) > 0
          and d_bad == 0 and w_bad > 0)
    say("  W3 %s   (cap: agree on the published grid, disagree on "
        "this one," % ("hold" if w3 else "REFUTED"))
    say("        and the field sides with the derived rule)")
    say("COUNT conjL_branch_missed %d" % w_bad)

    # ------------------------------------------------------------- W4
    say()
    say("W4  how large the invisible branch is")
    extra = int((kill2 & ~weak2).sum())
    frac = extra / n2
    w4 = frac >= 0.10
    say("  the branch alone annihilates %d of %d pairs, %.4f of the "
        "grid" % (extra, n2, frac))
    say("  W4 %s   (cap: a tenth)" % ("hold" if w4 else "REFUTED"))
    say("SHARE conjL_branch %.4f" % frac)

    say()
    say("=" * 70)
    say("W1 %s  W2 %s  W3 %s  W4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (w1, w2, w3, w4)))
    say()
    if w1 and w2:
        say("the stamp is re-verified: a second implementation, "
            "written from")
        say("conj:L's statement and not from the repository's code, "
            "lands on")
        say("the same %d pairs and the same exact agreement with the "
            "field." % n_kill)
    if w3 and w4:
        say("and it carries further than the stamp did. The mask has "
            "a branch")
        say("that the published grid holds fixed, so the stamp "
            "verified only")
        say("part of the rule -- correctly, but not the part that "
            "distinguishes")
        say("the odd-p field from an all-m one.")

    head = [
        "STATEMENT: Conjecture conj:L, paper/wall_v3.md -- M is the",
        "           deterministic local mask, computed exactly by",
        "           finite modular enumeration from the v_q-data of",
        "           (N,k,k'). The stamp under test is the blind mask",
        "           prediction, which the paper lists as not",
        "           independently re-verified.",
        "METHOD HERE: the annihilation rule derived from the",
        "           statement and from p being an odd prime, not read",
        "           from v1/code/; mu by a least-prime-factor sieve",
        "           written in this file; the field rebuilt over odd",
        "           primes p <= %d." % PMAX,
        "REPOSITORY'S NUMBER: %d annihilated of %d, exact in both"
        % (PUBLISHED_ZEROS, PUBLISHED_PAIRS),
        "           directions.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
