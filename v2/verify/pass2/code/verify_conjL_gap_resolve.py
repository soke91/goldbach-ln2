# -*- coding: utf-8 -*-
r"""
pass2, resolution: the 50 pairs that verify_conjL_mask_zeros.py found

WHAT HAPPENED

verify_conjL_mask_zeros.py re-derived conj:L's annihilation rule from
the conjecture's own statement and got **116000** annihilated pairs on
the grid the repository's stamp describes, against the stamp's
**115950**.  Its W1 was pre-registered to gate on that integer and it
is REFUTED.  The protocol for this tree says a disagreement is the
finding and neither side is assumed right until it is resolved.  This
script resolves it.

WHY THE FIELD CAN DECIDE

One direction of the rule is not a guess.  If q^2 divides N and q^2
divides k then q^2 divides N - pk for every p whatever, so mu(N - pk)
is zero for every p and the support is empty.  That is a proof, not a
measurement, and the field must obey it.  So either

  (a) the 116000 pairs really do all have empty support, and the
      stamp's 115950 undercounts on the grid its own note describes;
      or
  (b) the grid reconstructed here is not the grid the stamp used, and
      the 50 is the difference between two grids rather than an error.

Only (b) leaves the stamp standing, so (b) is the one that has to be
tried hardest.  The reconstruction used the only two facts the stamp
publishes about its grid -- "401 values of N, 1000 values of k" and
the note "every N here is divisible by 4 (3e6 and the step 2500 both
are)".  Neither pins the offset of either range, so this script
enumerates the grids consistent with those two facts and asks whether
any of them gives 115950.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  Every one of the 116000 pairs the derived rule predicts has
      empty support in the field built here.  This is the theorem
      above and a single counterexample refutes the derivation
      outright.
  R2  The gap is a whole number of contributions from single moduli,
      so listing the annihilated count per prime q, and the count
      that q alone contributes, will show where 50 could come from
      rather than leaving it unaccounted.
  R3  **No grid consistent with the stamp's two published facts gives
      115950.**  The scan covers k-ranges [k0, k0+1000) for k0 in
      1..64 and N-grids 3*10^6 + 2500 j offset by j0 in 0..64, both
      with the divisibility note respected.
  R4  So the resolution is (a): the derived rule is right, the stamp
      undercounts by 50, and because the stamp also reports
      "predicted-but-nonzero 0", those 50 pairs are ones its own
      prediction missed and its own check therefore could not see.

REFUTATION RULE (fixed before the run)

  R1  REFUTED by one pair with nonempty support.  Then the derivation
      is wrong, the stamp stands, and this pass reports its own error
      as the finding.
  R2  REFUTED if the per-q table cannot be made to sum to the totals,
      which would mean the accounting here is wrong.  This one is
      bookkeeping and gates nothing.
  R3  **REFUTED if any scanned grid gives exactly 115950.**  Then the
      answer is (b), the stamp stands, and the finding is that the
      stamp's grid is under-specified rather than wrong -- which is
      still worth recording, because a stamp that cannot be
      reconstructed from what it publishes is not independently
      checkable.
  R4  REFUTED if R1 or R3 is.  It has no content of its own; it is
      the reading, and it is written down so that it cannot be
      chosen after the fact.

  The scan in R3 is finite and cannot prove that no grid whatever
  gives 115950.  What it can do is fail to find one, and the honest
  statement of that outcome is "no grid consistent with the published
  facts and inside this scan", not "no grid".  That limit is stated
  before the run and is repeated in whatever the result turns out to
  be.

  NO NULL IS RUN and none applies: every quantity is an exact integer
  count of a deterministic set.
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
OUT = os.path.join(HERE, "..", "results", "verify_conjL_gap_resolve.txt")

PUBLISHED_ZEROS = 115950
PUBLISHED_PAIRS = 401000
NN, NK = 401, 1000
PMAX = 2000
SCAN = 64


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

    ks = np.arange(1, NK + 1, dtype=np.int64)
    Ns = np.arange(0, NN, dtype=np.int64) * 2500 + 3_000_000
    qs = qlist(NK + SCAN)
    kill = derived_rule(Ns, ks, qs)
    n_kill = int(kill.sum())
    say("the grid reconstructed from the stamp's two published facts")
    say("  N = 3e6 + 2500 j, j = 0..%d and k = 1..%d" % (NN - 1, NK))
    say("  derived %d annihilated, published %d, gap %d"
        % (n_kill, PUBLISHED_ZEROS, n_kill - PUBLISHED_ZEROS))
    say("COUNT gap_pairs %d" % (n_kill - PUBLISHED_ZEROS))

    # -------------------------------------------------------------- R1
    say()
    say("R1  does the field obey the derived rule?")
    ps = primes_upto(PMAX)
    ps = ps[ps > 2]
    mu = mobius_upto(int(Ns.max()))
    amu = np.abs(mu).astype(np.int32)
    sup = np.zeros((NN, NK), dtype=np.int32)
    for i, N in enumerate(Ns):
        idx = int(N) - np.outer(ps, ks)
        ok = idx >= 1
        sup[i] = (amu[np.where(ok, idx, 0)] * ok).sum(axis=0)
    obs = sup == 0
    bad = int((kill & ~obs).sum())
    extra = int((~kill & obs).sum())
    say("  over odd primes p <= %d (%d of them), a range this pass "
        "chose" % (PMAX, len(ps)))
    say("  predicted %d, observed empty %d" % (n_kill, int(obs.sum())))
    say("  predicted-but-nonzero %d, unpredicted zeros %d"
        % (bad, extra))
    r1 = bad == 0
    say("  R1 %s   (cap: a single pair)" % ("hold" if r1 else "REFUTED"))
    say("COUNT resolve_predicted_nonzero %d" % bad)

    # -------------------------------------------------------------- R2
    say()
    say("R2  where the annihilations come from, modulus by modulus")
    say("      q   q^2 | k   q^2 | N   pairs by q   q alone")
    tot = 0
    for q in qs:
        a = vq(Ns, q)
        b = vq(ks, q)
        m = (np.minimum(a[:, None], b[None, :]) >= 2)
        if q == 2:
            m = m | ((a[:, None] == b[None, :]) & (a[:, None] >= 1))
        if not m.any():
            continue
        others = np.zeros_like(m)
        for r in qs:
            if r == q:
                continue
            ar, br = vq(Ns, r), vq(ks, r)
            mr = (np.minimum(ar[:, None], br[None, :]) >= 2)
            if r == 2:
                mr = mr | ((ar[:, None] == br[None, :])
                           & (ar[:, None] >= 1))
            others |= mr
        alone = int((m & ~others).sum())
        tot += alone
        say("  %5d   %6d   %6d   %10d   %7d"
            % (q, int((b >= 2).sum()), int((a >= 2).sum()),
               int(m.sum()), alone))
    say("  sole-cause pairs %d of %d; the rest are annihilated by "
        "more than" % (tot, n_kill))
    say("  one modulus, which is why the columns do not add to the "
        "total")
    r2 = tot <= n_kill
    say("  R2 %s   (cap: the accounting closes)"
        % ("hold" if r2 else "REFUTED"))

    # -------------------------------------------------------------- R3
    say()
    say("R3  is there a grid consistent with the published facts "
        "that gives")
    say("    %d?  scanning k0 in 1..%d and j0 in 0..%d"
        % (PUBLISHED_ZEROS, SCAN, SCAN))
    hits = []
    seen = set()
    for j0 in range(SCAN + 1):
        Nsx = np.arange(j0, j0 + NN, dtype=np.int64) * 2500 + 3_000_000
        if not (vq(Nsx, 2) >= 2).all():
            continue                       # the note must hold
        for k0 in range(1, SCAN + 1):
            ksx = np.arange(k0, k0 + NK, dtype=np.int64)
            c = int(derived_rule(Nsx, ksx, qs).sum())
            seen.add(c)
            if c == PUBLISHED_ZEROS:
                hits.append((j0, k0, c))
    say("  %d grids scanned, %d distinct counts, range %d to %d"
        % (len(seen) and (SCAN + 1) * SCAN, len(seen),
           min(seen), max(seen)))
    say("  grids giving exactly %d: %d" % (PUBLISHED_ZEROS, len(hits)))
    for j0, k0, c in hits[:8]:
        say("    j0 = %d, k0 = %d" % (j0, k0))
    r3 = len(hits) == 0
    say("  R3 %s   (cap: one hit refutes)"
        % ("hold" if r3 else "REFUTED"))
    say("  the scan is finite: it can fail to find a grid, it cannot "
        "prove")
    say("  none exists, and the finding is stated at that strength")
    say("COUNT resolve_grid_hits %d" % len(hits))

    # -------------------------------------------------------------- R4
    say()
    say("R4  which resolution, then")
    r4 = r1 and r3
    if r4:
        say("  (a). The derived rule is a theorem in the direction "
            "that matters,")
        say("  the field obeys it on every one of the %d pairs, and "
            "no grid" % n_kill)
        say("  consistent with what the stamp publishes gives its "
            "number. So the")
        say("  stamp undercounts by %d on the grid its own note "
            "describes."
            % (n_kill - PUBLISHED_ZEROS))
        say("  Its 'predicted-but-nonzero 0' is consistent with "
            "that: a pair its")
        say("  rule never predicted cannot show up as a prediction "
            "that failed.")
    elif not r1:
        say("  Neither. The derivation here is wrong -- the field "
            "has a pair it")
        say("  predicted and did not get -- so this pass reports its "
            "own error")
        say("  and the stamp stands.")
    else:
        say("  (b). A grid consistent with the published facts gives "
            "the stamp's")
        say("  number, so the stamp is right and under-specified "
            "rather than")
        say("  wrong. A stamp that cannot be reconstructed from what "
            "it prints")
        say("  is not independently checkable, and that is the "
            "finding.")
    say("  R4 %s" % ("hold" if r4 else "REFUTED"))

    say()
    say("=" * 70)
    say("R1 %s  R2 %s  R3 %s  R4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (r1, r2, r3, r4)))

    head = [
        "STATEMENT: the gap between verify_conjL_mask_zeros.py's",
        "           derived %d and conj:L's blind-mask stamp of %d"
        % (n_kill, PUBLISHED_ZEROS),
        "           on the grid the stamp's own note describes.",
        "METHOD HERE: the field rebuilt over odd primes p <= %d and"
        % PMAX,
        "           compared pair by pair against the derived rule;",
        "           then a finite scan over the grids consistent with",
        "           the two facts the stamp publishes about its own.",
        "REPOSITORY'S NUMBER: %d of %d, reported exact in both"
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
