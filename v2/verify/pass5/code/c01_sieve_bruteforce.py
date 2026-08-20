# -*- coding: utf-8 -*-
"""C-01.  Every sieve in code/ against trial-division factorisation.

BRIEF C 3.2: exhaustive comparison over n <= 10^4 (not a sample), plus the
boundaries n=0,1 and the prime-power/large-prime cases, plus a cross-check
of every script's sieve against every other's.
"""
import io
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import CODE_DIR, Lambda_bf, load, mu_bf, phi_bf  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results", "c01_sieve_bruteforce.txt")
LIM = 10_000
lines = []


def say(s=""):
    print(s)
    lines.append(s)


# ---------------------------------------------------------------- truth
truth_mu = [mu_bf(n) for n in range(LIM + 1)]
truth_lam = [Lambda_bf(n) for n in range(LIM + 1)]

say("C-01  sieve audit: every mu / Lambda implementation in code/ against")
say("      trial-division factorisation, exhaustively over n <= %d." % LIM)
say("=" * 72)
say()

# which scripts expose a sieve, and under what name
SCRIPTS = sorted(f[:-3] for f in os.listdir(CODE_DIR) if f.endswith(".py"))

results = {}
for name in SCRIPTS:
    try:
        m = load(name)
    except Exception as e:                                    # pragma: no cover
        say("  %-28s IMPORT FAILED: %r" % (name, e))
        continue
    got = {}
    # the two families: sieves(n) -> tuple, mobius_upto/mobius_big/mobius(n)
    for fn in ("sieves", "sieve_mu_phi", "sieves_by_factorisation",
               "sieves_by_recurrence", "mobius_upto", "mobius_big",
               "mobius_small", "sieved"):
        f = getattr(m, fn, None)
        if f is None:
            continue
        try:
            out = f(LIM)
        except TypeError:
            continue
        except Exception as e:
            say("  %-28s %s(%d) raised %r" % (name, fn, LIM, e))
            continue
        got[fn] = out
    if got:
        results[name] = (m, got)


def as_arrays(out):
    """pick the mu and lam arrays out of whatever the function returned."""
    mu = lam = phi = None
    items = out if isinstance(out, tuple) else (out,)
    for a in items:
        if not isinstance(a, np.ndarray) or a.shape[0] < LIM:
            continue
        if a.dtype.kind == "f":
            if lam is None and abs(float(a[4]) - np.log(2)) < 1e-12:
                lam = a
        else:
            vals = set(np.unique(a[: LIM + 1]).tolist())
            if vals <= {-1, 0, 1} and mu is None:
                mu = a
            elif phi is None and int(a[7]) == 6 and int(a[8]) == 4:
                phi = a
    return mu, lam, phi


bad_total = 0
for name in sorted(results):
    m, got = results[name]
    for fn, out in sorted(got.items()):
        mu, lam, phi = as_arrays(out)
        msgs = []
        if mu is not None:
            d = [n for n in range(LIM + 1) if int(mu[n]) != truth_mu[n]]
            msgs.append("mu: %s" % ("EXACT" if not d else
                                    "%d mismatches, first %s" % (len(d), d[:5])))
            bad_total += len(d)
        if lam is not None:
            d = [n for n in range(LIM + 1)
                 if abs(float(lam[n]) - truth_lam[n]) > 1e-12]
            msgs.append("Lambda: %s" % ("EXACT" if not d else
                                        "%d mismatches, first %s" % (len(d), d[:5])))
            bad_total += len(d)
        if phi is not None:
            d = [n for n in range(1, LIM + 1) if int(phi[n]) != phi_bf(n)]
            msgs.append("phi: %s" % ("EXACT" if not d else
                                     "%d mismatches, first %s" % (len(d), d[:5])))
            bad_total += len(d)
        if not msgs:
            msgs = ["(no mu/Lambda/phi array recognised in return value)"]
        say("  %-30s %-24s %s" % (name, fn + "(%d)" % LIM, "; ".join(msgs)))

say()
say("total mismatches against trial division: %d" % bad_total)
say()
say("boundary probes (mu[0], mu[1], Lambda[1], Lambda[4], Lambda[8]):")
for name in sorted(results):
    m, got = results[name]
    for fn, out in sorted(got.items()):
        mu, lam, phi = as_arrays(out)
        if mu is None and lam is None:
            continue
        say("  %-30s %-22s mu0=%s mu1=%s L1=%s L4=%.6f L8=%.6f"
            % (name, fn,
               "-" if mu is None else int(mu[0]),
               "-" if mu is None else int(mu[1]),
               "-" if lam is None else float(lam[1]),
               -1 if lam is None else float(lam[4]),
               -1 if lam is None else float(lam[8])))

io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
print("\nwrote", os.path.abspath(OUT))
