# -*- coding: utf-8 -*-
"""C-04.  The FFT route against direct summation.

Nine of the thirty scripts form V, C, U or r~ by zero-padded FFT
convolution at length 2^25 over arrays of 1.6e7 float64 entries, and the
papers quote figures from them to six and seven decimals (V/(W A) = 1.0000002,
sd 0.0001659, cell floors to 1e-3 relative).  The quantities are heavily
cancelled, so the question BRIEF C 3.3 asks is whether float64 FFT at that
length supports the printed precision.

This script recomputes V, W, C, U and r~ at selected N by direct summation
in a different order (prime-power enumeration for V; a single vector product
for C, U, r~) and compares.  It also checks the FFT padding length actually
used against the 2X-1 a linear convolution needs.
"""
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import Lambda_bf, mu_bf  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results", "c04_fft_vs_direct.txt")
X = 16_000_000
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def my_sieves(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    pr = np.flatnonzero(s).astype(np.int64)
    omega = np.zeros(n + 1, dtype=np.int8)
    sqf = np.ones(n + 1, dtype=bool)
    for p in pr:
        p = int(p)
        omega[p::p] += 1
        if p * p <= n:
            sqf[p * p::p * p] = False
    mu = np.where(sqf, np.where(omega & 1, -1, 1), 0).astype(np.int8)
    mu[0] = 0
    lam = np.zeros(n + 1, dtype=np.float64)
    lg = np.log(pr.astype(np.float64))
    lam[pr] = lg
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lg[i]
            if q > n // p:
                break
            q *= p
    return pr, lam, mu


say("C-04  FFT convolution against direct summation")
say("=" * 78)
pr, lam, mu = my_sieves(X)
bad = [n for n in range(10_001)
       if int(mu[n]) != mu_bf(n) or abs(float(lam[n]) - Lambda_bf(n)) > 1e-12]
say("  sieve self-check to 10^4: %d mismatches" % len(bad))
sqf = (mu != 0).astype(np.float64)
lam2 = lam * lam
W = np.cumsum(lam2)

# ---- direct: V(N) = sum_{v<N} mu^2(v) Lambda(N-v)^2, by vector product
def V_direct(N):
    return float((sqf[1:N] * lam2[N - np.arange(1, N, dtype=np.int64)]).sum())


def V_direct_fsum(N):
    """same sum, math.fsum over the nonzero terms -- exactly rounded."""
    idx = np.flatnonzero(lam2[1:N]) + 1          # w = N - v ranges over these
    return math.fsum(float(lam2[w]) for w in idx if mu[N - w] != 0)


def C_direct(N):
    return float((lam[1:N] * mu[N - np.arange(1, N, dtype=np.int64)]).sum())


def r_direct(N):
    return float((lam[1:N] * lam[N - np.arange(1, N, dtype=np.int64)]).sum())


# ---- the FFT route, exactly as the target scripts build it
n = 1
while n < 2 * (X + 1):
    n <<= 1
say("  FFT length used by the targets: 2^%d = %d;  a linear convolution of"
    % (int(math.log2(n)), n))
say("  two length-%d arrays needs %d.  padding is %s"
    % (X + 1, 2 * (X + 1) - 1, "sufficient" if n >= 2 * (X + 1) - 1 else "SHORT"))
a = np.zeros(n, dtype=np.float64)
a[:X + 1] = lam2
F2 = np.fft.rfft(a)
a[:] = 0.0
a[:X + 1] = sqf
V_fft = np.fft.irfft(F2 * np.fft.rfft(a), n)[:X + 1]
del F2
a[:] = 0.0
a[:X + 1] = lam
FL = np.fft.rfft(a)
a[:] = 0.0
a[:X + 1] = mu
C_fft = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
a[:] = 0.0
a[:X + 1] = lam
r_fft = np.fft.irfft(FL * FL, n)[:X + 1]
del a, FL

say()
say("  N            V direct          V fft             rel diff    "
    "V fsum rel diff")
say("  " + "-" * 76)
for N in (1_000_000, 4_000_000, 8_000_000, 15_999_998):
    vd = V_direct(N)
    vf = float(V_fft[N])
    vs = V_direct_fsum(N)
    say("  %-12d %.6f  %.6f  %9.2e  %9.2e"
        % (N, vd, vf, abs(vd - vf) / vd, abs(vs - vf) / vs))

say()
say("  N            C direct          C fft             abs diff    "
    "|C| relative")
say("  " + "-" * 76)
for N in (1_000_000, 4_000_000, 8_000_000, 15_999_998):
    cd = C_direct(N)
    cf = float(C_fft[N])
    say("  %-12d %+.6f  %+.6f  %9.2e  %9.2e"
        % (N, cd, cf, abs(cd - cf), abs(cd - cf) / max(abs(cd), 1e-30)))

say()
say("  N            r~ direct         r~ fft            rel diff")
say("  " + "-" * 76)
for N in (1_000_000, 4_000_000, 15_999_998):
    rd = r_direct(N)
    rf = float(r_fft[N])
    say("  %-12d %.6f  %.6f  %9.2e" % (N, rd, rf, abs(rd - rf) / rd))

say()
say("  cross-check between two scripts in the corpus, on the same number:")
say("    audit_amplification.txt prints V(4e6) = 44684177.8625 (direct,")
say("      prime-power enumeration, its own mu sieve)")
say("    lab_second_moment.txt prints W(4e6)/V(4e6) = 1.270800 (FFT)")
say("    this script, direct : V(4e6) = %.4f" % V_direct(4_000_000))
say("    this script, fsum   : V(4e6) = %.4f" % V_direct_fsum(4_000_000))
say("    this script, FFT    : V(4e6) = %.4f" % float(V_fft[4_000_000]))
say("    W(4e6) = %.4f, so W/V = %.6f  (published 1.270800)"
    % (float(W[3_999_999]), float(W[3_999_999]) / float(V_fft[4_000_000])))
say("    note W(N) is sum_{w<N}, i.e. W[N-1]: using W[N] instead gives "
    "%.6f" % (float(W[4_000_000]) / float(V_fft[4_000_000])))

io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
print("\nwrote", os.path.abspath(OUT))
