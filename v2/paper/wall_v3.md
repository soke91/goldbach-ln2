# The demand side of the Elliott–Halberstam route to binary Goldbach,  and the wall's exact second moment

```latex
% 수식이 쓰는 매크로 — 렌더러/역변환용
\renewcommand{\SS}{\mathfrak{S}}
\newcommand{\AAA}{\mathfrak{A}}
```

## Abstract

Huang and Li proved that binary Goldbach for large even $N$ follows
from $EH_\mu(N^{\theta'})$ for a single $\theta' > 1/2$. This paper
maps that hypothesis from both of its sides and establishes what is
exactly computable about the single scalar it reduces to.

On the *demand* side we prove one unconditional theorem
(Theorem [thm:A]): the Möbius-weighted, fixed-class correlation
sum with weight $w_k = 1$ is $\ll_A N(\log N)^{-A}$, so Huang–Li's
$E_4$ consumption of $EH_\mu$ is unnecessary and the whole demand
collapses to the single scalar $E_3$. An unconditional identity
(Theorem [thm:C]) then evaluates that scalar: $E_3$ is exactly the
defect in Huang–Li's equation (22), so $E_3 \ll_A N(\log N)^{-A}$ is
equivalent to (22) itself, yields binary Goldbach through their
Theorem 1, and yields the asymptotic $\tilde r(N)\sim\SS(N)N$ precisely
when $C(N) = o(N)$ — which is the wall. The demand side is closed at
the level of identities, though not at the level of strength: Goldbach
itself needs only a one-sided bound on $E_3$, at a threshold
$\asymp N$ for almost all $N$ (Proposition [prop:onesided]). We then
close the interior of the design space those two endpoints open
(Theorems [thm:D] and [thm:Dprime],
Propositions [prop:E] and [prop:Dpp]). We also report a defect
in equation (18) of the source paper, and its repair.

The demand side reduces everything to
$C(N)=\sum_{n<N}\Lambda(n)\mu(N-n) = o(N)$.
Sections [sec:wall]–[sec:floor] give the exact facts about
that quantity: its second moment in closed form
(Proposition [prop:V]), an exact identity for the aggregate second
moment (Lemma [lem:MP]), an exact closed form for the fluctuation
of any cell mean (Lemma [lem:cellmom]), and the observation that
this error bar falls like $(\log N)^{-1/2}$ and not like $n_c^{-1/2}$
(Proposition [prop:coh]), which changes the width of every interval
in this subject built from a count. Two controls
(Lemmas [lem:coin] and [lem:placebo]) delimit what can be
measured at all: in particular the excess of $\mathrm{Var}\,C$ over its
random-sign value cannot be estimated by any centred statistic, because
such a statistic returns the same value on a coin.

We also quantify the obstruction to the natural conditional route:
Chowla's conjecture does *not* give $\mathrm{Var}\,C \sim V$,
because the coefficient amplifying $S(h)$ grows like $N/\log N$
(Proposition [prop:W]); and since that coefficient is nonnegative
while the true $S$ already overspends the budget it defines, no bound
on $|S(h)|$ suffices.

The deliverable is the map and the exact facts, not a theorem toward
Goldbach. Net progress toward Goldbach: zero. We state throughout what
is and is not claimed, and \S[sec:notclaimed] states what has been
deliberately withheld from this version.


## Introduction


### The conditional frame


Let $N$ be a large even integer,
$\tilde r(N) = \sum_{n<N}\Lambda(n)\Lambda(N-n)$, and $\SS(N)$ the
singular series. Write $EH_\mu(Q)$ for the assertion that for every
$A>0$

$$
\sum_{q\le Q}\max_{y<N}\max_{(a,q)=1}
  \Bigl|\sum_{\substack{n\le y\\ n\equiv a\,(q)}}\Lambda(n)\mu(N-n)
   -\frac{1}{\varphi(q)}\sum_{n\le y}\Lambda(n)\mu(N-n)\Bigr|
  \ll_A \frac{N}{(\log N)^A}.
$$

Huang and Li [HL], continuing Pan [Pan], proved that for each
fixed $A>0$, $EH(N^{\theta}(\log N)^{2A+8})$ together with
$EH_\mu(N^{1-\theta})$ implies
$\tilde r(N)\ge \SS(N)\bigl(1-\AAA(N)\bigr)N + O(N(\log N)^{-A})$, and
deduced (their Corollary 1) that in view of Bombieri–Vinogradov,
binary Goldbach for all sufficiently large even integers follows from
$EH_\mu(N^{\theta'})$ alone for any single $\theta' > 1/2$.

Here $\AAA(N)$ is their equation (7). We have given it the symbol
$\AAA$, rather than the $A$ of [HL], both because $A$ is already
the free exponent in the display above and because $\AAA(N)$ is not an
incidental constant: for even $N$ it is *identically* the local
factor of Proposition [prop:V], the exact second moment of the
scalar this paper reduces to. That is not a coincidence. Huang–Li
reach it by the triangle inequality
$|C(N)| \le \sum_{n<N}\Lambda(n)\mu^2(N-n)$, whose right-hand side
counts squarefree shifted primes — which is what
Proposition [prop:V] evaluates. So the deficiency $1-\AAA(N)$ in
their lower bound is one minus the wall's own local density, and the
two halves of this paper meet at a constant neither side names.

This is the sharpest conditional frame we are aware of: one sentence
stands between classical technology and Goldbach. This paper asks what
that sentence actually costs.


### The organizing principle: two couplings


The hypothesis contains a product of two Möbius factors, and there
are exactly two ways their arguments can be coupled. This turns out to
organize everything.


- **Divisibility-coupled** (the *demand* side). Inside
  Huang–Li's proof, $EH_\mu$ is consumed only in the fixed residue
  class $n \equiv N \pmod k$, i.e. $k \mid N-n$: one Möbius
  argument divides the other. Section [sec:demand].

- **Difference-coupled** (the *supply* side). Any attempt
  to *prove* $EH_\mu$ by dispersion arrives at the dilate field
  $D(k) = \sum_m \mu(m)\mu(N-mk)$, where the two arguments are related
  by a difference, not a divisibility. Section [sec:supply].


The two behave completely differently, and each is closed for its own
reason:


```latex
\begin{center}
```latex
\begin{tabular}{p{0.20\textwidth}p{0.34\textwidth}p{0.34\textwidth}}
\hline
 & \textbf{Divisibility-coupled} & \textbf{Difference-coupled}\\
\hline
Divisor switch & exact, and free & exact, but useless\\
M\"obius lands on & the \emph{short} variable, $m\le N^{1/2-\delta}$ &
 the \emph{long} variable, $m \asymp N/K \ge N^{2/3}$\\
Consequence & Bombieri--Vinogradov closes it: Theorem~\ref{thm:A} &
 no known machine applies\\
But & the weighted version is \emph{equivalent} to Goldbach
 (Theorem~\ref{thm:C}) & no named coupling surface found, internal or
 external (\S\ref{sec:closures})\\
\hline
\end{tabular}

```

\end{center}
```


### What this paper claims {#sec:claims}


- **Claimed**: Theorem [thm:A], Corollary [cor:B],
  Theorem [thm:C], Theorems [thm:D] and [thm:Dprime],
  Propositions [prop:E], [prop:Dpp], [prop:V],
  [prop:W], [prop:scaleinv] and [prop:coh],
  Lemmas [lem:MP], [lem:cellmom], [lem:coin]
  and [lem:placebo]; the defect in [HL]'s equation (18) and
  its repair; the measurements reported here, which are reproducible.

- **Conjectured**: Conjecture [conj:L].

- **Not claimed**: any theorem toward Goldbach.
  Theorem [thm:A] is unconditional but Goldbach-neutral: it
  removes exactly the half of the demand that carries no Goldbach
  content. See also \S[sec:notclaimed].


## The demand side {#sec:demand}


### Where $EH_\mu$ is actually consumed


For $(k,N)=1$ and $1\le t<N$ put

$$
E_\mu(t;k) = \sum_{\substack{n\le t\\ n\equiv N\,(k)}}
    \Lambda(n)\mu(N-n) - \frac{1}{\varphi(k)}\sum_{n\le t}
    \Lambda(n)\mu(N-n).
$$

With $K = (N-1)/\alpha$, Huang–Li's two consumption sites are

$$
E_3(\alpha)=\!\!\sum_{\substack{k<K\\(k,N)=1}}\!\!\mu(k)\log k\,
  E_\mu(N;k),
  \qquad
  E_4(\alpha)=\!\!\sum_{\substack{k<K\\(k,N)=1}}\!\!\mu(k)\,
  \bigl[\text{$E_\mu$ with the extra weight }\log(N-n)\bigr].
$$

They differ in exactly one respect: the weight $w_k$ is $\log k$ in
$E_3$ and $1$ in $E_4$ (the $\log(N-n)$ inside $E_4$ is
$k$-independent and comes off by partial summation). In both cases
Huang–Li discard the signs $\mu(k)$ on the first line by the triangle
inequality, and only then appeal to $EH_\mu$. Keeping the signs is
what the following results exploit.


### Results imported from the companion note {#sec:imported}

The demand-side results are stated and proved once, in
`paper/theorem_A.md`. They are not restated here; this section lists
what is imported and what each one is used for downstream. References
of the form [thm:A] resolve there.

| Imported | What it gives this paper |
|---|---|
| Theorem [thm:A] | $w_k=1$ is unconditional, so the $E_4$ half of the demand carries no hypothesis |
| Corollary [cor:B] | hence the whole $EH_\mu$ demand is the single scalar $E_3(\alpha)$ |
| Theorem [thm:C] | that scalar is Huang–Li's (22); its asymptotic half needs $C(N)=o(N)$, which is §[sec:wall] onward |
| Theorem [thm:D] | the interior of the weight space is empty, so $E_3$ is not one candidate among a family |
| Theorem [thm:Dprime] | and stays empty granting $EH$, so the missing ingredient is a mechanism and not a level |
| Proposition [prop:E] | the circle method has no margin on $C(N)$ either |
| Proposition [prop:Dpp] | nor does any polynomial weight |

Only one feature of the mechanism is needed below, so it is recorded
here rather than looked up. The class $n\equiv N\pmod k$ is $k\mid N-u$
with $u=N-n$, so the $k$-sum is an incomplete divisor sum; completing
it costs $N^{o(1)}$ and moves the work onto $k\ge K$, where $u=mk$
forces $m<N/K=N^{1-\theta'}\le N^{1/2-\delta}$ and, for squarefree
$u$, $\mu(u)\mu(k)=\mu(m)\mu^2(k)$. **The surviving Möbius sits on the
short variable.** That single assignment is what
Bombieri–Vinogradov consumes, and §[sec:R4] shows it is exactly what
the supply side's type-II cut forbids — which is why the switch
transfers to the supply side as a technique and not as a result.

The two weights differ in one respect only, $w_k=\log k$ against
$w_k=1$, and the difference decides everything: for $w_k=\log k$ the
complete divisor sum is $-\Lambda$ rather than an indicator, because
$\mu*\log=\Lambda$ — the identity Huang–Li *start* from (their
(10)). Closure on the demand side is therefore at the level of
identities, not of estimates, and no choice of $\theta'$, truncation
or smoothing evades it.

The loss exponent $\sqrt{(1/2)\log N}$ of Theorem [thm:D] is the
$\sqrt N$ barrier itself: it is the gap between where a weight must
live to see $C(N)$ and where Bombieri–Vinogradov permits it to live.
That is a no-go for one precisely specified method — divisor
switching with BV as the only input — over that method's entire
weight space; it is not an obstruction to other methods.


### A defect in the published paper


Equation (18) of [HL] replaces the $n$-dependent constraint
$k < (N-n)/\alpha$, present three displays earlier in the derivation of
$S_2(\alpha)$, by the $n$-free $k < (N-1)/\alpha$. The terms thereby
included are exactly those with $d=(N-n)/k \le \alpha$, and they sum to

$$
\Delta = -\sum_{2\le m\le\alpha}\mu(m)\log m
    \sum_{\substack{k<(N-1)/\alpha\\ (k,m)=1}}\mu^2(k)\Lambda(N-mk),
$$

whose trivial bound $\ll N(\log N)^2$ exceeds the target. It closes by
the same machinery — Möbius again on the short variable — under
hypotheses already assumed: unconditionally by Bombieri–Vinogradov in
the Corollary-1 regime, and under the assumed $EH$ in the Theorem-1
regime. [HL]'s Theorem 1 and their Corollary 1 stand as stated.


## The supply side {#sec:supply}


### The consumable


The chain's supply-side consumable, distilled, is

$$
\textbf{E1 (weak form).}\quad
  D(k) = \!\!\sum_{\sqrt N < m \le N/k}\!\! \mu(m)\mu(N-mk),
  \qquad
  \sum_{k\sim K}|D(k)|^2 \ll (\log N)^{-2A-2}\sum_{k\sim K} M_k^2,
$$

for $K \le N^{1/3}$ dyadic, where $M_k$ is the length of the $m$-range.
A fixed log-power saving is the currency; $o(1)$ and $\log\log$ savings
are not consumable.

The normalisation is $M_k^2$, i.e. the *trivial* bound, not
$M_k$, the square-root scale: the wall is
$T_{II} = \sum_{k\sim K}b_kD(k)$ with $|b_k|\ll\log N$, needed at
$N(\log N)^{-A}$, and Cauchy–Schwarz in $k$ gives
$\sum_k|D(k)|^2 \ll N^2/(K(\log N)^{2A+2})$ with
$\sum_{k\sim K}M_k^2 \asymp N^2/K$. The distinction is load-bearing:
at the square-root normalisation the demand would be
$\sum|D|^2/\sum M_k \ll 8.8\cdot10^{-6}$ at $N=10^8$, $A=1$, which our
own measurements refute directly — the measured ratio is
$0.34,\,0.39,\,0.32$ over three dyadic bands, which is the density of
the surviving support, i.e. of the $m$ for which $\mu(m)$ and
$\mu(N-mk)$ are both nonzero. Normalised by that support, measured
band by band, the ratio is $1.04,\,1.17,\,0.96$: the field exhibits
*exact* square-root cancellation on the terms that survive.


#### Remark (this row is not reproduced) {#rem:e1row}
<!-- evidence: audit_support_density.py -->

The two rows just quoted are the audit's rule Z6 and it fails. On the
bands that Remark [rem:supp] pins — the ones under which the support
density reproduces to four decimals —
$\sum_k|D(k)|^2/\sum_k M_k$ measures
$0.4795,\,0.3393,\,0.3480$ against the printed
$0.34,\,0.39,\,0.32$, and normalised by the support
$1.4518,\,1.0286,\,1.0482$ against the printed
$1.04,\,1.17,\,0.96$. The gap is not a convention: the row was
recomputed under all nine band offsets from $-4$ to $+4$ and all four
aggregations (ratio of sums, mean of per-modulus ratios, each with and
without the dead moduli). Over those offsets band $1$ ranges in
$[0.4199,\,0.4841]$, band $2$ in $[0.3230,\,0.3714]$ and band $3$ in
$[0.3441,\,0.3524]$: the printed band-2 value $0.39$ lies above its
whole range and the printed band-1 value $0.34$ below its whole range.
Band $1$ is the discrepancy that matters — $0.48$ against $0.34$ is
not square-root cancellation on the surviving support but half again
as much.

What survives is weaker and is what the section actually needs: the
ratio is of order the support density and not of order
$(\log N)^{-2A-2}\approx 8.8\cdot10^{-6}$, so the square-root
normalisation is refuted by four orders of magnitude either way. The
claim of *exact* square-root cancellation is withdrawn pending a
statement of what was computed.


#### Remark (the support density is not $\prod_p(1-2/p^2)$) {#rem:supp}
<!-- evidence: audit_support_density.py -->

The three bands are $K<k\le 2K$ for $K=58,116,232$, i.e.
$[59,116],\,[117,232],\,[233,464]$ with
$464=\lfloor N^{1/3}\rfloor$. The convention has to be stated: the
other reading of "dyadic", $K\le k<2K$, holds the same $58$, $116$,
$232$ values of $k$ and so is not excluded by the counts, but it gives
$0.3345,\,0.3278,\,0.3310$ for the aggregate below instead of
$0.3303,\,0.3298,\,0.3320$ — a shift larger than the effect the
paragraph is about. That was the audit's rule Z5, which fails under
the $K\le k<2K$ reading and holds under this one; the convention is
recorded here because the counts do not carry it.

It is tempting to read that density as $\prod_p(1-2/p^2)=0.3226341$,
the value obtained by assuming $p^2\mid m$ and $p^2\mid N-mk$ cut out
two distinct classes mod $p^2$. At $N=10^8=2^8\cdot5^8$ that assumption
fails at $p=2$ and $p=5$, where $p^2\mid N$ collapses the second
condition onto the first, and the true density depends on $v_p(k)$:
measured exactly over all $k$ in the three bands it ranges from $0$ to
$0.594330$, and lands within $0.01$ of $0.32263$ for only $1$, $1$ and
$3$ of the $58$, $116$ and $232$ values of $k$. In particular $D(k)$
vanishes identically whenever $4\mid k$ or $25\mid k$ — $2799$ of the
$9999$ values $k<\sqrt N$, or $28.0\%$, and that set is exactly the
vanishing set, with no other $k$ in the bands giving $D(k)=0$.
Aggregated the coincidence is close but not exact:
$\sum_k\mathrm{supp}(k)/\sum_k M_k$ is $0.3303,\,0.3298,\,0.3320$ on
the three bands, high by $2.2$ to $2.9\%$. The figures above use the
measured support throughout; the closed form is quoted here only to
say that it is not the right one.


At the
correct normalisation the target is cleared with margin
$(N/K)(\log N)^{-2A-2}$, which is asymptotic and is below $1$
everywhere accessible, so no computation is offered as evidence for it.


#### Remark (the mandatory comparison) {#rem:scale}

A target must be checked against one's own measurement of the same
quantity before it is used to adjudicate anything. Writing a target at
the square-root scale when the chain consumes it at the trivial scale
inflates every demand by a power of $N$, and no verdict passed on a
budget basis survives it. We record this as a standing rule because
the same slip occurred twice in the course of this work, and both times
it was caught only by that comparison.


### Conjecture L


#### Conjecture (the factorization law) {#conj:L}
<!-- evidence: audit_support_density.py -->

Each Möbius family probed in this program — prime-indexed dilate
pairs $C_{k,k'} = \sum_{p\sim P}\mu(N-pk)\mu(N-pk')$, integer-indexed
dilates and their pairs, and Möbius sums over thin progressions —
factorizes as

$$
\text{field} \;=\; \mathbf M \times \mathbf G,
$$

where $\mathbf M$ is the deterministic local mask, computed exactly by
finite modular enumeration from the $v_q$-data of $(N,k,k')$, and
$\mathbf G$ is, on the surviving support, a fluctuation with no mean
field and no class structure.


The evidence is of two kinds and we separate them. The E1 arm is stated
above and has been reproduced independently, and so, now, is the blind
mask prediction — see [rem:maskstamp], which also records what
re-verifying it cost and what its grid could not show. The exact
$(v_2,v_3)$ cells were attempted next and the attempt did not reach
them — [rem:cellstamp] records what it would take. The remaining two
stamps — pair statistics and the Wishart pair spectrum — are recorded
in the repository and have *not* been independently re-verified; they are cited here as the reason the
conjecture is stated, not as measurements this paper vouches for.

#### Remark (one of the four stamps re-verified, and what re-verifying it showed) {#rem:maskstamp}
<!-- evidence: verify/pass2/code/verify_conjL_mask_zeros.py -->
<!-- evidence: verify/pass2/code/verify_conjL_gap_resolve.py -->
<!-- evidence: verify/pass2/code/verify_conjL_mask_branch.py -->

Of the four stamps the paragraph above sets aside, the blind mask
prediction was the one worth taking first, because it is not a
statistic. It asserts a set equality — the field's support is empty on
exactly the pairs the mask names — so a second implementation either
lands on the same set or does not, and the failure mode that has
limited nearly every other measurement in this program, *too noisy to
tell*, does not exist for it. `verify/pass2` re-derived the rule from
this conjecture's own statement, without reading the original code.

**The rule holds.** On the grid re-derived here the mask names $116000$
pairs of $401000$, and the field rebuilt independently has empty
support on exactly those: zero predicted-but-nonzero, zero unpredicted
zeros. What `M.3 EXACT` was asserting is true, and it is now true on a
second implementation.

**The number cannot be checked.** The stamp reports $115950$, and the
re-derivation got $116000$. The gap is not a disagreement about the
rule; it is a disagreement about the grid. The stamp publishes exactly
two facts about its own grid — $401$ values of $N$, $1000$ values of
$k$, and a note that $3\cdot10^6$ and the step $2500$ are both
divisible by $4$ — and those two do not fix the offsets. Scanning the
$4160$ grids consistent with them gives counts from $115950$ to
$116101$, and **$100$ of them give exactly $115950$**. So the stamp is
right and unreconstructible: a reader who does what this pass did
cannot land on its number except by luck.

**And its grid cannot see part of the rule it stamps.** Because $p$
runs over odd primes there is a $2$-adic branch: if $v_2(N)=v_2(k)=a\ge
1$ then $N/2^a - p\,k/2^a$ is even for every odd $p$, so
$v_2(N-pk)\ge a+1\ge 2$ and the support is empty. At $a=1$ that is not
of the form $q^2\mid\gcd(N,k)$ for any $q$. Every $N$ on the stamp's
grid is divisible by $4$, which forces $a\ge2$ and folds the branch
into $\min(v_q(N),v_q(k))\ge 2$ — the two rules agree on all $401000$
pairs there, so no stamp on that grid can tell whether a mask contains
the branch. On a grid with $2\,\|\,N$ they separate sharply:

| | annihilated | field disagrees |
|---|---|---|
| with the branch | $25.94\%$ | $0$ of $101000$ |
| without it | $1.26\%$ | $24930$ of $101000$ |

and the branch alone carries $0.2468$ of the grid — the same quarter
the stamp's own note records one valuation higher.

This does not say the original mask omits the branch; the pass was
blind and did not read that code. It says the single witness is a
witness to less than the whole rule. Two things follow for the other
three stamps, and they pull in opposite directions: the substance of
this one survived independent re-derivation, which is evidence for the
conjecture; and the stamp's *form* — a number without the grid that
produced it — is what made re-verification cost three scripts instead
of one.


#### Remark (the second stamp, and the second parameter it does not print) {#rem:cellstamp}
<!-- evidence: verify/pass3/code/verify_conjL_exact_cells.py -->
<!-- evidence: verify/pass3/code/verify_conjL_cells_range.py -->

[rem:maskstamp] took the blind mask prediction and found it right and
unreconstructible. `verify/pass3` took the next of the four, the exact
cells, whose published form is a single interval: every viable
$(v_2,v_3)$ cell between $0.99$ and $1.06$. The quantity is not in
doubt — it is the variance ratio $E[D^2]/\mathrm{supp}$, which is $1$
under exact square-root cancellation on the surviving terms, so a cell
away from $1$ is the class structure [conj:L] forbids.

**What is in doubt is the cell.** $(v_2,v_3)$ does not say of what.
Three readings are natural — the valuations of $k$, of $\gcd(N,k)$, of
$N$ — and on the field rebuilt here none of them reproduces the
interval:

| reading | viable cells | outside $[0.99,1.06]$ | mean |
|---|---|---|---|
| $k$ | $13$ | $8$ | $0.9844$ |
| $\gcd(N,k)$ | $4$ | $3$ | $0.9750$ |
| $N$ | $31$ | $20$ | $0.9674$ |

**And a second parameter is missing.** The pooled ratio — which no
cell choice can affect — runs from $0.99701$ at $p\le250$ to $0.98179$
at $p\le2900$, a drift of $0.01522$, a fifth of the stamp's own
interval. So even a correct reading of the cell leaves the number
undetermined until the prime range is named, and it is not named. Over
the twenty-one configurations tried, three readings by seven ranges,
the best is $9$ viable cells of $13$ inside the interval, at the
smallest range; never all.

This is not a claim that the stamp is wrong, and the pre-registration
fixed that reading before the run: three readings is a guess at a
space whose size the stamp does not publish, and failing to find a
measurement is not finding it false. Two of this pass's own
predictions broke and are recorded broken — the band's asymmetry about
$1$ is not a property of the quantity (every reading's mean is
*below* $1$), and the drift is not monotone, one step of $+0.00016$
reversing it. A third, that the repository's pooled $0.99781$ would
lie inside the swept range, missed by $0.0008$ at the top: the
direction points below $p\le250$ and the scan did not go there.

**What the two passes together establish is about form, not truth.**
Both stamps examined are right as far as anything here can tell, and
neither can be reconstructed from what it prints — one missing its
grid, this one missing its cell index and its prime range. That is a
property of how the stamps were published, and it is the reason
[conj:L]'s evidence paragraph separates what has been re-verified from
what has been cited. The remaining two stamps are still cited only.


#### Remark (the same question turned on this repository, and what paying it costs) {#rem:fieldpinned}
<!-- evidence: audit_field_pinned.py -->

[rem:maskstamp] and [rem:cellstamp] found the same defect in two of
[conj:L]'s stamps: the number is right as far as anything here could
tell, and cannot be rebuilt from what was printed. Recording that
about another tree and not asking it here would be incoherent, so it
was asked. The test is one-sided and was stated that way before the
run: a script fixes its field in module-level constants, and if a
constant's value appears nowhere in that script's own result file, a
reader holding only the file cannot rebuild the run. A match may be
coincidence; only an absence carries information.

Of $154$ script–result pairs, $135$ are clean and $19$ are not — so
the discipline mostly holds (F1, F2). **But F3 was pre-registered as
"no absence is a field bound", and it is REFUTED.** As first measured
it failed fifteen times, seven of them the seed of a null:

| what is missing | at first | now |
|---|---|---|
| `SEED` | $7$ | $2$ |
| `UMAX = 400` | $4$ | $4$ |
| `WIN = 0.95` | $3$ | $3$ |
| `SWEEP = 0.1` | $1$ | $1$ |

The seed is the one that matters. In all seven it drives
`np.random.default_rng(SEED)` — the null, the bootstrap, the
permutation — and a null whose seed is unrecorded is not reproducible
from the file that reports it. That is the same thing
[rem:cellstamp] recorded about another tree's stamp, and this
repository's own G4, which requires every result file to open with
`STATISTIC:` and `FIELD:`, had already decided that a result file must
stand alone. G4 only checks that the line is present, so **G78 now
checks that the constants reach the file**, with the remaining
absences on a list that may only shrink.

**Five of the seven seeds were closed and two were not, and the reason
is worth more than the fix.** Adding the seed to a header changes the
script, so G18 requires the script be re-run, and G22 then requires
every consumer of its result to be re-run in turn. For five of them
that closure is minutes. `audit_residue_arithmetic` and
`lab_primorial_ladder` sit upstream of the rung ladder —
`lab_primorial_ladder` → `audit_primorial_reach` →
`audit_primorial_rung10` → rungs $11$, $15$, $16$ and the uniform-cap
sweep — so a one-line documentation fix on either of them costs a
recomputation of the ladder at $N$ up to $1.97\cdot10^9$, hours of it,
for numbers that would come back identical. **The debt is not
unpaid because nobody noticed; it is unpaid because the gate prices a
documentation change at the cost of the computation beneath it.**
That is the correct price for a *content* change and the wrong one
here, and the two files stay on the list carrying that reason.

Every result file re-run in closing the five came back byte-identical
apart from the header line — which is a second, unlooked-for check
that the recorded seeds do determine their nulls exactly.

Three limits, all real. The census reads constants written as
`NAME = value` and misses tuple assignments like
`NN, NK, PMAX = 401, 1000, 2000`, so the true count is higher than
$19$ and the clean share is optimistic. The bound family is a list of
name fragments fixed before the run, so it can miss a bound named
otherwise and can flag a constant that is cosmetic in its own script —
which is why every absence is printed by file and name, to be read
rather than counted. And G18 and G22 compare modification times, which
git does not carry: **the property that this repository's gate passes
is not reproducible from the repository**, only from a working copy
that has run everything in order. That is the same class of defect as
the two stamps, one level up, and it is not fixed here.

Conjecture [conj:L] is *stronger* than the chain needs: E1
requires only square-root cancellation on average. What is missing is
therefore not knowledge of the structure — the mask is an algorithm
and $\mathbf G$ is featureless in every measurement — but a proof
technique for the amplitude of a featureless object.


## The wall, exactly {#sec:wall}


The chain of Section [sec:demand] reduces the whole problem to one
scalar,

$$
C(N) \;=\; \sum_{n<N}\Lambda(n)\mu(N-n) \;=\; o(N).
$$

This section gives that scalar its exact second moment, an exact
identity for its aggregate second moment, and the two controls that
delimit what can be measured about it at all.


### The exact scale
#### Proposition (The exact scale) {#prop:V}
<!-- evidence: lab_second_moment.py -->

Let

$$
V(N) \;=\; \sum_{v<N}\mu^2(v)\,\Lambda(N-v)^2,
  \qquad
  W(N) \;=\; \sum_{w<N}\Lambda(w)^2,
$$

and let

$$
\AAA(N) \;=\; \prod_{q\,\nmid\,N}\Bigl(1-\frac{1}{q(q-1)}\Bigr).
$$

Then $V(N) = W(N)\,\AAA(N)\,(1+o(1))$, and consequently
$V(N) \sim \AAA(N)\,N\log N$.


**Proof (Derivation).** 
$\Lambda(w)^2$ is supported on prime powers with weight $(\log p)^2$,
so with $w=N-v$,

$$
V(N) \;=\; \sum_{p^k<N}(\log p)^2\mu^2(N-p^k)
        \;=\; \sum_{p<N}(\log p)^2\mu^2(N-p) + O(\sqrt N\log^2 N),
$$

a $(\log p)^2$-weighted count of *squarefree shifted primes*. Its
local density at a prime $q$ is $1$ when $q\mid N$ — there
$q^2\mid N-p$ forces $q\mid p$, hence $p=q$, a single term — and
$1-1/\varphi(q^2) = 1-1/(q(q-1))$ when $q\nmid N$, since the bad class
$p\equiv N \pmod{q^2}$ is then a unit class. The statement is Mirsky's
theorem on squarefree values of shifted primes [Mir49] together
with partial summation; it is recalled here, not claimed.
 ∎


**The local factor is $\AAA$, not $\SS$..**
The singular series $\SS(N)$ is Hardy–Littlewood's and is correct for
the Goldbach *count*; the right object for the second moment of
$C(N)$ is $\AAA(N)$. Both depend on $N$ only through the set of primes
dividing $N$ — though in opposite ways, $\SS$ through a product
*over* those primes and $\AAA$ through a product over the primes
*not* among them — which is why the substitution is easy to make
and hard to see.
The statistic, stated. For a candidate $c\in\{\AAA,\SS\}$ put
$z_c(N)=\bigl(V(N)/W(N)\bigr)/c(N)$; the figure quoted is
$\mathrm{sd}\bigl(z_c/\overline{z_c}\bigr)$, i.e. the candidate is
rescaled to the measured mean so that only its shape in $N$ is judged.
Over even $N$ in $[10^5,\,1.6\cdot10^7]$ it is $0.000323$ for $\AAA$
against $0.245235$ for $\SS$ — a factor of $759.3$. The reading has to
be stated: the other natural one, subtracting a least-mean multiple
rather than dividing, gives $0.000317$ and $0.299040$ instead, and the
two disagree in the second figure by more than a fifth.

Dividing through by $W(N)$ removes the analytic factor exactly, so no
asymptotic for $\sum(\log p)^2$ is assumed: the ratio $V(N)/W(N)$
against $\AAA(N)$ has mean $1.0000002$ with standard deviation
$0.0001659$ in the top octave. At $N=4\cdot10^6$ the ratio $W/V$ is
$1.270800$ against the predicted $1/\AAA(N) = 1.270204$.

By cells, with cell $j$ the even $N$ divisible by $2$ and by the first
$j-1$ of $3,5,7,11,13$:

$$
\begin{array}{r|cccccc}
 \text{modulus} & 2 & 6 & 30 & 210 & 2310 & 30030\\
 \text{count} & 7950000 & 2650000 & 530000 & 75714 & 6883 & 529\\\hline
 \overline{V/(W\AAA)} & 1.000000475 & 1.000000147 & 1.000000010
   & 0.999999570 & 0.999999037 & 0.999999646
\end{array}
$$

Version 3 said this "agrees to six decimals in each of the six radical
cells". Five of the six do; the $2310$ cell reads $0.999999$, missing
by $1$ in the sixth place on $6883$ values of $N$. The audit
pre-registered that claim as its rule X5 and it fails there.

The lower cutoff $10^5$ is not cosmetic: the figure for $\AAA$ rises as
the cutoff falls, reading $0.000346,\,0.000398,\,0.000470,\,0.000529$
at cutoffs $5\cdot10^4,\,10^4,\,10^3,\,10^2$, while the figure for
$\SS$ is unaffected because it is dominated by $\SS$'s wrong shape
rather than by noise. Any comparison of the two must state its range.
Version 3 put the degraded figure at $0.000582$; under the reading
above no lower cutoff reaches it — the sweep saturates near
$0.000555$ as the cutoff approaches the bottom of the range — so that
one number is withdrawn rather than corrected. It was the audit's rule
X6.


#### Remark (the second-order form) {#rem:secondorder}
<!-- evidence: lab_second_moment.py -->

Proposition [prop:V]'s $V\sim\AAA N\log N$ is correct but converges
slowly, and the audit's rule X7 — a $5\%$ band on
$V/(\AAA N\log N)$ at $N=1.6\cdot10^7$ — fails on it: the measured
ratio is $0.927137,\,0.933409,\,0.939701$ at
$N=10^6,\,4\cdot10^6,\,1.6\cdot10^7$. The band was set as an effect
size and not from the asymptotic, which is the same error as
Remark [rem:cap] in the companion note. Partial summation from
$\theta(x)\sim x$ gives $\sum_{p\le x}(\log p)^2 \sim x\log x-x$, so the
second-order form is $\AAA(N)\,(N\log N-N)$, a factor $1-1/\log N$
smaller — $0.939716$ at the top, against the measured $0.939701$.
Against that form the ratio reads
$0.999482,\,0.999134,\,0.999984$ at the same three $N$.


### The aggregate second moment, exactly


#### Lemma (the second-moment identity) {#lem:MP}
<!-- evidence: analytic -->

Fix $X$ and let

$$
\widehat C(N) \;=\; \sum_{\substack{n+v=N\\ n\le X,\ v\le X}}
    \Lambda(n)\,\mu(v)
$$

be the truncated convolution, whose support runs to $2X$. Then

$$
\sum_{N\le 2X} \widehat C(N)^2 \;=\; \sum_{|h|<X} M(h)\,P(h),
  \qquad
  M(h) = \sum_{v,\,v+h\le X} \mu(v)\mu(v+h),
  \quad
  P(h) = \sum_{w,\,w+h\le X} \Lambda(w)\Lambda(w+h).
$$


**Proof.** 
Expanding the square as a double sum over $(n,v)$ and $(n',v')$ with
$n+v=n'+v'$ and summing over *all* $N$ leaves the single
constraint $n-n' = v'-v =: h$, with $(n,n')$ and $(v,v')$ otherwise
ranging independently over $[1,X]^2$. Collecting the $\Lambda$-pairs
and the $\mu$-pairs separately gives the two factors.
 ∎


#### Remark

The truncation on the left is necessary and it is easy to lose. If one
writes instead $\sum_{N\le X}C(N)^2$ with $C$ the untruncated
convolution, the three surviving indices $(n,v,h)$ are constrained by
$n+v\le X$ — a *simplex* — while $M(h)P(h)$ ranges them over
the *box* $[1,X]^2$ and so also counts every pair with
$X< n+v\le 2X$. The two sides then differ by a factor near $1.57$ at
$X$ from $800$ to $3200$, and the ratio does not tend to $1$. The form
above is exact to machine precision at every $X$ tested, to
$X=1.6\cdot10^7$.


The identity is exact and unconditional, and it makes the aggregate
size of the wall a statement about two shifted correlations, one of
which ($P$) is the Hardy–Littlewood prime-pair sum and the other of
which ($M$) is the binary Chowla sum.


### What the excess is, and why Chowla does not control it


Write $\rho(N) = \mathrm{Var}\,C(N)/V(N)$, which is $1$ if the $\mu(v)$
on the surviving support are replaced by independent signs. Pointwise,
$C(N)^2 = V(N) + \mathrm{OffDiag}(N)$ with

$$
\mathrm{OffDiag}(N)
  = \sum_{h\ne0}\ \sum_{\substack{p,\,p+h\ \mathrm{prime}}}
    (\log p)(\log(p+h))\,\mu(N-p)\mu(N-p-h),
$$

which is algebra. Passing from there to $\sum_{h\ne0}c(h)S(h)$ with
$c(h)=\sum_{p'-p=h}(\log p)(\log p')$ and
$S(h)=\langle\mu(u)\mu(u-h)\rangle$ is *not* algebra: it replaces
the values $\mu(N-p)\mu(N-p-h)$, sampled on the sparse set of $p$ with
$p$ and $p+h$ both prime, by an average over all $u$. We do not assume
that decoupling. What we do record is that even granting it, the
conditional route does not close.


#### Proposition (the amplification) {#prop:W}
<!-- evidence: audit_amplification.py -->

Let $\Gamma(N) = \bigl(\sum_{h\ne0}c(h)\bigr)/V(N)$. Then
$\Gamma(N) \sim N/(\AAA(N)\log N)$, and consequently a hypothesis
$|S(h)|\le\varepsilon$ yields nothing better than
$|\rho-1| \le \varepsilon\,\Gamma(N)$. Because $c(h)\ge0$ that
inequality is sharp, so such a hypothesis can give $\rho\to1$ only for
$\varepsilon = o(1/\Gamma(N)) = o(\log N/N)$; and no bound of that
strength is true, because the measured $S$ already exceeds the budget
— display [eq:budget] below. Hence no bound on $|S(h)|$ alone
can give $\rho\to1$.


Measured with $c$ as defined above — over primes, not prime powers —
$\Gamma = 1.5128\cdot10^{3},\ 1.8412\cdot10^{4},\ 3.5759\cdot10^{5}$ at
$N=10^4,\,1.6\cdot10^5,\,4\cdot10^6$, with
$\Gamma\log N/N = 1.3933,\,1.3789,\,1.3590$ against the predicted
$1/\AAA(N) = 1.270204$. The numerator needs no correlation: summing
$c(h)$ over all $h\ne0$ is exactly $\theta(N)^2-\sum_{p<N}(\log p)^2$.


#### Remark (one symbol, two meanings, in one paragraph) {#rem:cdef}
<!-- evidence: audit_amplification.py -->
<!-- symbol: c = sum_{p'-p=h} (log p)(log p'), over primes only -->

Version 3 printed the row above as
$1.5489\cdot10^{3},\ 1.8517\cdot10^{4},\ 3.5798\cdot10^{5}$. Those are
the audit's rules Y1 and Y2 and they fail against the definition of $c$
given here, by $-2.33\%,\,-0.57\%,\,-0.11\%$ — a gap shrinking like
$N^{-1/2}$, which is the signature of prime powers, since
$\psi(x)-\theta(x)\asymp\sqrt x$ and the numerator is squared. Computed
instead with $c(h)=\sum_n\Lambda(n)\Lambda(n+h)$ the row reproduces to
every printed digit. The budget $\mathcal B(X)$ below, however,
reproduces to every printed digit with $c$ over *primes*: so the
two displays of this one paragraph were computed with two different
$c$, and the text defines only one of them.

Nothing in Proposition [prop:W] turns on it — the two differ by
$1+O(N^{-1/2})$ and the claim is the asymptotic $N/(\AAA\log N)$, which
holds for either. The declaration above fixes the meaning to the one
the text states, and the row is restated in it.

Two consequences. First, Chowla's conjecture asserts $S(h)=o(1)$ for
fixed $h$; the strength the display needs is $S(h) = o(\log N/N)$, so
the gap is a factor $N/\log N$.

Second — and this is what closes the absolute-value route rather than
merely costing it a power of $\log$ — the true $S$ overspends that
budget already. Nonnegativity of $c$ says the triangle inequality is
sharp; the measurement says it is also lost. Define the *absolute
budget*, which is exactly what the triangle inequality spends,

$$
\begin{equation}\label{eq:budget}
  \mathcal B(X) = \frac{1}{V(X)}\sum_{0<|h|<X} c(h)\,\bigl|S(h)\bigr|,
  \qquad
  S(h) = \frac{1}{X-|h|}\!\!\sum_{|h|<u\le X}\!\!\mu(u)\,\mu(u-|h|),
\end{equation}
$$

with $c(h) = \sum_{p'-p=h,\ p,p'\le X}(\log p)(\log p')$ as above.
Measured, $\mathcal B(X) = 13.3,\,17.3,\,23.1,\,30.5$ at
$X = 2\cdot10^4,\ 4\cdot10^4,\ 8\cdot10^4,\ 1.6\cdot10^5$: above $1$ by
an order of magnitude, and growing. (Normalising $S$ by $X$ rather than
by the number $X-|h|$ of terms gives $7.9,\,10.3,\,13.8,\,18.2$ —
smaller, same conclusion. Both are reported because the phrase "binary
Chowla sum" does not fix the normalisation, and a budget is not
interpretable until it does.)

So the only route to $\rho\to1$ is *signed cancellation across
$h$*, which is a different object from Chowla's smallness and from the
averaged absolute bound of [MRT15]. Naming what supplies that
cancellation is an open question and is stated as such.


### Two controls, and what they invalidate


#### Remark (the withdrawn law, re-measured: one statistic survives the coin) {#rem:cnlaw}
<!-- evidence: audit_cn_law.py -->
<!-- evidence: audit_cn_coin_spread.py -->

OPEN.md's wall item 1 is the first on that list and had never moved.
The draft's $C(N)=m(N)+\sqrt{V(N)}\,G(N)$ was withdrawn for reasons
about evidence rather than truth — the bulk and tail did not reproduce
under the cell index the text specified, the phase content was
reproduced by the coin and fell to [lem:coin], the mask's significance
was overstated — and the closing condition written there was one
thing: define those statistics and measure them again.

They can now be defined without fitting anything. [prop:V] gives
$V(N)=\sum_{v<N}\mu^2(v)\Lambda(N-v)^2$ exactly, so
$G(N)=C(N)/\sqrt{V(N)}$ is a computed quantity, not an estimated one.
Both $C$ and $V$ are linear convolutions, so one pair of FFTs gives
every $N$ in a band at once; over even $N$ in $(2\cdot10^6,
4\cdot10^6]$ that is $10^6$ values, and the convolution reproduces a
direct sum at eight sampled $N$ to $10^{-11}$ (P1).

**Two of the four predictions broke, and the reason they broke is the
finding.** $\mathrm{sd}(G)=0.92953$, outside the registered
$[0.95,1.05]$ (P2), and a coin arm — $\mu$'s support kept, its signs
replaced, so $V$ is unchanged by construction — missed the real arm by
$8.81$ blocked standard errors in $\mathrm{sd}$ and $15.79$ in excess
kurtosis (P4). Both readings assumed a single coin draw was an error
bar. It is not: one sign pattern is fixed and then used at every $N$,
so all $10^6$ values share it and blocking over $N$ cannot see an
offset common to every block. Sixty-four independent draws give the
error bar the blocking could not, and the two halves of the result
part company.

| | real | coin, over $64$ draws | $z$ |
|---|---|---|---|
| $\mathrm{sd}(G)$ | $0.92953$ | $0.93420 \pm 0.04596$ | $-0.10$ |
| excess kurtosis | $+0.26422$ | $-0.01825 \pm 0.03957$ | $+7.1$ |

The draw-to-draw spread of $\mathrm{sd}$ is $55.4$ times the blocked
standard error (Q2), so **P2's refutation is about the estimator and
not about $\mu$**: $V$ is the second moment *in expectation over
signs*, and no single pattern — $\mu$ among them — has mean square
exactly $V$ across a band. $0.92953$ is where a sign pattern lands
(Q3). Both predictions stay refuted as registered; what they refute is
the reading.

**The kurtosis does not go that way.** No draw among the $64$ comes
near $+0.264$; the largest is $+0.042$ (Q4). $G$ is not Gaussian, and
its non-Gaussianity is not something a sign pattern on $\mu$'s support
produces. That is one statistic of the withdrawn law that the coin
cannot make — the first thing in this branch to survive [lem:coin],
and it is a property of the object binary Goldbach's demand side
reduces to.

The tail is printed beside it and is not symmetric: at the $0.1\%$ and
$99.9\%$ quantiles the real arm reads $-3.1965$ and $+2.8074$ against
a normal's $\mp3.0902$ and a coin draw's $-2.9631$ and $+2.6990$ —
heavier on the left than the normal, lighter on the right. That
asymmetry was not pre-registered and is recorded as an observation.

Three limits. This is one band; nothing here measures how any of these
quantities drifts with $N$, and a law needs that drift. Sixty-four
draws make $z=7.1$ a bound rather than a measurement. And a single
non-Gaussian moment is not the withdrawn law — it is the first term of
it that has evidence, which is exactly what item 1 asked for and less
than what the draft claimed.

#### Remark (the non-Gaussianity is real at every scale reached, and it decays slower than independence) {#rem:cnkurt}
<!-- evidence: audit_cn_kurt_drift.py -->

[rem:cnlaw] found the first quantity in this branch to survive
[lem:coin], and said in the same breath what it could not do: one band
is not a law. Seven octaves, $(2^b,2^{b+1}]$ for $b=17\ldots23$, on a
single sieve and a single pair of FFTs, with a fresh ensemble of $32$
sign patterns recomputed over the whole ladder, are that measurement.
The gate is the overlap: re-measured at a different sieve top,
[rem:cnlaw]'s band returns $+0.26422$ exactly (D1).

| $b$ | even $N$ | excess kurtosis | $z$ against the coin |
|---|---|---|---|
| $17$ | $65\,536$ | $+1.57158$ | $33.4$ |
| $18$ | $131\,072$ | $+1.06049$ | $39.8$ |
| $19$ | $262\,144$ | $+0.69806$ | $20.7$ |
| $20$ | $524\,288$ | $+0.41549$ | $16.5$ |
| $21$ | $1\,048\,576$ | $+0.25305$ | $9.2$ |
| $22$ | $2\,097\,152$ | $+0.14257$ | $6.2$ |
| $23$ | $4\,194\,304$ | $+0.07411$ | $5.0$ |

**The separation is not local to one band** (D2): every octave sits
outside its coin ensemble, the smallest margin being $5.0$ standard
deviations at the top. **And it decays** (D3): the slope against
$\log_2$ of the band midpoint is $-0.24190\pm0.03457$, $t=-7.00$, and
in log-log the seven points lie on a line — $\text{kurtosis}\sim
N^{-0.7312}$ with a standard error of $0.019$ on the fitted slope. So
this is a finite-$N$ effect and not yet a law, which is the weaker of
the two readings the single band could not separate.

**WITHDRAWN.** This paragraph read that the decay is slower than an
independent-term model's $N^{-1}$, taking that rate from the excess
kurtosis of $\sum w_v\varepsilon_v$ as a random variable over the signs.
What is measured here is the excess kurtosis of the empirical
distribution over $N$ of one fixed sign pattern, which is a different
object, and the computed independent-sign null — the coin ensemble in
the table above — does not follow $N^{-1}$ at all. [rem:cnkurtlimit]
carries the correction and what replaces it: the real arm decays
toward the coin arm's level and reaches it at $b\approx25$.

D4 broke and its breaking matters. The coin's own mean excess
kurtosis is not zero but $-0.02103$ to $-0.01376$ across the bands,
$4.43$ of its own standard errors from zero at the top. The $z$
column above already scores against that mean so the separation is
unaffected, but the fitted exponent is not: refitting the power law
on real minus the control's mean gives $N^{-0.6953}$ against the raw
$N^{-0.7312}$. Both are reported and the pre-registered fit is the raw one.
[rem:cnkurtlimit] measures what that control does over two more
octaves: it drifts toward zero rather than holding, so the
baseline these bands are scored against is not a constant.

Two limits, both taken up by [rem:cnkurtlimit]. Seven octaves cannot
tell decay to zero from decay to a positive limit, and that limit is
the question this one becomes. Thirty-two draws make the large $z$
values bounds rather than measurements. A third limit written here —
that no ensemble of independent-term fields had been built — was
wrong when written: the coin ensemble in the table above is exactly
that, which is why the paragraph it justified is withdrawn.

#### Remark (the forecast held, the floor did not, and this branch has a noise floor) {#rem:cnkurtlimit}
<!-- evidence: audit_cn_kurt_limit.py -->

[rem:cnkurt] left one question: seven octaves cannot tell decay to
zero from decay to a positive limit, and a positive limit is what
would make the non-Gaussianity of $C(N)$ a property of $C$ rather than
of the range. Fitting $A N^{-a}+L$ on seven points is a question about
functional form, and this repository has failed at those four times
over. So the limit was not approached by fitting. The floorless power
law, fitted on $b=17\ldots23$ alone, was made to forecast two octaves
it had never seen, and then they were measured — $N$ to
$6.71\cdot10^7$, the shared band $b=23$ reproducing $0.07411$ exactly
as the gate (E1).

| $b$ | forecast | measured | $z$ |
|---|---|---|---|
| $24$ | $0.03913$–$0.06627$ | $0.04847$ | $-0.37$ |
| $25$ | $0.02302$–$0.04088$ | $0.02605$ | $-1.14$ |

**Both landed** (E2), and both slightly low — the decay is not
flattening, which is what approaching a positive floor would look
like from below. No floor is needed to predict where the kurtosis
goes.

**E3 broke, and the pre-registration had already fixed how to read
it.** Fitting $A N^{-a}+L$ on all nine octaves gives
$L=-0.06114\pm0.01015$, $t=-6.03$ — resolved, and *negative*. The
refutation rule written before the run says a resolved negative $L$
means the three-parameter family is fitting shape and not a limit,
since a quantity positive at every octave measured cannot have a
negative limit. The fitted $A=+1870.7$ with $a=0.58$ confirms it: a
large amplitude with a small exponent minus a constant is a degenerate
reparametrisation of a slowly varying function, not a floor. E2 is the
one to believe, and the pre-registration says so in advance.

One comparison in that output is mis-stated and is corrected here
rather than left standing: the three-parameter fit's residual r.m.s.
$0.02085$ is in kurtosis units while the floorless fit's $0.10062$ is
in log units, so the two are not comparable and the printed line
reads as though the floor fit were five times better. It is not a
comparison at all.

**E4 broke too, and it is the more consequential break.** The coin's
own mean excess kurtosis at the two new bands is $-0.00848$ and
$-0.01019$, outside the $-0.02103$ to $-0.01376$ that [rem:cnkurt]
measured. The control drifts toward zero rather than sitting at a
constant offset. And with it goes the separation: at $b=25$ the real
point scores $z=1.98$ against its coin ensemble, where every band up
to $b=23$ scored at least $5.0$.

That is this branch's floor, and it is not a budget problem. The
coin's *draw-to-draw* spread in excess kurtosis is $0.018$ at
$b=23$ and $0.018$ at $b=25$ — it does not shrink as $N$ grows,
because it is set by the field's correlation across $N$ and not by
sample size. The real signal decays like $N^{-0.73}$. **The two meet
at $N\sim3\cdot10^7$, and pushing $N$ further cannot separate them
because only one of the two is falling.** More draws sharpen the
estimate of where the coin's centre is; they do not narrow the
ensemble the real point has to escape.

**And one sentence of [rem:cnkurt] is withdrawn.** It read that the
decay is slower than an independent-term model's $N^{-1}$. That
$N^{-1}$ comes from the excess kurtosis of $\sum w_v\varepsilon_v$ as
a random variable *over the signs*, and what is measured here is the
excess kurtosis of the empirical distribution *over $N$* of one fixed
sign pattern. They are different objects, and the computed
independent-sign null — the coin ensemble, which was there all along —
does not behave like $N^{-1}$ at all: it sits near $-0.01$ and drifts
slowly toward zero. The comparison was against a formula that does not
govern the measured statistic. What replaces it is the paragraph
above: the real arm decays toward the coin arm's level and reaches it
at $b\approx25$. The sentence in
`results/audit_cn_kurt_drift.txt` that makes the same comparison is
superseded by this one.

So the reading after nine octaves is narrower than after seven, and
better founded. The excess kurtosis of $C(N)$ is real and outside its
coin ensemble from $b=17$ to $b=23$; it follows a power law well
enough to forecast two octaves out of sample; it has no floor this
design can see, and it reaches the coin's own noise at $b\approx25$,
which is where this instrument stops.

#### Remark (the control does shrink, and it does not help) {#rem:cnskew}
<!-- evidence: audit_cn_skew.py -->

[rem:cnkurtlimit] closed the kurtosis route by finding its floor: the
coin's draw-to-draw spread in excess kurtosis does not shrink with $N$
— it is set by the field's correlation across $N$, not by sample size
— while the real signal decays like $N^{-0.73}$, so the two meet at
$b\approx25$ and more $N$ cannot separate them. What was needed next
was not more $N$ but a statistic whose control *does* shrink.

The skewness was the named candidate, and there is a structural reason
it should behave differently: sending $\varepsilon\to-\varepsilon$
sends $G_{\text{coin}}\to-G_{\text{coin}}$, and skewness is odd, so the
coin's skewness distribution is exactly symmetric about zero — its
centre is a construction, not a measurement. S5 confirms the
instrument: over seven octaves the coin's mean skewness is never more
than $1.46$ of its own standard error from zero.

**The reasoning was right and the conclusion is negative.**

| $b$ | real skew | coin sd | $z$ |
|---|---|---|---|
| $17$ | $-0.47174$ | $0.03691$ | $-12.73$ |
| $19$ | $-0.21666$ | $0.05067$ | $-4.28$ |
| $21$ | $-0.09094$ | $0.03870$ | $-2.49$ |
| $23$ | $-0.03244$ | $0.02489$ | $-1.20$ |

The asymmetry is negative at every octave, so [rem:cnlaw]'s heavier
left tail is a property of $C$ and not of its band. **And the control
shrinks** (S3): the coin's skew spread goes like $N^{-0.1198}$,
resolved at $t=-2.32$, against a kurtosis spread that was flat. That
is exactly the property the kurtosis lacked.

**It buys nothing** (S4). The signal shrinks like $N^{-0.6367}$, five
times faster than the control, and $|z|$ falls from $12.73$ to $1.20$
across the seven octaves, a slope resolved at $t=-4.68$. S2 breaks
with it: the separation is already below $3$ at $b=21$ and gone by
$b=23$. **The skewness closes earlier than the kurtosis did**, at
$N\sim10^7$ rather than $3\cdot10^7$, despite having the better
control. The pre-registration wrote this case out in advance — a
control that shrinks with a signal that shrinks faster is the same
wall in a different place — and that is where it landed.

Two cases are not a law, but they are the same shape twice, and it is
worth stating as the thing to test rather than leaving implicit: in
both, $G$'s departure from Gaussianity decays like $N^{-0.6}$ to
$N^{-0.7}$ while the coin ensemble's own fluctuation in the same
statistic decays like $N^{-0.12}$ or not at all. **If that is general,
no standardised moment of $G$ separates $\mu$ from a coin
asymptotically, and the finite-$N$ separations this branch measured
are all of them.**

What follows is a change of object rather than of statistic. Every
quantity in this branch so far is a functional of the *marginal*
collection $\{G(N)\}$, and the coin matches those because it has
$\mu$'s support and independent signs, which is enough for a marginal.
What it does not have is $\mu$'s multiplicative structure linking
different $N$ — $G(N)$ against $G(2N)$, or against $G(N')$ with the
same radical. Those are correlations of the field, not moments of its
margin, and nothing here has measured one.

#### Remark (the coin reproduces the two-point function too, and the test was the wrong kind) {#rem:cnshift}
<!-- evidence: audit_cn_shiftcorr.py -->

[rem:cnskew] closed the second marginal statistic and named what
should replace it: every quantity measured in this branch is a
functional of the marginal collection $\{G(N)\}$, and the coin has
$\mu$'s support and independent signs, which is all a marginal needs.
What the coin lacks is $\mu$'s structure linking different $N$. The
shift correlation
$\rho(h)=\operatorname{mean}_N G(N)G(N+h)$ is the natural place to
look, and for a reason: averaging the coin over its signs kills every
off-diagonal $\mu(v)\mu(w)$ term, leaving the diagonal
$\sum_v\mu^2(v)\Lambda(N-v)\Lambda(N+h-v)$ exactly. **The off-diagonal
is where $\mu$'s own correlations live.**

They do not show. Over $32$ even shifts $h=2\ldots64$ at $b=23$, the
real arm's largest $|z|$ against the coin ensemble is $2.23$, while
the coin draws' own largest — each scored leave-one-out against the
other $31$ — run from $0.73$ to $4.86$, and **four of the $32$ draws
reach or exceed the real arm's value** (T2 refuted). At $b=21$ the
same comparison gives $1.56$ against a coin range of $0.74$ to $3.72$.
The real field sits inside its control at both octaves.

T3 held — the largest $|z|$ rose from $1.56$ to $2.23$ — and holding
inside a refuted frame carries no information: the coin's own range
rose with it, from a ceiling of $3.72$ to $4.86$. It is recorded as
held because that is what was registered, and read as nothing.

**T4 broke, and it says why.** The coin's spread in $\rho(h)$ falls by
a factor of $1.3608$ between the two octaves, where four times as many
$N$ would give $2$ if the spread were sample-limited. It is not: like
the excess kurtosis, this control is set by the field's correlation
across $N$, not by how many $N$ are averaged. So even had the real arm
escaped, the escape would not have widened with $N$.

That is three routes closed, and together they say something sharper
than any of them alone. **What the coin reproduces is not merely the
margin.** It has $\mu^2$ — the support — and the shift correlation is
dominated by a diagonal that $\mu^2$ fixes completely; the actual
signs contribute nothing that rises above the control's own
fluctuation at these ranges. [lem:coin]'s reach is wider than it was
stated for.

**And the test asked the wrong question.** [rem:cnskew] named
$G(N)$ against $G(2N)$, and against $G(N')$ of the same radical —
*multiplicative* relations, because $\mu$ is multiplicative and
$\mu(2v)$ is tied to $\mu(v)$. What was measured here is the
*additive* two-point function, $N$ against $N+h$. The diagonal that
the coin fixes is exactly the additive overlap of two $\Lambda$-shifts,
so an additive test was always the one the coin was best placed to
survive. Nothing here bears on the multiplicative question, which
remains as [rem:cnskew] left it, and the gap between what that remark
named and what this one ran is the author's and is recorded as such.

#### Remark (the multiplicative question, asked properly, closes too) {#rem:cndilation}
<!-- evidence: audit_cn_dilation.py -->

[rem:cnshift] recorded that it had asked the wrong question: the
additive two-point function is the one the coin is best placed to
survive, because the diagonal the coin fixes exactly *is* the additive
overlap of two $\Lambda$-shifts. What [rem:cnskew] had named was
multiplicative — $G(N)$ against $G(dN)$ — because $\mu$ is
multiplicative and the coin's $\varepsilon$ links nothing to anything.
This is that question, with
$\rho(d)=\operatorname{mean}_N G(N)G(dN)$ over even $N$ in
$(2^{20},2^{21}]$ and $d=2,\ldots,9$.

The algebra says the structure is there. Among the off-diagonal
$\mu(v)\mu(w)$ terms are the pairs $w=dv$, and at $d=2$ they are
exact: $\mu(2v)=-\mu(v)$ for $v$ odd, $0$ for $v$ even. The coin has
zero of that by construction. What the algebra does not say is
whether it is visible after normalisation and averaging, and it is
not.

**U2 broke, and it broke the wrong way round.** The real arm's largest
$|z|$ over the eight dilations is $0.72$, at $d=6$. Each coin draw's
own largest, scored leave-one-out against the other $31$, runs from
$0.80$ to $3.98$ — **all $32$ of them exceed the real arm's.** The
field is not merely inside its control; it is further inside than
every draw of the control. U3 broke with it ($d=6$, not $2$) and is
recorded without being read, since a location means nothing once the
escape is absent.

That is four routes closed, and this one was the route the algebra
most favoured. The pre-registration wrote down what that would mean
before the run: **the coin's agreement with $C(N)$ at these ranges is
not an artefact of which statistics were chosen.** Marginal third
moment, marginal fourth moment, additive two-point, multiplicative
two-point — every one is reproduced by independent signs on $\mu$'s
support.

U4 was an instrument check and it broke, so its power is reported
rather than left implicit. It asked whether the coin's mean $\rho(d)$
is positive at every $d$ and falling; the means run $+0.004225$ to
$+0.035083$ with one negative at $-0.006715$, against standard errors
of $0.003923$ to $0.012706$. **Only three of the eight are resolved
away from zero at all**, and the negative one is $1.51$ of its own
standard error. The check asked whether eight numbers are ordered when
five of them are not distinguishable from zero; its failure is
evidence that $32$ draws cannot say, not that $\rho(d)$ is other than
the overlap. U2 does not rest on it — scoring each draw leave-one-out
against the rest is self-calibrating whatever $\rho(d)$ is.

What this does **not** say. It does not say $\mu$ is a coin, and it
does not license the coin as a model in any proof. It says that at
$N\sim10^6$, across these four families of statistic, no measurement
in this repository distinguishes $C(N)$ from a field with $\mu$'s
support and independent signs. The distinguishing structure exists —
$\mu(2v)=-\mu(v)$ is not a conjecture — and the finding is that it
does not surface in these quantities at this size. [lem:coin]'s reach
is wider than it was stated for, and the demand-side measurement
programme on $C(N)$ has met a wall that more $N$ does not move.

#### Remark (the non-Gaussianity is not made anywhere) {#rem:cnkurtwhere}
<!-- evidence: audit_cn_kurt_where.py -->

Four routes closed, and one thing did not: the finite-$N$ separation
in excess kurtosis is real and large — $z=16.5$ at $b=20$, $33.4$ at
$b=17$ — and had no mechanism. Writing
$C(N)=\sum_{p<N}\log p\,\mu(N-p)$ and splitting the sum by $v=N-p$
dyadically splits it by how far the prime sits below $N$: small $v$ is
few terms of weight about $\log N$, large $v$ is many terms of small
weight. If one part makes the non-Gaussianity, that is a mechanism.

The decomposition is exact — the pieces sum to $C$ to $1.455\cdot
10^{-11}$ across the band, and the whole band reproduces
$0.41549$ (W1). Each piece is scored against its own ensemble of $32$
sign patterns restricted to the same $v$-range, so term counts and
weights are matched.

| $v$-piece | terms | real | coin mean | $z$ |
|---|---|---|---|---|
| $(2^{10},2^{11}]$ | $621$ | $+0.09297$ | $-0.04258$ | $4.96$ |
| $(2^{12},2^{13}]$ | $2491$ | $+0.13139$ | $-0.01907$ | $6.27$ |
| $(2^{14},2^{15}]$ | $9958$ | $+0.04683$ | $-0.01686$ | $1.68$ |
| $(2^{17},2^{18}]$ | $79672$ | $+0.01047$ | $-0.02350$ | $1.17$ |
| $(2^{20},2^{21}]$ | $637457$ | $+0.01005$ | $-0.01465$ | $0.54$ |

**W2 broke: it is not concentrated.** The largest per-piece $|z|$ is
$6.27$, well below the whole field's $16.5$. A quantity made in one
place would show a piece at least as separated as the whole; instead
fourteen pieces each carry a little and the whole carries more than
any of them. The pre-registration fixed how to read that before the
run — there is no mechanism of this kind to find, and that answers the
question rather than leaving it open.

**W4 broke, and for a reason worth recording as an error.** It
predicted that some piece would show a coin excess kurtosis above
$0.5$, on the thought that a piece with few terms is non-Gaussian for
any signs. The sign is wrong: a weighted sum of independent $\pm1$ is
*platykurtic*, with excess kurtosis about $-2\sum w^4/(\sum w^2)^2$,
and the measured coin means run $-0.11969$ to $-0.01102$ — negative
everywhere and largest in magnitude at the smallest piece, exactly as
that formula says and opposite to what was predicted. No piece of the
fourteen reaches $0.5$.

W3 held and is not read, since a location means nothing once the
concentration is absent.

**WITHDRAWN.** This paragraph observed that the per-piece $|z|$
is not flat — $2.90$ to $6.27$ over $v\le2^{14}$ and $0.54$ to
$1.74$ above it — and said the way to make it evidence is to
predict the same profile at a band this run never saw.
[rem:cnwherereach] did that and the profile is not a profile: a
piece's $|z|$ is a function of how many terms it holds, fitted
here and forecasting two unseen bands with an out-of-sample
r.m.s. error of $0.8194$ against a constant's $1.5392$. The
upper half of the ladder looks empty because those pieces hold
the most terms, and nothing about where the prime sits below
$N$ is involved.

#### Remark (there is no profile in v, only a count — and the coin's own law is not the textbook one) {#rem:cnwherereach}
<!-- evidence: audit_cn_where_reach.py -->

[rem:cnkurtwhere] left an observation it refused to call evidence: the
per-piece $|z|$ was not flat, running $2.90$ to $6.27$ over $v\le2^{14}$
and $0.54$ to $1.74$ above it. Two readings fit that shape — a profile
in $v$, meaning where the prime sits below $N$ matters; or no profile
at all, the $|z|$ being a function of how many terms the piece holds,
which for a dyadic piece is about $2^j/\log N$ and so nearly the same
at every band. The second explains the shape with no arithmetic, so it
is the one to try to confirm.

Fitted on $b=20$'s fourteen pieces alone, $|z| = 6.92713 - 0.32369
\log_2(\text{terms})$ with residual r.m.s. $1.36790$. That fit then
forecast the pieces of $b=18$ and $b=19$, which it had never seen, and
both bands first reproduced their published whole-band kurtosis
exactly — $1.06049$ and $0.69806$ (Y1).

**All $25$ pieces landed inside the two-standard-error forecast**
(Y2). That alone is a weak pass and is reported as one: the interval
is about $\pm2.7$ while the new $|z|$ span $0.38$ to $5.33$, so a
constant would catch most of them too. Against a constant at $b=20$'s
mean $|z|$ of $2.7662$, on the same $25$ unseen pieces:

| forecast | out-of-sample r.m.s. error |
|---|---|
| from the term count | $0.8194$ |
| constant | $1.5392$ |

a ratio of $1.8784$. The count does not merely fail to be contradicted
— it predicts the unseen pieces nearly twice as well as knowing
nothing.

**So [rem:cnkurtwhere]'s observation is withdrawn.** There is no
profile in $v$. A piece's separation from its own coin is set by how
many terms it holds, and the apparent emptiness of the upper half of
the ladder is that those pieces hold the most terms. Nothing about
where the prime sits below $N$ is involved. The non-concentration
reproduces out of sample too (Y3): largest piece $5.33$ against the
whole field's $39.8$ at $b=18$, and $4.61$ against $20.7$ at $b=19$.

**Y4 tested the formula that W4's failure identified, and the formula
is wrong in the same way twice.** The coin's per-piece excess kurtosis
is negative at all $25$ pieces, as $-2\sum w^4/(\sum w^2)^2$ requires,
and it falls with the term count at $t=-5.62$. But that expression
predicts a slope of $-1$ against $\log_2(\text{terms})$ and the
measured slope is $-0.17584\pm0.03130$. The reason is the one
[rem:cnkurtlimit] recorded when withdrawing a different comparison:
$-2\sum w^4/(\sum w^2)^2$ is the excess kurtosis of $\sum
w_v\varepsilon_v$ *as a random variable over the signs*, and what is
measured is the excess kurtosis of the empirical distribution *over
$N$* of one fixed sign pattern. **The same confusion has now produced
one withdrawn paragraph and one broken prediction**, and the size of
the gap is on the record: the empirical-over-$N$ kurtosis dies with
the term count about six times more slowly than the marginal formula
says.

#### Remark (the last axis is not closed: C(N) has class structure V does not carry) {#rem:cnclass}
<!-- evidence: audit_cn_class.py -->

Five routes in this branch closed, and every one used $N$ only as a
band. The axis left was the one the coin cannot follow: its
$\varepsilon$ does not know $N$. Cutting the even $N$ of
$(2^{20},2^{21}]$ by which of $3,5,7$ divide them gives eight classes,
and within each the first and second moments of $G=C/\sqrt{V}$ are
measured against $32$ sign patterns. Averaged over signs the coin
gives $E[G^2]=1$ in *every* class exactly, since $V(N)$ is the exact
second moment at each $N$ — so the null here is no structure at all,
not merely no difference between classes.

| $3\mid N$ | $5\mid N$ | $7\mid N$ | $N$ | $E[G^2]$ | coin | $z$ |
|---|---|---|---|---|---|---|
| — | — | — | $239675$ | $0.80803$ | $1.02249$ | $-0.90$ |
| 3 | — | — | $119837$ | $0.83798$ | $1.03152$ | $-0.78$ |
| 3 | 5 | — | $29960$ | $1.16445$ | $1.05787$ | $0.31$ |
| 3 | — | 7 | $19973$ | $0.91322$ | $1.03594$ | $-0.49$ |
| — | 5 | 7 | $9987$ | $0.95176$ | $1.07534$ | $-0.27$ |
| **3** | **5** | **7** | $4993$ | $\mathbf{4.37497}$ | $1.02257$ | $\mathbf{9.01}$ |

**Z2 holds against the iid coin.**  It was read as this branch's
first escape and that reading is superseded: [rem:cnmultdeep]
shows a resolved random-multiplicative ensemble reaches it, so
what the coin could not make, multiplicativity can. The real arm's
largest $|z|$ over classes is $9.01$; the coin draws' own largest, each
taken leave-one-out against the other $31$, run $0.54$ to $6.46$, and
none of the $32$ reaches it. The bookkeeping gate is exact: the
size-weighted average of the class second moments equals the whole
band's $0.875667582379$ to twelve places (Z1).

**The departure is one class, $105\mid N$, and its arithmetic is not
in doubt.** If $3,5,7$ all divide $N$ then for every prime $p>7$ the
argument $N-p$ is coprime to $105$: $N-p\equiv-p$, and $p$ is a unit
modulo each. So $C(N)=\sum_p\log p\,\mu(N-p)$ runs $\mu$ over
integers restricted to a class of $105$, while $V(N)$ counts
$\mu^2(N-p)\log^2p$ without regard to it. That restriction is a
deterministic local factor of exactly the kind [conj:L] calls the mask
$\mathbf M$. **So this is not a refutation of that conjecture; it is a
demonstration that its $\mathbf M$ is not absorbed into
[prop:V]'s $V$.** $V$ is the coin's second moment and the field's is
larger by a local factor at highly composite $N$ — here by a factor of
$5.41$ against the class where none of the three divides.

The size is measured and not attributed. Nothing here checks it
against a singular series or any other formula, and the obvious next
question is whether it is one.

Two cautions, both against this remark's own result. **Z3 holds by
$0.03$**: the first moment's largest $|z|$ is $3.20$ against a coin
ceiling of $3.17$, so the escape in the mean is a hair and should be
read as untested rather than shown. And **Z4 is REFUTED** — the median
coin draw's largest $|z|$ over classes is $1.09$ on the second moment
against a predicted floor of $1.5$, so the control was less necessary
than the design assumed; the $105\mid N$ class would have stood out in
a raw table. That does not weaken Z2, which is scored against the
ensemble either way, but it does mean this design was more cautious
than it needed to be.

One band and three primes. Whether the factor persists at other $N$,
whether it grows with the number of small primes dividing $N$ — the
lower rows hint at $1.16$ at $15\mid N$ and $0.95$ at $35\mid N$
against $4.37$ at $105\mid N$, which is not a pattern any of this
tested — and whether it matches a computable local density, are three
separate measurements and none is made here.

#### Remark (the escape is a shift, it lives at every band, and it fades slowest of all) {#rem:cnclassreach}
<!-- evidence: audit_cn_class_reach.py -->

[rem:cnclass] found what was read as this branch's first escape --
[rem:cnmultdeep] later showed the reading was against too weak a
null -- and left the two questions that decide what the
measurement is worth.  Those two are answered below and their
answers do not depend on the withdrawn reading. The first is answered by
arithmetic on its own numbers: of the $E[G^2]=4.37497$ at $105\mid N$,
the class mean $-1.70728$ contributes $2.9148$, leaving a variance of
$1.4602$. **It is not a few enormous $|G|$ but $4993$ values of $N$
sitting about $1.7$ standard deviations low together** — a shift, not
a spread.

Six bands, $b=17\ldots22$, answer the rest. The gate reproduces
$4.37497$ exactly (A1).

| $b$ | class $N$ | mean $G$ | $E[G^2]$ | $z$ | coin's own largest |
|---|---|---|---|---|---|
| $17$ | $624$ | $-2.91115$ | $10.20413$ | $18.01$ | $0.47$–$8.24$ |
| $18$ | $1248$ | $-2.42550$ | $7.69269$ | $11.80$ | $0.45$–$5.48$ |
| $19$ | $2497$ | $-2.03591$ | $5.79633$ | $18.05$ | $0.63$–$5.19$ |
| $20$ | $4993$ | $-1.70728$ | $4.37497$ | $17.94$ | $0.64$–$6.04$ |
| $21$ | $9986$ | $-1.42183$ | $3.42435$ | $12.68$ | $0.57$–$4.25$ |
| $22$ | $19973$ | $-1.16347$ | $2.67490$ | $5.38$ | $0.50$–$5.06$ |

**The escape is not one band's** (A2): at every band, none of the $32$
coin draws' own largest $|z|$ over classes reaches the real arm's.
**And the shift keeps its sign** (A3), negative at all six — at
$b=17$ it is $-2.91$, so $C(N)$ sits nearly three of its own standard
deviations below zero, systematically, for $N$ divisible by $105$.

**A4 broke.** The shift fades: $\log|{\rm mean}\,G|$ against $\log_2N$
has slope $-0.18183\pm0.00203$, so the shift goes like $N^{-0.2623}$.
The pre-registration fixed the reading — a sixth finite-$N$
separation, and the exponent is what this run adds.

Two things about that exponent are worth separating from the $t$ of
$89.38$, which is not evidence and is not offered as any. Six bands
cut from one sieve are not independent samples, so that standard error
is a lower bound. What the fit does support is straightness: the
residuals in $\log|{\rm mean}|$ are $-0.00332$, $-0.00400$, $+0.00273$,
$+0.00852$, $+0.00739$, $-0.01131$, an r.m.s. of $0.00695$ against a
mean $|\log|$ of $0.6173$. **Six octaves lie on a line to one part in
ninety.**

And the exponent is the slow one. Every other separation this branch
measured decayed like $N^{-0.64}$ ([rem:cnskew]) or $N^{-0.73}$
([rem:cnkurtlimit]). This one decays at $N^{-0.26}$, two to three
times more slowly, and the separation itself — the $|z|$, which
carries the control's own shrinking with it — goes like $N^{-0.2402}$.
From $5.38$ at $b=22$ that reaches $3$ near $b=25.5$. **That is an
extrapolation and not a measurement**, printed only so the next band
to try has a name.

So the sixth separation is still a finite-$N$ separation, and this
branch has now produced six of them and no seventh kind. What is
different here is only the rate, and what remains untouched is the
same thing [rem:cnclass] left: nothing has checked the shift against a
computable local density, at this class or any other.

#### Remark (the shift is not made anywhere: it is what survives a 198-fold cancellation) {#rem:cnclassomega}
<!-- evidence: audit_cn_class_omega.py -->

OPEN.md item 1 asks whether the $105\mid N$ shift is a computable local
density. It cannot be asked from [rem:cnclass]'s table, because
fitting a multiplicative law needs resolved classes and seven of the
eight are noise — the first-moment $z$ runs $-0.00,-0.09,+0.38,-0.94,
+0.17,-0.54,+0.87$, with only $-3.20$ at $105\mid N$. What has to come
first is where in the sum the shift is made.

The sum splits exactly. For squarefree $m$, $\mu(m)=(-1)^{\omega(m)}$,
so $C(N)=\sum_j(-1)^jS_j(N)$ with
$S_j(N)=\sum_{\omega(m)=j,\ \mu^2(m)=1}\Lambda(N-m)\ge0$. The
alternating sum rebuilds $C$ to $7.194\cdot10^{-10}$ across the band
and the class mean returns $-1.70728$ (B1).

| $j$ | $m$ in bucket | contribution to the class mean |
|---|---|---|
| $1$ | $155611$ | $-115.48990$ |
| $2$ | $425951$ | $+164.60003$ |
| $3$ | $436477$ | $-54.12597$ |
| $4$ | $207617$ | $+3.31153$ |
| $5$ | $45500$ | $-0.00396$ |

**B2 passed its cap and the cap does not measure what it was written
for.** It asked whether one bucket carries more than half the shift;
four of them clear that, at $67.65$, $96.41$, $31.70$ and $1.94$ times
the total, because any bucket larger than the residue clears it. The
buckets sum in absolute value to $337.53238$ and cancel to
$-1.70728$: **a cancellation of $197.70$ to one.** B3 passed on the
same vacuous footing and is recorded without being read.

So the answer is the one B2's refutation clause was written for even
though B2 did not formally break: the shift is not made in a place. It
is the residue of a near-total cancellation between $\omega$-buckets
each one to two orders of magnitude larger than it. The
$\omega$-distribution of the shifted primes $N-p$ peaks at $j=2$ with
a mean near $2.65$, the alternating sum of that distribution is small
by parity, and **the class shift is a $0.5$ per cent asymmetry in how
completely it cancels.**

That sharpens item 1's question and makes it much harder. A local
density describing this shift would have to predict the parity balance
of $\omega(N-p)$ for $N$ divisible by $105$ to about one part in two
hundred. Nothing in this repository predicts anything about
$\omega(N-p)$ at all.

B4 was the instrument check and it broke, at $2.13$ standard errors in
the $j=1$ bucket. **That is a threshold crossed, not an instrument
shown broken**: eight buckets at a two-standard-error threshold fire by
chance about one time in three. The refutation is recorded and read as
the unresolved case — the wording of B4's rule was amended after the
run to say so, disclosed in the script, with the cap and the verdict
unchanged.

#### Remark (a multiplicative null, and three things it showed that were not the question) {#rem:cnmultnull}
<!-- evidence: audit_cn_multnull.py -->

Every control in this branch had been iid signs on $\mu$'s support,
which is a weak null for a multiplicative function. The natural one
keeps the structure and changes the values: $f(p)=\pm1$ iid,
$f(m)=\prod_{p\mid m}f(p)$ on squarefree $m$, zero otherwise, so
$|f|=\mu^2$ exactly and $f$ is multiplicative exactly as $\mu$ is.
Anything the escape owes to multiplicativity plus the coprimality
constraint at $105\mid N$, this ensemble has too.

| $3\mid N$ | $5\mid N$ | $7\mid N$ | mean $G$ | ensemble mean | ensemble sd | $z$ |
|---|---|---|---|---|---|---|
| — | — | — | $+0.01376$ | $+0.02068$ | $0.33637$ | $-0.02$ |
| 3 | 5 | — | $-0.46174$ | $+0.06086$ | $0.44825$ | $-1.17$ |
| — | 5 | 7 | $+0.43997$ | $+0.01167$ | $0.41757$ | $1.03$ |
| **3** | **5** | **7** | $-1.70728$ | $+0.03422$ | $0.50217$ | $-3.47$ |

**C2 held by one draw of thirty-two.** The real arm's largest $|z|$ is
$3.47$; the draws' own largest run $0.34$ to $3.57$, and exactly one
reaches it. The pre-registration said a C2 that holds by one draw is
not the same as one that holds by twenty and that the count would be
reported, so: $\mu$ sits at about the $97$th percentile of this
ensemble, inside it and near the edge.

**C3 broke and the direction is the interesting part.** The
multiplicative ensemble's spread at $105\mid N$ is $0.50217$ against
the iid coin's $0.5646$ — a factor of $0.89$, *smaller*, where three
times larger was predicted. Making the null multiplicative did not
make it harder to escape. The two nulls behave almost identically on
this statistic.

**And the reading written in advance for "C2 holds" was wrong.** It
said a random multiplicative sign function "shows the same class
structure". It does not: the ensemble's class means are $+0.012$ to
$+0.064$, flat and at zero, and C4 confirms it is centred. What the
ensemble has is width, not structure — it covers $\mu$ by wandering,
not by going where $\mu$ goes. The error was assuming that a null
which covers a point must do so by reproducing what makes the point
extreme. That sentence stands in
`results/audit_cn_multnull.txt` as written and is superseded here.

**What this run did not test is the escape itself.** [rem:cnclass]'s
$z=9.01$ was on the *second* moment $E[G^2]$; the first moment was
already marginal there, holding by $3.20$ against a coin ceiling of
$3.17$, and that remark said so and said to read it as untested. This
script compared class *means* only. So it has re-run the marginal
comparison against a second null and found it still marginal, and has
left the $9.01$ untouched. That is a shortfall in the design and not
a finding, and the second-moment comparison against this ensemble is
the measurement that should have been made.

#### Remark (mu is outside the multiplicative null — and the same run says not to build on it yet) {#rem:cnmultm2}
<!-- evidence: audit_cn_multnull_m2.py -->

[rem:cnmultnull] built the right null and compared the wrong
statistic. This is the comparison that was missing: the second moment,
where [rem:cnclass]'s escape actually lives, against the random
multiplicative ensemble.

| $3\mid N$ | $5\mid N$ | $7\mid N$ | $E[G^2]$ | ensemble mean | ensemble sd | $z$ |
|---|---|---|---|---|---|---|
| — | — | — | $0.80803$ | $0.90346$ | $0.16415$ | $-0.58$ |
| 3 | 5 | — | $1.16445$ | $0.93021$ | $0.29866$ | $0.78$ |
| — | 5 | 7 | $0.95176$ | $0.85394$ | $0.14589$ | $0.67$ |
| **3** | **5** | **7** | $\mathbf{4.37497}$ | $0.90114$ | $0.29369$ | $\mathbf{11.83}$ |

**D2 is REFUTED, and the caveat below turned out to matter more
than the refutation.**  [rem:cnmultdeep] resolved the ensemble to
$512$ draws and the $11.83$ became $5.30$ with five draws
reaching it; what follows is recorded as measured and its
reading is superseded there. The real arm's
largest $|z|$ is $11.83$; the draws' own largest, leave-one-out, run
$0.70$ to $6.63$, and none of the $32$ reaches it. A random
multiplicative sign function has $\mu$'s support, $\mu$'s
multiplicativity, and $\mu$'s coprimality constraint at $105\mid N$,
and it does not go where $\mu$ goes. D3 held: the two nulls are of
similar strength here, the ensemble's spread being $0.79$ of the iid
coin's.

**And D4 is REFUTED in the way its own rule said would carry a
reading.** For a random multiplicative $f$, $E[C_f(N)^2]=V(N)$ is
exact — $f(m)f(m')$ is a product over the symmetric difference of two
prime sets, so it has mean zero unless $m=m'$ — hence $E[G_f^2]=1$ in
every class at every $N$. The measured ensemble means are $0.85394$ to
$0.93021$, **every one of the eight below $1$**, the furthest at
$5.66$ of its own standard error. D4's rule named the
too-noisy-to-tell case and named what would escape it: a crossing
repeated across classes. Eight of eight is that.

So the ensemble does not sit where the algebra puts it, and the
explanation available is that $G_f^2$ is heavy-tailed, so a $32$-draw
sample mean lands systematically below a mean carried by rare large
draws. **That is a hypothesis, not something measured here**, and it
bears directly on D2: the $z$ is computed against a sample mean that
is too low and a sample standard deviation that under-resolves a heavy
tail, so $11.83$ is an overstatement of an unknown size, and "none of
$32$" says less about a tail than about $32$.

The honest position is therefore narrower than the headline. What is
measured: on the statistic where the escape lives, $\mu$ is outside
what $32$ draws of the strongest available null produce, by a margin
the same $32$ draws cannot calibrate. What is not established: that it
would remain outside a properly resolved ensemble. **The fix is more
draws**, until the ensemble's class means come to $1$ as the identity
requires, and that is the measurement this remark asks for rather than
a conclusion it offers.

#### Remark (the escape was the tail: seven routes, none of them open) {#rem:cnmultdeep}
<!-- evidence: audit_cn_multnull_deep.py -->

[rem:cnmultm2] measured $\mu$ outside the multiplicative null at
$|z|=11.83$ and, in the same run, found the ensemble sitting at
$0.854$–$0.930$ where the identity $E[C_f(N)^2]=V(N)$ puts it at $1$.
It wrote the honest position — what is measured is that $\mu$ is
outside what $32$ draws produce, by a margin those $32$ draws cannot
calibrate — and asked for this run. $512$ draws answer it.

**E2 caught a bug before it caught anything else, which is what it was
for.** The first execution put the ensemble at $1387$ to $1795$
instead of $1$. The cause was in this script: masks drawn from
`integers(0, 1 << 63)` never set bit $63$, so one draw per pass had
every sign $+1$ and $C_f=\sum_m\mu^2(m)\Lambda(N-m)$, eight degenerate
draws swamping $512$. E2's refutation rule had named exactly that
possibility — "the construction of $f$ is not what it is taken to be"
— and E3's jump between $32$ and $64$ draws located the first one.
The masks are now built from two $32$-bit halves; the numbers below
are the corrected execution of the same pre-registration, and the
earlier two multiplicative runs are unaffected, their masks having
been drawn over the full $2^{32}$.

| $3\mid N$ | $5\mid N$ | $7\mid N$ | $E[G^2]$ | ensemble mean | ensemble sd | $z$ |
|---|---|---|---|---|---|---|
| — | — | — | $0.80803$ | $1.04608$ | $1.27572$ | — |
| 3 | 5 | — | $1.16445$ | $1.04262$ | $0.85337$ | — |
| **3** | **5** | **7** | $4.37497$ | $1.00609$ | $0.63553$ | $\mathbf{5.30}$ |

**The ensemble now sits where the identity puts it** (E2): the class
means run $1.00290$ to $1.06459$, the furthest $1.30$ standard errors
from $1$. **And it got there by converging** (E3): the mean deviation
falls $0.13591$, $0.10868$, $0.03961$, $0.02859$, $0.02769$ as the
draws go $32$ to $512$. E2 on its own would have been the
too-weak-to-tell case — a wide ensemble makes the standard error large
and any centre passes — and E3 is what makes it mean anything.

**With the instrument fixed, the escape closes** (E4). The real arm's
largest $|z|$ is $5.30$, not $11.83$: at $512$ draws the ensemble's
spread at $105\mid N$ is $0.63553$ against $0.29369$ at $32$, and
$5$ of the $512$ draws reach or exceed the real arm, their own largest
running $0.12$ to $30.88$. $\mu$ sits near the $99$th percentile of
the resolved multiplicative ensemble — inside it, near the edge,
which is where [rem:cnmultnull] had already put it on the first
moment at the $97$th.

So the seventh route closes with the other six, and
[rem:cnmultm2]'s caveat was right about its own headline. **What
remains true, and is worth separating from the withdrawn reading, is
everything that was measured rather than compared.** The class
structure at $105\mid N$ is real: $E[G^2]=4.37497$ against $0.80803$
where none of $3,5,7$ divides, a shift of $-1.70728$ in the mean,
present at all six bands [rem:cnclassreach] measured and fading like
$N^{-0.2623}$. $V(N)$ is not the right scale for highly composite $N$.
**What has changed is why**: it is what multiplicativity and the
coprimality constraint produce, not something particular to $\mu$ —
a random multiplicative sign function does the same. Read forward,
that is a statement about [conj:L]'s mask $\mathbf M$: the mask is a
consequence of multiplicativity, and $\mu$'s own values add nothing to
it that $512$ draws can see.

#### Remark (it was the null, not the depth: the iid coin is narrower than the right one) {#rem:cncoindeep}
<!-- evidence: audit_cn_coin_deep.py -->

[rem:cnmultdeep] changed two things at once — $32$ draws to $512$, and
iid signs to multiplicative ones — and nothing there separated which
did the work. The separation is one run: the same statistic, the same
band, the same $512$ draws, iid signs.

| | iid, $512$ | multiplicative, $512$ |
|---|---|---|
| ensemble mean at $105\mid N$ | $1.00642$ | $1.00609$ |
| ensemble spread there | $0.41025$ | $0.63553$ |
| real arm's largest $\lvert z\rvert$ | $8.21$ | $5.30$ |
| draws reaching it | $0$ of $512$ | $5$ of $512$ |

Both ensembles satisfy their own identity (F2): the iid class means
run $0.98707$ to $1.01534$, the furthest $1.01$ standard errors from
$1$, and the deviation falls $0.06330$, $0.02551$, $0.02300$,
$0.00409$, $0.00638$ as the draws go $32$ to $512$. So both are
resolved and they still disagree.

**F3 is REFUTED and F4 with it.** At full resolution the iid coin
*still* excludes $\mu$ — none of $512$ draws reaches $8.21$, their own
largest running $0.44$ to $7.73$ — while the multiplicative ensemble
covers it. The iid spread is $0.6455$ of the multiplicative one. **The
null type did the work, not the depth.**

That has a consequence outside this branch, and the pre-registration
said to state it plainly. [lem:coin] is an iid coin, and it is the
control that has sunk claim after claim here. Where the object under
test is multiplicative — $\mu$ is, and most of what this repository
measures is built from it — **an iid coin is narrower than the right
null**, so every significance calibrated against one is larger than it
should be.

The damage is bounded and the bound is worth stating exactly. A
narrower null is *conservative* for killing: if the coin's ensemble
already covers a measured value, a wider null covers it too, so every
claim [lem:coin] killed stays killed. It is the other direction that
fails. A value that escaped the coin has not been shown to escape the
right null, and this remark is an instance — [rem:cnclass]'s $z=9.01$
was a real escape from the coin at $32$ draws, is still an escape at
$512$, and is not an escape from the multiplicative ensemble at all.

One refinement to [rem:cnmultdeep], which read the closure as the
multiplicative ensemble's tail. That reading is right about why $32$
draws were not enough — the multiplicative spread more than doubles
from $32$ to $512$ while the iid spread moves about a tenth — but it
is not the whole cause. What covers $\mu$ is the multiplicative
structure, and the depth is what was needed to see that the structure
has a tail. Both were necessary; only one of them is the answer.

One band, one class family, one statistic. What this settles is not
[lem:coin] in general but that the question is worth asking of it,
which nothing before this had established.

#### Remark (counting what stands on a coin: eleven blocks, and the two that matter sit on eight draws) {#rem:coinsurface}
<!-- evidence: audit_coin_surface.py -->

[rem:cncoindeep] left a caveat with a fixed direction — a narrower
null is conservative for killing, so what [lem:coin] killed stays
killed, and only measurements recorded as *escaping* a coin are at
risk — and no count. A caveat touching three remarks is a footnote and
one touching thirty is a rewrite. This counts.

Of $237$ labelled blocks across the two documents, $62$ mention a coin
or an ensemble (H2) and **$11$ also contain one of ten phrases, fixed
before the run, that state an escape** (H3, REFUTED against a cap of
fewer than ten). $8$ of the $11$ are in this branch, labels beginning
`rem:cn` (H4).

**H1 caught the census before the census caught anything.** The
positive control `rem:cncoindeep` was missed on the first execution,
because the paper writes its numbers inside dollar signs and the
search read the markup literally. The fix is to the reader and not to
the rule: dollars are stripped, whitespace collapsed, and each phrase
matched at a word boundary so that "0 of" does not fire inside
"10 of". The phrase list is unchanged, and the disclosure sits in the
result file.

The count is a lower bound and every hit is a pointer to read, which
the script says it cannot do for itself. Reading the three outside
this branch:

| block | what caught it | reading |
|---|---|---|
| [rem:maskrivals] | "none of the three is separated" | a false hit — it says the rivals are *not* separated |
| [rem:identitynull] | "negative throughout where $0$ of $8$ draws are" | a real escape |
| [rem:decaynull] | "$0$ of $8$ draws fall monotonically at all" | a real escape |

So the surface outside the $C(N)$ branch is two blocks, not three, and
**both rest on ensembles of eight draws.** That is the number worth
carrying out of this run. [rem:cnmultdeep] measured a multiplicative
ensemble's spread more than doubling between $32$ draws and $512$;
eight is far below even the smaller of those. Neither block is
impeached here — the census cannot say whether their objects are
multiplicative or whether eight draws suffice for what they measure —
but both now have a specific question against them that did not exist
before, and it is the cheapest kind to answer.

H3's refutation is the honest headline: eleven is not a footnote. What
softens it is not the count but the reading. Eight of the eleven are
the chain this branch already worked through and corrected in place,
one is a false hit, and the two that remain need draws rather than a
rewrite.

#### Remark (the second branch says the same thing: the control was wrong, not small) {#rem:decaydeep}
<!-- evidence: audit_decayfamily_deep.py -->

[rem:coinsurface] narrowed the surface outside the $C(N)$ branch to two
blocks, both resting on eight draws. This asks the question of one of
them. [rem:decaynull]'s separation in size — $\mu$'s $|1/2-f|$ at
$0.2727$ and $0.2772$ against a coin maximum of $0.0826$ — is
re-measured at the two smallest $N$ with the field, weight, $k$-range
and $\theta'$ taken from that script unchanged, against $512$ iid
draws and $512$ random multiplicative ones. The gate reproduces both
published values exactly (J1).

| null | $N$ | $\mu$ | draw max | draws reaching $\mu$ |
|---|---|---|---|---|
| iid | $200000$ | $0.2727$ | $0.1213$ | $0$ |
| iid | $400000$ | $0.2772$ | $0.1197$ | $0$ |
| multiplicative | $200000$ | $0.2727$ | $0.4959$ | $67$ |
| multiplicative | $400000$ | $0.2772$ | $0.4937$ | $54$ |

**J2 holds and J3 is REFUTED, and the pair is the finding.** Depth
alone leaves the separation exactly where [rem:decaynull] left it —
none of $512$ iid draws comes near $\mu$ at either $N$. Against a
random multiplicative ensemble, about one draw in eight reaches it at
the smaller $N$ and one in ten at the larger. **Eight draws were not
hiding a tail; the control was the wrong kind.**

J4 holds and is resolved rather than nominal: the multiplicative
maximum is $0.4959$ against the iid $0.1213$, a gap of $0.3746$
against ensemble spreads of $0.1127$ and $0.0241$, so the widths
genuinely differ here by about a factor of four.

That is the shape [rem:cncoindeep] found for $C(N)$, now in a second
and independent place, and it does more here than there. In the $C(N)$
branch the iid coin's narrowness cost a reading that this repository
had already flagged as its own. **Here it flips a published verdict**:
[rem:decaynull] wrote that every draw sits an order below $\mu$ at
every $N$, and that is true of iid draws and false of multiplicative
ones. The sentence is flagged there.

Two things this does not touch, and both matter. It covers $N=2\cdot
10^5$ and $4\cdot10^5$ only — the statistic costs a strided sum for
every squarefree $k<N^{0.56}$ at every draw, and $512$ draws over the
published six-point sweep is hours rather than minutes — so the four
larger $N$ are untested. And [rem:decaynull] makes a second, separate
claim, that the $\alpha$-sweep discriminates: $\mu$'s interior
$\alpha^*=1.45$ at residual $0.011556$ against draws pinning at the
grid ends with residuals $2.991540$ to $10.078929$. **Nothing here
tests that**, and a multiplicative ensemble might behave differently
there too. It is the next thing to ask of this block, not something
this remark has answered.

#### Remark (both halves go the same way, and the cost objection was mine) {#rem:decaysweep}
<!-- evidence: audit_decayfamily_sweepdeep.py -->

[rem:decaydeep] flipped one of [rem:decaynull]'s two claims and said
what it had not done: two of six $N$, and nothing about the
$\alpha$-sweep, which is the sharper claim — $\mu$ gives an interior
minimiser at $\alpha^*=1.45$ with residual $0.011556$ while not one of
eight published draws has an interior minimiser at all. It also gave a
cost reason for stopping at two $N$, and **that reason was wrong.**
One draw over the whole six-point sweep takes under half a second: the
$k$-range is squarefree and coprime to $N$, which leaves $2186$ of
them at the top $N$ rather than the $6570$ that $N^{0.56}$ suggests.
So this run does all six $N$ and both claims, at $512$ draws of each
null, and the gate reproduces $\mu$'s published sweep exactly (K1).

**The size separation fails at every $N$, not only the two reached
before** (K2), and it gets worse as $N$ grows:

| $N$ | $\mu$ | iid draws reaching | multiplicative draws reaching |
|---|---|---|---|
| $200000$ | $0.2727$ | $0$ | $83$ |
| $800000$ | $0.2265$ | $0$ | $99$ |
| $3200000$ | $0.1793$ | $0$ | $125$ |
| $6400000$ | $0.1624$ | $0$ | $139$ |

**And the sweep goes the same way** (K3). Six multiplicative draws
have both an interior minimiser and a residual at or below $\mu$'s,
and the best residual any draw achieves is $0.000003$ against $\mu$'s
$0.011556$ — better by more than three orders. The iid ensemble does
not manage it once in $512$: twelve draws have an interior minimiser
and none of them gets under $\mu$'s residual, the best being
$0.089978$ (K4). **So the published claim survives depth and fails the
null type, exactly as the size claim did.**

Both halves of that block therefore rest on an iid coin, and
[rem:cncoindeep]'s pattern now holds in every place it has been looked
for: twice in the $C(N)$ branch and twice here.

One objection has to be raised against this run rather than left for a
reader. The multiplicative draws reach $0.4990$ to $0.5000$ in
$|1/2-f|$, which is the largest the statistic can be — a single draw
can make most $f(k)$ agree in sign, because the $k$ share small prime
factors and $f$ is multiplicative. A null that saturates a statistic's
range can be too wide to say anything. **It is not too wide here**, and
the counts are how one tells: only $83$ to $139$ of $512$ draws reach
$\mu$, so $\mu$ sits between roughly the seventy-third and
eighty-fourth percentile rather than in the bulk. The correlation that
makes the tail is the same correlation $\mu$ has, which is the whole
reason this is the right null and not merely a wider one.

What this does not settle. One seed per ensemble. And the counts
rising with $N$ is a trend across six points that nothing here tests
as a trend; it is printed because it points the other way from the
block's extrapolation, not because this run measured a slope.

#### Remark (the one that survives, and it survives for a reason that can be written down) {#rem:identitydeep}
<!-- evidence: audit_directidentity_deep.py -->

[rem:coinsurface] left three blocks outside the $C(N)$ branch standing
on eight draws; one was a false hit and [rem:decaydeep] and
[rem:decaysweep] took the second, both halves of which failed a random
multiplicative ensemble. This is the third and last, and it goes the
other way.

| null | $N=200000$ | $N=400000$ | $N=800000$ |
|---|---|---|---|
| iid, smallest of $512$ | $4.8387$ | $5.2286$ | $5.5530$ |
| multiplicative, smallest of $512$ | $1.9149$ | $2.0030$ | $2.1084$ |
| multiplicative, median | $4.2771$ | $4.5269$ | $4.7875$ |
| $\mu$ | $1.0039$ | $0.9865$ | $0.9893$ |

**Neither ensemble reaches $\mu$: $0$ of $512$ in both** (L2 holds, L3
REFUTED). The gate is exact — $|T-R|/R$ comes out $6.752\cdot10^{-15}$
and the three published ratios reproduce to four decimals (L1).

**And the two nulls are not the same, which is where the content is.**
Writing $v=N-i$, the statistic is
$T=\sum_v\Lambda(N-v)\,s(v)\sum_{k\mid v,\,k\ge2}\log k\;s(k)$, and a
random multiplicative $f$ has a cancellation an iid draw has none of:
$\sum_{d\mid n}f(d)=\prod_{p\mid n}(1+f(p))$ vanishes as soon as one
$f(p)=-1$. That cancellation is worth a factor of about two and a half
— it takes the best draw from $4.8387$ down to $1.9149$ (L4) — and it
is not worth the rest. Only $\mu$ has
$\sum_{k\mid v}\mu(k)\log k=-\Lambda(v)$, which collapses $T$ to $R$
exactly.

So there are three tiers and each is explained: independent signs give
about five, multiplicative signs about two at best, the Möbius
identity gives one. **This is the only thing in this repository that a
resolved multiplicative ensemble has failed to reach**, and every
other block [rem:coinsurface] named has now been tested.

Two things keep it in proportion. **The survival is not evidence of
hidden structure** — the identity is a theorem, so of course no other
sign function has it; what this run adds is the size of the gap that
the identity is responsible for, which was not known and is not
small. And the bound is the smallest ratio seen and not zero: $1.9149$
at the smallest $N$, one seed per ensemble, $512$ draws, and three of
the five $N$ that block publishes. L4 holds but not by much at the
smaller $N$ — the gap of $2.9238$ against a multiplicative spread of
$2.0400$ is a factor of $1.43$, which is a direction rather than a
measurement, as its rule said in advance.

#### Lemma (the coin control) {#lem:coin}
<!-- evidence: analytic -->

Let $\varepsilon(v)=\pm1$ be arbitrary signs on $\{v : \mu(v)\neq0\}$
and zero elsewhere. Then $\varepsilon^2 = \mu^2$ pointwise, so
$V(N)$ is unchanged, and the field
$C_\varepsilon(N)=\sum_v\varepsilon(v)\Lambda(N-v)$ has the same exact
second moment as $C(N)$ for every $N$. Consequently any estimator whose
output is reproduced when $\mu$ is replaced by $\varepsilon$ is not
measuring $\mu$.

**What this licenses, and what it does not.** The implication runs one
way: *reproduced under $\varepsilon$* $\Rightarrow$ *not measuring
$\mu$*. Its converse does not follow and is not proved anywhere here,
so a statistic that a coin ensemble fails to reach has not thereby
been shown to measure $\mu$ — it has been shown not to be killed by
this particular control. Four measurements now say how large that
distinction is.

[rem:cncoindeep] resolved an iid ensemble to $512$ draws beside a
random multiplicative one on the same statistic and found the iid
spread $0.6455$ of the multiplicative: **a coin is narrower than the
right null wherever the object is multiplicative**, so significance
read against one is larger than it should be. Because narrowness is
conservative for killing, every claim this lemma has killed stays
killed; it is the other direction that fails.

[rem:coinsurface] counted the blocks in the affected direction —
eleven contain a phrase claiming an escape, of which one is a false
hit — and the three outside the $C(N)$ branch have since been tested.
[rem:decaydeep] and [rem:decaysweep] took both halves of
[rem:decaynull]: each survives $512$ iid draws untouched and each
fails a multiplicative ensemble, $83$ to $139$ draws of $512$ reaching
$\mu$ in size and six matching it on the sweep. [rem:identitydeep]
took [rem:identitynull] and it survives both, with the reason written
down: no sign function but $\mu$ has
$\sum_{k\mid v}\mu(k)\log k=-\Lambda(v)$.

So the lemma is correct as stated and was read for more than it says.
**An estimator that escapes an iid coin has cleared the weakest
control in this paper**, and until it has been put against an ensemble
that shares $\mu$'s multiplicativity, the escape is a fact about the
control.


#### Remark (the seed is paid, and the tail is steadier than the bulk) {#rem:identityseeds}
<!-- evidence: audit_directidentity_seeds.py -->

[rem:identitydeep] found the one separation a resolved multiplicative
ensemble fails to reach and wrote its own debt: one seed, and another
owed before the result is leaned on. Two fresh seeds, $512$
multiplicative draws each, at the same three $N$.

| seed | $N=200000$ min | median | $N=800000$ min | median |
|---|---|---|---|---|
| first | $1.9149$ | $4.2771$ | $2.1084$ | $4.7875$ |
| $20260907$ | $2.0049$ | $4.4979$ | $2.1938$ | $5.0514$ |
| $20260908$ | $2.0631$ | $5.3443$ | $2.2752$ | $5.9581$ |

**M2 holds: no draw of $512$ at either fresh seed reaches $\mu$ at any
of the three $N$** — $0$ everywhere, as at the first seed. The debt is
paid and [rem:identitydeep]'s refutation is not one seed's luck.

**M4 is REFUTED and it broke the wrong way round.** It predicted the
bulk would agree more tightly than the tail; the median moves by as
much as $1.1706$ across seeds while the minimum moves by at most
$0.1702$. **The tail is the steady one, by a factor of seven.** M4's
own rule said a median moving that far would put M3 and the design in
question rather than just the tail — and read literally that is what
has happened, so it is recorded that way. But the tension resolves in
a direction the rule did not anticipate, and the resolution is worth
stating rather than hiding behind the verdict: M2 is a statement about
the *minimum*, and the minimum is precisely the quantity that
reproduced. M3 holds with room to spare, its gaps of $0.0854$ to
$0.1702$ sitting an order below the ensemble spreads of $2.0400$ to
$2.3376$.

What the design assumed, and got wrong, is that a heavy-tailed
ensemble's extremes are its least reliable part. Here they are its
most reliable. **That is an observation and not a measurement** — this
run did not ask why — but the shape it suggests is that the low end of
this ensemble is pinned by configurations in which many $f(p)=-1$ and
the divisor sums $\prod_{p\mid n}(1+f(p))$ vanish wholesale, which is
a structural floor rather than a fluctuation, while the bulk floats
with whatever the small primes happen to be.

So the standing claim narrows in one place and firms in another. It
firms because three seeds now agree that nothing reaches $\mu$. It
narrows because the bound quoted is the smallest of three minima and
those differ by $0.15$, so the honest bound is about $1.9$ and not
$1.9149$. Three of five $N$ remain untested, as before.

#### Remark (the last two N, and a floor the draws cannot get under) {#rem:identitybig}
<!-- evidence: audit_directidentity_bigN.py -->

[rem:identitydeep] and [rem:identityseeds] both ended on the same
sentence: the two larger $N$ of [rem:identitynull] are untested. This
tests them, and with them item 8's work list is finished.

| $N$ | $\mu$ | multiplicative min of $256$ | median | draws reaching $\mu$ |
|---|---|---|---|---|
| $1600000$ | $1.0017$ | $2.4598$ | $4.7259$ | $0$ |
| $3200000$ | $0.9968$ | $2.5636$ | $4.9806$ | $0$ |

**N2 holds where nothing had looked**, and N3 with it: across all five
$N$ the multiplicative minimum runs $1.9149$, $2.0030$, $2.1084$,
$2.4598$, $2.5636$ — monotone, so **the separation widens with $N$**
rather than closing. $\mu$ meanwhile sits at $1.00$ throughout.

**N4 is the one worth pausing on.** It was put in to price a shortcut:
these two $N$ cost about ten seconds a draw against under two at the
smaller ones, so this run used $256$ draws where the others used
$512$, and a shallower ensemble gives a higher minimum for free. N4
asked how much of N3's rise that buys. The answer is none that can be
seen: **the minimum over the first $64$ draws, over $128$, and over
$256$ is the same to four decimals — a drift of $0.0000$.** One draw
inside the first sixty-four reached $2.4598$ and the next hundred and
ninety-two never improved on it.

That is the second time this ensemble's low end has behaved like a
floor rather than a tail. [rem:identityseeds] found the minimum moving
by at most $0.1702$ across three seeds while the median moved by
$1.1706$; here it does not move with depth at all. The two together
say the same thing from different directions, and the reading offered
there — that configurations with many $f(p)=-1$ make the divisor sums
$\prod_{p\mid n}(1+f(p))$ vanish wholesale and pin the low end — now
has a second observation behind it. **It is still an observation**:
nothing here has measured the mechanism, only its signature.

The practical consequence is that the bias N3's rule named against
itself is not there. A rise from $2.1084$ to $2.4598$ across a factor
of two in $N$ is not a depth artefact, because depth does not move
this quantity over a fourfold range. **N3 holds for the reason it
intended.**

So the identity separation holds across the whole range that block
publishes, at three seeds where it was measured deeply and at one seed
here, and the honest bound is the smallest ratio seen at each $N$.
G19 caught two of those minima typed by hand into a sentence rather
than read from the file that produced them; they are read now.

#### Remark (the mechanism is measured, the word "floor" was wrong, and the bound is looser than it looked) {#rem:identityfloor}
<!-- evidence: audit_identity_floor.py -->

[rem:identityseeds] and [rem:identitybig] both saw the multiplicative
ensemble's low end behave oddly and both offered the same reading as
an observation. It is not a guess. For multiplicative $s$,
$\log k=\sum_{p\mid k}\log p$ gives

$$
W(v)=\sum_{k\mid v}s(k)\log k
 =\sum_{p\mid v}\log(p)\,s(p)\!\!\prod_{q\mid v,\,q\neq p}\!\!(1+s(q)),
$$

so **if two primes dividing $v$ carry $s(q)=-1$ every term has a zero
factor and $W(v)=0$**: $W$ lives on the $v$ with at most one negative
prime. Scoring each draw by $m$, how many of the ten smallest primes
it sends to $-1$, tests whether that is what drives the measured
minimum. The draws are the same draws — same seed, sieve top and pass
order as [rem:identitydeep], and the minima reproduce at $1.9149$ and
$2.0030$ (P1).

| $m$ | draws | ratio min | ratio median |
|---|---|---|---|
| $1$ | $6$ | $8.3324$ | $9.4563$ |
| $4$ | $104$ | $2.5903$ | $5.2067$ |
| $7$ | $60$ | $1.9149$ | $3.0533$ |
| $9$ | $5$ | $2.0899$ | $2.2391$ |

**The mechanism is real and it dominates** (P2): the slope of
$\log(\text{ratio})$ on $m$ is $-0.15802\pm0.00925$, $t=-17.08$, and
the median falls monotonically from $9.4563$ to $2.2391$ across the
range. The minimising draw has $m=7$, in the top decile (P3).

**P4 is REFUTED and the sentence this run had waiting for that is
wrong.** P4 asked whether the top decile of $m$ is tight — at most
half the spread of the bottom decile — and got $1.2152$ against
$1.8725$, a ratio of $0.6490$. The message written in advance for that
outcome said the low end is a tail and the earlier reading was wrong.
It does not follow, and P2 and P3 are why: the ratio is driven by $m$
with a $t$ of seventeen, and the minimum is where $m$ is largest. What
P4 actually shows is that the contraction with $m$ is **gradual**. The
top decile is $m\ge7$ and is mostly $m=7$, where the ratio still runs
$1.9149$ to $6.7781$; the tightness appears only at $m=9$, where five
draws span $2.0899$ to $2.7002$. The cap was set without knowing that
granularity.

**And the run found something neither remark suspected.** The $m$
counts are $6,27,59,104,137,94,60,20,5$ — no draw reached $m=10$.
Ten independent signs give $m=10$ with probability $1/1024$, so $512$
draws are expected to miss it, and they did. **The observed minimum is
therefore the best of $m\le9$ and not the ensemble's infimum.** What a
draw with every one of the ten smallest primes negative would give is
not measured here, and it is the natural place for the bound to move.
That does not touch [rem:identitydeep]'s refutation — no draw came
within a factor of $1.9$ of $\mu$ — but it does mean the number
$1.9149$ is a sample minimum in a stronger sense than "smallest of
$512$": it is the smallest of $512$ draws none of which realised the
configuration the mechanism says is best.

#### Remark (forced, the mechanism's best still stops short — and a prediction that could be omitted) {#rem:identityforced}
<!-- evidence: audit_identity_forced.py -->

[rem:identityfloor] found that no draw of $512$ reached $m=10$, so
the observed minimum was the best of $m\le9$ and the configuration the
mechanism calls best had never been seen. This run forces it: the
first $J$ primes are set to $-1$ and the rest left iid, for
$J=0,5,\dots,30$, $64$ draws each. The construction reproduces $\mu$'s
own ratios $1.0039$ and $0.9865$ exactly (Q1).

| $J$ | forced to | min | median |
|---|---|---|---|
| $0$ | — | $2.4168$ | $4.5493$ |
| $10$ | $29$ | $1.8018$ | $2.1264$ |
| $20$ | $71$ | $1.6716$ | $1.8172$ |
| $30$ | $113$ | $1.6100$ | $1.6894$ |

The ratio falls with $J$ at every step (Q2). **But it saturates, and
not at $\mu$** (Q3, refuted). The medians run $4.5493$, $2.6273$,
$2.1264$, $1.9401$, $1.8172$, $1.7391$, $1.6894$ — increments
$1.92,\,0.50,\,0.19,\,0.12,\,0.08,\,0.05$, a sequence settling near
$1.6$ while $\mu$ sits at $1.0039$. The closest any forced draw comes
is $1.6100$.

**The prediction registered here was that forcing would reach $\mu$**,
on the algebra that $W$ dies wherever two primes dividing $v$ are
negative, so emptying the support should drive $|T|$ down. The support
does empty and the ratio does fall; it stops. Emptying $W$'s support
removes mass from the numerator and the normalisation together, and
past $J\approx20$ the two move at the same rate. So the reasoning was
right about the direction and wrong about the destination, and the
result is the *stronger* one for [rem:identitydeep]: the separation
between $\mu$ and a multiplicative sign function is not a matter of
how deep the ensemble is sampled. It survives the mechanism's own best
configuration, forced rather than waited for. **Bounded at $J\le30$,
which is a choice this run made and not a limit it established.**

**And a registered prediction turned out not to be one** (Q4). Q4
priced the rarity of "the smallest $J$ that reaches $\mu$", and its
rule said it would be *omitted* if Q3 fell. Q3 fell, the object is
empty, and the gate refused a pre-registered item with no verdict —
correctly. Q4 is refuted, not omitted: it asserted an existence and
the existence fails. The rule is not rewritten here. What is recorded
is that **a prediction whose rule allows it to be omitted is not a
prediction**, which is a defect of the same family as M9 and its
fourth instance in this repository.

#### Lemma (the placebo key) {#lem:placebo}
<!-- evidence: lab_mask_placebo.py -->

Let the cells be indexed by a labelling $\ell(N)$, and let $\pi$ be a
permutation of the even integers in the band. Replacing $\ell$ by
$\ell\circ\pi$ preserves every cell size and leaves the field
$Z(N)=C(N)/\sqrt{V(N)}$ byte-identical, while destroying the
correspondence between cells and arithmetic. Any cell-indexed statistic
that survives this replacement is a property of the cell sizes and not
of that correspondence.


These are trivial and they are load-bearing. Lemma [lem:coin]
invalidates the natural measurement of $\rho$: the centred estimator
returns the same value on $\varepsilon$ as on $\mu$, so a quoted level
for $\rho$ calibrates nothing. **Any claim about
$\mathrm{Var}\,C$ relative to $V$ must first exhibit an estimator that
distinguishes $\mu$ from a coin, and show that it does.**
Lemma [lem:placebo] is what confirms, in Section [sec:floor], that
the cell floor is a property of the arithmetic and not of the
partition.


#### Proposition (the coin obstruction is arithmetic, not informational) {#prop:coindisc}
<!-- evidence: lab_coin_discriminator.py -->

An earlier version added that with one realisation of $\mu$ such an
estimator "may not be possible at all, and saying so would itself be a
result". It is possible, on both counts:

**(i)** Estimators that separate $\mu$ from a coin are cheap,
provided they read *multiplicativity* rather than variability.
$\mu(2v)=-\mu(v)$ for odd $v$ is an exact identity, while
$\varepsilon(2v)$ and $\varepsilon(v)$ are independent signs, so

$$
T(x) \;=\; \frac1x\sum_{v\le x}\mu(v)\mu(2v)
  \;\longrightarrow\; -\frac{4}{\pi^2}
$$

— minus the density of odd squarefree integers — while its coin
analogue is $O(x^{-1/2})$. Measured, $T(x)$ reads
$-0.405600,\,-0.405270,\,-0.405286,\,-0.405295$ at
$x=10^4,10^5,10^6,2\cdot10^6$ against $-4/\pi^2=-0.405285$, while the
worst of twenty coin draws on the same support is
$0.020600,\,0.004570,\,0.001214,\,0.000935$. At the top the separation
is about $430$ coin standard deviations.

**(ii)** The field itself determines $\mu$. Writing $M=N-2$ and
$\widehat\Lambda(m)=\Lambda(m+2)$, one has
$C(N)=(\widehat\Lambda*\mu)(M)$ as an additive convolution with
$\widehat\Lambda(0)=\Lambda(2)=\log2\neq0$, so $\widehat\Lambda$ is
invertible and $\mu(M)=\sum_{m<M}a(m)\,C(M-m+2)$ for the inverse
filter $a$. Lemma [lem:coin] is therefore not an
information-theoretic obstruction: two different sign patterns give
two different fields.




#### Remark (but the inverse filter is numerically dead) {#rem:filter}
<!-- evidence: lab_coin_discriminator.py -->

What kills the constructive route is amplification, not information.
The inverse filter grows geometrically: $|a(m)|$ runs
$1.443,\,2.287,\,2.182,\,20.91,\,512.3,\,3.443\cdot10^{5},\,
1.073\cdot10^{14},\,1.535\cdot10^{28},\,3.146\cdot10^{56}$ at
$m=0,1,2,5,10,20,50,100,200$, with
$|a(200)|^{1/200}=1.916413$. The recursion itself is exact — the
residual of $a*\widehat\Lambda$ against $\delta_0$ is
$7.018\cdot10^{-17}$ relative to the mass of the sum, i.e. machine
epsilon — but in double precision the absolute residual passes
$10^{-9}$ already at $m=31$, where $|a(m)|=4.479\cdot10^{8}$. That
absolute tolerance was the audit's rule H1 and it fails; the failure
is the amplification and not the recursion, which is why the relative
figure is the one quoted. A rule stated as an absolute tolerance on a
sum whose terms reach $10^{56}$ tests double precision, not the
mathematics — the same mis-specification as
Remark [rem:cap], and the third of its kind in this work.

The reconstruction behaves accordingly: it returns $\mu(M)$ to
$1.5\cdot10^{-5}$ at $M=40$ and is wrong by $4$ at $M=60$ and by
$1.1\cdot10^{12}$ at $M=100$. So the field carries $\mu$ exactly and
releases it only through a filter that amplifies by $1.92^{m}$.

That is the sharper form of Lemma [lem:coin], and it is the one to
apply to a new approach: the pass condition is not "exhibit any
discriminating estimator" — those exist — but "exhibit one whose
gain against the coin is not paid for by an amplification of the same
order".


#### Proposition (the demand-side discrepancy is a dilate) {#prop:dilate}
<!-- evidence: lab_dilate_identity.py -->

Write $A(N;k)=\sum_{n<N,\ n\equiv N\,(k)}\Lambda(n)\mu(N-n)$ for the
progression sum inside $\Emu(N;k)$. Then, unconditionally,

$$
\begin{equation}\label{eq:dilate}
  A(N;k) \;=\; \mu(k)\,H(N;k),
  \qquad
  H(N;k) \;=\; \sum_{\substack{m<N/k\\ (m,k)=1}}\Lambda(N-mk)\,\mu(m).
\end{equation}
$$




**Proof.** 
The class $n\equiv N\pmod k$ is $k\mid N-n$, so writing $N-n=mk$ turns
the sum into $\sum_{m<N/k}\Lambda(N-mk)\mu(mk)$; and
$\mu(mk)=\mu(m)\mu(k)$ when $(m,k)=1$, zero otherwise.
 ∎




Two things follow. First, the pass condition above is met, by the
ratio $r_s(N;k)=A_s(N;k)/\bigl(s(k)H_s(N;k)\bigr)$: it is identically
$1$ for $s=\mu$ and nothing of the kind for a coin, since
$\varepsilon(mk)$ is independent of $\varepsilon(m)$ and
$\varepsilon(k)$. Measured over every squarefree $k<N^{0.56}$ coprime
to $N$ at $N=2\cdot10^5,\,4\cdot10^5,\,8\cdot10^5$, the worst relative
error of [eq:dilate] is $2.294\cdot10^{-14}$,
$1.197\cdot10^{-13}$, $4.591\cdot10^{-14}$; the fraction of $k$ at
which the ratio is $1$ is $1.000000$ for $\mu$ and $0.000000$ for
every one of eight coin draws, whose median $|r_\varepsilon-1|$ is
$1.4366$. **Both sides of the ratio are the same size, so the gain
costs no amplification at all** — unlike the inverse filter of
Remark [rem:filter], which pays $1.92^{m}$.

Second, [eq:dilate] qualifies the organising principle of
\S[sec:claims]. The demand side is divisibility-coupled where
$EH_\mu$ is *consumed*, but its per-modulus discrepancy is
$\mu(k)$ times a *difference-coupled* sum at scale $N/k$ — the
supply side's own object, one dilation down. The two couplings are not
two fields; they are the same field seen before and after the switch.


#### Remark (why the coin beat $\mu$ on the level) {#rem:whycoinwins}
<!-- evidence: lab_dilate_identity.py -->

Remark [rem:levelmeas] withdrew a level measurement because a coin
reached a higher $K^*$ than $\mu$ at every $N$, and left the reason
open. [eq:dilate] supplies it. For $\mu$,
$|A(N;k)|=|H(N;k)|$ is a Möbius–prime correlation of length $N/k$
— the wall, dilated — and nothing makes it small. For a coin,
$\sum_m\Lambda(N-mk)\varepsilon(mk)$ is a sum of $N/k$ independent
signs, and square-root cancellation is both the best it can do and
what it does. Measured, $\langle|A_\mu|\rangle/
\langle|A_\varepsilon|\rangle$ reads

$$
\begin{array}{r|ccc}
 k\text{-band} & [2,32) & [32,256) & [256,2048)\\\hline
 N=2\cdot10^5 & 1.0761 & 1.5306 & 1.4820\\
 N=4\cdot10^5 & 1.2570 & 1.2671 & 1.4496\\
 N=8\cdot10^5 & 1.1733 & 1.5013 & 1.4830
\end{array}
$$

So $\mu$ is the noisier of the two in progressions, by half again at
the larger moduli, and the withdrawal of Remark [rem:levelmeas] has
a structural cause rather than merely a failed null. It also sharpens
what the demand asks: $B(N)=\sum_k(\log k)|H(N;k)|$ is a weighted sum
of dilated walls, so Proposition [prop:nolog]'s constant-factor
bound is a constant-factor bound on the wall itself, summed over
dilations.


#### Proposition (the weights are nonnegative) {#prop:posweights}
<!-- evidence: lab_positive_weights.py -->

Substituting [eq:dilate] into $E_3$ and using $\mu(k)^2=1$ on the
squarefree $k$,

$$
\begin{equation}\label{eq:posweights}
  E_3(\alpha) \;=\; \sum_{\substack{k<K\\(k,N)=1}}
    \mu^2(k)\,(\log k)\,H(N;k) \;-\; C(N)\,B_{\log}(K).
\end{equation}
$$

**The signs $\mu(k)$ cancel exactly.** They cancel against the
$\mu(k)$ that [eq:dilate] extracts from the progression sum, and what
is left is a sum of dilated walls whose weights $\mu^2(k)\log k$ are
nonnegative. Verified to machine zero — the relative error of
[eq:posweights] against a direct evaluation of $E_3$ is at most
$1.8\cdot10^{-16}$ over $N=2\cdot10^5$ to $3.2\cdot10^6$.

This is worth pausing on. Huang–Li discard the $\mu(k)$ by the
triangle inequality, and \S[sec:demand] says that keeping them is
what the results here exploit. Keeping them is what
Theorem [thm:A] exploits, on the $w_k=1$ branch. On the $\log k$
branch they are not there to keep.


#### Remark (there is no cancellation across dilations) {#rem:nocrossk}
<!-- evidence: lab_positive_weights.py -->

Nonnegative weights mean that whatever smallness $E_3$ has must come
from the signs of $H(N;k)$ across $k$. It does not come from there.
Write $G=\sum_k(\log k)|H|\big/\bigl|\sum_k(\log k)H\bigr|$ for the
gain that cancellation across dilations buys:

$$
\begin{array}{r|ccccc}
 N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6
   & 3.2\cdot10^6\\\hline
 \#k & 313 & 462 & 682 & 1004 & 1485\\
 G & 1.834 & 1.804 & 2.207 & 2.588 & 2.789\\
 \sqrt{\#k} & 17.7 & 21.5 & 26.1 & 31.7 & 38.5
\end{array}
$$

Independent signs would give $\sqrt{\#k}$. The measured gain is $2$ to
$3$. Equivalently the sum behaves as though it had
$n_{\mathrm{eff}}=G^2=3.4,\,3.3,\,4.9,\,6.7,\,7.8$ independent signs,
where it has hundreds to thousands of terms: **the dilated walls move
together.** For a coin, which keeps signed weights because
$\varepsilon(k)\varepsilon(mk)$ does not collapse, the same gain reads
$965.6,\,13.3,\,25.0,\,1000.4,\,96.0$ — erratic, since its
denominator is a near-cancellation, but always more than five times
$\mu$'s.

The audit's rule T4 guessed the wrong mechanism, and its failure is
informative: it proposed that a heavy tail explains the missing
$\sqrt{\#k}$, but the top decile of $k$ carries only
$0.3486$ to $0.3587$ of $\sum(\log k)|H|$ — concentrated, not
dominant — and the $(\log k)$-weighted fraction of $k$ with $H>0$ is
$0.3763$ to $0.4809$, near balanced. Neither domination nor a sign
bias: a correlation.

**What this costs the one-sided route, at accessible $N$.**
Proposition [prop:onesided] weakened the demand from two-sided to
one-sided, and the value of that depends on how much smaller a signed
sum can be than its absolute counterpart. Over the range measured here
the answer is a factor of about two, so [eq:onesided] and [eq:nolog]
sit within a small constant of each other. Remark [rem:leandecay]
shows that constant is not a constant — it grows — so this paragraph
describes the accessible range and not the asymptotics.


#### Remark (a Mertens reduction, closed) {#rem:mertens}
<!-- evidence: lab_mertens_reduction.py -->

Why the $H(N;k)$ move together has an obvious candidate answer, and it
is wrong. If the weight $\Lambda(N-mk)$ behaved like its mean over
$m$ — which is what the prime number theorem in progressions gives on
average — then $H(N;k)$ would be a constant times the
coprimality-restricted Mertens sum
$M_k(M)=\sum_{m<M,(m,k)=1}\mu(m)$, whose values at different arguments
are correlated over long ranges. The demand
$B(N)\le\SS(N)(1-A(N))N$ would then be a statement about Mertens sums
over dilations, a far more classical object.

All four pre-registered rules fail. The correlation
$\mathrm{corr}(H,M_k)$ across $k$ reads
$0.5066,\,0.5473,\,0.5411,\,0.3242,\,0.2724$ at
$N=2\cdot10^5$ to $3.2\cdot10^6$ — moderate, and *falling* with
$N$, where U1 asked for $0.7$. The slope of $H$ against $M_k$ is
$6.5$ to $11.2$ where U3 asked for $1$, and
$\sum(\log k)|M_k|$ is $0.058$ to $0.073$ of $\sum(\log k)|H|$ where
U4 asked for agreement to $30\%$: $|H|$ runs $14$ to $17$ times
$|M_k|$, against the $\sqrt{\log N}=3.5$ to $3.9$ that a
sparse-weight heuristic would give.

The control settles what is left. U2 predicted that a coin would
correlate just as well, so that the reduction would be a fact about
the weight; measured, the coin gives $0.25$ to $0.33$, and at the
largest $N$ $\mu$'s own $0.2724$ is *below* the coin's $0.3151$. So
the reduction fails for both and there is nothing $\mu$-specific in
what little correlation there is.

Decomposing $H=cM_k+R$, the residual carries
$0.7434,\,0.7004,\,0.7072,\,0.8949,\,0.9258$ of the variance, and its
share grows with $N$. **The $\Lambda$ weight is not a constant
multiplier; it supplies the bulk of $H$'s own fluctuation.** The line
is closed: the correlation across $k$ that Remark [rem:nocrossk]
found is not a common Mertens factor, and its mechanism stays open.


#### Remark (the correlation sits in the large terms) {#rem:signmass}
<!-- evidence: lab_sign_structure.py -->

Two explanations of Remark [rem:nocrossk] remained. Either a few $k$
carry the mass, so the sum fails to cancel for arithmetic reasons with
no correlation at all; or the signs of the large terms agree while the
small ones are balanced. They are separated exactly by the null that
keeps the magnitudes $|H(N;k)|$ and re-signs them at random, and it is
the second:

$$
\begin{array}{r|ccccc}
 N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6
   & 3.2\cdot10^6\\\hline
 G & 1.834 & 1.804 & 2.207 & 2.588 & 2.789\\
 G\ \text{re-signed} & 49.01 & 67.21 & 66.59 & 82.14 & 339.03\\
 \text{mass fraction } H>0 & 0.2273 & 0.2228 & 0.2735 & 0.3068
   & 0.3207\\
 \text{count fraction } H>0 & 0.4121 & 0.3831 & 0.4296 & 0.4392
   & 0.4808\\
 \text{low half} & 0.5096 & 0.4762 & 0.5396 & 0.5000 & 0.5774\\
 \text{high half} & 0.3141 & 0.2900 & 0.3196 & 0.3785 & 0.3841
\end{array}
$$

Re-signing buys $27$ to $122$ times more cancellation, so
heterogeneity is not the explanation. Counting moduli the signs are
near balanced; weighting by contribution, $68\%$ to $78\%$ of
$\sum(\log k)|H|$ is negative; and split at the median of $|H|$, the
small dilated walls are balanced while the large ones lean negative.
**That is why $E_3$ is negative at every $N$ measured and why the sum
does not cancel.**

The control is less clean than the null and is reported as it came.
Rebuilding the whole field from coin signs, ten draws give mass
fractions $0.3033,\,0.4977,\,0.5222,\,0.5012,\,0.5165,\,0.3773,\,
0.4544,\,0.5254,\,0.5085,\,0.4296$ — scattered about $\tfrac12$, one
of ten below $0.33$ — and half-split gaps that never reach $0.13$ in
absolute value, where every one of $\mu$'s five exceeds $0.12$ and the
largest is above $0.21$. So **part of the half-split is
mechanical**: conditioning on whatever sum was realised makes its
large terms share its sign, and a coin shows that too. What $\mu$ adds
is a stronger gap and a direction that is the same at all five $N$,
where the coin's sign is a fresh flip each draw. That said, the five
$N$ share one realisation of $\mu$, so the consistency of direction is
one observation and not five — the same caution
Remark [rem:mcratios] records for a different statistic.

Why the large dilated walls lean negative is answered next.


#### Remark (the lean is not the Mertens function) {#rem:leanmertens}
<!-- evidence: lab_lean_mechanism.py -->

Remark [rem:leandecay] calls the sign lean a finite-$N$ effect, which
without a mechanism is a label. [eq:dilate] offers one to test:
$H(N;k)$ sums $m$ over $1\le m<N/k$, so large $k$ means a *short*
inner sum, and $\sum_{m\le M}\mu(m)=M(M)$ is negative at most small
$M$. If $\Lambda(N-mk)$ were flat in $m$ the sign of $H(N;k)$ would be
the sign of $M(\lfloor N/k\rfloor)$.

**All four pre-registered rules fail.** The correspondence is not
there: over $2\le N/k\le10^3$ the agreement between
$\operatorname{sign}H(N;k)$ and $\operatorname{sign}M(\lfloor
N/k\rfloor)$ reads $0.5472,\,0.5483,\,0.5513,\,0.5522,\,0.5535$, while
shuffling the Mertens signs among the distinct values of $\lfloor
N/k\rfloor$ — which preserves both marginal sign distributions
exactly — reaches $0.6607,\,0.5610,\,0.6343,\,0.6297,\,0.6371$ at its
best draw. **P3** therefore fails at every $N$: the real alignment
does no better than a shuffle, and usually worse. **P4** fails too,
and in the opposite direction to the prediction: $|1/2-f_+|$ is
$0.1446$ where $N/k<10^3$ and $0.1667$ where $N/k\ge10^3$ at the
largest $N$, so the lean is not concentrated in the short sums at all.
The premise was only half true in any case — the fraction of
$2\le x\le10^4$ with $M(x)<0$ is $0.5630$, not one.

What the run does establish is a profile nobody had looked at. Binning
$k$ by the length $N/k$ of its inner sum, the mass-weighted $f_+$ runs

$$
\begin{array}{lccccc}
 N/k & [2,4) & [8,16) & [16,32) & [256,2^{10}) & [2^{12},\infty)\\\hline
 f_+ & 0.7325 & 0.1897 & 0.0936 & 0.2196 & 0.3849
\end{array}
$$

at $N=3.2\cdot10^6$ — **positive at the shortest sums, a sharp minimum
at $N/k\in[16,32)$, and a slow return toward $1/2$**. That is not
monotone, so **P1** and **P2** fail as stated, and the minimum sits at
the same octave at every $N$, drifting only from $0.0769$ to $0.0936$
across a factor $16$ in $N$.

The positive end is elementary and carries no information: $N$ is even,
so $N-mk$ is even whenever $mk$ is, and $\Lambda$ vanishes on even
numbers above $2$; for odd $k$ only odd $m$ survive, leaving the
nonnegative $m=1$ term almost alone. The minimum is not explained
here. What is closed is the Mertens reading of Remark
[rem:signmass]'s lean, and the closure is worth the cycle: it was the
only elementary mechanism on offer.


#### Remark (it is the Mertens function the parity leaves) {#rem:leanodd}
<!-- evidence: lab_lean_oddmertens.py -->

Remark [rem:leanmertens]'s own diagnostic supplies the repair. $N$ is
even, so $\Lambda(N-mk)$ vanishes whenever $mk$ is even; every *even*
$m$ contributes nothing to $H(N;k)$, while $M(x)=\sum_{m\le x}\mu(m)$
counts them. The predictor was summing over a set half of which cannot
appear. The object the parity leaves standing is
$M_{\rm odd}(x):=\sum_{m\le x,\ m\ \text{odd}}\mu(m)$.

With that one change the correspondence appears. Over
$2\le N/k\le10^3$ the agreement between $\operatorname{sign}H(N;k)$ and
$\operatorname{sign}M_{\rm odd}(\lfloor N/k\rfloor)$ is $0.7657$ at
$N=3.2\cdot10^6$, against $0.5738$ for the best of sixteen
permutations of the predictor's signs and $0.5563$ for the full
Mertens function — the same test, the same null, and the parity
restriction is the whole difference. Requiring $m$ coprime to $k$ as
well, which is the exact summation range of [eq:dilate], changes the
agreement by $+0.0000$ at every $N$: the coprimality is redundant
once the parity is imposed.

**So the sign lean of the dilated walls is $M_{\rm odd}$ read at the
length of the inner sum**, and Remark [rem:leandecay]'s "finite-$N$
effect" has a name. It also explains the decay without appeal:
$M_{\rm odd}(x)/\#\{m\le x\ \text{odd squarefree}\}$ runs
$-0.3333,\,-0.5610,\,-0.4756,\,-0.3407,\,-0.2262,\,-0.1352,\,-0.0535,
\,-0.0175$ across the octaves from $[4,8)$ up, so a fixed $k$-range
sees a shallower lean as $N$ grows and $N/k$ with it.

The fourth pre-registered rule fails and the failure is instructive.
**Q4** asked the predictor's deepest octave to be $[16,32)$, where the
measurement's minimum sits; it is $[8,16)$. The model behind the
prediction — that a contributing $k$ has one surviving $m$, so
$\operatorname{sign}H=\mu(m)$ and $f_+\approx\tfrac12(1+M_{\rm
odd}/{\rm count})$ — tracks the shape but under-predicts the depth,
by $-0.0298$ at $[8,16)$ and $-0.1686$ at $[16,32)$, widening to
$-0.2699$ at $[64,128)$. The sign of $H$ term by term is $M_{\rm
odd}$; how the surviving terms are *selected* is not one term drawn
uniformly, and that is what the miss measures. Q1 to Q3 are about the
first and hold; Q4 was about the second and does not.


#### Remark (the weights are free; the survivor set is everything) {#rem:levelweighted}
<!-- evidence: audit_level_weighted.py -->

Remark [rem:leveldemand] leaves $0.2271$–$0.2525$ of the demand
standing at $\alpha=1/2$ and OPEN asked whether that is a floor. It is
not, and the sieve family has no room above $1/2$ either: past
$\sqrt N$ nothing composite is left to remove and raising $Q$ only
strikes the small primes, which Remark [rem:sievedepth] measured as a
fall. So there is no $0.5<\alpha\le0.56$ inside this family.

**The $0.23$ is the cost of counting.** At $\alpha=1/2$ the survivors
are the primes, so $\log(N-mk)$ *is* $\log p$ on them — the
log-weighted predictor that Remark [rem:logweightpredictor] found
unhelpful at level $29$ becomes essentially exact at level $\sqrt N$.
**V1, V2, V3 and V4 all hold.** The log-weighted residual share at
$\alpha=1/2$ is

$$
0.0090,\ 0.0079,\ 0.0070,\ 0.0080,\ 0.0064 ,
$$

against $0.2525$ down to $0.2271$ unweighted — a factor of about
$28$. And what is left is identified rather than tolerated: the sieve
strikes each small prime along with its multiples, so the terms whose
$N-mk$ is a prime at or below $\sqrt N$ are true contributors it
removes, and their share of the demand measured directly is
$0.009755$ to $0.007108$ — the same interval.

**So the whole difficulty is the survivor set and none of it is the
weight**, which reverses how Remark [rem:logweightpredictor]'s failure
reads: the weight was useless there because the survivor set was
wrong, not because weights do not matter.

**And V4 is the part that costs.** At $\alpha=0.3$, the level at which
the sign stops slipping, the free weights take the share only from
$0.4866$–$0.6022$ to $0.4592$–$0.5858$ — a factor $0.93$ to $0.97$.
The demand's threshold does not move at all: it stays at
$\alpha=0.5$, well above the sign's $0.3$. What the weights change is
the size of what is left at the top, from $0.2525$ to $0.008981$.


#### Remark (holding the sign is not cutting the demand) {#rem:leveldemand}
<!-- evidence: audit_level_demand.py -->

Remark [rem:levelthreshold] put the level at which the wall's sign
stops slipping between $N^{0.2}$ and $N^{0.3}$, below the
$\theta'=0.56$ the reduction consumes, and left the question of
whether the shortfall there is usable. The currency is not the lean
but the demand: Remark [rem:predictable]'s residual share
$\sum(\log k)|H-\beta P_Q|\big/\sum(\log k)|H|$, which at level $29$
is $0.6310$ down to $0.5307$. Its far end is fixed by arithmetic —
at $\alpha=1/2$ the survivors are the primes, so $H$ and $P_{1/2}$ run
over the same terms and differ only by the weights $\log p$.

**U1 holds** — $\beta$ and all five shares reproduce exactly — and
**U3 holds**: at $\alpha=1/2$ the share is $0.2271$ to $0.2525$, so
the spread of $\log p$ alone accounts for under a quarter of the
demand and the rest is reachable in principle.

**U4 is refuted, and that is the answer.** At $\alpha=0.3$ the share
is $0.6022,\,0.5526,\,0.5147,\,0.4986,\,0.4866$ — against $0.6310$
down to $0.5307$ at the fixed level $29$. The level that holds the
sign takes almost nothing off the demand:

$$
\frac{\text{share}(29)-\text{share}(0.3)}
     {\text{share}(29)-\text{share}(0.5)}
= 0.0760,\,0.1288,\,0.1261,\,0.1383,\,0.1475 .
$$

**Eight to fifteen per cent of what an unbounded level would buy.**
Taken at matching $N$, the affordable level still leaves $2.10$ to
$2.39$ times what primality leaves. So the shortfall of Remark
[rem:levelthreshold] is the residue: the sign is held cheaply and the
demand is not.

**U2 is refuted by the ladder, not the data.** The fixed $Q=29$ is a
deeper sieve than $N^{0.1}$ and $N^{0.2}$ at every $N$ here — $29$
against $3$–$4$ and $11$–$20$ — so the first steps of the ladder go
backwards in depth. On the $\alpha$ ladder proper the share is
monotone at every $N$.

The two thresholds on this one axis are therefore different, and the
cheap one may not be carried to the other statistic: the sign's lean
stops slipping at $\alpha=0.3$, the demand halves only at
$\alpha=0.5$.


#### Remark (the level at which the sign stops slipping) {#rem:levelthreshold}
<!-- evidence: audit_level_threshold.py -->

Remark [rem:sievedepth] leaves the sign of the wall between two
extremes: held exactly at level $\sqrt N$, and degrading at the fixed
level $29$. Between them lies the only quantity this programme trades
in — the level as a power of $N$ — so sieving to $Q=N^\alpha$ and
asking where the degradation stops turns "not bounded" into a number
to set beside $\theta'=0.56$.

**T1 holds**, reproducing the top rung exactly. **T3 holds**:
$\alpha=0.4$ is flat, $+0.001564$ at $0.34$ standard errors.

**T2 and T4 are refuted, and both by want of power rather than by
direction.** Every level's agreement slope is negative and **none**
reaches two standard errors, the fixed level included at $1.10$. Six
points of a ratio of counts do not carry that much. The lean does, and
it is the statistic the route cares about:

$$
\begin{array}{l|ccc}
\text{level} & \text{slope} & \text{s.e.} & t\\\hline
Q=29 & -0.048145 & 0.013664 & 3.52\\
Q=N^{0.1} & -0.088359 & 0.012460 & 7.09\\
Q=N^{0.2} & -0.066398 & 0.015393 & 4.31\\
Q=N^{0.3} & -0.015265 & 0.010024 & 1.52\\
Q=N^{0.4} & +0.001937 & 0.008989 & 0.22\\
Q=N^{0.5} & +0.001303 & 0.000546 & 2.39
\end{array}
$$

**So the slipping stops between $\alpha=0.2$ and $\alpha=0.3$** — well
below the $\theta'=0.56$ the reduction consumes. What is held there is
not everything: at $\alpha=0.3$ the predicted lean is $0.87$–$0.91$ of
$\mu$'s, flat but short, and it climbs to $0.92$–$0.99$ at $0.4$ and
$0.99$ at $0.5$. The degradation and the level are two different
statements, and only the first has a threshold.

*Correcting Remark [rem:sievedepth].* It says of the fixed level that
the agreement and the lean ratio are "both getting worse as $N$
grows". The lean ratio is, at $3.52$ standard errors. The agreement's
slope is $-0.010101$ at $1.10$, and on six points that is not
resolved.


#### Remark (the missing piece is sieve depth, and that is the bad news) {#rem:sievedepth}
<!-- evidence: audit_sieve_depth.py -->

Remark [rem:logweightpredictor] ruled out the $\Lambda$ weight. What
had never been varied is the other half of $P$'s definition: the sieve
runs over the odd primes below $30$, a number that was inherited.
Sweeping it is not a search over models — a term survives $P_Q$ when
$N-mk$ has no odd prime factor at or below $Q$, and at
$Q=\lceil\sqrt N\rceil$ that set is the primes, so the sweep is a path
from the published predictor to $H$'s own support. (The published
weight's skip of primes dividing $k$ is vacuous: $k$ is coprime to
$N$, so such a prime never divides $N-mk$.)

**R1 holds** at $Q=29$ on both the agreements and the lean ratios.
**R3 and R4 hold, decisively.** At $Q=\lceil\sqrt N\rceil$ the
agreement is $0.9884$ to $0.9969$, the predicted lean is $0.9885$ to
$0.9942$ of $\mu$'s, and the decay slope is $-0.165946$ against
$\mu$'s $-0.167260$ — **$0.05$ standard errors.** The sign of the
dilated wall is carried, essentially exactly, by which $m$ make
$N-mk$ prime.

**And that is the bad news.** The agreement climbs
$0.7579\to0.7830\to0.7998\to0.8408\to0.8688\to0.9124\to0.9962$ across
$Q=29,\,53,\,101,\,211,\,503,\,1009,\,\lceil\sqrt N\rceil$ at the top
$N$: it reaches $\mu$ only as the level reaches $\sqrt N$. At the
fixed level $29$ it is $0.7367$–$0.8129$ and the lean ratio
$0.6571$–$0.8676$, **both getting worse as $N$ grows**. So the lean is
not a bounded-modulus object. Remark [rem:provablehalf] calls $P$
elementary because every condition defining its summation set is
multiplicative or a residue condition to a *bounded* modulus; the
predictor that carries the lean is exactly the one that gives that up.

**R2 is refuted, by arithmetic and not by data.** The agreement ties
at $N=8\cdot10^5$ and falls by $0.0001$ at $N=4\cdot10^5$, both past
$\lceil\sqrt N\rceil$, where the ladder's rungs overshoot: beyond the
square root no composite has a factor left to remove, and raising $Q$
further only strikes the small primes themselves, which are true
contributors. That is why the top rung is $\lceil\sqrt N\rceil$.

So OPEN's item — what carries the lean where the lean is measured —
closes, in the direction that costs the programme rather than helps
it. Nothing of bounded level carries it, and the thing that does is
primality itself.


#### Remark (the missing weight is not log(N-mk)) {#rem:logweightpredictor}
<!-- evidence: audit_logweight_predictor.py -->

Remarks [rem:oddmertensrange] and [rem:survivorrange] leave the sign
lean without an elementary carrier on its own $k$, and they fail in
opposite directions — $M_{\rm odd}$ overshoots with a flat decay, the
sieve-weighted $P$ undershoots and decays too fast. Opposite failures
say something is missing, and one thing is: $H(N;k)$ weights each
surviving $m$ by $\Lambda(N-mk)=\log p$ while $P$ counts it. Over a
short inner sum that hardly matters; over a long one $N-mk$ runs from
about $N$ down to about $k$, a factor of seven in the logarithm here.
So the natural candidate is
$P_{\log}=\sum_m\mu(m)w(m,k)\log(N-mk)$, elementary in the same sense
$P$ is.

**K1 holds** — $P$'s agreements and lean ratios reproduce exactly —
and **K2, K3 and K4 are all refuted.** The log weight *lowers* the
sign agreement at every $N$, $0.7955$ down to $0.7484$ against $P$'s
$0.8129$ down to $0.7579$. It improves the lean, $0.7018$–$0.9212$ of
$\mu$'s against $P$'s $0.6571$–$0.8676$, but not to within the
registered factor $1.2$. And its decay is $-0.231781$ against $\mu$'s
$-0.167260$, $2.48$ standard errors out — no better than $P$'s $2.11$.

**So the missing piece is not the $\Lambda$ weight**, and the two
scores move in opposite directions, which is worth recording as a
methodological point: a candidate that is best of three on the lean
ratio and worst on the agreement can be called an improvement by
choosing the score afterwards. The criterion registered here is the
agreement, the score these predictors were built for, and on it $P$
remains the best of the three.

The ledger, on $k<N^{0.56}$:

$$
\begin{array}{l|ccc}
 & \text{agreement} & \text{lean ratio} & \text{decay slope}\\\hline
M_{\rm odd} & 0.52\text{–}0.62 & 1.73\text{–}2.96 & -0.004298\\
P & 0.74\text{–}0.81 & 0.66\text{–}0.87 & -0.230633\\
P_{\log} & 0.72\text{–}0.80 & 0.70\text{–}0.92 & -0.231781\\
\mu & 1 & 1 & -0.167260
\end{array}
$$


#### Remark (which half of the predictor does the cutting) {#rem:predictablenull}
<!-- evidence: audit_predictable_null.py -->

Remark [rem:predictable] declines a control: "a randomisation would
move both parts together". That is true of one that moves both. Two do
not — keep $\operatorname{sign}P$ and permute $|P|$ across $k$, or
keep $|P|$ and redraw the sign — and each breaks exactly half the
predictor while leaving $H$ untouched. Given Remark
[rem:survivorrange], which half is doing the work now matters: the
sign field is the part validated only on short inner sums.

**D1 holds** — $\beta$ and all five residual shares reproduce exactly
— and **D2 holds**: the fall is $-0.037002$ per unit $\log N$ at
$6.61$ standard errors, so "the share it cuts is growing with $N$"
survives the test that withdrew rule U4.

**D4 holds and D3 is refuted, narrowly and informatively.** Scored on
the residual share, the sign-only null leaves $0.9252$ to $0.9543$,
just under the $0.95$ registered. Scored on the *cut*, which is what
the remark claims, the picture is unambiguous:

$$
\begin{array}{r|ccccc}
N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6 &
  3.2\cdot10^6\\\hline
\text{true cut} & 0.3690 & 0.4021 & 0.4446 & 0.4579 & 0.4693\\
\text{signs only} & 0.1012 & 0.1465 & 0.1512 & 0.1158 & 0.0984\\
\text{magnitudes only} & -0.0325 & -0.0173 & -0.0065 & -0.0032 &
  -0.0014
\end{array}
$$

The sign field alone reaches a tenth to a seventh of the cut; the
magnitudes alone reach none of it and are slightly worse than nothing,
as a random relative sign must be. **The $37$ to $47$ per cent needs
the per-$k$ pairing of the two.**

So the declined null was not empty, and what it settles is favourable
to the split as a statement about $P$ and unfavourable to reading it
through the sign agreement alone. The split does not inherit Remark
[rem:survivors]'s window restriction, because it is not the sign
agreement that carries it.


#### Remark (the sieve-weighted predictor on the lean's own k) {#rem:survivorrange}
<!-- evidence: audit_survivor_range.py -->

Remark [rem:oddmertensrange] withdrew $M_{\rm odd}$'s claim on the
lean and left one predictor standing. Remark [rem:survivors]'s
sieve-weighted $P(N;k)=\sum_m\mu(m)w(m,k)$ agrees with
$\operatorname{sign}H$ at $0.9274$ down to $0.9080$ — but on the same
window, $2\le N/k\le10^3$, which Remark [rem:predictable] states
plainly. The lean runs to $N/k=2.1\cdot10^6$.

**S1 holds.** Over every admissible $k$ in the published window rather
than its $60000$-$k$ subsample, the agreement comes back at $0.9106$
against $0.9080$ at the top $N$ and within $0.005$ throughout, and the
$M_{\rm odd}$ column likewise.

**S2, S3 and S4 are refuted, and narrowly.** On $k<N^{0.56}$ the
agreement is $0.8129,\,0.7632,\,0.7446,\,0.7367,\,0.7759,\,0.7579$ —
well below the $0.91$ of the demonstration. The lean $P$ predicts on
$\mu$'s own magnitudes is $0.6571$ to $0.8676$ of $\mu$'s, and its
slope is $-0.230633$ against $\mu$'s $-0.167260$, $2.11$ standard
errors apart.

Set beside the other predictor on the same window:

$$
\begin{array}{l|ccc}
 & \text{ratio to }\mu\text{'s lean} & \text{slope} &
   \text{agreement}\\\hline
 M_{\rm odd} & 1.73\text{–}2.96 & -0.004298 & 0.52\text{–}0.62\\
 P & 0.66\text{–}0.87 & -0.230633 & 0.74\text{–}0.81\\
 \mu & 1 & -0.167260 &
\end{array}
$$

**So no elementary predictor carries the lean where the lean is
measured**, and the two fail in opposite directions: $M_{\rm odd}$
overshoots by $73$ to $196$ per cent with no decay at all, $P$
undershoots by $13$ to $34$ per cent and decays too fast. $P$ is much
the closer of the two and its misses are marginal — a factor
$0.6571$ against a tolerance of $0.667$, a slope $2.11$ standard
errors out — but they are misses, and they are one-directional.

What this leaves of Remark [rem:survivors] is its own sentence with
the window attached: nine tenths of the sign of the dilated wall is a
sieve-weighted Möbius sum **over a short inner range**. On the range
the lean is measured over it is three quarters, and the quarter it
misses is where the lean's own trend lives.


#### Remark (the mechanism was shown on the wrong window) {#rem:oddmertensrange}
<!-- evidence: audit_oddmertens_range.py -->

Remark [rem:leanodd] establishes the agreement over $2\le N/k\le10^3$
with $k<N/2$ — **short inner sums** — and then credits it with
explaining Remark [rem:leandecay]'s decay. That lean is measured on
$k<N^{0.56}$, whose inner lengths run from $215$ to $2.1\cdot10^6$.
The two windows do not nest: the statistic overruns the demonstration
by a factor $2133$.

**P1 holds.** The agreement reproduces to four decimals at all five
$N$, $0.7704$ down to $0.7657$, and the full Mertens at $0.5472$ to
$0.5563$. Nothing in the mechanism as measured is in doubt.

**P2 holds, P3 and P4 are refuted.** On the lean's own $k$-range the
agreement is $0.5815,\,0.6161,\,0.5718,\,0.5618,\,0.5209,\,0.5201$ —
weaker than the demonstration and falling. Replacing each
$\operatorname{sign}H$ by $\operatorname{sign}M_{\rm odd}$ on $\mu$'s
own magnitudes gives a predicted lean of $0.48$–$0.50$ at every $N$,
which is **$1.8$ to $3.0$ times** $\mu$'s, and its slope is
$-0.004298$ at $0.94$ standard errors — flat, against $\mu$'s
$-0.167257$, a gap of $8.51$ standard errors. Against the floor of
Remark [rem:leanfloor] the two also part: $\mu$ rises at $+0.126950$
and the predictor at $+0.289908$, $3.40$ standard errors apart.

So the sentence to withdraw is "it also explains the decay without
appeal". **Where the lean is measured, $M_{\rm odd}$ predicts a lean
nearly twice too deep and predicts no decay at all** — it tracks the
floor's own decay instead, which is what a predictor blind to the
magnitudes would do. What stands is Remark [rem:leanodd]'s first
half, which is what P1 confirms: on short inner sums the sign of $H$
is the sign of $M_{\rm odd}$, and the parity restriction is the whole
difference. That is a fact about individual terms and Remark
[rem:survivors] is where its limits were already being measured.


#### Remark (which $m$ survive, and why the lean is deeper) {#rem:survivors}
<!-- evidence: lab_survivor_selection.py -->

Remark [rem:leanodd] settled the sign of each term and left the
selection open: the uniform single-survivor model under-predicts the
depth. The candidate is elementary. A term survives when
$\Lambda(N-mk)\ne0$, and $N-mk$ *shrinks* as $m$ grows — at the top of
the range it is of order $k$, where primes are far denser than near
$N$. So survivors should lean towards large $m$, and over a short
range the large odd squarefree $m$ are mostly primes, with $\mu=-1$.
The same argument makes $m=1$, the one guaranteed $\mu=+1$, the least
likely survivor, since $N-k$ is the largest value in the range.

Both consequences are there. Among sampled $k$ with exactly one
surviving $m$, the fraction whose survivor is $m=1$ reads
$0.4710,\,0.4678,\,0.4670,\,0.4601,\,0.4550$ against the uniform
expectation on the *same* $k$, $0.5420,\,0.5381,\,0.5305,\,0.5209,\,
0.5159$ — a deficit of $12$ to $13$ percent. And the survivors' mean
relative position $m/(N/k)$ is $0.5298,\,0.5293,\,0.5290,\,0.5286,\,
0.5268$ against the $0.5$ a uniform draw would give.

Replacing the count by the density each term actually has settles it.
With the standard sieve weight for $N-mk$ to be prime,
$w(m,k)=\prod_{q\nmid k}\bigl[N\not\equiv mk\ (q)\bigr]\,q/(q-1)$ over
the odd $q\le30$, the predictor $P(N;k)=\sum_m\mu(m)w(m,k)$ agrees
with $\operatorname{sign}H(N;k)$ at
$0.9274,\,0.9224,\,0.9204,\,0.9154,\,0.9080$ — against $0.8377$ down
to $0.8238$ for the odd Mertens function of Remark [rem:leanodd], and
against $0.5372$ to $0.5414$ for the best of sixteen permutations.
**Nine tenths of the sign of the dilated wall is a sieve-weighted
Möbius sum over its own inner range.**

That closes the chain [rem:signmass] opened. The lean is not a
property of $\mu$ that cancellation across $k$ might undo: it is
$\sum_m\mu(m)$ over a short range, biased by the prime density that
selects which $m$ appear. What it does not do is help. The
$9$ percent that $P$ misses is where any saving would have to live,
and it is the part no elementary weight predicts.


#### Remark (how much of the demand is elementary) {#rem:predictable}
<!-- evidence: lab_predictable_part.py -->

Remark [rem:survivors] validated $P(N;k)=\sum_m\mu(m)w(m,k)$ on short
inner sums, $N/k\le10^3$, and called the $9\%$ it misses the place a
saving would live. The demand $B_H(N;K)=\sum_{k<K}(\log k)|H(N;k)|$
runs down to $k=2$, where the inner sum has $N/2$ terms, so the
question is whether the predictor survives there and whether that is
where the mass is. Measured over the operative range $2\le k<3\cdot10^4$,
which covers $K^*$ at every $N$:

**All four pre-registered rules fail, and every one of them fails in
the direction that makes the elementary part matter more.** **S1**
predicted the correlation of $H$ with $P$ would collapse at long inner
sums, below $0.2$; it is $0.6149,\,0.8637,\,0.8599,\,0.8621,\,0.8096$
there against $0.8583$ down to $0.8492$ at the short end. **S2**
predicted the mass would sit at the long end; the share of
$\sum(\log k)|H|$ coming from $N/k\ge10^4$ is
$0.0171$ to $0.0384$, while $N/k\le10^3$ carries
$0.7705$ to $0.8991$ — there are simply far more $k$ at the short end.
**S4** predicted the sign agreement would be void where the mass is;
at $N/k\ge10^4$ it is $0.8462,\,0.8148,\,0.7925,\,0.8148$ at the four
largest $N$.

**S3** is the consequential one. With $\beta$ the least-squares scale
through the origin — $\beta=2.6992,\,2.8588,\,2.9952,\,3.1437,\,
3.0473$, stable — the elementary part is not a correction but most of
the object: $\sum(\log k)|\beta P|$ is $0.8595$ to $0.9553$ of
$\sum(\log k)|H|$, and subtracting it leaves
$\sum(\log k)|H-\beta P|/\sum(\log k)|H| = 0.6310,\,0.5979,\,0.5554,\,
0.5421,\,0.5307$. Removing everything an elementary sieve weight
predicts cuts the demand by $37$ to $47$ percent, and the share it
cuts is *growing* with $N$.

What that does not do is reduce the difficulty. The residue is still
of the same order as $B_H$ itself, so no bound follows; and $P$ is a
Möbius sum over the same short range, not a quantity with a known
bound. What it establishes is the split: **to within a factor $\beta$
of about $3$, the dilated wall is an elementary sieve-weighted Möbius
sum plus a residue of half its size**, and the residue is the only
part that is a Möbius–prime correlation in any essential sense. That
is a sharper localisation than Remark [rem:survivors]'s $9\%$, which
was measured where the mass is thinnest.


#### Remark (a copied floor, and the control M3 was owed) {#rem:copiedfloor}
<!-- evidence: lab_cell_singular.py -->

The fourth and last pointer-style decline. Proposition
[prop:scaleinv]'s evidence declines a control because $D_c$ "is a
deterministic arithmetic functional" and "the detection it explains
has its own control" elsewhere. Both halves are true and neither
covers **M3**, the surviving claim of that script: a correlation of
$0.9805$ across depths between $D_c$ and the exact floor $se_c$. The
control it points at permutes depth labels to test $z_c$, a different
statistic.

Auditing it turned up something first. The floor $se_c$ was not
computed there at all — six numbers were typed in, with a comment
naming the script they came from. **A copied value is a dependency no
check can see**: the gate compares a script with its own result, and a
script with the files it reads, and a number that was typed is
neither. It is now read from that script's result file, which puts it
under the ordering check, and the program's gate refuses a value whose
comment names a source it does not read.

**M3 survives the control it was owed, and by the widest margin
available.** With six depths the permutation null is exact: over all
$720$ orderings of the depth labels, exactly $1$ reaches
$|r|=0.9805$ — the identity — against a null median of $0.3475$ and a
$95$th percentile of $0.8084$. So the floor tracks $D_c$ at
$p=1/720$, and the tracking is not an artefact of both being monotone
in depth.

The other three rules of that script were already refuted and remain
so; what this adds is that its one surviving claim now has the null
its evidence declined, and that the number it correlates against is
computed rather than copied.


#### Remark (the decay sweep, controlled) {#rem:decaynull}
<!-- evidence: audit_decayfamily_null.py -->

Remark [rem:decayfamily]'s evidence declines a control by pointing at
the scripts that measured $f$, on the ground that it "only fits a
functional form to numbers already controlled there". Fitting a
functional form is not nothing: the $\alpha$-sweep of
$\log|1/2-f| = a - c(\log N)^\alpha$ reports a minimiser, a
one-percent band and a nine-order extrapolation, and none of those is
controlled by having measured $f$. Eight noisy points will produce a
minimiser and a band whatever they contain — or so it seemed worth
checking.

**The pointer missed, and half of what follows is superseded.**
[rem:decaydeep] re-measured the size separation below at $512$ iid
draws and $512$ random multiplicative ones: no iid draw reaches
$\mu$ at either of the two smallest $N$, and $67$ and $54$
multiplicative draws do. The claim that every draw sits an order
below $\mu$ is true of the control run here and false of the
wider one, so what follows about the *size* is a statement about
an iid coin and not about $\mu$. The $\alpha$-sweep half was untouched then and is
superseded now: [rem:decaysweep] finds six multiplicative draws of
512 with both an interior minimiser and a residual at or below
0.011556, the best at 0.000003, while no iid draw manages it in
512. Both halves of this block rest on an iid coin. With the field,
weight, $k$-range and $\theta'$ identical and only the sign pattern
changed, $\mu$'s $|1/2-f|$ reads
$0.2727,\,0.2772,\,0.2265,\,0.1932,\,0.1793,\,0.1624$ against a coin
maximum of $0.0826$ — every draw an order below $\mu$ at every $N$,
so the lean the pointed-to control established is confirmed here too.
And the sweep discriminates decisively, against the prediction:
**not one** of the eight draws has an interior minimiser, every one
pinning at $0.05$ or $1.50$, with residual sums of squares from
$2.991540$ to $10.078929$ against $\mu$'s interior $\alpha^*=1.45$ at
$0.011556$ — a separation of two to three orders in fit quality.

Two rules fail. **Y3** predicted the sweep would be uninformative and
it is not, which is the better outcome and the reason the rule was
written to be falsifiable in that direction. **Y2** asked $\mu$'s
$|1/2-f|$ to fall monotonically; it rises once, from $0.2727$ to
$0.2772$ at the first step, and falls cleanly thereafter — while
$0$ of $8$ draws fall monotonically at all, which is the half of Y2
that carried the content. **Y4** holds: the median draw's band is
$1.02$ against $\mu$'s $0.60$, so band width alone remains a weak
discriminator and Remark [rem:decayfamily]'s own note — that its
one-percent band is a threshold and not a null — stands, now measured
rather than conceded.


#### Remark (a pointer that missed, and the control it owed) {#rem:forecastnull}
<!-- evidence: audit_forecast_null.py -->

Remark [rem:forecast]'s evidence declines a control by pointing
elsewhere: $\gamma$ carries what the sign pattern contributed, and
that was measured for $\mu$ and for a coin in the dilation
extrapolation. The pointer is real and it misses. What was measured
there is the *median* of $|A(N;k)|/\sqrt{N/k}$; the forecast's own
diagnostic then rejects the median — "$\gamma$ was calibrated on the
median, but $B$ is a sum and needs the mean" — and rebuilds on a
mean-based $\gamma$. **The coin was never measured for the mean**, and
the two differ by half again: $1.6429,\,1.6259,\,1.4888,\,1.4725,\,
1.5033$ across the range. This is the failure Remark [rem:splitnull]
found in the split — a named control testing a different statistic
from the claim.

By the criterion of Remarks [rem:weightgapnull] and [rem:extendnull]
the control is usable: a mean of absolute values cannot be driven near
zero, and measured, the eight draws' $\gamma$ spread is $0.0434$ of
their median. So it is run, and **the calibration is $\mu$'s.** The
mean of $|A|/\sqrt{N/k}$ is $3.5647$ to $3.7497$ for $\mu$ against a
draw median of $2.4882$ to $2.7481$, every draw below $\mu$ at every
$N$ as Remark [rem:whycoinwins] requires; $\mu$'s fitted
$\gamma=1.0029$ sits outside the draws' range $[0.6963,\,0.7271]$; and
feeding each $\gamma$ through the same model, the $N$ at which $K^*$
first reaches $\sqrt N$ is $10^{6.20}$ for $\mu$ against
$[10^{4.46},\,10^{4.83}]$ for the draws — outside by more than an
order of magnitude.

Two things are worth stating exactly. The crossing audited here is the
level-one-half one, $K^*=\sqrt N$, not the $\theta'=0.56$ crossing
Remark [rem:forecast] quotes; they share the calibration, so what
transfers is that $\gamma$ is $\mu$'s and therefore so is any crossing
the model produces from it. And an earlier version of this audit
returned every crossing at the top of its bracket, because the
bisection moved the wrong way below the crossing — the numbers above
are from the corrected solve, and the broken ones were an instrument
fault, not a measurement.


#### Remark (the extrapolation's declined null, run) {#rem:extendnull}
<!-- evidence: audit_extendrange_null.py -->

Remark [rem:extendrange]'s evidence declines a control — "no new
detection is claimed, only a longer lever" — and its content is an
extrapolation, which is exactly the kind of claim a control reaches.
Remark [rem:weightgapnull] settled when a declined control is worth
running: when the statistic is well conditioned under it. $B(N)/N$ is
a sum of absolute values, so a coin cannot drive it near zero; across
eight draws its spread is $0.1510$ of the median at the smallest $N$
and $0.0436$ at the largest. This is the well-conditioned case, so the
control is run.

**The extrapolation is $\mu$'s.** With the field, the weight $\log k$,
the $k$-range and $\theta'$ identical and only the sign pattern
changed, the coin's $B/N$ is below $\mu$'s at every $N$ and every
draw — median $0.5399$ against $0.8086$ at $2\cdot10^5$, $0.3984$
against $0.5526$ at $6.4\cdot10^6$ — which is what Remark
[rem:whycoinwins] requires. The fitted log-law exponent is $-1.5274$
for $\mu$ against a draw range of $[-1.3570,\,-0.8243]$, outside it;
and the crossing of the Goldbach threshold sits at $10^{8.86}$ for
$\mu$ against $[10^{7.17},\,10^{7.87}]$ for the draws, outside it by
more than an order of magnitude. So the bracket Remark
[rem:extendrange] quotes is a statement about $\mu$ and not about the
fitting exercise.

One rule fails, and its failure was flagged in advance as the good
outcome. **W2** predicted that a good fit would be uninformative, every
coin fitting past $0.98$ in magnitude as $\mu$ does. Three draws do
not, running down to $0.91992$, so fit quality carries *some*
information — but not much, since $\mu$'s $0.98333$ is beaten by one
draw's $0.99494$. **The quality of the fit is a weak discriminator and
the rate is a strong one**, which is the right way round for an
extrapolation and was worth establishing rather than assuming.


#### Remark (the control the split never had) {#rem:splitnull}
<!-- evidence: lab_split_null.py -->

Remarks [rem:predictable] and [rem:residue] were built without one.
lab_survivor_selection.py permuted signs and showed $P$ predicts the
*sign* of $H$; nothing tested the *size* claim, and the size claim is
what the split rests on. It could easily have been empty: $P$ and $H$
are sums of the same length with bounded weights, so a least-squares
fit will absorb some fixed share of the mass whatever $P$ contains.

Two controls settle it, and each breaks $P$ without changing its
shape. Replacing $\mu(m)$ by a fixed $\pm1$ on the odd squarefree $m$,
sieve weights and summation range untouched, the residual share
$\sum(\log k)|H-\beta P|/\sum(\log k)|H|$ rises from $\mu$'s
$0.6310,\,0.5979,\,0.5554,\,0.5421,\,0.5307$ to a median over eight
coins of $0.9879,\,0.9566,\,0.9711,\,0.9905,\,0.9986$ — **the coin
absorbs essentially nothing**, and the gap $0.3569$ to $0.4680$ is
twice what the rule required. Replacing the sieve weight by $1$ and
keeping $\mu$, the plain Möbius sum reaches correlation
$0.6570$ down to $0.5532$ against $P$'s $0.8575$ to $0.8463$, and
leaves $0.8710$ to $0.8272$ of the mass against $P$'s $0.53$ to
$0.63$. So both ingredients earn their place: the arithmetic does the
work, and the sieve weight roughly halves what is left.

One rule fails, on its band rather than its content. **U2** asked the
coin correlations to sit in $[-0.1,0.1]$; they run $-0.2918$ to
$+0.2188$. The naive sampling error at $10132$ points would be
$0.0099$, but the inner sums are *nested* — consecutive $k$ share
almost all their $m$ — so the effective degrees of freedom are the
number of octaves, not the number of $k$, and $\pm0.29$ is the honest
noise scale. Read against the band the coins themselves span, $\mu$
clears the far edge by $+0.6300$ to $+0.8291$. The band should have
been computed from the controls rather than assumed, which is what
[rem:cap]'s rule says and what U1 did correctly.


#### Remark (the residue is square-root) {#rem:residue}
<!-- evidence: lab_residue_size.py -->

Remark [rem:predictable] split $H=\beta P+R$ and left the question the
whole route turns on. Remark [rem:directlevel] showed that the measured
level $K^*_H$ is exactly what square-root cancellation in $H$ would
predict; the split makes that testable on the part that is actually a
Möbius–prime correlation.

**It is square-root, and the lean is not in it.** Fitting the octave
means of $|R|$ against $N/k$ — with each octave's abscissa the mean of
$N/k$ inside it, not a nominal midpoint; bins closed at both ends and
required to hold at least ten $k$, as Remark [rem:elemreach] forced —
gives exponents

$$
0.4722,\ 0.4923,\ 0.5079,\ 0.5048,\ 0.4869
$$

at all five $N$, with correlations $0.99582$ to $0.99961$ and
leave-one-out spreads $0.0106$ to $0.0419$. The floor is a threshold
and is swept: over $5$, $10$ and $20$ $k$ per octave the exponent
moves by at most $0.0374$ (rule T5, added with the correction and
disclosed as such). The thinnest octave any of these five fits stands
on holds $25,\,13,\,25,\,13,\,25$ values of $k$ — declared to the gate
now, along with the correlations, rather than left to be discovered,
as in Remark [rem:elemsize]. Meanwhile the mass-weighted $f_+$ of $P$ is
$0.0796,\,0.0649,\,0.0941,\,0.1326,\,0.1836$, an order below the
$0.4854$–$0.5155$ band that sixteen sign draws on $|R|$ span. **The
elementary sieve part carries the lean; the residue carries the
square root.**

**This remark first reported something weaker, and the correction is
the point.** As it stood, its fit ran through an unbounded top bin
$[32768,\infty)$; T2 and T4 then failed at $N=2\cdot10^5$, and the
remark diagnosed the cause correctly — at that $N$ the top bin means
$k\le6$ and holds a handful of them — but treated it as one bad $N$
rather than as a defect in the estimator. Remark [rem:elemreach]
later measured how large that defect is on a longer lever. With every
bin closed and populated, **T2 and T4 hold at all five $N$**, the
exponents sit on $\tfrac12$ instead of running from $0.47$ to $0.49$,
and the old numbers are not reprinted here because nothing computes
them any more.

**T1** still fails everywhere: the residue is not perfectly
centred, its $f_+$ reading
$0.5516,\,0.5573,\,0.5484,\,0.5379,\,0.4832$ against draw bands of
width about $0.02$. So $\beta P$ slightly over-corrects, and the
residual lean shrinks and changes sign across the range — an order of
magnitude smaller than $H$'s own $0.09$ to $0.38$, but not zero.

What this gives is a clean statement of the remaining difficulty.
$\sum_{k<K}(\log k)|H|$ is, to within the factor $\beta\approx3$, an
elementary sieve-weighted Möbius sum plus a residue of half its size
that obeys $|R|\asymp(N/k)^{1/2}$ over the whole operative range. The
heuristic of Remark [rem:directlevel] — $K^*\asymp\SS^2N/(4\log^2K)$
from square-root cancellation — is therefore a statement about $R$
alone, and what is unproved is square-root cancellation for $R$, not
for $H$. That is a smaller object than the one the program started
with, and no easier: it is still a Möbius–prime correlation of length
$N/k$ with the elementary part removed. Remark [rem:residuecancel]
asks whether that size is bought or given.


#### Remark (the head is not a range, but its sign is the lean) {#rem:headidentity}
<!-- evidence: audit_head_identity.py -->

Remark [rem:gainsplit] left the shortfall in the top tenth by mass and
OPEN item 5 turning on where that tenth's single-sign fraction goes.
Before fitting a shape to it, two structural questions. Since
$|a_k|=(\log k)|H(N;k)|$ and $|H|$ is largest where the inner sum is
longest, the head ought to be the small-$k$ end; and its sign ought to
be the lean of Remark [rem:signmass].

**Z2 is refuted.** The overlap between the top decile by $|a|$ and the
smallest decile of $k$ is $0.2174$ to $0.3263$ — the head is spread
across the range, with its median member at the $0.22$–$0.32$ point of
the $k$-order. The $\log k$ weight grows with $k$ while $|H|$ falls,
and the two pull against each other. Dropping the weight recovers the
expected picture: the head overlaps the top decile by $|H|$ alone at
$0.7794$ to $0.8863$. **So a mass-ranked split is not a range
restriction**, and looking for an arithmetic description of "the
head" as a set of small $k$ is looking for something that is not
there.

**Z1, Z3 and Z4 hold.** The head's majority sign is negative at every
$N$, and its single-sign fraction minus the whole range's
mass-weighted negative share is

$$
0.2273,\ 0.2010,\ 0.1853,\ 0.1868,\ 0.1924,\ 0.2097,\ 0.2078,\
0.1882 ,
$$

a spread of $0.0421$. The head's fraction is a proportion over $31$
to $475$ terms, so it carries a binomial error of its own; the
expected span of eight such draws is $0.0626$, **larger than the
spread observed**. The difference is constant to within its own
estimation error.

So the head's fraction is not a free quantity with a limit of its
own: it is the lean of Remark [rem:signmass], offset by a constant and
read over the heaviest tenth. OPEN item 5 is the lean question again,
and Remarks [rem:sievedepth] and [rem:levelweighted] have already
found what carries the lean — the primality of $N-mk$, at a level no
bounded modulus reaches.


#### Remark (half of the lean is in the head, and half is not) {#rem:signmasshead}
<!-- evidence: audit_signmass_head.py -->

Remark [rem:signmass] is where the lean was first located, and it is
titled *the correlation sits in the large terms*. Remark
[rem:gainsplit] later put the cross-$k$ shortfall in the top tenth by
mass and Remark [rem:headidentity] locked that tenth's sign to the
same lean. Those readings ought to be one object. They are not, and
the title is half true.

The control reproduces Remark [rem:signmass] exactly: over the
squarefree $k<N^{0.56}$ coprime to $N$, the share of $\sum(\log k)|H|$
carried by $H>0$ is $0.2273$ to $0.3207$, the share of $k$ with $H>0$
is $0.4121$ to $0.4808$, and split at the median of $|H|$ the small
walls sit at $0.5096$–$0.5774$ against the large ones' $0.3141$–$0.3841$.

**H2 is refuted: the two splits are different objects.** The high half's
positive-mass share minus the top decile's is $0.3141$, $0.2711$,
$0.2440$, $0.2641$, $0.2704$ — a span of $0.0701$, not a constant
offset. The median of $|H|$ and the top decile of $(\log k)|H|$ are
not two names for one cut, as Remark [rem:headidentity] should have
warned: mass ranking is not range ranking, and it is not magnitude
ranking either.

**H3 is refuted: the head does not carry the whole correlation.** The gap
between counting and weighting is $0.1848$, $0.1603$, $0.1561$,
$0.1325$, $0.1601$; with the top tenth by $|a|$ removed it is
$0.1038$, $0.0863$, $0.0881$, $0.0637$, $0.0884$. Taken at the same
$N$ the ratio is $0.4810$ to $0.5642$: **the heaviest tenth carries
$44$ to $52$ per cent of the correlation, and the other half is spread
over the remaining nine tenths.** What is left is not noise — under
$256$ re-signings that hold the magnitudes and destroy only the sign
pattern, the gap is $0.0040$, $0.0044$, $0.0004$, $0.0027$, $0.0004$,
so the residual $0.06$–$0.10$ is $\mu$'s doing at a scale twenty times
the null's.

So Remark [rem:signmass]'s title is exact only in the weak sense that
the large terms carry more than their share. **Removing them does not
remove the lean**, and an argument that controls the head alone
controls about half of what has to be controlled. This does not move
OPEN item 5's conclusion — the residual half is still $\mu$'s and
still carried by the primality of $N-mk$ (Remarks [rem:sievedepth],
[rem:levelweighted]) — but it removes the hope that the obstruction
lives on a thin set.


#### Remark (the tail cancels better than randomly; the head not at all) {#rem:gainsplit}
<!-- evidence: audit_gain_split.py -->

Remark [rem:flatnessshape] left OPEN item 5 as one question: has
$e(G)=0.153911$ any reason to reach $\theta'/2$? Remark
[rem:nocrossk] read the shortfall as the dilated walls moving
together, and its rule T4 measured the top decile's *mass* share at
$0.3486$–$0.3587$ — but never asked how much that decile cancels.
Splitting the $k$ by mass rank, in a fixed fraction so that each part
keeps $\#S\asymp N^{\theta'}$ and the same reference $\theta'/2$:

$$
\begin{array}{l|ccc}
\text{part} & \text{exponent} & \text{s.e.} & t\\\hline
\text{whole} & +0.153911 & 0.011253 & 13.68\\
\text{top tenth} & +0.077963 & 0.009252 & 8.43\\
\text{bottom nine tenths} & +0.340006 & 0.025785 & 13.19
\end{array}
$$

**Y1, Y2 and Y3 hold** — the control reproduces $G$ and its exponent
exactly, the head is below the whole and the tail above it. **Y4 is
refuted, and in the strong direction**: the tail is not below
$\theta'/2$ but $2.33$ standard errors **above** it. The many small
terms cancel *better* than random signs would; the whole shortfall
lives in the largest tenth, whose gain runs $1.0000$ to $1.4639$
where $\sqrt{\#\text{head}}$ is $5.6$ to $21.8$.

And they move together in the plain sense. The fraction of the head
carrying a single sign is

$$
1.0000,\ 0.9783,\ 0.9118,\ 0.8800,\ 0.8716,\ 0.8721,\ 0.8545,\
0.8274 ,
$$

above four fifths at every $N$ — at $N=2\cdot10^5$ every term in the
top decile has the same sign, which is why its gain is exactly one.

So Remark [rem:nocrossk]'s reading survives and sharpens: it is not
that every dilation fails to cancel, but that **a positive proportion
of them — a tenth, growing like $N^{\theta'}$ — carries one sign**,
and that alone holds $e(G)$ down. Reading the whole range's $0.154$
without the split hides a part that already exceeds square root and a
part that does not cancel at all.

*Added later.* Every number above is eight doublings. The same split
on the field to $1.024\cdot10^8$ is Remark [rem:splitreach]; the
readings survive and two of them change size.


#### Remark (the same split on eighty-one points) {#rem:splitreach}
<!-- evidence: audit_split_reach.py -->

Remark [rem:gainsplit] is the only place this programme says *where*
in the $k$ the deficit of item 4(b) sits, and all of it stands on
eight doublings spanning a factor of $128$. Its sharpest sentence is
the one about signs: the head's one-sign fraction runs $1.0000$ down
to $0.8274$ across those eight, falling monotonically, and whether it
falls to a positive limit or to a coin is the difference between an
obstruction and a small-$N$ effect. Remark [rem:fieldreach] made the
field affordable to $1.024\cdot10^8$, so the split runs on $81$
on-field $N$ over a spread of $6.2383$ in $\log N$. The eight
doublings reproduce to $0.000050$ (Z1).

**The ordering survives and the tail's margin triples.** On the $81$,
$e(\text{head})=+0.072900$, $e(\text{whole})=+0.149567$,
$e(\text{tail})=+0.457005$ (Z2), against $+0.077963$, $+0.153911$,
$+0.340006$ on eight. The tail's excess over $\theta'/2$ is
$+0.177005$ — **$9.25$ standard errors**, after $2.33$ (Z3). The small
terms do not merely beat square-root cancellation; the margin by which
they beat it grows with the reach. The tail is also the noisy part:
its r.m.s. is $0.310423$ against $0.040867$ for the head.

**The head is falling further behind, not catching up.** Its gain
against the square-root reference $\sqrt{\#\text{head}}$ runs
$0.1796$, $0.0927$, $0.0558$ at the bottom, middle and top of the
sweep, an exponent of $-0.207342$ at $t=82.40$. A head of $1034$
dilations at the top $N$ has gain $1.7957$ where independent signs
would give about $32$.

**The one sign is falling, and is nowhere near a coin.** The fraction
runs $1.0000$ to $0.7786$ with a slope of $-0.027628$ in $\log N$ at
$t=24.17$, so the fall [rem:gainsplit] saw across eight points is
real and resolved. But at the top $N$ the head has $1034$ members, a
coin's standard deviation is $0.015549$, and the fraction is $0.7795$
— **$17.97$ standard deviations above one half** (Z4). Both halves of
the reading hold at once: the head's alignment is decaying, and at
every $N$ this programme can compute it is not decay to randomness.
No arrival at one half is published from that slope; the shape
carrying it is undetermined, exactly as [rem:fieldreach] found for
$F$.

**One prediction failed: the head is not a fixed share of the mass**
(Z5). Its share of $\sum|a|$ runs $0.3628$ down to $0.3337$ with
exponent $-0.006738$ at $t=7.75$ — small, and resolved. This does not
contradict [rem:gainprofile], whose flat share was the top block of a
ten-way split of $|w_ds_d|$ and a different object; what it costs is
the assumption that exponents taken on "the head" are taken on the
same fraction of the mass at every $N$. They are not, by about a
thirtieth over the sweep.

**What this does to item 4(b).** The demand is that $e(G)=+0.149567$
close $+0.134019$ onto $e(\ell^1/\ell^2)=+0.283586$. The split says
that closing has to happen in the head: the tail is already past
$\theta'/2$ by $9.25$ standard errors and pushing it further buys
nothing, while the head sits at $+0.072900$ and its distance to the
square-root reference is widening at $-0.207342$. **So item 4(b) is
not a demand about $\mu$'s cancellation in general; it is a demand
about a tenth of the dilations whose signs agree four times in five
and whose agreement decays too slowly to reach a coin anywhere this
can be computed.**


#### Remark (which axis holds the head, and for how long) {#rem:headaxis}
<!-- evidence: audit_headaxis_reach.py -->

If item 4(b) is a demand about the head, what the head IS matters, and
the two remarks that answer stand on eight doublings.
[rem:headidentity] says it is not the small-$k$ end.
[rem:headsign] says the alignment does not live on the $k$ axis
either: cutting the same deciles on each factor of
$|a_k|=(\log k)\,T_k|I_k|$ separately, the top-minus-bottom spread of
the negative share is $+0.3389$ on $|I|$, $+0.1537$ on $T$ and
$+0.0126$ on $k$. On the field to $1.024\cdot10^8$ that control
reproduces to $0.000047$ (Y1), and two of the three readings survive
while the third — the one nobody had asked — does not.

**Alignment: the ordering holds, and $k$ carries none of it.** At the
top $N$ the spreads are $+0.3182$ on $|I|$, $+0.1654$ on $T$ and
$+0.0077$ on $k$ (Y2). Each decile there holds about $1033$ of the
$k$, so a top-minus-bottom difference of negative shares has a
binomial deviation of $0.022001$; $k$'s spread is $0.35$ of those and
$|I|$'s is $14.46$ (Y5). Four octaves past the doublings, the
arithmetic of the dilation still does not select the sign.

**Selection: the same axis, and the $k$ factor works against it.** The
variance of $\log|a|=\log\log k+\log T+\log|I|$ splits by covariance
into shares summing to one, and $\log|I|$ carries the largest at every
$N$ — $0.7703$ to $0.8903$, against $0.1302$ to $0.2673$ for $\log T$
(Y4). The $\log\log k$ share is **negative**, $-0.0383$ to $-0.0280$:
a larger $k$ shortens the inner sum, and the factor that multiplies
$a_k$ up is anticorrelated with the object it multiplies. So the head
is not a compromise between axes. It is the imbalance's, twice over —
what makes a term large is the same thing that makes it negative.

**Y3 is refuted, and that is the measurement.** The imbalance axis is
losing its grip. The $|I|$ spread runs $+0.5854$ down to $+0.2752$
with a slope of $-0.032884$ at $t=11.62$, while the $T$ spread rises
at $+0.033468$ at $t=12.99$ — opposite signs, near-equal rates. The
selection share moves the same way, $\log|I|$'s falling at
$-0.007297$ ($t=6.94$) with $\log T$'s taking it up. **The head's
description is correct now and is being handed over**: the alignment
that [rem:headsign] attributed to imbalance is migrating to mass as
$N$ grows.

No crossing point is published. Both slopes are least-squares fits
over one window and the shapes carrying them are undetermined, exactly
as [rem:fieldreach] found for $F$ and [rem:splitreach] for the
one-sign fraction; this repository has now three quantities whose
direction is resolved and whose destination is not. What the
refutation costs is a sentence, not a number: **"the head's signs live
on the imbalance axis" may not be written without the reach it was
measured at**, and the one statement that survives every reach here is
the negative one — $k$ carries none of it.


#### Remark (no prefix of the inner sum owns the head's sign) {#rem:headbounded}
<!-- evidence: audit_head_bounded.py -->

$H(N;k)=\sum_m\Lambda(N-mk)\mu(m)$ is a sum over an ordered index, and
nothing in this repository had asked which $m$ make it negative. The
question has a consequence: [rem:provablehalf] calls a condition
elementary when it is multiplicative or of bounded modulus, so if the
head's sign were carried by a bounded number of leading terms — a
signed sum of at most $M_0$ values of $\Lambda$ — the obstruction
would be an object of the kind the elementary half already handles.
Splitting $H=A+B$ at $m\le M_0$ for the fixed cutoffs $29$ and $1000$
answers it. The split is exact to $5.482\cdot10^{-15}$ relative and
the head's one-sign fraction reproduces [rem:splitreach] within the
bound its printing forces (R1).

**Twenty-nine leading terms know nothing** (R2, refuted). Their sign
agrees with $H$'s on the head at $0.5995$ to $0.9167$, but against an
arm drawn at their own marginal rate the excess runs $-0.78$ to
$+2.23$ — unresolved, and at the top $N$ it is $-0.13$ on an agreement
of $0.6015$. This is what G73 was added to catch: on a head that is
four-fifths one sign, an agreement of $0.60$ is below chance.

**A thousand leading terms do know, and increasingly** (R2, the other
cutoff). The excess runs $0.00$, $4.98$, $6.70$, $11.04$, $16.85$ at
$N=2\cdot10^5$, $1.6\cdot10^6$, $6.4\cdot10^6$, $2.56\cdot10^7$,
$1.024\cdot10^8$ — the agreement at the top is $0.8859$ where the
matched arm gives about two thirds. The registered R2 asked for
resolution at *every* $N$ and both cutoffs, so it is refuted; what
refutes it is $M_0=29$ having no skill and $M_0=1000$ having none at
the small $N$ where it is not predicting but containing (its share of
$\lvert H\rvert$ there is $0.8091$ and its majority share is
$1.0000$, which the DEGENERATE marker declares).

**But the prefix is not what makes the sign** (R3, refuted). Delete it
and nothing moves: at the top $N$ the remainder's one-sign fraction is
$0.7795$ for both cutoffs — the head's own value to four decimals,
$17.97$ binomial deviations above one half. And the prefix's share of
the magnitude is going away: exponent $-0.238967$ at $t=53.62$ for
$M_0=29$ and $-0.221725$ at $t=63.36$ for $M_0=1000$, the shares
falling $0.0698\to0.0136$ and $0.8091\to0.1897$ across the field (R4).

**What the two refutations say together is one thing, and it is not
what either asked.** A prefix carrying a fifth of the magnitude agrees
in sign with the whole $88.59$ per cent of the time where chance gives
about $68$, and deleting it changes the whole's sign essentially
never. **The alignment is coherent across the $m$-decomposition**: the
first thousand terms and the remaining ones are aligned with each
other, not merely both present. That is [rem:nocrossk]'s "the dilated
walls move together" one level down — across $m$ inside a single
dilation, where it had never been looked for.

**And it closes the route it was opened for.** No truncation isolates
the head's sign: the bounded part neither carries it nor is needed for
it, and its weight vanishes like $N^{-0.22}$. [rem:sievedepth]'s
verdict for the slope — nothing of bounded modulus carries it — now
holds for the head's sign as well, by a different argument on a
different object.

*Added later.* The coherence paragraph above is corrected, not
withdrawn: most of that $0.8859$ is the selection, and what survives
it is smaller and appears only above a threshold in $N$. See Remark
[rem:coherencenull]. The route-closing paragraph is untouched.


#### Remark (most of the coherence was the conditioning) {#rem:coherencenull}
<!-- evidence: audit_coherence_null.py -->

[rem:headbounded] ended on a claim neither of its predictions had
asked for, and the claim had a hole. The head is the top tenth of $k$
by $\lvert a_k\rvert=(\log k)\lvert A+B\rvert$, and selecting on the
magnitude of a *sum* preferentially selects pairs whose parts agree in
sign, because agreeing parts add and disagreeing parts cancel. M4 is
in this repository for exactly that, and the control had not been run.
It is run here: $\lvert A\rvert$ and $\lvert B\rvert$ are kept exactly,
only the two signs are redrawn, independently, at the marginal rates
observed at that $N$, and the selection is redone on the surrogate.
The control reproduces [rem:headbounded]'s head agreement to
$0.000016$ (S1).

**The worry was justified and large** (S2). The coin arm's head
agreement exceeds its own all-$k$ agreement at every $N$, by $+0.0724$
to $+0.3335$; at the top $N$ it gives $0.7440$ on the head against
$0.4985$ on all $k$. **So most of the $0.8859$ that remark reported is
the conditioning.** A pair of independent coins with those magnitudes,
filtered the same way, looks three-quarters coherent.

**S3 and S4 are refuted, and where they fail is the finding.** Both
registered "at every $N$", and both fail on the small ones: the head
excess runs $-2.42$ to $13.09$ and the all-$k$ excess $-14.15$ to
$+19.62$. At the bottom of the field the effect is not weak but
*reversed* — at $N=2\cdot10^5$ the observed all-$k$ agreement is
$0.1182$ against the arm's $0.4106$, because $\operatorname{sign}B$ is
negative on only $0.1086$ of the $k$ there and the remainder is nearly
a constant. That is a degenerate split, not a measurement of $\mu$.

**Above a threshold it is real, and it grows.** The all-$k$ excess is
under two deviations at $37$ of the $81$ $N$, the largest being
$5242880$, and is resolved above two at every $N$ from $5000000$
upward — $44$ of them. As an effect size rather than a ratio,
observed minus arm on all $k$ runs $-0.2924$ at the bottom to
$+0.0916$ at the top, and the excess in deviations rises with slope
$+5.904067$ at $t=28.60$. At the top $N$: all $k$, $0.5901$ against
$0.4985$ ($17.29$ deviations); on the head, $0.8859$ against $0.7440$
($12.14$).

**What the corrected statement is.** The two halves of the inner sum
are not independent in sign given their magnitudes — the all-$k$
version cannot be blamed on any selection, and it clears its arm over
the upper half of the field by nine points of agreement. But the
coherence is a fraction of what [rem:headbounded] read off the head,
it is absent below $N\approx5\cdot10^6$, and its direction there is
opposite. **The honest version of "the alignment is coherent across
the $m$-decomposition" is that it becomes so**, somewhere inside the
range this repository can compute, having been anti-coherent before —
which is one more quantity whose direction is resolved and whose
destination this repository has no way to reach.


#### Remark (the gain's denominator is a truncated Chebyshev correlation) {#rem:denominator}
<!-- evidence: audit_denominator.py -->

Everything this repository has done to item 4(b) treats
$\lvert\sum a\rvert$ as a cancelling sum whose smallness is to be
improved. Nobody asked what that sum is. Because $\mu*\log=\Lambda$,
the unrestricted double sum collapses:
$$\sum_k\log k\sum_m\mu(m)\Lambda(N-mk)
 =\sum_j\Lambda(N-j)\!\!\sum_{mk=j}\!\!\mu(m)\log k
 =\sum_j\Lambda(j)\Lambda(N-j),$$
the Chebyshev–Goldbach correlation. Computed both ways at
$N=20000$, $50000$, $100000$ the two agree to $4.945\cdot10^{-15}$
relative (T1). **So $\sum a$ is that correlation with four
restrictions on the index — $k<N^{\theta'}$, $k$ squarefree, $k$
coprime to $N$, $m$ coprime to $k$ — and nothing else.**

**One restriction does the damage, and it is the truncation** (T4).
Imposed alone it turns the sum negative: $-8657.0$, $-17462.3$,
$-23198.0$ against the unrestricted $+36539.9$, $+86993.4$,
$+176573.6$. Each of the other three alone leaves it positive —
squarefree $+27661.8$, $+61575.0$, $+131917.7$; coprime to $N$
$+36539.9$, $+86993.4$, $+176566.2$; $m$ coprime to $k$ $+51233.9$,
$+125171.4$, $+248504.2$. On the field the restricted sum is negative
at $81$ of $81$ $N$ (T2). The truncation does not shrink the main
term; it removes it and overshoots.

**What survives is below main-term order, and by how much is the
whole question** (T3). On the field, $\lvert\sum a\rvert$ has exponent
$\alpha=+0.717916$ with a standard error of $0.002376$, while
$\sum_j\Lambda(j)\Lambda(N-j)$ divided by $N$ is flat to $0.72$
standard errors — it is a main term, as it must be. So
$\lvert\sum a\rvert$ falls short of order $N$ by $-0.282084$, which is
$118.71$ standard errors. Meanwhile $e(\ell^1)=+0.867483$ and
$e(G)=+0.149567$, and the identity $e(G)=e(\ell^1)-\alpha$ closes to
$1.1\cdot10^{-16}$.

**In these terms the demand is arithmetic, not fitted.** $e(G)\to
\theta'/2$ is exactly $\alpha\to e(\ell^1)-\theta'/2=+0.587483$: the
truncated correlation must be smaller than it is by a factor growing
like $N^{+0.130433}$.

**And the shortfall is already a recognisable number** (X1, written
after T3 and not pre-registered). $1-\alpha=+0.282084$ against
$\theta'/2=0.2800$ — a gap of $+0.002084$, which is $0.88$ standard
errors of $\alpha$; against half the measured $\#k$ exponent,
$+0.280178$, the gap is $+0.001906$. **The truncation leaves the
correlation divided by the square root of the index set it kept:**
$\lvert\sum a\rvert\sim(\Lambda*\Lambda)/\sqrt{\#k}$, which is
square-root cancellation over the $k$ that were thrown away, measured
and not assumed.

Granting that, item 4(b) is one line. $e(G)=e(\ell^1)-\alpha$ and
$\alpha=1-e(\#k)/2$, so $e(G)\to e(\#k)/2$ **is** $e(\ell^1)\to1$:
the demand is that the $\ell^1$ norm of the weighted dilation sums
reach main-term order. It is measured at $+0.867483$, short by
$+0.132517$. That is the same distance item 4(b) has always had, now
attached to a quantity with a size rather than to a ratio — and it
says the obstruction is not that $\mu$ cancels too little, but that
$\sum_k(\log k)\lvert H(N;k)\rvert$ is too small to be a main term
while $\sum_k(\log k)H(N;k)$ is already as small as chance allows.

The equality $\alpha=1-\theta'/2$ is measured over this range at
$0.88$ standard errors and is not proved here; every consequence above
inherits that.

*Added later.* It is a coincidence at one level, not a law: swept over
$\theta'$, $\alpha$ *rises* where the square-root reading needs it to
fall. See Remark [rem:thetalaw]. What survives is T1 to T4 and the
restatement **at $\theta'=0.56$**; the words "square-root cancellation
over the $k$ that were thrown away" are withdrawn.


#### Remark (raising the level moves away from the demand) {#rem:thetalaw}
<!-- evidence: audit_theta_law.py -->

[rem:denominator] measured $1-\alpha=+0.282084$ against $\theta'/2$ at
$0.88$ standard errors and read it as a law — the truncation leaving
the correlation divided by $\sqrt{\#k}$ — which with
[rem:headsign]'s square-root inner cancellation predicts
$e(\ell^1)=(1+\theta')/2$, $\alpha=1-\theta'/2$ and
$e(G)=\theta'-1/2$, three lines in $\theta'$ testable by sweeping it.
The levels $0.40$ to $0.64$ were swept on the whole field; at
$\theta'=0.56$ the three exponents reproduce that remark to
$0.00000035$ against a printing bound of $0.00000050$ (U1).

**The law is wrong, and $\alpha$ is wrong in the opposite direction**
(U3, U6 refuted). $\alpha$ runs $+0.581914$, $+0.609324$, $+0.670962$,
$+0.717916$, $+0.764301$, $+0.808475$ — a slope of $+0.978885$ in
$\theta'$ where the model wants $-1/2$, off by $21.91$ standard
errors. In hindsight it must rise: as $\theta'\to1$ the truncation
becomes no truncation and $\lvert\sum a\rvert$ returns to the
main-term order of $\sum_j\Lambda(j)\Lambda(N-j)$, and $\alpha\to1$ is
what the numbers are doing. **So $1-\alpha=\theta'/2$ holds at one
swept level out of six and is a coincidence there**; at $\theta'=0.40$
it is $0.418$ against $0.200$ and at $0.64$ it is $0.192$ against
$0.320$. [rem:denominator]'s X1 sentence is withdrawn as a law and
kept as an arithmetic fact at $\theta'=0.56$.

$e(\ell^1)$ does rise (U2), at $+0.625088$ against the model's $1/2$,
which is also refuted at $8.34$ standard errors. Only its sign
survives the sweep.

**And the direction that matters is the pessimistic one** (U4, U5
refuted). $e(G)$ *falls* with the level — $+0.190221$, $+0.193499$,
$+0.169761$, $+0.149567$, $+0.130229$, $+0.111700$ — while the demand
$\theta'/2$ rises. The gap $e(G)-\theta'/2$ is negative at every swept
level and widens monotonically: $-0.009779$, $-0.036501$,
$-0.090239$, $-0.130433$, $-0.169771$, $-0.208300$, a slope of
$-0.853798$ at $t=16.25$. **Raising the level — the one axis the
reduction lets anyone choose — moves away from the demand, and
quickly.**

The two lines do cross, at $\theta'=0.4041$ with bracket
$[0.3017,\,0.5351]$, but on the wrong side: below the $\theta'>1/2$
the reduction requires, and the point estimate is below $1/2$
outright. The extrapolation should not be leaned on in any case — the
slope refitted on the lower three levels is $-0.170500$ and on the
upper three $-0.473335$, a drift of $0.302835$ that exceeds the
fitted slope's own magnitude, so $e(G)$ against $\theta'$ is curved
and the crossing is a linear reading of a curve.

**A registered rule was misspecified and it has to be said.** U4 was
written as "the gap is increasing in $\theta'$" with the gloss that
its refutation — a flat or shrinking gap — would be "the outcome
worth having". That gloss is backwards: the gap is negative, so a
*shrinking* gap is a widening shortfall, which is what was found. The
prediction is refuted by its letter and its substance is confirmed in
the pessimistic direction. The rule should have been written on
$\lvert e(G)-\theta'/2\rvert$.


#### Remark (the surviving correlation is a composite-index object) {#rem:support}
<!-- evidence: audit_support.py -->

[rem:denominator] showed $\sum a=\sum_j\Lambda(N-j)\Lambda_K(j)$ and
[rem:thetalaw] showed the truncation is what removes the main term.
Neither says what the survivor sits on, and here the index decides it
before any measurement. Write $j=mk$ with $k$ squarefree,
$2\le k<N^{\theta'}$, $k$ coprime to $N$, $m$ coprime to $k$:

* $j$ prime leaves only $(m,k)=(1,j)$, since the other factorisation
  has $\log1=0$. So the prime part is
  $\sum_{p<N^{\theta'}}\log p\,\Lambda(N-p)$ — **a Goldbach count with
  one prime forced below the level**.
* $j=p^e$ with $e\ge2$ leaves only $k=p$ among squarefree divisors
  above $1$, and then $m=p^{e-1}$ is not coprime to $k$. **Higher
  prime powers contribute nothing at all.**
* Everything else has at least two distinct prime factors.

The split reproduces $\sum a$ to $5.346\cdot10^{-15}$ relative (V1),
and the prime-power bucket is $0.0$ — exactly, not approximately (V2).
The signs separate completely: the prime part is positive at $81$ of
$81$ $N$ and the composite part negative at $81$ of $81$ (V3).

**The composite part is the whole object.** Its exponent is
$+0.716454$ against the total's $+0.717916$, while the prime part is
$+0.579228$ — below the total by $-0.138688$, which is $22.39$
standard errors of the difference. In absolute terms at the top $N$
the prime piece is $+57085.5$ against a total of $-8013222.0$, and at
the bottom $+1417.9$ against $-87895.3$: a Goldbach count is in there,
positive, and it is a hundredth of what it is inside.

**So item 4(b)'s denominator has almost no Goldbach content.** What
has to be made small is not the main term of the binary problem —
that is removed by the truncation, exactly, for every prime index
above the level — but a correlation of $\Lambda$ against a
convolution supported on indices with two or more distinct prime
factors. Whatever the obstruction is, it is not the Goldbach main
term standing in the way.

**One prediction failed, and modestly** (V4). The prime piece should
be of size $N^{\theta'}$: the Hardy–Littlewood count of $p<K$ with
$N-p$ prime is about $\mathfrak S K/(\log K\log N)$ and the summand
$\log p\log(N-p)$ cancels both logs, leaving $\mathfrak S N^{\theta'}$
with $\mathfrak S$ constant along this family. Measured, the exponent
is $+0.579228\pm0.005719$, above $\theta'=0.56$ by $+0.019228$ —
$3.36$ standard errors. It is the noisy part, with r.m.s. $0.092817$
against the composite's $0.037934$, but the departure is resolved and
is recorded rather than explained: at this reach the explicit
Goldbach piece grows slightly faster than its own heuristic says.


#### Remark (it is a correlation and not a mean) {#rem:meanonly}
<!-- evidence: audit_meanonly.py -->

Reordered, $\sum a=\sum_q(\log q)\Lambda_K(N-q)$: the truncation
defect sampled at the shifted primes. If $\Lambda_K$ had a nonzero
mean and the shifted primes saw only that, the whole denominator would
be $(\psi(N)/N)\sum_j\Lambda_K(j)$ — a divisor-sum average with no
primes in it, and item 4(b) would become an elementary question. The
control M4 asks for is to recompute with $\Lambda(N-j)$ replaced by
its own mean over the same range, keeping the index set, the weights
$\mu(m)\log k$ and the truncation exactly.

The mean does point the same way — the two agree in sign at $81$ of
$81$ $N$ (W2) — and it is not most of the sum (W3, refuted). The
residual runs $0.8555$ to $0.8919$ of $\lvert\sum a\rvert$; at the top
$N$ the sum is $-8013222.0$, the mean-only value $-1065195.2$ and the
residual $-6948026.9$. Nor is the residual of lower order (W4,
refuted): its exponent is $+0.714942$ against the sum's $+0.717916$
and the mean-only arm's $+0.739002$, a difference of $-0.002973$ at
$0.85$ standard errors.

**So the denominator is a genuine correlation with the primes**, about
seven eighths of it, and it grows at the same rate as the whole. The
simplification is not available: item 4(b) keeps its arithmetic
content, which is the outcome the registered rule named as leaving the
problem where it was.

The control W1 passes at $0.049264$ against a bound of $0.05$ — the
totals it checks against are printed to one decimal on numbers of size
$10^7$, so that control is coarse by construction and says only that
the same sum was computed.


#### Remark (the signed sum is the head's, and increasingly) {#rem:sumhead}
<!-- evidence: audit_sum_head.py -->

Two lines had not been joined. [rem:splitreach] localised the gain
deficit in the head — the top tenth of $k$ by $\lvert a_k\rvert$ —
and [rem:denominator] showed the demand is about $\lvert\sum
a\rvert$, the signed total. Every head measurement was about gains and
mass; nobody had measured how much of the *signed sum* the head
carries, which is the quantity item 4(b) is about. A ratio of two
published gains suggests an answer but cannot settle it, because the
two parts' sums may oppose and partly cancel.

They do not (X2): head and tail each agree in sign with the total at
$81$ of $81$ $N$. **And the head is nearly all of it** (X3). Its share
of $\lvert\sum a\rvert$ runs $0.5467$ to $0.9710$, rising with slope
$+0.053778$ at $t=25.60$; at the top $N$ it carries $-7082988.7$ of
$-8013222.0$ while holding $0.3421$ of the mass and a tenth of the
$k$.

**The tail is already below the target** (X4). Its exponent is
$+0.414033$ against the total's $+0.717916$ — lower by $15.42$
standard errors — and the demand's target, $\ell^2$ order, is
$+0.583897$. So the bottom nine tenths of the $k$ already cancel to
better than item 4(b) asks. The head's exponent is $+0.787845$, above
the total's, which is what a share rising from a half to nearly one
means.

**So item 4(b) is one statement about one tenth of the $k$.** The
whole needs $+0.134019$ of improvement to reach $\ell^2$ order; on
this split every bit of that belongs to the head, whose own signed sum
must come down from $+0.787845$ to below $+0.583897$ — further than
the whole has to travel, because the part that is failing is
travelling the other way. The tail may be dropped from the problem
entirely.

Together with [rem:support] this says what the object is with no
slack left in the description: **a correlation of $\Lambda$ against a
convolution supported on indices with two or more distinct prime
factors, restricted to the tenth of the dilations whose sums are
largest, whose signs agree four times in five and whose share of the
total is heading to one.**

*Added later.* "Further than the whole has to travel" is withdrawn:
the head's excess exponent is a transient of a rising share, and over
the top octave the two exponents already coincide. See Remark
[rem:headfraction]. The rest of this remark stands.


#### Remark (the excess was transient, and the set is not thin) {#rem:headfraction}
<!-- evidence: audit_head_fraction.py -->

[rem:sumhead] ended by saying the head's signed sum must come down
*further* than the whole. It must not. Head and tail agree in sign at
every $N$, so the head's share of $\lvert\sum a\rvert$ is bounded by
one, and an exponent above the total's can only mean a share still
climbing. Whether that is a transient is measurable two ways, and
this cycle got one of them wrong before the other settled it.

**The ceiling argument, as this script's own setup stated it, is
refuted** (Y3). The setup read [rem:sumhead]'s share range
$0.5467$–$0.9710$ and took the upper end for the value at the top
$N$. It is not: the share at the top $N$ is $0.8839$, the range
maximum falls elsewhere in the field, and at the fitted slope
$+0.053778$ the bound is $2.1586$ units of $\log N$ away — a factor
of nearly nine in $N$, not the sliver the setup claimed.

**The transience is real anyway, on the other measurement** (Y4, but
see [rem:alphalocal] — an octave fit of these exponents carries a
standard error near $0.034$, so the gap Y4 compares is not resolved
either way and the transience rests on the bound alone).
Refitted on the top octave alone, twelve points, the head's exponent
is $+0.667854$ and the total's $+0.664374$ — a gap of $0.003480$
against $0.069930$ over the whole field. The convergence the bound
forces has already begun inside the computed range, so the excess is
a within-range effect and **the head's signed sum has the same
$+0.134019$ to travel as the whole, not more.** Both octave values
sit below their whole-range fits; the script prints no error for a
twelve-point fit and nothing is read from that here.

**And the failing set is not thinner than a tenth** (Y2, refuted).
Sweeping the fraction, the top one per cent carries $0.0969$ to
$0.2342$ of $\lvert\sum a\rvert$ — below a half at $81$ of $81$ $N$.
The tenth of [rem:gainsplit]'s convention is about the right
description and not a loose one.

**The partial sums overshoot.** At fraction $0.20$ the share is
$1.1740$ at the top $N$ and at $0.50$ it is $1.1060$: the largest
fifth of the dilations sums to more than the whole does, so the
remaining four fifths oppose it. Ordering by magnitude produces
opposing tiers rather than a head that simply dominates, and the
exponents fall monotonically across the sweep — $+0.805134$,
$+0.806838$, $+0.796624$, $+0.787845$, $+0.764263$, $+0.737664$
against the whole's $+0.717916$ — which is the same convergence Y4
found, seen along the fraction axis instead of along $N$.


#### Remark (alpha is a constant of the range, and the octaves cannot see) {#rem:alphalocal}
<!-- evidence: audit_alpha_local.py -->

Six remarks have quoted $\alpha=+0.717916$ as a number the field has.
[rem:headfraction] refitted the top octave alone, got $+0.664374$, and
declined to read it for want of a printed error. Here the errors are
printed, and they settle three things at once.

**$\alpha$ is a constant of this range** (Z3, refuted). The nine
octave fits run $+0.626273$ to $+0.823194$ — a spread of $0.196921$ —
and their slope against mid $\log N$ is $-0.000807$ with a standard
error of $0.013688$, $t=0.06$. Flat, as flat as a fit gets. The
quantity that six cycles treated as constant is one. *(Put in question by
[rem:alphareach]: on sixteen octaves the same slope reads
$-0.012011\pm0.005910$, $t=-2.03$ — but dropping one octave moves it
to $-0.005791$, so neither reading is robust.)*

**And the top octave is not low** (Z2, refuted). Its $+0.664374$ comes
with a standard error of $0.034400$, so it sits $1.56$ of its own
standard errors below the global $\alpha$. [rem:headfraction]'s
refusal to read it was right; this remark supplies the number that
makes the refusal a measurement rather than a caution. It also costs
that remark's Y4 its evidential weight: a gap between two octave fits
each carrying $0.034$ cannot distinguish $0.003480$ from $0.069930$.
The head's excess exponent is still transient — the share is bounded
by one and must stop rising — but that follows from the bound, not
from the octave.

**The deficit is out of the octaves' reach** (Z4, refuted). Locally
$\alpha-e(\ell^2)$ runs $-0.035856$ to $+0.250926$ and fails to
resolve positive in four octaves of nine, with standard errors of
$0.033$ to $0.105$ against a global deficit of $+0.134019$ measured to
$0.003$. Its own slope is $+0.014678\pm0.017249$, $t=0.85$: flat.
**Nothing is closing.** Four octaves fail to resolve because seven to
ten points cannot resolve a tenth, not because the demand is met
anywhere.

**A registered rule was misspecified, for the second time.** Z4 was
written "the deficit stays positive by more than two standard errors
in every octave", with the gloss that its refutation would be "the
demand met, at some scale, by the field itself". Unresolved is not
zero, and the gloss named the wrong event — the same error
[rem:thetalaw] recorded for its U4. Twice is a rule, and it is now
M9 in the README: **a refutation rule is written on the quantity, and
what its refutation would mean must be checked against the ways the
condition can fail** — of which "too noisy to tell" is always one.


#### Remark (the sign axis, seven octaves further, and what one octave does to it) {#rem:alphareach}
<!-- evidence: audit_alpha_reach.py -->

This is the axis the proof is stuck on. [rem:leanidentity] leaves one
requirement — $|\sum a|$ must reach the order of its own $\ell^2$
norm, a gap of $+0.134019$ in exponent — and every constructive route
to it is closed: [rem:levelmagnitude] in principle, [rem:filter]
because the inverse filter grows like $1.916413^m$, [rem:meanonly]
because the denominator is a correlation and not a mean. One
computational question is left: does the gap close as $N$ grows?

[rem:alphalocal] answered "nothing is closing" from nine octaves,
with the deficit's slope at $+0.014678\pm0.017249$. That is honest
and weak: a drift two standard errors inside it would close the gap
in about two decades, so the measurement could not exclude the budget
route, it could only fail to see it.

**The field stopped at $1.024\cdot10^8$ because of the sieve, not the
arithmetic.** The published route holds $\Lambda$ at eight bytes an
index with an int32 cofactor beside it, thirteen bytes in all; at
$8\cdot10^9$ that is $96.86$ GB. The kind byte of [rem:rung18] holds
the same information in two — $14.90$ GB — and the half-index trick
of [rem:rung17] does not apply here because $m$ runs over every
integer coprime to $k$, so $N-mk$ takes both parities. Recomputed
through this packing with the per-$k$ sums split into blocks, the top
published octave returns $+0.664374$ and $+0.108604$ to every digit
printed (R1).

**The extension buys the resolution it was run for** (R5). Seven more
octaves, seventy-five more $N$, and the deficit slope's standard
error falls from $0.017249$ to $0.007837$.

**And then one octave decides everything.** Over sixteen octaves the
$\alpha$ slope is $-0.012011\pm0.005910$, $t=-2.03$ — resolved, which
refutes R3 and would overturn [rem:alphalocal]'s Z3, the finding that
$\alpha$ is a constant of the range. The deficit slope is
$-0.007519\pm0.007837$, $t=-0.96$: unresolved, refuting R4 by the
route its rule named (b), and **with its sign flipped** from the
$+0.014678$ nine octaves gave.

Both of those are one octave's doing. Dropping each octave in turn
moves the $\alpha$ slope across $-0.015033$ to $-0.005791$, a spread
of $0.009242$ against its own standard error of $0.005910$; and the
deficit slope across $-0.012726$ to $+0.001402$, a spread of
$0.014129$ **that crosses zero**. The octave responsible sits at mid
$\log N=22.7030$ and holds four points where the others hold eleven
to fourteen. **It is short because the field was cut at
$8\cdot10^9$, inside an octave, and that was a choice made when the
bound was registered.** The bound is not moved here: choosing it
after seeing which side of zero it puts the answer on is the error
this repository exists to avoid.

So the honest report is narrower than either refutation sounds.
$\alpha$ is not shown to drift; a sixteen-octave fit says so at
$t=-2.03$ and a fifteen-octave fit says $-1.13$. The deficit is not
shown to close; it is not shown to grow either. **What the extension
did establish is the size of the blind spot.** A drift smaller than
$0.015673$ is invisible in this field, and a drift that size would
close the gap in $3.7$ decades — where nine octaves could only say
two. The route is not excluded, and it is excluded less loosely than
before by about two decades.

Two things follow for the next measurement. The first is that a field
cut inside an octave buys a partial octave that dominates the
regression, and the next reach should end on a boundary — $LO\cdot
2^j$ — rather than at a round number. The second is that this
analysis summarises each octave into one number and then fits sixteen
numbers, which throws away most of the data: the $87$ per-$N$ values
this run prints as `POINT` markers exist precisely so that a fit over
the points themselves can be done without measuring them again.


#### Remark (the deficit is closing, and where it closes is not quotable) {#rem:deficitdirect}
<!-- evidence: audit_deficit_direct.py -->

[rem:alphareach] left the sign axis with the deficit's drift at
$-0.007519\pm0.007837$ — unresolved, and moved across zero by dropping
any one octave. It also named the reason, and the reason was not the
field's size: **the analysis summarises each octave into one number
and then fits sixteen numbers.** A hundred and fifty-six measurements
become sixteen, and the drift is estimated with thirteen degrees of
freedom. The drift is the coefficient of the quadratic term in
$$\log(|\textstyle\sum a|/\ell^2) = c + \alpha x + \tfrac12\beta x^2,
\qquad x=\log N,$$
and fitting that directly uses every point.

Nothing is measured at large $N$ here. The seventy-five points above
$1.024\cdot10^8$ are read from the `POINT` markers [rem:alphareach]
printed for this purpose; only the eighty-one below are recomputed,
which costs minutes. Twelve are both, and they agree to a relative
$1.376\cdot10^{-11}$ (A1); a line on the published eighty-one returns
the whole-field deficit $0.134019$ exactly (A2).

**The direct fit is twelve times sharper** (A3). Its standard error on
$\beta$ is $0.000632$ against the octave route's $0.007837$ — a ratio
of $12.4$, where the trial on seventy-five points had suggested $5.5$.

**And the drift is resolved** (A4). $\beta=-0.007380$, $t=-11.67$.
**The deficit is closing.** This is the first resolved statement the
sign axis has produced about whether the requirement of
[rem:leanidentity] can be met by pushing $N$, and it says the
direction is the favourable one.

**The two routes agree** (A5): $-0.007380$ against $-0.007519$, a gap
of $0.000139$ inside the octave route's own $0.007837$. So the
sharpening is not a different answer, and the linear-drift model that
buys it is not contradicted by the model-free route.

**Where it closes is a different matter, and the run measures why.**
The fitted local deficit is $+0.247479-0.007380x$, which reaches zero
at $\log_{10}N=14.5635$ with a bracket $[13.5846,\,15.9638]$ from all
$4000$ draws — $4.7$ decades outside the field. And the same fit on
the published eighty-one points alone puts it at $10^{60.4543}$.
**Adding seven octaves moved the forecast by $45.8908$ decades.**
That is not a bracket that has settled; it is a quantity with no
stability at all, and it is exactly what [rem:shapepower] exists to
forbid — a level read off an underived shape extrapolated past its
data. **No closure $N$ is published.**

So the axis has moved, and less far than the resolved $t$ makes it
sound. What is measured, inside the field and by two routes that
agree, is that the deficit shrinks at $0.007380$ per log unit. What
is not measured is that it reaches zero anywhere: a linear drift in
$\log N$ is an assumption the data cannot distinguish from a drift
that flattens, and the forecast's $45.8908$-decade movement is the
evidence that it cannot. The budget route is no longer invisible; it
is visible and it is not thereby open.

Two things this leaves. The eighty-one published points were
recomputed here and their per-$N$ values still live only in prose —
the whole field should carry `POINT` markers, not part of it. And the
partial octave that [rem:alphareach] found dominating its regression
has no effect here, because the direct fit does not weight octaves;
that is a second reason the summary route was costing more than its
robustness was worth.


#### Remark (the drift steepens, and two of this run's own rules were misspecified) {#rem:deficitshape}
<!-- evidence: audit_deficit_shape.py -->

[rem:deficitdirect] resolved the deficit's drift at $t=-11.67$ and
declined to say where it reaches zero, because that forecast had
moved $45.8908$ decades on the last seven octaves. The reason is a
shape assumption: a quadratic in log-log makes the drift constant, so
the deficit must cross zero; a cubic lets the drift die, and then the
deficit approaches a positive limit and **never closes at any $N$**.
Both fit a shrinking deficit; only the first implies the budget route
exists. Nothing is measured here — all $156$ points are read from
`POINT` markers.

**The drift is not constant, and it steepens** (B2). The cubic
coefficient is $-0.00201923\pm0.00068284$, $t=-2.96$: resolved, and
negative, which means the deficit is closing *faster* than the
quadratic said rather than slowing. The r.m.s. residual falls from
$0.032663$ to $0.031762$, and walking forward from the fortieth point
the cubic's out-of-sample departure is $0.028534$ against the
quadratic's $0.029139$ (B3) — so it predicts as well as it fits,
which is the test [rem:signrun] applied and a term can fail.

**And two registered rules named the wrong event.** B4 asked that a
nominal five per cent test come in *below* five per cent; it measured
$0.0515$, which is $0.44$ Monte Carlo standard errors from $0.05$ on
$4000$ draws. That is a coin flip, and a rule whose outcome is a coin
flip tests the coin. B5 asked whether the cubic *resolves* under the
alternative — but the alternative's own curvature is $+0.00022008$
against a standard error of $0.00068284$, $0.32$ of it. Nothing could
resolve that, and the rule should not have asked it to. **Both stand
refuted as written.**

The question B5 should have asked is not whether the alternative
would be detected but whether the measurement can be told apart from
it, and it can: the gap between the measured cubic term and the
alternative's is $-0.00223931$, **$3.28$ standard errors**. So the
particular way of never closing that this run constructed — a drift
dying exactly where the quadratic put the closure — is excluded, and
excluded by the direct estimate rather than by a power calculation
that was never able to speak.

**This is the third time.** [rem:thetalaw]'s U4 and [rem:alphalocal]'s
Z4 made the same error and it became M9 in the README: *a refutation
rule is written on the quantity, and what its refutation would mean
must be checked against the ways the condition can fail.* M9 has now
failed to prevent its own third instance, which is worth recording as
a fact about rules rather than about these fits.

What this settles and what it does not. The drift is measured, it
steepens, and one construction of "never closes" is out by $3.3$
standard errors. That is not "the deficit closes": a drift that dies
later than the quadratic's crossing, or dies in some shape other than
a cubic's, is untouched, and the field's $10.5966$ in $\log N$ is not
long enough to enumerate those. **The budget route is narrower than
[rem:deficitdirect] left it and still not open**, and what would open
or close it is a functional form for the drift, which no measurement
in this repository has ever supplied.


#### Remark (the never-closes region has a boundary, and the family has none) {#rem:deficitregion}
<!-- evidence: audit_deficit_region.py -->

[rem:deficitshape] excluded one point of the "never closes" region at
$3.28$ standard errors. The region has a boundary and it is simpler
than that point. Writing the deficit as the cubic's derivative,
$$\mathrm{deficit}(x) = c_1 + c_2 x + \tfrac12 c_3 x^2,$$
a parabola in $x=\log N$: if $c_3<0$ it opens downward and the deficit
reaches zero at some finite $N$ whatever else is true; if $c_3>0$ its
vertex sits at $x=-c_2/c_3$, which with $c_2>0$ is at negative $x$,
below the field entirely, so past the top the deficit only rises and
never reaches zero. **Inside the cubic family the whole region is the
single inequality $c_3>0$.**

The premises hold (C2): the deficit at the top of the field is
$+0.057813$, $c_2$ is $+0.02820735$, and a positive $c_3$ would put
the vertex at $x=-13.97$ against the field's lowest point $12.2061$.
And the inequality is excluded (C3): $c_3=-0.00201923\pm0.00068284$,
one-sided $t=-2.96$.

**And then the family fails** (C4, refuted). A quartic term added to
the same fit is $+0.0025976491\pm0.0009936619$, $t=2.61$ — resolved,
and positive, which is the direction that bends the deficit back up
beyond the field and puts the region straight back. The refutation
rule named this as the outcome that would cost the most and said it
was not the one predicted. **So C3's exclusion is void, and the
parabola argument does not describe the deficit past the field.**

**A diagnostic, run after that refutation and predicted by nothing.**
C4 asked whether one more degree resolves; it does, and the question
that raises is whether the next ones do. Fitting degrees two through
eight, the newest coefficient's $t$ reads $-11.67$, $-2.96$, $+2.61$,
$-1.34$, $-1.16$, $+2.67$, $+7.41$, while the r.m.s. residual moves
only from $0.032663$ to $0.030092$ — eight per cent across six added
parameters. **Degree eight resolves at $t=7.41$ and buys almost
nothing.** That is not a deficit that is a degree-eight polynomial;
it is a family flexible enough to chase a shape it does not contain,
and it means no polynomial statement about this deficit past the
field is stable at any degree.

So the computational branch of the sign axis closes, and it closes
without an answer. [rem:deficitdirect] measured the drift and it is
real: inside the field the deficit shrinks, resolved at $t=-11.67$,
by two routes that agree. [rem:deficitshape] measured that the drift
steepens. Neither licenses the next step, because every attempt to
say what happens past $\log N=22.7030$ has to name a family, and the
degrees say this field will resolve a new coefficient in whatever
family is offered without the residual conceding that the shape was
found. **[rem:shapepower] said this for the level axis. It is now
measured for the sign axis, with the degree at which it bites.**

What is left is not more $N$. Extending the field lengthens the
lever arm for every degree at once, so a longer field resolves higher
degrees rather than settling lower ones — the degree-eight $t$ of
$7.41$ on ten and a half decades is what that looks like already.
What would settle it is a derivation that says which family the
deficit belongs to, and no measurement in this repository has ever
supplied one.


#### Remark (the target of the sign axis is derived, not fitted) {#rem:targetderived}
<!-- evidence: audit_target_derived.py -->

[rem:deficitregion] closed the computational branch of the sign axis
without an answer and named what would open it: *a derivation that
says which family the deficit belongs to, and no measurement in this
repository has ever supplied one.* The deficit is a difference of two
exponents, $e(\ell^1/\ell^2)-e(G)$, and [rem:leanidentity]
quotes the second as a measured $+0.283586$.

[rem:leanidentity]'s W4 wrote the ingredient down and read it
backwards. It observed $\ell^1\le\sqrt{\#k}\,\ell^2$ — Cauchy–Schwarz,
an identity of norms — and then called W4 refuted because the measured
exponent stood $3.15$ standard errors *above* the ceiling $\theta/2$.
**A ratio bounded by one cannot exceed its ceiling asymptotically.**
The excess is a transient, which is the argument [rem:headfraction]
made about the head's share and which was never applied here.

**And $\#k$ is not empirical.** The range is the squarefree $k<N^
\theta$ coprime to $N$, so
$$\#k=\frac{6}{\pi^2}\prod_{p\mid N}\frac{p}{p+1}\;N^{\theta}\,
(1+o(1)),$$
leading order derived. Measured against that formula at ten $N$ the
worst error is **$0.370$ per cent** (E2), so the count carries no free
parameter. And $\theta=0.56$ is a parameter of the *construction* —
it defines the $k$-range — **not the level exponent the program is
proving**, so the ceiling $\theta/2=0.28$ is an exactly known number
that does not move with the target.

| $N$ | $\#k$ | derived | $\ell^1/\ell^2$ | $r$ |
|---|---|---|---|---|
| $2\cdot10^5$ | $313$ | $314.2$ | $11.9596$ | $0.6760$ |
| $1.28\cdot10^7$ | $3226$ | $3225.6$ | $39.6783$ | $0.6986$ |
| $1.024\cdot10^8$ | $10338$ | $10335.9$ | $69.7186$ | $0.6857$ |

$r=(\ell^1/\ell^2)/\sqrt{\#k}$ stays in $[0.6622,\,0.6986]$, below one
everywhere (E3), and it has stopped climbing: the bottom half of the
field rises $+0.01120$ and the top half **falls** $-0.00150$ (E4). The
local exponent of $\ell^1/\ell^2$ on sliding windows of five runs
$+0.285443$, $+0.289948$, $+0.288106$, $+0.290289$, $+0.284047$,
$+0.277325$ — **the top window is below the ceiling**, at $-0.002675$
(E5). The remaining rise available above $(\theta/2)\log N$ is
$-\log r=0.3773$ (E6).

**So the demand of [rem:leanidentity] is a fixed number.** It stops
being "$e(G)$ must reach a fitted $+0.283586$" and becomes
$$e(G)\longrightarrow \theta/2=0.28 ,$$
and the target side of the deficit needs no shape at all. What happens
past the field is a question about $G$ alone. That is one of the two
free terms [rem:deficitregion] complained of, removed by derivation
rather than by fitting.

*Corrected.* This remark was written saying "six cycles have treated
both as fitted", and that is wrong: [rem:denominator] had already
written the arithmetic form, "$e(G)\to\theta'/2$ is exactly
$\alpha\to e(\ell^1)-\theta'/2$", and called the demand arithmetic
rather than fitted. Its X1 even measured the gap $1-\alpha=+0.282084$
against $\theta'/2=0.2800$. **The form was in the repository and this
run did not find it.** What is this run's own is narrower and is what
the rest of the remark reports: that $\#k$ matches its derived leading
order to $0.370$ per cent so the ceiling carries no fitted quantity,
that W4's "refuted" reading of a bounded ratio was the wrong reading,
that the local exponent is measured falling below the ceiling, and
that $\theta$ is a parameter of the $k$-range rather than the level —
so the ceiling does not move with what the program is proving, which
[rem:denominator]'s $\theta'$ notation leaves open.

**What this does not do, stated at its own strength.** E5's registered
bar was that the fall exceed the window error; it does, $+0.008119$
against $0.006187$, but that is $t=1.31$ and **would not be resolved
at the bar this repository normally uses**. The bar was set weak and
is recorded as weak. The unconditional half is the inequality: $r\le1$
gives $e(\ell^1/\ell^2)\le\theta/2$ with no assumption, while equality
needs $r$ bounded away from zero, which is measured over ten $N$ and
not proved. And $r$ turning over in the top half means the picture is
not the monotone rise the setup described — the ceiling holds either
way, but "bounded and rising to a limit" is not what the field shows.
**Nothing here measures $G$**, and the shape of $G$ past the field is
exactly the thing [rem:deficitregion] said no measurement supplies.

#### Remark (this field cannot tell a logarithm from a line) {#rem:deficitlog}
<!-- evidence: audit_deficit_log.py -->

[rem:deficitregion] closed the sign axis' computational branch on one
pathology — the polynomial family resolves a new coefficient at any
degree offered, degree eight at $t=+7.41$, while the r.m.s. residual
moves only $0.032663\to0.030092$ — and named the cure: *a derivation
that says which family the deficit belongs to.* A power of $\log N$
is exactly what does that to a polynomial fit, and there is a reason
to expect one. The diagonal of $\ell^2$'s second moment is the shape
prop:V evaluates exactly; if it carries $\ell^2$, then $\ell^2$ is
$\sqrt N$ times a power of $\log N$ and **the $+0.583897$ six remarks
quote is a logarithm absorbed into a power**. If $\lvert\sum a\rvert$
is too, the $\sqrt N$ cancels and
$$y(x)=\log\bigl(\lvert\textstyle\sum a\rvert/\ell^2\bigr)
 = c + C\log x$$
with no $x$ term at all — two parameters where the cubic has four.
The form even carries a check available before fitting: it forces
$y''/y'=-1/x$, and [rem:deficitdirect]'s $-0.007380$ and $0.134019$
give $x=18.16$, inside the field.

**It was the coincidence.** The run reads all $156$ points from POINT
markers and measures nothing new; the gate reproduces $0.134019$ and
$\beta=-0.007380$ exactly (F1). Then:

| | |
|---|---|
| log form, 2 parameters | r.m.s. $0.033651$ |
| line, 2 parameters | $0.044911$ |
| cubic, 4 parameters | $0.031762$ |

**F2 is refuted**: the two-parameter log form does not match the
cubic's residual. **F3 is refuted and by more**: the form implies a
deficit of $C/\bar x=2.030302/18.0518=0.112471$ against the field's
own $0.134019$, low by $16.08$ per cent. The pure $c+C\log x$ does not
reproduce the slope it has to reproduce.

**But the verdict this run is entitled to is narrower than either
refutation, and it was registered in advance.** Over this field $x$
runs $12.2061$ to $22.8027$ while $\log x$ runs $2.5019$ to $3.1269$,
and their correlation is
$$\mathrm{corr}(x,\log x)=0.996737 .$$
G69's threshold is $0.99$, the marker `COEFF NOT SEPARABLE` is
printed, and the rule fixed before the run says **F2, F4 and F5 may
not be read as a win for either family**. So the log family's $t=-0.43$
on an added $x$ term — against the polynomial family's $+7.41$ at
degree eight — is not evidence that it contains the shape, and its
walk-forward departure of $0.023586$ against the cubic's published
$0.028534$ is not either. Both are recorded and neither is read.

**What this establishes is a number for the reach, and it sharpens
[rem:deficitregion] rather than answering it.** That remark said a
longer field resolves higher degrees rather than settling lower ones.
This one says something worse about the same field: at
$\mathrm{corr}=0.996737$ **a logarithm and a line are the same
regressor here**, so no fit performed on these $10.6$ decades can
distinguish a family that closes from a family that does not, whatever
its motivation. The discrimination is not in the data, and picking a
better-motivated shape does not put it there — the same limit
[rem:curvereach] measured on the level axis, now measured on this one.
[rem:shapepower] is not repealed by a derivation.

#### Remark (the object has a name: a truncated Möbius sieve weight) {#rem:sieveweight}
<!-- evidence: audit_sieve_weight.py -->

OPEN item 5 ends its description of what the sign axis is left with by
saying the object **has no name**, and [rem:deficitlog] then measured
that this field cannot separate a logarithm from a line, so the shape
of $G$ cannot be settled by fitting and what remains is an
unconditional statement about $\lvert\sum a\rvert$. An unconditional
statement needs an object such statements exist about. The algebra
supplies one.

[rem:denominator] has $\sum a=\sum_j\Lambda(N-j)\Lambda_K(j)$ and
[rem:support] splits it. For **squarefree** $j$ coprime to $N$ every
condition on $k$ is automatic and the untruncated $\mu*\log$ is
$\Lambda(j)=0$ once $\omega(j)\ge2$, so $\Lambda_K(j)$ is exactly
minus what the truncation threw away; rewriting those divisors by
their cofactors $d=j/k$,
$$\Lambda_K(j)=-\sum_{d\mid j,\ d\le j/K}\mu(d)\log(j/d),
\qquad K=\lfloor N^{\theta}\rfloor,$$
the Eratosthenes–Legendre sieve weight cut at $D_j=j/K$.

**It holds, at machine precision** (H6): the worst relative departure
is $3.234\cdot10^{-15}$ at each of $N=20000$, $50000$, $100000$, over
$4494$, $11756$ and $24181$ such $j$. The rebuild is exact (H8) — the
prime part, the sieve weight where it applies, and $\Lambda_K$ itself
on the handful of $j$ the derivation does not cover, sum to $\sum a$
at relative $0$, $1.29\cdot10^{-16}$, $0$. **And the name covers the
object** (H9): the squarefree composite $j$ carry $1.0001$, $1.0000$,
$1.0001$ of the composite part, which is itself the whole of $\sum a$
(H5). The gate reproduces $\lvert\sum a\rvert=87895.3236$ at
$N=200000$ (H1).

**Three predictions were refuted first and stay refuted.** H2 wrote
the level as the real $N^{\theta}$; the repository's $k$-range is
`range(2, int(N**theta))`, so the truncation is at $\lfloor
N^{\theta}\rfloor$. At $N=50000$ those are $427$ and $427.979732$, and
$j=16653=3\cdot7\cdot13\cdot61$ has the divisor $427$ in the gap — one
$j$ in $11756$, and the identity failed by a relative $5.513$ there.
**The level being an integer is a fact about the construction, not a
rounding convenience**, and H2 as written did not know it. H3 and H4
ranged over every $j$ with $\Lambda(j)=0$, which includes $j=12$ and
its kind, where $\mu*\log$ is not $\Lambda$ and the rewriting was
never claimed: the test set was wrong, not the object.

**H7's cap could not be met by any computation.** It asked that a sum
of floating-point magnitudes be *exactly* zero; the sums are
$6.2172\cdot10^{-15}$, $1.3323\cdot10^{-14}$, $2.1316\cdot10^{-14}$,
which is zero to machine precision and is not zero. It stands refuted
as written and the defect is the cap. That is the same family as
[rem:identityforced]'s Q4 — a rule whose outcome no run could change —
and its second instance here.

**And this run's own closing sentence is wrong.** It was written to
fire when H7 failed and says the object "is not named here and item 5
keeps the description it had". H6, H8 and H9 say otherwise and H7's
failure is a floating-point cap. The support claim H7 was testing is
true and visible in the numbers it printed: below the level a
squarefree $j$ has no divisor above the level to throw away, so
$\Lambda_K$ vanishes there, and the composite part is supported on
$j\ge K$.

**What this buys, and no more.** $\lvert\sum a\rvert$ is a correlation
of $\Lambda$ against a classical sieve weight, supported on $j\ge K$ —
a family that unconditional bounds exist for, which is what item 5
now needs and what "no name" was blocking. **An identity is not a
bound.** Nothing here supplies an estimate, nothing says the classical
bounds are strong enough to decide item 5, and no exponent or forecast
is measured. [rem:shapepower] and [rem:deficitlog] are untouched.

#### Remark (the sum has two square-root barriers and they differ by a factor of two) {#rem:jbarrier}
<!-- evidence: audit_jbarrier.py -->

[rem:sieveweight] named the object. Item 5's demand is
$\lvert\sum a\rvert\lesssim\ell^2$, and $\ell^2$ is a norm **over
$k$** while $\sum a$ is a sum **over $j$**. One number, two groupings,
two square-root barriers — and only one had ever been computed:
$$\ell^2=\Bigl(\sum_k(\log k)^2H(N;k)^2\Bigr)^{1/2},\qquad
D=\Bigl(\sum_j\Lambda(N-j)^2\Lambda_K(j)^2\Bigr)^{1/2}.$$
[rem:denominator]'s "already as small as chance allows" is a statement
about the $k$ grouping; nothing here had said what chance allows in
the $j$ one.

Nine $N$ over $2.41$ decades, gate reproducing $87895.3236$ and
$11.9596$ (J1):

| quantity | exponent | s.e. |
|---|---|---|
| $\lvert\sum a\rvert$ | $+0.763327$ | $0.009864$ |
| $\ell^2$ over $k$ | $+0.625318$ | $0.005519$ |
| $D$ over $j$ | $+0.625057$ | $0.003179$ |

The sum stands above both by the same amount — $+0.138271$ against
$+0.138010$ (J3) — and above them at every $N$, not on average (J4).

**J2's reading is barred by J2's own rule.** The gap $-0.000261$ sits
inside the larger individual error $0.005519$, and the rule fixed
before the run says that in that case neither consequence may be
drawn. It is not drawn.

**The correction this run then registered failed, and failed
informatively** (K1). The rule had compared a difference to an error
bar the difference does not have; refitting $\log(D/\ell^2)$ directly
gives $-0.000261\pm0.003787$, $t=-0.07$ — unresolved, the same
verdict — but the error is *above* $D$'s own $0.003179$. **So the two
fits' errors are not strongly correlated**, which was the premise of
the correction, and K1 is refuted on that clause. The exponent route
cannot decide this question and neither of its two statistics does.

**The ratio says what the exponents could not** (K2, K3). $D/\ell^2$
runs

$$0.5090,\ 0.4959,\ 0.4942,\ 0.4900,\ 0.5106,\ 0.4933,\ 0.4853,\
0.5114,\ 0.4998$$

— mean $0.498823$, range $5.23$ per cent, flat across the field. **The
$j$-side barrier is half the $k$-side barrier**, and a flat ratio is a
stronger statement than equal exponents.

**Its evidential weight is low and the reason is disclosed.** K2 and
K3 were registered after the first run printed those ratios, so they
test a reading of data already seen; K3's cap of $0.02$ around a half
was nearly certain to be met once the range $0.485$–$0.511$ was on the
page. **They are confirmations, not blind predictions.** No mechanism
is offered for the constant, as K3's rule requires — the diagonal
sketch that suggests one gives $6/\pi^2$ in the wrong direction, not
a half.

*Withdrawn.* A blind second radical family refutes the constant and the
drift within that family resolves at $t=-7.90$: the two barriers do not
share a scale. See [rem:jbarrierreach]. The paragraph below is what this
remark suggested and it is no longer available; J1 to J4 and the
exponents at family A stand.

**What is suggested and not established.** If $D/\ell^2$ really is a
constant, then $\ell^2$ is not the wrong floor for a sum indexed by
$j$ — it is the right floor up to a factor of two, and item 5's demand
is exactly "achieve square-root cancellation", a named barrier rather
than a gap between two fitted exponents. This field is consistent with
that and does not establish it. What would is a second radical family
and a longer field, on the ratio rather than on the exponents. And
$D$ is a heuristic floor, not a theorem: it is what the sum would be
under random signs in $j$, and no lower bound on $\lvert\sum a\rvert$
follows from measuring it.

#### Remark (the half was the radical's, and the barriers do not share a scale) {#rem:jbarrierreach}
<!-- evidence: audit_jbarrier_reach.py -->

[rem:jbarrier] found $D/\ell^2$ flat at $0.498823$ across $2.41$
decades, said its own evidence was weak — K2 and K3 were registered
after that ratio was on the page — and named the test that would
settle it: a second radical family and a longer field, on the ratio.
This is that test, blind. **It refutes the constant.**

| family | radical | $\prod_{p\mid N}\frac{p}{p+1}$ | mean $D/\ell^2$ |
|---|---|---|---|
| A | $\{2,5\}$ | $0.555556$ | $0.503442$ |
| B | $\{2,3,5,7,11,13\}$ | $0.310330$ | $0.178268$ |

The gate reproduces all nine published ratios to four decimals (L1)
and family A extended two octaves keeps its mean at $0.503442$ (L2).
**L3 and L4 are refuted**, the second by $0.325173$ against a cap of
$0.02$. L4's rule named this outcome in advance as the one that would
matter most, and it fires: the relation between the two barriers
carries arithmetic, and a single number is the wrong shape for it.
It is not the arithmetic factor either — the means stand in ratio
$2.82$ while the factors stand in $1.79$ — and **no mechanism is
offered.**

**And it is not a constant within family B at all.** The ratio there
runs $0.2303$, $0.2055$, $0.1899$, $0.1753$, $0.1711$, $0.1642$,
$0.1577$, $0.1549$, $0.1555$ — falling throughout.

**The half's proposed mechanism fails too.** The diagonals of $D^2$
and $\ell^{2\,2}$ are the same quantity by algebra, so a ratio of a
half would mean $D^2=\mathrm{DIAG}$ and $\ell^{2\,2}=4\,\mathrm{DIAG}$.
L6 is refuted by $475$ per cent: $\ell^{2\,2}/\mathrm{DIAG}$ has mean
$23.00$ and runs from $1.76$ to $53.96$. **L5 "holds" and must not be
read as support**: its cap was on the mean of a series running $0.48$
to $1.80$, and a mean of $1.01065$ from that spread tests the cap and
not the claim. That is the same defect as [rem:identityforced]'s Q4
and [rem:sieveweight]'s H7 — a rule that no outcome could have
failed — and its third instance.

**The statistic the caps should have been written on** (a diagnostic,
run after the verdicts and predicted by nothing): the drift of
$\log(D/\ell^2)$ against $\log N$ within each family.

| family | drift | s.e. | $t$ |
|---|---|---|---|
| A | $+0.005736$ | $0.003313$ | $+1.73$ |
| B | $-0.068713$ | $0.008699$ | $-7.90$ |

**The two barriers do not share a scale.** Family B resolves the
drift at $t=-7.90$, and family A's agreement — the whole of
[rem:jbarrier]'s finding — is one radical's coincidence over a field
too short to see its own drift.

**So [rem:jbarrier]'s suggested reading is withdrawn.** It offered
that $\ell^2$ might be the right floor for a $j$-indexed sum up to a
factor of two, making item 5's demand a named barrier. It is not: the
$j$-side barrier is a different size in a different family and falls
away from $\ell^2$ as $N$ grows. What survives from that remark is
what it measured and not what it suggested — J1 through J4, the
exponents at family A, and $\lvert\sum a\rvert$ standing above both
barriers at every $N$ there. **Which floor item 5 is measured against
is open, and this run closes one candidate answer rather than
supplying another.**

#### Remark (the deficit carries the radical) {#rem:whichfloor}
<!-- evidence: audit_which_floor.py -->

[rem:jbarrierreach] left the question of which floor item 5's demand
is measured against, and handed over the tool to answer it. On family
A the two barriers sit on top of each other — $D/\ell^2$ drift
$+0.005736$ at $t=+1.73$ — so that field cannot tell. **Family B
separates them** at $t=-7.90$, and a sum can track at most one of two
floors that are pulling apart.

$\lvert\sum a\rvert$ on family B had never been computed. It is the
half of this run that decides, and it was blind.

| family | radical | $\lvert\sum a\rvert/\ell^2$ drift | s.e. | $t$ |
|---|---|---|---|---|
| A | $\{2,5\}$ | $+0.141100$ | $0.005708$ | $+24.72$ |
| B | $\{2,3,5,7,11,13\}$ | $+0.284941$ | $0.001449$ | $+196.70$ |

**M2 is refuted, and it is the outcome its own rule named as costing
the most.** The difference is $-0.143841$ against an error of the
difference of $0.005889$ — twenty-four times it, not a near miss and
not the unresolved case this run named in advance. **The deficit
carries the radical.**

**The ranges do not explain it** (a diagnostic, run after the verdicts
and predicted by nothing). The two families span different $\log N$
and [rem:deficitdirect] measured the deficit's drift as itself
drifting, so the windows were the obvious confound. Refitted on the
$\log N$ the two share, $10.3100$ to $15.8551$, family A gives
$+0.138586\pm0.010916$ against family B's $+0.284941\pm0.001449$: a
difference of $-0.146355$ at $t=-13.29$. It is larger on the shared
window, not smaller.

**What this costs.** The published field is $168$ $N$ and **every one
of them has radical $\{2,5\}$** — one family, family A's. So item 5's
$+0.134019$ is that family's number. Family A measured here gives
$+0.141100$ over eleven doubling points, consistent with it; family B
gives twice as much. Every remark that quotes $+0.134019$ as *the*
distance the demand has to travel is quoting a number of the
primorial-free family, and on a primorial family the sum runs away
from $\ell^2$ at twice the rate. **The demand is not one number.**

M3 holds — the $\lvert\sum a\rvert/D$ drifts differ by $-0.218290$ —
and it was registered only so the pair would be read together, since
the two differ by the $D/\ell^2$ drift already published. M4 holds:
$\lvert\sum a\rvert$ stands above both barriers at every $N$ of both
families. **Neither floor is the one the deficit is flat against**,
so the question [rem:jbarrierreach] left is not answered here — it is
made harder, because the answer now has to explain a rate that
depends on which primes divide $N$.

**And the gate failed first, on a defect of its own.** M1 read
$\lvert\sum a\rvert$ from a table printed to two decimals and judged
it at four; it came out refuted by the print bound rather than by any
disagreement — `TOL BELOW PRINT`, which G75 exists for. The script
exits at M1, so no verdict on M2 to M4 existed when this was found,
and the gate now reads markers carrying the digits it asks for with
the tolerance unchanged. Worth recording beside it: this run's care
went into naming M2's unresolved case *before* the fact, which is
where the previous three failures had been, and the defect moved to
the gate instead.

*Reading withdrawn.* A blind, pre-registered control at a radical this
branch had never used finds five bases sharing one radical disagreeing
beyond their own errors, at a chi-square of 14.233 against a cap of 9.49.
The drift is therefore not a function of the radical, and this remark read it
from one base per radical. **Its numbers stand and the reading built on
them does not.** See [rem:radicalblind].

#### Remark (which primes, not how many) {#rem:radicallaw}
<!-- evidence: audit_radical_law.py -->

[rem:whichfloor] found the deficit's drift $24$ errors apart on two
radicals and ended by saying two points do not determine a function.
Six doubling families here, each with a fixed radical, and **three of
them share $\omega=2$ and differ in which primes** — the one thing a
two-family measurement cannot ask. The gate reproduces family B's
drift to six decimals and $\lvert\sum a\rvert$ at $N=200000$ to four
(N1).

| radical | $\omega$ | $\prod\frac{p}{p+1}$ | drift of $\log(\lvert\sum a\rvert/\ell^2)$ |
|---|---|---|---|
| $\{2\}$ | $1$ | $0.666667$ | $+0.068281$ |
| $\{2,5\}$ | $2$ | $0.555556$ | $+0.138010$ |
| $\{2,7\}$ | $2$ | $0.583333$ | $+0.141444$ |
| $\{2,3\}$ | $2$ | $0.500000$ | $+0.212384$ |
| $\{2,3,5\}$ | $3$ | $0.416667$ | $+0.265357$ |
| $\{2,3,5,7,11,13\}$ | $6$ | $0.310330$ | $+0.284941$ |

**N3 is refuted and it decides the shape of the answer.** The three
$\omega=2$ families span $0.074374$ — against the largest single drift
error of $0.019289$, so the unresolved clause this run named in
advance does *not* fire and the disagreement is real. $\{2,3\}$ sits
half again above $\{2,5\}$ and $\{2,7\}$, which agree with each other.
**No function of $\omega$ can be the answer. It is which primes divide
$N$, not how many.**

N2 and N4 are refuted too, and **their readings are barred by the rule
that named this case**: a line through a variable already known not to
be the cause is a fit, not a law. For the record the residuals were
$0.044185$ in $\omega$, $0.019112$ in the arithmetic factor and
$0.052354$ in $\log\mathrm{rad}$, and none of the three is read.

**And the primes do not act one at a time** (a diagnostic, after the
verdicts and predicted by nothing). Taking the two-prime families as
each odd prime's own contribution — $3$ gives $+0.144102$, $5$ gives
$+0.069729$, $7$ gives $+0.073163$, on a base of $+0.068281$ — the
multi-prime families are then predictions:

$$\{2,3,5\}:\ +0.282112 \text{ against } +0.265357,\quad
\{2,3,5,7,11,13\}:\ +0.355275 \text{ against } +0.284941 .$$

Both fall short, by $-0.016756$ and $-0.070333$, and the second
shortfall is understated because the additive prediction there omits
$11$ and $13$ entirely. **The contributions are sub-additive and this
run does not say what they are instead.** What it does say is that
$3$ carries about twice what $5$ does while $5$ and $7$ are not
separated by their own errors — a shape, not a law, and no mechanism
is offered for it.

**What this costs item 5.** The demand's distance is not one number
and it is not a number of $\omega$ either. Any statement of the form
"the sum must travel $+0.134019$" is a statement about $N$ with
radical $\{2,5\}$; on $\{2,3\}$ it is half again as far and the
primes' contributions do not combine by any rule measured here.
[rem:shapepower] applies to a law in the primes exactly as it does to
one in $\log N$, and none is claimed.

*Reading withdrawn.* A blind, pre-registered control at a radical this
branch had never used finds five bases sharing one radical disagreeing
beyond their own errors, at a chi-square of 14.233 against a cap of 9.49.
The drift is therefore not a function of the radical, and this remark read it
from one base per radical. **Its numbers stand and the reading built on
them does not.** See [rem:radicalblind].

#### Remark (three is the whole of it, and there is no decay) {#rem:primecontrib}
<!-- evidence: audit_prime_contrib.py -->

[rem:radicallaw] left $f(p)$, the amount a single odd prime dividing
$N$ adds to the drift of $\log(\lvert\sum a\rvert/\ell^2)$, with three
points and no shape. Neighbouring primes cannot be ordered at this
precision — $f(5)$ and $f(7)$ differ by $0.0034$ against errors near
$0.010$ — so this run does not try. It spreads $p$ from $3$ to $101$,
where a $1/p$ decay would put $f(101)$ at a thirtieth of $f(3)$, and
asks for the exponent. The gate reproduces all four shared drifts to
six decimals (O1).

| $p$ | $3$ | $5$ | $7$ | $11$ | $13$ | $17$ | $23$ | $47$ | $101$ |
|---|---|---|---|---|---|---|---|---|---|
| $f(p)$ | $.1441$ | $.0697$ | $.0732$ | $.0745$ | $.0392$ | $.0339$ | $.0090$ | $.0560$ | $.0512$ |

**O2, O3 and O4 are all refuted.** $f(101)/f(3)=0.3555$ against a cap
of a quarter; the slope of $\log f$ on $\log p$ is $-0.317171\pm
0.234902$, which stands $2.91$ errors from $-1$; and the r.m.s.
residual is $0.642640$ against a cap of $0.35$. **The cost of
excluding $p$ from the $k$-range is not of order $1/p$**, so the
natural derivation — that the excluded $k$ are a $1/(p+1)$ share of
the squarefree ones — does not describe it. The sequence is not even
monotone: it bottoms at $p=23$ and comes back up.

**What the numbers say instead** (a diagnostic, after the verdicts and
predicted by nothing). Every $f(p)$ is a drift minus the *same* base
drift, whose error $0.019289$ is the largest here and is common to all
of them, so it cancels in differences and must not be used when
comparing them. Compared properly, **$f(3)$ stands apart from every
other prime**: the differences run $+0.0696$ to $+0.1351$ at
$t=+3.94$ to $+9.12$. The eight primes from $5$ to $101$ have weighted
mean $+0.054986$ with $\chi^2=19.04$ on seven degrees and a largest
deviation of $2.68$ of its own error — **higher than a constant
comfortably allows, but the departure is not a decay**: the low point
is $p=23$ and $47$ and $101$ are above it.

**And this run's closing sentence is wrong.** It was written to fire
when O2 failed and says that reading the two-prime families as one
contribution each "is wrong, not merely unshaped", withdrawing
[rem:radicallaw]'s contributions as a description. That does not
follow: a contribution that fails to decay is still a contribution,
and nothing here tests additivity — [rem:radicallaw] tested it and
found it sub-additive, which stands. **What is refuted is the decay
and the $1/(p+1)$ derivation, not the reading.**

**Where this leaves item 5.** The demand's distance depends on which
primes divide $N$; the dependence is carried almost entirely by
whether $3$ divides $N$, and among the primes above $3$ this field
measures scatter and no shape. A quantity that jumps at one small
prime and is flat-to-ragged afterwards is not a function anyone here
can extrapolate, and [rem:shapepower] applies to it as to the rest.
The design cannot separate $p$ from the base's own size or its
$2$-adic valuation, one base per prime being all it has.

*Reading withdrawn.* A blind, pre-registered control at a radical this
branch had never used finds five bases sharing one radical disagreeing
beyond their own errors, at a chi-square of 14.233 against a cap of 9.49.
The drift is therefore not a function of the radical, and this remark read it
from one base per radical. **Its numbers stand and the reading built on
them does not.** See [rem:radicalblind].

#### Remark (the control this branch never ran, and the clause that barred it) {#rem:basecontrol}
<!-- evidence: audit_base_control.py -->

Writing this run made a larger gap plain than the one
[rem:primecontrib] named. **Every measurement in this branch has used
one base per radical.** [rem:whichfloor] compared $25000$ against
$30030$; [rem:radicallaw] six bases with six radicals;
[rem:primecontrib] nine more. In none of them was the same radical
measured twice, so "the drift depends on which primes divide $N$" was
never separated from "the drift depends on the base", and three
remarks read a radical dependence off a design that could not tell.

Six bases on $\{2,3\}$ spanning $2$-adic valuation $2$ to $12$, three
on $\{2,5\}$, ten $N$ each. The gate reproduces both published drifts
to six decimals on the nine $N$ that produced them (P1).

| radical | drifts |
|---|---|
| $\{2,3\}$ | $.211888\ .213528\ .211865\ .205821\ .204441\ .208169$ |
| $\{2,5\}$ | $.144190\ .138572\ .153405$ |

P2, P3 and P4 all hold — ranges $0.009087$ and $0.014833$ against a
cap of $0.02$, and the two radicals stay $+0.063896$ apart against a
cap of $0.05$.

**And P2's reading is barred by P2's own clause.** The clause said
that a range below the largest single drift error means the control
is underpowered rather than passed. The range is $0.009087$ and that
error is $0.012752$, so it fires. **The control is not claimed to
have passed.**

**The clause used the wrong error, which is the sixth time in this
branch.** It said "the largest single drift error" without saying
within which radical, and the largest — $0.012752$, at base $31250$ —
is in the *other* group. The $\{2,3\}$ group's own largest error is
$0.007815$, below its range. The rule is not rewritten; it fires as
written and the reading stays barred.

**The statistic it meant to make** (a diagnostic, after the verdicts
and predicted by nothing): each group against one constant, on its own
errors. $\{2,3\}$ gives a weighted mean of $+0.209772$ with
$\chi^2=2.512$ on five degrees; $\{2,5\}$ gives $+0.142637$ with
$\chi^2=1.110$ on two. Both are comfortably consistent with one
constant per radical, and neither shows any trend against $2$-adic
valuation — the six $\{2,3\}$ bases run $v_2=6,8,2,10,4,12$ and their
drifts do not order with it.

**What this leaves.** The diagnostic says what the control was built
to ask and answers it the favourable way, but it is a statistic chosen
after seeing the numbers and **does not replace one chosen before
them**. So the position is exactly this: the drift being a function of
the radical is now *consistent with* a direct test at two radicals and
nine bases, and is still not *established* by a pre-registered one.
[rem:whichfloor], [rem:radicallaw] and [rem:primecontrib] keep their
numbers and their readings remain uncontrolled until a run registers
the $\chi^2$ test in advance — on a radical none of these used, so
that it is blind as well as correct.

#### Remark (the control, blind and registered, and it fails) {#rem:radicalblind}
<!-- evidence: audit_radical_blind.py -->

[rem:basecontrol] ran the control this branch had never run and then
could not claim it: its clause fired on an error from the other
radical group, and the statistic that answered the question was chosen
after the numbers were on the page. It named the fix — **register the
$\chi^2$ in advance and run it on radicals none of these used.** Two
fresh ones here, five bases each, $2$-adic valuation spanning $1$ to
$11$, with base $20736$ riding along as the gate (Q1, exact).

| $\{2,19\}$, $v_2$ | $1$ | $2$ | $6$ | $10$ | $11$ |
|---|---|---|---|---|---|
| drift | $.066885$ | $.069391$ | $.106752$ | $.104068$ | $.098167$ |

| $\{2,3,7\}$, $v_2$ | $1$ | $2$ | $4$ | $8$ | $10$ |
|---|---|---|---|---|---|
| drift | $.249138$ | $.249852$ | $.251525$ | $.253167$ | $.248982$ |

**Q2 is refuted at $\chi^2=14.233$ against a cap of $9.49$, and it is
the outcome its own rule named as the most expensive this branch can
produce.** Five bases sharing the radical $\{2,19\}$ disagree beyond
their own errors. **So the drift is not a function of the radical**,
and [rem:whichfloor], [rem:radicallaw] and [rem:primecontrib] each
measured a quantity that moves with the base. **Their numbers stand
and their readings are withdrawn**, as registered — including "the
deficit carries the radical", "which primes, not how many", and the
contribution function $f(p)$ itself, every one of which was read off
one base per radical.

**Q3 holds and makes the failure sharper, not softer.** The
$\{2,3,7\}$ bases give $\chi^2=2.688$ with errors five times smaller,
and that test could have detected a base scatter of $0.002093$ — while
the effect visible in $\{2,19\}$ runs across $0.040$. So the same
design finds base dependence at one radical and excludes anything a
twentieth its size at another. The two groups disagree about whether
the base matters, and this run does not say why.

**What it moves along is not established.** Regressed on $v_2$ the
$\{2,19\}$ drifts give $+0.003599\pm0.001334$, $t=+2.70$; on $\log$
base, $t=+0.64$. But a base of fixed size cannot raise $v_2$ without
lowering the odd primes' valuations, and the two are correlated at
$-0.993884$ in $\{2,19\}$ and $-0.996116$ in $\{2,3,7\}$ — `COEFF NOT
SEPARABLE` on both. **The $v_2$ slope names a direction, not a cause.
Q2 establishes that the bases disagree, not what they disagree
along.**

**Q4 holds at $t=+33.89$** — the two radicals sit at $+0.088453$ and
$+0.250223$ — so radicals are not irrelevant either. Both things are
true at once: the radical moves the drift by a great deal, and within
one radical the base moves it by more than the errors allow.

**Where this leaves the branch.** Six ticks of radical structure
reduce to two controlled statements: the drift differs enormously
between radicals, and it is not determined by the radical alone. Every
functional reading built on top — the count-versus-primes question,
$f(p)$, the sub-additivity, item 5's distance carrying the radical —
rests on a variable now known to be incomplete. The honest next step
is a design that varies $v_2$ and $v_p$ independently, which requires
bases of different sizes and therefore a control for size as well.

#### Remark (the exponent is not it either, and the slope is not a number) {#rem:valuation}
<!-- evidence: audit_valuation.py -->

[rem:radicalblind] refuted the radical and barred reading its own
$v_2$ slope. Writing this run made the reason plainer: **in a doubling
family $v_2$ is not a property of the family at all** — for base
$2^ap^b$ the family is $N_j=2^{a+j}p^b$, so $v_2$ runs through ten
values inside every family while $b$ stays fixed. Read by $b$ that
remark's numbers group exactly, so $b$ was registered here as the
variable and tested on bases built to separate it from size: sizes
held inside a factor of two while $b$ runs one to five, giving
correlations of $+0.126$ and $-0.055$ between $b$ and $\log$ base.

| $p=5$, $v_5$ | $1$ | $2$ | $3$ | $4$ | $5$ |
|---|---|---|---|---|---|
| drift | $.146109$ | $.147112$ | $.144190$ | $.146704$ | $.138572$ |

**R2 and R3 are both refuted** — slopes $-0.001548\pm0.000922$
($t=-1.68$) and $-0.005679\pm0.004998$ ($t=-1.14$) — with resolvable
sizes of $0.002767$ and $0.014995$ printed beside them. R4 holds: the
size slope is unresolved at both primes. **The branch has now refuted
the radical and the exponent and has no candidate left**, which its
own rule required be said in those words.

**And then the diagnostic found why there may be nothing to explain**
(run after the verdicts, predicted by nothing). The fitted quantity is
the slope of a curved function — [rem:deficitdirect] put the curvature
at $-0.007380$ and [rem:deficitshape] found it steepening — so a
ten-point slope depends on where its window sits. Refitting each
family on its own first six and last six points:

| | between bases | **within one family** |
|---|---|---|
| $p=5$ | $0.008540$ | $\mathbf{0.041057}$ |
| $p=7$ | $0.035607$ | $\mathbf{0.031688}$ |

**One family's own drift moves by $0.03$ to $0.08$ depending on which
six of its ten points are used** — base $25600$ gives $+0.184332$ and
$+0.114863$, base $25088$ gives $+0.180851$ and $+0.104773$. That is
larger than the $\{2,3\}$ against $\{2,5\}$ gap of $0.074374$ this
branch has been chasing. **And none of it is in the printed standard
errors**, which measure scatter about a line, not the movement of the
line under its own window.

**Every family in this branch sits on a different window**, because a
different base gives a different set of $N$. So every cross-family
comparison here has used an uncertainty that omits its dominant term.
The pattern across runs is what that predicts: ordered by how far
their bases span in size, $\{2,19\}$ at $2.837\times$ gave $\chi^2=
14.233$ and failed; $p=7$ at $2.000\times$ spread $0.035607$; $p=5$ at
$1.600\times$ spread $0.008540$; $\{2,3,7\}$ at $1.500\times$ gave
$\chi^2=2.688$ and passed. **The disagreement tracks the size span,
not the arithmetic.**

**What this does and does not do.** It is a diagnostic, run after the
verdicts and predicted by nothing, so it does not reinstate what
[rem:radicalblind] withdrew — a post-hoc argument does not undo a
registered refutation. What it establishes is that the statistic this
branch has been comparing is unstable under a change no arithmetic
variable controls, and that the next test must not be a comparison of
slopes at all. Comparing $\log(\lvert\sum a\rvert/\ell^2)$ itself at
matched $N$ removes the window entirely, and that is the design this
question now needs.

#### Remark (the level, at matched N, with no window in it) {#rem:levelmatched}
<!-- evidence: audit_level_matched.py -->

[rem:valuation] found the statistic this branch had compared for six
ticks unstable under its own window and said the fix was to stop
comparing slopes. $L(N)=\log(\lvert\sum a\rvert/\ell^2)$ is measured
at one $N$ with no fitting at all. Twelve $N$ inside a band of factor
$1.20$, ten radicals, three of them appearing twice as the control.

| $N$ | radical | $L$ |
|---|---|---|
| $2097152$ | $\{2\}$ | $+1.341935$ |
| $2085136$ | $\{2,19\}$ | $+1.727299$ |
| $1874048$ | $\{2,11\}$ | $+1.812152$ |
| $1827904$ | $\{2,13\}$ | $+1.816345$ |
| $1882384,\ 2151296$ | $\{2,7\}$ | $+2.008906,\ +2.075538$ |
| $2000000$ | $\{2,5\}$ | $+2.166023$ |
| $1889568,\ 2125764$ | $\{2,3\}$ | $+2.995958,\ +3.032323$ |
| $1921920$ | $\{2,3,5,7,11,13\}$ | $+3.071761$ |
| $2000376$ | $\{2,3,7\}$ | $+3.141281$ |
| $1800000,\ 2025000$ | $\{2,3,5\}$ | $+3.187272,\ +3.232112$ |

**S3 holds by seventy-five times its own margin.** $L$ spreads
$1.890177$ across the radicals against the $0.024960$ the band's width
can contribute. **This is the first statement about the radical in
this branch that no window can be blamed for**, and it is about the
level rather than the drift.

**S2 is refuted, on one pair of three.** The $\{2,7\}$ pair differs by
$0.066632$ against a cap of $0.05$; $\{2,3\}$ gives $0.036365$ and
$\{2,3,5\}$ gives $0.044840$. **And the bound the rule compared them
against was itself one number where the branch had measured a range.**
The registration used a single trend of $0.14$ per log unit; the
drifts this branch has published run $0.066885$ to $0.284941$, a
factor of four, so a pair's own bound is its $\log N$ gap times its
own radical's drift — $0.007878$ to $0.033561$ for the first two,
$0.008931$ to $0.038049$ for $\{2,7\}$. All three $\lvert\Delta
L\rvert$ sit above the low end of their own bounds and the $\{2,7\}$
one sits above the high end. **The rule is not rewritten and S2 stays
refuted; the defect is that its bound was a single number, the ninth
time this branch has capped on the wrong quantity.**

**S4 is refuted, and informatively.** The primorial $N$ is not the
largest — $\{2,3,5\}$ is, at $+3.232112$ against $+3.071761$, with
$\{2,3,7\}$ between them. Adding $7$, $11$ and $13$ to $\{2,3,5\}$
*lowers* $L$. That is the sub-additivity [rem:radicallaw] measured in
the drift, now visible in a quantity with no fitting in it.

**And this run's closing sentence is wrong.** It was written to fire
when S2 failed and says the radical does not determine the level
either, so that everything this branch has said about radicals
"is describing something that is not there". S3 says otherwise by a
factor of seventy-five, and two of the three control pairs sit inside
the range their own radicals' drifts allow. What S2's failure
establishes is narrower: **at $\{2,7\}$, two $N$ of one radical differ
by more than that radical's own drift accounts for**, in a measurement
with no window to blame. One pair out of three is a finding; it is not
the collapse the sentence claims.

#### Remark (packed close: the level carries scatter, and how much) {#rem:leveldense}
<!-- evidence: audit_level_dense.py -->

[rem:levelmatched]'s control failed on one pair of three: two $N$ of
radical $\{2,7\}$, $0.133531$ apart in $\log N$, with $L$ differing by
$0.066632$ where a smooth $L$ at that radical's drift gives $0.018$.
Packing $N$ close inside one radical decides what that was — halve the
gap and a smooth $L$ halves the difference, while scatter stays put.
Fourteen $N$ of each of $\{2,7\}$ and $\{2,5\}$, all $2^ap^b$;
$128/125=1.024$ makes $\{2,5\}$ supply seven pairs at a $\log$ gap of
$0.0237$, where the two readings are $0.0034$ and $0.05$.

**T3 is refuted: seven of twelve close pairs exceed the cap.**

| radical | $\log$ gap | $\lvert\Delta L\rvert$ | smooth |
|---|---|---|---|
| $\{2,5\}$ | $0.023717$ | $0.064590$ | $0.003428$ |
| $\{2,5\}$ | $0.023717$ | $0.028794$ | $0.003428$ |
| $\{2,7\}$ | $0.025490$ | $0.062849$ | $0.003310$ |
| $\{2,5\}$ | $0.023717$ | $0.002779$ | $0.003428$ |

The largest is $0.064590$ — **nineteen times what a smooth $L$ gives
at that spacing** — and the median of the twelve is $0.017737$, five
times it. **The differences do not shrink with the gap**, which is
what scatter at fixed radical means.

**T2 and T4 hold, and the first of them is the scatter's size, not
its absence.** A line through each radical's fourteen points has
r.m.s. residual $0.029924$ and $0.028850$ against a cap of $0.03$ — a
margin of one part in three hundred — and those residuals *are* the
scatter T3 detects. The slopes, $+0.163956$ and $+0.127948$ against
published drifts of $+0.144537$ and $+0.129861$, agree within $0.05$:
**the level's own slope is the drift**, measured without any window.

**So the size is the finding, as the rule required.** Every statement
about the radical in this branch carries an unmodelled term of about
$0.03$ in $L$. Against [rem:levelmatched]'s radical separation of
$1.890177$ that is a sixtieth: **the radical effect is real and sixty
times the noise, and the noise is the precision floor for any radical
claim here.**

**This run's closing sentence is wrong,** as written to fire on T3's
failure: it says the radical does not determine the level. It
determines it to within a sixtieth of its own range. What T3
establishes is the size of what the radical does *not* determine, not
that it determines nothing.

**And the first execution had two wrong $N$**, disclosed in the result
file. Its $p=7$ list carried $1834496$ and $3668992$ as $2^{18}\cdot7$
and $2^{19}\cdot7$; they are $2^9\cdot3583$ and $2^{10}\cdot3583$, and
$2^{18}\cdot7$ is $1835008$. Those two produced that run's largest
$\lvert\Delta L\rvert$, both near one, by comparing across radicals —
the radical effect, not scatter. The list is corrected, every verdict
above is on the corrected field, and the script now asserts each $N$
has the radical it claims. T3 was refuted in that execution and is
refuted here on pairs the slip never touched.

#### Remark (the scatter is not a scale effect, and it is smaller here) {#rem:levelfine}
<!-- evidence: audit_level_fine.py -->

[rem:leveldense] measured the scatter at about $0.03$ and could not
say what it is, because $2^ap^b$ cannot be packed tighter than a
$\log$ gap of $0.0237$ at these $N$. Radical $\{2,3,5\}$ is a finer
ruler: enumerating every $N=2^a3^b5^c$ with all three exponents at
least one, $983040=2^{16}\cdot3\cdot5$ and $984150=2\cdot3^9\cdot5^2$
sit $0.001129$ apart, a twenty-first of that. The band
$[700000,\,1400000]$ holds forty such $N$, one coprimality class, and
the gaps are whatever the arithmetic gives.

**U2 is refuted and it is the well-powered verdict here.** Over
thirty-nine adjacent pairs the correlation of $\lvert\Delta L\rvert$
with the $\log$ gap is $-0.310186$, against a cap of $+0.5$ — not
merely below it but the wrong sign. **$\lvert\Delta L\rvert$ does not
follow the spacing, so the scatter is not a scale effect**, which is
the reading U2's rule attached to exactly this outcome.

**U3 is refuted and its reading is barred**, by the clause registered
before the run: the band supplied two pairs under a $\log$ gap of
$0.005$ where three were required. For the record they are gap
$0.004683$ with $\lvert\Delta L\rvert=0.006131$, and gap $0.001129$
with $0.017552$ against a smooth prediction of $0.000299$ — fifty-nine
times it, at the finest spacing this ruler reaches. **Two pairs are
two pairs and the conclusion is declined.**

**U4 holds, and its residual is the news.** The slope is $+0.251027$
against the published drift $+0.265357$, and the r.m.s. residual about
the line is $\mathbf{0.007363}$ — **four times smaller than
[rem:leveldense]'s $0.029924$ and $0.028850$.** So the scatter is not
one number across radicals any more than the level is.

**What is established and what is not.** The scatter does not shrink
with spacing (U2, thirty-nine pairs), so refining the $N$ does not
refine $L$; and its size at $\{2,3,5\}$ in this band is $0.0074$, not
the $0.03$ [rem:leveldense] found at $\{2,5\}$ and $\{2,7\}$. **The
comparison is not clean**: this band spans $0.6931$ in $\log N$ against
their $1.1632$ and $1.2528$, and a shorter span through a curved $L$
returns a smaller residual whatever the noise is. So the $0.0074$ is
an upper bound on this radical's scatter and not a measurement of a
difference between radicals.

**And this run's closing sentence is incomplete.** It fires on U3's
underpowering and says the band did not supply enough fine pairs to
decide — true of U3, and it ignores that U2 decided the same question
with thirty-nine pairs and the opposite kind of evidence. The
spacing question is answered; what is not answered is whether the
answer holds at the finest spacing, which is a smaller claim than the
sentence implies.

#### Remark (the residual is noise, on a closed list of candidates) {#rem:levelresidual}
<!-- evidence: audit_level_residual.py -->

[rem:levelfine] showed the scatter in $L(N)=\log(\lvert\sum
a\rvert/\ell^2)$ is not a scale effect and left what it *is* open.
Writing this run closed the list. Inside one radical every
construction-level quantity is fixed — the singular series depends on
$N$ only through its radical, and so does which $k$ the range excludes
— so what varies from $N$ to $N$ is only how $\log N$ is split among
the three primes. And that split is not three free numbers:
$$v_2\log2+v_3\log3+v_5\log5=\log N,$$
so after the $\log N$ trend **exactly two directions remain**.
Regressing the residual on the shares $s_2$ and $s_5$ tests the whole
list. **This is the first question in this branch whose candidate set
is complete rather than chosen.** Nothing is measured here; all forty
$L$ are read from markers, and the identity closes to $2.22\cdot
10^{-16}$.

The regressors have room: $s_2$ runs $0.0502$ to $0.8037$ and $s_5$
runs $0.1137$ to $0.8193$, so this is not a test starved of variation.

**V2, V3 and V4 are all refuted.**

| | coefficient | s.e. | $t$ |
|---|---|---|---|
| on $s_2$ | $-0.007474$ | $0.006655$ | $-1.12$ |
| on $s_5$ | $-0.011157$ | $0.006788$ | $-1.64$ |

$R^2=0.071081$, against a permutation null over $4000$ draws with
median $0.037793$ and ninety-fifth percentile $0.153099$. **The fit
sits between the null's median and its percentile** — indistinguishable
from shuffling the residuals against the same shares. The
underpowering clause did not fire: the percentile is below V3's cap,
so V3's threshold was the conservative one and the answer does not
turn on it.

**So the residual is not a function of $N$'s arithmetic in any form
this construction can express**, which is the reading V2's rule
attached to this outcome. The split is the only thing that varies
inside a radical once the trend is removed, and the residual does not
see it. **The $0.007363$ [rem:levelfine] measured is noise of the
object, and it is irreducible for a reason rather than by
observation.**

**What that costs and what it buys.** It costs any hope of sharpening
a radical statement in this branch below its own scatter — there is
nothing left to condition on. It buys the reason: the scatter is not
a variable this branch failed to find, it is the absence of one. The
list is closed *for this construction*; a quantity reaching $N$
through something both the $k$-range and the singular series ignore
would not be on it, and this repository knows of none.

#### Remark (the object is a Type II bilinear form with disjoint ranges) {#rem:bilinear}
<!-- evidence: audit_bilinear.py -->

[rem:levelresidual] closed the radical branch — the scatter above the
radical is noise of the object on a candidate list that is complete —
and item 5's demand was untouched by any of it. [rem:sieveweight] had
said what kind of object would admit an unconditional statement. Push
its form one step. For squarefree $j$, $\Lambda_K(j)=-\sum_{d\mid j,\,
d\le j/K}\mu(d)\log(j/d)$; writing $j=dm$, the condition $d\le j/K$ is
*exactly* $m\ge K$, so exchanging the order of summation gives
$$\sum a\Big|_{\text{covered}}
 = -\sum_{d\le N/K}\mu(d)\!\!\sum_{K\le m<N/d}\!\!\mu^2(m)\,
   \Lambda(N-dm)\log m ,$$
$d$ and $m$ squarefree and coprime to each other and to $N$. **That is
a Type II bilinear form**, the shape Vaughan's identity and its
descendants are written for.

**The exchange is exact** (X1): the form reproduces the part it covers
at relative $0$, $1.00\cdot10^{-15}$ and $1.63\cdot10^{-16}$ at
$N=20000$, $50000$, $200000$, with the gate reproducing
$\lvert\sum a\rvert=87895.3236$ (W1).

**And the ranges do not meet** (W3). Since $\theta=0.56>1/2$,
$N/K<K$:

| $N$ | $d$ up to | $m$ from | gap |
|---|---|---|---|
| $20000$ | $77$ | $257$ | $180$ |
| $50000$ | $113$ | $427$ | $314$ |
| $200000$ | $213$ | $933$ | $720$ |

so every contributing pair has $d<m$ with $N^{0.12}$ of clear air
between the variables. **The short side is the derived length** (W4):
the contributing $d$ number $27$, $40$, $72$ against
$(6/\pi^2)\prod_{p\mid N}\frac{p}{p+1}N^{1-\theta}$ giving $26.4$,
$39.5$, $72.6$ — within $2.4$ per cent, the same derived count
[rem:targetderived] measured for $\#k$.

**W2 is refuted and its target was the error, not the exchange.** It
compared the form against $\sum a$'s whole composite part, which also
contains $j$ like $12=2^2\cdot3$ where $\mu*\log$ is not $\Lambda$ and
the cofactor rewriting was never claimed. The discrepancies are
$1.768148$ at $N=20000$ and the *same* $1.768148$ at $200000$, and
exactly zero at $50000$ — a fixed set of terms, not a growing error —
and X2 confirms they are exactly [rem:sieveweight]'s "uncovered"
column to relative $10^{-15}$. The rule is not rewritten; W2 stays
refuted and the second block is registered as second.

**What this buys.** Item 5's remaining requirement is a bound on a
bilinear sum whose Möbius-weighted variable runs to $N^{0.44}$ and
whose $\Lambda$-weighted variable starts at $N^{0.56}$, with the
ranges separated. That is a question with a literature rather than a
question with no name — which is what [rem:deficitlog] left the branch
needing when it measured that no fit on this field can decide the
shape.

**And what it does not.** A bilinear form is a shape, not an estimate.
Nothing here supplies a bound, nothing here claims the known bounds
reach $\ell^2$ order, and no exponent is measured or forecast.
[rem:deficitlog] and [rem:shapepower] stand exactly where they stood.

#### Remark (the smallness is cancellation between moduli, and better than square-root) {#rem:bilinearcancel}
<!-- evidence: audit_bilinear_cancel.py -->

[rem:bilinear] named the object; it supplied no bound. What the
classical tools give can be read off the dimensions before anything is
computed, and the arithmetic is not encouraging:

| | exponent |
|---|---|
| trivial | $1$ |
| Bombieri–Vinogradov, error | $1$ |
| GRH applied to each modulus | $1.5-\theta=0.94$ |
| measured $\alpha$ ([rem:denominator]) | $+0.717916$ |
| the demand ([rem:denominator]) | $+0.587483$ |

Two things follow from the ranges alone. **BV covers $d\le
N^{1-\theta}$ unconditionally exactly when $\theta>1/2$** — the regime
the program needs, so the $d$-range is inside reach for the same
reason the reduction is interesting — and **the demand sits below what
GRH gives applied per modulus**, so no per-modulus input of any
strength reaches item 5. The measured $0.717916$ is far below $0.94$,
so cancellation over $d$ is already happening. This run measures how
much. (Nothing here reopens [rem:thetalaw]: $\alpha=1-\theta'/2$ is a
coincidence at one level and stays withdrawn.)

| $N$ | $\#d$ | $\sum_d\lvert I(d)\rvert/\lvert S\rvert$ | $\sqrt{\#d}$ | ratio |
|---|---|---|---|---|
| $20000$ | $27$ | $15.3979$ | $5.1962$ | $2.9633$ |
| $50000$ | $40$ | $21.9558$ | $6.3246$ | $3.4715$ |
| $200000$ | $72$ | $36.5215$ | $8.4853$ | $4.3041$ |

**Y2 holds: the smallness is cancellation over $d$**, by a factor of
$36.5$ at the largest $N$. **Y4 holds and locates the problem**: the
largest piece is $d=1$ — $I(1)=1409338.1646$, which is $104$ times
$\ell^2$ — so a single term of the $d$-sum is two orders above the
target and the whole of the smallness is other moduli cancelling it.
That is [rem:denominator]'s "the truncation does not shrink the main
term; it removes it and overshoots", now visible term by term.

**Y3 is refuted, and in the direction the rule called unexplained.**
The cancellation is not square-root but $2.96$, $3.47$, $4.30$ times
better, and **the excess grows with $N$**. A $d$-sum cancelling better
than random signs is not what any bound in evidence here would give,
and this run measures it without explaining it. Three points is three
points.

**What this establishes for item 5.** Its requirement needs
cancellation *between* moduli, which is not what level-of-distribution
input supplies: BV covers this $d$-range only because $\theta>1/2$,
and covering a range is not bounding a sum over it below the target.
**No bound is established here and none is claimed** — the four
exponents are arithmetic on the ranges with BV and GRH cited for what
they give, and the reading is that the demand lies below the
per-modulus route, not that the demand is impossible.

**The gate failed twice on the same defect before holding.** Y1 first
read the covered part from a six-decimal table and judged it at a
relative $10^{-12}$; the print bound there is $3.7\cdot10^{-11}$, so it
was refuted by rounding and not by disagreement. Emitting the marker
at ten digits did not fix it — that is still $5\cdot10^{-7}$ absolute,
the same wall one decimal out — and only full double precision did,
after which the three relatives are $2.67$, $3.75$ and
$4.89\cdot10^{-16}$. **`TOL BELOW PRINT` twice in one run**, and the
second instance of the defect [rem:whichfloor]'s M1 had. The tolerance
was never changed; only the source of the digits.

#### Remark (the demand is a main term removed to a relative accuracy) {#rem:maintermremoval}
<!-- evidence: audit_mainterm_removal.py -->

[rem:bilinearcancel] located the whole problem at one modulus: the
largest piece of the Type II sum is $d=1$, at $104$ times $\ell^2$. So
split there. With $A=I(1)$, $B=\sum_{d\ge2}\mu(d)I(d)$ and $S=-(A+B)$,
the gate reproduces the published covered part at relative
$8.15\cdot10^{-16}$ (Z1).

**$B$ fights $A$ at every $N$, and ever harder** (Z2):

| $N$ | $25000$ | $200000$ | $800000$ | $3200000$ |
|---|---|---|---|---|
| $\lvert B\rvert/\lvert A\rvert$ | $0.874261$ | $0.936626$ | $0.959641$ | $0.976934$ |

monotone to $0.977$ — the second variable cancels the first to within
$2.3$ per cent at the top and improving. **That is the structure of
item 5**: a term of main-term order is removed by the rest of the
$d$-sum, and what the demand asks is how exact the removal is.

**And the model-free number needs no exponent at all.**
$\lvert S\rvert/\ell^2$ runs

$$5.0350,\ 5.3385,\ 5.8927,\ 6.6274,\ 8.0213,\ 8.2110,\ 8.3955,\
9.5118$$

across the field. **The demand needs it bounded; it grows by $89$ per
cent over $2.1$ decades.**

**Z4 is refuted and the error is mine, two lines above its own cap.**
It asked whether $\lvert A\rvert$ has exponent $1$; it has
$+1.107864\pm0.002950$. $I(1)=\sum_m\mu^2(m)\Lambda(N-m)\log m$
carries a $\log m$ weight — written in this run's own formula — so its
main term is $N\log N$, not $N$. Fitted over the same eight $N$,
$N\log N$ gives $+1.080429$, and $\lvert A\rvert$ stands
$+0.027436$ from it. **The cap was written on the wrong quantity, not
the measurement on the wrong object.**

**Z3 is refuted for the same reason.** It compared the residue's
exponent against $\alpha-1=-0.282084$, which presumes
$e(\lvert A\rvert)=1$. Measured, the residue runs $-0.342848\pm
0.009184$ — six errors from the target, so this is a real
disagreement and not an unresolved one, and it is what
$e(\lvert S\rvert)-e(\lvert A\rvert)$ on this short field gives once
$\lvert A\rvert$ is allowed its logarithm.

**And this run's closing sentence is wrong.** It fires on Z4 and says
that because $A$ is not of main-term order, calling the requirement a
main-term removal is the wrong description. $A$ *is* a main term —
times $\log N$. **The description stands; the exponent in the cap did
not.**

**What the restatement is worth, stated with its own limits.** Item 5
is that a main term of order $N\log N$ is cancelled to leave something
of order $\ell^2$. The exponent form of that — a relative $N^{-1/2}$ —
rests on $\ell^2$ being $\sqrt N$ times logs, which is an asymptotic
reading and not measured here; this repository's fitted $e(\ell^2)$ is
$+0.583897$, and **the two statements of the demand differ by $0.084$
in exponent and this run does not choose between them.** The bounded
ratio above is the part that needs no choice.

#### Remark (the two rates, and the deficit reached without fitting the sum) {#rem:tworates}
<!-- evidence: audit_two_rates.py -->

[rem:maintermremoval] left a puzzle in its own numbers: the
cancellation improves — $\lvert B\rvert/\lvert A\rvert$ runs
$0.874261$ to $0.976934$ — and yet $\lvert S\rvert/\ell^2$ grows,
$5.0350$ to $9.5118$. Both hold because two rates are racing. The
arithmetic fixes what to measure: $\lvert S\rvert=\lvert A\rvert r$
exactly with $r=\lvert A+B\rvert/\lvert A\rvert$, so $e(\lvert
S\rvert)=e(\lvert A\rvert)+e(r)$ — verified to $2.22\cdot10^{-16}$
(AA2), which is what makes the labels trustworthy — and the demand
$\lvert S\rvert\lesssim\ell^2$ is
$$e(r)\ \le\ e(\ell^2)-e(\lvert A\rvert).$$

| quantity | exponent | s.e. |
|---|---|---|
| $\lvert A\rvert$ | $+1.107864$ | $0.002950$ |
| $r=\lvert A+B\rvert/\lvert A\rvert$ | $-0.342848$ | $0.009184$ |
| $\lvert S\rvert$ | $+0.765016$ | $0.011771$ |
| $\ell^2$ | $+0.629287$ | $0.006482$ |

**The required residue exponent is $-0.478577$ and the achieved one is
$-0.342848$.** The gap is $+0.135730$, fitted directly on
$\log(\lvert S\rvert/\ell^2)$ so its error is the error of the
difference and not the sum of two — $\pm0.010646$, the mistake
[rem:jbarrier]'s K1 made and not repeated here.

**That gap is the published deficit, reached by a route that never
fits $\lvert\sum a\rvert$ against anything** (AA3): $+0.135730$
against $+0.134019$, agreeing to $0.0017$. It is a cross-check of the
published number and not a new one, and the two fields differ —
$e(\ell^2)$ here is $+0.629287$ on $2.1$ decades against the published
$+0.583897$ on ten and a half.

**And in plain terms** (AA4). At $N=3200000$ the cancellation between
$A$ and the rest of the $d$-sum is exact to $0.023066$. To meet item
5's demand there it would have to be exact to $0.002425$ — **a factor
of $9.5118$ closer, half a per cent instead of two and a third.**

**The gate failed on the same defect for a third tick, and this time
it was anticipated in words and missed in fact.** AA1 read its markers
at full double precision "because ten digits was not enough last
time" — and that fixed the *reader* while both *emitters* still
carried seven significant figures, giving relatives of
$1.17\cdot10^{-7}$ and $5.33\cdot10^{-10}$ against a $10^{-12}$ cap.
`TOL BELOW PRINT` again. The `mainA` marker is this repository's own
with exactly one consumer, so it was widened and re-run with no
cascade; `audit_jbarrier.py` was not, because re-running it would make
[rem:jbarrierreach]'s and [rem:whichfloor]'s runs older than what they
read and G22 would demand those be re-run too — **a print-width fix is
not worth invalidating two measured runs.** So BB1, registered second,
judges the $\ell^2$ leg at that marker's own print bound of
$5\cdot10^{-7}$ and keeps $10^{-12}$ on the other; both pass. AA1
stays refuted as written.

#### Remark (the primes do not drop out) {#rem:residuemodel}
<!-- evidence: audit_residue_model.py -->

[rem:tworates] ended by asking where the factor of $9.5118$ lives
among the $d$. Writing this run produced a candidate that would have
removed the primes entirely. $I(d)$ counts prime powers $N-dm$, which
lie in the class $N\bmod d$, so its size should be governed by
$1/\varphi(d)$; and $m$ runs over $[K,N/d)$, a range emptying as
$d\to D$. With $w(d)=I(d)/A$ that gives
$$A+B=A\sum_d\mu(d)w(d),\qquad
w(d)\approx\frac{1}{\varphi(d)}\cdot\frac{1-d/D}{1-1/D},$$
and if it held, the whole of item 5's remaining difficulty would be
the rate at which the elementary Fejér-weighted Möbius sum
$\sum_{d\le D}\mu(d)(1-d/D)/\varphi(d)$ approaches zero — a classical
question with no $\Lambda$ in it.

**It does not hold, on all three tests.**

| $N$ | residue | model | ratio |
|---|---|---|---|
| $25000$ | $+0.125739$ | $+0.091290$ | $0.7260$ |
| $200000$ | $+0.063374$ | $+0.026120$ | $0.4122$ |
| $800000$ | $+0.040359$ | $+0.007945$ | $0.1969$ |
| $3200000$ | $+0.023066$ | $-0.001965$ | $-0.0852$ |

**CC3 is refuted and it is the one that decides.** The model's sum
falls away far faster than the residue and **changes sign at the top
$N$** — where, as this run's own rule says, a magnitude within a factor
of two of a quantity of the opposite sign is not agreement. CC4 puts
the rates at $-0.342848\pm0.009184$ measured against
$-0.876211\pm0.090961$ modelled; its verdict is barred by the power
clause, the model's own error exceeding the cap, but CC3 had already
decided.

**CC2 is refuted too, and its pattern is the informative part.**
At $N=3200000$ the top ten weights miss the model by

| $d$ | $3$ | $7$ | $11$ | $13$ | $17$ | $19$ | $23$ | $21$ | $33$ |
|---|---|---|---|---|---|---|---|---|---|
| per cent | $-12.20$ | $-9.68$ | $-9.10$ | $-9.59$ | $-9.52$ | $-10.30$ | $-11.28$ | $\mathbf{-34.12}$ | $\mathbf{-33.60}$ |

— **around ten per cent for prime $d$ and around a third for the two
composite ones**, $21=3\cdot7$ and $33=3\cdot11$. The shape is wrong
in a way that depends on how many primes divide $d$, and this run
measures that without proposing a correction: a factor fitted to two
composite values would be a law drawn from two points.

**So the primes do not drop out**, which is the reading CC3's rule
attached to this outcome. Item 5's difficulty stays exactly where
[rem:bilinear] put it — in a correlation of $\Lambda$ against a sieve
weight — and the reduction to an elementary Möbius sum is not
available. That is worth the run: a tempting reduction is closed
rather than left as an unexamined hope, and the closing is by a
pre-registered test that named its own failure mode before it fired.

#### Remark (the correction is not a product over the primes of d) {#rem:weightshape}
<!-- evidence: audit_weight_shape.py -->

[rem:residuemodel] left a pattern: the weights miss the model
$w_{\text{mod}}(d)=(1/\varphi(d))(1-d/D)/(1-1/D)$ by about ten per
cent for prime $d$ and about a third for composite. Writing
$c(d)=w(d)/w_{\text{mod}}(d)$, this run asks whether $c$ is a product
over the primes dividing $d$ — a test with no free parameters, since
the $\omega=1$ values determine every other. The test set was fixed
before the run: every squarefree $d\le50$ coprime to $N$.

**DD2 is refuted, and the failure is systematic rather than
scattered:**

| $d$ | $c(d)$ | $\prod_{p\mid d}c(p)$ | per cent |
|---|---|---|---|
| $21=3\cdot7$ | $+0.658776$ | $+0.793053$ | $-16.93$ |
| $33=3\cdot11$ | $+0.664015$ | $+0.798162$ | $-16.81$ |
| $39=3\cdot13$ | $+0.656835$ | $+0.793825$ | $-17.26$ |

Three two-prime $d$, all short of the product by $16.8$ to $17.3$ per
cent — a spread of $0.45$ points. **So the correction is not a product
over the primes, and what it is instead is not scattered noise.** This
run measures that and proposes nothing: a factor fitted to three
points would be a law drawn from three points, and this branch has
made that mistake before.

**DD3 is refuted, narrowly, and $c(p)$ is not flat either.** The
spread is $0.053187$ against a cap of $0.05$, and the values are not
random about their mean of $0.886295$ — they rise from $0.8780$ at
$p=3$ to $0.9090$ at $p=11$ and fall monotonically to $0.8559$ at
$p=47$.

**DD4 is refuted and could not have tested what it was for.** Its
ratios run $1.3964$, $1.4461$, $1.4707$, $1.4356$, $1.3780$,
$1.4267$, $1.4767$ and — **at the very $N$ where $c(p)$ was fitted** —
$1.4754$. A test of transport that fails identically at its own
fitting point has not measured transport. The cause is in the
construction: $c(p)$ exists only for $p\le50$, the fixed test set, so
every $d$ with a larger prime factor was corrected by $1$ where the
measured corrections sit near $0.886$. **The verdict stands as
written and the defect is the test's**, one more in this branch's list
of rules no outcome could pass.

**And what DD4 does show, read for what it is, is the opposite of its
verdict.** The eight ratios span $0.0987$ about a mean of $1.4382$ —
nearly constant across two decades. **The correction transports; the
level is uniformly wrong** for the reason above.

**So this run's closing sentence overstates.** It fires on DD2 and
says the weights carry structure no multiplicative model expresses.
What is refuted is *this* multiplicative form, a product over the
primes of $d$; a systematic $-17$ per cent at every $\omega=2$ is
exactly what a model with an $\omega$-dependent factor would produce,
and such a model is not tested here. **The shape of $w$ remains
undescribed, and nothing here revives the reduction
[rem:residuemodel] closed.**

#### Remark (kappa predicts three primes from two, and the transport is real) {#rem:omegafactor}
<!-- evidence: audit_omega_factor.py -->

[rem:weightshape] refuted the product form for $c(d)=w(d)/
w_{\text{mod}}(d)$ and refuted it *systematically* — its three
two-prime $d$ came in $16.8$ to $17.3$ per cent below
$\prod_{p\mid d}c(p)$ — which is what an $\omega$-dependent factor
produces. It also found its own transport test void: DD4 failed by
$1.4754$ at the very $N$ where $c(p)$ was fitted, because $c(p)$
existed only for $p\le50$. Both are fixed here, and the second fix is
**written into the rule** rather than found afterwards.

**EE3 holds, blind.** With $\kappa$ taken from $\omega=2$ alone,
$c(d)=\prod_{p\mid d}c(p)\,\kappa^{2}$ predicts the three-prime $d$:

| $d$ | $231$ | $273$ | $357$ | $399$ | $429$ | $483$ |
|---|---|---|---|---|---|---|
| per cent | $-4.50$ | $-0.34$ | $+4.20$ | $+2.93$ | $+1.94$ | $-2.65$ |

six of them, all inside $4.5$ per cent, with no free parameter. **A
correction fitted on two primes carries to three.**

**EE4 is refuted on transport, and this time that means something.**
The fitted $N$ comes in at $1.0448$ against its cap of $1.05$ — so the
hole DD4 fell into is closed and the seven others are readable. They
run $1.2747$, $1.2853$, $1.2743$, $1.2096$, $1.1226$, $1.1175$,
$1.1085$ — the four smallest $N$ outside the $1.2$ cap. **The factors
fall monotonically toward $1$ as $N$ grows**, so the model is not
wrong so much as incomplete at small $N$; what it is missing is not
measured here.

**EE2 is refuted, and its statistic was the wrong one.** The cap was
written on $\max-\min$ over however many $d$ the range supplies — $73$
here against the $3$ [rem:weightshape] could test — and a range grows
with the count. The same $73$ ratios have s.d. $0.025290$,
interquartile $0.022543$ and median $0.829575$ against
$\kappa=0.832904$. **The rule is not rewritten and EE2 stands
refuted**; the defect in its statistic is recorded as one, and it is
the reason EE3 could still hold on a $\kappa$ its own test called
scattered.

**And this run's closing sentence is wrong.** It fires when EE2 fails
and says the blind set is too small or $\kappa$ is not a number, so
"the omega model is not tested here and nothing is claimed for it."
The blind set was six against a registered minimum of three, and EE3
held. **The $\omega$ model was tested and it passed the test it was
built for**; what failed is a spread statistic and a transport cap at
small $N$.

**What this describes and does not.** $w$ is $1/\varphi(d)$ times a
taper, times a product over the primes of $d$, times
$\kappa^{\omega-1}$ with $\kappa\approx0.833$ — a description that
survives a blind extension in $\omega$ and improves with $N$. It is a
description of the weights and not a bound on anything, and
[rem:residuemodel]'s closure of the reduction to an elementary Möbius
sum stands untouched.

#### Remark (the corrections drift, most for the smallest primes) {#rem:localcorrections}
<!-- evidence: audit_local_corrections.py -->

[rem:omegafactor] left the weight model right in shape and wrong in
transport, and asked what shape the convergence has. The shape of a
symptom is the wrong thing to fit when the cause is one step away: the
constants were measured at one $N$ and used at seven others, so the
question is whether $c(p)$ and $\kappa$ are themselves functions of
$N$. Fitting them locally at every $N$ answers it.

**They are, and decisively for the smallest primes.**

| | values across the eight $N$ | slope | $t$ |
|---|---|---|---|
| $c(3)$ | $0.9325\to0.8780$ | $-0.011551\pm0.000876$ | $\mathbf{-13.19}$ |
| $c(7)$ | | $-0.005152\pm0.000543$ | $-9.49$ |
| $c(11)$ | | $-0.001589\pm0.001154$ | $-1.38$ |
| $\kappa$ | $0.7775\to0.8329$ | $+0.010340\pm0.003562$ | $+2.90$ |

**FF3 holds** — $\kappa$'s range is $0.055367$ against a largest
standard error of $0.034522$, so the unresolved clause does not fire
and the movement is real. **FF4 is refuted**, at $t=+2.90$ against a
cap of $3$: $\kappa$ moves without a trend resolved at this bar, and
no direction may be read into it. The $c(p)$ carry no such doubt, and
their drifts fall away with $p$ — the correction moves fastest where
the prime is smallest.

**FF2 is refuted, and by how much matters.** The local factors are
$1.0230$, $1.0427$, $1.0728$, $0.9396$, $1.0435$, $1.0649$, $1.0672$,
$1.0448$, four of them outside the $1.05$ cap. But
[rem:omegafactor]'s transported factors ran to $1.2853$: **local
constants absorb most of the failure, from $28$ per cent down to $7$,
and not all of it.** The registered reading — that the model is
incomplete in a way locally fitted coefficients cannot absorb — stands
as written, and the size of what is left is $7$ per cent rather than
the $28$ it started from. The residual is not a systematic offset
either: $N=200000$ comes in *below* one at $0.9396$ while the rest sit
above.

**A crash was fixed before any verdict existed, and it is a fact about
the model.** The first execution divided by zero at $N=400000$, where
$D=291=3\cdot97$ is squarefree and coprime to $N$ so $d=D$
contributes — and the model's taper $(1-d/D)$ is exactly zero there.
At $N=200000$ it did not arise because $D=215=5\cdot43$ shares a
factor with $N$. **The model predicts $w(D)=0$ exactly and the
measurement gives $+0.00001449$**, so $c(D)$ is not a ratio that
exists; $d=D$ is excluded from the fits and its weight is printed
rather than papered over.

**What this leaves.** The description is $w=(1/\varphi(d))$ times a
taper, times $\prod_{p\mid d}c(p,N)$, times $\kappa(N)^{\omega-1}$ —
with the constants now known to move, most at $p=3$. Reading them
locally accounts for three quarters of the transport failure and not
the rest, so the question has moved to those constants **and has not
been answered by moving it**: that is progress only if they are
simpler than what they replace, and nothing here claims they are.
[rem:residuemodel]'s closure stands.

#### Remark (the drift belongs to two primes, not to the family) {#rem:driftbyprime}
<!-- evidence: audit_drift_by_prime.py -->

[rem:localcorrections] saw the weight corrections drift with $N$ and
drift most at the smallest prime, and said in its own words that three
primes is not a shape. Twenty primes contribute above the noise
threshold at all eight $N$ here, and they settle it.

**GG2 is refuted: two of the twenty resolve.**

| $p$ | $3$ | $7$ | $11$ | $17$ | $37$ | $73$ |
|---|---|---|---|---|---|---|
| drift | $-0.011551$ | $-0.005152$ | $-0.001589$ | $+0.004238$ | $-0.019126$ | $+0.018483$ |
| $t$ | $-13.19$ | $-9.49$ | $-1.38$ | $+1.27$ | $-1.80$ | $+1.18$ |

Beyond $p=7$ nothing reaches $\lvert t\rvert=2$, and the signs are
mixed — $+0.004238$ at $17$ against $-0.019126$ at $37$ and
$+0.018483$ at $73$, with errors three to twenty-six times the
$p=3$ signal. **The drift is a fact about $c(3)$ and $c(7)$ and not
about the corrections as a family**, which is the reading GG2's rule
attached to this outcome.

**GG3 and GG4 both say "hold" and neither is read**, by the rule fixed
before the run: with fewer than four resolved primes they have nothing
to run on. GG3's ordering is automatic on two points. GG4's fit
returns a slope of $-0.952938$ — temptingly close to $-1$ — with a
standard error that is **not finite**: two points give zero degrees of
freedom, the line passes through them exactly, and no $t$ is emitted,
because a $t$ against an infinite error is not a $t$. The gate caught
that as G39 when the run first tried to print one.

**And the separability the design named in advance did fire.** Over
the resolved set $\log p$ and $\log(p-1)$ correlate at $1.000000$;
`COEFF NOT SEPARABLE` is emitted, so even had four primes resolved,
$1/p$, $1/(p-1)$ and $1/(p+1)$ would not have been told apart —
only an exponent, and not a form.

**What this costs and what it leaves.** It costs
[rem:localcorrections]'s picture of a drift dying away with $p$: that
was drawn from three primes of which one was already unresolved, and
at twenty it does not survive. What stands is narrower and firmer —
$c(3)$ and $c(7)$ move with $N$ at $t=-13.19$ and $-9.49$, and
nothing else measurably does. Whether that is two primes being
special or eighteen being too noisy to tell, this field does not say:
the errors at $p\ge11$ are larger than the effect at $p=3$, so the
absence is an absence of power as much as an absence of drift, and
the remark claims no more.

#### Remark (power was bought, and the slope is not a slope) {#rem:driftpower}
<!-- evidence: audit_drift_power.py -->

[rem:driftbyprime] found two of twenty primes with a resolved drift
and refused to say whether $3$ and $7$ are special or the rest too
noisy, because the errors at $p\ge11$ exceeded the effect at $p=3$.
Power is buyable here: $\lvert w(p)\rvert$ is small partly because $p$
sits high in a $d$-range ending at $D=\lfloor(N-1)/K\rfloor$, so
raising $N$ raises $D$ and moves the same $p$ lower in it. Three more
$N$ take $D$ from $728$ to $1815$.

**HH3 holds and that is what makes the rest readable.** The $p=11$
standard error falls from $0.001154$ to $0.000587$, **$49.10$ per
cent** — power was genuinely bought, so a failure to resolve is about
the primes and not about the extension. That guard was written into
the rule before the run, in the shape that made
[rem:omegafactor]'s transport failure readable.

**HH2 is refuted, and the two watched primes fail differently.**

| $p$ | eight-$N$ $t$ | eleven-$N$ $t$ | drift |
|---|---|---|---|
| $11$ | $-1.38$ | $\mathbf{-2.37}$ | $-0.001589\to-0.001391$ |
| $13$ | $-0.50$ | $-0.49$ | $-0.001948\to-0.000977$ |

At $p=13$ the error halved and $t$ did not move: that is a prime with
no drift to find. At $p=11$ the error halved and $t$ rose from $1.38$
to $2.37$ — **still short of $3$ and heading for it**, so the run's
closing sentence, that $3$ and $7$ are special rather than merely
loudest, is right about $13$ and premature about $11$.

**HH4 holds at ninety-six per cent of its cap, and the diagnostic says
why that matters.** It capped how far the drifts moved, not by what
factor. Relatively:

$$p=3:\ 0.83298,\qquad p=7:\ 0.83429,\qquad p=11:\ 0.87551 .$$

**The two resolved primes shrink by the same factor, differing by
$0.00130$.** A slope that shrinks when the field lengthens is not a
slope, and two slopes shrinking together say the linear model is
failing the same way at both. This run measures that and proposes
nothing to replace it.

**And it is the same disease twice.** [rem:valuation] found that a
ten-point slope of the deficit moves $0.03$ to $0.08$ under its own
window, because the underlying function is curved; the corrections'
drifts now do the identical thing, by seventeen per cent over one
extra decade. **Every "drift" this branch has quoted for $c(p,N)$ —
including [rem:localcorrections]'s $t=-13.19$ and this run's own
numbers — is a window-dependent quantity**, and the resolution of
$c(3)$ and $c(7)$ says they move, not how fast.

#### Remark (dm/se, as v1 said to fit and nobody had) {#rem:maskdmse}
<!-- evidence: audit_mask_dmse.py -->

v1 handed one instruction forward with its open item on the mask's
decay exponent and it has sat unexecuted ever since: **fit $dm/se$,
not $dm$.** The reason is visible in the table it left. At depth 5
the cell population runs $2,3,4,5,8,11,16,22$ across the bands while
$se$ falls from $0.7136$ to $0.4727$; over the same bands $|dm|$ falls
from $7.0004$ to $4.9979$. **The error shrinks faster than the
amplitude**, and $|z|=|dm|/se$ — the quantity that says whether the
mask is there at all — *grows*, from $9.81$ to $10.57$.

Nothing is measured here; the fit v1 used is identified rather than
assumed. Weighted least squares of $\log|dm|$ on $\log N$ with
weights $(dm/se)^2$ and the covariance scaled by the residual
variance returns all six published exponents exactly and four of the
six published standard errors to the four decimals printed, the other
two in the last digit (E1).

**There is something to subtract** (E2). The exponent of $se$ itself
is positive and resolved at every depth, $t$ running $11.18$ to
$69.92$. So part of what was fitted as the mask decaying is the
error decaying, and the split is exact: the $|dm/se|$ exponent is the
$|dm|$ exponent minus the $se$ one (E3, smaller at every depth).

**And at the deepest cell the decay does not survive it** (E4). Depth
5 goes from $0.1434\pm0.0155$ — quoted at $9.2$ standard errors — to
$0.0317\pm0.0231$, $t=1.37$: **unresolved**. What that depth was
measuring was mostly its own error shrinking as the cell filled up.
Depth 1, already unmeasurable in v1's table, goes to $t=0.10$.

**But the reading v1 built on it survives** (E5, refuted). Excluding
depth 1, the $|dm|$ exponents deepest-first run $+0.1434$, $+0.2152$,
$+0.2713$, $+0.3686$, $+0.6289$ and the $|dm/se|$ exponents run
$+0.0317$, $+0.1655$, $+0.2309$, $+0.3283$, $+0.5869$ — **monotone
both times**. The order is a property of the mask and not of how the
cell populations grow with depth, and *"the mask decays faster where
fewer small primes divide $N$"* stands as written.

So the instruction was worth executing and it does not overturn the
paper. It removes one number: depth 5's exponent should not be quoted
as a measured decay. It leaves the rest smaller and intact, and it
leaves the open item where it was — [rem:shapepower]'s point applies
here too, since these are still exponents of an assumed form over a
factor $160$ in $N$, and nothing in this re-analysis touches that.

One thing about the run itself is worth recording. E5's comparison
was implemented in the wrong direction — v1 reports the exponent
rising as the cell gets shallower, which listed deepest-first is an
increase, and the code tested for a decrease. The verdict was
unaffected, but the line it printed said the $|dm|$ exponents were
*not* monotone as v1 reports, which contradicts a published table.
**That contradiction is what caught it**, and it is an argument for
printing the control's intermediate rather than only its verdict.


#### Remark (the floor test and the fit answer different questions) {#rem:maskfloornull}
<!-- evidence: audit_mask_floornull.py -->

OPEN.md item 3 has said the mask has no decay exponent because at the
three shallowest depths the amplitude does not clear the exact floor
at any scale measured, **so nothing can be fitted there.** The premise
is exact and F2 confirms it band by band: no band reaches $|z|=2$ at
depth 0 or 1, only two of fifteen do at depth 2, and every band does
at depths 3, 4 and 5.

**The conclusion does not follow, and the null says so.** Under dm
redrawn as $N(0,se^2)$ at each band and fitted the same way, the
exponent comes out at $0.0437$ in the median — near the exponent of
$se$ itself, because the magnitude of noise tracks its own scale —
with a 95 per cent range of $[-0.2076,\,0.2945]$. Depth 0's observed
$0.6289$ sits outside it at $p=0.0003$ (F3) and depth 2's $0.3686$ at
$p=0.0043$ (F4). Depth 1's $0.0437$ sits inside at $p=0.4755$ (F5),
so v1's "not measurable" stands exactly where v1 said it.

The reason the two tests disagree is not subtle once stated. The
floor test asks each band separately whether the amplitude is there;
the fit asks fifteen bands together whether it is *falling*. At depth
0 the amplitude runs $0.1005$, $0.0847$, $0.0632$, $0.0512$, $0.0411$
across the first five bands — a smooth fall whose every term is
inside its own error. **A quantity can decay systematically without
any single measurement of it being significant**, and OPEN.md's
sentence read one test's silence as the other's.

**But the same null corrects the exponents in the other direction,
and by more.** The fit's own standard error at depth 0 is $0.0121$,
which puts the exponent at $52.0$ standard errors. The null's spread
implies $0.1281$ — **eleven times larger** — and against that the
exponent is $4.9$. At depth 2 the quoted $0.0052$ becomes $0.1298$,
twenty-five times larger, and $70.9$ standard errors become $2.8$.
The fit's errors assume its weights are the true inverse variances of
$\log|dm|$, and that breaks exactly where $|dm|$ is comparable to
$se$ — which is the definition of these depths.

So item 3's sentence is wrong and its caution was right. Something
*can* be fitted at depths 0 and 2, and what is fitted is much less
certain than the fit reports: two of the six exponents survive as
measurements at about $5$ and $3$ standard errors rather than $52$
and $71$, one (depth 1) is noise as v1 said, and [rem:maskdmse]
already removed depth 5's by a different route. **Three of six
remain, at errors an order of magnitude wider than published.**

Nothing here touches the form. These are still exponents of an
assumed $N^{-a}$ over a factor $160$ in $N$, which
[rem:ladderdegree] and [rem:deficitregion] have now measured to be
the binding limitation on both other axes, and it is the binding one
here too.


#### Remark (a design that could not ask its question, and what that shows) {#rem:maskformreach}
<!-- evidence: audit_mask_formreach.py -->

v1 stated a limitation and OPEN.md item 3 carried it unquantified:
over a factor $160$ in $N$ the data do not separate $N^{-a}$ from
$(\log N)^{-b}$, so the mask's exponents are exponents of an assumed
form. That says the forms do not separate. It does not say what
would, and the answer decides whether the limitation is this
program's budget or the question's nature.

This run set out to measure it and **failed, in a way worth
recording.** The exponents reproduce (H1) and v1's statement holds
where it was made: drawing from the power law at the observed range,
the power law wins on weighted r.m.s. in $0.2440$ of draws at depth
0, so the forms are not separated (H2). Extending the range with the
band density and the error law fixed, the fraction does not rise —
$0.1970$, $0.2110$, $0.2155$, $0.2360$, $0.1995$, $0.2215$, $0.2185$
across extensions to a factor $10^{266}$ in $N$. **H3 and H4 are
refuted, and neither refutation means anything about the two forms.**

The diagnostic says why, and it is one number. The effective sample
size, $(\sum w)^2/\sum w^2$, reads $4.99$ at the observed range and
$5.01$ at every extension after it, while the band count runs $15$,
$29$, $57$, $113$, $225$, $449$, $897$, $1793$. **Every band added
past the observed range carries no weight.** The weights are
$(dm/se)^2$, and [rem:maskdmse] measured the amplitude's exponent at
$0.6289$ against the error's at $0.0420$; the weight therefore falls
like $N$ to twice their difference, and the design is asking for
information from a region where the mask has already sunk below its
own error.

So the registered question was malformed, and its malformation is
the finding. **The informative window is bounded.** Extending $N$
does not buy range for this measurement — past a point it buys bands
that measure nothing, and the point is set by how fast the signal
falls relative to the error, not by any budget. Whatever separates
$N^{-a}$ from $(\log N)^{-b}$ for this mask, it is not more $N$.

That is the third axis to reach the same place by a different route.
[rem:ladderdegree] found the level axis unable to challenge its own
family for want of points; [rem:deficitregion] found the sign axis
resolving a new degree at almost every order while the residual
conceded nothing; and here the wall's mask cannot be given more
information at all. **On none of the three does the form follow from
the data, and on none of the three does more computation change
that.**

One caution on this remark itself. What is shown is that *this*
design — v1's weights, v1's error law, the band density held fixed —
cannot separate the forms at any range. A design that reweighted, or
that measured the mask where it has not yet sunk, is untouched by
this and unexplored. The claim is about the reach of a measurement
already made, not about the mask.


#### Remark (where the mask has not sunk, the form is measurable) {#rem:maskdeepform}
<!-- evidence: audit_mask_deepform.py -->

[rem:maskformreach] could not measure what range separates $N^{-a}$
from $(\log N)^{-b}$, because at the depths it used the weights
collapsed — the effective sample froze at $5.01$ while the band count
went to $1793$. It named the untried design: measure where the mask
has not yet sunk below its own error. That design is in the same
table, and it works.

**The mechanism is visible in the exponents** (I1). The weight is
$(dm/se)^2$ and falls like $N$ to twice the gap between the amplitude
and error exponents. At depth 3 that is $N^{-0.4657}$ and at depth 4
$N^{-0.3375}$ — still collapsing. **At depth 5 it is $N^{-0.0659}$**,
almost flat, because [rem:maskdmse] found the amplitude and the error
decaying at nearly the same rate there.

**So the sample grows** (I2). At depth 5 the effective size runs
$13.89$ to $88.78$ as the range extends, a factor of $6.4$, against
the $4.99$ to $5.01$ [rem:maskformreach] recorded. Depths 3 and 4 sit
between, reaching $12.58$ and $17.34$ and then stopping.

**v1's statement holds where it was made** (I3). At the observed
range the power law wins in $0.5640$, $0.5705$ and $0.5515$ of draws
at the three depths — barely better than a coin, exactly as "the data
do not separate them" says.

**And extending the range does separate them** (I4). At depth 5 the
fraction runs $0.5515$, $0.7235$, $0.9900$, then $1.0000$ at every
larger extension. **The first extension to clear $95$ per cent
corresponds to a factor of $10^{7.72}$ in $N$.** Depths 3 and 4
plateau instead, at $0.34$ and $0.78$, for the reason their weight
exponents give.

So OPEN.md item 3's limitation is a reach and not a property of the
question, and the reach now has a number. It is not a small one (I5):
starting from the bands' floor at $N\sim1.189\cdot10^5$, a factor of
$10^{7.72}$ puts the top at $6.24\cdot10^{12}$. The kind-byte sieve
this session built carried $8\cdot10^9$ in $14.90$ GB; $6.24\cdot10^{12}$ is $780$ times that and would need about $12.5$ TB by the
same packing. **The form is measurable and this machine cannot
measure it.**

That is a different sentence from the three the other axes produced,
and worth keeping separate from them. [rem:ladderdegree],
[rem:deficitregion] and [rem:maskformreach] each found a place where
more computation buys nothing. **This one found a place where it
would buy the answer, and priced it.** A limitation with a price is a
different object from a limitation without one, and item 3 should
carry the price rather than the word "assumed".

Two cautions. The separation measured is between two named forms and
says nothing about a third. And depth 5 is the deepest cell — the one
with fewest bands and the smallest population — so the design that
succeeds here is the one with the least data behind it, which is not
a contradiction but is worth stating: what makes it work is the
*ratio* of the two exponents, not the amount of signal.


#### Remark (the price was against one rival, and there is no price of the form) {#rem:maskrivals}
<!-- evidence: audit_mask_rivals.py -->

[rem:maskdeepform] priced the mask's form question at a factor of
$10^{7.72}$ in $N$ and closed with a caution: the separation measured
is between two named forms and says nothing about a third. That
caution is the whole story, and this run measures it.

Three two-parameter rivals are put against the power law at depth 5,
on [rem:maskdeepform]'s design unchanged — the same bands, error law
and extensions — so the comparison of weighted r.m.s. is fair and the
$(\log N)^{-b}$ reach reproduces its $10^{7.72}$ exactly (J1). At the
observed range none of the three is separated: the power law wins in
$0.5560$, $0.5265$ and $0.5005$ of draws, which is a coin (J2).

**And the reaches are not one number.**

| rival | separates at |
|---|---|
| $A - b\log x$ | $10^{7.72}$ |
| $A - c\sqrt{x}$ | $10^{15.44}$ |
| $A - d\,x/\log x$ | $10^{30.88}$ |

with $x=\log N$. They span $23.16$ in $\log_{10}$ of the factor (J3),
and the hardest costs more than [rem:maskdeepform] quoted (J4).

**So $10^{7.72}$ is the price against $(\log N)^{-b}$ and not the
price of establishing the form**, and OPEN.md item 3 carried that
mistake for one commit on my writing of it. The correction is not
just a larger number. The ordering says what the number is doing:
the rivals separate in the order of how nearly linear in $x$ they
are, and $x/\log x$ — the shape a sieve bound takes — is nearly
linear over any finite range. **A two-parameter rival can be made
asymptotically as close to $x$ as one likes, so there is no supremum
to take.** "The price of the form" is not a large number; it is not a
number.

What survives is narrower and still worth having. Against a *named*
rival the question is decidable and the price is computable, and the
three prices are now on the record. That is the shape of every
statement this program can make about form: not *the mask decays like
$N^{-a}$*, but *the mask's decay is distinguishable from $(\log
N)^{-b}$ at $10^{7.72}$, from a stretched exponential at $10^{15.44}$,
and from $x/\log x$ at $10^{30.88}$* — each a fact, none of them the
fact that was wanted.

This closes the form branch of item 3 the way the other three axes
closed: not by answering the question but by showing what answering
it would take, and here the showing is exact. [rem:ladderdegree]
found a family the data cannot challenge, [rem:deficitregion] a
family the data refutes at every order, [rem:maskformreach] a design
that cannot be given more information. **This one found that even
where the design works and the price is finite, the price depends on
what you are pricing against, without bound.**


#### Remark (mu's share of the spectrum is structured, and by something Lambda's weights do not contain) {#rem:spectrummushare}
<!-- evidence: audit_spectrum_mushare.py -->

OPEN.md's wall item 6 stood untouched from the start. The spectral
measure is atomic at the rationals $j/q$ with weights
$\mu^2(q)/\varphi^2(q)$, and those weights are $\Lambda$'s structure;
$\mu$'s own contribution had been measured only through the
principal-arc deficit, at $q=3$ and $q=5$. Whether anything of $\mu$
lay beyond those two arcs was never asked.

**The table those two numbers come from already answered it.** v1's
`lab_atoms_perq.txt` records, for each of $31$ moduli, the share of
the spectrum the real field puts at $j/q$ and the share the coin puts
there. The coin replaces $\mu$ by random signs and keeps everything
else, so **the ratio of the two shares is $\mu$'s contribution by
construction**, not by inference. Nothing was computed here that was
not already on disk; the two quoted figures reproduce to the printed
decimals as the gate (K1).

And the ratio is not flat. It sorts by $\omega(q)$, the number of
prime factors:

| $\omega(q)$ | moduli | real/coin |
|---|---|---|
| $1$ | $3,5,7,11,13$ | $0.080$–$1.149$ |
| $2$ | $15,\dots,143$ | $1.404$–$21.627$ |
| $3$ | $105,\dots,1001$ | $2.692$–$6.951$ |
| $4$–$5$ | $1155,\dots,15015$ | $3.810$–$6.800$ |

**$\mu$ suppresses the prime moduli and enhances the composite ones.**
Regressing $\log(\text{real}/\text{coin})$ on $\omega$ gives slope
$+0.7473 \pm 0.1420$, $t = 5.26$ (K2). With $\log q$ added, $\omega$
keeps its sign and its resolution at $t = 2.49$ while $\log q$ does
not resolve at all ($t = -1.01$), so $\omega$ is not standing in for
the size of the modulus (K3) — and it survives a correlation of
$0.957$ between the two regressors to do it.

One check was not pre-registered and is reported because the figure
demands it: $q=15$ sits at $21.627$, three times the next largest, and
a point that size can manufacture a slope. It does not. **Dropping it
raises $t$ from $5.26$ to $6.97$** — the outlier was flattening the
fit, not making it — and over all $31$ leave-one-out fits the smallest
$|t|$ is $4.75$.

$\omega(q)$ is not in $\mu^2(q)/\varphi^2(q)$. So there is a component
of the spectrum that is $\mu$'s and not $\Lambda$'s, it is not
confined to the principal arcs, and it has a sign: **the deficit
OPEN.md quotes at $q=3$ and $q=5$ is the $\omega=1$ case of a
monotone law, not a fact about $3$ and $5$** (K4). Read back through
that law, $8.40$ and $15.16$ are the two most extreme suppressions in
the table because $3$ and $5$ are the two smallest primes on it — the
principal arcs are where $\mu$'s effect is largest, which is why they
were the first place it was seen.

Two limits. The direction is measured and the exponent is not: a slope
of $0.75$ per prime factor is what these $31$ moduli give, and no
claim is made that the true dependence is exponential in $\omega$
rather than, say, in $\log\varphi(q)/\log q$ — this program has been
shown four times over (rem:ladderdegree, rem:deficitregion,
rem:maskformreach, rem:maskrivals) that it cannot settle functional
form, and nothing here is an exception. And the field is one $n$;
whether the slope drifts with $n$ is a measurement v1's design would
have to be re-run to make.


#### Remark (the flatness cannot rise forever) {#rem:flatnessshape}
<!-- evidence: audit_flatness_shape.py -->

Remark [rem:leanidentity] left one exponent loose: $e(\ell^1/\ell^2)$
measures $0.287798\pm0.002472$ against its own ceiling $\theta'/2$,
and the excess is carried by the flatness
$F=(\ell^1/\ell^2)/\sqrt{\#k}$, which is still rising. **That rise
cannot continue**: $F\le1$ by Cauchy–Schwarz, so no power of $N$
describes it past the point where it reaches one.

**X1 and X2 hold.** The $\#k$ reproduce exactly and $F$ reproduces
Remark [rem:crosskreference]'s $0.6622$–$0.6854$ on the five $N$ it
published. And $\#k$ is a clean power: its exponent is $0.560706$ with
r.m.s. residual $0.000853$, so the ceiling of $e(\ell^1/\ell^2)$ is
$\theta'/2$ as a fact about the $k$-set and not a fitted quantity.

**X4 holds and X3 is refuted**, and the refutation is empty as a
discrimination. Fitted to eight points,

$$
\begin{array}{l|cc}
\text{shape} & \text{r.m.s.} & \\\hline
F\sim N^{e} & 0.006479 & e=+0.007445,\ F=1 \text{ at } 10^{28.68}\\
F=a+b/\log N & 0.006641 & a=0.754152
\end{array}
$$

The power fits better by $0.000163$, and an r.m.s. from eight points
with two parameters carries a standard error of $0.001870$ — $28.9$
per cent. **The gap is $0.09$ of that.** The two shapes are not
separated by the data at all; what separates them is that
Cauchy–Schwarz forbids one of them. At $\log_{10}N=12$ they still
disagree by only $0.0349$ in $F$.

So the flatness saturates — below $1$, at a fitted $0.754$ — and
asymptotically $e(\ell^1/\ell^2)=\theta'/2$. Item 5 then does reduce
to $e(G)$, and the lean-over-floor exponent tends to
$\theta'/2-e(G)=+0.126442$ against the $+0.159294$ measured over this
range. **The lean still grows; a fifth less fast than the accessible
range shows.**


#### Remark (three exponents, one identity) {#rem:leanidentity}
<!-- evidence: audit_lean_identity.py -->

Remark [rem:leanextended] reports the lean growing against its floor
with no mechanism attached, and it does not need one. Writing
$a_k=(\log k)H(N;k)$, the three quantities in play are ratios of the
same three norms: $G=\ell^1/|\!\sum a|$, the concentration
$\ell^1/\ell^2$, and $\text{lean}/\text{floor}=|\!\sum a|/(c\,\ell^2)$
with $c$ the constant a median sign sum sits on. So

$$
\frac{\text{lean}}{\text{floor}}
= \frac{\ell^1/\ell^2}{G}\cdot\frac1c
$$

identically. **The lean grows against its floor for exactly one
reason: the magnitude concentration outruns the cancellation.**

**W1 holds** — $G$ reproduces the published values to $0.00042$ — and
**W3 holds**: the exponents satisfy
$e(\ell^1/\ell^2)-e(G)-e(\text{lean}/\text{floor})=-0.025407$ against
two standard errors of $0.040473$. Measured over eight $N$,

$$
e(G)=+0.153911,\quad e(\ell^1/\ell^2)=+0.287798,\quad
e(\text{lean}/\text{floor})=+0.159294 .
$$

**So what would have to change is one number.** The lean stops growing
against its floor exactly when $e(G)$ reaches $e(\ell^1/\ell^2)$, and
that is a factor $1.87$ in the exponent — from $0.154$ to $0.288$.

*Added later.* Both exponents are eight-point values on the doubling
family. On the field to $1.024\cdot10^8$ they are $+0.149567$ and
$+0.283586$ — see [rem:fillfield] and [rem:fieldreach]. The identity
and the demand are unchanged; only the numbers move.

**W4 is refuted, narrowly and in the awkward direction.**
$\ell^1/\ell^2$ is bounded by $\sqrt{\#k}$ and $\#k\asymp N^{\theta'}$
here, so its exponent cannot exceed $\theta'/2=0.28$; it measures
$0.287798\pm0.002472$, **$3.15$ standard errors above**. The ratio
$(\ell^1/\ell^2)/\sqrt{\#k}$ is $0.6760,\,0.6622,\,0.6854,\,0.6764,\,
0.6802,\,0.6872,\,0.6986,\,0.6909$ — Remark [rem:crosskreference]
called it flat on five points and on eight it is mildly rising. So
the concentration is still gaining on its own ceiling and item 5 does
not reduce to $e(G)$ alone.

*Added later.* Eight points were a small sample and a mixed field.
On the clean field the excess over the ceiling is $4.98$ standard
errors and on the field to $1.024\cdot10^8$ it is $6.20$, with $F$
resolved rising at $t=5.92$ — [rem:fillfield] and [rem:fieldreach].
The reading here is right and its numbers are superseded.

**W2 is refuted by its cap.** The constant $c$ was asked to hold to
$0.02$ in the log; it runs over $0.219062$. But $c$ is a *median* over
$256$ draws. Splitting them into sixteen groups and measuring the
estimator's own scatter puts the sampling spread of eight such medians
at $0.237831$ — the same size. $c$ does not drift; it is estimated to
about a quarter by $256$ draws, and the cap was set without asking
that.


#### Remark (and it holds over a factor 128) {#rem:leanextended}
<!-- evidence: audit_lean_extended.py -->

Remark [rem:leanfloor] is six points, and this project's record is
that six points are where a trend can be a short-sweep artefact.
The longer lever already existed: Remark [rem:extendrange]'s table
carries the same $f$ to $N=2.56\cdot10^7$, eight octaves. **Its
control does not.** That file runs no null, on the grounds that "the
coin reference level for $f$ is $\tfrac12$ and was measured there" —
in Remark [rem:leandecay]'s evidence, which stops at $6.4\cdot10^6$.
The two largest $N$ in these papers had never had a floor.

**E1 and E2 hold**: the eight mass fractions reproduce to $0.00002$
with the same $\#k$ at every $N$, and the power exponent comes back
at $-0.153911$ against the published $-0.1539$.

**E3 and E4 hold, over the whole factor $128$.** The median lean of
$256$ sign vectors on the identical magnitudes falls at $-0.313205$
($24.21$ standard errors) against $\mu$'s $-0.153911$ — a difference
of $+0.159294$, **$9.29$ standard errors**. In units of that floor
the lean is

$$
8.49,\ 11.31,\ 11.92,\ 13.43,\ 13.25,\ 15.79,\ 16.84,\ 21.36 ,
$$

rising at $+0.159294$ per unit $\log N$, $9.58$ standard errors.

And the two sweeps agree on the thing that was in doubt: the
floor-relative slope is $+0.152263$ on six points and $+0.159294$ on
eight, a difference of $0.22$ standard errors. **The longer lever does
not soften it; it sharpens it.** Over every range measured, the sign
lean grows relative to what chance gives on the same magnitudes.


#### Remark (the floor moved faster than the lean) {#rem:leanfloor}
<!-- evidence: audit_lean_floor.py -->

Remark [rem:leandecay] reads the lean's decay as the thing that
rescues the asymptotic picture, and reads its coin arm — two draws per
$N$ — only for whether it "sits at $\tfrac12$ throughout, as it must".
Sitting at $\tfrac12$ is not the question. A random sign field on
these magnitudes leans by about $\ell^2/2\ell^1$, and Remark
[rem:crosskreference] has just measured $\ell^1/\ell^2$ growing, so
that floor **moves**.

**L1 and L2 hold.** The six mass fractions come back to $0.00042$ and
the raw decay is resolved: slope $-0.167257$ at $9.00$ standard
errors, reproducing the published $-0.1673$.

**But the floor falls faster.** With $256$ global sign vectors on the
identical magnitudes, the median coin lean runs
$0.0276,\,0.0237,\,0.0192,\,0.0162,\,0.0116,\,0.0095$ and its slope
is $-0.315933$ at $17.68$ standard errors — steeper than $\mu$'s by
$5.77$ standard errors. **L4 holds** and says the same thing
directly: $\mu$'s lean measured in units of its own floor is

$$
9.86,\ 11.71,\ 11.79,\ 11.95,\ 15.44,\ 17.14 ,
$$

rising at $+0.148676$ per unit $\log N$, $5.81$ standard errors.

**L3 is refuted by its instrument, not its direction.** It asked
$\mu$'s slope to sit above the $97.5$th percentile of the $256$
per-draw slopes; a single draw's six leans are so noisy that those
slopes run from $-1.63$ to $+0.63$, and no fixed quantity could clear
that percentile. The floor is not one draw but where the draws sit,
and that is what the median series estimates.

What survives and what does not. The raw table is right, $G=1/|2f-1|$
does rise, and Remark [rem:leandecay]'s first consequence — that
Remark [rem:nocrossk] is a statement about the accessible range — is
untouched. What does not follow is the second reading. **Over this
sweep the lean does not go away relative to chance; it grows.** That
the raw lean shrinks is the floor moving under it.


#### Remark (what independent signs would actually give) {#rem:crosskreference}
<!-- evidence: audit_crossk_reference.py -->

Remark [rem:nocrossk] measures the cross-$k$ gain at $1.834$ to
$2.789$ and says "independent signs would give $\sqrt{\#k}$" —
$17.7$ to $38.5$. **$\sqrt{\#k}$ is what independent signs give on
equal magnitudes**, and these are not equal: that remark's own rule
T4 put the top decile at $0.3486$–$0.3587$ of the mass. For unequal
magnitudes the reference is $\ell^1/\ell^2$, which is below
$\sqrt{\#k}$ by exactly the concentration and had never been computed.

**T1 holds**, the gain reproducing to $0.00042$ with the same $\#k$ at
every $N$. Measured,

$$
\frac{\ell^1/\ell^2}{\sqrt{\#k}}
= 0.6622\ \text{to}\ 0.6854 ,
$$

so the deficit factors, with no distributional constant needed on
either side:

$$
\frac{\sqrt{\#k}}{G}
= \underbrace{1.459\text{–}1.510}_{\text{concentration, flat}}
\ \times\
\underbrace{6.52\text{–}9.40}_{\text{correlation, rising}} .
$$

**T3 holds and is the one that matters**: random signs on $\mu$'s own
magnitudes give a median gain $9.94$ to $12.98$ times $\mu$'s, so a
correlation across dilations is there and Remark [rem:nocrossk]'s
conclusion stands. **T4 holds** — the coin's median gain over
$\ell^1/\ell^2$ runs $1.26$ to $1.59$, straddling the $1/\,\mathbb
E\text{-median}|Z|$ that a median-convention gain must sit on.

**T2 is refuted, and by its own mixing of conventions.** It asked the
median coin gain to stay below $\sqrt{\#k}$; a median gain carries a
factor the equal-magnitude $\ell^1/\ell^2$ does not, so at two $N$ it
sits just above. The factorisation above is the same comparison made
in one convention throughout, and it is the statement to keep.

What changes: $n_{\mathrm{eff}}=G^2$ should not be read against
$\#k$. Against the coin's own $n_{\mathrm{eff}}$ on the same
magnitudes the walls move together by $98.8$ to $168.5$, not by the
$313$ to $1485$ the published comparison suggests. The conclusion is
unchanged and its size is roughly halved in the gain.


#### Remark (mu ranks fourth of seventeen) {#rem:coinrank}
<!-- evidence: audit_residue_coin_rank.py -->

Remark [rem:residuecancel]'s conclusion — "nothing in $R$ is doing
better than random signs, so no argument that $\mu$ is special against
$\delta$ can help" — rests on a band, and a band is not a test. The
published ratio divides $\mu$'s octave mean by the mean of sixteen
coins, and the scatter of that mean was never computed; nor was the
sign counted, and reading down the column most entries are below $1$
with all of them below at the largest $N$.

Since the sixteen sign vectors are global, one per draw and held
across every $k$ exactly as $\mu$ is, the seventeen are exchangeable
under the hypothesis and $\mu$'s rank among them is uniform. That is
an exact test.

**V1 and V2 hold.** The $\ell^1$ exponents come back at
$0.9745,\,0.9919,\,0.9962,\,0.9963,\,0.9979$ and every octave mean of
$|R|/\ell^2$ to within $0.00005$; the band this measurement gives is
$[0.8195,\,1.2871]$, inside the published one.

**V3 holds.** Pooling each draw against the mean of the other sixteen
and averaging over the thirty octaves, $\mu$ scores $-0.033137$ and
three of the coins score lower: **rank four of seventeen.** $\mu$ is
inside the coins' own range, and the conclusion stands.

**V4 is refuted, and what it caught was its own null.** $\mu$ falls
below the coins' median at $21$ of $30$ octaves against an expected
$15.0$ — $2.2$ binomial sigma. But the octaves are not independent:
the same sixteen vectors serve every octave and every $N$, and
octaves at different $N$ run over overlapping $m$. Counting the same
statistic for each coin gives $3,\,9,\,9,\,10,\,11,\,11,\,12,\,12,\,
14,\,16,\,17,\,18,\,19,\,20,\,22,\,24$: **two coins lean at least as
much as $\mu$**, so $21$ is a three-in-seventeen event and not a
$2.2$-sigma one. The magnitude test and the sign test agree once the
dependence is taken out.

What survives is the sentence as written, now with a rank behind it
rather than a band. What does not survive is reading the column of
sub-unit ratios as a systematic advantage: at rank four of seventeen
$\mu$ leans a little, and so do several coins.


#### Remark (the residue buys exactly a coin's cancellation) {#rem:residuecancel}
<!-- evidence: lab_residue_cancellation.py -->

$|R|\asymp(N/k)^{1/2}$ is compatible with two opposite worlds and the
remark above does not separate them. Writing
$R(N;k)=\sum_m\mu(m)\delta(m,k)$ with
$\delta=\Lambda(N-mk)-\beta w(m,k)$ — the deviation of the von Mangoldt
weight from its own sieve prediction — either $\delta$ is already tiny
and no cancellation occurs, or $\sum_m|\delta|$ is of full size and
$\mu$ cancels it down to the square root. Only the second is hard, and
the two are told apart by comparing $|R|$ with the $\ell^1$ and
$\ell^2$ norms of its own summands.

**It is the hard world, and $\mu$ buys exactly a coin's worth.** The
$\ell^1$ norm is linear: fitted against $N/k$ its exponents are
$0.9745,\,0.9919,\,0.9962,\,0.9963,\,0.9979$ with correlations to
$1.00000$. Against that, $|R|/\ell^2$ has octave means running $0.91$
down to $0.67$ — and sixteen global sign vectors summed against the
*identical* $\delta$ run $0.93$ down to $0.68$, so $\mu/\text{coin}$
stays in $[0.78,\,1.33]$ and sits on $1$. The scale both are near is
$\sqrt{2/\pi}=0.7979$, which is what $\mathbb E|Z|/\sigma$ is for a
random sign sum. The control $\beta$ reproduces the published split to
six decimals.

The aggregate gain is therefore square-root on the nose:

$$
\frac{\overline{\ell^1}}{\overline{|R|}}\ \asymp\ (N/k)^{e},\qquad
e=0.5024,\,0.4997,\,0.4883,\,0.4915,\,0.5109
$$

with correlations $0.99610$ to $0.99970$. (Rule V4 asked this of
$\overline{\ell^1/|R|}$ and **is refuted**; its instrument is a mean
of a ratio whose denominator passes near zero, so a single $k$ at
which $R$ nearly vanishes sets the octave mean, and its correlations
wander in sign from $-0.90$ to $+0.98$. The aggregate has no such
denominator and is what the rule meant to ask.)

Two things follow, and both are negative in a useful way. Nothing in
$R$ is doing better than random signs, so **no argument that $\mu$ is
special against $\delta$ can help** — the route cannot be rescued by
finding extra structure here, because there is none to find. And
nothing is doing worse, so the square-root heuristic that
Remark [rem:directlevel] and Remark [rem:heuristic] price the level
with is not optimistic about $R$; it is exact. What remains is a
Möbius sum against a bounded, essentially two-valued deviation
sequence, achieving precisely square-root cancellation, at level $k$
up to $N^{\theta'}$ — and by Remark [rem:provablehalf] it is the half
for which no unconditional estimate exists at all.


#### Remark (the residue has no law, only a fit) {#rem:residueconstant}
<!-- evidence: audit_residue_constant.py -->

Remark [rem:heuristic] pinned $H$'s constant: $|H|\approx
c(N)\sqrt{N/k}$ with $c(N)/\sqrt{\log N}$ reading
$1.0138,\,0.9854,\,1.0039,\,0.9939,\,0.9844$ — a law to a percent and
a half, and enough to predict $K^*_H$ to $1.5\%$. Whether the
knife-edge of Remark [rem:residuelevel] is local or permanent is a
question about the same constant for $R$, which had never been
measured.

**$R$ has no such law.** Measured below the operative crossing, the
control reproducing all five $K^*_R$ exactly,

$$
c_R(N)=1.5925,\ 1.7303,\ 1.9439,\ 1.9894,\ 1.9188,
$$

and divided by $\sqrt{\log N}$ these are
$0.4558,\,0.4818,\,0.5273,\,0.5263,\,0.4958$ — **a spread of $14.4\%$
of the mean, non-monotone, where $H$'s was under $3\%$ and flat.**
Rule Q2 asked for $5\%$ and is refuted.

What survives is a fit and not a law, and the distinction is
operative. Q3 holds and holds well: with $c_R$ taken at each $N$, the
model $\sum_{k<K}(\log k)c_R\sqrt{N/k}=\SS(1-A)N$ reproduces the
crossing to $0.998,\,1.000,\,1.007,\,0.994,\,1.002$. Q4 holds too —
the model's own exponents $0.5652,\,0.5642,\,0.5604,\,0.5671,\,0.5800$
track the measured $0.5654,\,0.5642,\,0.5599,\,0.5675,\,0.5799$ to
five decimal places at four of five $N$, and put $N=8\cdot10^5$ at
$0.5604$ where the measurement reads $0.5599$: **the dip that refuted
U3 of Remark [rem:residuelevel] sits at the resolution of the model,
on the barrier rather than under it.**

**But no forecast is made, and the script refuses to make one.** The
extrapolation was pre-registered as conditional on Q2, and Q2 failed,
so the model may interpolate and may not be projected. That is the
finding. Remark [rem:residuelevel] established that the margin over
$\tfrac12$ is not closing over the accessible range; this establishes
that **nothing here entitles anyone to say where it goes.** The
$+0.0047$ slope of the exponent is measured over a factor $16$ in $N$
against a constant that itself wanders by $14\%$ over the same range.

One number is worth keeping. $c_R\approx0.50\sqrt{\log N}$ against
$H$'s $c\approx1.00\sqrt{\log N}$: **the split halves the constant and
buys no exponent**, which is Remark [rem:elemsize]'s conclusion seen
from the other side, and is why removing the elementary half moves
$\theta'$ by the $0.06$ that Remark [rem:splitbudget] measured and not
by more.


#### Remark (the lean is finite-$N$ and decays) {#rem:leandecay}
<!-- evidence: lab_lean_decay.py -->

The table above carries a hint it was not read for: the mass fraction
rises with $N$. Extending by one octave,

$$
\begin{array}{r|cccccc}
 N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6
   & 3.2\cdot10^6 & 6.4\cdot10^6\\\hline
 f & 0.2273 & 0.2228 & 0.2735 & 0.3068 & 0.3207 & 0.3376\\
 |{\tfrac12}-f| & 0.2727 & 0.2772 & 0.2265 & 0.1932 & 0.1793 & 0.1624\\
 G = 1/|2f-1| & 1.834 & 1.804 & 2.207 & 2.588 & 2.789 & 3.079\\
 \text{coin } f & 0.5053 & 0.5124 & 0.5137 & 0.4928 & 0.4619 & 0.4720
\end{array}
$$

The coin sits at $\tfrac12$ throughout, as it must; $\mu$'s lean
shrinks by $0.11$ over a factor $32$ in $N$, and $G$ rises from
$1.834$ to $3.079$. **So the lean is not a structural fact about
$\mu$: it is the finite-$N$ error of Theorem [thm:C], the same
error Remark [rem:thetasweep] measured at $0.17$ to $0.46$ of $N$,
and it decays.**

Two things follow, and the first corrects the paragraph above.
Remark [rem:nocrossk]'s "there is no cancellation across dilations"
is a statement about the accessible range, not an asymptotic one:
$G=1/|2f-1|\to\infty$ as $f\to\tfrac12$, so the cancellation does
arrive and the one-sidedness of Proposition [prop:onesided] buys a
*growing* factor rather than the constant near two that the accessible
range shows.

The second is that the rate is not determined here, and the two
candidates are far apart. A power law gives
$|\tfrac12-f|\asymp N^{-0.1673}$ with correlation $-0.97616$; a log
law gives $(\log N)^{-2.3166}$ with correlation $-0.97469$. Over a
factor $32$ in $N$ they are not separable, and they reach
$|\tfrac12-f|<0.01$ at $N=1.02\cdot10^{14}$ and at
$N=4.3\cdot10^{22}$ respectively. Since the lean *is* the error term
of Theorem [thm:C], and that error is a power of $\log$, the log law
is the more natural of the two — but nothing measured here chooses
between them, and the projected gains $G\approx3.4,\,7.3,\,23.1$ at
$N=10^7,10^9,10^{12}$ are the power law's and should be read as its
upper end. Remark [rem:leanbracket] measures how wide "not
separable" is.


#### Remark (the lean, by a second route, and its bracket) {#rem:leanbracket}
<!-- evidence: audit_lean_bracket.py -->

Two things are owed here. The lean above is computed from the inner
sum $H(N;k)=\sum_{m<N/k,(m,k)=1}\Lambda(N-mk)\mu(m)$ directly, and a
statistic that carries a theorem's error term should not rest on one
implementation; and the $1.02\cdot10^{14}$ is a forecast eight
million-fold beyond the computed range with no bracket on it.

The first settles cleanly. Recomputing $f$ from the dilation identity
[eq:dilate], $H(N;k)=\mu(k)A(N;k)$ with
$A(N;k)=\sum_{n\equiv N\ (k)}\Lambda(n)\mu(N-n)$ — a strided sum over
a different array, sharing no arithmetic with the inner sum beyond the
sieve — reproduces all six published values
$0.2273,\,0.2228,\,0.2735,\,0.3068,\,0.3207,\,0.3376$ **to every
printed digit**, and the refitted exponent agrees to $0.0000$. The
lean is not an artefact of how $H$ was summed.

The second does not. Refitting the power law on each leave-one-out
subset moves the forecast over $10^{13.17}$ to $10^{14.01}$ — only
$0.84$ decades, and *less* than the naive $1.99$ the $12.42\%$ spread
in the exponent would suggest, because the leave-one-out refits the
intercept alongside the exponent and a steeper decay comes with a
smaller constant. **The law choice has no such brake.** The log law
lands at $10^{22.64}$, $10^{18.98}$ and $10^{22.38}$ on the same three
subsets, and the total bracket is

$$
|{\tfrac12}-f|<0.01 \ \text{ at }\ N\in[10^{13.17},\,10^{22.64}],
$$

**nine and a half decades.** So $1.02\cdot10^{14}$ is not a forecast
and must not be quoted as one; what the data support is that the lean
decays and that $G\to\infty$, with the scale of the approach undecided
across nine decades.

Those nine decades are not an error bar, and the distinction matters
enough to state. Within a law the constant extrapolated is the decay
exponent $b$, and its own leave-one-out spread — $0.1673$ to
$0.1880$, a relative $0.1180$ — is what the narrow bracket is built
from, so bracket and drift are the same number here as in Remark
[rem:marginbracket]. The wide bracket is a **choice between shapes**:
the power law and the log law are not separable over the accessible
range, and no amount of care about a fitted constant addresses that.
Reporting it as one interval is a convenience for the reader and
should not be read as precision about anything. The projected $G$ at $10^7,10^9,10^{12}$ inherit
the same width and are the power law's end of it, as
Remark [rem:leandecay] already says.


#### Remark (two more octaves, and what they settle) {#rem:extendrange}
<!-- evidence: lab_extend_range.py -->

Three questions in this section were limited by the same thing, the
length of the $N$-range, so the computation was extended once to
$N=2.56\cdot10^7$ — a factor $128$ rather than $32$:

$$
\begin{array}{r|cccccccc}
 N & 2\!\cdot\!10^5 & 4\!\cdot\!10^5 & 8\!\cdot\!10^5
   & 1.6\!\cdot\!10^6 & 3.2\!\cdot\!10^6 & 6.4\!\cdot\!10^6
   & 1.28\!\cdot\!10^7 & 2.56\!\cdot\!10^7\\\hline
 B/N & 0.8086 & 0.7395 & 0.7303 & 0.6547 & 0.5916 & 0.5526 & 0.4992
   & 0.4527\\
 |E_3|/N & 0.4377 & 0.3837 & 0.3172 & 0.2608 & 0.2073 & 0.1855
   & 0.1459 & 0.1245\\
 f & 0.2273 & 0.2228 & 0.2735 & 0.3068 & 0.3207 & 0.3376 & 0.3533
   & 0.3608\\
 G & 1.834 & 1.804 & 2.207 & 2.588 & 2.789 & 3.079 & 3.407 & 3.592
\end{array}
$$

**The decay law is still not settled, and that is the finding.**
Over the doubled range the power fit is
$|\tfrac12-f|\asymp N^{-0.1539}$ with correlation $-0.98434$ and the
log fit $(\log N)^{-2.2380}$ with $-0.98619$ — a difference of
$0.00185$, with the log law now marginally ahead, as theory says it
should be. Quadrupling the range moved the discrimination by less than
$0.002$ in correlation. Extending $N$ is not the way to settle it, at
any scale this method reaches.

The out-of-sample check on $B(N)/N$ — the audit's rule X3 — fails,
and its direction is
useful. The $(\log N)^{-1.4526}$ law fitted on the five smallest $N$
over-predicts the three it never saw by $2.7\%$, $6.8\%$ and
$11.1\%$ — monotone, so systematic: $B/N$ falls *faster* than that
law. Refitting on all eight gives $(\log N)^{-1.7267}$ crossing the
Goldbach threshold $0.3745$ at $N=10^{8.44}$, or
$N^{-0.1196}$ crossing at $10^{8.17}$. So Remark [rem:relocate]'s
$10^{8.9}$ was the conservative end and **[eq:nolog] is projected to
become true between $10^{8.2}$ and $10^{9.1}$** — two orders above the
computable range, not twenty. Remark [rem:forecast], reaching the
same region by an entirely different route, put the $\theta'=0.56$
level crossing at $2.08\cdot10^8=10^{8.32}$, inside that bracket.

And the one-sided condition is no longer marginal: $|E_3|/N$ is
$0.1459$ and $0.1245$ at the two new $N$ against a threshold of
$0.3745$, so [eq:onesided] now holds with a factor of three to spare
and is still falling.


#### Remark (the decay exponent is not determined, and by how much) {#rem:decayfamily}
<!-- evidence: lab_decay_family.py -->

Extending $N$ will not settle the decay law, so the family is fitted
instead. The two candidates compared so far are special cases of

$$
\bigl|\tfrac12-f\bigr| \;=\; A\exp\!\bigl(-c\,(\log N)^{\alpha}\bigr),
$$

with $\alpha=1$ the power law and $\alpha\to0$ the log law — and
theory names a third. Remark [rem:thetasweep] found the finite-$N$
residual is dominated by the main-term cancellation over
$m<M=N^{1-\theta'}$, which Lemma [lem:mu] bounds by
$\exp(-c\sqrt{\log M})$: that is $\alpha=\tfrac12$, neither a power of
$N$ nor a power of $\log N$, and it had never been tested.

Swept over $\alpha\in[0.05,1.50]$ the residual sum of squares moves by
$23\%$ in total — $1.0096,\,1.0289,\,1.0554,\,1.0891,\,1.1300,\,
1.1780,\,1.2330$ relative to the minimum at
$\alpha=0.20,0.40,\dots,1.40$. **The theoretically named
$\alpha=\tfrac12$ costs $4.12\%$ over the best fit**, which is nothing
on eight points. And the same sweep applied to $B(N)/N$ prefers the
*opposite* end, $\alpha=1.50$, with $\alpha=\tfrac12$ costing $75\%$
there: the two quantities do not even agree on which way the
preference runs.

The audit's rule Y2 asked whether the data determines $\alpha$, using
a band of $1\%$ in residual sum of squares, and by that band it does
— the admissible set spans only $0.15$. **That band is my threshold
and not a confidence band**, and eight points with two fitted
parameters do not resolve $13\%$ of RSS. Under the proper $95\%$ band
for one parameter with six degrees of freedom, RSS within a factor
$2.00$, the admissible $\alpha$ is the whole sweep from $0.05$ to
$1.50$. The rule is recorded as refuted and the diagnosis with it;
this is the fourth threshold in this work stated as an effect size
rather than against a null, after Remarks [rem:cap], [rem:band]
and [rem:filter].

What that costs is quantified rather than waved at. The three
canonical laws put $|\tfrac12-f|<0.01$ at $N=10^{14.71}$,
$10^{17.46}$ and $10^{23.74}$ — a span of $9.02$ orders. **So no
crossing scale for the lean should be quoted**, and the projected
gains at the end of Remark [rem:leandecay] are one law's among
three. The bracket $10^{8.2}$ to $10^{9.1}$ of
Remark [rem:extendrange] is a different matter and survives: it
extrapolates $B/N$ by a factor of at most forty beyond the computed
range, where this one extrapolates by a factor $10^{7}$.


#### Remark (the dilated wall obeys a square-root law, and what that would buy) {#rem:dilateprofile}
<!-- evidence: lab_dilate_profile.py -->

With [eq:dilate] the demand has one governing profile, the relative
size of the dilated wall,

$$
\rho(k) \;:=\; \bigl|A(N;k)\bigr|\,\frac{k}{N}
  \;=\; \bigl|H(N;k)\bigr| \Big/ \frac{N}{k}.
$$

At $N=3.2\cdot10^6$, $\rho(k)/\sqrt k$ reads, as an octave median,

$$
\begin{array}{r|ccccccc}
 k & [64,128) & [256,512) & [1024,2048) & [4096,8192)
   & [8192,16384)\\\hline
 \rho/\sqrt k & 0.001234 & 0.001357 & 0.001268 & 0.001376 & 0.001512
\end{array}
$$

— flat to about $\pm10\%$ across eight octaves. Fitted, the exponent
in $\rho(k)\asymp k^{a}$ is $a=0.4431,\,0.4528,\,0.4686,\,0.4962,\,
0.5826$ for $\mu$ at $N=2\cdot10^5$ through $3.2\cdot10^6$, against
$0.3833,\,0.5284,\,0.5081,\,0.5210,\,0.4968$ for a coin on the same
support. **The two share the exponent $\tfrac12$ and differ only in
the constant** — by the factor $1.08$ to $1.53$ that
Remark [rem:whycoinwins] measures.

What $a=\tfrac12$ would buy is worth stating, because it is not small.
Writing $\rho(k)=c'\sqrt{k/N}$,

$$
\frac{B(N)}{N} \;=\; \sum_{k<K}(\log k)\,\frac{\rho(k)}{k}
  \;\approx\; 2c'\sqrt{K/N}\,\log K ,
$$

so $B(N)\le\SS(N)(1-A(N))N$ holds up to $K\asymp N/(\log N)^{2}$ if
$c'$ is a constant — an exponent of $1$ in $N$, where Huang–Li need
only $\theta'>1/2$. The constant is not constant, and
Remark [rem:extrap] carries the count of $\log$s through properly;
the exponent $1$ survives. So **the level requirement is not the
obstruction that a square-root law would leave.** The obstruction is
that no proof gives square-root cancellation for a dilated
$\Lambda$–$\mu$ correlation: that is the same shape as the $E1$
consumable of \S[sec:supply], and of the same difficulty class. The
demand side's remaining requirement has been carried onto the supply
side's ground.

Two of the audit's rules fail, and both to a mis-specified field. Q2
asked the coin's exponent to sit within $0.05$ of $\tfrac12$ and Q3
asked $\mu$'s to exceed the coin's; both fail at $N=2\cdot10^5$,
where the fit ran to $k<2\cdot10^4$ and so let the dilate length
$M=N/k$ fall to $10$ — a range in which a law in $k$ has nothing to
describe. Refitting over $M\ge1000$ gives the figures quoted above,
and there Q2 holds at every $N$ with enough moduli. Q3 does not
recover, and should not have been asked: the two exponents agree and
the difference between $\mu$ and a coin is in the constant, which is
what Remark [rem:whycoinwins] already showed. The extrapolation to
$K\asymp N/(\log N)^2$ is a heuristic from a fitted exponent and not a
theorem; the directly measured $K^*$ exponent over a factor $16$ in
$N$ was $0.7057$, which is what that polylog looks like from here.


#### Remark (the extrapolation, tested out of sample) {#rem:extrap}
<!-- evidence: lab_dilate_extrapolation.py -->

The step from $\rho(k)\asymp\sqrt k$ to $K\asymp N/(\log N)^2$ above
assumed that $c'=|A(N;k)|/\sqrt{N/k}$ does not depend on $N$, and
nothing tested it. It does depend on $N$, and the reason is visible in
the control: the summands $\Lambda(N-mk)\varepsilon(m)$ are nonzero at
about $M/\log N$ places with size $\log N$ each, so
$|H_\varepsilon|\asymp\sqrt{M\log N}$ and
$c'_\varepsilon\asymp\sqrt{\log N}$. Measured,
$c'_\varepsilon/\sqrt{\log N}$ reads
$0.5619,\,0.6034,\,0.6065,\,0.5861,\,0.6010$ — flat to $1.079$, which
is the calibration this remark rests on. Carrying that power through,

$$
\frac{B(N)}{N} \;\approx\; 2\gamma\sqrt{\log N}\,\sqrt{K/N}\,\log K ,
\qquad
K^* \;\asymp\; \frac{N}{\log N\,(\log K)^{2}} ,
$$

**three powers of $\log$, not two.** The exponent $1$ in $N$ is
unchanged, so nothing in the conclusion moves.

The model is then tested where it can be: $\gamma$ is fitted on
$N=2\cdot10^5,\,4\cdot10^5,\,8\cdot10^5$ alone and used to predict
$K^*$ at the two larger $N$ it never saw.


#### Remark (that prediction was an artifact; the corrected one is not) {#rem:artifact}
<!-- evidence: lab_level_forecast.py -->

The prediction just described was first reported as $1599$ against a
measured $1353$ and $3199$ against $2319$, ratios $1.182$ and
$1.379$. **It is withdrawn: those were not solutions.** The search
for $K^*$ ran over $k\le N/1000$, and $1599$ and $3199$ are exactly
the largest admissible $k$ below $1600$ and $3200$. The cumulative
sum never reached the threshold inside the range, and the routine
returned its last index. A prediction that lands on the endpoint of
its own search range is the endpoint, not a prediction. That was rule R4
of the extrapolation audit, and with the boundary now detected it
reads REFUTED there, as it always should have.

Searched properly, to $K=10^7$, the model with $\gamma=0.6520$ gives
$597,\,937,\,1483,\,2391,\,3903$ against the measured
$319,\,537,\,767,\,1353,\,2319$ — it overshoots by
$1.871,\,1.745,\,1.934,\,1.767,\,1.683$, a steady factor near $1.8$,
so it *understates* $B(N)$. That was the audit's rule S1 and it
fails.

The cause is a calibration error and it is instructive. $\gamma$ was
taken as the **median** of $|A(N;k)|/\bigl(\sqrt{N/k}\sqrt{\log N}\bigr)$,
but $B(N)$ is a **sum**, which needs the mean, and $|A|$ is
one-sided so the mean is the larger: measured, the ratio is
$1.4501,\,1.5899,\,1.4832,\,1.4473,\,1.4712$. With the mean-based
$\gamma=0.9803$ — still fitted on the three smallest $N$ alone — the
model gives

$$
\begin{array}{r|ccccc}
 N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6
   & 3.2\cdot10^6\\\hline
 K^*\ \text{model} & 329 & 569 & 779 & 1271 & 2171\\
 \text{model}/\text{measured} & 1.031 & 1.060 & 1.016 & 0.939 & 0.936
\end{array}
$$

within $6.4\%$ at every $N$, including the two it never saw. **That
is the out-of-sample validation**, and it is the corrected one that
Remark [rem:dilateprofile]'s conclusion rests on. Nothing in the
shape of that conclusion moves — the exponent in $N$ is still $1$ and
the polylog is still three — but the constant, and every number
extrapolated from it, does.


#### Remark (where the model puts the route's hypothesis) {#rem:forecast}
<!-- evidence: lab_level_forecast.py -->

With the corrected $\gamma$ the model can be asked the question no
computation reaches: at which $N$ does
$B(N)\le\SS(N)(1-A(N))N$ hold at the level Huang–Li need? Solving
$B(N)/N=\gamma\sqrt{\log N/N}\,S(K)$ against the threshold, with
$S(K)=\sum_{k<K}(\log k)k^{-1/2}$ over the admissible $k$ enumerated
exactly to $10^7$:

$$
\begin{array}{r|cccccc}
 N & 10^6 & 10^7 & 10^8 & 10^9 & 10^{12}\\\hline
 K^*/\sqrt N & 0.9570 & 1.4904 & 2.5967 & 4.9203 & 47.1443\\
 K^*/N^{0.56} & 0.4177 & 0.5666 & 0.8598 & 1.4190 & 8.9832
\end{array}
$$

so $K^*$ passes $\sqrt N$ at $N=1.299\cdot10^6$ and passes $N^{0.56}$
at $N=2.077\cdot10^8$ — each to be read with the bracket Remark
[rem:forecastbracket] computes, two-thirds of a decade wide. The first
of those is a check rather than a
forecast, and it passes: the measured $K^*/\sqrt N$ crossed $1$
between $8\cdot10^5$ and $1.6\cdot10^6$, which the model reproduces
without having been shown it. Rule S2 of the forecast audit asked
for that crossing under the *uncorrected* $\gamma$ and is refuted
there — under it $K^*/\sqrt N$ already exceeds $1$ at $N=10^5$,
which is the same overshoot rule S1 records.

Three things this is not. It is a forecast from a model fitted on
three values of $N$ spanning a factor $4$, extrapolated sixty-fold
beyond the largest of them. It is about $B(N)$, the aggregate
absolute sum that Proposition [prop:nolog] needs, and not about
$EH_\mu$, which asks for a maximum over residue classes as well. And
— the caution that Remark [rem:levelmeas] was withdrawn for not
taking — **$B(N)$ does not distinguish $\mu$ from a coin**: by
Remark [rem:whycoinwins] a coin gives a slightly *smaller* $B$ and
so a slightly earlier crossing, so this forecast is a statement about
the size of the object and about square-root cancellation, not about
any structure special to $\mu$. What it says is that the size is not
the obstruction, and where the size stops being one.

It also settles the gap Remark [rem:dilateprofile] left open between
an asymptotic exponent of $1$ and a measured $K^*$ exponent of
$0.7057$. The approach is slow: the effective exponent of
$N/(\log N)^3$ is $0.4298$ at $N=10^6$, $0.5612$ at $10^9$ and
$0.6397$ at $10^{12}$. **A local slope of $0.7$ over a factor $16$ in
$N$ is not evidence against an asymptote of $1$; it is what the
asymptote looks like from this range**, and no accessible computation
distinguishes the two.

Two pre-registered rules fail, and they fail because the range is too
short. R1 asked $c'(N)$ to increase; it reads
$2.3951,\,2.1060,\,2.5220,\,2.6162,\,2.5291$ and wobbles. R2 asked
$c'/\sqrt{\log N}$ to be flat to $15\%$; its spread is $1.1805$. Over
a factor $16$ in $N$, $\sqrt{\log N}$ moves by only $11\%$, which is
below the scatter in $c'$ — so $\mu$'s own data cannot separate
"$c'$ constant" from "$c'\asymp\sqrt{\log N}$", and the extra power of
$\log$ is carried by the coin calibration and by the count of
non-zero summands, not by $\mu$.


#### Remark (the bracket that forecast has to carry) {#rem:forecastbracket}
<!-- evidence: audit_forecast_bracket.py -->

Remark [rem:forecast] quotes one number, $N=2.077\cdot10^8$, and
Remark [rem:modeltransfer] has since shown that this family of models
fits its constant rather than deriving it — a $5\%$ drift in that
constant moved a measured crossing by $10\%$. A forecast extrapolated
sixtyfold past its calibration inherits that. Rebuilt independently
here, three things hold and the fourth, which was the point, does not.

The reconstruction is exact: solving the published model from an
independently sieved $S(K)$ reproduces every decade of the published
table from $10^5$ to $10^{20}$, in both $K^*/\sqrt N$ and
$K^*/N^{0.56}$, to $0.0001$. The constant is not a free parameter:
the forecast fits $\gamma=0.9803$ on the mean of $|A|/\sqrt{N/k}$,
and Remark [rem:heuristic] measures the same constant from a different
statistic over a different $k$-range at
$1.0138,\,0.9854,\,1.0039,\,0.9939,\,0.9844$ — **two independent
measurements agreeing to $3.4\%$**. And the fit is doing almost no
work: setting $\gamma=1$ exactly moves the two crossings only to
$1.5370\cdot10^6$ and $2.4181\cdot10^8$, by $18\%$ and $16\%$.

**But W4 fails, and by a factor of two and a half.** Perturbing the
constant by the $\pm10\%$ that [rem:modeltransfer] measured moves the
$\theta'=0.56$ crossing over $[9.0996\cdot10^7,\,4.2786\cdot10^8]$, a
span of $4.7019$; the $\sqrt N$ crossing spans $5.3909$. In logs:

$$
K^*=\sqrt N:\ 10^{6.11}\ [10^{5.72},10^{6.45}],\qquad
K^*=N^{0.56}:\ 10^{8.32}\ [10^{7.96},10^{8.63}].
$$

**Two-thirds of a decade, not the factor $2$ pre-registered.** The
reason is that the square root of $K$ doubles the exponent: at the
crossing $K^*=N^e$ the balance reads
$N^{(1-e)/2}\sim\gamma\sqrt{\log N}\,(e\log N-2)$, so a relative error
$d$ moves $N$ by $(1+d)^{2/(1-e)}$ — $1.4641$ and $1.5422$ at the two
exponents — and even that understates the measured $2.1814$ and
$2.0600$, because the right-hand side *grows* with $N$ and so $N$ must
move further still. **The logarithms amplify the uncertainty rather
than damping it.**

The bracket is honest where it can be checked: the measured
$K^*/\sqrt N$ crossing lies between $8\cdot10^5$ and $1.6\cdot10^6$,
and both the point estimate and the bracket contain it. **And it is
honest about itself**, which is the harder test: a bracket made by
wobbling a constant is worth nothing if the constant drifts by more
than the wobble. The $\pm10\%$ assumed here is against a measured
drift of $0.0295$ — the spread of the five independent readings of
$c(N)/\sqrt{\log N}$ — so the assumed wobble is three times the real
one. Remark [rem:residueconstant] is the case where that comparison
goes the other way and the forecast is therefore refused. What changes
is what may be said. Remark [rem:forecast]'s "$K^*$ passes $N^{0.56}$
at $N=2.077\cdot10^8$" is a statement about an order of magnitude and
must be written as one; the qualitative reading it draws — that the
size is not the obstruction, and roughly where it stops being one —
survives untouched, since two-thirds of a decade does not reach any
computable range.


## The cell floor {#sec:floor}


Throughout, cells are indexed by *depth* $d$, the number of
$3,5,7,11,13$ dividing $N$, and $Z(N)=C(N)/\sqrt{V(N)}$.


### The floor, in closed form


#### Lemma (exact cell moments) {#lem:cellmom}
<!-- evidence: lab_cellmom_montecarlo.py -->

Let a band of even $N$ be partitioned into cells, let $c$ be a cell of
size $n_c$ inside a band of size $n$, let $m_c$ and $\overline m$ be
the means of $Z$ over $c$ and over the band, and set

$$
u_c(v) \;=\; \sum_{N\in c}\frac{\Lambda(N-v)}{\sqrt{V(N)}},
  \qquad
  Q_{cd} \;=\; \sum_v \mu^2(v)\,u_c(v)\,u_d(v).
$$

Then, when the $\mu(v)$ on the surviving support are replaced by
independent signs,

$$
\mathrm{Var}\bigl(m_c-\overline m\bigr)
  \;=\; \frac{Q_{cc}}{n_c^{2}}
  \;-\; \frac{2\,Q_{ca}}{n_c\,n}
  \;+\; \frac{Q_{aa}}{n^{2}},
$$

where $a$ denotes the whole band. Every term is computable exactly by
one convolution; no simulation is involved.


**Proof.** 
Under independent signs $E[\varepsilon(v)\varepsilon(v')]
=\mu^2(v)\delta_{vv'}$, so
$\mathrm{Cov}(n_c m_c, n_d m_d) = Q_{cd}$ exactly. Expanding
$\mathrm{Var}(m_c-\overline m)$ bilinearly gives the three terms.
 ∎


The lemma does not need simulation, but the simulation is the only
thing that would catch an error in the formula, and every $z_c$ in this
paper is divided by it. Run in the band $(10^5,2\cdot10^5]$:

$$
\begin{array}{r|cccccc}
 \text{depth} & 0 & 1 & 2 & 3 & 4 & 5\\
 n_c & 19181 & 21097 & 8183 & 1421 & 115 & 3\\\hline
 \text{closed form} & 0.140080 & 0.052724 & 0.167579 & 0.273630
   & 0.378073 & 0.639330\\
 \text{MC}/\text{closed},\ 2000\ \text{draws} & 0.9920 & 0.9993
   & 0.9918 & 1.0044 & 1.0157 & 1.0320
\end{array}
$$

— worst deviation $0.0320$, at the depth whose cell holds three
elements, and the six ratios straddle $1$ three and three, which is
what an exact formula gives.


#### Remark (six ratios are one observation) {#rem:mcratios}
<!-- evidence: lab_cellmom_montecarlo.py -->

Version 3 reported this check at $60$ draws with "the ratios running
$0.88$ to $0.98$ against a Monte-Carlo precision of $\pm0.09$". The
precision is right — $1/\sqrt{2\cdot59}=0.092057$ — and the band is an
ordinary draw: $60$ draws here give $1.0192$ to $1.0647$, on the other
side of $1$ and equally unremarkable, both inside $1\pm0.276$.

What is worth saying is why neither run is evidence of a bias. All six
quoted ratios fell below $1$ and all six of mine fall above, and it is
tempting to read either as $2^{-6}$. It is not: the six depths are
estimated from the *same* draws, so their deviations are strongly
correlated and their common sign is one observation, not six. That is
why the check was repeated at $2000$ draws rather than re-run at
$60$ — the audit's rule N4 asked for the ratios to straddle $1$, and
at $2000$ they do.

**All three terms are needed.** The substitution
$u_c(v)\approx n_c/\sqrt V$, which gives the right scale for
$Q_{cc}/n_c^2$, is not specific to $c$; applied to all three terms it
returns the same value $S$ for each and hence $S-2S+S=0$. So it cannot
be used to estimate the variance of the difference. Exactly, the
variance is smaller than $Q_{cc}/n_c^2$ by a factor that is a
structural constant rather than noise: measured at the top octave of
$1.6\cdot10^7$, $\mathrm{Var}/(Q_{cc}/n_c^2)$ is
$0.113119,\,0.016302,\,0.288567,\,0.551516$ at depths $0,1,3,5$.

At a quarter of that $N$ the shallow ratios are the same to six
digits — $0.113118,\,0.016303,\,0.288571$ at depths $0,1,3$ — which
is what "structural constant rather than noise" means and is stronger
than stability to three digits. Depth $5$ is not: it reads $0.557654$
there against $0.551516$ at the top, differing in the third digit.
Version 3 said the ratios reproduce to three digits without excepting
it; that was the audit's rule C2 and it fails on that one cell, which
is also much the smallest — $266$ values of $N$ in the top octave
against $1{,}687{,}911$ at depth $1$.


### The floor's uncertainty does not fall like a count


#### Proposition (coherent sums) {#prop:coh}
<!-- evidence: lab_cell_floor.py -->

The error bar of Lemma [lem:cellmom] does *not* decay like
$n_c^{-1/2}$. For fixed $v$, about $n_c/\log N$ of the terms of
$u_c(v)$ are nonzero, each of size $\log N/\sqrt V$, so
$u_c(v)\approx n_c/\sqrt V$ and

$$
\frac{Q_{cc}}{n_c^2}
  \;\approx\; \frac{\sum_v\mu^2(v)}{V}
  \;\asymp\; \frac{1}{\log N},
$$

independently of $n_c$; and by the previous subsection the variance
itself is a fixed fraction of this. So the standard error falls like
$(\log N)^{-1/2}$.


The heuristic is quantitatively loose — at $N=1.6\cdot10^7$,
$Q_{cc}/n_c^2$ reads

$$
\begin{array}{r|cccccc}
 \text{depth} & 0 & 1 & 2 & 3 & 4 & 5\\\hline
 Q_{cc}/n_c^2 & 0.123577 & 0.121208 & 0.146026 & 0.184159 & 0.235171
   & 0.307074
\end{array}
$$

against the heuristic's $\bigl(\sum_{v<N}\mu^2(v)\bigr)/V(N) =
0.049540$ — so it is offered for the *form* and not the
constant. Version 3 said the row "runs $0.124$ to $0.307$"; the
audit's rule C3 fails on the low end, because the row is not monotone
in the depth and its minimum is $0.121208$ at depth $1$, not
$0.123577$ at depth $0$. C3 also fails on the heuristic, which is
$0.049540$ and so rounds to $0.050$ rather than to the printed
$0.049$.

The form is confirmed directly. Fitting $\mathrm{se}\propto N^{-b}$ to
the exact floor across eight octaves gives
$b = 0.039451,\,0.039671,\,0.039388$ at depths $2,1,0$ —
reproducing the printed $0.0395,\,0.0397,\,0.0394$ — against the
$0.5$ that a count would give. The comparison constant is a different
matter: version 3 put the apparent exponent at
$1/(2\langle\log N\rangle) = 0.0358$, but $\langle\log N\rangle$ is
nowhere defined, and the natural reading — the mean of $\log N$ over
the eight octave midpoints — gives $0.036038$. That was the rest of
rule C4. Nothing turns on the third decimal here, since the point is
that $0.039$ is near $0.036$ and nowhere near $0.5$; the constant is
restated rather than corrected because there is nothing to correct it
against.


#### Remark (a count is not an error bar)

The consequence is not confined to this paper. Over a factor $140$ in
$N$ — the sub-range of the field $10^5<N\le1.6\cdot10^7$ on which the
exact floor was fitted, and not that field itself, whose factor is
$160$ — $n_c^{-1/2}$ says the error bar shrinks by $11.8$; it shrinks
by $1.21$, against the $1.19$ that $(\log N)^{-1/2}$ predicts over the
same factor. **An interval built from a count is
therefore about ten times too narrow at the top of this range relative
to the bottom**, and in absolute terms the discrepancy is larger still:
at the top octave the exact floor exceeds
$\mathrm{sd}(Z)/\sqrt{n_c}$ by factors of $5.8$ to $160$, growing with
the cell, exactly as $Q_{cc}/n_c^2 \asymp 1/\log N$ predicts. The cause
is that $u_c(v)$ is a coherent sum, not a self-averaging one, and the
same is true of any cell mean of a field whose summands share a common
arithmetic input.


### The mask exists


Against the exactly computed floor of Lemma [lem:cellmom], the cell
means are not zero. Over every octave from $6.25\cdot10^4$ to
$1.6\cdot10^7$, $\max_c|z_c|$ runs $9.1$ to $13.0$ and clears
Bonferroni in each. The signal is carried by the deep cells: at the top
octave the depths $3,4,5$ sit at $z = -1.6,\,-4.5,\,-9.1$ while depths
$0,1,2$ sit at $+0.0,\,+0.7,\,-0.1$.


#### Proposition (the mask survives its placebo) {#prop:placebo}
<!-- evidence: lab_mask_placebo.py -->

The permutation of Lemma [lem:placebo] is the control this claim
rests on, and until now it had been cited rather than run. Run on the
octave $(2\cdot10^6,\,4\cdot10^6]$, with cell sizes preserved exactly
and the floor recomputed from scratch for each permutation:

$$
\begin{array}{r|cccccc}
 \text{depth} & 0 & 1 & 2 & 3 & 4 & 5\\
 n_c & 383617 & 421978 & 163568 & 28507 & 2263 & 67\\\hline
 z_c\ \text{(true)} & +0.1069 & +1.1133 & -0.2314 & -2.3990
   & -5.9997 & -11.0258
\end{array}
$$

against $\max_c|z_c|$ over ten label permutations reading
$0.9359,\,1.6073,\,1.9921,\,1.5687,\,2.1392,\,1.1216,\,1.7178,\,
1.1348,\,3.2040,\,1.3192$ — mean $1.6741$, never above $3.3$. The
correlation between depth and $z_c$ is $-0.9106$ for the true
labelling and $+0.0042$ on average under permutation. So what is
detected is the correspondence between cells and divisibility, not the
cell sizes.


#### Remark (the floor is itself signal) {#rem:floorsignal}
<!-- evidence: lab_mask_placebo.py -->

The placebo also corrects a reading of Lemma [lem:cellmom]. The
audit pre-registered, as its rule L3, that the floor $\mathrm{se}_c$
would be a property of the cell sizes and so move by less than $10\%$
under permutation. It fails, and badly: the floor **collapses**, by
factors from $3.8$ at depth $5$ to $105$ at depth $0$ —
$1.2400\cdot10^{-1}$ against $1.1794\cdot10^{-3}$ there.

The mechanism is the three-term structure. For a random subset of the
band, $u_c(v)\approx (n_c/n)\,u_a(v)$, and then
$Q_{cc}/n_c^2-2Q_{ca}/(n_cn)+Q_{aa}/n^2$ cancels to leading order,
leaving only the fluctuation around proportionality. An arithmetic
cell breaks that proportionality, and the floor is what is left over.
**The exact floor is therefore not a noise level but a second
measurement of the same correspondence**, and quoting $z_c$ against it
is conservative by about two orders of magnitude: with the permuted
floor in the denominator the depth-$5$ cell would read $z\approx-42$
rather than $-11$. The figure kept is the conservative one.

What predicts the floor's *size* is the excess of the shift's
singular series over same-cell pairs. Write $\SS_2(h)$ for the
Hardy–Littlewood singular series of the shift $h$ — the local
density of pairs $(p,p+h)$ — and, for a cell $c$ in the band, let
$E_{\mathrm{same},c}[\SS_2]$ be the mean of $\SS_2(N-N')$ over pairs
$N,N'$ both lying in $c$, and $E_{\mathrm{all}}[\SS_2]$ the same mean
over all pairs in the band. The floor's magnitude tracks
$D_c := E_{\mathrm{same},c}[\SS_2]-E_{\mathrm{all}}[\SS_2]$, and it
tracks it including its shape:

$$
\begin{array}{r|cccccc}
 \text{depth} & 0 & 1 & 2 & 3 & 4 & 5\\\hline
 D_c & 0.302138 & 0.046766 & 0.412504 & 1.052129 & 1.946109
   & 3.171298\\
 \mathrm{se}_c & 0.12400 & 0.04662 & 0.14834 & 0.24178 & 0.33253
   & 0.43685
\end{array}
$$

with correlation $0.9805$ across depths at the octave
$(2\cdot10^6,\,4\cdot10^6]$. Neither row is monotone: both dip at
depth $1$, and that is the mechanism made visible. Depth $0$ is the
cell on which *none* of $3,5,7,11,13$ divides $N$, and excluding a
class concentrates the shift just as fixing one does — if
$3\nmid N$ and $3\nmid N'$ then both lie in $\{1,2\}\bmod 3$ and
$3\mid h$ with probability $\tfrac12$ against $\tfrac13$ for a random
pair. Both ends of the depth index concentrate $h$; depth $1$ mixes
the patterns and washes out. **Depth is a coarse index for $D_c$, and
the floor follows the same coarse index**, which is why the
correlation is what it is. An earlier reading of this paragraph, that
$D_c$ increases with depth, was the audit's rule M1 and is wrong.


#### Proposition (the size mechanism is scale-invariant) {#prop:scaleinv}
<!-- evidence: lab_cell_singular.py -->

$D_c := E_{\mathrm{same},c}[\SS_2] - E_{\mathrm{all}}[\SS_2]$ depends
only on the distribution of the shift $h$ in the residue classes mod
$3,5,7,11,13$ that the cell fixes. That distribution is the same at
every scale, so $D_c$ is scale-invariant and predicts an exponent of
zero at every depth.

Measured over $(10^6,2\cdot10^6]$, $(2\cdot10^6,4\cdot10^6]$ and
$(4\cdot10^6,8\cdot10^6]$, the fitted exponents are

$$
\begin{array}{r|cccccc}
 \text{depth} & 0 & 1 & 2 & 3 & 4 & 5\\\hline
 e & -0.000879 & -0.016226 & +0.000936 & +0.002904 & +0.000241
   & -0.008428
\end{array}
$$

— zero to three decimals at every depth but one. The exception is
depth $1$, and it is the depth at which $D_c=0.046766$ is six times
smaller than any other while the sampling error of the pair average is
the same $0.0013$ throughout; its fitted exponent flips sign from
$+0.065546$ to $-0.016226$ when the sample is taken ten times larger,
which is the signature of noise rather than of a trend. Depth $1$ is
also the cell that carries no mask ($z=+1.11$), so nothing rests on
it.

The audit pre-registered a $2\%$ band on the spread across octaves and
$|e|<0.01$, as rules M2 and M4, and both fail at depths $0$ and $1$ at
the original sample size. Depth $0$ passes once resampled — spread
$0.0043$, exponent $-0.000879$ — so its failure was sampling; depth
$1$ does not resolve. Two further caveats belong with the verdict: the
pre-registration allowed three sampling standard errors and the code
applied the band with no allowance, so the code was the stricter of
the two, and the verdict quoted is the code's; and the resampling is
post hoc.


**How the mask decays is open, and this version does not answer
it.** The amplitudes at depths $3,4,5$ fall with $N$; at depths $0,1,2$
they are within the exact floor at every octave measured, so no decay
exponent can be fitted there at all. Whether the exponents differ by
depth is therefore not decided by the range accessible here. What
Proposition [prop:scaleinv] does settle is that whatever produces
the decay, it is not the arithmetic excess that produces the size.

None of this threatens $C(N)=o(N)$: the mask is lower order under every
parameterisation considered.


## The negative map {#sec:closures}


Every entry below was pre-registered: the decision rule was fixed in
writing before the computation ran. Route adjudications were done in
fresh context against the source papers' verbatim lemma hypotheses. The
verdicts are stated here without their supporting statistics; those
live in the repository, where each carries its own null and its own
error bar. Where a design's threshold was chosen as an effect size
rather than in standard errors of its own null, the repository records
that too, and the verdict rests on the measurement rather than on the
criterion.


### Adjudication of existing machinery (5)


```latex
\begin{longtable}{p{0.24\textwidth}p{0.10\textwidth}p{0.56\textwidth}}
\hline
Route & Verdict & Blocking coordinate\\
\hline
\endhead
MRT \cite{MRT15} / Lichtman \cite{Li20}, shift $\to$ dilate & Blocked &
 the orthogonality factorization needs a \emph{linear} pair
 constraint ($h = m-n$); the dilate constraint
 $m'u - mu' = N(m'-m)$ is bilinear. The $h$-average is a translation;
 the $k$-average is a dilation, with no diagonalizing character
 family.\\
Tao 2016 \cite{Tao} entropy decrement, $k$-averaged & Blocked $\times 3$ &
 (i) no $k$-analog of approximate affine invariance; (ii) the sampling
 prime enters the phase multiplicatively, so the sparsification fails
 and the surviving input is the target itself; (iii) the saving is
 triple-log, below spec.\\
Lichtman rerun in dilate coordinates \cite{Li23} & Blocked &
 the congruence coupling is a genuine isomorphism but powers only the
 typical-set restriction; the decoupling blocks at the same bilinear
 coordinate, and the uniformity slot would need an $EH_\mu$-grade
 input (circular).\\
Dirichlet-polynomial fourth moment / Perron & Blocked &
 every log-saving mean-value pillar assumes coefficient
 multiplicativity in its own variable; $\mu(N-u)$ has none. Unfolding
 by Perron costs $T \gtrsim N^{1-o(1)}$; opening the fourth moment
 reproduces the binary correlation.\\
Partial slices & Partial &
 the type-I slice is already consumed; the \emph{$N$-averaged}
 theorem is provable but lands in exceptional-set territory, which
 Huang--Li cannot consume; no nontrivial fixed-$N$ slice exists.\\
\hline
\end{longtable}
```


**The common obstruction.** $\mu(m)\mu(N-mk)$ couples its
variables simultaneously through the product $mk$ and the difference
$N-mk$; the pair constraint is bilinear and is diagonalized by no
single character family, additive or multiplicative, while each
route's decisive lemma consumes exactly that diagonalization as a
hypothesis.


### Technique designs, kill-tested (9)


```latex
\begin{longtable}{p{0.05\textwidth}p{0.30\textwidth}p{0.55\textwidth}}
\hline
\# & Design & Result (pre-registered rule applied)\\
\hline
\endhead
K1 & multiplicative Fej\'er kernel on the exact dilation ladder &
 \textbf{open}. On the type-II field the design is not computable: the
 $\sqrt N$ cut empties all but thirteen of the sixty-three orbit
 columns, so what was measured is a fifth of the orbit the design
 names. On the untruncated field, where the whole orbit is live, the
 statistic lands in the band the pre-registration reserved for
 ``repeat at a second $N$ before deciding''.\\
K2 & manufactured pair congruence (determinant/Kloosterman) &
 \textbf{dead}: congruent pairs are statistically $h$-blind, and the
 verdict is quantitative --- the design would need several orders of
 magnitude more than its own detection floor.\\
K3 & Wishart / operator moment method &
 \textbf{dead}: the second, third and fourth traces sit on the Wishart
 null. No sub-Wishart surplus to fund the exchange.\\
K4 & $N$-average descent &
 \textbf{dead}: the dual field is $\delta$-blind and
 $m$-uncorrelated.\\
R1 & zero-spectrum visibility (explicit formula) &
 \textbf{dead} on the measurement, which sits below its own null. The
 quoted precision of that null is not supportable --- its spread was
 estimated from six draws --- so no exclusion interval is claimed.\\
R2 & determinant / Kloosterman phase &
 \textbf{dead} on the measurement, which sits below its control. Its
 two criteria are uncalibrated and neither decides; and the
 coherent-gain arm is not blind but sits about $2.8$ standard errors
 above its own null, under a bar that was set at $5.4$.\\
R3 & character transform of the $k$-average &
 \textbf{dead by analysis}: for $K \le N^{1/3}$ the transform deposits
 every component into thin progressions (moduli $\ge N^{2/3}$);
 Parseval forbids a statistical gain.\\
R4 & divisor switch &
 \textbf{dead}: see \S\ref{sec:R4}.\\
R5 & the circle method applied directly to $C(N)$ &
 \textbf{dead, zero margin}: both bills lie at or above the trivial
 bound, and the cap on the pointwise route is Parseval
 (Proposition~\ref{prop:E}).\\
\hline
\end{longtable}
```


#### Remark (what a null verdict costs)

A kill-test that does not fire is evidence of absence only against a
stated floor, and a threshold chosen as an effect size is not a
threshold. Three of the rows above have been restated on that basis:
K1 is reopened because its design was not computable on the field it
was run on; R1's and R2's verdicts are retained but their quoted
precisions are withdrawn. In each case the DEAD direction, where it
stands, stands on the measurement — never systematically below the
null — and not on the count of flagged levels.


### Representation classes (3, plus one open)


```latex
\begin{longtable}{p{0.10\textwidth}p{0.25\textwidth}p{0.55\textwidth}}
\hline
Class & Test & Result\\
\hline
\endhead
C-I abelian & rational-peak energy vs mask-null &
 \textbf{closed}: the abelian spectrum is mask-exact.\\
C-II inverse domain & special-frequency excess in the
 modular-inverse re-indexing &
 \textbf{closed on verification}: fired under eight-draw nulls,
 collapsed under sixty-four-draw nulls, permutation null concurring,
 second-$N$ replication empty.\\
C-IV manufactured modularity & approximate Fricke law for
 $\Phi(z) = \sum t_m e(mz)$ &
 \textbf{closed}. A true modular form drives the defect to
 $\approx 0$, so absence is meaningful.\\
C-III Motohashi type & (no finite test) &
 \textbf{open}, needing exactly the three items of \S\ref{sec:c3}.\\
\hline
\end{longtable}
```


### R4 in detail: the divisor switch does not
localize {#sec:R4}


Applied to the dilate field over its *full* ranges the switch
gives an exact identity

$$
\sum_{k\ge1}\ \sum_{m:\,mk\le N-1}\mu(m)\mu(N-mk)
  = \sum_{u<N}\mu(N-u)\sum_{m\mid u}\mu(m)
  = \mu(N-1),
$$

verified by brute force at six values of $N$. This is perfect
cancellation: $O(1)$ for a double sum of $\sim N\log N$ terms, far
beyond square root, and it is the strongest cancellation found
anywhere in this work.

E1 imposes two restrictions that make the inner divisor sum
incomplete: the type-II cut $m > \sqrt N$ and the dyadic band
$k \sim K$. The kill-test asked whether any of the cancellation
survives, by measuring block sums $S_B(j) = \sum_{k\in\text{block}}D(k)$
against the $B$-independence that Conjecture [conj:L] predicts.

*The statistic, stated.* Two normalisations of a block sum both
equal $1$ under independence, and they disagree; we use the
unweighted mean of per-block ratios,

$$
r(B) \;=\; \Bigl\langle\,
    S_B(j)^2 \Big/ \sum_{k\in\text{block }j}\mathrm{supp}(k)
  \,\Bigr\rangle_j ,
$$

whose $B=1$ baseline is $1$ under exact square-root cancellation on the
surviving support. (The ratio-of-sums normalisation gives a $B=1$
baseline of $1.783927$, because the five largest moduli carry
$0.3285$ of $\sum_k\mathrm{supp}(k)$; it is reported in the repository
and is not used here.)

**A block is a range of $k$ of width $B$**, with the dead moduli
dropped from the sums inside it — not $B$ consecutive *surviving*
moduli. The distinction is not cosmetic and version 3 left it to be
inferred: blocking by consecutive survivors gives
$1.024067,\,1.034667,\,1.036551,\,1.034833$ where blocking by $k$-range
gives $1.024067,\,1.016910,\,1.015973,\,1.025319$. Both are flat, which
is the conclusion, but they are flat at different levels and only the
second is the printed row. That was the audit's rule D3, which fails
under the first reading.

*The dead moduli.* Over the entire band on which the type-II field
is non-empty — $k<\sqrt N$, since $m>\sqrt N$ forces it — a
substantial part of the band carries no field at all. At
$N=10^8=2^8\cdot5^8$ one has $\mathrm{supp}(k)=0$, hence $D(k)=0$
identically, exactly when $4\mid k$ or $25\mid k$: $2799$ of the $9999$
moduli, or $28.0\%$ (Remark [rem:supp]). Both statistics below are
therefore $0/0$ on more than a quarter of their domain, and we state
the convention rather than leave it to be inferred: **moduli with
$\mathrm{supp}(k)=0$ are excluded**, from the blocks and from the
autocorrelation alike.

*Result: none of the cancellation survives.* $r(B)$ reads
$1.024067,\,1.016910,\,1.015973,\,1.025319$ at $B=1,2,4,8$: flat. It
degrades only at block sizes where its own sampling spread swamps it.
The sharper diagnostic, the lag-1 autocorrelation of
$D(k)/\sqrt{\mathrm{supp}(k)}$ over the surviving moduli, reads
$+0.005462$ at $N=10^8$ against a $400$-draw permutation null of
standard deviation $0.011515$ — $+0.47$ standard errors, dead zero.

We claim nothing from its sign. Under the other convention, keeping the
dead moduli with $D(k)/\sqrt{\mathrm{supp}(k)}$ set to $0$, the same
statistic reads $-0.010763$ against its own null of $0.010606$, or
$-1.01$ standard errors: still dead, but of the opposite sign. Neither
reading is significant and that is the whole of what this diagnostic
establishes. (Note also that under the stated convention consecutive
surviving moduli are $1$, $2$ or $3$ apart — the maximum gap is
exactly $3$ — so "lag-1" names position in the surviving sequence and
not a fixed gap in $k$.)


#### Remark (a null that is not reproduced) {#rem:r4null}
<!-- evidence: audit_r4_blocks.py -->

Version 3 put the second null at $0.0094$ and the second reading at
$-1.2$ standard errors. That was the audit's rule D6 and it fails: the
point estimate $-0.010763$ reproduces $-0.0108$ exactly, but no reading
of "permutation null" reaches $0.0094$. Permuting the whole
$9999$-vector gives $0.010606$; permuting only the surviving entries
with the zeros held in place gives $0.010246$. Both put the statistic
at about one standard error rather than $1.2$, which does not change
the verdict — it is dead either way — but the quoted precision is
withdrawn. What a permutation null permutes has to be said when part of
the vector is structurally zero, and $28.0\%$ of this one is.

*The mirror.* Switching the banded $L^1$ sum gives
$\sum_{k\sim K}D_{\text{full}}(k)
 = \sum_{u<N}\mu(N-u)\sum_{m\mid u,\ u/2K<m\le u/K}\mu(m)$, where
$m \asymp N/K \ge N^{2/3}$: the surviving Möbius sits on the
*long* variable, the exact opposite of the assignment that makes
Theorem [thm:A] work. The switch is not a technique that happens
to fail on the supply side; its single requirement — Möbius on the
short variable — is precisely what the type-II cut forbids.


### C-III: what it needs {#sec:c3}


C-III is the one representation class with no finite test, and it is
open. Three requirements are identified; two are settled here against
the natural construction, and the third is external.


**(1) A legitimate transform — closed for changes of
variable..** Put $A = N-a$, $B = N-b$, so $A = mk$ and $B = m'k$; then
$m'A - mB = 0$ exactly. In centered coordinates the dilate family is
therefore the *pencil of lines through the origin*, while the
shift family $B = A - h$ is a family of *parallel* lines. A
pencil is characterised by its vertex, here a finite point; a parallel
family is a pencil whose vertex is at infinity. Affine maps send
finite points to finite points, and since $\mu$ lives on $\mathbb Z$
the available changes of variable are exactly the integral affine
ones. So no change of variable carries one family to the other. This
proves what a direct probe had only measured: the obvious re-indexing
is circular, the off-diagonal computed directly and via that
re-indexing agreeing exactly. What remains is the summation-formula
class, blocked in general by the roughness of the outer weight $\mu$,
with every named candidate in it already closed.


**(2) A classification covering the type-II region —
settled against the natural bookkeeping..** That bookkeeping needs the
Möbius-side $a$ to satisfy $a \le y^{O(1)} = M^{o(1)}$. Measuring the
absolute Heath–Brown weight of the identity of level $J$ with cut
$z = M^{1/J}$, the overwhelming majority of it sits outside any region
of the form $a \le M^{\eta}$ with $\eta$ small, and the fraction
increases with $M$: measured at $M=10^4,\,10^5,\,10^6$ the
$j\in\{6,7,8\}$ share at $J=8$ reads
$0.772590,\,0.833180,\,0.886081$.

The statistic has to be stated, because the section's own point is that
it is convention-dependent. Take the absolute weight of the $j$-th term
to be its $L^1$ mass,
$W_j = \binom{J}{j}\,\#\{(m_1..m_j,n_1..n_j):
m_i\le z\ \text{squarefree},\ \prod m_i\prod n_i\le x\}$, and the share
to be $W_j/\sum_i W_i$. Then at $M=10^6$ the top-$j$ share is
$0.848668$ at $J=3$ and the $j\in\{6,7,8\}$ share is $0.886081$ at
$J=8$.

The obstruction is parameter-independent, and that is what carries the
paragraph. The identity with cut $z$ needs $z^J \ge x$, and its $j$-th
term has $a \le z^j$, which at $j = J$ is $x$ for every admissible
$(z,J)$ — while the weight concentrates in exactly those high-$j$
terms, as the shares above show and as the trend in $M$ confirms.


#### Remark (the rounding is not a convention) {#rem:hbround}
<!-- evidence: audit_hb_weight.py -->

Version 3 wrote that "the rounding of $z$ alone moves the $J=8$ entry
by $0.017$", offering that as the reason not to quote three decimals.
Two things are wrong with it. First the size: under the statistic
stated above the three roundings of $z=M^{1/J}$ move the $J=8$ entry by
$0.004356$, a quarter of the quoted figure — the reason not to quote
three decimals is a real one but it was overstated. That was the
audit's rule E3 and it fails.

Second, and this is the substantive part: the rounding is not a free
convention at all. The identity requires $z^J\ge x$, and rounding
*down* violates it — at $J=8$ it gives $z=3,4,5$ against
$x=10^4,10^5,10^6$, so $z^J = 6561,\ 65536,\ 390625$, short at every
one; the same happens at $J=3$. Only rounding up is admissible at all
three $M$. So one of the three "conventions" is not a choice between
readings of the same object but a choice that breaks the identity, and
a sensitivity computed across all three is measuring partly that.
A repair would have to bound
$\sum_a\sum_b \alpha(a)\beta(b)\mu(N-abk)$ with $\alpha$ rough and $a$
long, where every individual piece is trivial and all the content is
cancellation across $a$ — a type-II estimate for $\mu(N-\cdot)$,
i.e. the wall. Every identity decomposing $\mu$ produces such a term;
one whose every piece had either a long free variable or a
divisor-structured rough coefficient would dispose of the parity
obstruction. **Completing the construction and breaking the wall
are the same task.**


**(3) Quantitative averaged Chowla..** Shift-averaged binary
$\mu\mu$ correlations at a *fixed* log-power saving, against a
best known of $(\log)^{1-c}$. This is a named external open problem
and is where C-III stands.

That third item sets the character of the obstruction. The
adjudication's central finding was that the **dilate** average
admits no diagonalizing character family; but the **shift**
average is precisely the home ground of the
Matomäki–Radziwiłł–Tao machinery. If a legitimate transform
converting one into the other exists, what remains is quantitative
strength inside an active area rather than the absence of any coupling
surface.


## What the closures constrain


Any transform or invariance $T$ that linearizes the pair constraint at
a cost of at most a log power must satisfy all of the following.


- **(K2)** $T$ cannot manufacture its congruence externally:
  collapse structure pays only when it arises inside an intrinsic
  average that is already present.

- **(K3)** $T$ cannot bootstrap from the field's own moments:
  every moment consumes higher $\mu$-correlations and the field
  carries no sub-Wishart surplus.

- **(K4)** $T$ cannot descend from the $N$-average: its
  linearization is load-bearing and the $k$-average holds no shadow
  of it.

- **(R4)** $T$ cannot be the divisor switch or a relative:
  the switch's cancellation is a property of *complete* divisor
  sums, and the type-II cut makes every divisor sum incomplete.


The constraint that K1 was to have supplied — that $T$ cannot act
through divisibility sub-sums — is *not* asserted here, because
K1 is reopened above. Together with the round-2 findings —
characters merely relocating the difficulty into thin progressions,
determinant phases not reaching their bar — the field offers no
coupling surface in any direction that has an existing mathematical
name. In the automorphic world, shifted convolutions
$\sum a(n)a(n+h)$ are controlled because the coefficients $a(n)$
*come with* a spectral representation; $\mu$ has none, and every
route above failed exactly where it tried to borrow one. The technique
to be created is a spectral (or equivalent structural) representation
of $\mu$-pairs themselves, not a projection onto existing spectra —
and its construction is a purely theoretical act, beyond what
kill-testing can reach.


## The margin, and where the difficulty is


The requirement constrains every $N$, so the figure to quote is the one
at the extreme, not at a typical $N$. With $\max|C|\approx
a_n\sqrt{V(N)}$ and $a_n$ the Gumbel location for the number of even
$N$ below the point,

$$
\frac{N}{\max_{N\le X}|C(N)|}\;\approx\;
  \frac{\sqrt N}{a_n\sqrt{\AAA\log N}},
$$

which is $10^{4.466}$ at $N=10^{12}$ and $10^{22.842}$ at $N=10^{50}$
with $a_n=\sqrt{2\log(N/2)}$ and $\AAA=0.787275$.

Measured on the octave grid anchored at $1.6\cdot10^7$ and halving —
the grid \S[sec:floor] uses — $\max|C|/N$ falls from $0.113524$ on
$(31250,\,62500]$ to $0.010068$ on
$(8\cdot10^6,\,1.6\cdot10^7]$, so the margin at the top of the computed
range is a factor $99.325$, and extrapolating by $N^{-0.43}$ gives
$99.325\cdot 6.25^{0.43}=218.42$ at $N=10^8$. At $N=10^8$ the factor
is not extrapolated any more: Remark [rem:marginoos] measures it at
$278.4734$, which no exponent this grid supports would have given. The
requirement is not remotely tight either way.

The maximum is attained at primorials — the argmax runs
$30030,\,66990,\,139230,\,300300,\,510510,\,1021020,\,2042040,\,
4084080,\,9699690$ up the grid — which is what
Proposition [prop:V] predicts, since $\AAA(N)$ is largest exactly
where $N$ has the most prime factors and $\max|C|\approx
a_n\sqrt{V(N)}$. Above the grid it is attained at near-primorials
instead, for the same reason and with the same consequence: Remark
[rem:marginoos] is what happens when a fit is run across the steps of
that staircase and pushed past the last one.


#### Remark (the octave grid, and the fitted exponent) {#rem:grid}
<!-- evidence: audit_margin.py -->

Two of the figures above were the audit's rules M1 and M2 and both
fail as version 3 stated them. The top figure $0.0101$ and the factor
$99$ reproduce on any grid; the bottom figure $0.114$ reproduces only
on the grid anchored at $1.6\cdot10^7$, where it is $0.113524$ on the
octave $(31250,\,62500]$. On a grid anchored at
$3\cdot10^4$ and doubling, the same octave reads $0.169339$. Note also
that version 3 named one endpoint by the octave's bottom
($3\cdot10^4$) and the other by its top ($1.6\cdot10^7$).

The exponent is worse. Fitting $\max|C|/N\propto N^{-b}$ gives
$b=0.440619$ on the anchored grid and $b=0.477527$ on the doubling
one; the printed $0.43$ is neither. It is a range-specific number and
is quoted here only as the input to the extrapolation, which is what
the paragraph uses it for. **The claim that followed — that the
extrapolated factor $218.42$ is insensitive to the third decimal of
$b$ — is withdrawn**; Remark [rem:marginbracket] measures the
sensitivity and $218.42$ falls outside the range the data support.
That range is itself withdrawn by Remark [rem:marginoos], which
measures the octave rather than extrapolating to it: no exponent this
fit supports reaches the answer.

Placing the two standard estimates in the same units makes the shape of
the difficulty explicit. The target is $N(\log N)^{-A}$. The trivial
bound sits a factor $(\log N)^{A}$ above it; Cauchy–Schwarz gives
$N\sqrt{6(\log N-1)/\pi^2}$, a factor $(\log N)^{A+1/2}$ above it; and
the truth sits a *power of $N$ below* it. **The entire
difficulty is a power of $\log N$**, and it is not the size but the
proof that is missing — by Theorems [thm:D] and [thm:Dprime]
the one classical route to proving it is closed over its whole design
space.


#### Remark (the base checked without the transform, and the bracket) {#rem:marginbracket}
<!-- evidence: audit_margin_bracket.py -->

Two things were owed. The octave maxima come from a single
length-$2^{26}$ real FFT, which is exactly where a silent error would
live and where rerunning the same transform would not find one; and
the extrapolation to $10^8$ carries no bracket.

The base survives an independent route. Summing
$C(N)=\sum_{n<N}\Lambda(n)\mu(N-n)$ directly, term by term, at each
octave's argmax — one pass over the sieve, sharing no arithmetic with
the transform — reproduces every published $\max|C|/N$ to between
$1.6\cdot10^{-7}$ and $3.4\cdot10^{-7}$, which is the rounding of the
sixth printed decimal. **The two routes agree to everything that was
printed.**

The bracket does not vindicate the quoted figure. Refitting over both
grid anchors and every leave-one-out subset gives $b$ from $0.433213$
to $0.494841$ and so

$$
\text{the factor at } N=10^8 \in [219.71,\ 245.98],
$$

a span of $1.1196$. The paper's "near $220$" is inside it — barely, at
the very bottom edge — but **the $218.42$ that
audit\_margin.py prints is outside**, because it is built on the
published $b=0.43$ that the same audit refutes as M2. The right
reading is a factor between $220$ and $246$, not a number.

What makes this bracket narrow, where the two others this repository
has computed span two-thirds of a decade and nine and a half, is
neither the fit's quality nor its pre-registration record — by both
measures this is the worst of the three. The extrapolation multiplies
by $(\text{reach})^b$, so an absolute error $d$ in the exponent costs
exactly $(\text{reach})^d$: $6.25^{0.061628}=1.1196$, matching the
measured span to four decimals. **Reach sets the bracket, and this one
reaches a factor $6.25$ where the others reached millions.**

This bracket is also of the safest of the three kinds the repository
now carries. A bracket made by *assuming* a wobble is worth nothing if
the constant drifts by more than the wobble, so the two must be
compared. Here they are the same number by construction: the bracket
is the measured spread of $b$ over both grid anchors and every
leave-one-out subset, $0.433213$ to $0.494841$, a relative $0.1339$,
and nothing is assumed. Remark [rem:forecastbracket] is the second
kind — an assumed $\pm10\%$ against a measured drift of $0.0295$, so
the assumption is three times too generous and therefore safe. Remark
[rem:residueconstant] is the third: measured drift $0.1436$ against
any $\pm10\%$ one might assume, so **no bracket was published and the
forecast was refused.**

**That grading is withdrawn.** Remark [rem:marginoos] measures the
octave this bracket forecasts and finds it outside.


#### Remark (the bracket, tested on the octave it forecasts) {#rem:marginoos}
<!-- evidence: audit_margin_oos.py -->

The paragraph above graded this the safest of the three brackets
because it assumes no wobble. That grade was about how the bracket was
built, and it is correct about that. It was never about whether the
fitted shape continues past the data, and the shape is the whole of
what an extrapolation rests on.

The forecast is testable. $1/\max|C|/N$ on the octave
$(5\cdot10^7,\,10^8]$ is a quantity and not a projection, and a
length-$2^{28}$ transform reaches it. Two controls come first: the new
transform reproduces the published maximum on
$(8\cdot10^6,\,1.6\cdot10^7]$ at the published argmax to
$3.91\cdot10^{-8}$ (W1), and direct summation
$C(N)=\sum_{n<N}\Lambda(n)\mu(N-n)$ at each new argmax reproduces the
transform to $4.34\cdot10^{-18}$, a few units in the last place (W2).

**W3 is refuted.** The measured factor at $N=10^8$ is $278.4734$
against the published $[219.71,\,245.98]$ — outside it, and $1.1321$
times its top. No exponent the grids support reaches it: landing
$278.4734$ takes $b=0.562554$ where the two grids gave
$[0.433213,\,0.494841]$.

**W4 is refuted as well.** The three new octave-to-octave exponents are
$0.568360$, $0.616593$ and $0.302358$, none inside the published
spread and not on one side of it either. Locally there is no power law
to have a bracket about.

Why it had to fail is already stated two remarks above. That paragraph
records that the maximum is attained at primorials, which is what
Proposition [prop:V] predicts, since $\AAA(N)$ is largest where $N$ has
the most prime factors. A smooth $N^{-b}$ was then fitted across the
steps of that staircase and pushed a factor $6.25$ past the last one.
The published grids contain three octave steps on which $\max|C|/N$
*rises* rather than falls — $b_{\mathrm{local}}=-0.028471$ where
$30030$ first becomes available and $-0.005500$ where $510510$ does —
so the fitted $b$ averages over jumps, and the stretch it was
extrapolated across contains none: the last primorial below $10^8$ is
$9699690$, which entered at the top of the fitted range, and the next
one is past $1.28\cdot10^8$. What the maximum does instead is move onto
near-primorials that trade one prime for a larger one,
$16546530=2\cdot3\cdot5\cdot7\cdot11\cdot13\cdot19\cdot29$ and
$70450380=2^2\cdot3^2\cdot5\cdot7\cdot11\cdot13\cdot17\cdot23$, and it
falls faster than every fitted $b$ while doing so — $0.495770$ over the
whole new stretch against a largest fitted $0.494841$. That was W5,
registered before the run, and it holds.

One caveat on the top octave. $(6.4\cdot10^7,\,1.28\cdot10^8]$ and
$(5\cdot10^7,\,10^8]$ both peak at $N=70450380$, so the last local
exponent is not an independent measurement of the decay — it records
where the maximum stopped moving. The whole-stretch $0.495770$ is
therefore a lower bound on the true decay, which is the direction that
makes W5's reading safe rather than the one that makes it lucky.

**What is left is smaller and is enough.** A bracket whose wobble
equals its drift is honest about its own construction and says nothing
about its reach: the drift of $b$ is measured *inside* the fitted range
and cannot see a step that lies outside it. The extrapolation is
withdrawn and the measurement replaces it — at $N=10^8$ the margin is a
factor $278.4734$, and the requirement $C(N)=o(N)$ is not remotely
tight. **The general lesson is the one this section can least afford to
forget: a self-consistent bracket is not a tested one, and all three
brackets graded above were graded on construction alone.** This is the
first of the three to be tested against a measurement on the other
side of it. The comparison to draw is with Remark
[rem:primorialreach], the only other bracket in these papers that has
been run out of sample: there the point estimate failed and the
*interval survived*, which is what an interval is for. Here the
interval failed too — and the difference is not care but shape. That
bracket extrapolated a slope along a ladder built to be smooth; this
one extrapolated across a staircase whose next step was known, from
Proposition [prop:V], to lie outside the range.


#### Remark

Everything measured here is at $N\le1.6\cdot10^7$, with two arms at
$N\approx10^8$. The no-go theorems of Section [sec:demand] are
asymptotic and acquire content only near $N\approx10^{480}$; nothing
measured here constrains that range, and nothing there is contradicted
by these measurements.


## What this version does not claim {#sec:notclaimed}


Stated so the gap is visible rather than inferred, and read against the
claim list of \S[sec:claims].


- **No law for $C(N)$ is conjectured here.** An earlier draft
  stated a conjecture of the form
  $C(N) = m(N) + \sqrt{V(N)}\,G(N)$ with $G$ Gaussian, supported by
  four measurements of the bulk, the tail, and the phase content in
  $\log N$. Re-verification found that the bulk and tail figures are
  not reproducible under the cell index and pooling the text specifies,
  that the phase-content measurement is reproduced by a coin and so by
  Lemma [lem:coin] is not measuring $\mu$, and that the mask's
  quoted significance was overstated. The conjecture may well be true;
  the evidence offered for it was not sound, and it is withheld until
  the statistics it rests on are defined and re-measured.

- **No level or rate for $\rho$ is quoted.** By
  Lemma [lem:coin] the centred estimator cannot distinguish $\mu$
  from a coin, so any such figure calibrates nothing.

- **No decay exponent for the mask is quoted.** At the three
  shallowest depths the amplitude does not clear the exact floor at any
  scale measured, so nothing can be fitted there.

- Conjecture [conj:L]'s non-E1 stamps are single-witness and are
  cited as motivation, not as verified measurement.


Each of these is an open question rather than a retraction of a
mechanism, and each is carried forward in the repository's continuation
notes.


## Methodology


Four rules were followed throughout, and they are the reason we believe
the negative results.


- **Pre-registration.** Every design's decision rule was
  written before its computation ran, including the rule that would
  refute the hypothesis under test.

- **Adversarial review in fresh context**, against the source
  papers' verbatim lemma hypotheses rather than against summaries of
  them, and repeated: three independent reviews have been run, each
  from the statements alone and without access to its predecessors'
  findings or code. Their coverage is uneven and we state that rather
  than round it up. Material inherited from Version 1 has been through
  all three. Statements first written for Version 2 — among them
  Proposition [prop:W] as stated, Lemma [lem:MP]'s truncated
  form, Proposition [prop:coh]'s derivation,
  Lemma [lem:placebo]'s statement and \S[sec:R4]'s block-ratio
  definition — had one reviewer when Version 2 was written and have
  had one since. The measured detection rate of a single review of this
  kind is about one half, so nothing here should be read as certified.

- **Nulls before thresholds.** A threshold means nothing until
  the spread of the quantity it judges has been measured, and a null
  must preserve the structure of the field it is a null for. A null
  built by shuffling a correlated field is not a null for a statistic
  that reads its correlation.

- **Weights and fields before comparisons.** Two summaries of
  one object are not comparable until each one's weight is stated
  (Remark [rem:scale]), and a measured figure is not interpretable
  until the range it was measured on is stated. Both rules were adopted after the
  corresponding mistakes were made and caught.


Three traps are worth recording because they are traps and not slips.


- *The $(q,N)>1$ main-term trap.* In the proof of
  Theorem [thm:A], assigning a main term $T_m/\varphi(q)$ to a
  residue class with $(q,N)>1$ — which is degenerate, its true
  contribution being $O(\log^2 N)$ — shifts the density by exactly
  $N/\varphi(N)$. Carried into the $\log k$ branch, that error
  produces an apparent *refutation of $EH_\mu$*. This is the most
  plausible false positive we met.

- *Nulls estimated from too few draws* inflate maxima across a
  family of tests; one representation-class "hit" survived
  eight-draw nulls and died under sixty-four-draw nulls plus
  replication, and one kill-test's quoted precision rests on a
  six-draw spread.

- *Local factors evaluated at the wrong modulus.* A
  reconstruction of the wall's excess divided by $W(N)$ where the
  definition calls for $V(N)$; the two differ by exactly $\AAA(N)$, the
  object Proposition [prop:V] is about.


## Reproducibility


All measurements run on a laptop with Python and numpy. The
repository's `PROVENANCE.md` maps every numbered statement in
this paper to the code that verifies it and the result file it was
read from. One-shot verification of the core corpus is
`python code/verify/verify\_all.py`. The consistency gate, which
runs every consistency checker and exits nonzero if any fails, belongs
to the program's own record rather than to this paper.


#### Remark (the two constants, and how far they were trusted) {#rem:constants}
<!-- evidence: audit_constants.py -->

Every threshold in this program is $\SS(N)(1-A(N))$ or $\SS(N)$, and
both are built from $2C_2=2\prod_{p>2}(1-1/(p-1)^2)$ and
$\prod_p(1-1/(p(p-1)))$, which nine scripts compute inline. They were
checked here by a route that enumerates no primes at all: with
$u=1/p$ each Euler factor is a power series in $u$, so the log of the
product is $\sum_{n\ge2}g_nP(n)$ with $P$ the prime zeta function,
obtained from $\zeta$ by Möbius inversion. That route returns
$1.320323631694$ and $0.373955813619$, missing the published constants
by $6.261\cdot10^{-12}$ and $1.920\cdot10^{-11}$. The sieved
production route at bound $10^7$ returns $1.320323639431$ and
$0.373955815811$, and its miss tracks $1/(P\log P)$ across every
bound — it is the truncation tail and nothing else.

Three of the four pre-registered rules fail and they fail together,
for one reason. **W1**, **W2** and **W3** asked a *truncated* Euler
product to match an *exact* constant to $10^{-9}$, when at the largest
bound used the tail is already $6.204\cdot10^{-9}$. The tolerance was
set below a floor the mathematics imposes, which is Remark
[rem:cap]'s mistake again: an effect size in place of a null. Judged
against the floor, both routes are correct and the prime-free one is
three orders inside it.

What the run found instead is worth the cycle. **W4** asked whether
the bound reaches the printed precision, and the answer is that it
does. Recomputing the Goldbach threshold at $N=2^a5^b$ from each
bound gives $0.374486$ at $10^5$ and $0.374487$ from $2\cdot10^5$
upward: the last digit this program prints depends on how far the
script that printed it happened to have sieved. Seven of the nine
implementations took the product over the measurement's own prime
list. None was low enough to be wrong — the smallest was $10^6$,
twelve times inside the margin — but the quantity was not pinned, and
pinning it moved a printed number in this note's companion: $R$ at
$\theta'=0.51$, $N=8\cdot10^5$ now reads $0.170167$, one unit in the
last place from what had been recorded. A constant that shifts a
published digit according to which script computed it is not a
constant, and it is now fixed at a single bound by the program's
gate.


#### Remark (the arithmetic underneath, audited) {#rem:sieve}
<!-- evidence: audit_sieve.py -->

Every measurement in this paper and its companion is a sum over
$\Lambda$ and $\mu$ built inside the measuring script itself. The
dominant construction's $\mu$ half is not a plain sieve: an `int32`
cofactor array carries the single prime factor that may exceed
$\sqrt n$, and a sign flip at the end accounts for it. A fault there
would be invisible to every consistency check in the repository,
because every script would be wrong in the same direction and would
agree with every other script. That is the one failure mode internal
agreement cannot detect, so it is checked from outside.

Auditing it turned up something first. There is not one construction
but three: fifteen scripts use the cofactor trick, three build $\mu$
from the recurrence $\mu(v)=-\mu(v/p)$ off a smallest-prime-factor
table, and two more differ only in what they return. An audit that
pins one and reports on "the sieve" claims more than it checked, so
all three are compared, and the audit carries a manifest of the
implementations it has seen; the program's gate holds that manifest
against the repository, so a fourth variant cannot enter unaudited.

$\mu$ and $\Lambda$ were rebuilt on $1\le n\le2\cdot10^6$ by a
structurally different route — a smallest-prime-factor sieve followed
by explicit factorisation of each $n$ — and compared elementwise. The
$\mu$ arrays disagree at $0$ indices, and the worst relative
difference on $\Lambda$ is $0$: not within tolerance, identical. The
recurrence construction, compared the same way on $1\le n\le5\cdot10^5$,
also disagrees at $0$ indices with worst relative $0$ on $\Lambda$. The
defining identities were then swept directly on the production arrays
to $2\cdot10^5$: $\sum_{d\mid n}\mu(d)=[n=1]$ fails at $0$ values of
$n$, and $\sum_{d\mid n}\Lambda(d)=\log n$ holds to a worst relative
$2.911\cdot10^{-16}$. Rebuilding $B(N)/N$ from the *independent*
arrays returns $0.808567150$, $0.739493602$, $0.730264309$ at
$N=2\cdot10^5,\,4\cdot10^5,\,8\cdot10^5$, reproducing the published
table at its printed width. At the top of the range actually used,
$n\le2.56\cdot10^7$, the two global consequences a sign error would
destroy both hold: $|M(x)|/\sqrt x$ stays at $0.3240$ or below over
eight checkpoints, and $|\psi(x)-x|/x=0.000050$.

This settles nothing mathematical. It removes one way for everything
else to be wrong at once.


## Summary


- The demand side of the Huang–Li reduction is closed: one half
  is now an unconditional theorem, the other half is equivalent to
  Huang–Li's equation (22) — hence already gives binary Goldbach,
  and gives the asymptotic $\tilde r(N)\sim\SS(N)N$ exactly when
  $C(N)=o(N)$ — and the interior of the design space between them is
  empty.

- The supply side is closed to measurement: the object is
  featureless in every direction with a name, and the pre-registered
  closures say precisely which directions were checked and why each
  fails. One of them, K1, is reopened here.

- The wall has an exact second moment, an exact aggregate identity,
  an exact closed form for the fluctuation of any cell mean, and a
  deterministic mask that is detected against that floor at every scale
  measured. Its error bar falls like $(\log N)^{-1/2}$ and not like a
  count.

- The natural conditional route does not close: the coefficient
  amplifying $S(h)$ grows like $N/\log N$ and is nonnegative, and the
  true $S$ already overspends the budget that coefficient defines, so
  no bound on $|S(h)|$ gives $\mathrm{Var}\,C\sim V$; only signed
  cancellation across $h$ would, and nothing in the literature supplies
  it.

- Net progress toward Goldbach is zero. The contribution is the
  map, the exact facts, one unconditional theorem that removes a
  Goldbach-neutral demand, one conjecture that any future construction
  must reproduce, and one defect report for the source paper.


## References

- **[HL]**  Jing-Jing Huang and Huixi Li, *On the connection
between the Goldbach conjecture and the Elliott–Halberstam
conjecture*, arXiv:2005.03811v2 [math.NT], 2022.

- **[ThmA]**  *An unconditional bound for a Möbius-weighted
correlation sum in fixed residue classes*, companion note,
`paper/theorem\_A.tex` in this repository.

- **[Pan]**  Cheng-Dong Pan, *A new attempt on Goldbach
conjecture*, Chinese Ann. Math. **3** (1982), 555–560.

- **[Tao]**  Terence Tao, *The logarithmically averaged Chowla
and Elliott conjectures for two-point correlations*, Forum Math. Pi
**4** (2016).

- **[Li20]**  Jared Duker Lichtman, *Averages of the Möbius
function on shifted primes*, arXiv:2009.08969.

- **[Mir49]**  L. Mirsky, *The number of representations
of an integer as the sum of a prime and a $k$-th power free
integer*, Amer. Math. Monthly **56** (1949), 17–19.

- **[MRT15]**  Kaisa Matomäki, Maksym Radziwiłł
and Terence Tao, *An averaged form of Chowla's conjecture*,
Algebra Number Theory **9** (2015), 2167–2196.

- **[Vau88]**  R. C. Vaughan, *The $L^1$ mean of exponential
sums over primes*, Bull. London Math. Soc. **20** (1988),
121–123.

- **[Li23]**  Jared Duker Lichtman (with an appendix by Sary
Drappeau), *Primes in arithmetic progressions to large moduli, and
Goldbach beyond the square-root barrier*, arXiv:2309.08522
[math.NT], 2023.

- **[MV17]**  M. Ram Murty and Akshaa Vatwani, *Twin primes and
the parity problem*, J. Number Theory **180** (2017), 643–659.
