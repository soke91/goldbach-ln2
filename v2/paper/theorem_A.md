# An unconditional bound for a Möbius-weighted correlation sum  in fixed residue classes, with two consequences for the  Huang–Li reduction of Goldbach to Elliott–Halberstam

```latex
% 수식이 쓰는 매크로 — 렌더러/역변환용
\renewcommand{\SS}{\mathfrak{S}}
\newcommand{\Emu}{E_\mu}
\DeclareMathOperator{\rad}{rad}
\DeclareMathOperator{\lcm}{lcm}
```

## Abstract

Huang and Li [arXiv:2005.03811] show that the binary Goldbach conjecture
for large even $N$ follows from $EH$ together with a Möbius-twisted
variant $EH_\mu$ whose levels sum to more than $1$; by Bombieri–Vinogradov
their Corollary 1 reduces this to $EH_\mu(N^{\theta'})$ alone for a single
$\theta' > 1/2$. Their proof consumes $EH_\mu$ at exactly two places, the
functionals they call $E_3(\alpha)$ and $E_4(\alpha)$. We prove
unconditionally (Theorem [thm:A]) that the Möbius-weighted
correlation sum underlying $E_4$ is
$\ll_A N(\log N)^{-A}$ for every $A>0$, uniformly in the truncation
point; consequently (Corollary [cor:B]) the $E_4$ consumption of
$EH_\mu$ is unnecessary and the entire $EH_\mu$ demand of the Huang–Li
argument collapses to the single scalar $E_3(\alpha)$. We then show
(Theorem [thm:C]) that the corresponding statement for $E_3$ — the
same sum with the weight $\log k$ — is not merely hard but, by an
unconditional identity, *equivalent* to Huang–Li's own
equation (22): it therefore already yields binary Goldbach for large
even $N$, and yields the asymptotic $\tilde r(N)\sim\SS(N)N$ exactly
when $\sum_{n<N}\Lambda(n)\mu(N-n) = o(N)$. The root cause is
$\mu * \log = \Lambda$. Finally we record a defect in the published
paper, first observed by S. Zheleznov: equation (18) drops an
$n$-dependent constraint present in the definition of $S_2(\alpha)$.
The missing term is exhibited and shown harmless under the hypotheses
already assumed (Section [sec:delta]); the authors' own correction
takes a different route, on which that term does not arise
([rem:movingswitch]). The two routes differ by exactly that term, and
Theorem [thm:A] together with it gives the bound against the corrected
formulation (Proposition [prop:movingcut]).

We state plainly what is *not* claimed: Theorem [thm:A] removes
only the part of the demand that carries no Goldbach content, and the
net progress toward the Goldbach conjecture is zero. We also state
plainly what is *standard*. The ingredients of Theorem [thm:A]
are classical and its mechanism — a divisor switch that moves the work
onto a short variable, followed by Bombieri–Vinogradov — is ordinary
practice; what is offered is the application, together with the
corrections it turns out to require. Theorem [thm:D] is, in genre, a
precise form for this particular reduction of the phenomenon Bombieri's
asymptotic sieve records, and a reader who finds its conclusion
unsurprising is right to. Section [sec:lit] places both against the
record, including Lichtman's level of distribution beyond the
square-root barrier [Li23], which postdates [HL].


## Introduction


### The Huang–Li reduction


Let $N$ be a large even integer, $\Lambda$ the von Mangoldt function,
$\mu$ the Möbius function, and

$$
\tilde r(N) \;=\; \sum_{n<N} \Lambda(n)\Lambda(N-n),
  \qquad
  \SS(N) \;=\; 2\prod_{p>2}\Bigl(1-\tfrac{1}{(p-1)^2}\Bigr)
                \prod_{\substack{p\mid N\\ p>2}}\Bigl(1+\tfrac{1}{p-2}\Bigr).
$$

Write $EH_\mu(Q)$ for the assertion that for every $A>0$

$$
\begin{equation}\label{eq:EHmu}
  \sum_{q\le Q}\ \max_{y<N}\ \max_{(a,q)=1}
  \Bigl|\sum_{\substack{n\le y\\ n\equiv a\ (q)}}\Lambda(n)\mu(N-n)
        -\frac{1}{\varphi(q)}\sum_{n\le y}\Lambda(n)\mu(N-n)\Bigr|
  \;\ll_A\; \frac{N}{(\log N)^A}.
\end{equation}
$$

Huang and Li [HL] prove that $EH(N^\theta(\log N)^{2A+8})$ together
with $EH_\mu(N^{1-\theta})$ implies
$\tilde r(N)\ge \SS(N)(1-A(N))N + O(N(\log N)^{-A})$, where

$$
\begin{equation}\label{eq:AN}
  A(N)=\prod_{\substack{p\nmid N\\ p>2}}\Bigl(1-\frac{1}{p(p-1)}\Bigr),
\end{equation}
$$

and deduce (their Corollary 1) that, in view of Bombieri–Vinogradov,
binary Goldbach for all sufficiently large even integers follows from
$EH_\mu(N^{\theta'})$ alone, for any single $\theta'>1/2$.

Throughout this note we work in that *Corollary-1 regime*: a fixed
$\theta'\in(1/2,1)$ is given, we put

$$
\delta := \theta' - \tfrac12 > 0, \qquad K := N^{\theta'},
  \qquad M := N/K = N^{1-\theta'} \le N^{1/2-\delta}.
$$


### Where $EH_\mu$ is consumed


The hypothesis [eq:EHmu] enters the Huang–Li proof at exactly two
places. With $\alpha$ the Vaughan-type parameter of their §3 and
$K=(N-1)/\alpha$, these are (their notation, and note the residue class
is the *fixed* class $a\equiv N$):

$$
\begin{align}
  E_3(\alpha) &= \sum_{\substack{k<K\\ (k,N)=1}}\mu(k)\log k
    \Bigl[\sum_{\substack{n<N\\ n\equiv N\ (k)}}\Lambda(n)\mu(N-n)
      -\frac{1}{\varphi(k)}\sum_{n<N}\Lambda(n)\mu(N-n)\Bigr],
      \label{eq:E3}\\
  E_4(\alpha) &= \sum_{\substack{k<K\\ (k,N)=1}}\mu(k)
    \Bigl[\sum_{\substack{n<N\\ n\equiv N\ (k)}}\Lambda(n)\mu(N-n)\log(N-n)
      -\frac{1}{\varphi(k)}\sum_{n<N}\Lambda(n)\mu(N-n)\log(N-n)\Bigr].
      \label{eq:E4}
\end{align}
$$

In both cases Huang–Li discard the signs $\mu(k)$ by the triangle
inequality on the first line and then appeal to [eq:EHmu] (their
Lemma 4 for $E_4$). The two functionals differ in exactly one respect:
$E_3$ carries the weight $w_k=\log k$ and $E_4$ carries $w_k=1$ (the
$\log(N-n)$ inside $E_4$ is $k$-independent and is removed by partial
summation). This note shows that this one difference decides everything.


### Results


For $(k,N)=1$ and $1\le t<N$ put

$$
\begin{equation}\label{eq:Emu}
  \Emu(t;k) := \sum_{\substack{n\le t\\ n\equiv N\ (k)}}\Lambda(n)\mu(N-n)
             - \frac{1}{\varphi(k)}\sum_{n\le t}\Lambda(n)\mu(N-n),
  \qquad
  C(t) := \sum_{n\le t}\Lambda(n)\mu(N-n),
\end{equation}
$$

and, for a weight $w$,

$$
T_w(t) \;:=\; \sum_{\substack{k<K\\ (k,N)=1}} \mu(k)\,w(k)\,\Emu(t;k).
$$

Thus $E_4$ is obtained from $T_1$ and $E_3$ from $T_{\log}$.


#### Theorem {#thm:A}
<!-- evidence: analytic -->

Fix $\theta'\in(1/2,1)$ and put $K=N^{\theta'}$. Then for every $A>0$,

$$
\sup_{1\le t<N}\ \bigl|T_1(t)\bigr|
  \;=\;\sup_{1\le t<N}\
  \Bigl|\sum_{\substack{k<K\\(k,N)=1}}\mu(k)\,\Emu(t;k)\Bigr|
  \;\ll_{A,\theta'}\; \frac{N}{(\log N)^{A}} .
$$

Moreover every ingredient of the proof except Bombieri–Vinogradov
contributes only $O\!\left(N e^{-c\sqrt{\log N}}\right)$: the power of
$\log$ in the statement is imposed by Bombieri–Vinogradov alone.


#### Remark

Theorem [thm:A] is a statement about a *signed* sum in a
*fixed* residue class. It is not an instance of, and does not
imply, $EH_\mu(N^{\theta'})$, which asks for absolute values and a
maximum over residue classes. The mechanism is described in
Section [sec:mechanism]: after divisor switching one *completes*
the divisor sum, and on the complementary range $k\ge K$ the Möbius
factor on the long variable squares itself away while the surviving
Möbius sits on a variable of size $\le N^{1/2-\delta}$. The completion
is not a formality — without it the cofactor is not short and the
configuration is the one with no known machine.


#### Corollary {#cor:B}
<!-- evidence: analytic -->

In the Corollary-1 regime, $E_4(\alpha)\ll_A N(\log N)^{-A}$ for every
$A>0$, unconditionally. Consequently Huang–Li's estimate
$S_4(\alpha)=O(N(\log N)^{-A})$ requires no hypothesis, their Lemma 4 is
not needed, and the entire $EH_\mu$ demand of their argument collapses
to the single scalar $E_3(\alpha)$ (together with the term $\Delta$ of
Section [sec:delta] where the published form of (18) is used; it is
likewise unconditional here, and on the corrected form it is absent).


#### Theorem {#thm:C}
<!-- evidence: audit_E3_constant.py -->

In the Corollary-1 regime one has the unconditional identity

$$
E_3(\alpha) \;=\; \sum_{n<N}\Lambda(n)\Lambda(N-n)
     \;-\; \SS(N)\Bigl(N-\sum_{n<N}\Lambda(n)\mu(N-n)\Bigr)
     \;+\; O_A\!\left(\frac{N}{(\log N)^{A}}\right),
$$

In particular the bound $E_3(\alpha)\ll_A N(\log N)^{-A}$ — the form
of $EH_\mu$ their argument consumes, though not the weakest that
suffices; see Proposition [prop:onesided] — is *equivalent* to
Huang–Li's equation (22). It therefore yields
binary Goldbach for large even $N$ through their Theorem 1; and granted
it, the asymptotic $\tilde r(N)\sim\SS(N)N$ holds if and only if
$C(N)=o(N)$.


#### Remark (three statements, not two) {#rem:threeway}

The identity says that the difference displayed above *is*
$E_3(\alpha)$, so the hypothesis $E_3\ll_A N(\log N)^{-A}$ and
Huang–Li's (22) are the same assertion. What [HL] record
immediately after their (22) is a different equivalence, between
$\tilde r(N)\sim\SS(N)N$ and $C(N)=o(N)$. Chaining the two,
$E_3\ll_A N(\log N)^{-A}$ delivers binary Goldbach — Huang–Li's
Theorem 1 supplies $|C(N)|\le A(N)N(1+o(1))$ by the triangle
inequality, and $A(N)<1$ — but it does not on its own deliver the
asymptotic, which needs $C(N)=o(N)$ in addition. An earlier version of
this note stated the last two as one equivalence; the missing
hypothesis is exactly the quantity the companion paper's
Sections 4–5 are about.


Theorem [thm:C] is the more consequential half. It says that the
demand side is closed at the level of *identities* rather than of
estimates: no choice of $\theta'$, of truncation, or of smoothing can
evade it, because the root cause is the identity $\mu*\log=\Lambda$
from which Huang–Li start. The weight $\log k$ is the carrier of the
Goldbach content, so every divisor switch must hand it back.

Closure at the level of identities is not closure at the level of
*strength*, and the next proposition separates the two. Every
condition on $E_3$ unwinds, through Theorem [thm:C], to a condition
on $\tilde r(N)$ and $C(N)$ — that is the identity part, and it is
permanent. But Huang–Li consume $E_3$ two-sidedly and at a saving of
every power of $\log$, and Goldbach does not need either.


#### Proposition (the demand is one-sided, and its threshold is not a log power) {#prop:onesided}
<!-- evidence: lab_onesided_margin.py -->

In the Corollary-1 regime, binary Goldbach for all sufficiently large
even $N$ follows from

$$
\begin{equation}\label{eq:onesided}
  E_3(\alpha) \;>\; -\,\SS(N)\bigl(1-A(N)\bigr)N\,\bigl(1+o(1)\bigr),
\end{equation}
$$

with $A(N)$ as in [eq:AN]. This is strictly weaker than
$E_3\ll_A N(\log N)^{-A}$ in two independent ways: it constrains
$E_3$ from one side only, and its threshold is
$\SS(N)(1-A(N))N$, which is $\asymp N$ for almost all even $N$ and
never smaller than $c\,N/(\log N\log\log N)$.




**Proof.** 
By Theorem [thm:C],
$\tilde r(N)=\SS(N)(N-C(N))+E_3(\alpha)+O_A(N(\log N)^{-A})$. The
triangle inequality gives
$|C(N)|\le U(N):=\sum_{n<N}\Lambda(n)\mu^2(N-n)$, and $U(N)=A(N)N(1+o(1))$
by Mirsky's theorem on squarefree shifted primes together with partial
summation — this is the count the companion paper's
Proposition \ref{prop:V} evaluates. Hence
$\SS(N)(N-C(N))\ge \SS(N)(1-A(N))N(1+o(1))$, and $A(N)<1$ for every
$N$, so the right-hand side is positive. Substituting and using that
the $O_A$ term may be taken below any fixed power of $N/\log N$ gives
$\tilde r(N)>0$ under [eq:onesided]. For the size of the threshold,
$A(N)=\prod_{q\nmid N}(1-\tfrac1{q(q-1)})$ is largest when $N$ carries
the most small primes, so $1-A(N)$ is smallest at primorials, where
$1-A(N)\asymp\sum_{q>p}\tfrac1{q(q-1)}\asymp 1/(p\log p)$ with
$p\asymp\log N$; and $\SS(N)\gg1$.
 ∎




#### Remark (what this does and does not open) {#rem:onesided}
<!-- evidence: lab_onesided_margin.py -->

Version 3 of this note wrote that "the demand side admits no weaker
sufficient condition" and called $E_3\ll_A N(\log N)^{-A}$ "the weakest
form of $EH_\mu$ that their argument actually needs". Both are too
strong: [eq:onesided] is weaker on both counts, and the gap in the
threshold alone is a factor $(\log N)^{A}$ for every $A$.

The margin is measured rather than asserted. Over even
$N\le1.6\cdot10^7$ the quantity $\SS(N)(1-A(N))$ has median
$0.333459$ and $0.1$-percentile $0.119639$; its minimum is
$0.060890$, attained at $N=9699690=2\cdot3\cdots19$, the largest
primorial in range, and the ten smallest values are all
primorial-like. Multiplied by $\log N\log\log N$ the minimum over
$N\ge10^5$ is $2.482019$, attained at $N=510510$ — the mechanism of
the proposition, visible directly. The step that produces the margin
is exact to five decimals: $U(N)/(A(N)N)$ has mean $0.999996$ with
standard deviation $0.000170$ over the top octave. Checked against the
truth: no even $N$ in $[10^5,1.6\cdot10^7]$ has
$\tilde r(N)<\SS(N)(1-A(N))N$, and the slack between them has minimum
$3.7038$, median $4.2902$ and maximum $95.0889$ — the maximum at
$N=9699690$, so the margin is thinnest exactly where the Goldbach
count is largest.

One pre-registered rule of that audit, F3, is refuted, and by a
mis-specified field rather than by the mathematics: it asked for
$\SS(N)(1-A(N))\log N\log\log N>1$ over *all* even $N$, and
$\log\log N<0$ below $N=e^{e}=15.15$. On $N\ge16$ the minimum is
$0.810202$, at $N=30$; on $N\ge10^3$ it is $1.737799$, at $N=2310$;
on $N\ge10^5$ it is $2.482019$. Every one of those argmins is a
primorial, which is the content. The constant $c$ in the proposition
is therefore stated as a constant and not as $1$.

What it does not open is the divisor switch. Theorem [thm:D] costs
$\exp(c_1\sqrt{\tfrac12\log N})$, which exceeds $\log N\log\log N$ as
it exceeds every power of $\log N$, so the weaker threshold is still
out of reach of that design space. By the test this program applies to
any new sufficient condition — does it imply the target, and is it
*easier* than the target — the first half is proved here. The
second is what the next proposition is for.


#### Proposition (the demand needs no saving in $\log N$) {#prop:nolog}
<!-- evidence: lab_onesided_demand.py -->

Put

$$
\begin{equation}\label{eq:Bsum}
  B(N) \;:=\; \sum_{\substack{k<K\\ (k,N)=1}} (\log k)\,
  \bigl|\Emu(N;k)\bigr| ,
\end{equation}
$$

the quantity Huang–Li reach by the triangle inequality before appealing
to $EH_\mu$. Then binary Goldbach for all sufficiently large even $N$
follows from

$$
\begin{equation}\label{eq:nolog}
  B(N) \;\le\; (1-\varepsilon)\,\SS(N)\bigl(1-A(N)\bigr)N
\end{equation}
$$

for a single fixed $\varepsilon>0$. Since $\SS(N)(1-A(N))\asymp1$ for
almost all even $N$, [eq:nolog] is a *constant-factor* bound on
exactly the object $EH_\mu$ bounds — absolute values kept, and no
saving in $\log N$ whatsoever — where $EH_\mu(N^{\theta'})$ asks for
$\ll_A N(\log N)^{-A}$ for every $A$.




**Proof.** 
Immediate from Theorem [thm:C], $|E_3(\alpha)|\le B(N)$, and
$|C(N)|\le A(N)N(1+o(1))$ as in Proposition [prop:onesided]:
$\tilde r(N)\ge\SS(N)(1-A(N))N(1+o(1))-B(N)+O_A(N(\log N)^{-A})$.
 ∎




#### Remark (what the weakening relocates) {#rem:relocate}
<!-- evidence: lab_onesided_demand.py -->

Measured at $\theta'=0.56$ and $N=2\cdot10^5$ through $3.2\cdot10^6$,
$B(N)/N$ reads
$0.8086,\,0.7395,\,0.7303,\,0.6547,\,0.5916$ against a threshold of
$0.3745N$: the ratio $B(N)/\bigl(\SS(1-A)N\bigr)$ falls
$2.159,\,1.975,\,1.950,\,1.748,\,1.580$. So the object is already
within a factor $1.6$ of sufficing, and the factor is falling. For
contrast, the consumed bound needs $B$ beaten by $9$, $133$ and $1988$
at $A=1,2,3$ at the same $N$.

**The weakening therefore takes the demand entirely off the
saving axis and leaves it on the level axis.** Bombieri–Vinogradov
supplies every power of $\log$ at level $N^{1/2-\delta}$ and nothing at
$N^{\theta'}$ with $\theta'>1/2$; moving what is asked from "every
power of $\log$" to "a constant factor" does not move the level, and
the level is the whole obstruction. That is the honest answer to the
second half of the test: [eq:nolog] is weaker than $EH_\mu$ by every
power of $\log N$ and is not thereby easier, because the difficulty was
never on that axis. What has been gained is that this is now visible
rather than conjectured — and that the level axis can then be
measured, which is what Remark [rem:levelmeas] does.

Three pre-registered rules of that audit fail and each one is
informative. G1 asked $B(N)\asymp N(\log N)^2$, the size
Brun–Titchmarsh gives through $|\Emu(N;k)|\ll N/\varphi(k)$ and
$\sum_{k<K}1/\varphi(k)\ll\log K$; measured, $B(N)\asymp N$, so the
truth already carries both powers of $\log$ and the provable bound is
a gross overestimate. G2 asked whether keeping the signs $\mu(k)$ is
what makes $E_3$ small: $|E_3|/B(N)$ reads
$0.54132,\,0.51889,\,0.43431,\,0.39839,\,0.35041$, so the signs buy a
factor between $2$ and $3$ and no more — the smallness is inside
$\Emu$ itself, not in the cancellation across $k$. G3 asked for the
required saving to grow like $(\log N)^2$; it does not grow at all,
which is the proposition.

Finally, [eq:nolog] is asymptotic and the accessible range straddles
it, which rule G4 was written to detect and which it confirms:
$E_3(\alpha)/N$ reads $-0.4377,\,-0.3837,\,-0.3172,\,-0.2608,\,-0.2073$
against the threshold $-0.3745$, so [eq:onesided] *fails* at
$N=2\cdot10^5$ and $4\cdot10^5$ and holds from $8\cdot10^5$ upward. A
sufficient condition that is false at the bottom of the computable
range is not thereby suspect, but it cannot be confirmed there either.


#### Remark (the threshold is not a constant, and the sweep held it fixed) {#rem:threshfam}
<!-- evidence: audit_threshold_arithmetic.py -->

The number $0.3745$ above, and everywhere else in this program, is
$\SS(N)(1-A(N))$ — and both factors depend on which primes divide
$N$. The $N$ of every sweep here are $2\cdot10^5\cdot2^{\,j}$, and
$2\cdot10^5=2^6\cdot5^5$, so all eight have odd radical $5$ and hence
*one* threshold. The sweeps move the size of $N$ by a factor $128$ and
its arithmetic not at all. Recomputing the constant from Euler
products truncated at $10^7$ gives $0.374487$, so the published value
is right; what was never asked is whether the verdict survives
changing the arithmetic.

It does not survive it, and by a wide margin. At seven $N$ in
$[1.40\cdot10^6,\,1.63\cdot10^6]$ chosen for their odd radicals, the
threshold runs from $0.073312$ at $N=2\cdot3^2\cdot5\cdot7\cdot11
\cdot13\cdot17$ to $0.374487$ at the published family — a factor of
five, and the published family sits at the *maximum*, not merely above
the median $0.270682$. The measured $|E_3|/N$ moves the other way over
the same seven, from $0.1289$ to $2.7763$; the two effects compound,
and [eq:onesided] holds at three of the seven and fails at four.
Against the constant $0.3745$ the verdict column is identical, at $0$
of $7$ — not because the confound is harmless but because the failures
are gross rather than marginal.

Two things follow, and they point opposite ways. The first is a
mechanism for the failures. When $N$ has several small odd prime
factors, $k$ must avoid them, and by [eq:dilate] so must $m$: for
$N=2\cdot3^2\cdot5\cdot7\cdot11\cdot13\cdot17$ the smallest admissible
$k$ is $19$, and the $m$ below $N/k$ that are coprime to $\rad(N)$ are
almost all primes, so $\mu(m)=-1$ dominates and $H(N;k)$ acquires a
sign. Measured, all $513$ terms of $E_3$ are negative and
$|E_3|=B(N)$ exactly: at these $N$ there is no cancellation across $k$
whatsoever. Compare $442$ against $562$ at the published $N$.

The second says the damage is not where it looks. Proposition
[prop:onesided] reaches its threshold by bounding the wall,
$|C(N)|\le A(N)N$, and that bound is slack: $|C|/(A N)$ measures
$0.0012$ to $0.0208$ across the seven. Putting the measured wall in
its place, Theorem [thm:C] asks for $|E_3|/N<\SS(N)(1-|C(N)|/N)$, and
*all seven* satisfy it — $2.7763$ against $5.3505$ at the worst one.
So the arithmetic dependence of the verdict lives entirely in the
slack of $|C|\le AN$, not in $E_3$. This is a measurement and not a
proof: $C(N)=o(N)$ is the open problem and nothing here supplies it.
What it locates is which of [prop:onesided]'s two inequalities is
doing the damage, and it is not the one this program has been
measuring.


#### Proposition (the wall cancels out of the count) {#prop:direct}
<!-- evidence: lab_direct_route.py -->

Substitute Proposition [prop:posweights],
$E_3=\sum_{k<K}\mu^2(k)(\log k)H(N;k)-C(N)B_{\log}(K)$, into the identity of
Theorem [thm:C]. The two occurrences of $C(N)$ carry opposite signs
and cancel:

$$
\begin{equation}\label{eq:direct}
  \tilde r(N) \;=\; \SS(N)N
   \;+\!\!\sum_{\substack{k<K\\(k,N)=1}}\!\!\mu^2(k)(\log k)H(N;k)
   \;-\; C(N)\bigl(B_{\log}(K)+\SS(N)\bigr)
   \;+\; O_A\!\left(\frac{N}{(\log N)^{A}}\right).
\end{equation}
$$

Since $B_{\log}(K)\to-\SS(N)$, the wall enters the Goldbach count only
through the vanishing factor $B_{\log}(K)+\SS(N)$. Consequently

$$
\begin{equation}\label{eq:directcond}
  \sum_{\substack{k<K\\(k,N)=1}}\mu^2(k)(\log k)\bigl|H(N;k)\bigr|
   \;\le\;(1-\varepsilon)\,\SS(N)\,N
\end{equation}
$$

suffices for $\tilde r(N)>0$, for all large even $N$.


**Proof.** 
Immediate from [eq:direct] and $|C(N)|\le A(N)N$, the last bound now
multiplied by $B_{\log}(K)+\SS(N)=o(1)$ rather than by $\SS(N)$.
 ∎


#### Remark (what the cancellation buys) {#rem:directmargin}
<!-- evidence: lab_direct_route.py -->

[eq:directcond] is the same shape as [eq:nolog] — a constant-factor
bound at level $N^{\theta'}$, so Remark [rem:relocate]'s verdict
stands and nothing has moved off the level axis. What has changed is
the constant, and it has changed by more than a constant's worth.

The old threshold is $\SS(N)(1-A(N))$ and the new one is $\SS(N)$.
Measured at $\theta'=0.56$ over $N=2\cdot10^5$ to $3.2\cdot10^6$,
$B_H(N)/(\SS(N)N)$ with
$B_H=\sum_{k<K}\mu^2(k)(\log k)|H(N;k)|$ reads
$0.4578,\,0.4064,\,0.4079,\,0.3769,\,0.3338$ — under $1$ throughout
and falling — against $2.1591,\,1.9747,\,1.9500,\,1.7483,\,1.5798$ for
the old ratio, which is above $1$ throughout. So [eq:directcond] is
*already satisfied* at every accessible $N$ where [eq:nolog] is not.

The gain is largest exactly where Remark [rem:threshfam] found the
collapse. Across the seven $N$ of that audit the new ratio reads
$0.2462,\,0.3970,\,0.3769,\,0.3124,\,0.4920,\,0.4716,\,0.4073$ against
the old $2.3920,\,1.7080,\,1.7483,\,5.6100,\,37.8696,\,28.2303,\,
1.6290$: the new ratio spans a factor $2.00$ where the old spans
$23.25$, a spread ratio of $0.0860$. Dropping $A(N)$ removes the
arithmetic sensitivity along with the slack, because it was the same
thing. The residual coupling is negligible at these $N$:
$|B_{\log}+\SS|/\SS$ stays under $0.0383$ and
$|C(N)(B_{\log}+\SS)|/N$ never exceeds $3.315\cdot10^{-4}$, four
orders below $\SS(N)$.

Two pre-registered rules of that run fail, and neither touches the
above. **Y1** asked $|B_{\log}+\SS|/\SS$ to fall monotonically; it
falls, rises at $N=1.6\cdot10^6$, then falls. Asking a
Möbius-weighted partial sum to converge monotonically was the error;
the level is what [eq:direct] uses and the level holds. **Y2** asked
the residual of [eq:direct] to be under $0.20$ at the largest $N$; it
is $0.2076$. The cap was an effect size and not a null, which is the
mistake Remark [rem:cap] recorded before, and the substantive half —
that the residual falls, $0.4559,\,0.3732,\,0.3110,\,0.2596,\,0.2076$
— held.

The refutation of Y2 is worth more than the rule was, because the
run's null came out inverted and explains it. A coin on
$\operatorname{supp}\mu^2$ gives a residual of median $0.0157$, a
factor $13.2$ *smaller* than $\mu$'s: for a coin both
$\sum(\log k)H_\varepsilon$ and $\tilde r-\SS N$ are separately near
zero, so the null as designed tests nothing about the cancellation.
What it exposes is that $\mu$'s residual is not noise but tracks
$\sum(\log k)H$ itself, at ratios
$1.0373,\,0.9409,\,0.9559,\,1.0124,\,0.9853$. At accessible $N$ the
unspecified $O_A$ term of Theorem [thm:C] is numerically the very
quantity [eq:directcond] must bound. That is a sharper form of the
known fact that [thm:C] cannot be checked numerically here: the
residual decays as $N^{-0.2794}$ with correlation $0.99930$ and
reaches $2\%$ of $\SS$ only near $N=10^{10.16}$, so no computation
below that separates the identity from its error term. [eq:directcond]
is therefore offered as a *weaker sufficient condition*, proved from
two identities, and not as something these measurements confirm.


#### Remark (where the direct condition's level actually sits) {#rem:directlevel}
<!-- evidence: lab_direct_level.py -->

[eq:directcond] is a demand at level $K$, so the question it leaves is
the only one the reduction cares about: where does
$B_H(N;K)=\sum_{k<K}\mu^2(k)(\log k)|H(N;k)|$ cross $\SS(N)N$? Write
$K^*_H(N)$ for that crossing and measure it in units of $\sqrt N$.

Over $N=2\cdot10^5$ to $3.2\cdot10^6$, walking $k$ upward until the
threshold is passed, $K^*_H$ reads $2973,\,5109,\,8021,\,13557,\,23397$
and $K^*_H/\sqrt N$ reads
$6.6478,\,8.0780,\,8.9678,\,10.7177,\,13.0793$ — above $1$ throughout
and growing. Fitted, $K^*_H\sim N^{0.7361}$ with correlation
$0.99954$. **At accessible $N$ the direct condition holds well past the
square-root barrier**, at a measured level near $\theta'=0.74$ where
the reduction needs any $\theta'>1/2$. Remark [rem:budget] audits that
exponent against the free parameter in it and it survives, which two
other exponents this program reported did not.

The control decides whether that is a fact about $\mu$. Replacing
$\mu$ by $\mu^2$ — same support, every sign $+1$, no cancellation —
gives $K^*=39,\,37,\,39,\,39,\,39$: bounded, independent of $N$, and
$K^*/\sqrt N$ falls from $0.0872$ to $0.0218$. The level is bought
entirely by cancellation and not at all by the size of the terms.

It is bought by *exactly* square-root cancellation, and that is the
sting. If $|H(N;k)|\asymp\sqrt{N/k}$ then
$B_H\asymp2\sqrt{NK}\log K$ and $K^*\asymp\SS^2N/(4\log^2K)$;
evaluated at the measured $K^*$ this predicts
$2423,\,4251,\,7669,\,13693,\,24496$ against the measured values, a
miss of $23\%$ at the bottom closing to $4\%$ at the top. Remark
[rem:heuristic] audits that agreement and finds it a coincidence of two
omitted factors; with both restored the prediction is within $1.5\%$ at
every $N$ and does not drift, so the law is confirmed more sharply than
stated here. So the
measurement does not reduce the difficulty — it identifies it exactly.
Everything [eq:directcond] needs, and nothing less, is square-root
cancellation for the dilated Möbius–prime correlation $H(N;k)$, which
is the object Remark [rem:dilateprofile] already identified as the
supply side's own consumable.

One thing this does correct. Remark [rem:levelmeas] withdrew an
earlier level measurement because a coin reached a higher $K^*$ than
$\mu$ at every $N$. That happens here too — the coin's $K^*/\sqrt N$
runs $11.5448$ to $22.7246$ against $\mu$'s $6.6478$ to $13.0793$ —
but Remark [rem:whycoinwins] has since explained it: by [eq:dilate]
the coin's progression sum is a sum of independent signs and gets
square-root cancellation for free, while $\mu$'s *is* the dilated
wall. The coin is better than $\mu$ by construction, so it is not a
null for a level measurement; it is a competitor. The withdrawal was
right about its control and wrong to be read as saying the level axis
cannot be measured. It can, against a control that is worse than $\mu$
rather than better.


#### Proposition (the untruncated sum is the count) {#prop:untrunc}
<!-- evidence: lab_direct_identity.py -->

Sum the direct route's terms over *all* $k$, with no truncation and no
coprimality restriction. Since $A(N;k)$ sums $n$ over $k\mid N-n$,
exchanging the order and applying $\sum_{d\mid u}\mu(d)\log d=-\Lambda(u)$
gives

$$
\begin{equation}\label{eq:untrunc}
  \sum_{k\ge1}(\log k)\,\mu(k)\,A(N;k)
   \;=\; -\sum_{n<N}\Lambda(n)\,\mu(N-n)\Lambda(N-n)
   \;=\; \sum_{p<N}\Lambda(N-p)\log p ,
\end{equation}
$$

the last step because $\mu(u)\Lambda(u)=-\log p$ at $u=p$ prime and
vanishes at every higher prime power. So the untruncated object *is*
$\tilde r(N)$, up to its prime-power part.


**Proof.** 
Both exchanges are finite sums. $\mu*\log=-\Lambda$ is Möbius
inversion of $\log = \mathbf 1 * \Lambda$.
 ∎


#### Remark (where the count actually lives) {#rem:whereitlives}
<!-- evidence: lab_direct_identity.py -->

[eq:untrunc] is not a decoration. It says [eq:direct] is a
*rearrangement* of $\tilde r(N)$ rather than an approximation waiting
for an estimate, so all of its content is in the truncation at $K$ —
and it lets one ask where the count is assembled. Verified at
$N=2\cdot10^5$ through $3.2\cdot10^6$, [eq:untrunc] holds to a worst
relative error of $3.136\cdot10^{-15}$, and the quantity is
$\SS(N)N$ to within a percent: the ratio reads
$1.0039,\,0.9865,\,0.9893,\,1.0017,\,0.9968$. Restricting $k$ to
$(k,N)=1$, as the Huang–Li setup does, changes the total by $0.0000$
at every $N$ — the moduli sharing a factor with $N$ contribute
nothing.

The pre-registered **Z3** fails and its failure is the finding. It
predicted that the truncated range carries under $10\%$ of the total.
It carries far more than $100\%$, with the wrong sign: the partial
sums over $k<N^{0.90}$ read
$-4.3412,\,-4.6951,\,-5.0055,\,-5.3075,\,-5.6291$ in units of $N$,
$2.4564$ to $3.2078$ times the total in magnitude and negative where
the total is positive. The band $k\in(N^{0.9},N/2]$ then supplies
$0.8034$ to $0.9146$ of the total, and everything above $N/2$ supplies
$2.6529$ to $3.2932$ of it. The count is not assembled gradually; it
is a difference of two large one-signed masses, and the truncation
cuts through the middle of the cancellation.

That is the sharpest statement so far of why $\theta'$ cannot simply
be pushed. Above $N/2$ only $m=1$ survives in
$H(N;k)=\sum_{m<N/k,(m,k)=1}\Lambda(N-mk)\mu(m)$, so that tail is
$\sum(\log k)\Lambda(N-k)$ — nonnegative, elementary, and carrying
three times the answer. The reduction truncates it away and asks the
remaining negative mass to be small, which it is not: it is $-5.6N$ at
$N=3.2\cdot10^6$ and growing. What [eq:directcond] demands is that the
*difference* between two quantities of size several $N$ be smaller
than $\SS(N)N$, and [eq:untrunc] says the difference is exactly
$\tilde r(N)$. The demand is therefore not merely as hard as binary
Goldbach at $K=N$; at $K=N$ it *is* binary Goldbach.


#### Proposition (the count as signed layers) {#prop:layers}
<!-- evidence: lab_layer_decomposition.py -->

Cut the same double sum along the inner variable of [eq:dilate]
instead. Exchanging the order in
$\sum_k\mu^2(k)(\log k)H(N;k)=\sum_k\mu^2(k)(\log k)\sum_{m<N/k,(m,k)=1}
\Lambda(N-mk)\mu(m)$,

$$
\begin{equation}\label{eq:layers}
  \sum_{k}\mu^2(k)(\log k)H(N;k) \;=\; \sum_{m}\mu(m)\,L(N;m),
  \qquad
  L(N;m) := \!\!\sum_{\substack{k<N/m,\ (k,m)=1\\ \mu(k)\neq0}}\!\!
            (\log k)\,\Lambda(N-mk) .
\end{equation}
$$

Every $L(N;m)$ is nonnegative. With [eq:untrunc], the Goldbach count
is therefore an alternating sum, signed by $\mu$, of nonnegative
layers.


**Proof.** 
A finite exchange; the coprimality conditions on the pair $(k,m)$ are
symmetric.
 ∎


#### Remark (what the layers weigh) {#rem:layerdecay}
<!-- evidence: lab_layer_decomposition.py -->

[eq:layers] holds to a worst relative error of $2.641\cdot10^{-15}$
over $N=2\cdot10^5$ to $1.6\cdot10^6$, reproducing [eq:untrunc]'s
totals digit for digit. What it exposes is the size of what cancels.
The first layer alone is $8.8266$ to $10.4611$ times $N$, which is
$5.0139$ to $5.9423$ times the whole answer; the layers together weigh
$24.6N$ to $33.0N$, so the surviving total is
$7.198\cdot10^{-2}$ down to $5.348\cdot10^{-2}$ of the mass. Growing
like $\log N$ on the way in and staying at $\SS(N)N$ on the way out,
this is a Mertens cancellation in $m$ and nothing gentler.

**Here $\mu$ beats a coin, and that is the point.** Remark
[rem:whycoinwins] showed a coin *beating* $\mu$ when the signs sit on
$k$, because there a coin buys square-root cancellation inside the
progression sum while $\mu$'s is the dilated wall. In [eq:layers] the
layers are fixed and nonnegative, so a sign draw cannot buy anything
inside them; it can only rearrange the pattern across $m$. Drawing
$\varepsilon(m)=\pm1$ on the squarefree $m\ge2$, sixteen times,
$|\sum\varepsilon(m)L|/N$ has minimum $1.8609,\,4.4112,\,4.6742,\,
7.0547$ against $\mu$'s $1.7673,\,1.7366,\,1.7416,\,1.7634$ — $\mu$ is
below every draw at every $N$, and the gap widens. The two remarks are
not in tension: they say the demand's difficulty lives in the outer
variable and its cancellation lives in the inner one.

Two further readings. The layers with $(m,N)>1$ are empty — if $q\mid
N$ and $q\mid m$ then $q\mid N-mk$, so $\Lambda(N-mk)$ survives only at
the single point $N-mk=q$ — and measured, the coprime layers carry
$0.999986$ to $0.999999$ of the mass. That is the dual of
[prop:untrunc]'s finding that restricting $k$ to $(k,N)=1$ changes the
total by nothing. And the pre-registered **Z3** fails: it asked the
partial sum at $M=N^{0.25}$ to still exceed twice the total, and the
ratio is $1.5163,\,1.6593,\,1.7090,\,1.6143$. The convergence in $m$
is faster than predicted — by $m\approx30$ the alternating sum is
already within a factor $1.7$ of its limit — so the cancellation is
not spread thinly over all $m$ but is essentially finished at small
$m$. The threshold $2$ was an effect size and not a null, which is
Remark [rem:cap]'s mistake; what the run measured is worth more than
what the rule asked.


#### Remark (the truncation in $m$ is cheap, and where the barrier goes) {#rem:layertail}
<!-- evidence: lab_layer_tail.py -->

$\Lambda(N-mk)$ is supported where $N-mk$ is a prime power, so writing
$p=N-mk$, the layer $L(N;m)$ of [eq:layers] is a sum over primes
$p<N$ with $p\equiv N\pmod m$, weighted by $\log\bigl((N-p)/m\bigr)$
and restricted to $(N-p)/m$ squarefree and coprime to $m$. The layers
are prime counts in progressions — the *supply* side's object — and
[eq:layers] turns the demand into a $\mu$-signed average of them. Two
numbers then decide what that is worth.

**The truncation in $m$ is cheap.** The tail
$\bigl|\sum_{m\ge M}\mu(m)L(N;m)\bigr|/N$ at
$M=10,\,30,\,100,\,300,\,1000,\,3000$ reads, at $N=1.6\cdot10^6$,
$3.5477,\,1.0484,\,0.4088,\,0.0660,\,0.0067,\,0.0015$; at $M=3000$ it
is between $0.0001$ and $0.0047$ of $\SS(N)$. Taking $M^*(N)$ to be
the least $M$ past which the tail stays under $0.01N$,
$M^*=8264,\,8624,\,9332,\,10352$ — a factor $1.25$ while $N$ grows by
$8$, so $M^*/\sqrt N$ falls from $18.4789$ to $8.1840$ across the
range.

**That last sentence is as far as the numbers reach, and an earlier
version of this remark went further than it should have.** It read a
fitted exponent off those four values and concluded that the
truncation in $m$ sits strictly inside the Bombieri–Vinogradov range
and moves further inside as $N$ grows. Remark [rem:tolerance]
withdraws the exponent: $M^*$ is defined by a tolerance chosen by
hand, and at other tolerances the fit gives a different answer with no
correlation to speak of. What survives is what is printed above —
$M^*$ is of order $10^4$ at these $N$ and grows slowly — together
with the threshold-free statement that the tail at a *fixed* cut
grows with $N$, so the truncation cannot be bounded.

That is as far as the good news goes, and the reason is worth stating
precisely. The cofactor $k=(N-p)/m$ is required *squarefree* — the
condition $\mu(k)\neq0$ that [eq:layers] inherits from
$H(N;k)=\mu(k)A(N;k)$, and it is not decoration. Recomputing the same
alternating sum over $m<2000$ without it gives
$2.5644,\,2.5168,\,2.5388,\,2.5479$ against
$1.7882,\,1.7690,\,1.7730,\,1.7819$ — a factor $1.4340,\,1.4227,\,
1.4319,\,1.4299$, flat in $N$. Detecting it costs
$\mu^2(k)=\sum_{d^2\mid k}\mu(d)$, which turns each layer into a sum
over moduli $d^2m$ with $d$ up to $\sqrt{N/m}$. **The cheap truncation
in $m$ does not remove the square-root barrier; it relocates it inside
the layer.** That is the same negative as Remark [rem:relocate], now
in the dual variable, and it is the first place in this program where
the two truncations can be compared on the same object.

Two pre-registered rules fail. **Z3** normalised the layers as
$R(m)=\varphi(m)L(N;m)\log N/\bigl(N\log(N/m)\bigr)$ and asked for
$R(m)\in[0.3,1.5]$; measured, $R$ runs $6.0127$ to $10.4611$. The band
was mis-specified by a constant — the substantive content, that the
$1/\varphi(m)$ law captures the layers up to a bounded factor, holds,
with $R$ spanning $7.2404$ to $10.4611$ at the largest $N$, a spread of
$1.44$. **Z4** predicted that the tail's smallness is a cancellation
only $\mu$ can supply, and it is not: $\mu$'s tail beats the minimum
of sixteen sign draws at two of four $N$ and not at the others. The
diagnosis is that past $M=3000$ the layers still carry $\ell^1$ mass
$1.4885$ to $3.9771$ but the largest single layer is $0.0020N$, so the
$\ell^2$ norm governing a random sign sum is only $0.0240$ to $0.0398$
— a coin cancels that tail as well as $\mu$ does. The cancellation
$\mu$ is genuinely needed for is spent at small $m$, which is what
Remark [rem:layerdecay] measured.


#### Proposition (the count over a combined modulus) {#prop:combined}
<!-- evidence: lab_combined_modulus.py -->

Expand the squarefree condition inside [eq:layers] by
$\mu^2(k)=\sum_{d^2\mid k}\mu(d)$ and write $k=d^2j$. Then

$$
\begin{equation}\label{eq:combined}
  \sum_{p<N}\Lambda(N-p)\log p
   \;=\!\!\sum_{\substack{m,d\ \text{squarefree}\\ (d,m)=1,\ md^2<N}}\!\!
     \mu(m)\mu(d)\,L_2(N;m,d),
\end{equation}
$$

$$
L_2(N;m,d) \;=\!\!\sum_{\substack{j<N/(md^2)\\ (j,m)=1}}\!\!
   \log(d^2j)\,\Lambda\bigl(N-md^2j\bigr),
$$

and $L_2$ carries no squarefree condition: it is a prime count in the
progression $p\equiv N\pmod{md^2}$, an unadorned Bombieri–Vinogradov
object of modulus $q=md^2$.


**Proof.** 
Finite exchange, and $(k,m)=1$ splits as $(d,m)=(j,m)=1$.
 ∎


#### Remark (the modulus the demand really needs) {#rem:combmod}
<!-- evidence: lab_combined_modulus.py -->

[eq:combined] makes the question a single number: how large must $q$
be allowed to grow? Remark [rem:layertail] found the truncation in
$m$ alone of order $10^4$ at these $N$ and argued that the squarefree
condition puts the barrier back. [eq:combined] lets that be measured
instead of argued.

Verified over $N=2\cdot10^5$ to $1.6\cdot10^6$ — $166384$ to $1331048$
pairs — [eq:combined] holds to a worst relative error of
$4.525\cdot10^{-15}$ and reproduces the totals of [eq:untrunc] and
[eq:layers] digit for digit. Truncating to $q<Q$ and taking $Q^*(N)$
to be the least $Q$ past which the deficit stays under $0.01N$,

$$
Q^*=160112,\ 200848,\ 482589,\ 589203,
\qquad
Q^*/\sqrt N = 358.0213,\ 317.5686,\ 539.5509,\ 465.8059,
$$

so $Q^*$ exceeds $\sqrt N$ by two orders at every $N$ tested, and
exceeds $M^*$ by a factor between $19$ and $57$. **Expanding the
squarefree condition costs one to two orders in the modulus**, which
is Remark [rem:layertail]'s relocation, measured.

An earlier version of this remark said more: it fitted
$Q^*\sim N^{0.6904}$, concluded that the combined modulus must reach a
power of $N$ above $1/2$, and observed that the $k$-side crossings of
Remark [rem:signedlevel] land near the same $N^{0.7}$, so that the
exponent looked like a property of the count rather than of the
dissection. Remark [rem:tolerance] withdraws all of that. The
exponent is not stable in the tolerance that defines $Q^*$, the ratio
$Q^*/\sqrt N$ has no consistent trend in $N$ across tolerances, and
the agreement with $0.6716$ and $0.7250$ was a coincidence of one
choice of cut-off. What is left is the comparison at a fixed
tolerance, which is what the paragraph above now claims and no more.

Two pre-registered rules fail. **Y2** asked the deficit to fall
monotonically over $Q=10^2$ to $10^6$; on the finer grid it reads
$0.22560,\,0.01496,\,0.01058,\,0.00000$ at the largest $N$ — three
orders of fall and then a wander at the $10^{-3}$ level, which is the
alternating sum finishing rather than failing to converge. Asking an
alternating sum for monotone convergence was the error. **Y4**
predicted that the convergence in $Q$ is bought by $\mu(d)$: it is
not. Holding every $L_2$ fixed and drawing the sign on $d$ at random,
sixteen draws bracket $\mu$'s deficit rather than losing to it —
$\mu$ gives $0.0733$ at the largest $N$ against a draw range of
$0.0038$ to $0.0825$. What makes the far pairs negligible is their
size, not their signs. The Möbius structure that does earn its keep is
the one on $m$, measured in Remark [rem:layerdecay], where $\mu$ beats
every draw.


#### Remark (the withdrawn level measurement, re-examined) {#rem:levelaudit}
<!-- evidence: audit_levelmeas_budget.py -->

Remark [rem:levelmeas] withdrew $K^*(N)$ on the strength of a coin
control. Two things have since made that verdict worth re-testing:
Remark [rem:whycoinwins] showed the coin wins here by construction,
and Remarks [rem:tolerance] and [rem:budget] established that a
crossing of a *monotone* sum survives a sweep of its free parameter
where a crossing of an alternating tail does not.
$B(N;K)=\sum_{k<K}(\log k)|\Emu(N;k)|$ is a sum of nonnegative terms,
so this crossing is in the surviving class.

It survives. Scaling the budget to $c\,\SS(N)(1-A(N))N$ and refitting
at $c=0.3,\,0.5,\,1,\,2,\,3$ gives $K^*$ of
$69,\,111,\,157,\,273,\,457$ at the loosest through
$1517,\,2463,\,3877,\,6701,\,11377$ at the tightest, increasing in $c$
at every $N$, with exponents

$$
e(c)=0.6753,\ 0.6813,\ 0.7057,\ 0.7299,\ 0.7258,
$$

spread $0.0545$ and correlations $0.99298$ to $0.99937$ — the
signature of the surviving class, against the $1.3033$ spread and
$0.00931$ correlations of the two withdrawn exponents. Every $e(c)$
exceeds $1/2$, the least being $0.6753$. The parameter-free check
agrees: $B(N;K)/\bigl(\SS(1-A)N\bigr)$ at $K<N^{0.50}$ reads
$1.3120,\,1.1459,\,1.1258,\,0.9509,\,0.8369$, falling through $1$
between $8\cdot10^5$ and $1.6\cdot10^6$, which is the same crossing
the fit reports and has no free parameter in it. Independently
recomputed here, $K^*$ at $c=1$ reproduces
$319,\,537,\,767,\,1353,\,2319$ exactly.

**The control the withdrawal should have used.** Replacing $\mu$ by
$\mu^2$ — same support, every sign $+1$, main term still subtracted,
so no cancellation available anywhere — gives $K^*$ of
$73,\,77,\,69,\,77,\,77$: bounded, independent of $N$, against $\mu$'s
$319$ to $2319$. The pre-registered **V4** asked its
$K^*/\sqrt N$ to be under $0.1$ at every $N$ and it reads
$0.1632,\,0.1217,\,0.0771,\,0.0609,\,0.0430$, so V4 is refuted at the
two smallest $N$ — the cap was an effect size and not a null, Remark
[rem:cap]'s mistake once more. What it was written to test holds
plainly: the control does not grow at all, so the level is bought by
cancellation and not by the size of the terms.

So Remark [rem:levelmeas]'s numbers stand and its verdict does not.
The withdrawal was right to demand a null and wrong to accept one that
had to win. What is claimed is narrow: at accessible $N$, and against
a control that cannot cancel, $\sum_{k<K}(\log k)|\Emu(N;k)|$ stays
under the Goldbach budget out to $K$ past $\sqrt N$, at a level whose
exponent is stable in the budget. Nothing about $\theta'$ asymptotic
follows; five values of $N$ cannot supply that, and Remark
[rem:relocate]'s obstruction is untouched.


#### Remark (every exponent, swept) {#rem:allswept}

Remarks [rem:tolerance], [rem:budget] and [rem:levelaudit] swept the
free parameter of the three crossings. The other six exponents these
notes quote are direct fits of measured values, whose free parameter
is not a budget but the $N$-range, and the corresponding check is a
leave-one-out refit — computed from data each script already holds, so
it costs nothing to run. Every results file that prints a fitted
exponent now emits the spread on a `SWEPT` line of its own, and the
gate refuses a file that prints an exponent without one. Each figure
below is read from the script that produced the exponent it belongs
to; there is no single source for the table.

The spreads are

$$
\begin{array}{lc}
\text{residual decay} & 0.0100\\
\text{flat sum decay} & 0.0094\\
|E_3| \text{ decay} & 0.0058\\
\text{lean decay} & 0.0208\\
B/N \text{ power} & 0.0084\\
\text{wall max decay} & 0.0173\\
K^*_H & 0.0143\\
M^* & 0.0441\\
Q^* & 0.1055
\end{array}
$$

and they rank the same way the budget sweeps did, without being told
to. $Q^*$ — the exponent withdrawn in Remark [rem:tolerance] — is the
least stable by a factor of two and a half over the next worst;
$M^*$, the other withdrawal, is second. The six direct fits sit
between $0.0058$ and $0.0208$, and $K^*_H$, the crossing that survived
its budget sweep, sits with them at $0.0143$. Two independent
robustness tests — varying the defining tolerance, and dropping an end
of the $N$-range — agree on which exponents are worth quoting.


#### Remark (the surviving exponent, audited) {#rem:budget}
<!-- evidence: audit_directlevel_budget.py -->

Remark [rem:tolerance] withdrew two exponents by sweeping the free
parameter that defined them. The same audit is owed to the third.
$K^*_H$ is the crossing of $B_H(N;K)=\sum_{k<K}\mu^2(k)(\log k)|H(N;k)|$ above
$\SS(N)N$, and [eq:directcond] actually asks for
$B_H\le(1-\varepsilon)\SS(N)N$, so the budget is a free parameter that
was set to $\varepsilon=0$ and never moved. Scaling it to $c\,\SS(N)N$
and refitting at $c=0.3,\,0.5,\,1,\,2,\,3$ gives

$$
e(c)=0.7133,\ 0.7126,\ 0.7361,\ 0.7257,\ 0.6932,
$$

a spread of $0.0429$ against the $1.3033$ that sank $Q^*$, with fit
correlations $0.99733$ to $0.99969$ against correlations as low as
$0.00931$ for $M^*$. All five exceed $1/2$, the smallest being
$0.6932$, and $K^*_H/\sqrt N$ at the largest $N$ runs $2.1013$ to
$93.0646$ across the budgets. **The exponent survives.**

The reason it survives is structural and worth stating, because it is
what distinguishes this crossing from the two that failed. $B_H$ is a
sum of *nonnegative* terms, so it is monotone in $K$ and its crossing
is unique; $Q^*$ and $M^*$ were crossings of an alternating tail,
which oscillates while decaying, so they read the last excursion
rather than a trend. The parameter-free check confirms it: at a fixed
cut, $B_H(N;K)/(\SS(N)N)$ reads
$0.4601,\,0.4081,\,0.4079,\,0.3772,\,0.3341$ at $K<N^{0.56}$ and
$1.3754,\,1.3566,\,1.3703,\,1.3452,\,1.2953$ at $K<N^{0.70}$ — falling
with $N$ at every fixed exponent, so the crossing moves outward, which
is the fit's conclusion reached with no free parameter in it. In the
withdrawn cases the parameter-free tables and the fitted exponents
disagreed; here they agree.


#### Remark (the exponents were the tolerance, not the count) {#rem:tolerance}
<!-- evidence: audit_truncation_exponent.py -->

$M^*$ and $Q^*$ are both defined as "the least truncation past which
the deficit stays under $0.01N$", and that $0.01$ was chosen by hand.
Sweeping it over $\varepsilon=0.3,\,0.1,\,0.03,\,0.01,\,0.003$ and
refitting at each value — from one enumeration of the pairs $(m,d)$,
since grouping them by $m$ reconstitutes $\mu(m)L(N;m)$ exactly — the
fitted exponent for $Q^*$ reads

$$
0.2658,\ 0.2922,\ -0.4303,\ 0.6904,\ 0.8730,
$$

a spread of $1.3033$, with fit correlations of $0.98006$, $0.99066$,
$-0.56215$, $0.96357$, $0.99474$. For $M^*$ it reads
$0.2238,\,0.0358,\,0.0022,\,0.1086,\,0.8116$, spread $0.8094$, with
correlations as low as $0.00931$. **All four pre-registered rules fail — X1 and X2 on the spreads, X3 because the exponent for $Q^*$ drops to $-0.4303$, X4 because the gap between the two goes negative with it — and that is the finding.** The exponents this
program reported are properties of the tolerance. The one value that
produced a clean-looking fit, $\varepsilon=0.01$, did so by luck: its
neighbours give noise.

The reason is visible in the mechanism. Fitting the deficit itself
against the truncation as a power law gives correlations of
$-0.42134$ to $-0.54921$ in $Q$: the deficit does not decay like a
power, it *oscillates* while decaying, because it is the tail of an
alternating sum. $Q^*(\varepsilon)$ therefore reads the last excursion
above $\varepsilon N$ rather than a trend, and inverting an erratic
envelope at four values of $N$ produces a slope that means nothing.

What is threshold-free, and is what should have been reported, is the
envelope at a *fixed* truncation — the largest deficit at or beyond
$Q$, in units of $N$:

| cut | $2\cdot10^5$ | $4\cdot10^5$ | $8\cdot10^5$ | $1.6\cdot10^6$ |
|---|---|---|---|---|
| $q<10^4$ | $0.03318$ | $0.04266$ | $0.06288$ | $0.08131$ |
| $q<10^5$ | $0.02248$ | $0.01444$ | $0.01953$ | $0.01172$ |
| $m<300$ | $0.11679$ | $0.13028$ | $0.16745$ | $0.16424$ |
| $m<3000$ | $0.01247$ | $0.01890$ | $0.02081$ | $0.02104$ |

Read across, that is the decay at fixed $N$; read down, it is how the
error at a fixed cut moves with $N$, with no free parameter in it. It
rises with $N$ at $q<10^4$, at $m<300$ and at $m<3000$, so both
truncations must be unbounded — which is the qualitative conclusion
those remarks were entitled to. Neither table supports a clean power
of $N$ over this range, and saying so is the correct report.


#### Remark (what the signs across $k$ are worth) {#rem:signedlevel}
<!-- evidence: lab_signed_level.py -->

[eq:direct] contains the signed sum; [eq:directcond] discards the
signs because a bound on $|H(N;k)|$ is what an estimate supplies.
Walking both first crossings upward over the squarefree $k<N/2$
coprime to $N$: the absolute sum crosses $\SS(N)N$ at
$K^*_H=2973,\,5109,\,8021,\,13557,\,23397$, exponents
$\log K^*/\log N$ of $0.6552$ to $0.6716$; the signed sum crosses
$-\SS(N)N$ at $5331,\,9077,\,16327,\,29373,\,52017$, exponents
$0.7030$ to $0.7250$. Keeping the signs is worth a factor $1.793$ to
$2.223$ in $K^*$ and about $0.053$ in $\theta'$ — a real gain, and a
bounded one.

Two qualifications, both now declared by the result file. Every $N$
here is $2^a5^b$ — **one odd radical** — and so is every $N$ in Remark
[rem:residuesigned], which runs the identical family. What the signs
across $k$ are worth is therefore known at one arithmetic type, both
for $H$ ($0.053$) and for the residue ($0.21$ to $0.29$, the largest
single discard in the chain); Remark [rem:residuearithmetic] shows the
underlying level moving across $\tfrac12$ as the radical changes, and
nothing here says the signs are worth the same there.

**And bounded, and smaller than the choice of budget.** Both columns
are crossed against $\SS(N)N$; Proposition [prop:nolog] asks for
$\SS(N)(1-A(N))N$, and Remark [rem:modeltransfer] measures that
difference at $0.1677$ in the exponent — **three times the whole gap
between the two columns here.** So the signed-versus-absolute question
that [eq:direct] and [eq:directcond] separate is a smaller effect than
which budget the crossing is read against, and no exponent from this
remark may be compared with one from Remark [rem:residuelevel] without
that correction.

The control says the gain is not $\mu$'s doing. Holding every
$|H(N;k)|$ fixed and drawing the sign of each term at random, all
$16$ draws fail to cross anywhere in the walk, at every $N$. So random
signs cancel and $\mu$'s do not: $\mu$ is *worse* than random here,
which is Remark [rem:signmass]'s lean seen at the level. Measured
directly, the signed sum reaches $3.3348$ times its own threshold by
the end of the walk while the absolute sum reaches $7.9$ times its
own, a ratio flat in $N$ at $0.2887$ to $0.2918$.

**A correction.** The first version of this measurement summed
$(\log k)A(N;k)$ where $(\log k)\mu(k)A(N;k)$ was meant — by
[eq:dilate] the missing factor is exactly $\mu(k)$ — and reported that
the signed sum never crosses at all. It does. The error was found by
[prop:untrunc]'s independent recomputation of the same object, which
is what that check is for; nothing in [prop:direct] or Remark
[rem:directlevel] depended on it, since both use $|H|$.


#### Remark (a level measurement, and its withdrawal) {#rem:levelmeas}
<!-- evidence: lab_level_coin_null.py -->

Once the demand is a constant-factor bound, the level at which it holds
becomes a *measurable* quantity rather than a hypothesis. Define

$$
\begin{equation}\label{eq:Kstar}
  K^*(N) \;:=\; \max\Bigl\{K \;:\;
    \sum_{\substack{k<K\\(k,N)=1}} (\log k)\bigl|\Emu(N;k)\bigr|
    \;\le\; \SS(N)\bigl(1-A(N)\bigr)N \Bigr\},
\end{equation}
$$

the largest truncation at which [eq:nolog] still holds. Huang–Li need
$K=N^{\theta'}$ for a single $\theta'>1/2$; measured,

$$
\begin{array}{r|ccccc}
 N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6 & 3.2\cdot10^6\\\hline
 K^*(N) & 319 & 537 & 767 & 1353 & 2319\\
 K^*(N)/\sqrt N & 0.7133 & 0.8491 & 0.8575 & 1.0696 & 1.2964
\end{array}
$$

so $K^*(N)\asymp N^{0.7057}$ over this range, and $K^*/\sqrt N$ rises
monotonically, crossing $1$ between $8\cdot10^5$ and $1.6\cdot10^6$.

**That reads as an empirical level past the square-root barrier,
and it is withdrawn.** It was published with no null, which is the
one thing Lemma [lem:coin] exists to forbid. Running the control —
$\varepsilon(v)=\pm1$ on $\{v:\mu(v)\ne0\}$ and zero elsewhere, with
the threshold, support, weight and $k$-range held identical, so that
the sign pattern is the only difference — gives, averaged over eight
draws,

$$
\begin{array}{r|ccccc}
 N & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5 & 1.6\cdot10^6
   & 3.2\cdot10^6\\\hline
 K^*_\mu & 319 & 537 & 767 & 1353 & 2319\\
 K^*_\varepsilon & 553.0 & 887.5 & 1386.8 & 2145.5 & 3423.2\\
 K^*_\mu/K^*_\varepsilon & 0.58 & 0.61 & 0.55 & 0.63 & 0.68
\end{array}
$$

The coin survives to a *higher* level than $\mu$ at every $N$
tested, with $\beta_\varepsilon=0.6534$ against
$\beta_\mu=0.7057$ — two fits over five points that do not separate.
So $K^*$ is a statement about the support of $\mu^2$ and about
square-root cancellation, not about $\mu$: the control's discrepancy
obeys $|E_\varepsilon(N;k)|\asymp k^{-a}$ with
$a=0.5082,\,0.4982,\,0.4978,\,0.5053,\,0.5056$, exactly the
square-root law, and $\mu$ does not beat it.

Four of the null's five pre-registered rules therefore fail — K1
($K^*_\mu>K^*_\varepsilon$), K2 (ratio $\ge3$), K3
($\beta_\varepsilon>\beta_\mu$) and K5 — and only K4, the coin's own
$k^{-1/2}$ scaling, holds. The refutation is of this note's own
previous claim and it is total: **no level statement is claimed
here.** What survives is negative and worth keeping: whatever
carries the Huang–Li route past $\sqrt N$, it is not visible in
$\sum_{k<K}(\log k)|\Emu(N;k)|$ at accessible $N$, because that sum
does not distinguish $\mu$ from a coin.

Remark [rem:levelaudit] revisits that verdict. The withdrawal was
right that the measurement had no null and right to insist on one; it
was wrong about which null. The coin is better than $\mu$ here by
construction — Remark [rem:whycoinwins] — so K1 and K2 could only
come out as they did, and a control that must win is not a control.
Against the reference that *cannot* win, $\mu$ replaced by $\mu^2$,
the numbers above survive.


The net progress toward the Goldbach conjecture in this note is
*zero*. Theorem [thm:A] removes precisely the half of the
$EH_\mu$ demand that carries no Goldbach content, and
Theorem [thm:C] shows the other half is the conclusion itself. What
is gained is structural: the demand collapses to one scalar, that
scalar is proved equivalent to the conclusion, and one defect of the
published paper is identified with its repair.


## Notation and the mechanism {#sec:mechanism}


$p$ always denotes a prime. $\varphi$ is Euler's function, $\tau_3$ the
$3$-fold divisor function, $\rad(u)=\prod_{p\mid u}p$. Constants
$c,c_1,\dots$ are positive and absolute unless subscripted;
implied constants may depend on $A$ and $\theta'$. We write
$P(N)$ for the set of primes dividing $N$.

The mechanism of Theorem [thm:A] is the following. The residue class
in [eq:Emu] is the fixed class $n\equiv N\pmod k$, i.e.

$$
n \equiv N \pmod{k}\quad\Longleftrightarrow\quad k \mid N-n .
$$

So the $k$-sum is a divisor sum over $u=N-n$, namely the incomplete sum
$\sigma_K(u)=\sum_{k\mid u,\ k<K,\ (k,N)=1}\mu(k)$, and the mechanism
has three steps rather than one.

*First* one **completes** it. By Lemma [lem:complete], for
squarefree $u$ the complete sum $\sum_{k\mid u,\,(k,N)=1}\mu(k)$ equals
$\mathbf 1_{\rad(u)\mid N}$, which forces $u\mid\rad(N)$ and so leaves
only $N^{o(1)}$ terms. *Then* the work sits entirely in the
complementary sum, over $k\ge K$ — and it is only on that range that
writing $u=mk$ gives

$$
m \;=\; u/k \;<\; N/K \;=\; N^{1-\theta'} \;\le\; N^{1/2-\delta}.
$$

*Finally*, for squarefree $u=mk$ the factorisation forces
$(m,k)=1$ and $\mu(u)\mu(k)=\mu(m)\mu^2(k)$: the Möbius factor on the
*long* variable $k$ cancels against itself and leaves the
nonnegative $\mu^2$, while the surviving Möbius sits on the
*short* variable $m$. That is the classical
"Möbius on the short variable plus Bombieri–Vinogradov"
configuration.

The completion step is load-bearing and easy to omit. Over the range
$k<K$ as it stands the cofactor $u/k$ runs up to $u$ itself, so the
bound $m<N^{1-\theta'}$ is simply false there and the surviving
Möbius would sit on the *long* variable — the configuration
this program records as having no known machine. Steps 1–3 of
\S[sec:proof] carry out the argument in that order.


## Proof of Theorem [thm:A
]{#sec:proof}

Fix $A>0$ and set

$$
\begin{equation}\label{eq:trunc}
  D_0 := (\log N)^{A+2}, \qquad E_0 := (\log N)^{2A+4} .
\end{equation}
$$

(The truncations must be allowed to depend on $A$; a fixed choice covers
only bounded $A$.)


### Step 0: the subtracted mean term


By [eq:Emu],

$$
\begin{equation}\label{eq:split0}
  T_1(t) \;=\; D(t) \;-\; C(t)\,B(K),
  \qquad
  D(t):=\sum_{\substack{k<K\\(k,N)=1}}\mu(k)
        \sum_{\substack{n\le t\\ n\equiv N\ (k)}}\Lambda(n)\mu(N-n),
\end{equation}
$$

where $B(K)=\sum_{k<K,(k,N)=1}\mu(k)/\varphi(k)$. Huang–Li's Lemma 1
(a form of Goldston–Y{ı}ld{ı}r{ı}m) gives
$B(K)\ll e^{-c_1\sqrt{\log K}}$, and trivially
$|C(t)|\le\sum_{n<N}\Lambda(n)\ll N$. Hence

$$
\begin{equation}\label{eq:meanterm}
  |C(t)B(K)| \;\ll\; N e^{-c\sqrt{\log N}}
\end{equation}
$$

uniformly in $t$.


### Step 1: divisor switching (an exact identity)


Since $n\equiv N\pmod k$ with $n\le t<N$ is the same as $u:=N-n$ being a
multiple of $k$ in $[N-t,N)$,

$$
\begin{equation}\label{eq:switch}
  D(t) \;=\; \sum_{N-t\le u<N}\Lambda(N-u)\,\mu(u)\,\sigma_K(u),
  \qquad
  \sigma_K(u):=\sum_{\substack{k\mid u,\ k<K\\ (k,N)=1}}\mu(k) .
\end{equation}
$$

This is a finite rearrangement, hence exact. (Machine-verified in
`code/audit\_switch\_identity.py`, which accumulates the two sides in
genuinely different index orders — the $k$-side by modular slices, the
$u$-side against a divisor-sum array — and finds them equal to
$10^{-16}$ relative to $N$.)


### Step 2: the complete divisor sum


#### Lemma {#lem:complete}
<!-- evidence: audit_switch_identity.py -->

For squarefree $u\ge1$,
$\displaystyle\sum_{k\mid u,\ (k,N)=1}\mu(k)=\mathbf 1_{\rad(u)\mid N}$.


**Proof.** 
The sum is multiplicative in $u$ with local factor $1$ at $p\mid N$ and
$1-1=0$ at $p\nmid N$. Hence it vanishes unless every prime of $u$
divides $N$.
 ∎


Splitting $\sigma_K(u)$ into the complete sum minus the tail $k\ge K$,
[eq:switch] becomes $D(t)=P(t)-R(t)$ with

$$
\begin{equation}\label{eq:PR}
  P(t)=\!\!\sum_{\substack{N-t\le u<N\\ \rad(u)\mid N}}\!\!\Lambda(N-u)\mu(u),
  \qquad
  R(t)=\!\!\sum_{N-t\le u<N}\!\!\Lambda(N-u)\mu(u)
       \sum_{\substack{k\mid u,\ k\ge K\\ (k,N)=1}}\mu(k).
\end{equation}
$$

In $P(t)$ the conditions $\mu(u)\ne0$ and $\rad(u)\mid N$ force
$u\mid\rad(N)$, so there are at most $2^{\omega(N)}=N^{o(1)}$ terms, each
$\ll\log N$:

$$
\begin{equation}\label{eq:P}
  P(t)\ll N^{o(1)} .
\end{equation}
$$


### Step 3: the residual, and the degeneracy lemma


Write $u=mk$ with $k\ge K$, so $m=u/k<N/K=M$. For squarefree $u$ the
factorisation forces $(m,k)=1$ and $\mu(u)\mu(k)=\mu(m)\mu^2(k)$, so

$$
\begin{equation}\label{eq:R1}
  R(t)=\sum_{m<M}\mu(m)
      \sum_{\substack{k\ge K,\ (k,m)=1,\ (k,N)=1\\ N-t\le mk<N}}
      \mu^2(k)\,\Lambda(N-mk).
\end{equation}
$$


#### Lemma (degeneracy) {#lem:degen}
<!-- evidence: analytic -->

The contribution to [eq:R1] of the terms with $(k,N)>1$ is
$O(N^{o(1)})$. The same bound holds for the terms in which the modulus
constructed in Step 4 is not coprime to $N$.


**Proof.** 
Suppose $p\mid(k,N)$ and $\Lambda(N-mk)\ne0$, say $N-mk=q^\ell$. Since
$p\mid k$ and $p\mid N$ we get $p\mid N-mk$, hence $q=p$ and
$n:=N-mk=p^\ell$ with $p\mid N$. The number of such $n<N$ is
$\le\sum_{p\mid N}\log N/\log p\ll(\log N)^2$; for each, the pairs
$(m,k)$ with $mk=N-n$ number $\le\tau(N-n)\ll N^{o(1)}$, and each term is
$\ll\log N$.
 ∎


By Lemma [lem:degen] we may drop the condition $(k,N)=1$ from
[eq:R1] at a cost of $O(N^{o(1)})$. This is essential: expanding
$\mathbf 1_{(k,N)=1}$ by Möbius over $e\mid N$ would multiply the number
of Bombieri–Vinogradov calls by $2^{\omega(N)}=N^{o(1)}$ and destroy the
budget.


### Step 4: unfolding $\mu^2$ and the coprimality


Insert $\mu^2(k)=\sum_{d^2\mid k}\mu(d)$ and
$\mathbf 1_{(k,m)=1}=\sum_{e\mid(k,m)}\mu(e)$ and truncate at $D_0,E_0$
of [eq:trunc].

*Tail $d>D_0$.* Such terms have $d^2\mid k$, hence
$n=N-mk\equiv N\pmod{md^2}$; the number of $n\le N$ in that progression
is $\ll N/(md^2)+1$ and each $\Lambda\ll\log N$. Summing,

$$
\ll \log N\sum_{m<M}\sum_{d>D_0}\Bigl(\frac{N}{md^2}+1\Bigr)
  \ll \frac{N(\log N)^2}{D_0} + \sqrt{NM}
  \ll \frac{N}{(\log N)^{A}} + N^{1-\theta'/2}.
$$


*Tail $e>E_0$.* Such terms have $e\mid k$ and $e\mid m$, hence
$n\equiv N\pmod{me}$, and

$$
\ll \log N\sum_{m<M}\sum_{\substack{e\mid m\\ e>E_0}}
      \Bigl(\frac{N}{me}+1\Bigr)
  \ll \frac{N\log N}{E_0}\sum_{m<M}\frac{\tau(m)}{m} + M
  \ll \frac{N(\log N)^{3}}{E_0}+M .
$$


Both are $\ll N(\log N)^{-A}$. In the remaining range put
$L:=\lcm(d^2,e)$ and $q:=mL$; the conditions $d^2\mid k$, $e\mid k$ and
$n=N-mk$ give $n\equiv N\pmod q$, and $mk<N$, $mk\ge\max(N-t,mK)$ give
$1\le n\le T_m(t)$ where

$$
\begin{equation}\label{eq:Tm}
  T_m(t) := \min\bigl(t,\ N-mK\bigr).
\end{equation}
$$

Therefore, with $\psi(y;q,a)=\sum_{n\le y,\,n\equiv a\,(q)}\Lambda(n)$,

$$
\begin{equation}\label{eq:R2}
  R(t) = \sum_{m<M}\mu(m)\sum_{d\le D_0}\mu(d)
         \sum_{\substack{e\mid m\\ e\le E_0}}\mu(e)\,
         \psi\bigl(T_m(t);\,q,\,N\bigr)
         \;+\;O\!\left(\frac{N}{(\log N)^{A}}\right).
\end{equation}
$$

By Lemma [lem:degen] we may further restrict to $(q,N)=1$: if
$g=(q,N)>1$ then $g\mid n$ forces $n$ to be a power of a prime dividing
$N$, and the total over the $\ll MD_0E_0\tau$-many triples is
$\ll N^{1-\theta'+o(1)}$.


#### Remark (the false-positive trap) {#rem:trap}

Restricting the *main terms* to classes with $(q,N)=1$ is not
cosmetic. Assigning a main term $T_m/\varphi(q)$ to a degenerate class
shifts the density of the whole computation by the factor
$N/\varphi(N)$. Carried into the $\log k$ branch of
Section [sec:C], that error produces an apparent *refutation*
of $EH_\mu$. It is recorded here because it is the most plausible
false-positive we have encountered.


Note the size of the modulus: $q=m\lcm(d^2,e)\le MD_0^2E_0$, so by
[eq:trunc]

$$
\begin{equation}\label{eq:level}
  q \;\le\; N^{1-\theta'}(\log N)^{4A+8} \;=\; N^{1/2-\delta}(\log N)^{4A+8}
  \;\le\; N^{1/2-\delta/2}
\end{equation}
$$

for $N$ large. This is the level at which Bombieri–Vinogradov is
applied.


### Step 5: Bombieri–Vinogradov


#### Lemma ($\tau$-weighted Bombieri–Vinogradov) {#lem:BV}
<!-- evidence: analytic -->

For all $A,B>0$ there is $C=C(A,B)$ such that for
$Q\le N^{1/2}(\log N)^{-C}$,

$$
\sum_{q\le Q}\tau_3(q)^{B}\ \max_{y\le N}\ \max_{(a,q)=1}
  \Bigl|\psi(y;q,a)-\frac{y}{\varphi(q)}\Bigr|
  \;\ll_{A,B}\; \frac{N}{(\log N)^{A}} .
$$


**Proof.** 
Write $\mathcal E_q$ for the inner maximum. By Cauchy–Schwarz,
$\sum_q\tau_3^B\mathcal E_q\le
(\sum_q\tau_3^{2B}\mathcal E_q)^{1/2}(\sum_q\mathcal E_q)^{1/2}$.
For the first factor use the trivial bound
$\mathcal E_q\ll (N/\varphi(q))\log N$, giving
$\ll N(\log N)^{C_1(B)}$; for the second use Bombieri–Vinogradov with
exponent $2A+C_1(B)$.
 ∎


A modulus $q$ arises in [eq:R2] from at most
$\sum_{m\mid q}\tau(q/m)^2 \le \tau(q)^3 \le \tau_3(q)^3$ triples
$(m,d,e)$: given $m\mid q$ there are at most $\tau(q/m)$ choices of $d$
with $d^2\mid q/m$ and at most $\tau(q/m)$ of $e$ with
$\lcm(d^2,e)=q/m$. Writing $\psi(T_m;q,N)=T_m/\varphi(q)+\mathcal E$
and applying Lemma [lem:BV] with $B=3$ and [eq:level],

$$
\begin{equation}\label{eq:R3}
  R(t) \;=\; \mathrm{MT}(t) \;+\; O\!\left(\frac{N}{(\log N)^{A}}\right),
  \qquad
  \mathrm{MT}(t)=\sum_{\substack{m<M\\ (m,N)=1}}\mu(m)\,T_m(t)\,c_{D_0,E_0}(m),
\end{equation}
$$

where
$c_{D_0,E_0}(m)=\sum_{d\le D_0}\sum_{e\mid m,\,e\le E_0}
\mu(d)\mu(e)\mathbf 1_{(d,N)=1}/\varphi(m\lcm(d^2,e))$.


### Step 6: the density


#### Lemma {#lem:density}
<!-- evidence: audit_density_identity.py -->

Let $m$ be squarefree with $(m,N)=1$ and put

$$
c(m):=\sum_{d\ge1}\ \sum_{e\mid m}
     \frac{\mu(d)\mu(e)\mathbf 1_{(d,N)=1}}{\varphi(m\lcm(d^2,e))} .
$$

Then $c(m)=A(N)\lambda(m)/m$ with $A(N)$ as in [eq:AN] and
$\lambda(m)=\prod_{p\mid m}\bigl(1-\tfrac{1}{p(p-1)}\bigr)^{-1}$.
Moreover $c_{D_0,E_0}(m)=c(m)+O(1/(mD_0))$.


**Proof.** 
The sum is multiplicative. If $p\nmid m$ then $p\nmid e$ and the local
factor is $1-1/\varphi(p^2)=1-1/(p(p-1))$ for $p\nmid N$, and $1$ for
$p\mid N$. If $p\mid m$ (so $p\nmid N$) the $p$-part of
$m\lcm(d^2,e)$ is $p^{1+\max(2v_p(d),v_p(e))}$, and the four choices
$(v_p(d),v_p(e))\in\{0,1\}^2$ contribute

$$
\frac{1}{p-1}-\frac{1}{p(p-1)}-\frac{1}{p^2(p-1)}+\frac{1}{p^2(p-1)}
  \;=\;\frac{1}{p}.
$$

Hence $c(m)=\prod_{p\mid m}p^{-1}\cdot
\prod_{p\nmid m,\ p\nmid N}(1-\tfrac{1}{p(p-1)})
= m^{-1}A(N)\lambda(m)$. The truncation error is the tail
$\sum_{d>D_0}1/\varphi(md^2)\ll 1/(mD_0)$ and similarly for $e$.
 ∎


#### Remark (the load-bearing line) {#rem:loadbearing}

Equivalently, for squarefree $m$,
$\sum_{g\mid m}\mu(g)/(\varphi(m/g)\,g\,\varphi(g))=1/m$: the local
factor is exactly $p^{-1}$, so the density exponent is exactly $1$ and
$1/\zeta$ occurs to the *first* power in Lemma [lem:mu]
below. A non-integral exponent would give, by Selberg–Delange, only
$(\log x)^{-c}$ for a fixed $c$ in Lemma [lem:mu], and
Theorem [thm:A] would be false. The identity is verified in exact
rational arithmetic for all squarefree $m<400$, with zero mismatches
(`code/audit\_density\_identity.py`).


#### Lemma {#lem:mu}
<!-- evidence: analytic -->

Let $f(m):=\mu(m)\lambda(m)\mathbf 1_{(m,N)=1}$ and
$G(x):=\sum_{m\le x}f(m)/m$. Then
$G(x)\ll \frac{N}{\varphi(N)}\,e^{-c\sqrt{\log x}}
 + x^{-1/4}e^{C\sqrt{\log N}}$
for $x\ge2$.


**Proof.** 
Write $F(s)=\sum_m f(m)m^{-s}=\prod_{p\nmid N}(1-\lambda(p)p^{-s})$ and
$F(s)=\zeta(s)^{-1}H(s)$, so $f=\mu*h$ with $h$ multiplicative,
$h(p^j)=1$ for $p\mid N$ and $h(p^j)=1-\lambda(p)=O(p^{-2})$ for
$p\nmid N$. Then
$G(x)=\sum_{b\le x}\frac{h(b)}{b}\sum_{a\le x/b}\frac{\mu(a)}{a}$.
For $b\le\sqrt x$ use
$\sum_{a\le y}\mu(a)/a\ll e^{-c\sqrt{\log y}}$ (the classical
zero-free region) together with
$\sum_b|h(b)|/b\ll\prod_{p\mid N}(1-1/p)^{-1}\ll N/\varphi(N)$. For
$b>\sqrt x$ bound $\sum_{b>y}|h(b)|/b\le
y^{-1/2}\sum_b|h(b)|b^{-1/2}\ll y^{-1/2}\prod_{p\mid N}(1-p^{-1/2})^{-1}
\ll y^{-1/2}e^{C\sqrt{\log N}}$.
 ∎


Since $M=N^{1-\theta'}$ is a fixed power of $N$, the second term of
Lemma [lem:mu] is negligible at $x\asymp M$, and
$N/\varphi(N)\ll\log\log N$; so $G(x)\ll e^{-c'\sqrt{\log x}}$ in the
relevant range.


### Step 7: the main term dies, uniformly in $t$


#### Proposition {#prop:MT}
<!-- evidence: analytic -->

$\displaystyle \sup_{1\le t<N}|\mathrm{MT}(t)|\ll N e^{-c\sqrt{\log N}}$.


**Proof.** 
By Lemma [lem:density],
$\mathrm{MT}(t)=A(N)\sum_{m<M}\frac{f(m)}{m}T_m(t)+O(M\log N/D_0)$
with $T_m(t)$ as in [eq:Tm]. Put $m_0:=(N-t)/K$, so
$T_m(t)=t$ for $m\le m_0$ and $T_m(t)=N-mK$ for $m>m_0$; in particular
$x\mapsto T_x(t)$ is non-increasing, is constant on $[1,m_0]$, has
derivative $-K$ on $(m_0,M)$, and $T_M(t)=0$. Abel summation against
$G$ of Lemma [lem:mu] gives

$$
\sum_{m<M}\frac{f(m)}{m}T_m(t)
  = G(M)\,T_M(t) + K\!\int_{m_0}^{M}\! G(x)\,dx
  = K\!\int_{m_0}^{M}\! G(x)\,dx .
$$

Hence, uniformly in $t$,

$$
\Bigl|\sum_{m<M}\frac{f(m)}{m}T_m(t)\Bigr|
  \le K\int_{1}^{M}|G(x)|\,dx
  \ll K\int_1^M e^{-c\sqrt{\log x}}dx
  \ll KM e^{-c'\sqrt{\log M}} \ll N e^{-c''\sqrt{\log N}},
$$

using $\log M\asymp\log N$. The truncation error is
$\ll M\log N/D_0=o(N(\log N)^{-A})$.
 ∎


Note that the cancellation is genuine and global: the term $m=1$ alone
contributes $T_1(t)\asymp N$, so the smallness of $\mathrm{MT}$ is
entirely due to the cancellation in $\sum f(m)/m$.


### Conclusion


Combining [eq:split0], [eq:meanterm], [eq:P],
[eq:R3] and Proposition [prop:MT]:

$$
T_1(t) = \underbrace{P(t)}_{N^{o(1)}}
         - \underbrace{\mathrm{MT}(t)}_{\ll Ne^{-c\sqrt{\log N}}}
         - \underbrace{O\!\left(N(\log N)^{-A}\right)}_{\text{BV}}
         - \underbrace{C(t)B(K)}_{\ll Ne^{-c\sqrt{\log N}}},
$$

uniformly in $t<N$. This proves Theorem [thm:A], and exhibits
Bombieri–Vinogradov as the only ingredient that does not give an
exponential saving. \qed


#### Remark {#rem:bound}

An earlier version of this statement claimed the bound
$Ne^{-c\sqrt{\log N}}$ outright. That is not justified: the
Bombieri–Vinogradov input yields $N(\log N)^{-A}$ for each fixed $A$
and no more, because the Siegel–Walfisz range in its proof is
$q\le(\log N)^{B}$ with $B$ fixed. The corrected statement is the one
given above; it is what Corollary [cor:B] needs.


## Proof of Corollary [cor:B
]

Fix $k$ and set
$a_n:=\Lambda(n)\mu(N-n)\bigl(\mathbf 1_{n\equiv N (k)}-1/\varphi(k)\bigr)$,
so that $\Emu(t;k)=\sum_{n\le t}a_n$ and the bracket in [eq:E4] is
$\sum_{n<N}a_n\log(N-n)$. Since $\log(N-n)$ vanishes at $n=N-1$ and has
derivative $-1/(N-t)$, Abel summation gives

$$
\sum_{n<N}a_n\log(N-n) \;=\; \int_{1}^{N-1}\frac{\Emu(t;k)}{N-t}\,dt .
$$

The weight $\log(N-n)$ does not depend on $k$, so multiplying by
$\mu(k)$ and summing over $k<K$ with $(k,N)=1$ (a finite sum) gives the
exact identity

$$
E_4(\alpha) \;=\; \int_{1}^{N-1}\frac{T_1(t)}{N-t}\,dt .
$$

Hence $|E_4(\alpha)|\le\bigl(\sup_{t<N}|T_1(t)|\bigr)\log N
\ll_A N(\log N)^{1-A}$ by Theorem [thm:A]. As $A$ is arbitrary this
is Corollary [cor:B]. Huang–Li's Lemma 4 is invoked only to bound
$E_4$; with the above it is not needed, and the only remaining
appearance of $EH_\mu$ in their §3 is through $E_3(\alpha)$. \qed


## Equation (18) as published, and the range it drops {#sec:delta}


The published form of equation (18) drops a range. What follows
records the defect and one way of closing it. Neither is claimed as
new: the observation is S. Zheleznov's, and the authors have since
corrected the manuscript by a different route, in which the dropped
range never appears — see [rem:movingswitch]. What is set down here is
kept because the mechanism it uses is the one Theorem [thm:A] rests on,
and because the two routes meet at the same configuration.

Huang–Li define
$S_2(\alpha)=\sum_{n<N}\Lambda(n)\,\mu^2(N-n)\,\tilde\Lambda_\alpha(N-n)$
with
$\tilde\Lambda_\alpha(u)=\sum_{d\mid u,\ d>\alpha}\mu(d)\log(1/d)$, and
substitute $k=u/d$, so that the constraint $d>\alpha$ becomes

$$
\begin{equation}\label{eq:constraint}
  k \;<\; \frac{N-n}{\alpha},
\end{equation}
$$

an *$n$-dependent* bound. Their displayed equation (18) reads

$$
S_2(\alpha) \;=\; \sum_{k<\frac{N-1}{\alpha}}\mu(k)
    \sum_{\substack{n<N\\ n\equiv N\ (k)}}
    \Lambda(n)\mu(N-n)\log\Bigl(\frac{k}{N-n}\Bigr),
$$

in which [eq:constraint] has been replaced by the $n$-free bound
$k<(N-1)/\alpha$. The two differ. Writing $N-n=mk$, the condition
[eq:constraint] is $m>\alpha$, so the right-hand side of (18)
contains, in addition to $S_2(\alpha)$, exactly the terms with
$m\le\alpha$. On those terms
$\mu(k)\mu(N-n)=\mu(m)\mu^2(k)\mathbf 1_{(k,m)=1}$ and
$\log(k/(N-n))=-\log m$, so

$$
\begin{equation}\label{eq:Delta}
  S_2(\alpha) \;=\; \text{RHS of (18)} \;+\; \Delta,
  \qquad
  \Delta \;=\; \sum_{2\le m\le\alpha}\mu(m)\log m
    \sum_{\substack{k<(N-1)/\alpha\\ (k,m)=1}}\mu^2(k)\,\Lambda(N-mk).
\end{equation}
$$

(The term $m=1$ drops out since $\log1=0$.) The trivial bound is
$\Delta\ll N(\log N)^2$, which exceeds the target
$O(N(\log N)^{-A})$, so $\Delta$ is not negligible and needs its own
lemma.

It does close, by the machinery already used above and under hypotheses
Huang–Li already assume. Indeed $\Delta$ has exactly the shape of the
residual [eq:R1]: the Möbius factor sits on the *short*
variable $m\le\alpha$, and the long variable $k$ carries only
$\mu^2\ge0$. Since $\mu(m)\ne0$ on the range of summation, the two
conditions on $k$ collapse into one,

$$
\begin{equation}\label{eq:mu2mk}
  \mu^2(k)\,\mathbf 1_{(k,m)=1} \;=\; \mu^2(mk),
\end{equation}
$$

and expanding $\mu^2(mk)=\sum_{b^2\mid mk}\mu(b)$ with the truncation
$b\le B:=(\log N)^{A+4}$ — the same one Huang–Li use in their §3.1 —
reduces $\Delta$ to sums of $\Lambda$ over progressions to moduli
$[m,b^2]\le\alpha B^2=\alpha(\log N)^{2A+8}$, plus a main term
$A(N)\sum_{m\le\alpha}\mu(m)\lambda(m)(\log m)\,T_m/m$ which is
$O(Ne^{-c\sqrt{\log N}})$ by Lemma [lem:mu] and partial summation.
The level is then exactly the one their hypothesis supplies, with no
enlargement of $A$. Hence:


- In the Corollary-1 regime, $\alpha=N^{1-\theta'}<N^{1/2}$ and
      Bombieri–Vinogradov closes $\Delta$ *unconditionally*.

- In the general Theorem-1 regime, $\alpha\asymp N^{\theta}$ and
      $\Delta$ closes under the hypothesis
      $EH(N^{\theta}(\log N)^{2A+8})$ which is assumed there.


So [HL]'s Theorem 1 and their Corollary 1 stand as stated. In the
published text the treatment of $\Delta$ is missing; in the corrected
text $\Delta$ does not arise, because the truncation is not moved off
the inner variable in the first place.




#### Remark (the cut has since moved, and the switch survives the move) {#rem:movingswitch}
<!-- evidence: audit_moving_switch.py -->

The omission recorded above was reported to the authors by
S. Zheleznov, and they have prepared a corrected version of the
manuscript in which the moving truncation is kept throughout. The
correction does not restore the missing range; it removes the need for
one. Writing $Y_k=\lceil N-\alpha k\rceil-1$ and
$R_n=\lceil (N-n)/\alpha\rceil-1$, so that

$$
n\le Y_k \iff n<N-\alpha k \iff k\le R_n ,
$$

their (18) becomes a sum over $k<K$ whose *inner* range stops at
$Y_k$. The terms with $m\le\alpha$ never enter, and $\Delta$ does not
arise. The corrected argument also carries $EH_\mu$ at the endpoint
$y=Y_k$, so the maximum over $y$ in [eq:EHmu] becomes load-bearing
where it was not before.

**Every statement here was proved against the fixed cut**, and is
therefore a statement about the superseded formulation unless the
mechanism survives the move. It does, and it arrives one step earlier.
Exchanging under the moving cut gives

$$
D^\ast=\sum_{u<N}\Lambda(N-u)\mu(u)\,\sigma^\ast(u),
\qquad
\sigma^\ast(u)=\!\!\sum_{\substack{k\mid u,\ (k,N)=1\\ u/k>\alpha}}\!\!\mu(k),
$$

the outer bound $k<K$ being implied by $u<N$. So the incomplete divisor
sum is now cut by the **cofactor**. Completing it by
Lemma [lem:complete] leaves

$$
D^\ast=P^\ast-R^\ast,\qquad
R^\ast=\sum_{m\le\alpha}\mu(m)\!\!\sum_{\substack{(k,mN)=1\\ mk<N}}\!\!
       \mu^2(k)\,\Lambda(N-mk),
$$

with $P^\ast\ll N^{o(1)}$ exactly as in [eq:P]. Under the fixed cut the
completion had to *produce* a short variable, and did so only because
$k\ge K$ forced $m<N/K$; under the moving cut the tail is indexed by
$m\le\alpha$ from the start. The favourable configuration — $\mu$ on
the short variable, $\mu^2\ge0$ on the long one — is reached without
that step.

$R^\ast$ is the shape of [eq:Delta], and it closes the same way: by
[eq:mu2mk] the two conditions on $k$ collapse to $\mu^2(mk)$, and
expanding with the truncation $b\le(\log N)^{A+4}$ gives moduli
$[m,b^2]\le\alpha(\log N)^{2A+8}$. In the Corollary-1 regime
$\alpha=N^{1-\theta'}<N^{1/2}$, so Bombieri--Vinogradov closes it
unconditionally.

The subtracted mean term moves the same way. Exchanging endpoints,

$$
M^\ast=\sum_{\substack{k<K\\(k,N)=1}}\frac{\mu(k)}{\varphi(k)}\,C(Y_k)
      =\sum_{n<N}\Lambda(n)\mu(N-n)\,\rho_N(R_n+1),
$$

with $\rho_N$ as in Lemma [lem:extract]. Splitting at $N-n=H$ with
$H=2\alpha N^{(1-\theta)/2}$ gives $R_n\ge N^{(1-\theta)/2}$ on the bulk,
where Huang--Li's Lemma 1 supplies $e^{-c\sqrt{\log N}}$, and $O(H\log N)$
on the boundary; since $H=2N^{(1+\theta)/2}$ and $\theta<1/2$ there, both
are below $N(\log N)^{-A}$.

**Verified as identities, not as estimates**
(`code/audit\_moving\_switch.py`; X1--X4 hold, X5 reported not judged).
Over $N=2\cdot10^5$ to $3.2\cdot10^6$ at $\alpha=N^{0.44}$, each
rearrangement is accumulated in two genuinely different index orders
and the pairs agree to $6.4\cdot10^{-17}$ (X1, the switch),
$2.5\cdot10^{-16}$ (X2, the completion) and $2.1\cdot10^{-18}$ (X4, the
mean term) relative to $N$; the largest $m$ carrying a surviving term
is $213$, $291$, $393$, $533$, $727$ against $\alpha=215$, $292$,
$396$, $537$, $728$, so the tail is on the short variable at every $N$
(X3); and $|M^\ast|/N$ reads $0.000461$, $0.000059$, $0.000026$,
$0.000045$, $0.000015$.

One of the four checks first failed at a tolerance of $10^{-12}$
relative to $N$ and passed at $10^{-19}$ once the two accumulations
were made exactly rounded rather than running float. Nothing in the
formulas changed. $M^\ast$ is a sum whose cancellation reaches
$10^{-5}$ of its own mass, so a tolerance set against $N$ was below
what the arithmetic could deliver — the same defect this repository
records as `TOL BELOW PRINT`.

**What this does not do.** It does not restate Theorem [thm:A] or
Theorem [thm:C] against the corrected formulation; it establishes that
the mechanism they rest on transfers. The statements are still written
for the fixed cut, and the endpoint $Y_k$ appears in none of them.


#### Remark (the two cuts are a $\Delta$ apart) {#rem:cutbridge}
<!-- evidence: audit_delta_bridge.py -->
<!-- evidence: audit_cut_bridge.py -->

[rem:movingswitch] shows the mechanism transfers, which leaves open the
question it does not answer: is the fixed-cut statement *superseded* by
the corrected formulation, or is it the stronger of the two? The answer
is the second, and the term that separates them is the one this section
is about.

**The fixed-cut object is the larger.** Write
$J=T_1(N-1)=\sum_{k<K,(k,N)=1}\mu(k)\Emu(N-1;k)$ for the fixed cut at
its extreme endpoint and $T_1^\ast=\sum_k\mu(k)\Emu(Y_k;k)$ for the
moving cut. Over $N=2\cdot10^5$ to $3.2\cdot10^6$ at $\alpha=N^{0.44}$,
$|T_1^\ast|$ is $0.2910$, $0.1054$, $0.2496$, $0.2290$, $0.2126$ of
$\sup_{t<N}|T_1(t)|$ — the supremum taken over every $t$, not a grid.
Moreover $|J|/N$ reads $0.078425$, $0.062218$, $0.048288$, $0.038048$,
$0.028801$ against a supremum of $0.078491$, $0.062218$, $0.048305$,
$0.038048$, $0.028801$: the supremum is essentially attained at the
extreme endpoint, so nothing below starts from a convenient point of
Theorem [thm:A].

Uniformity in the outer truncation is free where the mechanism works.
Maximising $|T_1(t;K')|$ over both $t$ and $K'\in[N^{1/2},K]$ returns
$1.0000$ times the value at $K'=K$ at every $N$ — the maximum is
attained at the full truncation. ($K'<N^{1/2}$ is printed and not
judged: there the completion leaves a cofactor that is not short, and
no computation at these $N$ can show that.)

**The bridge is not a telescoping.** Since $Y_k$ decreases in $k$,

$$
T_1^\ast=T_1(Y_J;K)+\sum_{j<J}\bigl[T_1(Y_j;j+1)-T_1(Y_{j+1};j+1)\bigr]
$$

is an exact identity in the two-parameter family (verified to $0$ and
$2.3\cdot10^{-18}$ relative to $N$), but its increments live on
intervals of length about $\alpha$ and there are $K$ of them. Their
total $\sum_j|\cdot|$ falls like $N^{-0.1707}$, against the
$N^{-0.2200}$ that square-root cancellation over an interval of length
$\alpha$ would give. **That agreement is the reason the route is
useless.** A power saving here *is* a short-interval estimate for
$\Lambda$ on intervals of length $N^{1-\theta'}$, which no
unconditional argument supplies. The numbers say the bridge is cheap
and no proof follows. The pre-registered rule Y3 of
`code/audit\_cut\_bridge.py` asked instead whether the cost exceeded
$N$; it does not, so **Y3 is refuted as registered**, and the rule is
left standing with its verdict because the question it asked was the
wrong one and the record should show what was asked before the numbers
were seen. Y5, added afterward and disclosed as such in that script,
carries the question Y3 meant to ask.

**The bridge is $\Delta$.** The difference between the two cuts is
carried by the $n$ the moving cut drops, and there

$$
n>Y_k \iff N-n\le\alpha k \iff m\le\alpha,
\qquad m=\frac{N-n}{k},
$$

so on exactly those terms $\mu(k)\mu(N-n)=\mu(k)\mu(mk)
=\mu(m)\mu^2(k)$ — the Möbius factor on the *short* variable, the long
one carrying only $\mu^2\ge0$. That is [eq:Delta]'s shape with weight
$1$ in place of $\log m$. Hence

$$
\begin{equation}\label{eq:cutbridge}
  T_1^\ast=J-\Xi,\qquad \Xi=\Xi_{\mathrm{main}}-\Xi_{\mathrm{mean}},
\end{equation}
$$

$$
\Xi_{\mathrm{main}}=\sum_{m\le\alpha}\mu(m)\!\!
  \sum_{\substack{k<K\\ (k,mN)=1}}\!\!\mu^2(k)\,\Lambda(N-mk),
\qquad
\Xi_{\mathrm{mean}}=\!\!\sum_{\substack{k<K\\ (k,N)=1}}\!\!
  \frac{\mu(k)}{\varphi(k)}\bigl(C(N-1)-C(Y_k)\bigr).
$$

$\Xi_{\mathrm{main}}$ closes by [eq:mu2mk] and the $b\le(\log N)^{A+4}$
expansion exactly as $\Delta$ does, and more easily, there being no
$\log m$ to remove by partial summation; $\Xi_{\mathrm{mean}}$ is
$W_K\,C(N-1)$ less the term treated in [rem:movingswitch], and
$W_K=\sum_{k<K,(k,N)=1}\mu(k)/\varphi(k)\ll e^{-c\sqrt{\log N}}$ by
Lemma [lem:mu], which absorbs the trivial $C(N-1)\ll N\log N$.

**So the fixed cut is not superseded — it is one of two ingredients,
and the omitted term is the other.** [eq:cutbridge] was verified as an
identity, not as an estimate (`code/audit\_delta\_bridge.py`; Z1--Z3
hold): the two sides share no index order, $T_1^\ast$ being built from
the cofactor cut and $J$ from the full truncation, and they agree to
$2.274\cdot10^{-18}$ relative to $N$; $\Xi$ accumulated over residue
classes agrees with $\Xi_{\mathrm{main}}-\Xi_{\mathrm{mean}}$
accumulated over $m$ to $4.547\cdot10^{-18}$; and the largest $m$
carrying a surviving term is $213$, $291$, $393$, $533$, $727$ against
$\alpha=215$, $292$, $396$, $537$, $728$. The bridge is not a small
correction to be waved through: $|\Xi|/N$ reads $0.055582$,
$0.055661$, $0.036230$, $0.029335$, $0.022678$, the same order as
$|J|/N$. It needs $\Delta$'s estimate, not a triangle inequality.

[eq:cutbridge] is what Proposition [prop:movingcut] is built on:
Theorem [thm:A] at $t=N-1$, plus the estimate of Section [sec:delta]
with weight $1$, gives the moving-cut bound. The correction the
published equation needed and the formulation that removes the need for
it are two readings of one identity.


#### Proposition (the bound survives the corrected cut) {#prop:movingcut}
<!-- evidence: analytic -->

Fix $\theta'\in(1/2,1)$, put $\alpha=N^{1-\theta'}$ and
$K=(N-1)/\alpha$, and let $Y_k=\lceil N-\alpha k\rceil-1$. Then for
every $A>0$,

$$
\begin{equation}\label{eq:movingcut}
  T_1^\ast \;:=\; \sum_{\substack{k<K\\(k,N)=1}}\mu(k)\,\Emu(Y_k;k)
  \;\ll_{A,\theta'}\; \frac{N}{(\log N)^{A}},
\end{equation}
$$

unconditionally.


**Proof.** 
Write $J=T_1(N-1)$ and $a_n=\Lambda(n)\mu(N-n)$. Since
$\Emu(N-1;k)-\Emu(Y_k;k)$ is carried by $Y_k<n<N$,

$$
J-T_1^\ast=\Xi_{\mathrm{main}}-\Xi_{\mathrm{mean}},
\qquad
\Xi_{\mathrm{main}}=\sum_{\substack{k<K\\(k,N)=1}}\mu(k)
  \!\!\sum_{\substack{Y_k<n<N\\ n\equiv N\ (k)}}\!\! a_n ,
$$

$\Xi_{\mathrm{mean}}=\sum_{k<K,(k,N)=1}\frac{\mu(k)}{\varphi(k)}
\bigl(C(N-1)-C(Y_k)\bigr)$. We bound the three pieces.

*The fixed-cut term.* $J\ll_{A,\theta'}N(\log N)^{-A}$ is
Theorem [thm:A] at $t=N-1$, which lies in the range $1\le t<N$ of its
supremum. This is the only place $\theta'>1/2$ is used, and it is used
exactly where Theorem [thm:A] needs it.

*The bridge term.* Put $N-n=mk$. Then $n>Y_k$ iff $n\ge N-\alpha k$
iff $m\le\alpha$, and $m\le\alpha$ together with $k<K=(N-1)/\alpha$
forces $mk<N-1$, so the constraint $n\ge1$ is automatic. On these
terms $\mu(k)\mu(N-n)=\mu(k)\mu(mk)$, which vanishes unless $mk$ is
squarefree, in which case it equals $\mu(m)\mu^2(k)$. Hence

$$
\Xi_{\mathrm{main}}=\sum_{m\le\alpha}\mu(m)
  \!\!\sum_{\substack{k<K\\ (k,mN)=1}}\!\!\mu^2(k)\,\Lambda(N-mk).
$$

The term $m=1$ is $\sum_{k<K,(k,N)=1}\mu^2(k)\Lambda(N-k)\le
\sum_{N-K<j<N}\Lambda(j)\ll K\log N=N^{\theta'}\log N$, which is
$\ll_A N(\log N)^{-A}$ because $\theta'<1$. The terms $2\le m\le\alpha$
are [eq:Delta] with the weight $\log m$ replaced by $1$ and are closed
by the argument of Section [sec:delta] verbatim: by [eq:mu2mk] the two
conditions on $k$ collapse to $\mu^2(mk)$, expanding
$\mu^2(mk)=\sum_{b^2\mid mk}\mu(b)$ with the truncation
$b\le B:=(\log N)^{A+4}$ leaves sums of $\Lambda$ over progressions to
moduli $[m,b^2]\le\alpha B^2=N^{1-\theta'}(\log N)^{2A+8}$, and
$\alpha<N^{1/2}$ since $\theta'>1/2$, so Bombieri–Vinogradov closes
them. The main term is
$A(N)\sum_{m\le\alpha}\mu(m)\lambda(m)T_m/m$, which is
$O\!\left(Ne^{-c\sqrt{\log N}}\right)$ by Lemma [lem:mu] applied
directly — the weight being $1$, no partial summation is needed to
remove a $\log m$, so this step is strictly shorter here than in
Section [sec:delta].

*The mean term.* Since $n>Y_k$ iff $k>R_n:=\lceil (N-n)/\alpha\rceil-1$,
exchanging gives
$\Xi_{\mathrm{mean}}=\sum_{n<N}a_n\bigl(\rho_N(K)-\rho_N(R_n+1)\bigr)$
with $\rho_N$ as in Lemma [lem:extract]. Its hypothesis is met at both
arguments used below — $K=N^{\theta'}$ and, on the bulk,
$R_n\ge N^{\theta'/2}$, against $n=N$ — so
$\rho_N(x)\ll e^{-c_1\sqrt{\log x}}$ applies there. Split at $N-n=H$ with
$H:=2N^{1-\theta'/2}$. On $N-n\le H$ use
$|\rho_N|\ll1$ and $\sum_{N-H\le n<N}|a_n|\ll H\log N$; since
$1-\theta'/2<1$ this is $\ll_A N(\log N)^{-A}$. On $N-n>H$ one has
$R_n\ge H/(2\alpha)=N^{\theta'/2}$, so both $\rho$ terms are
$\ll e^{-c_2\sqrt{\theta'\log N}}$, and $\sum_{n<N}|a_n|\ll N$ gives
$\ll N e^{-c_2\sqrt{\theta'\log N}}\ll_A N(\log N)^{-A}$.

Collecting, $T_1^\ast=J-\Xi_{\mathrm{main}}+\Xi_{\mathrm{mean}}
\ll_{A,\theta'}N(\log N)^{-A}$.
 ∎


#### Remark (what the proposition does not cover) {#rem:movingcutscope}

Proposition [prop:movingcut] is stated in the Corollary-1 regime only.
In the general Theorem-1 regime $\alpha\asymp N^{\theta}$ the outer
truncation is $K\asymp N^{1-\theta}$, which drops below $N^{1/2}$ as
soon as $\theta>1/2$, and Theorem [thm:A] is stated for
$K=N^{\theta'}$ with $\theta'\in(1/2,1)$. The bridge and the estimate
of $\Xi_{\mathrm{main}}$ are unaffected — the latter closes under the
$EH$ assumed there, as $\Delta$ does — but the fixed-cut input is not
available, so no claim is made outside $\theta'>1/2$. That is the same
threshold [thm:D] turns on, and it is not a coincidence: below it the
cofactor left by the completion is no longer short.



## Proof of Theorem [thm:C
: the permanent closure]{#sec:C}

Now take the weight $w_k=\log k$, i.e. the functional $E_3$ of
[eq:E3]. Steps 0–1 are unchanged, and the complete divisor sum of
Lemma [lem:complete] is replaced by the following.


#### Lemma {#lem:completelog}
<!-- evidence: analytic -->

For squarefree $u$, write $u=u_N u'$ where $u_N$ is the largest divisor
of $u$ composed of primes dividing $N$. Then
$\sum_{k\mid u,\ (k,N)=1}\mu(k)\log k=-\Lambda(u')$.


**Proof.** 
The sum is over $k\mid u'$ and equals $-\Lambda(u')$ by
$\mu*\log=\Lambda$.
 ∎


Hence the complete piece of $E_3$ is

$$
-\sum_{1\le u<N}\Lambda(N-u)\,\mu(u)\,\Lambda(u')
  \;=\; -\sum_{\substack{u<N,\ (u,N)=1}}\Lambda(N-u)\mu(u)\Lambda(u)
        + O(N^{o(1)}\log^2 N).
$$

Now $\mu(u)\Lambda(u)$ is supported on primes $u=p$, where it equals
$-\log p$; prime powers $u=p^\ell$, $\ell\ge2$, have $\mu(u)=0$. So the
complete piece equals

$$
\begin{equation}\label{eq:goldbachback}
  \sum_{p<N}\Lambda(N-p)\log p \;+\;O\bigl(N^{o(1)}\log^2N\bigr)
  \;=\;\sum_{n<N}\Lambda(n)\Lambda(N-n)+O\bigl(N^{1/2+\varepsilon}\bigr),
\end{equation}
$$

the binary Goldbach sum itself. The door is free: divisor switching
costs nothing, and immediately behind it stands the object one was
trying to avoid.

Two further terms have to be evaluated, and neither vanishes.

*(i) The residual main term.* Repeating Steps 3–6 with the extra
factor $\log k=\log((N-n)/m)$ produces the main term
$A(N)\sum_{m<M}\mu(m)\lambda(m)m^{-1}\int_{mK}^{N}\log(v/m)\,dv$.
Substituting $v=mt$ evaluates the integral as
$N\log(N/m)-N-mK\log K+mK$, of which only the $-N\log m$ piece survives
the cancellation: by Lemma [lem:mu] the coefficients $\sum f(m)/m$ and
$\sum f(m)$ both die, and what is left is
$-A(N)\,N\,\widetilde G(1)$ with
$\widetilde G(1)=\lim_x\sum_{m\le x}\mu(m)\lambda(m)\mathbf 1_{(m,N)=1}
\log m/m$. The constant is fixed by

$$
\begin{equation}\label{eq:Gtilde}
  A(N)\,\widetilde G(1) \;=\; -\SS(N),
\end{equation}
$$

so that the residual main term is $+\SS(N)\,N+O(N(\log N)^{-A})$.


#### Remark (the sign, and why it is the whole identity) {#rem:sign}
<!-- evidence: audit_E3_constant.py -->

Version 3 of this note printed [eq:Gtilde] as
$A(N)\widetilde G(1)=+\SS(N)$ and the residual main term as
$-\SS(N)N$. Since $T=P-\mathrm{MT}-C\,B$, that pair delivers
$E_3=\tilde r(N)+\SS(N)N+\SS(N)C(N)$, which differs from
Theorem [thm:C] by $2\SS(N)N$ — the size of the object itself. The
sign above is the one that makes the section agree with its own
theorem, and it is forced: with
$F(s)=\sum_m f(m)m^{-s}=\zeta(s)^{-1}H(s)$ as in Lemma [lem:mu] one
has $\widetilde G(1)=-F'(1)=-H(1)<0$, while
$A(N)H(1)=\SS(N)$ identically — the local factors pair as
$\bigl(1-\tfrac{1}{p(p-1)}\bigr)\bigl(1-\tfrac{1}{(p^2-p-1)(p-1)}\bigr)
=1-\tfrac1{(p-1)^2}$ at $p\nmid N$, and as $p/(p-1)$ at $p\mid N$.
Measured, $A(N)\widetilde G(x)=-1.760250$ at $x=4\cdot10^6$ against
$\SS(N)=1.760432$, and the brute-force $E_3$ sits within $0.26N$ of
$\tilde r-\SS(N-C)$ and $2\SS N$ away from the other candidate.


#### Remark (the identity cannot be tested numerically, at any $N$) {#rem:thetasweep}
<!-- evidence: lab_theta_sweep.py -->

That $0.26N$ is the finite-$N$ residual
$R(N,\theta')=\bigl|E_3(N;\theta')-(\tilde r(N)-\SS(N)(N-C(N)))\bigr|/N$,
and it decides only the sign of the identity because the right-hand
side it is compared against is $\asymp10^{-3}N$. Two ways out suggest
themselves and both are closed.

$\theta'$ does not help. Swept over $0.51$ to $0.95$, $R$ is
essentially monotone *increasing*: at $N=8\cdot10^5$ it runs
$0.170167$ at $\theta'=0.51$ to $4.991178$ at $\theta'=0.90$, a factor
$29$ worse. The finite-$N$ error is dominated by the main-term
cancellation over $m<M=N^{1-\theta'}$, which $\theta'\to1$ destroys,
and the gain in $B_{\log}(K)$ is second order against it. The best
$\theta'$ is the smallest admissible one, sitting just above the
$\tfrac12$ barrier. Pre-registered as rules I1, I2 and I3 — an
interior optimum at $\theta'\ge0.70$ giving $R<0.10$ — all three
fail, and they fail as a hypothesis rather than as a mis-stated
threshold. Rule I4 ties the sweep to the earlier audit and holds:
$R$ at $\theta'=0.56$ recomputes as
$0.4558,\,0.3729,\,0.3108$ against $0.456,\,0.373,\,0.311$.

$N$ does not help either — it makes matters worse. At the best
$\theta'$ the ratio of error to signal reads
$15.19,\,18.38,\,26.84$ at $N=2\cdot10^5,\,4\cdot10^5,\,8\cdot10^5$:
it *widens*. The signal falls by factors $0.598$ and $0.586$ per
doubling-and-a-bit while the error falls only by $0.723$ and $0.856$,
because $C(N)=o(N)$ being true kills the right-hand side fast while
the finite-$N$ error dies only at a slow power of $\log$. So there is
no $N$ at which a brute-force computation confirms
Theorem [thm:C]'s content, and waiting for larger $N$ points the
wrong way.

*(ii) The subtracted mean term.* It carries the factor

$$
B_{\log}(K):=\sum_{\substack{k<K\\ (k,N)=1}}\frac{\mu(k)\log k}{\varphi(k)}
  \;=\;-\SS(N)+O\bigl(e^{-c\sqrt{\log K}}\log K\bigr),
$$

which is exactly Huang–Li's Lemma 1: subtracting
$\sum_{k\le K}\mu(k)\varphi(k)^{-1}\log(K/k)=\SS(N)+O(e^{-c_1\sqrt{\log K}})$
from $\log K\cdot\sum_{k\le K}\mu(k)/\varphi(k)=O(\log K\,e^{-c_1\sqrt{\log K}})$
gives the claim. (Independently: with
$f(s)=\prod_{p\nmid N}\bigl(1-\frac{p^{-s}}{p-1}\bigr)$ one has
$f(s)=\zeta(s+1)^{-1}h(s)$ with $h(0)=\SS(N)$, so
$B_{\log}(\infty)=-f'(0)=-\SS(N)$.) Since the mean term is
$-C(N)B_{\log}(K)$ with $C(N)=\sum_{n<N}\Lambda(n)\mu(N-n)$, it
contributes $+\SS(N)\,C(N)$.

Collecting [eq:goldbachback], (i) and (ii):

$$
E_3(\alpha) = \sum_{n<N}\Lambda(n)\Lambda(N-n)
    -\SS(N)N+\SS(N)C(N)+O_A\!\left(\frac{N}{(\log N)^{A}}\right),
$$

which is the assertion of Theorem [thm:C]. Since $\SS(N)\asymp1$,
the bound $E_3(\alpha)\ll_A N(\log N)^{-A}$ is therefore equivalent to

$$
\tilde r(N)=\SS(N)\bigl(N-C(N)\bigr)+O_A\!\left(N(\log N)^{-A}\right),
$$

which is Huang–Li's equation (22). Given that, $\tilde r(N)$ is
positive for large even $N$ — binary Goldbach — because
$|C(N)|\le A(N)N(1+o(1))$ with $A(N)<1$; and
$\tilde r(N)\sim\SS(N)N$ holds if and only if $C(N)=o(N)$, which is not
implied and is the subject of the companion paper. \qed


#### Remark

The mechanism is structural, not accidental. The identity
$\mu*\log=\Lambda$ is the starting point of the Huang–Li argument
(their (10)); it is therefore also what any divisor switch inside their
$\log k$-weighted functional must return. No choice of $\theta'$,
truncation, or smoothing can circumvent an identity. What the identity
does not fix is the *strength* at which $E_3$ must be bounded, and
Proposition [prop:onesided] weakens that.


## The demand side is empty: a no-go over all weights {#sec:D}


Theorem [thm:A] and Theorem [thm:C] treat the two weights that
occur in [HL], $w_k = 1$ and $w_k = \log k$. It is natural to ask
whether some *other* weight sits between them: one whose complete
divisor sum is still cheap (as for $w=1$) but whose main-term
coefficient is still of size $1$ (as for $w=\log k$), and which would
therefore extract the scalar $C(N)=\sum_{n<N}\Lambda(n)\mu(N-n)$ from
Bombieri–Vinogradov alone. This section shows that no such weight
exists, and that the obstruction is quantitatively the $\sqrt N$
barrier.


### The design space


Let $w:\mathbb N\to\mathbb R$ be arbitrary and let
$b := \mu * w$, so that $w_k=\sum_{d\mid k}b_d$; every weight has this
form, and $b$ is determined by $w$. Put

$$
B_w := \sum_{\substack{k<K\\(k,N)=1}}\frac{\mu(k)w_k}{\varphi(k)},
  \qquad
  \|b\|_1 := \sum_{u<N}|b_u| .
$$

Running the switch of Section 3 with the weight $w$ gives the identity

$$
\begin{equation}\label{eq:extract}
  B_w\cdot C(N) \;=\; \underbrace{\sum_{u<N}\Lambda(N-u)\mu^2(u)\,b_u}
  _{\text{complete part}} \;-\; \underbrace{\mathcal R_w}_{\text{residual}}
  \;-\; \underbrace{T_w}_{\text{the object itself}}
  \;+\;O\!\left(N^{o(1)}\right),
\end{equation}
$$

the complete part being evaluated by the following lemma, which
identifies $b$ as exactly the complete divisor transform.


#### Remark (the third term is not a remainder, and it was missing) {#rem:extractTw}
<!-- evidence: analytic -->

$T_w$ is the object of this whole note — at $w=\log$ it *is*
$E_3(\alpha)$, which Theorem [thm:C] evaluates as
$\tilde r(N)-\SS(N)(N-C(N))+O_A(N(\log N)^{-A})$ — and it is of size
$N$, not $N^{o(1)}$. Earlier printings of [eq:extract] omitted it, and
that turns a finite rearrangement into a false statement. The
correction is algebra, not estimation: $T_w$ is defined as
$\sum_k\mu(k)w_k[A(N;k)-C(N)/\varphi(k)]$, so $B_wC(N)$ is
$\sum_k\mu(k)w_kA(N;k)$ minus $T_w$, and it is the first of those that
the switch and Lemma [lem:Gb] evaluate. A blind re-verification
computed the two sides separately and found the omitted form off by
$\asymp N$ — by exactly $|T_{\log}|$ at each $N$ — and the corrected
form exact to machine precision; the figures are recorded in
`verify/pass4/FINDINGS.md`.

The correction strengthens Theorem [thm:D] rather than weakening it,
but it changes what that theorem must say. [eq:extract] on its own
yields nothing about $C(N)$: an unconditional bound on $T_w$ is
precisely what the problem asks for, and the only weight for which one
is known is $w=1$, by Theorem [thm:A]. Theorem [thm:D] therefore
*grants* $T_w\ll_A N(\log N)^{-A}$ — the most $EH_\mu$ would give —
and shows the extraction fails anyway, on the separation of the two
thresholds alone.


#### Lemma {#lem:Gb}
<!-- evidence: analytic -->

For squarefree $u$,
$\sum_{k\mid u}\mu(k)w_k = \mu(u)\,b_u$.


**Proof.** 
$\sum_{k\mid u}\mu(k)\sum_{d\mid k}b_d
 = \sum_{d\mid u}b_d\,\mu(d)\sum_{j\mid (u/d)}\mu(j)
 = \sum_{d\mid u}b_d\,\mu(d)\,[\,u/d=1\,] = \mu(u)b_u$,
using $(d,u/d)=1$ for squarefree $u$.
 ∎


The two known cases are the two ends of this space: $w=1$ gives
$b=\delta_1$, $\|b\|_1 = 1$ and (Huang–Li's Lemma 1)
$B_w\ll e^{-c\sqrt{\log K}}$; while $w=\log$ gives $b=\Lambda$,
$\|b\|_1\asymp N$ — the complete part *is* the binary Goldbach
sum — and $B_w\to-\SS(N)\asymp1$.


#### Proposition (both ends are the same sum of dilated walls) {#prop:flatsum}
<!-- evidence: lab_weight_gap.py -->

Put [eq:dilate] into $T_w$. The transform carries $\mu(k)$, so terms
with $\mu(k)=0$ contribute nothing on the left; on the remaining $k$,
$\mu(k)A(N;k)=\mu(k)^2H(N;k)$. The restriction must therefore be
carried on the right, and

$$
\begin{equation}\label{eq:flatsum}
  T_w(N) \;=\; \sum_{\substack{k<K\\(k,N)=1}} \mu^2(k)\,w_k\,H(N;k)
           \;-\; B_w\,C(N) .
\end{equation}
$$

Dropping $\mu^2(k)$ does not give a weaker identity but a false one:
$H(N;k)$ does not vanish at squarefull $k$, and the terms it would add
are of the size of $T_w$ itself.

The two known cases are therefore

$$
T_1(N) = \!\!\sum_{\substack{k<K\\(k,N)=1}}\!\! \mu^2(k)H(N;k)
  \;-\; B_1 C(N),
\qquad
E_3(N) = \!\!\sum_{\substack{k<K\\(k,N)=1}}\!\! \mu^2(k)(\log k)H(N;k)
  \;-\; B_{\log}C(N),
$$

and both weights are nonnegative on the range. Since
$B_1\ll e^{-c\sqrt{\log K}}$ kills the second term while
$B_{\log}\to-\SS(N)$ does not, Theorem [thm:A] — unconditional —
says exactly that the *flat* sum of dilated walls is
$\ll_A N(\log N)^{-A}$, and the wall says the same sum weighted by
$\log k$ is not small. Everything between what is proved and what is
open is the factor $\log k$ inside a positively weighted sum of the
same terms.

The identity is not the point; the sizes are. Over
$N=2\cdot10^5$ to $2.56\cdot10^7$ by doubling at $\theta'=0.56$,
[eq:flatsum] holds to a worst relative error of
$1.875\cdot10^{-16}$, and

$$
\frac{\bigl|\sum_k \mu^2(k)H(N;k)\bigr|}
     {\bigl|\sum_k \mu^2(k)(\log k)H(N;k)\bigr|}
  \;=\; 0.1785,\ 0.1559,\ 0.1484,\ 0.1483,\ 0.1367,\ 0.1334,\
        0.1241,\ 0.1135 ,
$$

falling throughout; at the top $|T_1|/N=0.01425$ against
$|E_3|/N=0.1245$. Fitted against $\log N$ the flat sum decays as
$N^{-0.3505}$ and $|E_3|/N$ as $N^{-0.2658}$: the gap is widening, as
it must if one side is $\ll_A N(\log N)^{-A}$ and the other is
$\asymp N$.

The index set here runs from $k=1$, as [eq:flatsum] does. Dropping that
one term changes the ratios to
$0.1807,\ 0.1740,\ 0.1624,\ 0.1389,\ 0.1456,\ 0.1216,\ 0.1258,\ 0.1188$
and the exponent to $N^{-0.3620}$, and it is the second of these that
`lab\_weight\_gap.py` computes — its loop starts at $k=2$. The
difference is exactly $C(N)$, since $H(N;1)=C(N)$, so it moves the flat
sum and not $T_1$, whose $k=1$ term is $A(N;1)-C(N)/\varphi(1)=0$. The
figures printed above are from
`verify/pass5/code/c02\_flatsum\_k1.py`, which
recomputes both index sets from an independent sieve.


#### Remark (the last declined null, rule by rule) {#rem:identitynull}
<!-- evidence: audit_directidentity_null.py -->

Remark [rem:whereitlives]'s evidence declines a control by argument
rather than by pointer: Z1 is an identity, Z2 compares against
$\SS(N)N$ which is the reference, and Z3 and Z4 "compare two sums over
the same terms, so a sign control would move both sides together".
Remark [rem:weightgapnull] fixed how to judge such an argument — a
control is worth running exactly when the statistic stays well
conditioned under it — and applied rule by rule it is right in two
places and wrong in one.

Z3 and Z4 divide by the *total*, and for a coin the total is not a
fixed reference, so those ratios are ill conditioned and the decline
is correct. Z2 divides by $\SS(N)N$, a fixed nonzero constant: that is
well conditioned and the control was owed. Run, it is decisive.
$T/(\SS(N)N)$ is $1.0039,\,0.9865,\,0.9893,\,1.0017,\,0.9968$ for
$\mu$ and $4.9835$ to $6.2827$ across eight draws — no draw within a
factor of five of the reference.

**The pre-registered V1 fails, and its failure is the finding.** It
predicted the coin's untruncated sum would be a small fluctuation,
under $0.2$ of the count. It is $4.9641$ to $6.3027$ times the count,
and *rising* with $N$. So [eq:untrunc] is not a small correction that
$\mu$ happens to satisfy: without the Möbius identity the same double
sum is six times the answer and growing, and $\mu$ collapses it to the
answer exactly — $|T-R|/R$ of $3.294\cdot10^{-16}$ to
$3.136\cdot10^{-15}$.

Z3 admits a repair the criterion suggests: measured against $N$ rather
than against the total, its partial sums are well conditioned. So
measured, $\mu$'s partial over $k<N^{0.90}$ is
$-4.3412,\,-4.6951,\,-5.0055,\,-5.3075,\,-5.6291$ against draws that
run $+2.1498$ to $+2.4747$ — larger in magnitude than every draw at
every $N$, and negative throughout where $0$ of $8$ draws are. The
one-signed mass Remark [rem:whereitlives] found is $\mu$'s, not a
property of the truncation.

This was the last of the six results that declined a control. Four
declined by pointing at another script and all four pointers missed;
this one declined by argument and the argument was right for two rules
of four.


#### Remark (the declined null, run) {#rem:weightgapnull}
<!-- evidence: audit_weightgap_null.py -->

The evidence for the next two statements declines a control, on the
ground that its claims compare two weightings of the same numbers so a
sign control would move both sides alike. Remark [rem:splitnull]
refuted a reason of exactly that shape, so it is run here: the coin
$\varepsilon(v)=\pm1$ on $\operatorname{supp}\mu^2$, field, weights,
$k$-range and truncation identical, eight draws.

**It was worth running, and it splits the claims in two.** The
identity survives it — the worst relative error of
$T_1=\sum H-B_1C(N)$ over $\mu$ and all draws is
$3.198\cdot10^{-14}$, as algebra requires. And the *decay ordering*
is $\mu$'s: over $N=2\cdot10^5$ to $3.2\cdot10^6$ the flat sum decays
as $N^{-0.3698}$ against $|E_3|$'s $N^{-0.2713}$, and only $3$ of $8$
coins reproduce that ordering, their exponents scattering from
$-1.1977$ to $1.3741$.

The other two do not survive, and not because they are wrong.
**V1** asked $\mu$'s ratio $|\sum H|/|\sum(\log k)H|$ to be below every
coin's; it is not — the coin band is $[0.0735,\,2.4871]$,
$[0.0000,\,0.2564]$, $[0.0422,\,8.0955]$, $[0.0631,\,0.2236]$,
$[0.0285,\,1.2501]$, and $\mu$'s $0.1807$ to $0.1456$ sits inside it.
**V2** asked the same of the profile's spread and gets the same answer,
the coin reaching $1.0017$ against $\mu$'s $1.0146$ at the smallest
$N$ and $2217.7429$ at the largest excursion. The reason is visible in
the bands: a coin drives *both* sums to square-root size, so their
ratio is a quotient of two near-zero quantities and takes any value it
likes. **The coin is not a usable reference for a ratio of this kind.**

So the original refusal was right in outcome and wrong in reason. It
is not that a coin moves both sides alike; it is that a coin makes
both sides small and the ratio ill-conditioned. What follows is that
the ratio of Remark [rem:weightgap] and the effective modulus read off
its profile stand **uncontrolled** — not refuted, and not supported by
a control either. The decay ordering is the part of that remark a
control does reach, and it holds.


#### Remark (the heuristic, corrected) {#rem:heuristic}
<!-- evidence: audit_directlevel_heuristic.py -->

Remark [rem:directlevel] reports that the measured crossing is not
only of the right exponent but of the right size, the naive solution
of $K=\SS^2N/(4\log^2K)$ missing by $23\%$ at the bottom and $4\%$ at
the top and "closing as $N$ grows". Two later measurements say that
heuristic dropped two large factors. The constant in the square-root
law is not $1$: measured below $K^*$, $\overline{|H|/\sqrt{N/k}}$ is
$3.5421,\,3.5393,\,3.7012,\,3.7565,\,3.8100$, which divided by
$\sqrt{\log N}$ is $1.0138,\,0.9854,\,1.0039,\,0.9939,\,0.9844$ —
**$c(N)=\sqrt{\log N}$ to within a percent and a half**. And the $k$
are squarefree and coprime to $N$, of density $0.3374$ to $0.3378$,
where the heuristic integrated over every $k$.

The two nearly cancel: their product is
$1.1950,\,1.1950,\,1.2496,\,1.2682,\,1.2869$. Reconstructed here, the
naive prediction reproduces the published ratios exactly —
$1.2271,\,1.2020,\,1.0458,\,0.9900,\,0.9551$ — and those are not
converging, they are **drifting monotonically** through $1$ near
$N=1.6\cdot10^6$. The $4\%$ at the top is where the drift happens to
cross, not a limit.

**Putting both factors back makes the prediction far better, and the
pre-registered X3 fails at five of five $N$** — the outcome flagged in
advance as the good one. Solving
$\sum_{k<K}(\log k)\,c(N)\sqrt{N/k}=\SS(N)N$ over the *actual*
admissible $k$ gives $2981,\,5057,\,8061,\,13589,\,23059$ against the
measured $2973,\,5109,\,8021,\,13557,\,23397$ — ratios
$0.9973,\,1.0103,\,0.9950,\,0.9976,\,1.0147$, **within $1.5\%$ at
every $N$ and with no drift at all**.

So the square-root law is confirmed far more sharply than Remark
[rem:directlevel] claimed, and its stated agreement should be read as
a coincidence of two omissions rather than as evidence. What the
corrected form adds is the constant: the crossing is governed by
$|H|\approx\sqrt{\log N}\sqrt{N/k}$, so
$K^*\asymp\SS^2N/(4\,c^2\log^2K^*)$ with $c^2=\log N$, and
$K^*/\sqrt N$ still grows — but as $\sqrt N/\log^3$ rather than
$\sqrt N/\log^2$. Remark [rem:modeltransfer] asks whether that
corrected form is a model or a fit.


#### Remark (a model of the exponent, not of the constant) {#rem:modeltransfer}
<!-- evidence: audit_model_transfer.py -->

A prediction that reproduces the crossing it was calibrated on to
$1.5\%$ has shown nothing yet. The repository holds a second crossing
of the *same* sum $B_H(N;K)=\sum_{k<K}\mu^2(k)(\log k)|H(N;k)|$ against a
different budget: [eq:nolog] asks it against $\SS(N)(1-A(N))N$, a
factor $4.70$ smaller, and that crossing sits an order of magnitude
lower. Calibrating $c(N)$ **only below the $\SS N$ crossing** and then
predicting both is asking the model something it cannot have fitted.

It half survives. Pre-registered V2 holds — the in-sample crossing is
reproduced at $1.0027,\,0.9898,\,1.0050,\,1.0024,\,0.9856$, confirming
[rem:heuristic]. V3 holds: the calibration is not local, $c$ taken
below the small crossing agreeing with $c$ below the large one to
$5\%$ ($1.0157,\,0.9514,\,1.0136,\,0.9912,\,0.9467$). **But V1 fails
and V4 with it.** Out of sample the model predicts
$323,\,521,\,791,\,1283,\,2099$ against the measured
$321,\,579,\,789,\,1311,\,2357$ — ratios
$1.0062,\,0.8998,\,1.0025,\,0.9786,\,0.8905$, missing by $10\%$ at two
of the five $N$, and getting the *ratio* of the two crossings right at
only three of five.

The failure is not random and it is not small. **The two $N$ at which
the model misses are exactly the two at which $c$ falls by about $5\%$
between the two $k$-ranges.** Since $B_H\sim cK^{1/2}\log K$, an error
in $c$ must be amplified into the crossing by $1/(\tfrac12+1/\log K)$,
about $1.56$ here; over the three rows where $c$ drifts *down* the
measured amplification is $2.1814$ against that $1.5587$, some $40\%$
stronger. Over the rows where $c$ drifts *up* it is $0.2920$: the
model absorbs an upward drift and not a downward one. That asymmetry
is unexplained and this audit does not settle it.

What survives is worth stating precisely, because it is what the
$\theta'$ target actually needs. A budget factor of $4.7009$ moves the
crossing by $0.1677$ in the exponent — measured gaps
$0.1824,\,0.1688,\,0.1706,\,0.1635,\,0.1532$ between
$\log K^*(\SS N)/\log N$ and $\log K^*(\SS(1-A)N)/\log N$ — and the
single calibration reproduces that shift without being told it. The
two constants are $\SS(N)=1.760432$ and $\SS(N)(1-A(N))=0.374487$ over
this family, and this is the measurement that fixes what the
difference between them is worth; every $K^*$ exponent in either paper
belongs to one of the two and to no other.

**And this gap is measured at one odd radical.** Every $N$ here is
$2^a5^b$. Remark [rem:residuearithmetic] measures the same response
across seven arithmetic types by regressing the exponent on the
logarithm of the threshold, and gets a slope of $+0.0516$ against the
$+0.1084$ that $0.1677$ per factor $4.7009$ implies — **half as
strong**. So the $0.1677$ is the right correction *within this
family*, and applying it to another arithmetic would overstate the
cost of a smaller budget by about a factor of two. That is the one
place in the whole chain where the arithmetic works in the
programme's favour.
**So $|H|\approx c(N)\sqrt{N/k}$ is a model of the exponent and only a
fit of the constant.** A $5\%$ wobble in $c$ costs $10\%$ in $K^*$,
which is fatal to any argument that has to decide a crossing to better
than a factor, and harmless to one that only has to place a level.
Remark [rem:directlevel]'s $K^*/\sqrt N\to\infty$ is of the second
kind and stands; the forecasts in [rem:forecast] and the bracket in
[rem:extendrange] are of the first kind and inherit a $10\%$ error bar
on $K^*$ that neither of them carries. Remark [rem:forecastbracket]
pays that debt for the first: propagated, the $10\%$ becomes
two-thirds of a decade on the forecast $N$, because the square root of
$K$ doubles the exponent and the logarithms amplify rather than damp.


#### Remark (the binding half is square-root too) {#rem:elemsize}
<!-- evidence: lab_elementary_size.py -->

Remark [rem:splitbudget] moved the target to the elementary half, so
the first thing to know about $P$ is its size. Fitting the octave
means of $|P(N;k)|$ against $N/k$ — abscissa the mean of $N/k$ inside
each octave, as Remark [rem:residue] fixed; bins closed at both ends
and required to hold at least ten $k$, as Remark [rem:elemreach]
forced — gives exponents

$$
0.5793,\ 0.5147,\ 0.4972,\ 0.5012,\ 0.4678
$$

with correlations $0.99810$ to $0.99989$ and leave-one-out spreads
$0.0031$ to $0.0210$. Against the residue's $0.4869$, read from its
own result file, the gap at the largest $N$ is $0.0191$. The
population floor is a threshold and is swept: over $5$, $10$ and $20$
$k$ per octave the exponent moves by at most $0.0325$ (rule Y5, added
after the correction and disclosed as such). The thinnest octave any
of these five fits stands on holds $25,\,13,\,25,\,13,\,25$ values of
$k$ — a fit is only as good as its emptiest bin, and that number, like
the correlations above, is now declared to the gate rather than left
to be discovered.

**These are not the numbers this remark first carried, and the old
ones are not reprinted here because nothing computes them any more.**
Version one fitted through an unbounded top bin $[32768,\infty)$. Its
five exponents fell monotonically across the sweep, ending well below
$\tfrac12$; its correlations were an order of magnitude looser than
those above and its leave-one-out spreads several times wider. Remark
[rem:elemreach] predicted that the drift was the thin end of the fit
and not the object. Closing the bin and requiring it to be populated
removes the drift, tightens every correlation, and leaves the
exponents sitting on $\tfrac12$ — which is what the prediction said
would happen.

**So every part of $H$ is square-root in $N/k$: $H$ itself, the
elementary $\beta P$, and the residue $R$.** The split of Remark
[rem:predictable] buys a constant — the mass shares — and not an
exponent, and Remark [rem:splitbudget]'s finding that the elementary
half binds is a statement about that constant. What binds is
therefore hard in the same shape as what does not: square-root
cancellation in a Möbius sum, which is the thing that is unproved
either way. Remark [rem:elemreach] tries to push that exponent
further out than $H$ can be measured, fails, and finds the drift in
the numbers above to be an artefact of the thin end of the fit.


#### Remark (the reach that cannot be bought) {#rem:elemreach}
<!-- evidence: lab_elementary_reach.py -->

Every measurement of $H$ here stops at $N=3.2\cdot10^6$ because $H$
needs $\Lambda(N-mk)$ and so a sieve to $N$. $P$ needs no primes. The
sieve weight factors as $w(m,k)=C_k\cdot\mathbf1[m\not\equiv
Nk^{-1}\ (q)\ \forall q]$ with $C_k=\prod_{q\nmid k}q/(q-1)$ over the
odd $q\le29$, so $P$ is one constant times an integer sum of $\mu$
over an explicit sifted set — a boolean mask, no $\Lambda$ anywhere.
That reaches $N=10^8$ and inner lengths $3.3\cdot10^7$ against the
$1.6\cdot10^6$ above. The control Y1 confirms it is the same $P$:
$\beta=\sum HP/\sum P^2$ recomputed this way reproduces the published
cross-check at all five $N$.

**The reach does not buy the measurement, and the reason is
structural.** At fixed $N$ the longest inner lengths come from the
fewest $k$; inverting the sampling to a ladder of $N$ with $k<400$
each does not repair it, because the number of pairs landing at inner
length $L$ falls like $1/L$. The pooled octave counts are
$598,\,668,\,583,\,190,\,48,\,8$: **every factor of four further out
costs a factor of four in sample count.** Declared as gate check G30
now requires, the two fits here stand on thinnest bins of $3$ and $8$
pairs — which is the whole finding in two numbers.

**And their correlations are $0.99724$ and $0.98087$.** A fit through
bins of three pairs is as tight as anything in this paper and means
nothing at all, which is why the two declarations are needed together:
correlation measures whether the points lie on a line, and says
nothing about whether the points are worth anything. Neither number
alone would have caught this; the pair does. Pre-registered Y5 asked for
$200$ pairs in every octave and a coin band under $0.10$ wide and
**is refuted** at both, and Y2, Y3 and Y6 are refuted with it —
exponents $0.3622$ and $0.3674$, far below square-root. Those
refutations are not evidence about $\mu$: the coin arm, on the same
sifted set and the same bins, scatters over $[0.4244,\,0.5294]$, a
band four times wider than the question. This is the conditioning
criterion of Remark [rem:weightgapnull] applied to our own statistic.

What the data do support is sharp. Keeping only the octaves with at
least $200$ pairs — three of them, inner lengths $18265$ to
$265770$ — gives

$$
\text{$\mu$: } 0.5178, \qquad
\text{coins: } [0.4844,\ 0.5095],\ \text{width } 0.0251,
$$

so **$|P|$ is square-root and if anything marginally above it**, not
below. The same trade shows in the two fixed-$N$ sweeps: where both
reach the same inner length their octave means agree to
$0.985,\,0.928$ while the bins hold thousands of pairs, and diverge to
$2.201$ where one holds thirteen.

**This corrected the reading of the exponents in Remark
[rem:elemsize].** As that remark first stood, its exponents drifted
monotonically downwards across the sweep, in the same direction as
the artefact isolated here, and its top bin was the unbounded
$[32768,\infty)$ — precisely the thin end. The prediction made here
was that the drift belonged to the fit and not to the object. Closing
the bins and requiring ten $k$ in each confirmed it: the drift is
gone and the exponents sit on $\tfrac12$.

Two rules still fail, and one of them says something the first
version did not. **Y2** asked $\mu$ to match a coin on the same sum
and pinned the coins themselves to $[0.40,\,0.60]$; their exponents
run $0.3582$ to $0.6031$, which is the noise scale of a six-point
octave fit repeated eight times, so the band was too tight for the
control rather than wrong about it. $\mu$'s exponent lies inside the
coins' observed range at every $N$. **But it lies above their median
at all five**, by $0.1111,\,0.0576,\,0.0458,\,0.0629,\,0.0269$ — five
of five on one side, where the coins themselves are symmetric by
construction. Remark [rem:elemreach] found the same sign
independently on a longer lever, $0.5178$ against a coin band topping
out at $0.5095$. Two measurements, different ranges and different
estimators, both put $|P|$ *above* what random signs give. The offset
is inside the per-$N$ coin spread each time and so is not resolved at
any single $N$; what is not accidental is its sign.

Remark [rem:muvscoin] settles that sign with an instrument that does
not go through a fit, and the answer is far larger than either
estimate suggested.

**Y3** asked the sieve weight to inflate uniformly and it does not.
Against the unweighted Möbius sum $F$ over the same range,
$\overline{|P|}/\overline{|F|}$ climbs from $0.9364$ at
$N/k\in[32,128)$ to $1.9632$ at $[32768,131072)$ — a spread of
$2.2338$ at its worst across the octaves fitted. The weight is nearly
inert where the inner sum is short and doubles the sum where it is
long, which is the opposite end from where Remark [rem:leanodd] found
the sign structure. Past the population floor the ratio runs away
entirely — $3.1071$ and $32.2243$ in the two bins holding a handful of
$k$ — which is the same thin-end artefact seen from another angle.


#### Remark (the elementary sum beats a coin, by a constant) {#rem:muvscoin}
<!-- evidence: lab_mu_vs_coin_size.py -->

Two measurements had put $|P|$ above what random signs give and
neither could resolve it: the offset sat inside the per-$N$ coin
spread every time, and all that was not accidental was its sign. Both
went through a fitted exponent, which is the wrong instrument. An
exponent is a slope through six points carrying the noise of all six;
the question is a ratio of two magnitudes at the *same* inner length,
which needs no fit. Since the sieve weight is $C_k$ times an
indicator, $C_k$ cancels from that ratio and what is compared is
$\sum_{m\in S}\mu(m)$ against $\sum_{m\in S}\varepsilon(m)$ on the
identical sifted set.

One trap has to be avoided first, and this remark fell into it once.
$C_k$ cancels from the ratio **at a fixed $k$**; it does not cancel
from a ratio of *means over $k$*, because $C_k$ depends on which small
primes divide $k$. The budget carries $C_k$, so the means must, and
the figures below are the corrected ones. Remark [rem:provablehalf]
caught the discrepancy by measuring the same constant a second way.

**The effect is not marginal.** Against $32$ global sign vectors —
global, so the coins carry the same across-$k$ correlation $\mu$ does
— the ratio of $\mu$'s octave mean to the coins' median runs $1.3$ to
$1.8$ across the range, $\mu$ ranks $25$ to $32$ of $32$ in almost
every cell, and it is above the coin median in $29$ of $30$. Rule Z4
asked $\mu$ to stay inside the coins' range and **is refuted**: at one
cell $\mu$ beats all thirty-two draws. The control Z1 reproduces every
published octave mean to $10^{-4}$.

**And it is a constant, not an exponent.** Pooled over the five $N$,
the ratio by octave reads

$$
0.8877,\ 1.2454,\ 1.4487,\ 1.3969,\ 1.4879,\ 1.5586,\ 1.4634,\ 1.3172
$$

from $N/k\in[2,8)$ upwards — it climbs out of the short octaves and
then flattens. Fitting $\log(\text{ratio})$ on $\log(N/k)$ from
$N/k=128$ up, where the transient is over, gives a slope of
$-0.0097$ with a leave-one-out spread of $0.0443$: **indistinguishable
from zero, around a mean ratio of $1.4621$.** Rule Z3 now holds at
three of five $N$, but its instrument is still the wrong one and the
pooled profile says why — each $N$ starts and ends at a different
octave, so differencing the last against the first compares different
parts of a non-monotone curve at different $N$.

Three things follow. The exponent measurements of Remarks
[rem:elemsize] and [rem:elemreach] stand: $|P|$ is square-root, and so
is a coin, and the factor between them does not touch the exponent.
The *constant* does not: any heuristic that prices
$\sum_{k<K}(\log k)|P|$ at what square-root cancellation with a random
sign pattern would give understates it by about $1.46$, and Remark
[rem:heuristic]'s calibration $c(N)=\sqrt{\log N}$ is a measurement of
$|H|$ precisely because no such pricing is available a priori. And the
one cell where $\mu$ falls *below* every coin is the shortest octave,
$N/k\in[2,8)$ at ratio $0.8877$ — where parity leaves the
non-negative $m=1$ term nearly alone, the structure Remark
[rem:leanodd] identified.


#### Remark (size and provability point at opposite halves) {#rem:provablehalf}
<!-- evidence: lab_elementary_provable.py -->

Remark [rem:splitbudget] puts the target on the elementary half:
$\beta B_P$ takes $0.81$ to $0.87$ of the budget and $B_R$ about half,
so **the wall is the sieve-weighted Möbius sum and not the residue.**
That is a statement about size at accessible $N$. About what can be
*proved* the two halves are not symmetric at all, and the asymmetry
runs the other way.

Write $S$ out. $P(N;k)=C_k\sum_{m\in S}\mu(m)$ with

$$
S=\{\,m<N/k:\ m \text{ odd},\ \mu^2(m)=1,\ (m,k)=1,\
m\not\equiv Nk^{-1}\ (q)\ \forall q\le29\,\}.
$$

Every condition here is multiplicative or a residue condition to a
**bounded** modulus. Coprimality to $k$ is not a residue condition
modulo $k$ — in the Dirichlet series it deletes Euler factors and
contributes $L(k)=\prod_{p\mid k}(1-1/p)^{-1}$, nothing worse. So the
classical unconditional $|\sum_{m\le x}\mu(m)|\le Ax\exp(-c\sqrt{\log
x})$ is of the *right shape* for $P$ uniformly in $k$, whereas for
$R$ — a genuine Möbius–prime correlation of length $N/k$ at level $k$
— nothing of the kind is available past $k=N^{1/2}$. **If the
uniformity holds then $\sum_{k<K}(\log k)|P|=o(N)$ for every fixed
$\theta'<1$, the elementary half is asymptotically free, and the
entire obstruction is the half [rem:splitbudget] measured as the
smaller one.**

This does not prove the uniformity and does not claim to. What it does
is test the shape where the shape can be falsified, and then price it.

The shape survives. Rule W1 asked the ratio $|P|/[(N/k)\exp(-c\sqrt{\log
(N/k)})L(k)]$ at $c=0.2098$ to stay below $1$ everywhere and **is
refuted** — but the maximum is attained at $N/k=7,13,30,63,108$, which
is $N$ over the $k$-cap, and an asymptotic estimate has no content at
$x=7$. Restricted to $N/k\ge2,8,32,128$ the maxima are
$1.2119,\,1.0710,\,0.7309,\,0.3363$: the violation lives entirely
below $x=32$ and the constant forced anywhere an estimate could speak
is under $1$. W2 holds at five of five — the ratio falls monotonically
across the octaves, from $7.8$ to $174.9$ in looseness at the largest
$N$, because the truth is $(N/k)^{1/2}$ and the bound is $N/k$ damped
sub-logarithmically, so the gap is itself a power.

The null separates the two causes. Coins on the identical sifted set,
measured against the identical bound, sit at
$0.10836,\,0.08210,\,0.04182,\,0.02266,\,0.01152,\,0.00546$ across the
octaves at the largest $N$ against $\mu$'s
$0.12883,\,0.08699,\,0.04734,\,0.02639,\,0.01356,\,0.00572$: **a coin
falls away from the bound at the same rate**, so the looseness is the
shape and not $\mu$, and what $\mu$ adds is the constant of Remark
[rem:muvscoin]. Measuring that constant here a second way is what
caught the $C_k$ weighting error in it.

**And the price is the whole answer.** At accessible $N$ the bound
alone would spend $13.98,\,15.38,\,16.82,\,18.29,\,19.83$ times the
budget at $\theta'=0.56$ — not merely useless but *getting worse*
across the sweep, since $\exp(-c\sqrt{\log x})$ falls more slowly than
$(\log K)^2$ rises. Solving
$A\,d_L\int_0^{\theta'u}v\,e^{-c\sqrt{u-v}}\,dv=\SS(N)(1-A(N))$ with
$u=\log N$, the measured $d_L=0.3994$ and the conservative
$A=1.2119$ forced by the data gives

$$
N\approx10^{5475},\qquad
[\,10^{2093},\ 10^{13093}\,]\ \text{over } c\in[0.15,0.30].
$$

That bracket sweeps $c$ and nothing else, which is worth saying now
that gate check G33 asks what a bracket covers. The forecast also
rests on the measured $d_L$ and on the implied constant $A$. The first
does not drift at all here — every $N$ in the sweep has the same odd
radical, so the admissible $k$-set and $d_L$ are identical across it,
a fact about this sweep and not a general one. The second is not a
drifting constant either: $A$ falls monotonically because the maximum
of the ratio is always at the shortest inner sum, $N$ over the
$k$-cap, which grows with $N$ — it is a function of where one looks,
and the *largest* value is the one used, which is the conservative
choice for an upper bound. Scaling $A\,d_L$ by $2.1071$ either way
moves the forecast over $[10^{4839},\,10^{6140}]$ against the
$c$-sweep's $[10^{2093},\,10^{13093}]$. Both matter and neither is
within nine thousand orders of magnitude of anything computable.

So the reduction is real in shape and empty in size. What it changes
is where the program should push: **not at $P$.** Any effort spent
proving a bound for the elementary half is spent on the half that
already has one in principle; the half with no bound at all past
$k=N^{1/2}$ is $R$, and [rem:splitbudget]'s ranking is a fact about
finite $N$ that inverts in the limit. What [rem:residue] measures —
$|R|\asymp(N/k)^{1/2}$ with the lean removed — is therefore the
statement the whole route turns on, and it is smaller than $H$ by a
constant and no easier by an exponent.


#### Remark (what the conditional reduction actually buys) {#rem:residuelevel}
<!-- evidence: audit_residue_level.py -->

Grant Remark [rem:provablehalf]'s uniformity. Then
$|H|\le\beta|P|+|R|$ gives $B(N)\le B_R(N)+o(N)$ and the route's
condition is a condition on the residue alone, permitting the level
$K^*_R$. Remark [rem:splitbudget] prints one — $9191$ to $63399$, at
exponents $0.7477$ down to $0.7382$ — and $0.74$ would clear
$\theta'=0.56$ with room to spare.

**It is against the wrong budget.** That table crosses each half
against $\SS(N)N$; Proposition [prop:nolog] needs
$\SS(N)(1-A(N))N$, smaller by $4.7009$. Recomputed here — the control
reproduces all five published $K^*_R$ to $1.0000$ — the operative
crossings are $993,\,1447,\,2019,\,3319,\,5923$ and the exponents are

$$
0.5654,\ 0.5642,\ 0.5599,\ 0.5675,\ 0.5799 .
$$

The gap between the two budgets is $0.1823,\,0.1775,\,0.1794,\,0.1696,\,
0.1583$, matching the $0.1677$ that Remark [rem:modeltransfer]
measured for $K^*_H$ — so the effect is the budget and not the half.

**This is the sharpest statement of where the program stands, and it
is a knife-edge.** Rule U2 holds: the residue alone carries the level
past the square-root barrier at every $N$, by a margin of $0.06$ to
$0.08$. Rule U3 asked it to clear $\theta'=0.56$ as well and **is
refuted** — at $N=8\cdot10^5$ the exponent is $0.5599$, short by one
part in five thousand. Rule U4 holds under its own registered rule —
the least-squares slope of the exponent against $\log N$ is
$+0.004692$, which is not negative — **but its reading is withdrawn.**
That slope reaches $1.60$ standard errors and its two-sigma interval
contains zero, and the leave-one-out agreement offered for it is a
property of least squares rather than of the data (Remark
[rem:slopes]). Whether the margin is closing is not determined *here*
— **it is determined two octaves further out**, where the same slope
over seven points reaches $3.71$ standard errors and is positive
(Remark [rem:slopereach]). The conclusion U4 drew is right; the five
points it drew it from could not support it.

What this does and does not say. It does not prove anything: the
uniformity it is conditional on is unproved, and an exponent measured
over a factor $16$ in $N$ forecasts nothing (Remark
[rem:forecastbracket]). **And every $N$ here is $2^a5^b$: the sweep
has one odd radical, which its own result file now declares.** Remark
[rem:residuearithmetic] repeats the measurement across arithmetic
types and finds the barrier not cleared at the primorial-like ones, so
the figures below are about this family and not about even numbers. What it does is price the reduction. **The
elementary half, which spends most of the budget, is the half with a
classical estimate in the right shape; strip it and the half that is
left carries the level to $0.56$–$0.58$, not to $0.74$ and not to
$0.50$.** The programme's target $\theta'>1/2$ is cleared by the
conditional reduction, and cleared by less than a tenth.


#### Remark (what the signs cost once the elementary half is gone) {#rem:residuesigned}
<!-- evidence: lab_residue_signed.py -->

[eq:direct] carries the signed sum over $k$; [eq:directcond] discards
the signs. Remark [rem:signedlevel] priced that discard for $H$ at a
factor $1.793$–$2.223$ in $K^*$, about $0.053$ in $\theta'$ — a real
gain and a bounded one. The same question for $R$, which is what the
conditional reduction of Remark [rem:provablehalf] leaves behind, had
not been asked.

**The discard costs four to six times as much there.** Against the
operative budget, and with the control S1 reproducing the absolute
crossings $993,\,1447,\,2019,\,3319,\,5923$ exactly, the signed walk
$\sum_{k<K}(\log k)R(N;k)$ does not leave
$[-\SS(1-A)N,\,+\SS(1-A)N]$ until

$$
K=35597,\ 37623,\ 48957,\ 68669,\ \text{and not at all below }10^5,
$$

factors of $35.8,\,26.0,\,24.2,\,20.7$ and more, at exponents
$0.8586,\,0.8167,\,0.7945,\,0.7796$ — **a gain of
$+0.2932,\,+0.2526,\,+0.2346,\,+0.2121$ in $\theta'$**, where $H$'s
was $0.053$. Rules S2 and S4 hold at every $N$: the signed level
clears $\theta'=0.56$ with room the absolute level never had.

**Rule S3 is refuted, and its refutation is the caution.** It asked
$R$'s signs to be random-like — not worse than a redraw. They are not:
none of the sixteen draws that hold every $|R(N;k)|$ fixed and redraw
its sign crosses anywhere below $k=10^5$, while $\mu$'s walk crosses
at four of the five $N$. So the residual lean survives the split.
Remark [rem:residue] located the lean in the elementary half by mass
fraction — $P$'s $f_+$ an order below the sign band — but $R$'s own
$f_+$ of $0.5516$ to $0.4832$ against a band of width about $0.02$ was
already outside it, and that is what this measures at the level.
**$\mu$ is still worse than random signs; it is merely much less worse
than before the split.**

What this changes. It does not open a route: a bound on $|R(N;k)|$ is
what an estimate supplies, and the signed sum is not that. What it
does is size the single largest loss in the chain. Of everything the
programme discards, **the signs across $k$ in the residue are worth
more than the split itself** — Remark [rem:splitbudget] measured the
whole elementary/residue division at about $0.06$ in $\theta'$, and
this one step at $0.21$ to $0.29$. Any future estimate that could
retain even part of the cancellation across $k$ would buy more than
every other refinement in these papers put together.

*Added later.* The comparison in the paragraph above takes its two
sides at different $N$: $0.06$ is the split at the top of the sweep
and $0.29$ the signs at the bottom. Matched, the ratio is $2.98$ to
$3.28$ and flat, not $4.41$ — see Remark [rem:splitvalue]. The
ordering stands; the size of it was overstated.

*Added later.* The fifth $N$ above is censored, not absent: the walk
was truncated at the same $k=10^5$ at which $\beta$ is fitted, and
that cap hides exactly the largest crossings. Remark
[rem:signedgain] separates the two caps, locates the missing
crossing, and finds the gain **rising** there — so "those four numbers
fall" was a statement about what the cap left visible. The gain does
decline, at a rate measured there; it does not decline monotonically.


#### Remark (the knife-edge is a fact about one family) {#rem:residuearithmetic}
<!-- evidence: audit_residue_arithmetic.py -->

Every measurement of $R$ in these papers runs over $N=2\cdot10^5\cdot
2^j$. Those five $N$ are all $2^a5^b$ and share one odd radical —
which Remark [rem:provablehalf]'s evidence made visible by accident,
its density factor $d_L$ coming out identical to four decimals at
every $N$ because the admissible $k$-set never changes. So the
knife-edge of Remark [rem:residuelevel] is measured at **one
arithmetic type**, and the quantity it is measured against is the one
that varies most with type: $\SS(N)(1-A(N))$ runs from $0.073312$ to
$0.374487$ across the test set, a factor of five.

**It does not survive.** With the control P1 reproducing the family
member to $3\cdot10^{-5}$, the seven $N$ of comparable size give

$$
\begin{array}{r|ccccccc}
 \text{odd part} & 3 & 7 & 5 & 3\cdot5 & 3\!\cdots\!17 & 3\!\cdots\!13 & 17\cdot47059\\\hline
 \log K^*_R/\log N & 0.5422 & 0.5591 & 0.5675 & 0.5294 & \mathbf{0.4808} & \mathbf{0.4747} & 0.5424
\end{array}
$$

and rule P2 asked every one of them to clear $\tfrac12$. **At the two
primorial-like $N$ it fails**: $1531530$ gives $0.4808$ and $1621620$
gives $0.4747$, both *below* the square-root barrier. The conditional
reduction of Remark [rem:provablehalf] — strip the elementary half and
let the residue carry the level — **does not reach $\theta'>1/2$ at
the arithmetic where the budget is thinnest.**

That is exactly the place Proposition [prop:onesided] identified: the
threshold is $\asymp N$ for almost all even $N$ and sinks towards
$N/(\log N\log\log N)$ at the primorial-like ones. P3 and P4 confirm
the mechanism rather than complicating it — the spread across types is
$0.0928$, and the exponent regressed on $\log$ of the threshold has
correlation $0.97565$. It is the budget, and nothing else, that moves
the level.

**And that spread survives the test that killed the $k$-exponent's.**
The five $N$ of Remark [rem:residuelevel] are one radical, so their
scatter about their own trend — r.m.s. $0.0049$ — is what this
statistic does with the arithmetic held fixed; seven draws at that
width have an expected span of $0.0134$ by simulation, against the
measured $0.0928$: **a ratio of $6.92$**, where the $k$-exponent of
Remark [rem:kexponent] managed only $1.40$ and had to be withdrawn.
The two claims are now judged on the same footing and only one of them
is a measurement.

Two things to record precisely. The response is weaker than the model
expects: Remark [rem:modeltransfer] prices a budget factor at
$0.1677$ per factor $4.7009$, a slope of $+0.1084$ per natural log,
against a measured $+0.0516$ — **the level is half as sensitive to the
budget as the model says**, which is the only thing here that works in
the programme's favour. And $c_R$ across types spreads by $0.1510$,
matching the $0.1436$ Remark [rem:residueconstant] found across the
family: the constant has no more law across arithmetic than it has
across size.

**So Remark [rem:residuelevel]'s $0.5654$–$0.5799$ must be read as a
statement about $2^a5^b$.** The margin over $\tfrac12$ that Remark
[rem:betafree] showed could not be tuned away by the split constant is
nevertheless erased by changing the arithmetic of $N$, and erased at
the $N$ that were always going to be the hard ones.


#### Remark (the primorial failure is finite-$N$, and nearly reachable) {#rem:primorialladder}
<!-- evidence: lab_primorial_ladder.py -->

Remark [rem:residuearithmetic] found the conditional reduction failing
at primorial-like $N$ and could not say whether that is about the
arithmetic or about the size: seven $N$ of one magnitude cannot
separate them. The arithmetic is not obviously fatal. At $N$ primorial
to $y$ the one-sided margin collapses — $1-A(N)=\sum_{p>y}1/(p(p-1))
\approx1/(y\log y)$ — while $\SS(N)$ grows like $\log y$, so the
budget is of order $N/\log N$ rather than $N$; balancing
$\sum_{k<K}(\log k)c_R\sqrt{N/k}$ against that still gives
$K\asymp N/\log^4N$, whose exponent tends to one.

**It is finite-$N$.** Sweeping a second family — $N=30030\cdot2^j$,
the radical $3\cdot5\cdot7\cdot11\cdot13$ held fixed over a factor
$64$, four times the lever the main family has, with $N=1621620$
recomputed as a control and reproducing its published $0.4747$ to
$10^{-5}$ — the exponent reads

$$
0.4550,\ 0.4595,\ 0.4682,\ 0.4633,\ 0.4688,\ 0.4746,\ 0.4876 .
$$

Rule R2 holds: the slope against $\log N$ is $+0.006623$ with
correlation $0.92378$. R3 holds too — that slope exceeds the
$2^a5^b$ family's $+0.004692$ — but **the ordering is not a
measurement**: the family's slope has standard error $0.002914$
(Remark [rem:slopes]), so the two differ by well under one of them and
"catching up" is not shown. R4 holds:
nothing in the ladder reaches $\tfrac12$, the top being $0.4876$.

The null says what the rise is and is not. A coin on the identical
deviations rises too — eight global sign vectors give slopes from
$-0.003700$ to $+0.014290$, median $+0.008310$, and $\mu$'s
$+0.006623$ sits inside that band. **The rise itself is a fact about
magnitudes, not about $\mu$**, which is what Remark
[rem:residuecancel] would predict since $R$'s sizes are exactly a
coin's. What R3 compares is not $\mu$ against a coin but this radical
against the other, where the budget differs by five — and that
comparison is the content.

What makes this ladder the clean experiment is that the threshold does
not move along it at all. With the radical fixed, $\SS(N)$ and $A(N)$
are constants and $\SS(N)(1-A(N))=0.087306$ at every rung; the rise is
the level alone, against a budget held still.

**And the barrier is nearly within reach.** On the fitted slope the
exponent reaches $\tfrac12$ at $N=10^{7.10}$, with a bracket of
$[10^{7.07},\,10^{7.36}]$ over the leave-one-out extremes of the slope
— whose own drift is $0.2758$, declared. That is a factor of six to
twenty past the top rung, not the eight thousand orders of magnitude
of Remark [rem:provablehalf] or the nine decades of Remark
[rem:leanbracket]. **Extending this ladder three more doublings would
settle it**, and the bracket is narrow for the reason Remark
[rem:marginbracket] gives: the reach is short. Narrow is not the same
as right — Remark [rem:marginoos] measured the octave that bracket
forecast over exactly that short reach and found it outside — so what
the reach buys here is that the ladder can be extended and the bracket
tested, which Remark [rem:primorialreach] then does.


#### Remark (the first bracket to be tested, and it failed) {#rem:primorialreach}
<!-- evidence: audit_primorial_reach.py -->

Remark [rem:primorialladder] forecast that the primorial ladder's
exponent reaches $\tfrac12$ at $N=10^{7.10}$, bracket
$[10^{7.07},\,10^{7.36}]$, and said three more doublings would settle
it. **Every other forecast in these papers is out of reach; this one
was not, and it has now been run.**

The rungs $j=7,8,9$ put the ladder at $3843840$, $7687680$ and
$15375360$ — the last at $10^{7.19}$, past the point estimate and
inside the bracket. The control E1 reproduces all seven published
rungs, and E2 holds with the slope essentially unchanged: $+0.006623$
on seven rungs, $+0.006643$ on ten, correlation improving from
$0.92378$ to $0.95981$ and the leave-one-out spread collapsing from
$0.001827$ to $0.000195$. **The slope was right.**

**E3 is refuted.** The three new exponents are
$0.4824,\,0.4965,\,0.4941$ — none reaches $\tfrac12$, so there is no
crossing inside the ladder, and the point estimate $10^{7.10}$ is
**excluded**: the top rung sits at $10^{7.19}$ with the exponent still
short.

**E4 is not refuted, and the distinction is the whole value of having
had a bracket.** With no crossing to place, E4 cannot be evaluated as
written; what is decidable is exclusion, and the bracket's upper end
$10^{7.36}$ lies *above* the top rung, so a crossing anywhere in
$[10^{7.19},\,10^{7.36}]$ remains open. **The point estimate failed
its first live test and the interval survived it** — which is what an
interval is for, and the first time in these papers that the
difference has been observable. That interval has since been filled
and closed by confirmation: Remark [rem:primorialgap] measures four
$N$ inside it and finds the crossing there.

The point estimate's failure is worth more than its success would have
been. That bracket was built
from the leave-one-out spread of the *slope*, and the slope is now
known to ten times better precision than when it was built. What it
never covered is the scatter of the rungs about the line:

$$
+0.0007,\ +0.0006,\ +0.0047,\ -0.0048,\ -0.0039,\ -0.0027,\ +0.0057,\
-0.0041,\ +0.0053,\ -0.0016
$$

with r.m.s. $0.0039$ against a trend of $0.0046$ per doubling —
**the scatter is eight tenths of what the trend gains in a whole
rung.** A level is therefore crossed several rungs before or after the
line says it will be, and no amount of precision in the slope can know
that. Remark [rem:forecastbracket] required brackets; Remark
[rem:residueconstant] required the constant's drift to be checked
before extrapolating; neither asks the question this failure asks,
which is how far the data sit off their own fit.

Redone with the scatter carried, the line reaches $\tfrac12$ at
$N=10^{7.47}$ with a bracket $[10^{7.21},\,10^{7.72}]$ at one r.m.s.
residual either way. Remark [rem:primorialladder] now publishes the
same correction on its original seven rungs alone, and **that
interval, built from the data that existed before the test, is not
contradicted by the three rungs that refuted the point estimate.** The correction
is not merely a wider bracket: the widened bracket the original data
would have given is consistent with what the new data show, and agrees
with the ten-rung refit to three decimals.


#### Remark (a better instrument that turned out worse) {#rem:primorialshare}
<!-- evidence: lab_primorial_share.py -->

Remark [rem:primorialreach] left a bracket half a decade wide because
the ladder's rungs scatter about their line by $0.8$ of what the trend
gains in a rung. That looked like a fault of the instrument. $K^*_R$
is where a step function first exceeds a level — an integer, discrete,
sensitive to a single term — whereas the question it encodes is not
about a location at all:

$$
\frac{\log K^*_R}{\log N}>\tfrac12
\iff K^*_R>\sqrt N
\iff \rho(N):=\frac{\sum_{k<\sqrt N}(\log k)|R(N;k)|}
{\SS(N)(1-A(N))N}<1 .
$$

$\rho$ is a ratio of two smooth sums at a fixed abscissa. Remark
[rem:muvscoin] made exactly this move for a different question and
resolved in one step what a fitted exponent could not resolve at all.

**It is the worse instrument here, and F2 and F4 are refuted.** The
control F1 holds at all ten rungs — $\rho$ and the exponent agree on
which side of the barrier every rung lies — and F3 holds, $\rho$
falling with slope $-0.035113$. But $\rho$'s residual scatter is
$0.0689$ against a trend of $0.0243$ per doubling, a ratio of $2.829$
where the exponent's is $0.85$, and the forecast bracket **widens**
from $1.6824$ to $1.7034$ decades. (That factor does not survive. Remark
[rem:densenoise] recomputes both ratios where the noise is separable
from the shape and finds the two instruments tied; what was compared
here was a scatter that carried both.) (That comparison has moved since
this was first written. The exponent-based bracket it is measured
against covered only the line's scatter then and was a third as wide;
Remark [rem:laddershape] has since forced it to cover the choice of
shape as well, and the old figure is not reprinted because nothing
computes it any more. The ratio is still the wider instrument, by far
less than it was.)

The mechanism is worth keeping. The two statements are equivalent but
their *trends* are not. Near the crossing
$\log\rho\approx-c\log N\,(e-\tfrac12)$ with $e$ the exponent, so the
exponent carries a division by $\log N$ that $\rho$ does not. That
division damps the fluctuation of $K^*$ without damping the trend by
as much: the exponent's scatter times the mean $\log N$ of $13.4291$
is $0.0524$ against $\rho$'s own $0.0689$ — **the same noise seen
through two lenses** — while the exponent's trend times $\log N$ is
$0.0618$ against $\rho$'s $0.0243$, because $\log\rho$ carries a
second drift $-c(e-\tfrac12)\,d\log N$ that the exponent has divided
away. The null confirms the diagnosis is about the statistic: the coin
arm's own $\rho$ scatters by $0.749$ to $2.529$ of its trend, the same
range.

**And it loses twice over.** Gate check G36 has since forced both
instruments to declare how many functional forms their brackets rest
on. For the exponent on eleven rungs, three forms survive and put
$\tfrac12$ within $0.2146$ decades of each other — the scatter bracket
already covers them and widens by nothing. For $\rho$, **all four
forms survive** and put $\rho=1$ at
$8.3524,\,9.1856,\,9.7957,\,11.4859$: a spread of $3.1336$ decades,
widening its bracket by $2.2818$ to $[10^{7.5007},\,10^{11.4859}]$.
The ratio's scatter is $3.3$ times the exponent's relative to trend,
and its shape ambiguity sits on top of that.

**So the lesson of Remark [rem:muvscoin] does not generalise.** A
ratio at a fixed abscissa beats a fitted exponent when the exponent is
a slope through points; it loses to one when the exponent is itself a
ratio whose denominator grows. The bracket of Remark
[rem:primorialreach] stands as the sharper of the two, by a factor of
seven once both are made to carry the same uncertainties.


#### Remark (the barrier crossed, at the hard arithmetic) {#rem:primorialrung10}
<!-- evidence: audit_primorial_rung10.py -->

The corrected forecast of Remark [rem:primorialreach] put the
primorial ladder's crossing of $\tfrac12$ at $N=10^{7.4684}$. The next
rung is $N=30030\cdot2^{10}=30750720$, at $10^{7.4879}$ — just above
it, and on the fitted line worth $0.5004$ against a scatter of
$0.0039$, so a coin flip written down as one.

**It crossed.**

$$
N=30750720:\qquad K^*_R \text{ gives } \frac{\log K^*_R}{\log N}
= \mathbf{0.5023}.
$$

The control G1 reproduces the previous top rung to $10^{-5}$. G2 holds
— the measured $0.5023$ is $0.0020$ from the line fitted on the ten
published rungs, half the r.m.s. residual. G4 holds, and is the part
worth recording: refitting on all eleven puts the crossing at
$10^{7.4249}$ with bracket $[10^{7.1866},\,10^{7.6631}]$, **inside the
$[10^{7.2133},\,10^{7.7189}]$ that was published before this rung
existed.** That bracket is also robust to the choice of functional
form, which Remark [rem:laddershape] shows is not fixed by these data:
the three shapes surviving at one standard error put the crossing at
$7.4249,\,7.5618,\,7.6394$, a spread of $0.2146$ decades, and the
scatter bracket already contains all of them — so widening it for
shape changes it by $0.0000$. The scatter-corrected bracket has now survived the test
that killed its slope-only predecessor's point estimate, and the
eleven-rung fit is tighter than the ten: slope $+0.006650\to+0.006780$
with the leave-one-out spread at $0.000227$, scatter
$0.0039\to0.0037$, correlation $0.97008$.

What this is. **Conditional on Remark [rem:provablehalf]'s uniformity,
at $N=30750720$ the residue alone carries the truncation past
$\sqrt N$ — at a primorial-like radical, where Remark
[rem:residuearithmetic] found the conditional reduction failing and
Proposition [prop:onesided] locates the worst budget.** The failure
was finite-$N$, as Remark [rem:primorialladder] argued it should be,
and the crossing has been observed rather than extrapolated.

*Added later.* "Observed rather than extrapolated" was said of one
rung, and the margin it clears $\tfrac12$ by is $0.0023$ against this
ladder's own r.m.s. scatter of $0.0037$. **The margin is inside the
floor.** By the standard adopted in Remark [rem:kexponent] and gate
check G37 this rung on its own cannot separate crossed from
fluctuated. It is separated one rung further out — see Remark
[rem:primorialrung11], where the margin is $0.0099$ — so the sentence
above stands, but it did not stand on the evidence given for it here.

What it is not. The uniformity it is conditional on is unproved and
Remark [rem:provablehalf] prices its unconditional form at
$10^{5475}$. One $N$ is not all $N$, and this ladder is one radical
(declared). And $\tfrac12$ is not the $\theta'=0.56$ these papers use:
the same line puts $0.5$ at $\log_{10}N=7.4249$ and $0.56$ at
$\log_{10}N=11.2680$ — **a number Remark [rem:laddershape] then
withdraws, because it is the line's answer and the line is not
distinguishable from shapes that answer very differently.** **What is settled is that
the square-root barrier is not where the primorial arithmetic stops
the residue** — the thing seven $N$ of one size could not tell apart
from the arithmetic itself.


#### Remark (the shape decides, and the data do not fix the shape) {#rem:laddershape}
<!-- evidence: audit_ladder_shape.py -->

Remark [rem:primorialrung10] quoted $\log_{10}N=11.2680$ for
$\theta'=0.56$. That is a linear fit in $\log N$ extrapolated four
decades, and nothing has justified the shape. The heuristic that
governs the level says otherwise: balancing
$\sum_{k<K}(\log k)c_R\sqrt{N/k}$ against a budget of order $N/\log N$
gives $K\asymp N/\log^4N$, an exponent $1-c\log\log N/\log N$, rising
to one with a falling derivative. A saturating $a+b/\log N$ rises to
$a$ and stops. Over a factor $1024$ in $N$ these look alike.

Fitted to the eleven rungs, with the control H1 reproducing the
published line to $2\cdot10^{-6}$:

$$
\begin{array}{l|cc}
 \text{shape} & \text{r.m.s.} & \log_{10}N \text{ at } 0.56\\\hline
 a+b\log N & 0.00372 & 11.2700\\
 a+b\log\log N & 0.00412 & 14.6167\\
 a+b\log\log N/\log N & 0.00432 & 19.6207\\
 a+b/\log N & 0.00473 & 82.5771\\
 1-c\log\log N/\log N & 0.03764 & 7.6204
\end{array}
$$

Rules H2 and H3 asked for an alternative fitting *at least as well* as
the line and **are refuted**: the line is the best of the five. But
that criterion is too strict, and the numbers say by how much. With
eleven points and two parameters the r.m.s. is estimated from nine
degrees of freedom, so its own standard error is
$0.00372/\sqrt{18}=0.00088$, **$23.6\%$**. At one standard error,
three shapes survive — the line at $0.00$, $a+b\log\log N$ at $0.46$,
$a+b\log\log N/\log N$ at $0.69$ — and $a+b/\log N$ is excluded by
only $1.15$.

**They agree about the measurement and disagree about the
extrapolation.** The three surviving shapes put $\tfrac12$ at
$7.4255,\,7.5627,\,7.6404$, a spread of $0.2148$ decades inside the
$0.4765$-decade bracket already published — which is H4, and it holds.
Remark [rem:primorialrung10] now carries the same comparison and finds
its bracket **unchanged to four decimals** when the shape ambiguity is
added to the scatter: at the crossing that was observed, the choice of
form costs nothing.

**Earlier in the ladder it cost a great deal, and the cost collapses
monotonically as the data approach the crossing.** The same five
shapes, fitted to the seven rungs the original forecast had, to the
ten of Remark [rem:primorialreach], and to the eleven here:

$$
\begin{array}{r|ccc}
 \text{rungs} & \text{surviving} & \tfrac12 \text{ spread} &
 \text{bracket widened by}\\\hline
 7 & 4 & 1.4329 & 1.1833\\
 10 & 4 & 0.5960 & 0.3432\\
 11 & 3 & 0.2146 & 0.0000
\end{array}
$$

**The shape ambiguity was real when the forecast was made and
evaporated once the ladder reached the crossing.** That is the honest
shape of extrapolation: the further the reach, the more the form
matters, and the form is the last thing data fix. The line's answer at
seven rungs, $7.4684$, is the one the eleven went on to confirm; the
saturating shape's $8.9013$ is excluded by observation, five times its
own residual scale away.
They put $0.56$ at $11.2700,\,14.6167,\,19.6207$: **a spread of
$8.3508$ decades**, and the shape excluded at $1.15$ standard errors
puts it at $82.5771$ and caps the exponent at $0.5663$ for ever.

So $10^{11.2680}$ is withdrawn. What eleven rungs fix is where the
square-root barrier is crossed, to a fifth of a decade. What they do
not fix is anything about $\theta'=0.56$ — not when, and not whether.


#### Remark (no derived shape is available either) {#rem:laddermodel}
<!-- evidence: audit_ladder_model.py -->

Every shape Remark [rem:laddershape] compared was *fitted*. The
heuristic derives one: $K^*_R$ solves
$\sum_{k<K}(\log k)|R|=\SS(1-A)N$, and with $|R|\approx
c_R(N)\sqrt{N/k}$ and $c_R\approx\gamma\sqrt{\log N}$ the crossing
follows with no shape freedom at all. Remark [rem:residueconstant]
confirmed that model to $0.7\%$ at the $2^a5^b$ family. **It fails
here, and at its foundation.**

The control J1 reproduces all eleven exponents to $5\cdot10^{-5}$.
J2 is refuted, though only just: the model's crossings sit within
$0.9479$ to $1.0367$ of the measured, one rung past the $5\%$ asked.
J4 is refuted outright — the model puts $0.56$ at $10^{9.0269}$,
*below* every fitted shape's answer, so it is a fifth opinion and not
a resolution.

**J3 says why, and is the finding.** The constant it rests on is not a
constant:

$$
\frac{c_R}{\sqrt{\log N}}:\quad
0.3724,\ 0.4235,\ 0.4239,\ 0.4757,\ 0.5185,\ 0.5140,\ 0.4985,\
0.5645,\ 0.5456,\ 0.5968,\ 0.6139,
$$

rising monotonically, a spread of $0.4790$ of its mean against the
$0.1436$ that family showed. Fitted directly, $c_R\asymp(\log
N)^{1.3838}$ with correlation $0.98535$ and a leave-one-out spread of
$0.0740$ — **not $\tfrac12$, and not a wobble around $\tfrac12$ but a
different law**, growing $2.8$ times faster in the exponent. The
model's residuals are correspondingly systematic, positive at the
bottom of the ladder and negative at the top, with r.m.s. $0.0112$
against the line's $0.0037$.

So the primorial radical does not merely shift the constants of
Remark [rem:heuristic]. Remark [rem:cRwindow] tried to break that
reading by attributing the growth to the widening window $c_R$ is
averaged over, and failed to: the growth is real, and its cause is
sharper still.


#### Remark (the exponent moves, not the constant) {#rem:cRwindow}
<!-- evidence: audit_cR_window.py -->

Remark [rem:laddermodel]'s headline — $c_R\asymp(\log N)^{1.3838}$ at
the primorial radical, "a different law" — rests on one convention:
$c_R$ is averaged over $k<K^*_R$, and $K^*_R$ grows along the ladder
from $109$ to $5773$. Remark [rem:modeltransfer] found precisely this
failure for the other constant, $c$ moving by five per cent between
two $k$-ranges at the *same* $N$. If $|R|=c_R(N)\sqrt{N/k}$ exactly,
the window cannot matter; if it does, the growth is an artefact.

**The attempt to break it failed, and informatively.** The control K1
reproduces $1.3838$ exactly. K2 predicted a *smaller* exponent from a
fixed window and **is refuted**: $k<300$ at every rung gives
$1.5918$, steeper, and $k<N^{1/4}$ gives $3.5892$. K3 falls with it.
The growth is not the window.

**K4 finds the real cause and is itself refuted at one rung.** Fitting
the octave means of $|R|$ against $N/k$ at each rung gives

$$
0.4223,\ 0.3991,\ 0.4192,\ 0.4648,\ 0.4652,\ 0.4971,\ 0.4791,\
0.5029,\ 0.4627,\ 0.5168,\ 0.5074,
$$

**rising** with slope $+0.014647$ against $\log N$, correlation
$0.85888$, leave-one-out spread $0.000993$ — and $0.3991$ at
$N=60060$ is outside the $[0.40,0.60]$ the rule asked for. So $|R|$ is
**not** $(N/k)^{1/2}$ at this radical below $N\approx9.17\cdot10^6$;
the exponent climbs to $\tfrac12$ from below and only reaches it near
the top of the ladder.

That settles what $c_R$ was. Dividing $|R|$ by $\sqrt{N/k}$ when the
true power is $a(N)<\tfrac12$ leaves $(N/k)^{a-1/2}$, which falls with
$N/k$ — so the "constant" has no value independent of where the
average stops, and its apparent growth is the mismatch closing as
$a(N)\to\tfrac12$. **The right statement is not that the residue's
constant obeys a different law at the hard arithmetic, but that its
$k$-exponent is still short of a half there and rising.** Remark
[rem:laddermodel]'s conclusion stands — no derived shape is available
— for this reason rather than the one it gave.


#### Remark (the k-exponent is too noisy to carry an arithmetic) {#rem:kexponent}
<!-- evidence: audit_residue_kexponent.py -->

Remark [rem:cRwindow] found $|R|$'s $k$-exponent rising towards
$\tfrac12$ along the primorial ladder, and Remark
[rem:residuearithmetic] found the level falling with the radical and
attributed that to the budget. If the $k$-exponent also moved with the
radical there would be a second source, and the seven $N$ of the
arithmetic test set — all near $1.6\cdot10^6$ — isolate the radical
with $N$ nearly fixed.

The control L1 holds, the family member reproducing its published
$k$-exponent to $0.0061$. L2 holds: the seven span $0.0727$. **L3 and
L4 are refuted, and in the direction opposite to the one predicted** —
the exponent *rises* with the number of odd prime factors (slope
$+0.006207$, correlation $0.52062$) and *anti*-correlates with the
level ($-0.46098$).

**None of that is a measurement.** The noise floor for this statistic
is available: the eleven rungs of the primorial ladder sit at one
radical, so their scatter about their own trend — r.m.s. $0.0191$ — is
what the exponent does with the arithmetic held fixed. Seven draws at
that width have an expected range of $0.0519$ by simulation, against
the measured $0.0727$: **a ratio of $1.40$.** And the same radical
appears twice across the two files at nearly the same $N$ —
$1621620$ here at $0.5341$, $1921920$ on the ladder at $0.4791$ —
**apart by $0.0550$, three quarters of the whole seven-radical span,
at one radical.**

So L2's "hold" must not be read as arithmetic dependence, and L3's and
L4's signs must not be read at all. The level's dependence on the
arithmetic remains what Remark [rem:residuearithmetic] measured — the
budget, at correlation $0.97565$ — with no second source shown. What
survives from Remark [rem:cRwindow] is the ladder's own rise, which is
about six times its scatter and stands; what does not survive is any
comparison of the $k$-exponent between radicals at a single $N$.


#### Remark (repairing the model's constant does not derive the shape) {#rem:ladderderived}
<!-- evidence: audit_ladder_derived.py -->

Remark [rem:laddermodel] refuted the heuristic's derived crossing and
identified where: the constant it rests on is not $\gamma\sqrt{\log
N}$ but $c_R\asymp(\log N)^{1.3838}$. It was refuted with its own
replacement in hand. Since OPEN item 1's bottleneck is that every
surviving shape is *fitted* (Remarks [rem:laddershape12],
[rem:shapetrust]), the repair is worth making: with $|R|=\gamma(\log
N)^{\rho}\sqrt{N/k}$ and $\sum_{k<K}(\log k)k^{-1/2}\sim
2\sqrt{K}(\log K-2)$, writing $K=N^{e}$ gives

$$
e \;=\; 1 \;-\; 2(1+\rho)\,\frac{\log\log N}{\log N}
\;+\; \frac{d}{\log N},
$$

whose leading coefficient is *derived* and which has one free constant
against every rival's two. Both are computed here from a ladder
recomputed off its own sieve; the control M1 reproduces all eleven
published exponents and constants and Remark [rem:primorialrung11]'s
twelfth exponent to $5\cdot10^{-5}$, and adds that rung's constant,
$c_R=2.5790$. M2 holds: eleven rungs give $1.3838$ at correlation
$0.98535$, twelve give $1.3518$ at $0.98651$.

**M3 is refuted, and not marginally.** With $\gamma=0.053849$ and
$\rho=1.351815$ fitted to $c_R$ alone and nothing fitted to the
exponents, the predicted crossings run $0.9068$ to $1.0643$ of the
measured — *wider* than the $0.9479$ to $1.0367$ the uncorrected model
gave. Correcting the constant's law does not improve the crossing at
all, so what fails is the $\sqrt{N/k}$ profile underneath it, not the
constant sitting in front.

**M4 is refuted at $3.84$ times the fitted r.m.s.** The derived shape
fits at $0.01419$ against $0.00370$, with residuals $+0.0200$,
$+0.0167$, $+0.0169$ at the bottom of the ladder and $-0.0194$,
$-0.0196$, $-0.0201$ at the top — a systematic curvature, not scatter.
Freeing the coefficient recovers the fit exactly ($5.7691$ against the
derived $4.7036$, r.m.s. $0.00395$), so the functional form is right
and the derived *coefficient* is $23$ per cent low. But the ladder
cannot say so: the two regressors $\log\log N/\log N$ and $1/\log N$
correlate at $0.99883$ over $2.78$ decades, so the coefficient and the
constant are not separable here, and the extra parameter buys
$0.01024$ of r.m.s. without being measurable on its own.

The consequence for OPEN item 1 is worse than a third opinion. The
derived shape reaches $0.56$ at $10^{8.9794}$ — within $0.05$ decades
of Remark [rem:laddermodel]'s $10^{9.0269}$, so the answer is stable
under the repair — against the fitted survivors' $10^{11.0762}$ and
$10^{13.8607}$. **But it reaches $1/2$ at $10^{6.7141}$, while the
ladder measures $0.5023$ at $10^{7.4879}$ and $0.5099$ at
$10^{7.7889}$.** The derivation is wrong *inside* the measured range,
by seven tenths of a decade, so its forecast carries no authority
against the fits. The fitted shapes are all there is, and Remark
[rem:shapetrust]'s reading stands: statements about $1/2$ are the
data's, statements about $\theta'$ are the shape's.


#### Remark (how far the ladder can be read) {#rem:shapetrust}
<!-- evidence: audit_shape_trust.py -->

Remark [rem:laddershape12] leaves $\theta'=0.56$ between two shapes
$2.78$ decades apart, and Remark [rem:primorialrung11] reads the same
ladder confidently at $\tfrac12$. Both are true, and the number that
reconciles them had not been computed: the $N$ up to which the
surviving shapes are interchangeable — differing by less than the
ladder's own r.m.s. scatter.

**C1 holds** (the fits reproduce exactly) and **C2 is refuted, badly.**
With the twelve-rung scatter at $0.00370$ and the top rung at
$\log_{10}N=7.7889$, the line and $a+b\log\log N$ first differ by more
than that scatter at

$$
\log_{10}N = 8.1253 ,
$$

**$0.34$ decades above the top rung**, not the two decades registered.
The gap grows to $0.0080$ at $10^9$, $0.0140$ at $10^{10}$ and
$0.0288$ at $10^{12}$.

**C3 holds**, which is the consequence: both $\theta'$ forecasts,
$11.0762$ and $13.8607$, lie far outside that range, so neither is a
reading of the data. The crossings of $\tfrac12$, at $7.3605$ and
$7.4266$, lie inside it — which is why Remark [rem:primorialrung11]'s
statement about $\tfrac12$ stands while the $\theta'$ statement does
not.

**C4 holds**: the flatness adjudication of Remark
[rem:flatnessshape] has the same structure — its two shapes part at
$\log_{10}N=8.6994$ while the bound one of them would reach sits at
$28.6782$.

*Added later.* That arrival is withdrawn — see [rem:fillfield] and
[rem:fieldreach], which publish no arrival for $F$ at all. The
structural point C4 makes survives it.

So the operative statement about this ladder is narrow and exact:
**it can be read to about $N=10^{8.1}$, a factor $2.2$ past its top
rung, and no further.** Everything the papers say about $\tfrac12$ is
inside that; everything about $\theta'$ is outside it.


#### Remark (what one more rung buys the shape question) {#rem:laddershape12}
<!-- evidence: audit_ladder_shape12.py -->

Remark [rem:laddershape] adjudicated five shapes on **eleven** rungs
and left $\theta'=0.56$ undetermined: three shapes survive at one
standard error of the best r.m.s. and put the crossing $8.35$ decades
apart. Remark [rem:primorialrung11] then measured a twelfth rung and
that adjudication was never redone.

**J1 holds** — all five r.m.s. and all five crossings reproduce
exactly — and **J2 holds**: on twelve the line is still the best.
**J3 holds too.** The best r.m.s. moves from $0.00372$ to $0.00370$
and its own standard error from $23.6\%$ to $22.4\%$, and the field
narrows:

$$
\begin{array}{l|ccc}
\text{shape} & \text{r.m.s. on }12 & \text{s.e. from best} &
  \log_{10}N\text{ at }0.56\\\hline
a+b\log N & 0.00370 & 0.00 & 11.0762\\
a+b\log\log N & 0.00433 & 0.76 & 13.8607\\
a+b\log\log N/\log N & 0.00464 & 1.12 & 17.5614\\
a+b/\log N & 0.00519 & 1.80 & 41.0798\\
1-c\log\log N/\log N & 0.03934 & 43.02 & 7.7105
\end{array}
$$

Three survivors become **two**, and the spread of their forecasts
falls from $8.3508$ decades to $2.7845$. **J4 holds**: it is still
more than a decade, so $\theta'=0.56$ is still not located.

What this prices is the rung itself. One doubling of $N$ — the
computation of Remark [rem:primorialrung11] — bought a factor $3$ in
the shape ambiguity, and the two survivors now agree on $\tfrac12$ to
within $0.07$ of a decade ($7.3605$ and $7.4266$), which is behind the
ladder and so not an extrapolation at all. The remaining ambiguity is
entirely about where the ladder meets $0.56$, and closing it by rungs
costs a doubling of $N$ for each factor of about three.


#### Remark (the premise, checked along the axis it is about) {#rem:provableuniformity}
<!-- evidence: audit_provable_uniformity.py -->

Remark [rem:provablehalf]'s premise is uniformity **in $k$**. What it
tests is the shape, and it reports the constant maximised by inner
length: $1.2119,\,1.0710,\,0.7309,\,0.3363$ at $N/k\ge2,8,32,128$.
That is the axis the classical estimate already controls. Holding the
inner length in an octave and letting $N$ run, the modulus $k=N/x$
moves by a factor $16$ — and nobody had looked there.

**F1 holds** (all four restricted maxima reproduce to $0.00005$) and
**F3 holds**: a resolved rise at two octaves of six. **F2 and F4 are
refuted, and each says something.**

$$
\begin{array}{l|ccc}
 N/k \text{ octave} & \text{growth} & \text{slope in } \log k & t\\\hline
 [8,32) & 1.3786 & -0.555097 & 2.13\\
 [32,128) & 2.0964 & -0.302804 & 1.91\\
 [128,512) & 1.2865 & +0.039912 & 0.91\\
 [512,2048) & 1.0878 & +0.023475 & 1.84\\
 [2048,8192) & 1.5307 & +0.167112 & 3.37\\
 [8192,32768) & 1.6719 & +0.227768 & 2.18
\end{array}
$$

**The axis is not flat, and it changes sign.** At short inner sums the
ratio *falls* with the modulus; at long ones it *rises*, resolvably at
$3.37$ and $2.18$ standard errors. So the uniformity is not merely
unproved: along the axis it is about, the accessible range shows a
resolved upward drift.

**But not where the mass is.** The two rising octaves carry $0.0673$
of $\sum(\log k)|P|$ — Remark [rem:predictable] put $0.77$ to $0.90$
of it at $N/k\le10^3$, which is the falling side. The drift is real
and sits where the elementary sum is thin.

F4 failed on a badly chosen threshold rather than on the object: the
argmax's inner length is $7,\,13,\,30,\,63,\,108$ — rising with $N$
because the argmax sits at the $k$-cap ($28549$ to $29451$) at every
$N$, exactly as Remark [rem:provablehalf] says. Asking it to stay
below a fixed $32$ was asking the wrong question.

*Added later.* "The drift sits where the elementary sum is thin" is a
statement about a factor-$16$ lever. With $512$ it is the other way
round; see Remark [rem:uniformityreach]. Nothing else here changes.


#### Remark (the drift is real, and nowhere near the bound) {#rem:uniformityreach}
<!-- evidence: audit_uniformity_reach.py -->

[rem:provableuniformity] checked [rem:provablehalf]'s premise on the
axis the classical estimate does not control — fixed inner length,
varying modulus — over a modulus range of a factor $16$, and closed
by saying the drift sits where the mass is not. The statistic turns
out to cost seconds, so the same check runs over a factor $512$: ten
doublings to $N=1.024\cdot10^8$, the $k$-cap of $30000$ and everything
else imported unchanged. The five published maxima reproduce to
$0.000033$ (A1).

**The drift is not a short-lever effect** (A2). Four octaves are
resolved rising, not two, and the two that were rising survive with
better resolution: $[2048,8192)$ at $+0.005283$ ($t=2.53$),
$[8192,32768)$ at $+0.006071$ ($t=5.60$), and the two longer
inner-length octaves the extra reach opens, $[32768,131072)$ at
$+0.004537$ ($t=3.88$) and $[131072,524288)$ at $+0.001827$
($t=3.03$). The short inner lengths fall, as they did:
$-0.212164$ and $-0.137950$ at $t=2.59$ and $3.28$. **The axis is not
flat, and the sign of the drift is set by the length of the inner
sum.**

**And it no longer sits where the mass is not** (A3). The share of
$\sum(\log k)\lvert P\rvert$ carried by resolved-rising octaves runs
$0.0266$ to $0.9804$ across the doublings, against the published
$0.0673$. The reason is structural and has to be said with the
result: the $k$-cap is fixed, so as $N$ grows every admissible $k$ has
a long inner sum, and long inner sums are exactly the rising ones. At
the top $N$ an uncapped sweep would have added only $2$ more $k$, so
within *this* field — which is the field [rem:provablehalf]'s
elementary half is defined on — the drift is where the mass is.

**But it is nowhere near the bound** (A4). The overall maximum ratio
runs $1.2119$ down to $0.0843$ with slope $-0.446123$ at $t=37.31$.
The rising octaves rise at five thousandths per unit $\log N$ while
the maximum falls at four tenths. **So the premise is unproved and
measurably non-uniform, and at every reach that can be computed the
constant it needs is improving by two orders faster than the drift
threatens it.** That is not a proof of the uniformity and does not
become one by being repeated at larger $N$; it is the exact shape of
what a computation can say about a statement quantified over all $k$
and all $N$.


#### Remark (two answers to two questions, 5465 decades apart) {#rem:provableforecast}
<!-- evidence: audit_provable_forecast.py -->

Remark [rem:provableshare] showed the implied constant $A$ is not one:
it is each $N$'s maximum ratio and falls from $1.2119$ to $0.3487$.
Remark [rem:provablehalf]'s forecast solves $A\,d_L\,I(u)=\SS(1-A)$
with $A$ frozen at that maximum, so the question is what the forecast
becomes when $A$ is allowed to move.

**E1 holds** — $d_L$, the frozen $A$ and the forecast $5474.8$
reproduce exactly. **E2 holds**: $A$'s decay is $-0.451786$ per unit
$\log N$ at $10.30$ standard errors. **E3 holds**: its shape is not
pinned — fitting $\log A$ against $\log N$ and against $\log\log N$
gives r.m.s. $0.074485$ and $0.087846$, a gap of $0.44$ of the
r.m.s.'s own standard error.

**E4 holds, and the size of it is the point.** Extrapolating $A$ by
each shape and solving the same equation:

$$
\begin{array}{l|c}
A \text{ frozen at the sweep maximum} & 10^{5474.8}\\
A\sim N^{-0.451786} & 10^{8.96}\\
A\sim(\log N)^{-6.085957} & 10^{10.29}
\end{array}
$$

The two unfrozen answers differ from each other by $1.34$ decades and
from the frozen one by **$5465$**.

**They answer different questions and neither corrects the other.**
$A$ is a measured maximum ratio, not a constant any theorem supplies.
Freezing it at the largest observed value is the conservative choice
for an *upper bound* on the spend, which is what Remark
[rem:provablehalf] wants; letting it fall describes what the data do,
and extrapolated far enough would put the "bound" below $|P|$ itself
and stop being a bound. So $10^{5475}$ is where the bound *provably*
pays and $10^{9}$–$10^{10}$ is where the measured ratio *would* — and
the gap between them is the whole distance between a theorem's
constant and a measurement.

What this settles for OPEN item 2: the published figure is exact
about a frozen constant and says nothing about the bound's behaviour
once the constant moves, because $A$'s shape over four thousand
decades is an extrapolation this repository refuses on its own
standard (Remark [rem:forecastbracket]).


#### Remark (the overspend worsens only with the constant frozen) {#rem:provableshare}
<!-- evidence: audit_provable_share.py -->

Remark [rem:provablehalf]'s rule W3 reports the classical bound
spending $13.98$ to $19.83$ times the budget at accessible $N$ — "not
merely useless but *getting worse* across the sweep". Two things had
never been asked of that: whether the rise clears its own standard
error, and what it does when the implied constant $A$ is not held
fixed.

**D1 and D2 hold.** The five shares reproduce to $0.0043$, and with
$A=1$ the rise is $+0.125950$ per unit $\log N$ at **$50.35$ standard
errors** — as firm as anything in these papers.

**D3 is refuted, and it reverses the sentence.** $A$ is a per-$N$
maximum and it collapses from $1.2119$ to $0.3487$ across the sweep.
With each $N$'s own constant the slope is

$$
-0.325836 \quad\text{at } 7.03 \text{ standard errors} ,
$$

so the overspend **improves** with $N$ rather than worsening. Neither
convention is wrong — Remark [rem:provablehalf] freezes $A$ at the
sweep maximum because it wants an upper bound and says so, noting that
$A$ "is a function of where one looks". What may not be said is that
the overspend worsens, full stop: the direction is a fact about the
convention. The forecast at $10^{4785}$ is untouched, since it uses
the frozen constant deliberately.

**D4 is refuted the other way, and usefully.** Five points *can* tell
the mechanism: fitting the share against $\log N$ gives r.m.s.
$0.004248$ and against $\log\log N$ gives $0.000653$, a gap of
$13.47$ standard errors of the r.m.s. itself. The share grows like
$(\log N)^{1.707104}$, which is the shape W3's stated mechanism
predicts — $\exp(-c\sqrt{\log x})$ falling more slowly than a power of
$\log$ rises. **The mechanism is confirmed; only the direction of the
trend was convention-dependent.**


#### Remark (the exposure that remark names, tested) {#rem:provablearithmetic}
<!-- evidence: audit_provable_arithmetic.py -->

Remark [rem:provablehalf] says of its forecast's inputs: "the first
does not drift at all here — every $N$ in the sweep has the same odd
radical, so the admissible $k$-set and $d_L$ are identical across it,
**a fact about this sweep and not a general one**." That is an
exposure named and left untested, and since Remark
[rem:arithmeticreach] the arithmetic dependence is known to be
standing rather than closing.

Recomputed independently, **Z1 holds**: $d_L=0.3994$ to $0.00004$ and,
with the published $A$, the forecast to $0.02$ of a decade. Across the
seven arithmetic types of Remark [rem:residuearithmetic], **Z2
holds**:

$$
d_L = 0.3318,\,0.4288,\,0.3994,\,0.2629,\,0.1730,\,0.1848,\,0.4726,
$$

a relative spread of $0.9309$ against a floor of exactly $0.0000$ —
with the radical held, the admissible $k$ are literally the same set
and $d_L$ cannot move at all.

**Z3 holds and Z4 holds.** With every input measured per type the
forecast runs $4784.6$ to $5984.4$; holding $A$ at the published
sweep-maximum, which is the remark's own convention, it runs $5474.8$
to $6184.3$. Over both conventions: $[10^{4784.6},\,10^{6184.3}]$.
That is inside the $c$-sweep's $[10^{2092.7},\,10^{13093.3}]$ — the
analytic exponent still dominates — but **outside the $A\,d_L$
sensitivity the remark declares**, $[10^{4838.5},\,10^{6139.9}]$, at
both ends.

The mechanism is that the two channels pull against each other. A
primorial-like $N$ admits fewer $k$, which lowers $A\,d_L$ and brings
the crossing in; and it carries a thinner budget — $0.0733$ against
the family's $0.3745$ — which pushes it out. The budget wins: the two
primorial-like types are the two latest, at $5984.4$ and $5883.9$.

Nothing here changes a conclusion, and nothing could: the whole span
is nine thousand orders of magnitude from anything computable. What
changes is the bracket. **The declared sensitivity was not a bound on
the exposure the remark itself named**, and it is now measured.


#### Remark (the budget gap is not a constant, and both halves pay it alike) {#rem:budgetgap}
<!-- evidence: audit_budget_gap.py -->

"A budget factor of $4.7009$ costs $0.1677$ in the exponent" is quoted
across these papers as a constant of the method, and three scripts
read it out of Remark [rem:modeltransfer]'s results file. **It is a
mean.** The series is $0.1824,\,0.1688,\,0.1706,\,0.1635,\,0.1532$ —
a fall of more than a tenth of itself across the sweep.

Recomputed independently, **Y1 holds** at every $N$, and **Y2 holds**:
the gap's least-squares slope against $\log N$ is $-0.009165$ at
$5.04$ standard errors. It is a declining quantity, and $0.1677$ is
attained at no $N$ in the sweep it was averaged over.

**Y3 and Y4 hold, and they are the useful part.** The gap had only
ever been measured on $H$, while every use of it is about $R$.
Computed for $R$ it is $0.1823,\,0.1775,\,0.1794,\,0.1696,\,0.1583$ —
within $0.0088$ of $H$'s at every $N$, and declining at $-0.008079$,
a difference of $-0.001087$ against two standard errors of
$0.005525$. **What the operative budget costs is a property of the
budget ratio and not of the half being truncated**, which is what
using it for $R$ has been assuming without evidence.

One consequence for Remark [rem:residuearithmetic]. It divided the
quoted mean by $\log$ of the budget factor to predict a response of
$+0.1084$ per natural log of threshold, and measured $+0.0516$ across
seven radicals. At its own $N=1.6\cdot10^6$ the gap gives $0.1057$,
not $0.1084$. The disagreement is not the averaging: it survives at
roughly a factor two either way, and remains what that remark said it
was — the harder radicals admit fewer $k$, so their truncation moves
less than the budget alone would move it.


#### Remark (one variable, and it explains more as N grows) {#rem:arithmeticonevar}
<!-- evidence: audit_arithmetic_onevar.py -->

Remark [rem:arithmeticreach] leaves the arithmetic spread standing at
seven times its floor over a factor $8$ in $N$, and says nothing about
whether the spread is *one* thing. Remark [rem:residuearithmetic]
attributed it to the budget with a correlation of $0.97565$ — at one
scale, and Remark [rem:kexponent] is the standing warning about what
one scale can show.

The regressor does not move here. Doubling $N$ leaves the odd radical
alone, so $\SS(N)(1-A(N))$ is the same number at $1,2,4,8$ times each
of the seven $N$: the abscissae are fixed and only the exponents move.

**A1, A2 and A4 hold.** The base scale reproduces to $0.00004$ with
slope $0.0516$ and correlation $0.97565$ exactly; the correlation
stays above $0.9$; the slope's spread over the four scales is
$0.0063$.

**A3 is refuted, and in the direction that favours one variable.** The
residual spread after removing the budget is

$$
0.0218,\quad 0.0240,\quad 0.0129,\quad 0.0103
$$

against a floor of $0.0133$ — **it crosses the floor between $2N$ and
$4N$**, and the correlation rises $0.97565\to0.98390\to0.99433\to
0.99577$. So a second arithmetic variable is visible at the accessible
sizes and gone at four and eight times them.

Two of the seven types keep the sign of their residue at all four
scales — $\text{odd}=5$ above the line and $17\cdot47059$ below — which
is what a second variable would look like; both shrink, from $+0.0090$
to $+0.0056$ and from $-0.0128$ to $-0.0047$.

So OPEN item 3 stands as a limitation but not as a mystery: the
caveat does not close, and what it is made of is one measured
quantity, the budget, with the residue falling inside the noise by
$4N$.


#### Remark (the one-radical caveat does not close) {#rem:arithmeticreach}
<!-- evidence: audit_arithmetic_reach.py -->

Every level result in these papers carries the caveat that its sweep
has one odd radical, and Remark [rem:residuearithmetic] is why: a
spread of $0.0928$ across seven radicals against a floor of $0.0134$.
**That was measured at one scale.** One scale cannot say whether a
dependence is standing or is being worked off, and Remarks
[rem:slopereach] and [rem:primorialrung11] had made the second
possibility concrete — the ladder rising at $+0.007013$ against the
family's $+0.005112$.

Multiplying $N$ by $2$ leaves its odd radical alone, so the same seven
types follow up the scale exactly. At $N,\,2N,\,4N,\,8N$ — twenty-eight
measurements, with **Q1 reproducing all seven published exponents**
(worst discrepancy $0.00004$) — the spread is

$$
0.0928,\quad 0.1063,\quad 0.1011,\quad 0.0981 .
$$

**Q3 is refuted.** It does not close. Its own slope is $+0.001535$ at
$0.36$ standard errors: flat, over a factor $8$ in $N$. **Q4 holds** —
the spread is still $7.4$ times its floor at the top — and **Q5
holds**: the same two primorial-like radicals are the bottom two at
every scale. So the caveat is a property of the method at these
sizes, not an artefact the range removes.

**Q2 holds but must not be over-read, and the reason is the second
finding.** All seven slopes are positive, but the number of odd prime
factors does not order them: the slowest is $3\cdot5\cdot7\cdot11\cdot
13\cdot17$ at $+0.002239$ — six odd primes, and $0.59$ standard
errors, so not resolved at all — while the fastest is
$17\cdot47059$ at $+0.012871$ with two. Remark
[rem:primorialladder]'s rule R3 read one radical rising faster than
one other as the hard arithmetic catching up; across seven radicals
there is no such ordering. What survives is that the individual
ladders rise, which Remark [rem:primorialrung11] establishes for the
one it measures.


#### Remark (both losses decline together) {#rem:splitvalue}
<!-- evidence: audit_split_value.py -->

Remark [rem:residuesigned] closes with "Remark [rem:splitbudget]
measured the whole elementary/residue division at about $0.06$ in
$\theta'$, and this one step at $0.21$ to $0.29$" — a point against a
range. Both are series, and Remark [rem:signedgain] has since shown
one of them declining, so the other's shape is not optional.

The independent implementation of Remark [rem:slopereach] reproduces
Remark [rem:splitbudget] exactly: **V1 holds** on all fifteen
crossings and **V2 holds** with the worst exponent discrepancy
$0.000046$. What the split is worth — the residue's exponent minus
$H$'s, at the same $N$ — is

$$
0.0925,\ 0.0797,\ 0.0779,\ 0.0711,\ 0.0666 ,
$$

not $0.06$. **V3 holds**: it falls at $-0.008722$ per unit $\log N$,
$6.70$ standard errors.

**And V4 holds, which is the point.** Taken at matching $N$, the two
losses stand in the ratios
$3.1709,\,3.1681,\,3.0108,\,2.9828,\,3.2772$ — slope $+0.003946$ at
$0.06$ standard errors, a flat line. The signs are worth about three
times the split, and that factor is a property of the range and not
of a particular $N$: both losses are being worked off together.

What the published sentence did was divide the top of one range by the
bottom of the other, which takes the two sides at different $N$ and
gives $4.41$ — $1.13$ above the largest matched ratio. The claim
survives, smaller and better founded: **of everything the programme
discards, the signs across $k$ in the residue are worth about three
times the split, at every $N$ measured.**


#### Remark (the largest loss is a decline, and it is not monotone) {#rem:signedgain}
<!-- evidence: audit_signed_gain.py -->

Remark [rem:residuesigned] is the largest number in the chain, and it
had been computed once. Recomputed by the independent implementation
of Remark [rem:slopereach], **W1 and W2 hold exactly**: the absolute
crossings $993,\,1447,\,2019,\,3319,\,5923$ and the signed
$35597,\,37623,\,48957,\,68669$ come back digit for digit. Nothing in
that remark's arithmetic is in doubt.

Its *reading* is. $\beta$ is fitted on $k<10^5$ and the walk was
truncated at the same $k$, and the cap hides exactly the largest
crossings — the fifth $N$ was printed as "none", the factor as
$>16.9$. **Censoring that removes the big values is not a missing
data point, it is a bias in whatever trend the survivors show.** The
two caps need not be the same number: the fit is a property of the
split, the walk is a sum. Continuing the walk to $k<4\cdot10^5$ with
$\beta$ untouched, **W3 holds**:

$$
N=3.2\cdot10^6:\qquad K^*_{\text{signed}}=155333,\quad
\frac{\log K^*}{\log N}=0.7980,\quad \text{gain } +0.2181 .
$$

The factor is $26.2$, not "$>16.9$" — and the gain is **above** the
$+0.2121$ before it. The published sequence
$+0.2932,\,+0.2526,\,+0.2346,\,+0.2121$ looked monotone because the
one point that would have broken the pattern was the one the cap
removed.

**W4 holds all the same.** Over five points the gain's least-squares
slope against $\log N$ is $-0.027526$ at $4.20$ standard errors,
two-sigma interval $[-0.040634,\,-0.014419]$: about $0.0191$ per
octave, against the $0.0530$ that Remark [rem:signedlevel] measured
for $H$ as a whole. So the largest single discard in the chain is
being worked off by the range — **it is not a constant of the
method** — while remaining, at every $N$ measured, four times $H$'s.

No $N$ at which the gain would vanish is quoted. That is an
extrapolation over a factor $16$, and Remark [rem:forecastbracket] is
the standing reason not to make it.


#### Remark (the rung that separates crossed from fluctuated) {#rem:primorialrung11}
<!-- evidence: audit_primorial_rung11.py -->

Remark [rem:primorialrung10]'s crossing clears $\tfrac12$ by $0.0023$
against a floor of $0.0037$, so the next rung
$N=30030\cdot2^{11}=61501440$ decides whether it was a crossing or a
fluctuation. Recomputed by the independent implementation of Remark
[rem:slopereach] — **P1 reproduces the published rung to $0.000034$** —
the answer is:

$$
N=61501440:\qquad \frac{\log K^*_R}{\log N} = \mathbf{0.5099}.
$$

**P2 and P3 hold.** The barrier stays crossed, and now by $0.0099$,
which is $2.7$ times the ladder's floor. Two rungs above $\tfrac12$,
the second outside its own noise: the crossing is an observation.

**P4 is refuted**, and this is where the registered rule was the
wrong instrument rather than the ladder the wrong shape. The line on
eleven rungs predicts $0.5057$ and the measurement is $0.5099$, off by
$0.0042$ against an in-sample r.m.s. of $0.0037$ — but an
out-of-sample point is not compared with an in-sample r.m.s. The
prediction standard error at that abscissa,
$s\sqrt{1+\tfrac1n+(x_0-\bar x)^2/\sum(x-\bar x)^2}$, is $0.0049$,
and $0.0042$ is $0.87$ of it. **The confirming test is the refit:** on
twelve rungs the slope is $+0.007013$ at $14.32$ standard errors and
the scatter is $0.0037$, unchanged. A ladder that bent would have
raised it.

What is *not* touched. The shape question of Remark
[rem:laddershape] is about where the ladder meets $\theta'=0.56$,
which is ahead of it; $\tfrac12$ is now behind it, and no forecast is
made here.

*Added later.* A thirteenth rung has been computed, the first ever
measured after the line was fitted. See Remark [rem:rung12].


#### Remark (the thirteenth rung, out of sample) {#rem:rung12}
<!-- evidence: audit_primorial_rung12.py -->

OPEN item 1 was blocked on one sentence: the ladder is trusted only so
far and "the only way to push it is to raise $N$". The sign-axis
cycles have since shown $N=10^8$ costs seconds here, and the ladder's
rungs are $30030\cdot2^j$, so $30030\cdot2^{12}=123002880$ at
$\log_{10}N=8.0899$ is affordable — three tenths of a decade past the
top published rung and, more to the point, **the first rung computed
after the line was fitted**. Rung 11 reproduces to $0.000006$ with the
same $K^*_R=9367$ (B1).

**The barrier is cleared again, and by more than before** (B2). The
new exponent is $0.5178$, a margin of $0.0178$ over $\tfrac12$ against
the ladder's scatter $0.0037$ — $4.8$ times the floor, where rung 10
had $0.0023$ (inside it) and rung 11 $0.0099$. The margin has grown at
each of the three rungs, and the rise continues: $0.5099$ to $0.5178$
(B3).

**B4 is refuted, and the refutation is this repository's own error
repeated.** The new rung's residual from the twelve-rung line is
$+0.0061$ against that line's in-sample r.m.s. $0.0037$, so B4 fails
as registered. But an out-of-sample point is not judged against an
in-sample scatter — [rem:primorialrung11] diagnosed exactly that for
its own P4 and printed the right width — and B4 was written against
the wrong one anyway. Against the prediction standard error at the new
abscissa, $0.0048$, the departure is still outside, at $1.27$
prediction standard errors: high, and not resolved.

**And the ladder has not bent, by the test that remark named.** On
thirteen rungs the slope is $+0.007301$ against $+0.007013$ and the
scatter is $0.0038$ against $0.0037$ — a rise of $0.0001$. "If the
ladder bent, the scatter would rise"; it did not. The line stands,
now with one point that was never fitted to it sitting $1.27$
prediction errors above.

No forecast is made from this. $\theta'$ lives at $\log_{10}N=10.6180$
on [rem:primorialdense]'s shape, two and a half decades past this
rung, and [rem:shapepower] is why nothing is published there.

*Added later.* A fourteenth rung, past the point where the shapes
part, is Remark [rem:rung13]. The line does not survive it as cleanly
as it survived this one.


#### Remark (the fourteenth rung is above both shapes) {#rem:rung13}
<!-- evidence: audit_primorial_rung13.py -->

$30030\cdot2^{13}=246005760$ sits at $\log_{10}N=8.3909$, past
[rem:shapetrust]'s parting point $8.3256$ — the first rung ever
measured where [rem:primorialdense]'s surviving shape and its
runner-up differ by more than the surviving fit's own r.m.s. Both are
recoverable exactly from the crossings that remark prints, so what the
new value had to be compared against was fixed before it was computed.
Rung 12 reproduces to $0.000035$ with the same $K^*_R=15461$ (C1).

**The margin escalates again** (C2). The exponent is $0.5283$, a
margin of $0.0283$ over $\tfrac12$ — $7.6$ times the ladder's scatter
$0.0037$, after $0.0023$, $0.0099$ and $0.0178$ at the three rungs
before it. Four rungs, four margins, each larger than the last.

**And the line does not hold it** (C3, refuted). The thirteen-rung
line predicts $0.5185$; the measurement is $0.5283$, a departure of
$+0.0098$, which against the prediction standard error $0.0048$ at
that abscissa is $2.02$. Both out-of-sample points now sit above the
line — $1.27$ prediction errors at rung 12, $2.02$ at rung 13 — and
the refit scatter, the test [rem:primorialrung11] named, has begun to
move: $+0.0001$ when rung 12 was added, $+0.0005$ now, from $0.0038$
to $0.0043$. **The ladder is bending upward.**

**The shapes do part there, and the point still cannot choose between
them** (C4 holds, C5 holds) — **because it is above both.**
$a+b\log N$ predicts $0.519492$ and $a+b\log\log N$ predicts
$0.515223$, a separation of $0.004269$ against the r.m.s. $0.003978$,
so the parting point is where it was said to be. But the measured
$0.5283$ is $0.008768$ from the first and $0.013037$ from the second;
the difference between those distances is $0.004269$, inside the
prediction error $0.004838$. The surviving shape is the nearer one and
that is all that may be said. **What the rung actually reports is that
both candidates underpredict it.**

This cuts one way and it is not the comfortable way. A ladder rising
faster than $a+b\log N$ would meet $\theta'$ *earlier* than
$10^{10.6180}$, not later. No number is put on that: two points cannot
fit a third shape, [rem:shapepower] measured what this repository's
shape discriminator can and cannot do, and the honest content here is
a direction and a warning that the published forecast rests on a curve
the last two measurements have both exceeded.

*Added later.* "The ladder is bending upward" is corrected to "the
ladder sits above the line": a third out-of-sample rung is above it
too, but by less, and the departures are not growing. See Remark
[rem:rung14].


#### Remark (three rungs above the line, none of them growing) {#rem:rung14}
<!-- evidence: audit_primorial_rung14.py -->

$30030\cdot2^{14}=492011520$ at $\log_{10}N=8.6920$ is the third rung
measured after the line was fitted, and it decides between the two
readings [rem:rung13] left open — a curve the shape does not have, or
a noisy pair. Rung 13 reproduces to $0.000040$ with the same
$K^*_R=27077$ (C1).

**The margin grows a fifth time** (C2). The exponent is $0.5333$, a
margin of $0.0333$ over $\tfrac12$ — $9.0$ times the ladder's scatter
— after $0.0023$, $0.0099$, $0.0178$, $0.0283$. That escalation is
the one thing in this thread that has not wavered.

**The rung is above the line again, and again above both shapes** (C3,
C5). The fourteen-rung line predicts $0.5264$ against a measured
$0.5333$, a departure of $+0.0069$; the candidate shapes predict
$0.524968$ and $0.519270$, and the measurement clears both by
$+0.008305$ and $+0.014003$. Three out-of-sample rungs, three above
the line, three above both shapes.

**But it is not a bend** (C4, refuted). In units of its own prediction
standard error the departure is $1.29$, against rung 13's $2.02$ and
rung 12's $1.27$. The departures are not growing; they are sitting.
C3's cap was one prediction standard error rather than two, so
"exceeds" there is a weak test and is reported as such. What the three
support is a roughly constant positive offset from the line the lower
rungs fit, not an accelerating curve — and [rem:rung13]'s "the ladder
is bending upward" is corrected to that.

The refits agree with the milder reading. Adding the three rungs took
the slope from $+0.007013$ to $+0.007301$, $+0.007703$, $+0.007953$
and the scatter from $0.0037$ to $0.0038$, $0.0043$, $0.0045$ —
steepening and loosening a little at each step, as a line does when
new points sit consistently above it, and not more than that.

Three of three above is $p=\tfrac18$ under a symmetric null and is
quoted for what it is worth, which is not much on its own. Still no
forecast is made: the published $\theta'$ crossing rests on a curve
that three measurements have now cleared from above, and whether that
means an earlier crossing or a slightly wrong constant is not
something five rungs past the fit can say.

*Added later.* "Not an accelerating curve" was measured on the wrong
quantity — departures from one line, rather than the rungs' own
slopes. Measured on those, the curvature is resolved. See Remark
[rem:laddercurve].


#### Remark (measured on its own slopes, the ladder curves) {#rem:laddercurve}
<!-- evidence: audit_ladder_curve.py -->

[rem:rung14] read three out-of-sample departures of $1.27$, $2.02$ and
$1.29$ prediction standard errors, saw no growth in them, and
concluded a constant offset rather than a curve. That conclusion is
measured on the wrong quantity. A departure is taken from a line
fitted on *all* the rungs, low ones included; if the ladder curves,
that line is already tilted by the curvature it is being used to
detect, and the prediction error grows with distance besides. The
question is about the rungs' own slopes, and fifteen of them are
published. Nothing new is computed here; the whole-set slope
reproduces to $0.000035$ (E1).

**The slope rises with $N$** (E2). On the lower seven rungs it is
$+0.006626\pm0.001225$; on the upper eight, $+0.010202\pm0.000732$.
The difference is $+0.003576$ against a combined standard error of
$0.001427$ — $2.51$ of them. The ladder is steeper at the top than at
the bottom by half again.

**A quadratic term is resolved** (E3). Fitted on all fifteen, the
coefficient of $(\log N)^2$ is $+0.00038561\pm0.00011737$, $t=3.29$,
and the r.m.s. residual falls from $0.0045$ to $0.0032$. **So the
ladder curves upward, and [rem:rung14]'s reading is corrected: what
looked like a fixed offset from a straight line is the straight line
failing.**

**Consistently, a line fitted on the upper rungs predicts the new ones
well** (E4). On the five rungs between the halfway mark and the
out-of-sample points the slope is $+0.008772$, and the three new rungs
depart from it by $0.44$, $1.05$ and $0.77$ of their own prediction
standard errors. One of three outside is what E4 asked for and it is
thin evidence; what carries the reading is E2 and E3. The two are not
in tension — a line fitted where the ladder is already steep absorbs
most of the curvature, which is exactly why the departures from the
all-rung line were not growing.

Three cycles have now given three readings of these points — a bend, a
flat offset, and a resolved curve — and the third is the one measured
on the ladder's own slopes rather than on residuals from a fit that
the curvature contaminates. The direction it implies is unchanged and
still not published as a number: a ladder steepening with $N$ meets
$\theta'$ earlier than $10^{10.6180}$, and extrapolating a quadratic
two decades past its data is the thing [rem:shapepower] exists to
forbid.


#### Remark (the shape contest could not have seen it) {#rem:curvereach}
<!-- evidence: audit_curve_reach.py -->

[rem:laddercurve] resolved an upward curvature the surviving shape
$a+b\log N$ does not have. That is either a defect in
[rem:primorialdense]'s contest or a limit of its reach, and the
difference is settled by refitting on prefixes — no new arithmetic.
The fifteen-rung quadratic reproduces to $5.0\cdot10^{-9}$ (F1).

**It is a limit, not a defect** (F2). On the eleven published rungs
the quadratic coefficient is $+0.00021493$ with a standard error of
$0.00030034$ — $t=0.72$. The curvature was not in the data that
contest had, and could not have been read from it.

**Two of the three new rungs sufficed** (F3, refuted). The sequence of
$t$ as rungs are added is $0.72$, $1.14$, $1.74$, $2.71$, $3.29$, so
it first clears two at fourteen rungs rather than fifteen. The
prediction that all three were needed is refuted; the third raised
$t$ from $2.71$ to $3.29$ and confirmed rather than established.

**F4 holds and is worth nothing as registered.** It asked whether
$a+b\log N$ is still the best of the shapes ranked by r.m.s., and it
is not — but the shape that beats it is the quadratic, which has three
parameters against the others' two, and a rule that ranks across
parameter counts is not a ranking. Among the four two-parameter shapes
the order is unchanged: $a+b\log N$ leads on both sets, at $0.003718$
on eleven rungs and $0.004465$ on fifteen, against
$0.004119\to0.006113$ for $a+b\log\log N$ and worse for the other two.
**What the numbers do say is in the comparison the rule should have
made:** the quadratic's edge over the best two-parameter shape grows
from $0.003605$ against $0.003718$ on eleven rungs — three per cent —
to $0.003240$ against $0.004465$ on fifteen, twenty-seven per cent.
That growth is the curvature arriving, and it is the third badly
written rule this programme has caught in its own scripts.

So [rem:primorialdense] stands as a correct reading of its own reach
and is superseded above it. Its $\theta'$ crossing was computed from a
shape that the rungs beyond its range now exclude; what replaces it is
not another number, for the reason [rem:shapepower] gives.


#### Remark (the curvature predicts, and has an expiry) {#rem:curvebound}
<!-- evidence: audit_curve_bound.py -->

A curvature that fits fifteen points is not yet a shape. The line was
held to out-of-sample prediction and failed it three times; the
quadratic had not been asked. It is asked here, and it has a ceiling
to declare besides: the level exponent is $\log K^*_R/\log N$ with
$K^*_R<N$, so it cannot pass $1$, and an upward quadratic must
therefore expire.

**It predicts** (G2). Fitted on the first thirteen rungs and asked for
the last two, the quadratic gives $0.5237$ and $0.5310$ against
measured $0.5283$ and $0.5333$ — departures of $+0.0046$ and $+0.0023$
against its own prediction standard errors $0.0054$ and $0.0063$, both
inside, and the second smaller than the first. The line on the same
thirteen gives $+0.0098$ and $+0.0097$: outside, and not shrinking.
**The curvature is the first shape in this thread to predict points it
was not fitted to.**

**The ceiling is nowhere near** (G3). At the four rungs that print it,
$K^*_R/N$ is $0.000152$, $0.000126$, $0.000110$, $0.000088$ and the
largest exponent is $0.5333$. Nothing measured is near the bound; it
constrains the extrapolation, not the data.

**And the crossing moves a decade earlier** (G4). The quadratic
reaches $0.56$ at $\log_{10}N=9.6068$ with bracket
$[9.3291,\,10.0373]$, against the line's published $10.6180$. Refitted
on the lower nine rungs it lands at $9.0380$, so the drift is
$0.5688$ — comparable to the bracket's own width, and declared as
such. The ceiling $1$ arrives at $18.5824$, bracket
$[16.6790,\,22.7350]$, drift $2.2830$ (G5): the curvature has a
computable expiry and it is far above anything at issue here.

**None of this is a forecast and the distinction is the point.** The
bracket is the fit's parameter spread and nothing more;
[rem:shapepower] measured that this repository's shape discriminator
has no power at this reach, and a quadratic carried past its data is
the case that warning was written about. What may be said is
comparative and does not need the numbers: the shape that predicts
out of sample puts $\theta'$ **earlier** than the shape that does not,
and [rem:primorialdense]'s $10^{10.6180}$ was computed from the
latter.


#### Remark (the sixteenth rung, and an extrapolation that settles) {#rem:rung15}
<!-- evidence: audit_primorial_rung15.py -->

[rem:curvebound] left the quadratic predicting but its crossing
unstable: a drift of $0.5688$ against a bracket of width $0.708$ is an
extrapolation that moves as much as it is worth.
$30030\cdot2^{15}=984023040$ at $\log_{10}N=8.9930$ is the point that
settles it. Rung 14 reproduces to $0.000027$ (H1).

**A sixth growing margin** (H2). The exponent is $0.5407$, a margin of
$0.0407$ over $\tfrac12$ — $11.0$ times the ladder's scatter, after
$0.0023$, $0.0099$, $0.0178$, $0.0283$, $0.0333$.

**The curvature predicts a third time, and better** (H3, H4). Fitted
on the fifteen published rungs it gives $0.5419$ against a measured
$0.5407$: a departure of $-0.0012$ against its own prediction standard
error $0.0049$, which is $0.26$ of it. The line on the same fifteen
gives $+0.0071$. **The quadratic's out-of-sample errors now read
$+0.0046$, $+0.0023$, $-0.0012$ — shrinking, and the last one has
changed sign** — against the line's $+0.0098$, $+0.0097$, $+0.0071$,
which do not.

**And the crossing settles** (H5). Refitted on sixteen rungs the
quadratic reaches $0.56$ at $\log_{10}N=9.6358$ with bracket
$[9.4216,\,9.9340]$; the fifteen-rung value was $9.6068$, so it moved
$0.0290$ against a declared drift of $0.5688$ — twenty times less than
the last addition moved it. The bracket has narrowed from $0.708$ to
$0.512$ and the r.m.s. from $0.0032$ to $0.0031$. **This is the first
extrapolation in this repository whose drift has come in an order of
magnitude below its own bracket.** *(Withdrawn by [rem:rung16]: the
next rung moved the crossing by $0.0425$ against the $0.0290$
declared here, so this was a pause and not a settling.)*

Two things follow and one does not. What follows: the crossing sits
$0.64$ decades above the top rung — a factor of $4.4$ in $N$, two more
doublings — so it is no longer an extrapolation into the far distance
but one just past the data, and the ladder could in principle be
carried to it rather than extrapolated to it. What also follows is the
limit: the sieve at $2^{16}$ needs sixteen gigabytes for $\Lambda$
alone and $2^{17}$ twice that, so the next doubling is at the memory
wall of this machine and the one after is past it. *(Withdrawn by
[rem:rung16]: the memory is not the wall and the k-cap is.)*

What does not follow is a forecast for $\theta'$. The quadratic has no
derivation, it is a local description with a computable expiry at
$\log_{10}N=18.5824$, and [rem:shapepower] measured that shape
selection at this reach has no power — a shape that predicts three
times running is better evidence than a shape that does not, and it is
still not a reason to quote where an underived curve meets a level.
What has changed is narrower and worth saying exactly: **the number
$10^{10.6180}$ came from a shape the rungs now exclude, and every
shape that survives them puts the crossing lower.**


#### Remark (the seventeenth rung, and the wall that was named wrong) {#rem:rung16}
<!-- evidence: audit_primorial_rung16.py -->

[rem:rung15] closed by naming the limit: "the sieve at $2^{16}$ needs
sixteen gigabytes for $\Lambda$ alone and $2^{17}$ twice that, so the
next doubling is at the memory wall of this machine and the one after
is past it." Both halves are wrong, and they are wrong in opposite
directions.

**The memory was never the wall.** The routine the ladder has used
since rung 11 holds four arrays of order $N$ where the information
fits in half that. Storing $\Lambda$ as the *prime* rather than its
logarithm halves it — the logarithm is taken at the point of use, of
the same float64 integer, so the float that enters the sum is the
identical float; building $\mu$ blockwise makes the int32 cofactor
block-local and the prime list below $N$ unnecessary; and the residue
mask is addressable only on odd indices, because $N$ is even and $k$
and $m$ are odd, so $N-mk$ is odd at every index the statistic reads.
Rung 16 costs $11.00$ GB of resident array against the $21.99$ GB the
old route would hold and the $7.33$ GB cofactor it builds them with.
Nothing moves: $\mu$ and $\Lambda$ agree with the production route at
every one of $2\cdot10^7$ indices, elementwise and bit for bit (C1),
and rung 14 and rung 15 return the published $K^*_R$ exactly, $43171$
and $72857$, at departures $0.000027$ and $0.000005$ (C2).

**The k-cap was the wall, and it binds exactly here.** At the
published cap rung 16 has no crossing below $k=10^5$: $K^*_R$ has
left the range the ladder searches. It was always going to. $K^*_R$
reads $43171$ at rung 14 and $72857$ at rung 15 against a cap that
does not move, so rung 15 was already at three quarters of it. **The
ladder stopped one rung short of its own definition running out, and
the remark attributed the stop to the wrong resource.**

**And the cap is inside the statistic, not beside it.** $\beta$ is
fitted by least squares over the same $k$-range the search runs on,
so widening the cap moves $\beta$, moves $R=H-\beta P$, and moves
$K^*_R$ at every rung already printed. The extension used here is the
one that leaves every published integer alone: $\beta$ keeps the
published window $k<10^5$ and only the truncation search widens, to
$k<4\cdot10^5$. C2 is what says that leaves the old rungs untouched.

**The rung.** $30030\cdot2^{16}=1968046080$ at $\log_{10}N=9.2940$
gives $K^*_R=122873$ and an exponent of $0.5476$.

**A seventh growing margin** (H1). The margin over $\tfrac12$ is
$0.0476$ — $12.9$ times the ladder's scatter — after $0.0023$,
$0.0099$, $0.0178$, $0.0283$, $0.0333$, $0.0407$.

**The curvature predicts a fourth time, at the first rung that could
choose** (H2, H3). Fitted on the sixteen published rungs the
quadratic gives $0.5499$ and the line $0.5408$: the two shapes are
$0.0091$ apart against a prediction standard error of $0.0046$,
$1.97$ of them. Every earlier rung had the shapes within about one
standard error of each other, so this is the first out-of-sample
point that discriminates. The measurement is $0.5476$: the
quadratic's departure is $-0.0023$, $0.50$ of its own prediction
error; the line's is $+0.0068$. **The quadratic's out-of-sample
errors now read $+0.0046$, $+0.0023$, $-0.0012$, $-0.0023$, all
inside their prediction errors, against the line's $+0.0098$,
$+0.0097$, $+0.0071$, $+0.0068$, none of which are.**

**The extrapolation did not settle; it paused** (H4, refuted).
[rem:rung15] read a drift of $0.0290$ and called it "the first
extrapolation in this repository whose drift has come in an order of
magnitude below its own bracket." Refitted on seventeen rungs the
$0.56$ crossing moves to $\log_{10}N=9.6783$ from $9.6358$ —
$0.0425$, which is **larger** than the $0.0290$ the last rung
declared. The bracket does narrow, to $[9.5061,\,9.9014]$ from
$[9.4216,\,9.9340]$; the centre does not hold. So that sentence is
withdrawn: one rung of small movement was a pause, and the test that
would have confirmed settling is the one that failed.

**And the convention has a cost that grows** (C3, refuted). Refitting
$\beta$ on the widened window instead of the published one moves the
exponent by $0.000824$ at rung 14 and $0.001157$ at rung 15. Both sit
well inside the ladder's scatter $0.0037$, which was the first half
of C3; the second half asked that the cost shrink as the rungs rise,
and it grows. The rule fixed before the run says what follows: the
frozen-window extension is withdrawn in favour of recomputing every
rung at one cap, which is [rem:laddercap]. What does not depend on
the convention is everything before the rung — the memory finding,
the cap finding, and C1 and C2.


#### Remark (the whole ladder at one cap, and what survives it) {#rem:laddercap}
<!-- evidence: audit_ladder_cap.py -->

[rem:rung16] left the ladder's exponents resting on an arbitrary
$k$-cap that turned out to be inside the statistic rather than beside
it, and its own C3 withdrew the frozen-window repair. So all
seventeen rungs are recomputed here with the $\beta$ fit and the
truncation search tied to the same cap, at four caps
$10^5,\,2\cdot10^5,\,4\cdot10^5,\,10^6$. One pass to $k<10^6$ at each
$N$ yields every cap as a subset operation, because $H(N;k)$ and
$P(N;k)$ do not depend on the cap at all — it enters only through
which $k$ are fitted and how far the sum is read.

**The control holds** (U1). At the published cap every rung whose
$K^*_R$ is printed returns it exactly: $9367$, $15461$, $27077$,
$43171$, $72857$.

**The first nine rungs cannot see the cap at all.** Rungs 0 through 8
return the identical $K^*_R$ at all four caps. The dependence begins
at rung 9 and grows with $N$, which is the same direction
[rem:rung16]'s C3 measured and the reason it was refuted.

**The cap converges, and U2 is refuted on a technicality worth
stating.** At rungs 12 through 15 the second step is smaller than the
first every time — $0.000332$ then $0.000124$, $0.000461$ then
$0.000045$, $0.000824$ then $0.000218$, $0.001157$ then $0.000406$.
U2 named the top five rungs, and the fifth is rung 16, which has no
value at the base cap at all; that is precisely [rem:rung16]'s
finding. So the prediction as written is not satisfied, and it is
recorded refuted rather than reworded after the fact. **What the four
measurable rungs say is that the cap has a limit and the published
value is an approximation to it, not a point on a drift.**

**And the approximation is inside the noise** (U3). The largest move
any rung makes between the published cap and $10^6$ is $0.001563$, at
rung 15 — $0.42$ of the ladder's scatter $0.0037$.

**The escalation survives** (U4). On the uniform ladder at $10^6$ the
top six margins over $\tfrac12$ are $0.0106$, $0.0183$, $0.0288$,
$0.0343$, $0.0423$, $0.0488$, growing at every step as they do at the
published cap.

**The curvature survives, and is better resolved than at the
published cap** (U5). On seventeen rungs the coefficient of
$(\log N)^2$ is $+0.00036885\pm0.00008050$, $t=4.58$, against the
$t=3.29$ [rem:laddercurve] reported on fifteen rungs at the published
cap. The r.m.s. residual is $0.0031$.

**The crossing survives, and tightens** (U6). The uniform ladder
reaches $0.56$ at $\log_{10}N=9.6200$ with bracket
$[9.4574,\,9.8277]$, inside the $[9.4216,\,9.9340]$ [rem:rung15]
published and narrower at both ends.

So the branch is not an artefact of the cap: every reading it rests
on is reproduced on a ladder whose cap was generous from the start,
and the one reading that improves is the curvature. **The rung-16
exponent on this ladder is $0.548808$**, against the $0.5476$
[rem:rung16] measured under the frozen window — a difference of
about a thousandth, well inside the scatter, and the uniform figure
is the one that is comparable with the other sixteen.

What still does not follow is a forecast for $\theta'$.
[rem:shapepower] measured that shape selection at this reach has no
power, and a quadratic that now predicts four times running and
resolves at $t=4.58$ is better evidence than one that does not
without being a reason to quote where an underived curve meets a
level. What is worth recording is where the next two rungs fall:
rung 17 sits at $\log_{10}N=9.5951$ and rung 18 at $9.8961$, and the
crossing estimated here is between them. **If the ladder reaches
them, $0.56$ stops being a level extrapolated to and becomes one
bracketed by the sign of two measurements** — which is what
[rem:primorialgap] did for $\tfrac12$, and the only kind of statement
[rem:shapepower] does not forbid.


#### Remark (the eighteenth rung, where the fit could not call it) {#rem:rung17}
<!-- evidence: audit_primorial_rung17.py -->

[rem:laddercap] left a uniform ladder of seventeen rungs whose
quadratic reaches $0.56$ at $\log_{10}N=9.6200$, and observed that
rung 17 sits at $9.5951$ — just under it. Fitted on those seventeen
the quadratic predicts $0.5592$ here against a prediction standard
error of $0.0044$: **the level is $0.17$ prediction errors away, so
the shape puts the rung on $0.56$ and cannot say which side.** That
is the situation in which a measurement is worth more than a fit, and
it is the only situation about $\theta'$ that [rem:shapepower] does
not forbid — [rem:primorialgap] did the same for $\tfrac12$.

**The packing, again.** Rung 17's three arrays under [rem:rung16]'s
whole-index route would be $21.99$ GB, and the machine's commit
headroom was below that. The same observation that had already
halved the residue mask halves the rest: $N$ is even, $k$ is coprime
to $N$ and hence odd, and $m$ is odd, so $N-mk$ is odd at every index
the statistic reads, and the even halves of $\Lambda$'s support and
of $\mu$ are never addressed. Everything is kept on half-indices
$v\mapsto v\gg1$ and the block sieve runs there, the odd multiples of
an odd $p$ being $p(2s+1)$ — an arithmetic progression of step $p$
from $(p-1)/2$. Three arrays, $12.83$ GB. Nothing moves: $\mu$ and
$\Lambda$ agree with the production route at every one of $10^7$ odd
indices, bit for bit (C1), and rungs 15 and 16 return the uniform
ladder's $K^*_R$ exactly, $75253$ and $126079$, with exponents
$0.542257$ and $0.548808$ to every digit printed (C2).

**The rung.** $30030\cdot2^{17}=3936092160$ at $\log_{10}N=9.5951$
gives $K^*_R=215843$ and an exponent of $0.555925$.

**An eighth growing margin** (H1). The margin over $\tfrac12$ is
$0.0559$, against rung 16's $0.0488$ — $15.1$ times the ladder's
scatter $0.0037$.

**The curvature predicts a fifth time** (H2, H3), and this is its
first out-of-sample test on the uniform ladder. The departure from
the quadratic is $-0.0033$, $0.75$ of its own prediction error; the
line, which predicts $0.5491$, misses by $+0.0068$. The two shapes
are $0.0101$ apart here, $2.29$ prediction standard errors — the
sharpest discrimination the ladder has offered, against $1.97$ at
rung 16.

**And $0.56$ is not yet crossed** (H4). The exponent is $0.555925$,
short of the level by $0.004075$ against the ladder's scatter
$0.0037$. **So the answer is resolved, and it is resolved by a tenth
of the floor.** The registered rule named the other outcome as the
one worth having — an exponent reaching $0.56$ would have put the
crossing at or below $\log_{10}N=9.5951$ as an observation rather
than a fit — and the measurement went the other way. What it buys is
narrower and still not a forecast: **the crossing of $0.56$ is now
bounded below by a measurement**, at a rung the fit was not entitled
to call.

**The refit, and what it does not settle.** On eighteen rungs the
$(\log N)^2$ coefficient is $+0.00033853$ with r.m.s. residual
$0.0031$, and the $0.56$ crossing moves to $\log_{10}N=9.6667$ from
$9.6200$ — a drift of $0.0467$, the second in a row larger than the
last rung's own declared drift. The bracket $[9.6009,\,9.8386]$ is
narrower than before at both ends, but the centre keeps receding by
about the same amount each time a rung is added, and nothing here
resolves whether that is convergence from below or a systematic.
**The bracket is also censored**: $662$ of $4000$ draws put the
crossing at or below the top rung and are dropped, so the fit itself
gives about a sixth of its own weight to "already crossed" while the
measurement at that rung says otherwise. The interval quoted is
conditioned on the rest and should be read as such.

Rung 18 sits at $\log_{10}N=9.8961$, above the upper end of that
bracket. If it can be reached, $0.56$ is bracketed by the sign of two
measurements rather than by any shape — which is the whole of what
this branch can deliver, and it does not touch the sign axis, where
[rem:leanidentity] and [rem:destination] leave the actual obstruction.


#### Remark (the nineteenth rung: 0.56 crossed, and inside the floor) {#rem:rung18}
<!-- evidence: audit_primorial_rung18.py -->

[rem:rung17] left rung 17 at $0.555925$, below $0.56$ by $0.004075$
against a floor of $0.0037$, and the eighteen-rung crossing at
$\log_{10}N=9.6667$. Rung 18 sits at $9.8961$, past the upper end of
that bracket, and the two surviving shapes disagreed about it: fitted
on the eighteen uniform rungs the quadratic predicts $0.5668$ here
and the line $0.5565$, $0.0103$ apart — $2.40$ prediction standard
errors — so one says the level is cleared and the other says it is
not. **Every earlier rung asked the shapes for a number; this one
asked them a yes-or-no question.**

**Neither the prime nor the ceiling.** $N=7872184320$ exceeds $2^{32}$,
so [rem:rung17]'s uint32 could not carry the prime that carries
$\Lambda$, and its three half-index arrays would be $25.66$ GB here.
Both problems are the same problem: the prime is stored at all. On
the support of $\Lambda$ every $v$ is a prime power, and **when $v$ is
prime the prime is $v$** — nothing to store. What remains is the
powers $p^j$ with $j\ge2$, and those have $p\le\sqrt N$, so the whole
range holds $9023$ odd ones: a table, not an array. One byte per
half-index says which of the three cases $v$ is in and the logarithm
is taken at the point of use, of the same float64 integer. The mask
sheds its second byte for a reason particular to this ladder: $k$ is
coprime to $N=30030\cdot2^{18}$, so none of $3,5,7,11,13$ can divide
$k$ and their five bits are always required-zero together — one bit,
"coprime to $15015$" — leaving only $17,19,23,29$ needing their own.
Nine bits become five. Three arrays of one byte, $11.00$ GB.

Nothing moves. $\mu$ and $\Lambda$ agree with the production route at
every one of $10^7$ odd indices, bit for bit, with the support split
$1270606$ prime and $709$ prime power against a table of $709$ (C1);
the five-bit mask returns the nine-bit mask's keep decision at every
odd $v$ for all $184$ admissible $k$ below $1000$, with zero
disagreements (C2); and rungs 16 and 17 return $126079$ and $215843$
with exponents $0.548808$ and $0.555925$ to every digit printed (C3).

**The rung.** $30030\cdot2^{18}=7872184320$ at $\log_{10}N=9.8961$
gives $K^*_R=370859$ and an exponent of $0.562768$.

**A ninth growing margin** (H1). The margin over $\tfrac12$ is
$0.0628$ against rung 17's $0.0559$.

**A sixth out-of-sample hit, and the narrowest** (H2, H3). The
departure from the quadratic is $-0.0041$, $0.95$ of its own
prediction error — inside, but only just, and the closest any of the
six has come to failing. The line misses by $+0.0062$.

**$0.56$ is crossed** (H4) **and the crossing is inside the floor.**
The exponent is over the level by $0.002768$, against the ladder's
scatter $0.0037$: the result file declares `INSIDE FLOOR` on the same
line that declares the crossing. So the two facts must be said
together and in this order. The level is between rungs 17 and 18 —
$(10^{9.5951},\,10^{9.8961}]$, $0.3010$ decades — and that bracket is
made by the sign of two measurements rather than by any fitted shape,
which is the thing [rem:shapepower] does not forbid. But **only its
lower end is resolved**: rung 17 falls short by $1.10$ floors and
rung 18 clears by $0.75$ of one. By the standard this repository
adopted in G37 and G40, and applied to itself at
[rem:primorialrung10], a single point clearing a barrier by less than
the scatter has not cleared it. **This is a point-estimate bracket
with a resolved floor and an unresolved ceiling**, exactly as
[rem:primorialgap] had to say of its own.

What settles the upper end is rung 19, or a second $N$ of the same
odd radical placed between these two, which is what
[rem:primorialgap] did for $\tfrac12$ — and what it found there is
worth remembering: the interval closed by confirmation, and the
resolution of the individual points was worse than the bracket
suggested.

**And the departures have a sign.** The quadratic's six out-of-sample
errors read $+0.0046$, $+0.0023$, $-0.0012$, $-0.0023$, $-0.0033$,
$-0.0041$ — four negative in a row, growing in size each time, and
the last at $0.95$ of its prediction error. Each is individually
inside, so no registered test has failed; but a one-signed run that
grows is what a shape looks like when it begins to fail from above,
and it would mean the ladder is flattening against the quadratic and
the crossing of any level above $0.56$ is later than the quadratic
places it. Nothing here measures that. It is recorded as open.
*(Withdrawn by [rem:signrun], which measured it: the six were spliced
from two ladders, the run on one ladder is ordinary under the null at
$0.0636$, and the growth is unresolved at $t=-0.80$.)*


#### Remark (the run is real, ordinary, and not growing) {#rem:signrun}
<!-- evidence: audit_ladder_signrun.py -->

[rem:rung18] closed by pointing at the sign of the quadratic's
out-of-sample errors — four negative in a row and growing — and said
that is what a shape looks like when it begins to fail from above.
That paragraph is withdrawn. Two of its three components do not
survive being measured, and the one that does means less than it was
made to mean.

**First, the six numbers were not one series.** Four were computed on
the ladder at the published cap and two on the uniform ladder of
[rem:laddercap]. So the series is rebuilt here on the uniform ladder
alone by walking forward — for each $j$, fit the quadratic on rungs
$0..j-1$ only and predict rung $j$ — which is the construction the
rung remarks used at the top, applied everywhere it can be. It
reproduces those two rungs exactly: departure $-0.0033$ with error
$0.0044$ at rung 17 and $-0.0041$ with $0.0043$ at rung 18, to the
decimals they print (S1).

**The run is real, and longer than reported** (S2). Thirteen
walk-forward points end in a run of $5$, not four: rungs 14 through
18 are all negative.

**And it is ordinary** (S3, refuted). Under the null — the fitted
quadratic at the same abscissae with i.i.d. Gaussian errors at the
ladder's own r.m.s. $0.0031$, walked forward by the same code —
$1272$ of $20000$ ladders end in a run of $5$ or longer. That is
$0.0636$, above the $0.05$ the rule fixed. **A correct shape with
this much noise and this few points produces a run this long about
one time in sixteen, so the run is not evidence that the shape is
failing.** This was registered as the prediction most likely to be
refuted, and it was.

**And it is not growing** (S4, refuted). Regressed on $\log N$ the
departures have slope $-0.00051419\pm0.00063930$, $t=-0.80$. The sign
is not resolved. So "growing" was never established; only
one-signedness was, and that is the weaker statement.

The full record is worth stating too, because publishing only the top
six flattered it. Across all thirteen walk-forward points the ratio
of departure to its own prediction error reads $1.55$, $1.41$,
$0.70$, $1.09$, $0.05$, $0.35$, $0.38$, $0.81$, $0.07$, $0.16$,
$0.58$, $0.75$, $0.95$ — **three exceed one, all three at the bottom
of the ladder**, and the largest departures in absolute terms
($+0.0116$, $-0.0111$) are the earliest, not the latest. The
quadratic's out-of-sample behaviour is not deteriorating with $N$;
by this measure it is improving.

What this closes and what it does not. It closes the alarm: nothing
in the departures says the ladder is flattening against the
quadratic, so [rem:rung18]'s inference that levels above $0.56$ are
crossed later than the quadratic places them has no support and is
withdrawn. It does not promote the quadratic. [rem:shapepower]
measured that shape selection at this reach has no power, and a run
that is ordinary under the null is equally ordinary under shapes the
null was not drawn from. **The departures are simply not informative
about the shape, in either direction, and that is the finding.**


#### Remark (the floor was the wrong one: 0.56 resolved) {#rem:rung18fill}
<!-- evidence: audit_rung18_fill.py -->

[rem:rung18] crossed $0.56$ and declared `INSIDE FLOOR` on the same
line, because the clearance $0.002768$ sits below the ladder's
scatter $0.0037$. That verdict was correct under the rule it was
judged by. What it could not say is whether $0.0037$ is the right
yardstick for the question, and it is not.

**What the floor measures.** $K^*_R$ is a deterministic integer once
$N$ is fixed; there is no sampling noise in it. The ladder's scatter
is the spread of its exponents about a curve fitted across ten
decades, so it mixes arithmetic fluctuation with the shape's own
misfit. The question `INSIDE FLOOR` raises is narrower and local:
*could a neighbouring $N$ of the same odd radical give an exponent
low enough to fall under the level?* That is answered by measuring
neighbouring $N$, which is what [rem:primorialgap] did for
$\tfrac12$.

Four are placed here within $0.0208$ decades of rung 18 on either
side, each $30030m$ with $m$ composed only of $2,3,5,7,11,13$ so the
odd radical, and with it the threshold, is the ladder's. The window
spans $0.0413$ decades, over which the ladder's slope contributes far
less than the floor. Rung 16 recomputes to $126079$ and $0.548808$
through this harness (C1), and the threshold reads $0.087306$ at all
five $N$ with a spread of $0.00000049$, inside the printing bound
(C2) — which is what "the same odd radical" has to mean.

**Every one of them clears the level** (F1). The exponents are
$0.563206$, $0.562560$, $0.563119$, $0.563611$, against rung 18's
$0.562768$, clearing $0.56$ by $0.003206$, $0.002560$, $0.003119$ and
$0.003611$.

**And the local floor is twelve times smaller than the ladder's**
(F2). A line through the five has slope $+0.005751$ per log unit and
an r.m.s. residual of $0.000308$. **So the arithmetic fluctuation
between neighbouring $N$ is not $0.0037$; across this window it is
$0.000308$, and the ladder's scatter at that scale is almost entirely
the shape's misfit rather than the arithmetic's noise.**

**The crossing is resolved** (F3). The lowest of the five clears
$0.56$ by $0.002560$ against a local floor of $0.000308$ — a ratio of
$8.30$. Nothing in this window comes near falling under the level.

Two things follow, and they should be kept apart. The first is the
one that was asked: **the upper end of the $0.56$ bracket is
resolved**, and [rem:rung18]'s "point-estimate bracket with an
unresolved ceiling" is closed — not by a bigger margin at one point
but by the floor being measured where the question lives. The second
is unasked and larger: **every single-point verdict in this branch
judged against the ladder's $0.0037$ has been judged against a floor
that is too big at local scale**, including the $\tfrac12$ crossings
of [rem:primorialrung10] and [rem:primorialgap]. This does not
overturn them — a conservative floor makes a crossing harder to
claim, not easier — but it means the branch has been understating
what its own measurements resolve, and nothing has measured the local
floor anywhere else on the ladder.

The bracket also narrows on its own. The lowest of the four new
points sits at $\log_{10}N=9.8756$ and is already above the level, so
the crossing is below it: with rung 17 at $9.5951$ still short, $0.56$
is now located in $(10^{9.5951},\,10^{9.8756}]$ by the sign of two
measurements. **No shape is used to say it.** What remains open there
is the interval itself, which is wide and unfilled, and
[rem:shapepower] still forbids quoting a crossing inside it.


#### Remark (the floor is scale-dependent, and 1/2 is crossed more than once) {#rem:localfloor}
<!-- evidence: audit_local_floor.py -->

[rem:rung18fill] measured the arithmetic fluctuation between
neighbouring $N$ near rung 18 at $0.000308$ and observed that every
single-point verdict in this branch had been judged against the
ladder's $0.0037$ instead. Two of those verdicts are at $\tfrac12$
and both retreated. This measures the floor where they live, the same
way, and the answer changes three things.

Rungs 9 and 10 recompute to $3551$ with $0.494008$ and $5779$ with
$0.502394$ (C1) and the threshold reads $0.087306$ at all fifteen $N$
with zero spread (C2).

**The floor is not one number** (L1). Six same-radical $N$ around rung
10 give a local floor of $0.001883$; eight around
[rem:primorialgap]'s interval give $0.001662$. Both are under the
ladder's $0.0037$, so L1 holds — **but both are about six times the
$0.000308$ measured near rung 18.** The fluctuation shrinks with $N$,
as it should, and that means no floor may be carried from one part of
the ladder to another. [rem:rung18fill]'s number licenses nothing
outside its own window, and neither does this one.

**Rung 10's retreat was unnecessary, and only just** (L2). It clears
$\tfrac12$ by $0.002394$ against a local floor of $0.001883$, a ratio
of $1.27$. So [rem:primorialrung10]'s crossing survives being judged
where it lives — by a quarter, not by an order of magnitude, and it
would not have survived a floor twice as careful.

**One end of the gap bracket resolves and the other does not** (L3
hold, L4 refuted). At $N=20180160$ the exponent clears $\tfrac12$ by
$0.002978$, $1.79$ floors. At $N=18018000$ — the point
[rem:primorialgap] named "last short" — it falls short by
$0.000194$, which is $0.12$ of the floor. That was registered as the
prediction likely to fail and it failed.

**And the reason it failed is the finding.** Eight points across
$0.0582$ decades read, against $\tfrac12$: $-0.000194$, $-0.004142$,
$+0.000511$, $-0.000204$, $+0.001345$, $-0.001365$, $+0.002978$,
$+0.000886$. **The sign changes five times.** Six of the eight are
inside the local floor, so their side of $\tfrac12$ is not
determined at all; the two that are resolved are $18498480$ below and
$20180160$ above.

So [rem:primorialgap]'s reading has to be corrected, though not its
arithmetic. Its four interior points read $0.4966$, $0.4998$,
$0.5030$, $0.5040$ and looked monotone, and it concluded "the
crossing is *inside* the interval". With eight points in the same
place the exponent does not rise through $\tfrac12$ once; **it sits
on $\tfrac12$ within the fluctuation and changes sign repeatedly.**
There is no single crossing to bracket here, and the interval that
was reported as one is better described as the region where the
ladder is indistinguishable from $\tfrac12$. The four points that
suggested otherwise were a sparse sample of an oscillation.

What survives, stated as narrowly as it should be: **the exponent is
resolved below $\tfrac12$ at $10^{7.2671}$ and resolved above it at
$10^{7.3049}$**, so the ladder does pass from one side to the other
across that span; nothing measured says it does so once. And rung 10,
higher up at $10^{7.4879}$, is resolved above.

Two cautions follow for the level this branch actually cares about.
[rem:rung18fill] found five consecutive $N$ all above $0.56$ by eight
to twelve local floors, with no sign change among them, so nothing
here disturbs that. But the span from rung 17 to those points is
unfilled, and **the oscillation seen at $\tfrac12$ is a reason not to
assume a single crossing of $0.56$ either.** The bracket
$(10^{9.5951},\,10^{9.8756}]$ remains what it was said to be: two
resolved measurements with an unmeasured interior.


#### Remark (0.56 is crossed once, and the band says why 1/2 was not) {#rem:targetband}
<!-- evidence: audit_target_band.py -->

[rem:localfloor] found the sign of exponent $-\tfrac12$ turning over
five times across eight same-radical $N$ and read it as the ladder
sitting on the level rather than passing through it. That reading was
right about those points and incomplete about the reason, and the
missing piece decides where to look for any level.

**A point cannot be told from a level when $|{\rm exponent}-{\rm
level}|$ is under the local floor, so the indeterminate band is about
$2\,{\rm floor}/{\rm slope}$ wide.** At $\tfrac12$ that window's slope
was $0.024973$ per log unit against a floor of $0.001662$: a band of
$0.0578$ decades, and the window spanned $0.0582$. **The window was
the band.** Every point in it had to be indeterminate, and five sign
changes is what sampling a band from inside looks like. Nothing was
learned there about how the ladder meets $\tfrac12$ because nothing
could be.

So the same question is put to $0.56$ where the band can be crossed.
Seven same-radical $N$ span $0.0641$ decades about $10^{9.7521}$,
where the chord between rung 17 and the lowest point of
[rem:rung18fill] puts the level, spaced about $0.0107$. Rung 16
recomputes to $126079$ and $0.548808$ (C1) and the threshold reads
$0.087306$ at all eight $N$ with zero spread (C2).

**The floor keeps shrinking** (W1). The local floor here is
$0.000264$, against $0.001662$ at $\tfrac12$ and $0.000308$ near rung
18. Arithmetic fluctuation falls with $N$ across three measured
places now, and no floor may still be carried from one to another.

**The window straddles the level** (W2): four of the seven are below
it and three above.

**And the sign turns over exactly once** (W3). In $N$ order the
exponents read $0.559068$, $0.558602$, $0.559720$, $0.559837$,
$0.560264$, $0.560247$, $0.560576$ — signs $-\,-\,-\,-\,+\,+\,+$
against the level. Three of the seven sit inside the local floor, so
their side is undetermined, but they are contiguous and they surround
the turn. **That is what a crossing looks like**, and it is not what
$\tfrac12$ looked like.

**Because the band here is a third as wide** (W4). $2\,{\rm
floor}/{\rm slope}$ is $0.0190$ decades against $0.0578$ at
$\tfrac12$: the floor is six times smaller and the slope only twice,
so the window at $0.0641$ decades covers the band instead of sitting
inside it. The difference between the two levels is not that the
ladder behaves differently at them; it is that one was sampled across
its band and the other inside it.

**The bracket narrows by a factor of six, and still without a shape.**
The highest $N$ resolved below the level is $5516751240$ at
$10^{9.7417}$ with exponent $0.559720$, and the lowest resolved above
is $6090084000$ at $10^{9.7846}$ with $0.560576$. So
$$0.56 \text{ is crossed in } (10^{9.7417},\,10^{9.7846}],$$
a span of $0.0429$ decades, against the $(10^{9.5951},\,10^{9.8756}]$
[rem:rung18fill] left. Both ends are resolved measurements and no
fitted curve is used to say it.

**The upper end is one rounding from being tighter, and that should be
said.** The point at $10^{9.7632}$ clears the level by $0.000264$ and
the local floor is $0.000264$ — equal in every digit printed. The
comparison is strict, so that point counts as undetermined and the
resolved end falls back to $10^{9.7846}$; the other way it would be
$10^{9.7632}$. Nothing here is wrong, but a verdict that turns on
digits past the printing is one to quote with its fragility attached.

Two things this does not do. It does not make $\theta'$ quotable:
[rem:shapepower] forbids reading a level's crossing off an underived
shape, and nothing here is a shape — but nothing here is a forecast
either, which is the point. And it does not revive
[rem:primorialgap]: that window was inside its band, so it could not
have resolved a crossing whatever the ladder did, and
[rem:localfloor]'s correction of it stands with a reason attached
rather than a warning.

What it leaves is a rule the branch did not have. **To resolve a
crossing, sample a window wider than $2\,{\rm floor}/{\rm slope}$ at
that place** — and since the floor falls with $N$ while the slope
does not, that gets easier, not harder, the further the ladder is
carried.


#### Remark (the floor falls like a power, and 1/2 was understated by twelve) {#rem:floorlaw}
<!-- evidence: audit_floor_law.py -->

[rem:localfloor] found the local floor differing by about six between
$2\cdot10^7$ and $8\cdot10^9$ and concluded that none of them may be
carried from one place to another. Four have now been measured —
$0.001883$ around rung 10 and $0.001662$ across
[rem:primorialgap]'s interval, both near $10^{7.3}$; $0.000264$ across
the $0.56$ crossing and $0.000308$ around rung 18, both near
$10^{9.9}$. Nothing is measured in this remark; every value is read
from a whole marker line and declared with a `READ` line that G76
checks against its source, and the windows' abscissae are imported
from the scripts that defined them rather than parsed back out of
prose. Both halves of that are the lesson from the failure G76 was
built for.

**The floor falls as a power** (P4). Regressing $\log$ floor on
$\log N$ over the four windows gives a slope of
$-0.321639\pm0.033495$, nominally $t=-9.60$ over a span of $6.0098$
in $\log N$. **The $t$ should not be read as nine standard errors.**
The four windows sit in two clusters $2.6$ decades apart, so the
slope is in substance a two-point determination and the small
within-cluster scatter flatters the standard error; the script says
so on the line beneath it. What the two clusters do support is the
direction and the rough size: the floor falls roughly like
$N^{-1/3}$, and a window at ten times the $N$ buys about half the
floor.

That is worth having because it sizes a window in advance.
[rem:targetband] showed that a window narrower than
$2\,{\rm floor}/{\rm slope}$ cannot resolve a crossing however many
points it holds, and [rem:primorialgap] fell into exactly that. With
a floor that is predictable to a factor, the band can be estimated
before the compute is spent rather than after.

**And the branch's headline is understated by twelve** (P2, P3). The
target section has been reporting the crossing of $\tfrac12$ as
[rem:primorialrung11] left it — a margin clearing the ladder-wide
floor by a factor under three. At rung 18 the margin over $\tfrac12$
is $0.062800$. Against the ladder-wide $0.0037$ that is $17.0$
floors; **against the floor measured in rung 18's own window,
$0.000308$, it is $203.9$.** Rung 17's $0.0559$ is $211.7$ of the
floor measured next to it. The ladder-wide number was never wrong; it
was answering a different question, and it cost this branch an order
of magnitude in what it could say.

Rung 10 is the other end of that. Its margin $0.002394$ clears its
own window's floor $0.001883$ by $1.3$ — resolved, as
[rem:localfloor] found, and barely. **So the crossing of $\tfrac12$
is resolved at both ends of the measured ladder, by a factor of
$1.3$ at $10^{7.5}$ and by two hundred at $10^{9.9}$**, and the
growth between them is the branch's real evidence rather than any
single rung's clearance.

What this does not do is touch $\theta'$. The exponent crossing
$\tfrac12$ with room to spare is what the level axis was always
measuring, and [rem:shapepower] still forbids reading a value of
$\theta'$ off any shape fitted to it. Nothing here is a shape. The
statement is narrower and now quantitative: **at the top of the
ladder the level exponent stands two hundred local floors above
$\tfrac12$**, and the yardstick that said seventeen was measuring the
ladder's departure from a fitted curve rather than the arithmetic's
own noise.


#### Remark (the floor law predicts a window it was not fitted on) {#rem:floormid}
<!-- evidence: audit_floor_midband.py -->

[rem:floorlaw] fitted $\log$ floor on $\log N$ over four windows and
got $-0.321639$, and said in the same breath why its $t=-9.60$ should
not be read as nine standard errors: the four sat in two clusters
$2.6$ decades apart, so the slope was in substance a two-point
determination. A line through two clusters is a line, not a law, and
the way to tell them apart is a third place.

$10^{8.5}$ is almost exactly between them — between rung 13 at
$10^{8.3909}$ and rung 14 at $10^{8.6920}$ — and at $N$ near
$3\cdot10^8$ a window costs minutes. **So the law was asked for the
number before the number was measured.** Fitted on the four windows
alone it puts the floor at this window's geometric mean
$N=316799937$ at $0.000767$, and the tolerance was fixed at twice the
four-window fit's own r.m.s. residual, $0.2668$ in log — a factor of
$1.31$ — because that r.m.s. is also the scatter between two floors
measured at the same $N$.

Rung 13 recomputes to $27343$ and $0.528766$ (Q1) and the threshold
reads $0.087306$ at all eight $N$ with zero spread (Q2).

**The prediction holds** (Q3). Seven same-radical $N$ across
$0.0389$ decades give a local floor of $0.000897$ against the
predicted $0.000767$: a log ratio of $+0.1569$, a factor of $1.17$,
inside the registered $1.31$. It lies between the clusters' own
floors, $0.000308$ below and $0.001662$ above (Q4).

**And the slope survives the clustering being broken** (Q5). On five
windows it is $-0.322586\pm0.030876$, $t=-10.45$, against
$-0.321639$ on four; the r.m.s. residual moves only from $0.1334$ to
$0.1348$, which is another way of saying the new window fell on the
line rather than being absorbed by it. The $t$ is now read off a
spread that is no longer two points, which is the first time this
slope has meant what a $t$ ordinarily means — on five windows, which
is still few.

What this buys is not a fact about the ladder. It is the ability to
size a window before spending the compute on it. [rem:targetband]
showed that a window narrower than $2\,{\rm floor}/{\rm slope}$
cannot resolve a crossing however many points it holds, and that
[rem:primorialgap] had fallen into exactly that at $\tfrac12$ — its
window *was* its band, so five sign changes were the only thing it
could have seen. With a floor predictable to about a sixth, the band
at any $N$ can now be estimated in advance and the window chosen to
straddle it. **The branch stops discovering after the fact that it
sampled from inside.**

Two things it does not do. It does not extend to $N$ far outside
$[10^{7.3},\,10^{9.9}]$: five windows over six in $\log N$ constrain
a power there and say nothing about where the power comes from, and
no mechanism is offered for why $N^{-1/3}$ rather than anything else.
And it does not touch $\theta'$ — the floor is the resolution of a
measurement, not the thing measured, and [rem:shapepower] is
untouched by knowing how well one can see.


#### Remark (the ladder does not escalate, and that is not a vindication) {#rem:ladderdegree}
<!-- evidence: audit_ladder_degree.py -->

[rem:deficitregion] found the sign axis resolving a new polynomial
coefficient at almost every degree it was offered while the residual
barely moved, and read that as a family chasing a shape it does not
contain. The level axis rests on a quadratic — [rem:laddercurve]
resolved it at $t=3.29$, [rem:laddercap] at $t=4.58$, and three rung
remarks each read a rung as that curvature predicting out of sample.
So the same question is put to it. Nothing is measured; the nineteen
uniform rungs are read from three result files, and the seventeen
[rem:laddercap] fitted reproduce its $t=4.58$ exactly (D1).

**The ladder does not escalate** (D2, D3). On all nineteen rungs the
cubic coefficient is $-0.0001035913\pm0.0001091288$, $t=-0.95$, and
across degrees two to six the top coefficient reads $+5.13$, $-0.95$,
$-1.00$, $-0.15$, $+0.01$: **not one coefficient past the second
resolves.** Where the sign axis resolved five of seven, the ladder
resolves none. The quadratic also strengthens with the two rungs
added since [rem:laddercap], from $t=4.58$ to $t=5.13$.

**And that is not the quadratic being vindicated** (D4, refuted). The
r.m.s. residual runs $0.003070$ at degree two to $0.002878$ at degree
six — a fall of $0.0625$, against the $0.08$ the sign axis conceded
over the same span of added parameters. **The ladder buys even less
per parameter than the deficit did.** A family whose extra degrees
are silent *and* buy nothing is not a family the data has chosen
from; it is a family the data cannot argue with. Nineteen points
resolve fewer degrees than a hundred and fifty-six whatever the
shape, and D4 was registered precisely to keep that reading available
against D2's silence.

So the two axes fail differently and arrive together. The sign axis
has the points to challenge its family and they refute it; the level
axis has a family nothing challenges because there is not enough of
it. **Neither licenses a statement past its own field**, and the
level axis's licence was never the weaker of the two — it was the
less tested.

What this does not touch is what the rung remarks actually measured.
[rem:rung16], [rem:rung17] and [rem:rung18] each recorded an
out-of-sample prediction and each one landed; [rem:signrun] measured
that the run of departures is ordinary under a correct-shape null,
which is a fact about the departures and stands. A term that predicts
six times is a term that predicts six times. What this remark denies
is the step from there to *the shape is the quadratic* — the data has
never been in a position to refuse a different one, and now the size
of that inability is on the record.

The consequence for $\theta'$ is the one [rem:shapepower] already
carries and it is unchanged: no crossing is read off this curve. The
consequence for what to do next is new. **Adding rungs will not
settle the shape** — the sign axis shows what more points buy, and
what they bought there was higher degrees resolving, not a lower one
confirmed. The ladder's own residual, falling six per cent across
four extra parameters, says the same about its future.


#### Remark (the void between two rungs, filled) {#rem:primorialgap}
<!-- evidence: audit_primorial_gap.py -->

Remark [rem:primorialreach] left "a crossing anywhere in
$[10^{7.19},\,10^{7.36}]$" open, and nothing since has been measured
in it. The ladder is $N=30030\cdot2^j$, so it steps from rung 9 at
$10^{7.1868}$ straight to rung 10 at $10^{7.4879}$ and the open
interval sits in the void between them. **It was never excluded; it
was stepped over.**

Doubling is not what the ladder holds fixed. What it holds fixed is
the prime set $P(N)=\{2,3,5,7,11,13\}$ — that is what makes
$\SS(N)(1-A(N))$ constant along it, and it is the single odd radical
whose caveat Remark [rem:arithmeticreach] says cannot be computed
away. Any $N=30030m$ with $m$ composed of those primes holds it just
as well, and the interval contains more than twenty. Four are taken,
spread across it, and the controls say they are on the ladder: the
same code returns the published rung-9 exponent to $0.000010$ and the
published fitted column to $0.000027$ (Q1), and the threshold reads
$0.087306$ at all five, to every digit printed, with $18863$ admissible
$k$ at each (Q4).

**Q2 is refuted, and it is the refutation the rule named as the one
worth having.** The four exponents are $0.4966,\,0.4998,\,0.5030,\,
0.5040$: the crossing is *inside* the interval, not above it as the
eleven-rung line said. *(Corrected by [rem:localfloor]: eight points
in the same place change sign five times, so these four are a sparse
sample of an oscillation and there is no single crossing to bracket.)* So [rem:primorialreach]'s interval closes by
**confirmation**. Between the last point short of $\tfrac12$ and the
first point over it the crossing is bracketed to
$(10^{7.2557},\,10^{7.3049}]$, $0.0492$ decades against the $0.3010$
that separated the two rungs — a narrowing of $6.12$, and one made by
the sign of a measurement rather than by a fit.

**Q6 is refuted too, and it bounds how much of that may be said.**
The registered expectation was that nothing in the interval would
resolve against $\tfrac12$. One point does: $0.5040$ clears the
ladder's floor $0.0037$ by $0.0040$. But only one. The other three sit
inside the floor — $0.0034$, $0.0002$ and $0.0030$ — so **the bracket
above is a point-estimate bracket and not a resolved one**, and its
lower end in particular rests on a point $0.0002$ from $\tfrac12$.
What is resolved is a weaker and still new statement: the exponent is
observably over $\tfrac12$ already at $10^{7.3526}$, a third of a
decade below the rung where [rem:primorialrung10] first saw it.

**Q3 is refuted, and the r.m.s. is the least interesting part of it.**
The four residuals about the eleven-rung line are $+0.0001,\,+0.0025,\,
+0.0049,\,+0.0052$ — four of four positive and monotone increasing.
They do not scatter about the line, they leave it. The r.m.s.
$0.0038$ against the published $0.0037$ is a hair over and would be
easy to dismiss; the sign run is not, and no check in this repository
was reading it.

What that does *not* license is a claim that the ladder bends.
Refitting on all fifteen points **lowers** the scatter, $0.0037\to
0.0035$, with the slope moving $+0.006778\to+0.007139$ — by
[rem:primorialrung11]'s own test, a ladder that bent would have raised
it. The systematic residual is therefore a **slope** correction and
not a curvature: eleven rungs put the slope too shallow, and the
points in the gap tilt it. The local slope between rung 9 and the top
interior point is $+0.025989$ per unit $\log N$ against the ladder's
$+0.006778$, a factor $3.83$ — but it is a two-point slope over a short
lever, and propagating the ladder's own scatter through that lever puts
it under two standard errors of itself, so the factor is **not
resolved** and the script declares it so. Remark [rem:slopes] is the
standing warning against reading it. The honest version is the one the
published rungs already contained: rung 9 to rung 10 alone rises at
$+0.011815$, a factor $1.74$ over the global slope, and there was no
point between them to say whether that meant anything.

**Where this lands.** The fifteen-point crossing is $10^{7.3191}$,
inside the $[10^{7.1866},\,10^{7.6631}]$ published on eleven (Q5), so
no bracket is broken here. What moves is the instrument: any quantity
read off the eleven-rung slope past the end of the ladder —
$\theta'=0.56$ included — was read off a line now known to be too
shallow in the one region where it has been checked at finer than a
doubling. The fifteen-point refit does not repair that; it is still
one line through a region whose residuals have a sign.


#### Remark (fifty points to a doubling, and what they buy) {#rem:primorialdense}
<!-- evidence: audit_primorial_dense.py -->

Remark [rem:primorialgap] left a factor it could not read. The four
points it put between rung 9 and rung 10 rose at $3.83$ times the
ladder's global slope, but that was a two-point slope over a lever of
$0.3817$ in $\log N$, and propagating the ladder's own scatter through
the lever puts it at $1.90$ standard errors, which [rem:slopes]
forbids reading. The closing condition was named there and in the
repository's standing note: fill the space between the rungs until the
slope is separated by the **number of points** instead of by the length
of the lever. It is filled here. Every $N=30030m$ whose $m$ is composed
of $\{2,3,5,7,11,13\}$ with $m\in[2^6,2^{10}]$ — $202$ of them, fifty
to a doubling, the published rungs 6 through 10 among them — carries
the identical threshold $0.087306$ and the identical $18863$ admissible
$k$ (D2), and the imported measurement returns the five rungs the sweep
contains to $0.000048$ (D1).

**D3 is refuted, and it fixes what the rest can possibly say.** Fitted
inside each octave on that octave's own points, the exponent scatters
by $0.004084,\,0.004470,\,0.003975,\,0.003113$ — three of the four at
or above the ladder's published $0.0037$. That number was therefore
never an artefact of sampling one point per doubling. It is the
statistic's own noise: $K^*_R$ is where a step function first exceeds a
level, so it jumps, and the exponent at one $N$ of this family does not
predict it at the next. **The rungs' departures from their line are not
a shape waiting to be sampled**, and the hope that filling in would
resolve one is dead.

**D4 is refuted with it, and the closing condition is not met.** The
four octave slopes stand at $t=0.28,\,2.03,\,2.72,\,4.55$: three
resolve, the first does not. The failure is not one unlucky octave. A
per-point noise of the ladder's own size against a trend that gains
less than that across a whole doubling means an octave is too short a
lever whatever fills it — more points shrink the error of a slope, they
do not lengthen the base it is taken over. The factor $3.83$ is not
recoverable by this method and it is not recovered.

**What is recoverable is the slope over the whole sweep, and there D5
is refuted.** The four local slopes are
$-0.001039,\,+0.006957,\,+0.007174,\,+0.008091$, and regressed on the
$\log N$ of their own midpoints they rise at $+0.003983$, $2.17$
standard errors — a drift, not a constant. Taken as one fit the sweep
gives $+0.008688$ with standard error $0.000359$, $t=24.17$ over a
lever of $2.7726$ and correlation $0.86310$, against the eleven rungs'
$+0.006778$ with standard error $0.000565$: the two are $5.31$ of the
sweep's errors apart. And $64$ of the top octave's $78$ residuals about
the eleven-rung line are positive, so [rem:primorialgap]'s four of four
was not luck — though it is the octave's offset from the line that this
counts, not curvature within it. **The eleven-rung slope understates
the ladder in the region where the ladder has now been measured
densely, and anything read off it past the end — $\theta'$ included —
inherits that.**

**D6 holds, and it is worth more than the three that fell.** On the
$209$ distinct $N$ this ladder now has — the twelve published rungs and
the sweep — the best r.m.s. is $0.003978$ and carries a standard error
of $0.000196$, which is $4.9$ per cent of itself against the $22.4$ per
cent twelve rungs gave. **One shape survives at one standard error**,
$a+b\log N$, where twelve rungs left two. The margin is real and thin:
the gap to the runner-up is $0.000238$ against that $0.000196$, so it
separates, where on twelve rungs the same comparison was $0.000626$
against $0.000828$ and Remark [rem:laddershape12] had to call it tied.
The $2.7845$ decades that remark could not close are closed, not by an
argument about which shape is right, but by enough points to tell them
apart inside the range where both were fitted.

What the survivor says: $\tfrac12$ at $10^{7.3193}$ and
$\theta'=0.56$ at $10^{10.6180}$, bracketed at one r.m.s. residual
either way as [rem:primorialreach] brackets, $[10^{10.3993},\,
10^{10.8367}]$. Twelve rungs put it at $10^{11.0762}$, and the shape
that no longer survives put it at $10^{13.8607}$.

*Added later.* The survivor does not survive the rungs above this
range: see Remarks [rem:laddercurve], [rem:curvereach] and
[rem:rung15]. The $\tfrac12$ crossing is untouched and has since been
passed by direct measurement.

**And the bracket is a floor on the uncertainty, not the whole of it.**
Two things say so and neither is repaired here. The two best shapes
part by more than the fit's own r.m.s. at $10^{8.3256}$, only $0.5367$
decades above the top point $10^{7.7889}$ — Remark [rem:shapetrust]'s
boundary moves out by less than a fifth of a decade, so the forecast is
still far outside the range the data reach, and it is still the shape's
word and not the data's. What changed is that there is one shape to
take the word of instead of two. Second, D5: the slope this line
extrapolates is itself drifting, by $1.724191$ of its own mean across
the four octaves, and no shift of the level covers that.

The gate needed a repair to state this at all. Check G60 collected
`SHAPEGAP` by file order but `SHAPES TIED` as a union across files, so
the tie that twelve rungs correctly declared survived into a namespace
where $209$ points had separated the same two shapes, and contradicted
the current adjudication. It now reads, per target, only the file whose
`SHAPESURVIVE` point count is the largest — the one G53 already forces
to exist. An older file's marker stays as a statement about its own
point count and no longer binds the current one.


#### Remark (the noise is the field's, and both lenses see it) {#rem:densenoise}
<!-- evidence: audit_dense_noise.py -->

Remark [rem:primorialdense] showed that the ladder's scatter of
$0.0037$ is not a shape sampled once per doubling, and left one
question behind it with two opposite answers available. If the scatter
is the integrality of $K^*_R$ — which is a location, an integer, where
a step function first exceeds a level — then the ladder is being read
through a quantised instrument and the smooth equivalent
$\rho(N)<1\iff\log K^*_R/\log N>\tfrac12$ is the one to rebuild every
bracket on. If it is the field $|R(N;k)|$, nothing built on that field
escapes it. Both statistics are computed here on one walk of the same
$202$ points, so no two of them can disagree about the field: the
exponents return the sweep's own to $0.000050$, which is the width of
its printed column, and $\rho<1$ agrees with the exponent about which
side of the barrier every one of the $202$ falls on (N1), $\rho$ running
$0.9119$ to $1.2989$.

**The integrality is measurable and it is small.** The admissible $k$
are a fixed set of $18863$ values, so the crossing can only land on one
of them and the exponent moves in steps of
$(\log k'-\log k)/\log N$ — a quantity nobody had printed. Those steps
run $0.000019$ to $0.000618$ with median $0.000100$, and a step is a
bound on what quantisation can move, so its square bounds the share of
the variance it can account for. Per octave that share is
$0.0034,\,0.0010,\,0.0008,\,0.0003$: **under a hundredth of the
variance everywhere, and falling.** N2 is nonetheless refuted, and the
refutation is the rule's and not the reading's — the rule was written
against the largest step at any single $N$ and the largest one, in the
lowest octave, reaches $0.1595$ of the pooled scatter. On any variance
accounting the integrality is a rounding error and the rest is the
field.

**N3 is refuted, and it overturns a published comparison rather than
confirming it.** With the within-octave scatter as noise and the whole
sweep's slope as trend — two questions [rem:primorialdense] is the
reason for asking separately — the exponent gives noise $0.003876$
against $0.006022$ per doubling, a ratio of $0.6435$; $\log\rho$ gives
$0.039423$ against $0.063539$, a ratio of $0.6204$, on a slope of
$-0.091668$ at $24.91$ standard errors. Remark
[rem:primorialshare] measured the same pair on ten rungs at $2.829$
against $0.85$ and concluded $\rho$ was the worse instrument by a
factor of more than three. **The factor is gone.** Read against its own
error — an r.m.s. on $197$ degrees of freedom carries $1/\sqrt{2\,df}$,
a slope carries its own standard error — the two ratios differ by
$0.0231$ against a combined $0.0579$, which is $0.40$ of it. The
instruments are **tied**.

That was recoverable from the earlier remark's own mechanism paragraph,
which said the two are "the same noise seen through two lenses" and
computed the conversion. What its ratio comparison lacked was a scatter
that was only noise: $0.0689$ on ten rungs was noise and shape
together, and the shape half is what made $\rho$ look worse. **A ratio
of scatter to trend is a comparison between instruments only when the
scatter is the instrument's.**

**N4 is refuted too, and its refutation is the most useful thing here.**
Cutting the sum at $k=N^{1/4}$, the upper part holds $0.9257$ to
$0.9889$ of the mass, and the registered expectation was that it would
carry the fluctuation with it. It does not carry it in proportion. The
lower part's log scatters by $0.280725$ against the upper part's
$0.038733$ — **per unit of mass it fluctuates $7.25$ times as much** —
and weighting each by its mass share gives contributions $0.012212$ and
$0.037048$, so the $0.0435$ of the mass that sits below $N^{1/4}$
supplies $0.2479$ of what moves $\rho$. The fluctuation is not where
the mass is. A statistic that reweights toward the mass therefore does
not average it away, and neither does one anchored lower on the
$k$-axis, because three quarters of the contribution is still above the
cut.

**Where this lands.** There is no instrument left to switch to on this
ladder. The quantisation is a hundredth of the variance, the smooth
equivalent carries the same noise to within half its error, and the
noise is spread across the $k$-range rather than sitting on a thin set
that could be controlled separately. So the boundary
[rem:primorialdense] measured — the shapes parting at $10^{8.3256}$,
half a decade above the top point — is not moved by a better statistic
computed at the same $N$. It is moved by larger $N$, and that is a
budget.


#### Remark (part of the cancellation is between the blocks) {#rem:gainprofile}
<!-- evidence: audit_gain_profile.py -->

Remark [rem:leanidentity] reduced item 4(b) of the standing note to one
exponent: $\mathrm{slope}/\mathrm{floor}=(\ell^1/\ell^2)/G/c$ is exact,
so the slope stops growing against its floor exactly when $e(G)$
catches $e(\ell^1/\ell^2)$, and the measurement is $+0.153911$ against
a ceiling $\theta'/2=0.28$. Remark [rem:gainsplit] split the range by
magnitude and got the top tenth at $+0.077963$ against the remaining
nine tenths at $+0.340006$; Remark [rem:signmasshead] found the head
carries only half of what needs explaining. The note closed on the
sentence that the rise has to come from the whole range and not from
one decile.

That sentence rests on a two-way split, which cannot distinguish two
demands that cost a proof very different things. Order the $k$ by
$|a_k|$, cut into $B$ blocks of equal count, and write $w_d$ for block
$d$'s share of $\ell^1$ and $s_d=(\sum_d a)/(\sum_d|a|)\in[-1,1]$ for
its signed imbalance, so $G_d=1/|s_d|$. Then

$$
\frac1G=\Bigl|\sum_d w_d s_d\Bigr|
$$

identically, both sides being $|\sum a|/\ell^1$ — it holds here to
$1.110\cdot10^{-16}$ at every $N$, and the gains reproduce
[rem:gainsplit] to $0.000046$ with the same exponent to six places
(W1). The gain is therefore the reciprocal of a mass-weighted average
of ten signed imbalances, and it can be small either because each
$|s_d|$ is small — blocks cancelling internally — or because the $s_d$
have opposite signs, blocks cancelling against each other. Forbidding
the second gives $\sum_d w_d|s_d|$, and the two can be compared.

**W2 is refuted, and it names a mechanism nothing had named.** Of the
$0.153911$ at which $1/G$ decays, only $0.098386$ survives when the
blocks are forbidden to oppose one another; the difference
$0.055525$ stands at $3.67$ standard errors. **A third of the measured
cross-$k$ cancellation is between magnitude blocks, not inside them.**
That is a cheaper thing to ask a proof for than ten blocks each
cancelling at the square-root rate, and it was invisible to a two-way
split because a two-way split has almost no "between" to see.

**But the number is not one number, and saying so is the point.** With
$B$ blocks the forbidden-opposition value is $1/G$ at $B=1$ and exactly
$1$ at $B=\#k$, where it does not move with $N$ at all — so the share
attributed to opposition runs from $0$ to $1$ as the partition is
refined, and quoting it at one resolution reads as a canonical
decomposition when it is not. Measured at $B=2,5,10,20,50$ the share
is $0.2452,\,0.3684,\,0.3608,\,0.3216,\,0.2064$: **a factor $1.7848$
apart and not even monotone.** What is invariant is the direction — the
exponent falls at every resolution when opposition is forbidden — and
the coarsest reading, $B=2$, which is exactly the split
[rem:gainsplit] published, is the one that claims least at $0.2452$.
Gate check G72 now forces any such split to declare three resolutions
and to mark itself resolution-dependent when they disagree by more than
half again.

**W3 is refuted, and the profile is not a gradient.** The block gains
do not increase from top to bottom at any $N$, and four of the ten
blocks have exponents within their own standard error of zero. Only the
top two are both resolved and low — $+0.077963$ and $+0.112199$ at
$8.43$ and $8.98$ standard errors — and they hold $0.3400$ and $0.1867$
of the mass. Below them the blocks carry little mass and their gains
swing by two orders of magnitude between neighbouring $N$ (block 9
reaches $663.9864$ at $N=4\cdot10^5$ and $6.0390$ at
$6.4\cdot10^6$), so their individual exponents are not measurements of
anything. **The object that fails to cancel is the top fifth by
magnitude, not the top tenth and not a gradient**; every block is under
the ceiling and the mass-weighted shortfall is $0.227835$ in the
exponent.

**W4 holds, and only by the letter of its rule.** The top block
supplies the largest single term $w_d|s_d|$ at every one of the eight
$N$ and its share does not fall, which is what the rule asked. But the
share runs $0.5466$ to $0.6033$ with an exponent of $+0.007716$ at
$0.96$ standard errors, so the rise is unresolved and the honest
reading is that it is **flat**: over a factor $128$ in $N$ the top
block's grip on the whole neither tightens nor loosens. The script
declares it unresolved rather than reading the sign.

**Where this lands.** The demand on item 4(b) is now stated in the
right variables. It is not "ten blocks each cancelling to the
square-root rate", because a quarter to a third of the decay is already
opposition between blocks and that component is real at every
resolution. It is not "control the head", because the head's share of
the weighted sum is flat and the second block is resolved and low as
well. It is: **raise the internal cancellation of the top fifth by
magnitude, or find more opposition between blocks than the field
already supplies** — and the second alternative is new, has a measured
size, and has nothing said about it anywhere. (The second alternative does not
survive. Remark [rem:gainopposition] runs the two nulls it needs and
finds $\mu$'s opposition *below* chance at every $N$ with the chance
level itself flat — by an identity rather than a fit — so that route is
a bounded factor and not an exponent, and $0.2318$ of the $0.360763$
above is what a coin gets from the partition alone. What is left is the
first alternative.)


#### Remark (the opposition is below chance and chance does not grow) {#rem:gainopposition}
<!-- evidence: lab_gain_opposition.py -->

Remark [rem:gainprofile] measured that a third of the decay of $1/G$
sits between the magnitude blocks rather than inside them and closed by
naming a second route to $\theta'/2$: find more opposition between
blocks than the field supplies. **That route does not exist in the form
it was named, and this remark is the correction.** The sentence rested
on a difference of exponents, which is a rate; whether opposition can
be had is a question about a level, and it has an answer that is exact
rather than fitted.

Two one-sided nulls decide it. The coin arm randomises every sign on
$\mu$'s own magnitudes, so $\#k$, the block edges and every $w_d$ are
preserved and only $\mu$'s signs are broken — the arm
[rem:crosskreference] used. The block-sign arm randomises only the ten
block signs, keeping every $|w_d s_d|$ exactly as measured, so it breaks
the blocks' opposition and keeps their internal cancellation; it is
enumerated over all $1024$ patterns and needs no draws. The controls
reproduce [rem:gainprofile] to $3.983\cdot10^{-7}$ and the identity
$1/G=|\sum_d w_d s_d|$ to $1.110\cdot10^{-16}$ (O1).

**O2 is refuted: part of the split is the partition's.** The coin arm
has a between-block component of its own, $0.083641$ at $10$ blocks,
separated at $4.95$ standard errors — ordering by magnitude and cutting
manufactures some opposition all by itself, because the blocks are then
sorted by size and a sorted partition of a random-sign sequence does not
have independent block sums of equal scale. So $\mu$'s $0.360763$ is not
all $\mu$'s: $0.2318$ of it is available to a coin, and what is left,
a factor $4.31$ more, is the field's. The direction of
[rem:gainprofile] survives; its size is overstated by about a quarter.

**O3 holds, and it reverses the reading.** Write
$F=(\sum_d w_d|s_d|)/|\sum_d w_d s_d|\ge1$ for the factor by which
forbidding opposition raises $1/G$. Measured, $F$ runs $1.0316$ to
$1.3829$. Under the block-sign arm — same ten magnitudes, random
signs — it runs $1.6575$ to $1.8296$, and the ratio $F/F_{\text{null}}$
runs $0.6224$ to $0.8343$: **$\mu$'s blocks are more aligned than chance
would leave them, at every one of the eight $N$.** The ratio rises at
$+0.063241$, $8.42$ standard errors. The growth
[rem:gainprofile] attributed to opposition is $\mu$ closing a gap to
chance, not opening one past it.

**O4 holds, and it is an identity rather than a measurement.** If the
largest $|w_d s_d|$ exceeds the sum of the other nine then every signing
carries that term's sign, the enumerated mean of
$|\sum_d\pm|w_d s_d||$ is that term itself, and

$$
F_{\text{null}}=\frac{1}{\text{share of the dominant block}} .
$$

The largest term does exceed the rest at $8$ of $8$ $N$ — the shares are
$0.5659$ to $0.6033$, all above a half — and the identity holds to
$1.776\cdot10^{-15}$. So $e(F_{\text{null}})=-e(\text{top share})$
exactly, and [rem:gainprofile] already measured that share as flat at
$0.96$ standard errors. Numerically $F_{\text{null}}$'s own exponent is
$-0.007716$ at the same $0.96$, as it must be.

**Therefore the route is a bounded factor and not an exponent.**
Opposition among blocks of these relative sizes is worth a factor
$1/\!$(top share), that factor does not grow, and $\mu$ sits below it
with a factor $1.1987$ of it still unspent at the top $N$. The
$+0.055525$ that [rem:gainprofile] credited to opposition is what
closing that gap costs; when the gap closes it stops, and going further
means beating chance among ten numbers, for which nothing in this
programme names a mechanism. **Item 4(b) is harder than
[rem:gainprofile] left it, not cheaper**, and the honest exponent
available from within-block cancellation alone is the $0.098386$ that
remark measured — against a ceiling of $\theta'/2=0.28$.

**O5 holds and sharpens where the opposition lives.** The $\ell^1$ share
carried by blocks whose sign differs from the dominant block's runs
$0.0088$ to $0.0670$ with an exponent of $+0.216756$ at $1.73$ standard
errors, so it is not resolved as growing and the script declares it so.
Under $7$ per cent of the mass is doing the opposing at every $N$. A
mechanism that made the opposing mass grow would be a route even with
$F$ below chance, since $F$ is a ratio and mass is not; nothing here
shows one.

Both between-block shares are resolution-dependent, as gate check G72
requires them to declare: $\mu$'s runs $0.2452,\,0.3684,\,0.3608,\,
0.3216,\,0.2064$ over $B=2,5,10,20,50$ and the coin's
$0.0224,\,0.0664,\,0.0836,\,0.1003,\,0.1155$. The coin's rises with $B$
where $\mu$'s peaks and falls, which is the shape the argument above
predicts: refining the partition manufactures more of the artefact and
leaves less room for the real thing.


#### Remark (the inner sums already cancel; the sign is what does not) {#rem:headsign}
<!-- evidence: audit_head_sign.py -->

Remark [rem:gainopposition] left item 4(b) with one demand: raise the
internal cancellation of the top blocks by magnitude, since the
opposition route is a bounded factor and the honest exponent from
within-block cancellation is $0.098386$ against a ceiling
$\theta'/2=0.28$. Those blocks fail to cancel because
$0.8274$ to $1.0000$ of the top decile carries one sign
([rem:gainsplit]), and **why** had never been asked.

Split the inner sum of [eq:layers] where its sign is. With
$P_k=\sum_{\mu(m)=+1}\Lambda(N-mk)$, $M_k$ the same over
$\mu(m)=-1$, $H_k=P_k-M_k$, $T_k=P_k+M_k$ and $I_k=H_k/T_k$, the
magnitude factors as $|a_k|=(\log k)\,T_k\,|I_k|$ — a mass and an
imbalance. The split reproduces [rem:gainsplit]'s head fractions and
gains to $0.000046$ (Q1).

**Q2 is refuted: the head is both.** The top decile's average $|I_k|$
exceeds the whole range's by $2.0992$ to $2.5633$, as predicted, but so
does its average $T_k$, by $1.4707$ to $4.6521$ — and that mass
enrichment *grows* along the sweep where the imbalance enrichment does
not. The registered dichotomy, mass selection washes out and imbalance
selection does not, was the wrong frame.

**Q3 holds and is an identity, which is a warning about the check and
not a strength of it.** $T_k>0$, so
$\operatorname{sign}I_k=\operatorname{sign}H_k=\operatorname{sign}a_k$
and the share of the head with $I_k<0$ *is* its one-sign fraction
whenever the majority is negative; the difference is $0.0000$ at all
eight $N$ because the two columns are one column. What Q3 measures is
only that the majority is negative, which is [rem:headidentity]'s
finding.

**Q4 holds, and the reference it needed changes what it means.** The
average $|I_k|$ falls at $-0.191871$, $40.29$ standard errors. But
independent signs would not give zero: the inner sum has $n_k$
contributing terms and random cancellation leaves an imbalance of order
$1/\sqrt{n_k}$, a quantity computable term by term that had never been
put beside $|I_k|$. The average $n_k$ runs $100.33$ to $843.03$ and its
own reference falls at $-0.184058$; the ratio of the measured imbalance
to it runs $1.0310$ to $1.0995$ with a drift of $-0.007812$ at $1.78$
standard errors, **unresolved**. So:

> **the inner sums already cancel at the independent-sign rate, to
> within three to ten per cent and no resolved drift. There is nothing
> to win inside $k$, and the whole deficit against $\theta'/2$ is
> across $k$.**

**And the sign is not random even though the size is.** Under
independent signs every decile would sit at a half. Cut on $|a_k|$ the
share with $I_k<0$ runs $0.8274$ down to $0.4989$; the top-minus-bottom
spread is $+0.3284$. Cut on each factor alone it is $+0.3389$ on
$|I_k|$, $+0.1537$ on $T_k$ and $+0.0126$ on $k$. **The correlation
lives on the imbalance axis and essentially not on $k$** — so the
$|a_k|$ profile is inheriting $|I_k|$'s, and the head's alignment is
not the $k$-size effect the discussion around [eq:layers] proposed when
it observed that small $N/k$ leaves the $m$ almost all prime.

That discussion is not wrong about the bias, only about where to look
for it. The signed mean of $I_k$ is negative at every $N$,
$-0.106299$ to $-0.017828$, and by $k$-decile at the top $N$ it runs
$-0.0023,\,-0.0074,\,-0.0082,\,-0.0129,\,-0.0185,\,-0.0191,\,-0.0291,\,
-0.0273,\,-0.0240,\,-0.0295$ — growing in size with $k$, a spread of
$0.027147$ against a mean of $-0.017828$, which is the direction
[eq:layers] predicted. **The bias is $k$-dependent; the head is not
$k$-selected; those are consistent, and [rem:headidentity] is why** —
the head's median sits at the $0.22$ to $0.32$ point of $k$-order.

**The bias is also disappearing, measured against the only scale it can
be measured against.** In units of $1/\sqrt{n_k}$ it runs $0.6339$ down
to $0.2590$, falling at $-0.185280$ with $20.72$ standard errors, and
the whole range's share with $I_k<0$ runs $0.6169$ down to $0.5128$.
A bias that is a quarter of the noise and shrinking is not, on its own,
an obstruction.

**Which leaves one thing unexplained, and it is now the whole of item
4(b).** The top decile by $|I_k|$ is $0.8547$ negative at the top $N$
and falls only at $-0.032100$ ($6.94$ standard errors) — far slower
than the bias that is supposed to produce it. A bias of $0.2590$ of the
fluctuation scale does not obviously put $85$ per cent of the largest
deviations on one side; reconciling the two needs a null that signs the
observed $|I_k|$ at the observed rate and reads off the tail, and this
script did not run one. **So the demand on item 4(b) is no longer
"cancel the top fifth" — the inner sums are already at the
independent-sign rate — it is: account for a tail asymmetry that the
measured bias appears too small to explain.**


#### Remark (a null refuted before it could adjudicate) {#rem:tailasymmetry}
<!-- evidence: lab_tail_asymmetry.py -->

Remark [rem:headsign] left one thing between item 4(b) and a statement:
the inner sums of [eq:layers] already cancel at the independent-sign
rate, so the whole deficit against $\theta'/2$ is a sign correlation
across $k$, and that correlation shows as a top decile of the $k$ by
$|I_k|$ which is $0.8547$ negative at the top $N$ while the bias
supposed to produce it is only $0.2590$ of the fluctuation scale. Three
arms were registered to settle whether the bias accounts for the tail.
**T1 gates and T1 is refuted, so it does not settle it**, and the shape
of that failure is the result.

**Two things failed in T1 and they are not the same kind of thing.**
The calibration cap was $10^{-6}$ on the arm's average $|I_k|$, which is
a Monte Carlo mean over $400$ draws and therefore carries a sampling
error of its own; it came in at $1.938\cdot10^{-3}$. **That rule was
badly written and it is reported as written, not relaxed.** The second
failure is not a cap. The arm's bias profile is the measured one, so its
mean $I_k$ is the observed mean by construction — and its sign share is
not the observed sign share, standing above it by $+0.0615$ to
$+0.1076$ at all eight $N$. A shifted symmetric law with the right mean
gets the sign share wrong, in the direction of too many negatives; so
the observed distribution carries the same negative mean with **fewer
small negatives and larger ones**, which is to say the negativity is
concentrated in the tail rather than spread as a shift. The script exits
non-zero on this, as its own rule requires.

That blocks T3. Its numbers are computable — the observed tail stands
$+1.72$ to $+7.60$ of the structured arm's spreads above it, and the
arm reaches only $0.6976$ against the observed $0.8547$ at the top $N$ —
but they are a comparison against a model already known to be the wrong
shape, and the registered rule says they may not be read as evidence
about the tail. T2 and T4 are unaffected by the shape question and
stand: signing the observed magnitudes at one overall rate falls short
by $+4.66$ to $+14.91$ of its own binomial spreads (T2), and letting the
bias vary over $k$ rather than holding it flat is worth $+0.97$ to
$+2.06$ spreads (T4).

**What does settle the question is an arm with no distributional
assumption, and it was not registered.** Keep every $|I_k|$ exactly as
measured and draw only the signs, each $k$ at the negative rate its own
decile in $k$ exhibits. This assumes nothing about the law and keeps
whatever $k$-dependence the sign rate has; the top decile by $|I_k|$ is
then a fixed set, so the arm is closed form. It gives $0.5162$ to
$0.6395$ against an observed $0.8547$ to $1.0000$ — the observed tail
stands **$+4.28$ to $+14.77$ of that arm's spreads above it**, and the
sign rate by $k$-decile at the top $N$ is flat, $0.4853$ to $0.5420$.

So three ways of accounting for the tail have been tried and all three
miss it: a uniform sign rate on the observed magnitudes, a
$k$-dependent sign rate on the observed magnitudes, and a shifted
symmetric law with the measured bias — the last of which is also the
wrong shape independently. **The coupling between the size of $I_k$ and
its sign is direct: it is not mediated by $k$, it is not a shift, and it
is what item 4(b) is now about.** It has no name in this programme and
no measurement beyond the deciles of [rem:headsign].

One number is worth keeping for whoever measures it. The observed tail
falls at $-0.032100$ ($6.94$ standard errors) while the structured arm's
falls at $-0.054160$ ($18.28$): the thing that is unexplained is not
merely present, it is **decaying more slowly than the bias-driven part
of it**, so the gap between them widens along the sweep rather than
closing.


#### Remark (the one candidate, eliminated for being constant) {#rem:tailmertens}
<!-- evidence: audit_tail_mertens.py -->

Remark [rem:tailasymmetry] said the coupling between the size of $I_k$
and its sign has no name in this programme. It said so without testing
the one predictor of $\operatorname{sign}H(N;k)$ the repository has.
Remark [rem:leanmertens] found agreement with
$\operatorname{sign}M_{\mathrm{odd}}(\lfloor N/k\rfloor)$ at $0.7669$ to
$0.7704$ on inner lengths $2\le N/k\le1000$, and
[rem:oddmertensrange] found it does not transfer to $k<N^{\theta'}$,
where the agreement is $0.5201$ to $0.6161$. But the tail was never
looked at, and it is a different question: a predictor right half the
time on average can still be right almost always where $|I_k|$ is
largest. There was a reason to expect it might —
$M_{\mathrm{odd}}(\lfloor N/k\rfloor)$ is violently non-monotone in $k$,
so the $k$-decile arm of [rem:tailasymmetry] would have averaged
exactly this structure away.

**It is eliminated, and the way it fails is the finding.** The controls
hold to $0.000048$ on the published agreement and exactly on the
observed tail share (R1). R2 holds too — the tail agreement exceeds the
whole range's by $+0.3383$ to $+0.4453$ — but **R2 is empty and R3 and
R4 are what say so.**

The predictor is *constant* on this range. It is negative on $0.9829$
to $0.9970$ of the $k$, and on **every** $k$ of the tail at every one of
the eight $N$; of the $3936$ distinct values of $\lfloor N/k\rfloor$ at
the top $N$, $3843$ have $M_{\mathrm{odd}}<0$. A constant predictor
agrees with anything at the rate that thing takes its value and carries
no information, so:

* the tail agreement is numerically **identical** to the tail's negative
  share, $1.0000,\,1.0000,\,0.9853,\,0.9600,\,0.9662,\,0.9315,\,0.8854,\,
  0.8547$, and the agreement profile by decile of $|I_k|$ —
  $0.8547,\,0.5903,\,0.5021,\,0.4370,\ldots$ — is [rem:headsign]'s
  negative-share profile read back;
* the permutation baseline reaches it (R4 refuted), gaining $+0.0000$ at
  six of the eight $N$ and at most $+0.0066$, because permuting a
  constant changes nothing;
* and R3 is refuted in the only direction a degenerate predictor can
  fail: it puts every tail sign negative, so its top-decile share is
  $1.0000$ exactly, against an observed $0.8547$, and the gap **widens**
  from $+0.0000$ to $+0.1453$ along the sweep because the observed share
  falls while the prediction cannot.

**This strengthens [rem:oddmertensrange] rather than contradicting
it.** Its check P2 recorded that the agreement is below $0.70$ on
$k<N^{\theta'}$ and held. The truer statement is stronger: on that range
the predictor is not weak, it is **informationless**, and the $0.5201$
it reports is $H$'s own marginal negative rate reflected back. An
agreement fraction is only a measurement when the predictor has
variance, and nothing in the gate was checking that.

**Where this leaves item 4(b).** The coupling is still unnamed, but the
statement has changed from "no candidate has been tried" to "the one
candidate has been tried and is degenerate on the range that matters".
The two exponents are also worth keeping: the tail's negative share
falls at $-0.032100$ ($6.94$ standard errors) and the whole range's
agreement at $-0.035401$ ($4.88$), so both drift at the same slow rate
and neither is on a path to a half within reach of computation.

**And the gate now checks this everywhere, which is what says how far
the fault reaches.** Check G73 requires any result file declaring a sign
agreement to publish the predictor's own majority sign share — the
largest over every set it reports on, so that one degenerate window
cannot hide behind a well-conditioned one. The convention started as the
smallest, and [rem:oddmertensrange]'s two windows, $0.7928$ against
$0.9970$, are what forced the change. Nine files declare such an
agreement and all nine now publish the share. **Only two are
degenerate**: this remark's, at $0.9970$, and [rem:oddmertensrange]'s, at
the same $0.9970$ from its $k<N^{	heta'}$ window. Every other predictor
has real variance — $0.8316$ and $0.8081$ for the Mertens predictor on
$2\le N/k\le1000$, so [rem:leanmertens]'s agreement stands as a
measurement; and $0.7726$, $0.7091$, $0.6686$, $0.5879$, $0.5801$ for
the sieve predictors of the survivor, level-threshold, survivor-range,
log-weight and [rem:sievedepth] results. So the degeneracy is not a
fault of agreement measurements in this programme. **It is specific to
the Mertens predictor pushed onto the range where the gain lives**, and
every other agreement claim in the paper survives the check that caught
it.


#### Remark (a bounded modulus does reach the tail) {#rem:tailpredictors}
<!-- evidence: audit_tail_predictors.py -->

Two things were outstanding after [rem:tailmertens]. The Mertens
predictor was eliminated for being constant on the range that matters,
and G73's survey then said the fault was narrow: the sieve predictors
declare majority shares of $0.5801$ to $0.7726$, so they have real
variance, and none of them had been looked at on the tail. And a
boundary had never been checked. **Every sign predictor in this
repository sums over odd $m$; the gain of item 4(b) sums over all $m$.**
They ought to coincide, since $N$ is even and $k$ coprime to $N$ is odd,
so an even $m$ makes $N-mk$ even and $\Lambda$ vanishes unless $N-mk$ is
a power of two — but that is an argument, not a measurement.

**S1 settles the boundary exactly.** $\operatorname{sign}H$ over all $m$
and over the odd $m$ agree at $1.000000$ of the $k$, at every one of the
eight $N$ — not approximately, at every $k$. And the $Q=29$ agreement
reproduces [rem:sievedepth]'s published column to $0.000050$. So the
predictors were always predictors for the gain's field, and everything
compared across the two conventions in the preceding remarks was
compared inside one field.

**S2 is refuted, and it changes how any of this may be read.** The
sieve predictor at $Q=29$ has a majority sign share of $0.5170$ to
$0.6255$ on the whole range — genuine variance — but on the top decile
by $|I_k|$ that rises to $0.8484,\,0.8820,\,0.9309,\,0.9595,\,0.9600,\,
0.9851,\,1.0000,\,1.0000$. **Selecting the tail pushes the predictor
towards one sign by itself.** So the level of agreement on the tail is
not a measurement for this predictor either, and G73 marks this file
degenerate at $1.0000$ accordingly. What is left as a measurement is the
*excess* over what the predictor's own marginals give, and unlike the
Mertens case there is one, because the share stops short of $1$ at the
larger $N$ and falls as $N$ grows.

**S3 holds.** Replacing every sign by $\operatorname{sign}P_{29}$ on
$\mu$'s own magnitudes reproduces the observed top-decile negative share
at every $N$, worst difference $-0.0147$: $0.9706$ against $0.9853$,
$0.9224$ against $0.9315$, $0.8484$ against $0.8547$. The tail's sign
pattern is reproduced by a **bounded** modulus.

**S4 is refuted as registered, and where it fails is not where it
matters.** The matched-marginal arm draws signs independently at the
predictor's own negative rate on the tail, so it contains everything the
marginals can explain. The tail agreement beats it by $+3.06,\,+5.36,\,
+8.24,\,+11.94$ standard spreads at the top four $N$ — resolved, and
**growing with $N$** — but by only $+1.42$ and $+2.88$ at
$8\cdot10^5$ and $1.6\cdot10^6$, and the two smallest $N$ have no
spread at all because the arm is already at $1.0000$ there. The rule
asked for three spreads at every $N$ and does not get them. The honest
reading is that the excess is unmeasurable at the small $N$, where the
predictor has no variance on the tail to speak of, and strongly resolved
at the large ones, where it does.

**Where this leaves item 4(b), which is further than it has been.** The
tail's sign — the thing [rem:tailasymmetry] could not account for with a
uniform rate, a $k$-dependent rate or a shifted law, and
[rem:tailmertens] could not account for with the Mertens function — is
accounted for by the survivors of a sieve to $Q=29$. That is a bounded
modulus, which is the property [rem:provablehalf] needs and
[rem:sievedepth] found nothing had for the *slope*. The agreement does
not decay away either: on the tail its exponent is $-0.003755$ at
$4.06$ standard errors, essentially flat across a factor $128$ in $N$,
while the whole-range agreement's $-0.009494$ is unresolved at $1.49$.

What this does **not** say. It does not say $P_{29}$ carries the gain:
[rem:sievedepth] measured the slope ratio at a fixed level and found it
short, and [rem:leveldemand] measured that the residual demand barely
moves from level $29$ to $\alpha=0.3$. It says the *sign* of the largest
imbalances is a bounded-modulus object, which is one of the two factors
$|a_k|=(\log k)T_k|I_k|$ splits into, and the one item 4(b) was stuck
on. Whether a bounded level can also supply the *magnitudes* is a
different question and this remark does not touch it.


#### Remark (the level's signs, and why swapping them proves nothing) {#rem:levelmagnitude}
<!-- evidence: audit_level_magnitude.py -->

Remark [rem:tailpredictors] left one factor of
$|a_k|=(\log k)T_k|I_k|$ untested. [rem:leanidentity] gives the question
an exact form — the demand is $e(G)\to e(\ell^1/\ell^2)$, measured at
$+0.153911$ against $+0.287798$ — and both sides are computable for any
vector, so three surrogates separate the factors on the gain's own
field: $|a_k|$ with $\operatorname{sign}P_{29}$, $(\log k)|P_{29}|$ with
$\operatorname{sign}H$, and $(\log k)P_{29}$. No gain had been computed
for any of them; [rem:sievedepth] and [rem:leveldemand] measure
one-sided sums, and a gain is a ratio of two norms.

The control is exact: the gain and its exponent reproduce
[rem:gainsplit]'s to $0.000046$ and to six places (V1).

**V3 and V4 hold and are the durable part.** The level's *magnitudes*
are ruinous: on $\mu$'s own signs they give $e(G)=+0.066210$, below
$\mu$'s by $6.95$ standard errors, and a deficit of $+0.219441$. Taken
whole the level is worse than $\mu$ — $e(G)=+0.103219$, deficit
$+0.182432$ against $\mu$'s $+0.133887$. **A bounded modulus, used for
both factors, cancels no better than $\mu$ relative to its own
concentration.**

**V2 is refuted, and the refutation looked like the best news in this
programme until it was checked.** Putting the level's signs on $\mu$'s
magnitudes raises the exponent to $+0.260989$, $5.05$ standard errors
above $\mu$'s, and cuts the deficit to $+0.026809$ — a fifth of
$\mu$'s, and within a hair of the ceiling. Read naively that says a
bounded modulus almost supplies the cancellation item 4(b) needs.

**It does not, and the reason is that a small deficit is what noise
looks like.** A random sign vector has $|\sum|$ of order $\ell^2$, so
its gain *is* $\ell^1/\ell^2$ up to a constant and its deficit is zero
by construction. Measured on $\mu$'s own magnitudes the coin gives
$e(G)=+0.289456$ and a deficit of $-0.001658$ — zero, as the argument
says — with gains of $17.3688$ to $72.1485$ against $\mu$'s $1.8337$ to
$3.5925$, the factor [rem:crosskreference] already published. So
"deficit near zero" is not a property to aspire to; it is what a
sign vector carrying no information about $\mu$ achieves.

The arm that decides it is $\mu$'s own signs with **exactly as many
flipped at random as the predictor gets wrong** — the agreement is
$0.7367$ to $0.8129$, so $59$ to $1197$ flips. That arm gives
$e(G)=+0.180978$ and a deficit of $+0.106820$. The sign swap's
$+0.260989$ sits above it by $+0.080011$, $2.79$ standard errors, so the
level's signs are not merely $\mu$'s signs with noise added. **But the
separation is not where it would have to be.** At the top of the sweep
the two gains are $7.0377$ and $7.2326$ — the arm is *ahead* — and the
exponent gap is made at the small-$N$ end, where the sign swap gives
$2.0355$ against the arm's $2.9068$. The level's signs cancel worse than
random flips of the same count at small $N$ and no better at large $N$;
what the exponent records is the crossing, not an advantage.

**What is left of the question, stated in the right variables.** The
deficit scale now has both ends measured: $0$ for a coin, $+0.133887$
for $\mu$. Item 4(b) asks for $\mu$'s own sum to behave like the coin's
on $\mu$'s own magnitudes, and nothing measured here does that *for
$\mu$*. A sign swap cannot, in principle, say anything about it —
replacing $\mu$'s signs discards the object a proof has to bound, and
the surrogate's smallness is not $\mu$'s smallness. So the demand is
unchanged and now has a scale: close $+0.133887$, on the vector
$\mu$ actually gives, with the concentration $\ell^1/\ell^2$ fixed at
$+0.287798$ — which is itself only $+0.007798$ above the ceiling
$\theta'/2$ that [rem:flatnessshape] shows it cannot pass.

*Added later.* The scale is now $+0.134019$ onto a measured
$+0.283586$ on the field to $1.024\cdot10^8$ ([rem:fieldreach]), the
concentration's excess over the ceiling is $6.20$ standard errors
rather than $+0.007798$'s worth, and the deficit has since been
localised in the head and identified with a truncated Chebyshev
correlation — [rem:sumhead], [rem:denominator]. The reading above is
unchanged; every number in it has moved.


#### Remark (the excess over the ceiling was eight points) {#rem:flatnessfill}
<!-- evidence: audit_flatness_fill.py -->

Remark [rem:levelmagnitude] left item 4(b) with one computational axis:
the concentration exponent $e(\ell^1/\ell^2)=+0.287798$ sits $3.15$
standard errors **above** the ceiling $\theta'/2=0.28$ that
$\#k\asymp N^{\theta'}$ imposes, and [rem:flatnessshape] explained that
by $F=(\ell^1/\ell^2)/\sqrt{\#k}$ running $0.6622$ to $0.6986$ and
still rising towards the Cauchy–Schwarz bound $F\le1$. Every $N=2^a5^b$
has odd radical $5$, so the family can be filled without leaving the
field: there are $70$ such $N$ in $[2\cdot10^5,\,2.56\cdot10^7]$
against the $8$ doublings that had been used. The controls are exact —
$\#k$, $\ell^1/\ell^2$, $F$ and the gain reproduce the published to
$0.000050$ and $0.000046$ (U1).

**U3 is refuted, and that dissolves the puzzle instead of solving it.**
On $70$ points $e(\ell^1/\ell^2)=+0.280072$ with a standard error of
$0.004748$. The ceiling is $0.2800$. **The excess is $+0.000072$,
which is $0.02$ standard errors.** The concentration exponent sits on
the ceiling, where the arithmetic says it must; the $3.15$ standard
errors were eight points.

**U4 is refuted with it: $F$ does not rise.** Over the whole sweep its
slope is $+0.001803$ at $1.23$ standard errors and over the top octave
$-0.006772$ at $0.53$ — neither resolved, and the script declares both
so. $F$ runs $0.6622$ to $0.7388$ over the $70$. So the rise
[rem:flatnessshape] reported was the eight points wobbling, and **the
extrapolation built on it is void**: that remark's $e=+0.007445$ giving
$F=1$ at $10^{28.6782}$, its bracket
$[10^{24.2459},\,10^{28.6782}]$ and the parting point
[rem:shapetrust] reads off it all rest on a slope that is not resolved
on the denser set.

**U2 is refuted too, and it explains why more points do not settle the
shape.** Fitted inside each octave on its own points, $\log F$ scatters
by $0.022751,\,0.028360,\,0.011050,\,0.015638,\,0.012620,\,0.010021,\,
0.008713$ — every one above the published eight-point $0.006479$, which
was two parameters fitted to eight points and therefore an
underestimate. Refitting the two shapes on all $70$ gives $0.011916$
for the power law against $0.011929$ for $a+b/\log N$, a gap of
$0.000014$ against the r.m.s.'s own standard error $0.001022$. They
are **more** tied than before. More points in the same range do not
decide a shape; only a longer lever would.

**And the eight-point standard errors were too small on both sides.**
Refitted on the $70$, $e(G)=+0.142121$ with a standard error of
$0.033937$ against the published $+0.153911\pm0.011253$ — the value
holds, the error triples. $e(\ell^1/\ell^2)$'s error nearly doubles,
$0.002472\to0.004748$. The gain is the noisy one: its octave-wise
$\log$ scatter runs $0.278278$ to $0.491612$, so $G$ moves by half
again between neighbouring $N$ of the same family, where
$\ell^1/\ell^2$ moves by a twentieth of that. That is what a ratio with
a cancelling sum in the denominator does, and it is why the doublings
looked steadier than the family is.

**Where this leaves item 4(b): the target is now an arithmetic
constant.** The deficit on the dense set is $+0.137951$, against
$+0.133887$ on eight — unchanged within the errors. What changed is
what it has to close *to*. It is no longer a fitted $+0.287798$ sitting
unexplained above a bound; it is $\theta'/2=0.2800$ itself, which
follows from $\#k\asymp N^{\theta'}$ and not from any fit. So the
demand reads: **$e(G)$ must travel $+0.137879$, from $+0.142121$ to
$\theta'/2$, on the vector $\mu$ actually gives** — with the other side
of the identity pinned by arithmetic rather than by a curve whose shape
this repository has now twice failed to determine.

*Added later.* Every reading above that rests on the $70$ points is
withdrawn. The family is not one field: see Remark [rem:fillfield].
What survives untouched is U1, the controls, and the observation that
$G$ is the noisy side of the identity.


#### Remark (the seventy points were three fields) {#rem:fillfield}
<!-- evidence: audit_fill_field.py -->

Remark [rem:flatnessfill] filled the family with every $N=2^a5^b$ in
$[2\cdot10^5,\,2.56\cdot10^7]$ on the ground that "every $N=2^a5^b$ has
odd radical $5$, so the admissible $k$-set and the threshold are fixed
exactly as along the doublings". The sentence is false at both edges of
the family, because the enumeration starts at $a=0$ and at $b=0$:

* $N=2^a$ has odd radical $1$, and there $k$ ranges over the squarefree
  $k$ coprime to $2$ alone — every $k$ divisible by $5$, which the
  doublings exclude, is admitted. Seven of the $70$.
* $N=5^b$ is **odd**, and there $k$ ranges over the squarefree $k$
  coprime to $5$ alone — every **even** $k$ is admitted. Three of
  the $70$.

Three coprimality classes were fitted as one. Nothing rang because the
`RADICALS` line counts *odd* radicals and $2^a$ contributes the empty
one, so the file declared $2$ and G34 was satisfied; the partition by
what actually fixes the $k$-set gives $3$. The control is exact:
recomputing $\#k$ from the squarefree $k$ coprime to $N$ reproduces all
$70$ printed counts with a worst departure of $0$, and the ten
off-field $N$ are exactly the ten whose $k$-set is not the
coprime-to-$10$ one (V1).

**The ten are a different field, not a loud sample of this one.** Every
off-field gain lies outside the range of on-field gains within a factor
of two of it, both ways (V2). On-field $G$ runs $1.7228$ to $3.7589$
across the whole family; off-field, $3.7926$ to $12.6968$.

**On the field alone the excess over the ceiling returns, sharper than
on eight points.** With the $60$ on-field $N$, $e(\ell^1/\ell^2)
=+0.284380$ with a standard error of $0.000880$, so the excess over
$\theta'/2$ is $+0.004380$ — $4.98$ standard errors, against $3.15$ on
the doublings and the $0.02$ the mixture reported (V3). What the
mixture did was not add noise symmetrically: the three classes carry
almost the same $\#k$ exponent ($+0.560714$, $+0.560414$, $+0.560174$)
and quite different gains ($+0.298478$, $+0.146489$, $+0.081546$), and
the off-field points sit above the on-field concentration at their own
$N$ because their $k$-set is denser. The mixture's $r.m.s.$ of $\log G$
is $0.397647$ against $0.042144$ on the field.

**$F$ is resolved rising on the field.** The measured ceiling
$e(\#k)/2=+0.280207$, and $e(F)=+0.004166$ with a standard error of
$0.000874$, $t=4.76$ (V4). On the mixture the same slope was
$+0.001798$ at $t=1.23$ — which is exactly the unresolved sign
[rem:flatnessfill] read as "$F$ does not rise". So the rise
[rem:flatnessshape] reported is a fact about the field and not eight
points wobbling.

**But the extrapolation stays void, for the other reason.** Refitting
the two shapes on the $60$ gives $0.006480$ for the power law against
$0.006520$ for $a+b/\log N$, a gap of $0.000041$ against the r.m.s.'s
own standard error $0.000602$ (V5). A clean field does not decide the
shape either. And $F$'s octave-wise scatter is not below the published
eight-point $0.006479$ everywhere — $0.006673$, $0.009468$, $0.006251$,
$0.012181$, $0.004813$, $0.006741$, $0.008411$, four of seven above it
(V6). So $e=+0.007445$, the arrival $10^{28.6782}$, its bracket and the
parting point [rem:shapetrust] reads off it remain withdrawn; the
slope that replaces $+0.007445$ is $+0.004166$, and no arrival is
published from it here, because V5 says the shape carrying it is
undetermined.

**Where this leaves item 4(b): the target is not an arithmetic
constant.** On the field the deficit is $e(\ell^1/\ell^2)-e(G)
=+0.137891$, and the distance from $e(G)$ to $\theta'/2$ is
$+0.133511$. The two are not interchangeable, because the
concentration exponent is $4.98$ standard errors above $\theta'/2$ and
comes down to it only when $F$ saturates — which V5 leaves
undetermined. So the demand reads as [rem:leanidentity] wrote it and
not as [rem:flatnessfill] rewrote it: **$e(G)$ must close $+0.137891$
onto the measured $e(\ell^1/\ell^2)=+0.284380$, on the vector $\mu$
actually gives.**

*Added later.* The lever V5 asked for was affordable and has been
pulled; see Remark [rem:fieldreach]. The deficit and the exponents
above move within their errors; what changes is that the shape question
is retired rather than open.


#### Remark (two more octaves, and the shape question retires) {#rem:fieldreach}
<!-- evidence: audit_field_reach.py -->

Remark [rem:fillfield] refuted two of its own predictions in the same
sentence: on a clean field the two shapes of [rem:flatnessshape] are
still tied, and $F$'s octave-wise scatter is at or above the published
eight-point r.m.s. in most octaves. Both said that more points inside
$[2\cdot10^5,\,2.56\cdot10^7]$ cannot say where $F$ saturates — only a
longer lever can. The lever costs seconds per $N$: the field runs to
$1.024\cdot10^8$, which is $81$ on-field $N$ against $60$ and a spread
in $\log N$ of $6.2383$ against $4.8520$. The $60$ rows already
published reproduce with a worst departure of $0$ on a count and
$0.000050$ on a printed ratio (W1).

**The excess over the ceiling is not a small-sample artefact: it
sharpens with every lever.** $e(\ell^1/\ell^2)=+0.283586$ with a
standard error of $0.000578$, so the excess over $\theta'/2$ is
$+0.003586$ — **$6.20$ standard errors**, after $3.15$ on eight
doublings and $4.98$ on sixty (W2). The measured ceiling agrees:
$e(\#k)=+0.560357\pm0.000068$, half of it $+0.280178$.

**$F$ is resolved rising over the sweep and not at its top.** Globally
$e(F)=+0.003408$ with a standard error of $0.000576$, $t=5.92$ (W3),
and $F$ runs $0.6622$ to $0.7017$. But the topmost octave's local slope
is $+0.008112$ at $t=1.16$ (W6 refuted), and the two octaves below it
are $-0.006633$ and $-0.006279$. The rise is a fact about the range and
not about the neighbourhood of any $N$ in it.

**The shape question is retired, not left open.** Refitted on the $81$:
the power law gives r.m.s. $0.006380$ and $a+b/\log N$ gives
$0.006366$ — the bounded shape now marginally ahead, the reverse of
before and equally meaningless, because the gap is $0.000015$ against
the r.m.s.'s own standard error $0.000506$ (W4 refuted). Two octaves of
extra lever moved the ratio of gap to error from $0.07$ to $0.03$: the
shapes are not separating, they are converging, which is what two
curves do when they agree on the whole interval where data exist and
differ only where $F$ approaches its bound. **So where $F$ saturates is
not decidable by computation at any reach this programme can afford**,
and this repository publishes no arrival point for $F=1$ — not
[rem:flatnessshape]'s $10^{28.6782}$, and none from
$e=+0.003408$ either.

**What that costs item 4(b), and what it does not.** The deficit is
stable: $e(G)=+0.149567\pm0.002424$ and
$e(\ell^1/\ell^2)-e(G)=+0.134019$ against the $+0.137891$ measured on
$60$, a move of $-0.003872$ against a combined standard error of
$0.002492$, which is $1.55$ of it (W5). So the demand is unchanged in
size. What the retirement changes is the standing of the two ends: the
concentration exponent is $6.20$ standard errors above $\theta'/2$ and
the only thing that would bring it down is a saturation whose location
is now known to be uncomputable here. **$e(G)$ must therefore close
onto a measured $+0.283586$ and not onto an arithmetic $0.28$**, and no
future reach of this family will settle which of the two it is really
closing onto.

*Added later.* One sentence above is withdrawn — "the shapes are not
separating, they are converging". A shrinking ratio of gap to error
cannot be read that way from a test with no power, and the test has
none; see Remark [rem:shapepower]. The retirement itself survives and
is strengthened.


#### Remark (the tie was never a measurement) {#rem:shapepower}
<!-- evidence: audit_shape_power.py -->

Three cycles have ended on one sentence: direction resolved,
destination undetermined. [rem:fieldreach] on $F$, [rem:splitreach] on
the head's one-sign fraction, [rem:headaxis] on the imbalance axis's
grip. In each the trend clears its own error many times over and the
shape carrying it does not, so nothing is published about where it
goes. Rule 3 of this repository exists for exactly the question nobody
put to that: **a tie between two shapes is evidence of similarity only
if the test could have told them apart.**

It could not. Feeding the comparison $2000$ synthetic sweeps drawn
from the fitted power law itself, at the observed abscissae and the
observed residual r.m.s., it declares a separation in $0.0000$ of them
and picks the generating shape in $0.0000$ (P2). Drawing from the
bounded shape instead gives the same. **The discriminator this
repository has used four times — the gap in r.m.s. against
$\min(\mathrm{r.m.s.})/\sqrt{2(n-2)}$ — has no power at this design at
all**, and the scatter it is run at does not matter: at the loudest
octave the sweep shows and at the quietest, the answer is the same
$0.0000$.

**The shapes are distinguishable; the noise is thirteen times their
difference.** With the noise set to zero the bounded shape leaves
r.m.s. $0.00047608$ against the power law's exact fit (P4), and the
greatest pointwise gap between the two fitted curves is $0.001195$, at
$\log N=12.2061$ — $0.1877$ of a single residual. So the two curves
are not the same function on this window; they are nowhere near far
enough apart to be told apart through the scatter $F$ actually
carries.

**What span would do it.** Holding the family's density of $12.9843$
points per unit $\log N$ and extending upward, the power runs
$0.0000$ at $10^{8.0103}$ — the reach [rem:fieldreach] achieved —
$0.0000$ at $10^{9}$, $0.0075$ at $10^{10}$, $0.7075$ at $10^{12}$,
and first clears $0.95$ at $10^{15}$ (P3). Sieving $\Lambda$ and $\mu$
to $10^{15}$ is not a matter of patience.

**So the retirement in [rem:fieldreach] is right and its reasoning was
not.** That remark read the gap-to-error ratio falling from $0.09$ to
$0.07$ to $0.03$ as the shapes converging. With power $0.0000$ that
ratio is not a measurement of anything, and the sentence is withdrawn.
What replaces it is stronger: the question was never being asked. Four
statements of the form "the shapes are tied" in this repository —
[rem:flatnessshape], [rem:flatnessfill], [rem:fillfield],
[rem:fieldreach] — are all correct as refusals to conclude and none of
them is evidence that the two shapes describe $F$ equally well.

**And the control refuted its own tolerance, again.** P1 was
registered at $0.000001$ and the refit departs by $0.00000812$. The
table prints four decimals, so each $F$ is within $0.000050$ of the
value that produced it; an r.m.s. is nonexpansive under that
perturbation and the slope of $\log F$ is a fixed linear functional
with sensitivity $0.000036$. Every departure is inside its own bound.
This is [rem:slopes]'s M1 a second time — **a control's tolerance that
does not ask what the printed table carries refutes the tolerance, not
the fit** — and the script gates on the computed bound and prints that
it does.

One thing this does not measure. [rem:splitreach]'s one-sign fraction
and [rem:headaxis]'s spreads rest on the same kind of comparison and
have not been put through this test. What transfers to them is the
question, not the number.

*Added later.* They have now; see Remark [rem:destination]. The
numbers transfer, and asking for the destination instead of the shape
gets one answer the shape contest could not give.


#### Remark (asking for the destination instead of the shape) {#rem:destination}
<!-- evidence: audit_power_reach.py -->

[rem:shapepower] left two quantities untested and one method
unexamined. Both are fixed here. The controls reproduce inside the
bound each table's printing forces, computed and not assumed
(Q1) — the lesson that remark paid for, applied at the start this
time.

**The transfer holds.** Drawing from each fitted linear shape at its
own abscissae and residual r.m.s., the discriminator picks the
generating shape in $0.1375$ of trials for the head's one-sign
fraction and $0.0000$ for both axis spreads (Q2). Three quantities,
three shape questions, none answerable at the reach computed.

**But the shape contest is the wrong question.** Nobody needs to know
which of two curves fits; they need to know where the trend goes.
That is an asymptote, and an asymptote has a standard error whether or
not a contest has a winner: fit $y=L+b\,g(\log N)$ for bounded $g$ and
$L$ is the destination. The catch is that $L$ depends on $g$, and that
is the shape ambiguity in the one place it can be priced — the spread
of $L$ across bases against the error within one.

**Priced, it is the basis and not the noise** (Q3). On the bases
$1/x$, $1/\sqrt{x}$, $1/x^2$:

| quantity | $L$ across bases | spread | $2\times$ largest s.e. |
|---|---|---|---|
| one-sign fraction | $+0.462564,\ +0.042611,\ +0.672482$ | $0.629871$ | $0.073194$ |
| $\lvert I\rvert$ spread | $-0.061735,\ -0.557187,\ +0.185867$ | $0.743054$ | $0.177565$ |
| $T$ spread | $+0.578516,\ +1.094751,\ +0.320392$ | $0.774360$ | $0.151457$ |

Each destination is pinned to a few hundredths within a basis and
undetermined by six or seven tenths across them. **The noise is not
what leaves these open; the choice of curve is**, and no amount of
extra $N$ at this density changes that — it is the same wall
[rem:shapepower] measured, stated in the units the question was asked
in.

**Q4 is refuted, and it costs a reading.** Two of the three bases put
the head's one-sign asymptote at or below one half: $1/x$ gives
$+0.462564$ with two-sigma interval $[+0.425102,+0.500025]$ and
$1/\sqrt{x}$ gives $+0.042611$; only $1/x^2$ clears, at
$[+0.652438,+0.692527]$. **A coin is inside what the data allow.**
[rem:splitreach]'s sentence — at every $N$ this programme can compute,
the head's alignment is not decay to randomness — stands exactly as
written, because it was written about the measured range. What is
excluded is the extrapolation: nothing here entitles anyone to say the
head keeps a positive majority.

**One thing the destinations do agree on.** [rem:headaxis] published
no crossing for the handover from imbalance to mass. On all three
bases the $\lvert I\rvert$ destination lies below the $T$ destination
— gaps $-0.640250$, $-1.651938$, $-0.134524$, each clearing its own
two-sigma — so **the handover completes on every basis tried**, and
only its size is undetermined, by a factor of $12.28$. That is the
first destination this repository has been able to state, and it is a
direction rather than a place.


#### Remark (which of these slopes is above its own noise) {#rem:slopes}
<!-- evidence: audit_slope_significance.py -->

Remark [rem:kexponent] retired a *span* by measuring the scatter it
had to beat. The same question was never put to the **slopes**, and
four of them are load-bearing. A least-squares slope's standard error
is available from residuals the tables already print,
$\mathrm{s.e.}=\sqrt{\sum r^2/(n-2)\,/\sum(x-\bar x)^2}$, so no new
arithmetic is needed to put all four on one footing.

M2 holds and M3 holds; **M4 is refuted.** The primorial ladder's level
slope stands at $11.99$ standard errors and its $k$-exponent at
$5.03$ — the "six times its scatter" of Remark [rem:kexponent],
computed properly. The share's $\log\rho$ reaches $2.87$. **But the
family's level exponent — rule U4's entire basis — reaches $1.60$,
and its two-sigma interval is $[-0.001169,+0.010489]$, which contains
zero.** The margin over $\tfrac12$ may be opening at the rate U4
reports, or closing; five points cannot tell.

The argument U4 actually offered makes the failure sharper. Its three
leave-one-out refits span $0.000289$ to $0.007892$, which lies
*inside* that two-sigma interval — as it must, since dropping one
point of five is a bounded perturbation of a least-squares fit
whatever the noise. **Agreement among leave-one-out refits is a fact
about the arithmetic of least squares and not about the data**, and it
was read here as if it were evidence.

The control M1 is also refuted, at a tolerance of $10^{-5}$ set
without asking what the printed tables carry. That is answerable
rather than worrying: a slope is a fixed linear functional
$\sum c_iy_i$ of the ordinates, so rounding the $y$ column to $d$
decimals moves the refit by at most $\sum|c_i|\cdot\tfrac12 10^{-d}$.
Every M1 gap lies inside its own bound — the largest, $0.0000321$
against $0.0000433$ — so what M1 refutes is the tolerance, and no
published slope is in doubt from it.

What this costs: **rule U4 is withdrawn as a reading.** Its verdict
under its own registered rule was "hold", the slope being positive,
and that rule was the badly chosen thing; the sign it turns on is not
resolved. Remark [rem:residuelevel]'s $0.5654$–$0.5799$ is unaffected
— those are measurements, not a trend — but nothing may be inferred
from them about where the margin goes. The ladder's rise, which is
what Remark [rem:primorialrung10] rests on, is untouched and is the
only level trend this project has that clears its own noise.

*Added later.* That last sentence was true when written and is no
longer. What the table above does not say is over how much range each
$t$ was earned, and the family's five points span $2.7726$ in
$\log N$ against the ladder's $6.9315$. Two further octaves resolve
the family's sign as well; see Remark [rem:slopereach].


#### Remark (two more octaves settle it) {#rem:slopereach}
<!-- evidence: audit_level_slope_reach.py -->

Remark [rem:slopes] left one thing open rather than closed. The
family's level slope was unresolved at five points — but a slope's
standard error falls with the spread of the abscissae, and five
octaves of $N$ give a spread in $\log N$ of only $2.7726$ against the
primorial ladder's $6.9315$. **Unresolved because the statistic is
noisy and unresolved because the sweep is short are different
complaints, and only the second can be answered by computing more.**
It was the second.

An independent reimplementation — the sieve weight applied by a
precomputed residue bitmask over $N-mk$, with $C_k$ factored out of
the sum rather than carried in the weights — reproduces all five
published exponents. **N1 holds to $0.0000$ at every one of the
five**, which is the strongest control this project has: two
independently written codes agreeing on $K^*_R$ exactly.

Pushed to $N=6.4\cdot10^6$ and $N=1.28\cdot10^7$, **N2 holds**: the
exponents are $0.5767$ and $0.5834$, still clear of $\tfrac12$ and now
clear of $\theta'=0.56$ at both. Over all seven the spread in $\log N$
is $4.1589$ and **N3 and N4 hold**: the slope is $+0.005112$ at
$3.71$ standard errors, two-sigma interval $[+0.002353,+0.007871]$.

**So U4's reading is restored, on seven points and not on five.** The
margin over $\tfrac12$ is opening, at about $0.005$ per unit
$\log N$ across the accessible range. What restored it was range and
not argument: the slope barely moved, $+0.004692\to+0.005112$; the
standard error fell by the factor $0.4737$, of which $0.5976$ is the
extra spread and the rest is scatter, which fell from $0.0049$ to
$0.0043$.

Two things this does not say. It is still a measurement over a factor
$64$ in $N$ and forecasts nothing (Remark [rem:forecastbracket]).
And it is still one odd radical: Remark [rem:residuearithmetic]'s
primorial-like $N$ are not in this sweep, and Remark
[rem:primorialrung10] is the corresponding statement for those.


#### Remark (the split constant is free and buys nothing) {#rem:betafree}
<!-- evidence: audit_beta_optimal.py -->

Every exponent in Remark [rem:residuelevel] depends on $\beta$, and
$\beta$ is not given by the problem. $H=\beta P+R$ holds for any
$\beta$ and so does $|H|\le\beta|P|+|R|$; granting
Remark [rem:provablehalf] the conditional bound is
$B(N)\le B_R(\beta;N)+o(N)$ **for any $\beta$ whatever**, so the
argument is free to choose it. What it has instead is
$\sum HP/\sum P^2$, which minimises the $\ell^2$ distance from $H$ to
the ray $\beta P$ — and the budget is spent by
$\sum_{k<K}(\log k)|H-\beta P|$, an $\ell^1$ norm weighted by $\log k$
and truncated at the crossing. Two different problems, and nobody had
checked how far apart their answers are.

**They are not apart, and the freedom is worth nothing.** The control
T1 reproduces the operative exponents to $5\cdot10^{-5}$. Rule T2
asked the level to be sensitive to $\beta$ at all and **is refuted**:
over $\pm10\%$ the exponent moves by
$0.0092,\,0.0022,\,0.0026,\,0.0018,\,0.0016$ — at the four larger $N$
by less than a fifth of what the rule called detectable. Maximising
$K^*_R$ over $241$ values of $\beta$ from $0.25$ to $6.25$ moves it
from $993,\,1447,\,2019,\,3319,\,5923$ to
$1023,\,1447,\,2019,\,3323,\,5923$ — at three of the five $N$ not at
all. **T4 is refuted**: at $N=8\cdot10^5$ the exponent is still
$0.5599$, so the one $N$ that missed $\theta'=0.56$ misses it under
every $\beta$ the argument could have chosen.

T3 holds — the optimum is more than $1\%$ from the fit at three of
five — but the diagnostic shows why that is not a finding. At four of
the five $N$ the $\ell^2$ residual and the truncated budget sum both
read $1.000$ to three decimals at the optimum: **the two objectives
are not merely agreed, they are flat in the same place.** Only at the
smallest $N$ do they separate, where the optimal $\beta$ is worse by
$\ell^2$ and better by $\ell^1$ and buys $0.0024$ in the exponent.
(The optimum is also coarse in a way worth stating: $K^*$ is
integer-valued and jumps, so an interval of $\beta$ shares one $K^*$
and a $\beta$-ratio of $0.9829$ can mean the identical crossing.)

So the knife-edge of Remark [rem:residuelevel] is a property of the
object and not of a fitted constant. **That is the useful reading of a
negative result: the $0.5599$ cannot be tuned away, and neither can
the $0.06$ of margin over $\tfrac12$.**


#### Remark (which half exhausts the budget) {#rem:splitbudget}
<!-- evidence: lab_split_budget.py -->

Remark [rem:residue] found $|R|\asymp(N/k)^{1/2}$ and the program's
open list recorded the next task as proving square-root cancellation
for $R$. That task is worth its cost only if $R$ is what binds. By the
triangle inequality $B_H\le\beta B_P+B_R$, and each half has its own
crossing of the budget $\SS(N)N$.

**It is not $R$ that binds.** At the operative truncation $K^*_H$ the
elementary half already takes
$0.8112,\,0.8573,\,0.8484,\,0.8691,\,0.8456$ of the budget while the
residue takes $0.4737,\,0.5119,\,0.5062,\,0.5214,\,0.5275$. Removing
the residue *entirely* would move the truncation from
$K^*_H=2973,\ldots,23397$ only to $K^*_P=4119,\ldots,30369$ — a factor
$1.2402$ to $1.3855$, which in exponent is a move from
$0.6552$–$0.6716$ to $0.6799$–$0.6891$, about $0.02$ in $\theta'$.
Removing the elementary half instead would reach
$K^*_R=9191,\ldots,63399$, exponent $0.7371$–$0.7477$, about three
times as much. **So the wall is the sieve-weighted Möbius sum, not the
residue**, and a perfect result about $R$ moves $\theta'$ from $0.67$
to $0.69$.

Two qualifications, and each is now declared by the result file
itself. **Every figure in this remark is against the budget
$\SS(N)N$, which is not the one the route asks for.** Proposition [prop:nolog] needs
$\SS(N)(1-A(N))N$, smaller by $4.7009$ and worth about $0.18$ in the
exponent. The comparison *between* halves above is unaffected — both
are crossed against the same thing — but the exponents themselves are
not the ones $\theta'$ is measured against, and Remark
[rem:residuelevel] recomputes $K^*_R$ against the operative budget,
where it reads $0.56$–$0.58$ and not $0.74$.

**And every $N$ here is $2^a5^b$: one odd radical.** The *ordering* of
the two halves — that $\beta P$ takes the larger share — is what this
remark is for, and an ordering is more robust to arithmetic than a
level; but the exponents themselves belong to that radical, and Remark
[rem:residuearithmetic] shows the residue's own level crossing
$\tfrac12$ downward as the radical changes.

The pre-registered **Z4** capped $K^*_P/K^*_H$ at $1.30$ and it reads
$1.3855$ at the smallest $N$, so Z4 fails there and holds at the other
four ($1.2402$ to $1.2980$). The refutation was flagged in advance as
the outcome that would vindicate the recorded task; at one $N$ out of
five, and by four percent, it does not.

The control sharpens where the budget goes. Permuting each half's
magnitudes across $k$, weights left in place, the crossings move
*outward* by a factor of three to five — $K^*_P=4119$ against a band
of $[17627,\,20017]$ at the smallest $N$, and similarly throughout,
with narrow bands. So the crossings are strongly a property of the
pairing, and in the direction that matters: $|P|$ and $|R|$ are
largest at the *smallest* $k$, where the inner sum is longest, so the
measured cumulative front-loads and reaches the budget while a random
pairing is still climbing. **The budget is consumed at the bottom of
the $k$-range**, which is exactly where the reduction has no freedom
to move.


#### Remark (the control the ratio could take) {#rem:pairingnull}
<!-- evidence: audit_weightgap_pairing.py -->

Remark [rem:weightgapnull] found the coin unusable for two of the
weight-gap statistics — it drives both $|\sum H|$ and
$|\sum(\log k)H|$ to square-root size, so their ratio is a quotient of
two near-zero quantities — and left them standing uncontrolled. There
is a control that keeps them well conditioned: **permute $H(N;k)$
across $k$**, leaving each weight attached to its own $k$. A
permutation does not change a plain sum, so $\sum H$ is exactly
invariant — measured, it drifts by $2.291\cdot10^{-16}$ at worst — the
ratio's numerator is pinned, and the only thing destroyed is the
pairing between a modulus and its dilated wall.

**Under it the effective modulus does not survive.** The profile in
$j$ comes out *more* geometric for the permutations than for $\mu$:
their spread of consecutive ratios has median
$1.0054,\,1.0052,\,1.0053,\,1.0036,\,1.0033$ against $\mu$'s
$1.0146,\,1.0201,\,1.0132,\,1.0099,\,1.0220$. Near-geometry in $j$ is
therefore what an *unpaired* weighting gives, not evidence of
concentration — and $k^*$ read off that profile tracks the plain
geometric mean of the $k$-range, $\log k^*$ running
$5.5347$ to $6.8669$ against $\overline{\log k}$ of $5.8487$ to
$7.3930$, a gap of $0.0244$ to $0.0670$ of $\log K$.

The ratio itself does carry a little pairing information, in the
direction opposite to the withdrawn reading. $\mu$'s
$|\sum H|/|\sum(\log k)H|$ is $0.1807,\,0.1740,\,0.1624,\,0.1389,\,
0.1456$ against permutation maxima of
$0.1790,\,0.1695,\,0.1547,\,0.1481,\,0.1386$ — above the band at four
of the five $N$, by about four percent. So the $\log k$ weight
amplifies $\mu$'s sum slightly *less* than it would amplify the same
values in random order.

All four pre-registered rules fail, on their bands rather than their
content. **W2** and **W3** asked $\mu$ to sit *inside* the
permutations' range and it sits just outside at four of five $N$ —
the outcome flagged in advance as the informative one. **W4** capped
the gap at $0.05$ of $\log K$ and it reaches $0.0670$. **W1** capped
the permuted ratio's spread at $0.05$ of its median and it reaches
$0.0881$; the invariance half of W1 held exactly, and for comparison
the coin of Remark [rem:weightgapnull] spanned $[0.0000,\,8.0955]$ on
the same statistic, so this control is better conditioned by two
orders even where the cap was too tight.


#### Remark (where the weight puts its mass) {#rem:weightgap}
<!-- evidence: lab_weight_gap.py -->

[eq:flatsum] makes it cheap to turn the weight on continuously, through
$w_k=(\log k)^j$. At $N=2.56\cdot10^7$,
$|\sum_k(\log k)^jH|/N$ reads $0.01492$, $0.02545$, $0.04335$,
$0.07379$, $0.12565$ at $j=0,\,0.25,\,0.5,\,0.75,\,1$ — a factor
$1.7034$ per quarter step, with the five ratios spread by only
$1.0018$. Geometric in $j$ is what a sum concentrated at a *single*
effective modulus $k^*$ would give, $\sum(\log k)^jH\approx
(\log k^*)^j\sum H$; and $k^*$ read off that way sits at
$\log k^*/\log K=0.8814$ at the top $N$, between $0.7959$ and
$0.9368$ across the range.

**That reading is withdrawn.** An earlier version concluded from it
that the $\log k$ weight reweights the sum towards $k\approx K^{0.88}$,
the top of the range, and joined that to Lemma [lem:extract]'s
requirement on where a weight must put its mass. Remark
[rem:pairingnull] supplies the control that was missing, and it does
not survive: permuting $H$ across $k$ leaves $\sum H$ exactly fixed
and destroys only the pairing, and under that permutation the profile
comes out *more* geometric than $\mu$'s while $k^*$ tracks the
geometric mean of the $k$-range. What is left of this remark is the
table above and the fact that the log-weighted sum exceeds the flat
one by a factor near $8$ — not any statement about concentration.


### The two opposing constraints


#### Lemma (extraction) {#lem:extract}
<!-- evidence: audit_extraction_tradeoff.py -->

$\displaystyle B_w=\sum_{\substack{d<K\\(d,N)=1}}b_d\,
 \frac{\mu(d)}{\varphi(d)}\,\rho_{dN}\!\left(\frac{K}{d}\right)$,
where $\rho_{n}(x)=\sum_{j<x,\,(j,n)=1}\mu(j)/\varphi(j)$ satisfies
$\rho_n(x)\ll e^{-c_1\sqrt{\log x}}$ **whenever $\log n\ll\log x$**.


**Proof.** 
Write $k=dj$; for squarefree $k$ one has $(d,j)=1$,
$\mu(k)=\mu(d)\mu(j)$ and $\varphi(k)=\varphi(d)\varphi(j)$. The bound
on $\rho$ is Huang–Li's Lemma 1, whose hypothesis $\log n\ll\log x$ is
carried above rather than dropped. It is not decoration: the table of
§"Numerical confirmation of the load-bearing quantity" prints
$|B_w|\varphi(d_0)=1.000000$ at $\lfloor K/d_0\rfloor=1$ — no saving at
all — and that is the corner where $\log n\asymp\log N$ while
$\log x=O(1)$. Theorem [thm:D] is unaffected, because it maximises only
over $d\le N^{\theta'-1/2-\delta}$, where $x=K/d\ge N^{1/2+\delta}$ and
$n=dN\le N^{3/2}$, so $\log n\asymp\log x$ holds automatically.
 ∎


#### Lemma (BV-accessibility) {#lem:bv}
<!-- evidence: analytic -->

Expanding $w_k=\sum_{d\mid k}b_d$ inside the residual gives
$\mathcal R_w=\sum_d b_d\,\mathcal R^{(d)}$, where $\mathcal R^{(d)}$
is a sum of $\Lambda$ over progressions to moduli $md$ with
$m<N^{1-\theta'}$. Bombieri–Vinogradov bounds each
$\mathcal R^{(d)}\ll_A N(\log N)^{-A}$ provided
$d\le N^{1/2-\delta}/N^{1-\theta'}=N^{\theta'-1/2-\delta}$; hence, on
that support, $\mathcal R_w\ll_A \|b\|_1\,N(\log N)^{-A}$.


The two lemmas pull in opposite directions. Lemma [lem:extract]
says $B_w$ is controlled by the part of $b$ sitting *at* the
truncation point $K$: mass at $d\le K^{1-\varepsilon}$ is damped by
$e^{-c\sqrt{\varepsilon\log K}}$. Lemma [lem:bv] says $b$ may only
live *below* $N^{\theta'-1/2}$. The two thresholds are separated
by a factor $N^{1/2}$.


#### Theorem (no weight extracts $C(N)$) {#thm:D}
<!-- evidence: analytic -->

Fix $\theta'\in(1/2,1)$ and $\delta>0$, and let $w$ be any weight
whose transform $b=\mu*w$ is supported in
$[1,\,N^{\theta'-1/2-\delta}]$ — the range in which the residual is
accessible to Bombieri–Vinogradov. Then

$$
\frac{\|b\|_1}{|B_w|} \;\gg\;
   \exp\!\left(c_1\sqrt{(\tfrac12+\delta)\log N}\right),
$$

and consequently, *even granted* $T_w \ll_A N(\log N)^{-A}$ for every
$A>0$ — what Theorem [thm:A] supplies unconditionally at $w=1$, and
the most $EH_\mu$ would supply for a general $w$ (see
[rem:extractTw]) — [eq:extract] yields at best

$$
|C(N)| \;\ll_A\; \exp\!\left(c_1\sqrt{\tfrac12\log N}\right)
   \cdot \frac{N}{(\log N)^{A}},
$$

which is not a saving of any power of $\log N$. No weight in the
design space extracts $C(N)$ by divisor switching plus
Bombieri–Vinogradov.


#### Remark (relation to Bombieri's asymptotic sieve) {#rem:bombieri}

The genre of Theorem [thm:D] is classical. Bombieri's asymptotic
sieve [Bom76] shows that sieve weights alone cannot detect primes,
however the weights are chosen, and the parity obstruction it isolates
is the reason one expects a statement of this shape. What
Theorem [thm:D] adds is not the phenomenon but its precise form
here: the obstruction is quantified as a separation of $N^{1/2}$ between
two thresholds — the point at which a weight's transform must carry
mass in order to extract $C(N)$, and the point below which the residual
is accessible at all — it is specific to the divisor-switch route, and
it survives granting the full Elliott–Halberstam conjecture rather than
merely Bombieri–Vinogradov. The reason the statement is worth making is
narrow but real: without it, $E_3$ looks like one member of a family of
candidate consumptions, and Theorem [thm:C] then looks like an
accident of one choice of weight rather than a property of the route.


**Proof.** 
By Lemma [lem:extract] and $\varphi(d)\ge1$,

$$
|B_w| \le \sum_{d\le D}\frac{|b_d|}{\varphi(d)}
    \left|\rho_{dN}\!\left(\frac Kd\right)\right|
  \le \|b\|_1\max_{d\le D}\left|\rho_{dN}\!\left(\frac Kd\right)\right|
  \ll \|b\|_1\,e^{-c_1\sqrt{\log(K/D)}},
$$

with $D=N^{\theta'-1/2-\delta}$, so that
$K/D=N^{1/2+\delta}$. This is the displayed ratio.

For the second assertion, solve [eq:extract] for $C(N)$ and bound the
three terms on the right: the complete part is $\ll\|b\|_1\log N$ by
Lemma [lem:Gb] and $\Lambda\ll\log N$; the residual is
$\ll_A\|b\|_1N(\log N)^{-A}$ by Lemma [lem:bv], the support
hypothesis on $b$ being exactly what that lemma needs; and
$T_w\ll_AN(\log N)^{-A}$ is the granted hypothesis. Dividing by
$|B_w|$ and using the displayed ratio gives the bound on $|C(N)|$; and
$e^{c_1\sqrt{\log N/2}}$ exceeds every power of $\log N$.
 ∎


#### Remark (what the obstruction is)

The loss factor is $\exp(c_1\sqrt{(1/2)\log N})$, and the $1/2$ in the
exponent is literally the exponent of the $\sqrt N$ barrier: it is the
gap between where a weight must live to see $C(N)$ (at $N^{\theta'}$,
the truncation point) and where it may live for Bombieri–Vinogradov
to close its residual (below $N^{\theta'-1/2}$). Theorem [thm:C]
exhibited the two endpoints of this space; Theorem [thm:D] shows
the interior is empty, and shows why.


#### Remark (scope, stated honestly)

Theorem [thm:D] is a no-go for one precisely specified method —
divisor switching with Bombieri–Vinogradov as the only input — over
the whole weight space of that method. It is not, and does not claim
to be, an obstruction to other methods. Its value is that it closes
the design space that Theorems [thm:A] and [thm:C] opened,
rather than leaving it to be re-explored.


### The no-go survives Elliott–Halberstam


The obvious objection to Theorem [thm:D] is that
Bombieri–Vinogradov is not the only available input: Huang–Li's
Theorem 1 itself *assumes* $EH$ for $\Lambda$ at level
$N^{\theta}$. Raising the assumed level does not help.


#### Theorem {#thm:Dprime}
<!-- evidence: analytic -->

Suppose $\Lambda$ has level of distribution $\theta_E \in (0,1)$, i.e.
$EH(N^{\theta_E})$ holds. Then the residual is accessible only for $b$
supported on $d \le N^{\theta_E-(1-\theta')}$, while extraction still
requires mass at $d \asymp K = N^{\theta'}$; the two thresholds are
separated by $N^{1-\theta_E}$, and

$$
\frac{\|b\|_1}{|B_w|} \;\gg\;
  \exp\!\left(c_1\sqrt{(1-\theta_E)\log N}\right),
$$

which exceeds every power of $\log N$ for each fixed $\theta_E<1$.


**Proof.** 
Identical to Theorem [thm:D] with $N^{1/2-\delta}$ replaced by
$N^{\theta_E}$: the residual's moduli are $md$ with $m<N^{1-\theta'}$,
so $d \le N^{\theta_E}/N^{1-\theta'}$, whence
$K/D = N^{\theta'}/N^{\theta_E-1+\theta'} = N^{1-\theta_E}$, and
Lemma [lem:extract] gives the ratio.
 ∎


#### Remark

Closing the gap would require $\theta_E = 1$ exactly: equidistribution
of $\Lambda$ in progressions to moduli of size $N$ itself, where each
progression contains $O(1)$ terms and the statement carries no
information. So the demand side stays closed *even granting the
full Elliott–Halberstam conjecture* — the divisor-switch route needs
not a stronger level but a different mechanism.


### Smooth weights: $\mu * \log^D = \Lambda_D$ {#sec:Dpp}


Theorem [thm:D] closes the weights whose transform $b=\mu*w$ is
supported low enough for Bombieri–Vinogradov. That hypothesis
excludes the most natural family of all, $w_k = f(\log k)$ with $f$ a
polynomial, whose transform is spread over every scale. For those the
structure is explicit:

$$
b \;=\; \mu * \log^D \;=\; \Lambda_D,
$$

the generalised von Mangoldt function. $\Lambda_D$ vanishes on
integers with more than $D$ prime factors, and on a squarefree
$u=p_1\cdots p_r$ with $r\le D$ it is the $r$-th finite difference of
$x^D$ at the points $\log p_i$; in particular $\Lambda_D(u) =
D!\,\log p_1\cdots\log p_D$ when $r=D$. Hence for $f = x^D$ the
complete part splits by the number of prime factors,

$$
\mathrm{CP}_D(N) = \sum_{r=1}^{D}\ \sum_{\substack{u<N\ \text{squarefree}\\ \omega(u)=r}}
   \Lambda(N-u)\,\Lambda_D(u),
$$

the $r=1$ piece being a Goldbach-type sum and the $r\ge2$ pieces
Chen-type sums (prime plus $r$-almost-prime).


#### Proposition (polynomial weights) {#prop:Dpp}
<!-- evidence: analytic -->

Let $w_k = f(\log k)$ with $f=\sum_{j\le D}c_jx^j$ real.


- **(i)**  If $D=0$ then $B_w \ll e^{-c\sqrt{\log K}}$ (Huang–Li's
  Lemma 1): no extraction.

- **(ii)**  If $f$ is a single monomial $x^D$, $D\ge1$, then every
  term of $\mathrm{CP}_D$ is nonnegative, since $\Lambda\ge0$ and
  $\Lambda_D\ge0$; by classical sieve bounds
  $\mathrm{CP}_D \asymp N(\log N)^{D-1}$, with a fixed sign. It is
  never $o(N)$.

- **(iii)**  Cancellation therefore requires coefficients of opposite
  sign, i.e. tuning $c_1$ against $c_2,\dots$; but the ratio that
  would have to be matched involves the asymptotics of the $r=1$
  piece, which for $D=1$ *is* the binary Goldbach sum. The tuning
  is circular.


Consequently no polynomial weight makes the complete part a
classically known object.


#### Remark (a prediction of ours that the measurement refuted) {#rem:toprdom}
<!-- evidence: audit_polyweight.py -->

We first expected the top-$r$ piece to dominate the others by a power
of $\log N$, which would have given (ii) without any appeal to
nonnegativity. That is false: measured at $N = 10^6, 4\cdot10^6,
1.6\cdot10^7$ (`code/audit\_polyweight.py`), the two pieces of
$\mathrm{CP}_2$ are the *same* order, with

$$
\begin{array}{r|ccc}
 N & 10^6 & 4\cdot10^6 & 1.6\cdot10^7\\\hline
 r=1:\ \sum_p\Lambda(N-p)\log^2p\ /\,N & 22.5145 & 25.0438 & 27.4579\\
 r=2:\ 2\sum_{pq}\Lambda(N-pq)\log p\log q\ /\,N
   & 17.3648 & 19.7926 & 22.2512\\
 \text{ratio } r_2/r_1 & 0.7713 & 0.7903 & 0.8104\\
 \mathrm{CP}_2/(N\log N) & 2.8866 & 2.9494 & 2.9967
\end{array}
$$

— the ratio drifts towards $1$, not towards $0$ or $\infty$. The
closure stands, but on nonnegativity rather than on separation of
scales. Calibration: the $D=1$ column reproduces the Goldbach sum
$\sum_p\Lambda(N-p)\log p = 1.7565N,\,1.7633N,\,1.7614N$ against
$\SS(N)=1.76043$, which is the same at all three $N$ because
$P(N)=\{2,5\}$ for each.


#### Remark (the last table was truncated, not rounded) {#rem:trunc}
<!-- evidence: audit_polyweight.py -->

Version 3 printed the last row of the table above as
$2.886,\,2.949,\,2.997$. The audit pre-registered, as its rule W4,
that each entry should sit within half a unit in the last printed
place; the first entry fails, because $2.8866$ was truncated to
$2.886$ where rounding gives $2.887$. Nothing turns on it — the row
is a shape and the shape is unchanged — and it is recorded only
because it is the one entry in this section that a reader checking to
the printed precision would find wrong, and because the convention a
table rounds by is not visible from the table.


#### Remark (the canonical tuning moves the wrong way)

Since $b=\mu*w$ we have $w=1*b$ and hence $B(s)=W(s)/\zeta(s)$. That
series has a *double* pole at $s=1$ for every degree-$2$ $f$, so
no such $f$ removes it and $b$ cannot be made mean-zero; the most a
degree-$2$ tuning can do is kill the second-order term. With
$f=x^2+cx$ one has $W=\zeta''-c\zeta'$, and since neither $\zeta''$ nor
$\zeta'$ carries a $(s-1)^{-1}$ term,

$$
B(s) \;=\; \frac{2}{(s-1)^2} \;+\; \frac{c-2\gamma}{s-1} \;+\; O(1),
  \qquad
  \sum_{u\le x} b_u \;=\; 2x\log x + (c-2\gamma-2)x + o(x).
$$

So the $x$-term vanishes for $f(x)=x^2+(2+2\gamma)x$ and for no other
degree-$2$ choice. Measured to $x=1.6\cdot10^7$
(`code/audit\_polyweight.py`), the residual coefficient
$\bigl(\sum_{u\le x}b_u-2x\log x\bigr)/x$ is $+0.0001$ for that $f$,
$-3.1543$ for the untuned $x^2$, and $-4.3087$ for $x^2-2\gamma x$,
against the predicted $c-2\gamma-2$ of $0$, $-3.1544$ and $-4.3089$.
The computation is $\sum_{u\le x}b_u=\psi_2(x)+c\,\psi(x)$ with
$\psi_2(x)=\sum_{n\le x}\Lambda(n)\log n+\sum_d\Lambda(d)\psi(x/d)$,
which needs no sieve of $\Lambda_2$ itself.

Neither tuning buys smallness, and the canonical one is the worse:
against the untuned $\mathrm{CP} = 39.8793N,\,44.8364N,\,49.7091N$, the
tuning $x^2+(2+2\gamma)x$ *raises* it to
$45.4201N,\,50.3988N,\,55.2654N$ ($+13.89\%,\,+12.41\%,\,+11.18\%$),
while $x^2-2\gamma x$ lowers it to $37.8516N,\,42.8008N,\,47.6756N$
($-5.08\%,\,-4.54\%,\,-4.09\%$). Killing the average of $b$ is not
smallness: what must be small is $b$'s correlation with
$\Lambda(N-\cdot)$, not its mean, and the correlation is the whole
problem.

(Version 2 of the companion paper justified the tuning
$x^2-2\gamma x$ by a pole cancellation in the product $\zeta W$, and
concluded that it makes $b$ mean-zero. The series governing $b$ is the
quotient $W/\zeta$; its pole at $s=1$ is double and not removable by
any degree-$2$ $f$; and that particular $f$ moves the removable term
away from zero rather than to it.)


### The circle method has zero margin on $C(N)$ {#sec:circle}


Since the switch is closed over its whole design space, it is worth
recording the exact position of the other classical mechanism. Write

$$
C(N) = \int_0^1 S_\Lambda(\alpha)\,S_\mu(-\alpha)\,e(-N\alpha)\,
   d\alpha,\qquad
  S_\Lambda(\alpha)=\sum_{n\le N}\Lambda(n)e(n\alpha),\quad
  S_\mu(\alpha)=\sum_{m\le N}\mu(m)e(m\alpha).
$$


#### Proposition (zero margin) {#prop:E}
<!-- evidence: audit_circle_margin.py -->

Both standard ways of estimating this integral lie at or above the
trivial bound $|C(N)|\le\psi(N)\sim N$:


- **(i)**  $\|S_\Lambda\|_2\|S_\mu\|_2
   = \bigl(\sum_{n\le N}\Lambda(n)^2\bigr)^{1/2}
     \bigl(\sum_{m\le N}\mu^2(m)\bigr)^{1/2}
   \sim (6/\pi^2)^{1/2}\,N(\log N)^{1/2}$, exceeding the trivial bound
   by a factor $\asymp(\log N)^{1/2}$, which *grows*;

- **(ii)**  any bound of the shape
  $\sup_\alpha|S_\mu(\alpha)|\cdot\|S_\Lambda\|_1$ is at least
  $\|S_\mu\|_2\|S_\Lambda\|_1 \gg N^{1/2}\cdot N^{1/2}=N$.


The two factors in (ii) are bounded below for different reasons, and
only the first is Parseval's. Parseval gives
$\sup_\alpha|S_\mu| \ge \|S_\mu\|_2 = (6N/\pi^2)^{1/2}$, an
identity-level constraint: no improvement in Möbius exponential-sum
technology can lower $\sup_\alpha|S_\mu|$ below
$(6/\pi^2)^{1/2}N^{1/2}$. For the second factor Parseval runs the
*wrong way*: it gives $\|S_\Lambda\|_1\le\|S_\Lambda\|_2
\ll (N\log N)^{1/2}$, an upper bound. The trivial lower bound
$\|S_\Lambda\|_1 \ge \bigl|\int_0^1 S_\Lambda(\alpha)
e(-n\alpha)\,d\alpha\bigr| = \Lambda(n)$ reaches only $\log N$, short
of $N^{1/2}$ by a full power. The bound $\|S_\Lambda\|_1\gg N^{1/2}$
that (ii) uses is a theorem of Vaughan [Vau88].


Davenport's uniform bound $S_\mu(\alpha)\ll_A N(\log N)^{-A}$ is
therefore useless here: it is a saving over $N$, whereas the pairing
needs a bound at the scale $N^{1/2}$, which Parseval forbids. This is
the parity obstruction in circle-method language — *the binary
problem sits exactly at the trivial bound with no margin* — and it is
why the method that settles the ternary problem cannot be pushed to
the binary one by any sharpening of the Möbius input.

The constants are measured in `code/audit\_circle\_margin.py`
($N=2^{14}\ldots2^{20}$, exact FFT on a $4N$-point grid):

$$
\begin{array}{r|cccc}
 N & 2^{14} & 2^{16} & 2^{18} & 2^{20}\\\hline
 \|S_\mu\|_2/\sqrt N & 0.7798 & 0.7797 & 0.7797 & 0.7797\\
 \sup|S_\mu|/\sqrt N & 3.058 & 2.742 & 2.853 & 2.801\\
 \|S_\Lambda\|_1/\sqrt N & 1.946 & 2.084 & 2.219 & 2.346\\
 \text{(i) bound}/N & 2.297 & 2.473 & 2.639 & 2.795\\
 N/\bigl(\sup|S_\mu|\,\|S_\Lambda\|_1\bigr) & 0.168 & 0.175 & 0.158
  & 0.152
\end{array}
$$

The first row reproduces $\sqrt{6/\pi^2}=0.7797$ exactly, confirming
the computation; the last row — the margin route (ii) would need to
exceed $1$ — sits near $0.16$ and falls over the last three
abscissae. Route (i)
diverges from the trivial bound like $(\log N)^{1/2}$, as predicted.
For reference the object itself is small: $C(N)/N = -0.0105$, $0.0001$,
$0.0059$, $0.0032$ at these $N$.


### Numerical confirmation of the load-bearing quantity


Lemma [lem:extract] carries the theorem, so it is checked directly
(`code/audit\_extraction\_tradeoff.py`, $N = 99{,}999{,}998$,
$\theta'=0.56$, $K=\lfloor N^{\theta'}\rfloor=30199$). For the
single-divisor weights $w_k=[\,d_0\mid k\,]$ (i.e. $b=\delta_{d_0}$)
the lemma is an identity, and it is verified as one: over **all**
$10{,}262$ squarefree $d_0<K$ coprime to $N$, the brute-force $B_w$
over $k<K$ and the factorised $\mu(d_0)\varphi(d_0)^{-1}
\rho_{d_0N}(K/d_0)$ agree to $6.5\cdot10^{-18}$.

$$
\begin{array}{r|cccccc}
 d_0 & 15101 & 3777 & 2159 & 211 & 31 & 1\\
 \lfloor K/d_0\rfloor & 1 & 7 & 13 & 143 & 974 & 30199\\\hline
 |B_w|\,\varphi(d_0) & 1.000000 & 0.750000 & 0.066667 & 0.003384
  & 0.004166 & 0.000529
\end{array}
$$

— of order $1$ only when $K/d_0=O(1)$, i.e. only when the weight's
transform sits at the truncation point, and damped
throughout the range Bombieri–Vinogradov admits (for these parameters
$N^{\theta'-1/2}=3$, where the largest value attained is $0.004145$).


#### Remark (the ratio does not determine the entry) {#rem:ratiotable}
<!-- evidence: audit_extraction_tradeoff.py -->

Version 3 of this note printed the table above indexed by
$K/d_0$ alone. That index does not determine the entry:
$\rho_{d_0N}$ carries the condition $(j,d_0N)=1$, so the value depends
on which small primes divide $d_0$, not only on how far $d_0$ sits
from the truncation point. At $\lfloor K/d_0\rfloor=7$ the $184$
admissible $d_0$ produce exactly four values, split by $\gcd(d_0,15)$:

$$
\begin{array}{r|cccc}
 \gcd(d_0,15) & 1 & 3 & 5 & 15\\\hline
 |B_w|\,\varphi(d_0) & 0.250000 & 0.750000 & 0.500000 & 1.000000
\end{array}
$$

and at $\lfloor K/d_0\rfloor=13$ the $55$ admissible $d_0$ produce six
values spanning $0.066667$ to $0.566667$. The printed $0.750000$ and
$0.066667$ are therefore each one column of several, reproduced by
$37$ of $184$ and by $28$ of $55$ admissible $d_0$ respectively; the
$d_0$ row above is what makes them reproducible. The three right-hand
columns are unambiguous — one admissible $d_0$ each — and the
left-hand column is attained by all $5130$.

The audit also pre-registered a second explanation and **T5 is
refuted**: it predicted that the printed digits were unreachable under
the strict $j<K/d_0$ that Lemma [lem:extract] states, so that the
table had been computed under $j\le K/d_0$ instead. All six columns
are attained under the strict convention. The defect is the index and
nothing else.

This changes nothing in Theorem [thm:D], whose proof uses
$\max_{d\le D}|\rho_{dN}(K/d)|$ and not any particular $d_0$. It is
recorded because a table whose index does not determine its entries
cannot be checked by a reader, and this one was offered as the direct
check of the lemma that carries the theorem.


#### Remark (a pre-registered threshold of ours that failed) {#rem:cap}
<!-- evidence: audit_extraction_tradeoff.py -->

The audit above also pre-registered a cap as its rule T6: no $d_0$
with $K/d_0>1000$
should give $|B_w|\varphi(d_0)>0.05$. It is exceeded — the maximum
over that range is $0.062792$ — and the script exits non-zero on that
account. The cap was chosen as an effect size and not derived from any
null, which is the failure the companion paper's methodology section
names; $\rho\ll e^{-c_1\sqrt{\log x}}$ fixes no constant, so no
numerical cap follows from it. The damping itself is visible —
the octave maxima of $|B_w|\varphi(d_0)$ run
$1.000,\,1.000,\,0.754,\,0.511,\,0.114,\,0.063,\,0.0014,\,0.00053$
over $K/d_0\in[1,2),[4,8),\dots,[16384,32768)$ — but a properly formed
test of the *form* $e^{-c\sqrt{\log x}}$ is not offered here and
the pre-registered verdict stands as failed.


## Numerical verification


All computations are reproducible from the accompanying repository;
$\theta'=0.56$ throughout, $N$ ranging over $2.5\cdot10^4$ to
$8\cdot10^5$ (the ranges are small because the residual is computed
exactly, by enumeration).


- **Switching identity** (`audit\_switch\_identity.py`).
Both [eq:switch] and the split $D(t)=P(t)-R(t)$ of [eq:PR] hold to
$10^{-16}$ relative to $N$, at $N=2.5\cdot10^4$ through $4\cdot10^5$
and five truncations $t$ each. Lemma [lem:complete] is checked not
numerically but as an integer identity, on *every* squarefree
$u<N$: zero mismatches. The support of $P$ is exactly the set of
divisors of $\rad(N)$ below $N$, which is what [eq:P] asserts, and
$\max_t|P(t)|$ stays an order of magnitude under $2^{\omega(N)}\log N$.


- **The residual is its main term**
(`audit\_residual\_mainterm.py`). With the main terms restricted to
$(m,N)=1$ as required by Remark [rem:trap], and
$c(m)=A(N)\lambda(m)/m$ of Lemma [lem:density]:

$$
\begin{array}{r|ccccc}
 N & 5\cdot10^4 & 10^5 & 2\cdot10^5 & 4\cdot10^5 & 8\cdot10^5\\\hline
 R/N               & 0.1140 & 0.0965 & 0.0785 & 0.0618 & 0.0483\\
 \mathrm{MT}/N     & 0.1190 & 0.0964 & 0.0770 & 0.0616 & 0.0491\\
 |R-\mathrm{MT}|/R & 0.0437 & 0.0017 & 0.0191 & 0.0033 & 0.0160
\end{array}
$$

The residual *is* the predicted main term to between $0.2\%$ and
$4.4\%$, and that main term decays through the cancellation of
Lemma [lem:mu]. The observed decay is
$R\asymp N^{1-\theta'/2}$ (the ratio $R/N^{1-\theta'/2}$ reads
$2.358,\,2.425,\,2.394,\,2.290,\,2.172$, inside $[2.17,2.46]$ over the
whole range), not $\asymp N$; a least-squares fit of $\log R$ against
$\log N$ gives the exponent $0.6880$, against $1-\theta'/2=0.72$.

The $(m,N)=1$ restriction is not a refinement. Dropping it and
applying the same $c(m)$ to every $m<M$ takes $\mathrm{MT}/N$ to
$-0.0460,\,-0.0357,\,-0.0293,\,-0.0210,\,-0.0157$ — the wrong sign,
and $134\%$ to $140\%$ away from $R$. The local factor at $p\mid N$ is
$1$ and not $p^{-1}\lambda(p)$, so applying the coprime formula there
is not a small error in a constant.


#### Remark (a pre-registered check of ours that failed) {#rem:band}
<!-- evidence: audit_residual_mainterm.py -->

Version 3 of this note summarised the table above as
"$1$–$4\%$ agreement". The audit pre-registered that band as its
rule U3 and it fails: the worst entry is $0.0437$ at
$N=5\cdot10^4$, so the band is $0.2\%$ to $4.4\%$ and the script exits
non-zero on that account. Both printed rows themselves reproduce to
the last digit shown — the recomputed $R/N$ and $\mathrm{MT}/N$
differ from the printed values by at most $4.8\cdot10^{-5}$ and
$2.4\cdot10^{-5}$ — so what failed is the sentence describing the
table and not the table. It is recorded because a stated band is a
falsifiable claim and this one was false, and because a summary that
rounds its own worst case downward is the failure that a reader cannot
detect without recomputing.


- **The $\log k$ branch does not decay**
(`audit\_logweight\_branch.py`). The same quantity with $w_k=\log k$
has $R_{\log}/N=2.2842,\,2.2625,\,2.2068,\,2.1333,\,2.0669$ against the
predicted $2.3257,\,2.2589,\,2.1902,\,2.1302,\,2.0754$, agreeing to
between $0.14\%$ and $1.81\%$ — an object of size $\asymp N$,
consistent with Theorem [thm:C]. The contrast with the bullet above is
the whole point: on the same field and the same five $N$, the $w_k=1$
residual falls from $0.1140N$ to $0.0483N$ while this one stays above
$2N$, a ratio of $42.8$ at $N=8\cdot10^5$. Both rows fall towards
$\SS(N)=1.76043$ from above, which is what
Theorem [thm:C] requires and what fixes the constant.


- **The constants** (`audit\_E3\_constant.py`). $B_{\log}(K)$
approaches $-\SS(N)$, and $A(N)\widetilde G(x)=-1.760250$ at
$x=4\cdot10^6$ against $-\SS(N)=-1.760432$, as required by (i) and (ii)
of Section [sec:C]. **Both constants are negative**; an earlier
version of this bullet printed the second as $+\SS(N)$, which is the
error Remark [rem:sign] records.


- **The load-bearing identity**
(`audit\_density\_identity.py`). $\sum_{g\mid m}
\mu(g)/(\varphi(m/g)g\varphi(g))=1/m$ verified in exact rational
arithmetic for all $243$ squarefree $m<400$: zero mismatches. See
Remark [rem:loadbearing].


## Summary


- Theorem [thm:A]: the Möbius-weighted, fixed-class
correlation sum with weight $w_k=1$ is $\ll_A N(\log N)^{-A}$,
unconditionally, for any level $N^{\theta'}$ with $\theta'>1/2$.

- Corollary [cor:B]: hence Huang–Li's $E_4$ (their Lemma 4)
needs no hypothesis, and the whole $EH_\mu$ demand of their reduction
collapses to the scalar $E_3$.

- Theorem [thm:C]: that scalar is unconditionally equivalent to
Huang–Li's equation (22), hence already gives binary Goldbach for
large even $N$, and gives the asymptotic $\tilde r(N)\sim\SS(N)N$
exactly when $C(N)=o(N)$.

- Proposition [prop:onesided]: Goldbach itself needs only the
one-sided $E_3 > -\SS(N)(1-A(N))N(1+o(1))$, whose threshold is
$symp N$ for almost all even $N$ and never below
$cN/(\log N\log\log N)$ — weaker than the consumed bound by a factor
$(\log N)^A$ for every $A$. It does not reopen Theorem [thm:D].

- Section [sec:delta]: equation (18) of [HL] omits an
$n$-dependent constraint — S. Zheleznov's observation; the missing term
$\Delta$ is exhibited and closes under hypotheses already assumed. The
authors' correction keeps the constraint instead, and
[rem:movingswitch] shows the switch of Theorem [thm:A] survives that
move, reaching the favourable configuration one step earlier.

- Proposition [prop:movingcut]: the bound of Theorem [thm:A] holds
against the corrected formulation, unconditionally for $\theta'>1/2$.
The two formulations differ by a term of $\Delta$'s shape
([eq:cutbridge]), so the fixed cut is not superseded: it is one of the
two ingredients and the omitted term is the other.

- Section [sec:notclaimed]: the net progress toward Goldbach is
zero.


## Relation to the literature {#sec:lit}


This section is a placement, not a survey, and it is deliberately
conservative about what is offered here as new.


**Theorem [thm:A].** The mechanism is standard. After the
divisor switch the Möbius on the long variable squares away and the
surviving Möbius sits on a variable of length at most
$N^{1/2-\delta}$, which is the classical configuration in which
Bombieri–Vinogradov applies. Every ingredient
(Bombieri–Vinogradov; the Goldston–Yıldırım estimate
[GY] as used in [HL]; the elementary density identity of
Proposition [prop:MT]) is off the shelf. The statement is offered as
a *lemma about the Huang–Li reduction*, not as a contribution to
the theory of Bombieri–Vinogradov; its content is that this particular
consumption of $EH_\mu$ is unnecessary, and what difficulty it has lies
in the bookkeeping — in particular that the condition $(k,N)=1$ must
be discarded by the degeneracy argument rather than expanded over
$e\mid N$, and that main terms may be assigned only to classes with
$(q,N)=1$, the second of which shifts the density by exactly
$N/\varphi(N)$ if missed.


**Theorem [thm:D.** and Proposition \ref{prop:E].} See
Remark [rem:bombieri] for the relation to [Bom76].
Proposition [prop:E] is an observation — Parseval bounds any
pointwise route from below by the $L^2$ norm — and is recorded because
the circle method is the first thing one tries on $C(N)$, not because it
is difficult.


**Beyond the square-root barrier..** Since [HL] appeared,
Lichtman [Li23] has shown that the primes have level of
distribution $66/107 \approx 0.617$ using triply well-factorable
weights, and has used it to improve upper bounds for Goldbach
representations — the first use of a level beyond the square-root
barrier on that problem. That work does not collide with the present
one: it concerns the primes rather than
$\sum_{n<N}\Lambda(n)\mu(N-n)$, and it yields upper bounds rather than
the asymptotic that [HL] targets. It is nevertheless the state of
the art against which any claim of the form "beyond $1/2$" must now be
measured.

A level of distribution $3/5$ for the Möbius function with triply
well-factorable weights belongs to the earlier part of that line of
work, [Li20]; we have not verified its precise statement or
location against the published version, and cite it here as the object
the question below refers to rather than as a result we have checked.
It is *not* in [Li23], which does not treat the Möbius
function at all — the $3/5$ appearing there is Maynard's prior record
for the *primes*, which [Li23] improves to $66/107$. Since
[HL] require $EH_\mu(N^{\theta'})$ for a single $\theta'>1/2$, and
since Theorem [thm:C] identifies $w_k=\log k$ as the weight carrying
the Goldbach content, the following question is well posed and, so far
as the present authors know, open:

\begin{quote}
Does the weight that the Huang–Li consumption actually requires lie in
the triply well-factorable class for which level $3/5$ is known for
$\mu$?
\end{quote}

Two cautions attach to it. First, $EH_\mu$ as stated in [eq:EHmu]
carries $\max_a|\cdot|$, whereas well-factorable results bound signed
weighted sums $\sum_q \lambda_q E(x;q,a)$; the comparison must therefore
be made against what the argument *consumes*, which is a signed
sum, and not against the absolute-value form as written. Second,
Theorem [thm:D] closes extraction *by divisor switching*; the
technology of [Li23] is the Deshouillers–Iwaniec spectral large
sieve, a different mechanism, and Theorem [thm:D] says nothing about
it.


## References

- **[HL]**  Jing-Jing Huang and Huixi Li,
*On the connection between the Goldbach conjecture and the
Elliott–Halberstam conjecture*, arXiv:2005.03811v2 [math.NT], 2022.

- **[GY]**  Daniel Goldston and Cem Y{ı}ld{ı}r{ı}m,
*Higher correlations of divisor sums related to primes I*, Integers
**3** (2003), A5.

- **[Pan]**  Cheng-Dong Pan, *A new attempt on Goldbach
conjecture*, Chinese Ann. Math. **3** (1982), 555–560.

- **[MV]**  Hugh Montgomery and Robert Vaughan,
*Multiplicative number theory I: classical theory*, CUP, 2007.

- **[Bom76]**  Enrico Bombieri, *The asymptotic sieve*,
Rend. Accad. Naz. XL (5) **1/2** (1975/76), 243–269.

- **[Li23]**  Jared Duker Lichtman, *Primes in arithmetic
progressions to large moduli, and Goldbach beyond the square-root
barrier*, arXiv:2309.08522 [math.NT], 2023.

- **[Li20]**  Jared Duker Lichtman, *Primes in arithmetic
progressions to large moduli, II: well-factorable estimates*,
arXiv:2006.07088 [math.NT], 2020.

- **[Vau88]**  R. C. Vaughan, *The $L^1$ mean of exponential
sums over primes*, Bull. London Math. Soc. **20** (1988),
121–123.
