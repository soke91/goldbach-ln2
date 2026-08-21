<!-- v3. Restructured around Theorem A; v3/README.md records what
     changed and why. Proof sections are carried over from
     v2/paper/theorem_A.md unaltered. -->

# The flat branch of $EH_\mu$ is unconditional, and what follows for the Huang–Li reduction of Goldbach

```latex
% 수식이 쓰는 매크로 — 렌더러/역변환용
\renewcommand{\SS}{\mathfrak{S}}
\newcommand{\Emu}{E_\mu}
\DeclareMathOperator{\rad}{rad}
\DeclareMathOperator{\lcm}{lcm}
```


## Abstract

Murty and Vatwani **[MV17]** posed an Elliott–Halberstam conjecture for
the Möbius function on shifted primes and showed that, with $EH$, it
breaks the parity barrier for twin primes. Huang and Li **[HL]**
transposed the hypothesis to the Goldbach shift $n\mapsto N-n$ and
showed that binary Goldbach for large even $N$ follows from $EH$
together with $EH_\mu$ when the levels sum to more than $1$; by
Bombieri–Vinogradov their Corollary 1 reduces this to
$EH_\mu(N^{\theta'})$ alone for a single $\theta'>1/2$. In both papers
the hypothesis is consumed at exactly two places, distinguished by the
weight $w_k$ that the divisor decomposition attaches: $w_k=\log k$,
which is Huang–Li's $E_3(\alpha)$, and $w_k=1$, which is their
$E_4(\alpha)$.

**This note is about the second.** Theorem [thm:A] proves
unconditionally that the Möbius-weighted correlation sum in the fixed
class $n\equiv N\ (k)$ with $w_k=1$ is $\ll_A N(\log N)^{-A}$ for every
$A>0$, uniformly in the truncation point. Both **[MV17]** and **[HL]**
spend hypothesis on that branch — the first by partial summation on
$EH_{\mu_h}$, the second through its Lemma 4 — and Corollary [cor:B]
is that neither needs to: the $E_4$ consumption is unnecessary, and the
whole $EH_\mu$ demand of the reduction collapses to the single scalar
$E_3(\alpha)$.

What that leaves is then identified exactly. Corollary [thm:C] states
the identity behind the collapse: $E_3(\alpha)$ *is* the difference
between the prime-pair count and $\SS(N)\bigl(N-C(N)\bigr)$, up to
$O_A(N(\log N)^{-A})$, where $C(N)=\sum_{n<N}\Lambda(n)\mu(N-n)$. The
chain is **[MV17]**'s and its endpoint is **[HL]**'s equation (22);
what is added is that the residual is named rather than absorbed into
an $o(1)$, and that the statement is unconditional — which is possible
only because Theorem [thm:A] has discharged the other branch. The root
cause is $\mu*\log=\Lambda$. Propositions [prop:onesided] and
[prop:nolog] then separate closure at the level of identities from
closure at the level of strength: Goldbach needs the bound on one side
only, and at a constant factor rather than a saving of every power of
$\log N$.

Finally, Section [sec:delta] concerns a defect in the published form of
**[HL]**'s equation (18), which drops an $n$-dependent constraint. The
defect is not claimed here: it was reported to the authors, before this
note, by S. Zheleznov, and they have corrected the manuscript by
keeping the constraint so that the dropped range never arises. The
missing term $\Delta$ and its treatment are what this section sets
down, and Proposition [prop:movingcut] is what they are for: the
published form and the corrected one differ by exactly a term of
$\Delta$'s shape, so Theorem [thm:A], proved against the published cut,
holds against the corrected one as well.

**What is not claimed.** Theorem [thm:A] removes the part of the demand
that carries no Goldbach content, and the net progress toward the
Goldbach conjecture is zero. The ingredients are classical and the
mechanism — a divisor switch that moves the work onto a short variable,
followed by Bombieri–Vinogradov — is ordinary practice; what is offered
is the application. Section [sec:lit] places the whole against the
record, including the companion no-go over the design space of weights,
which is not reproduced here.

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

Throughout, $\theta'\in(1/2,1)$ is fixed, $K=N^{\theta'}$ and
$\alpha=N^{1-\theta'}$. Theorem [thm:A] is the result; Corollaries
[cor:B] and [thm:C] are what it buys, and the two propositions after
them say how much of the remaining demand Goldbach actually needs. The
proof of Theorem [thm:A] occupies Section [sec:mechanism] onward and is
self-contained.

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


Theorem [thm:C] says that the demand side is closed at the level of
*identities* rather than of estimates: no choice of $\theta'$, of
truncation, or of smoothing can evade it, because the root cause is the
identity $\mu*\log=\Lambda$ from which Huang–Li start. The weight
$\log k$ is the carrier of the Goldbach content, so every divisor
switch must hand it back.


#### Remark (whose chain this is) {#rem:thmCprior}
<!-- evidence: analytic -->

The chain behind Theorem [thm:C] is not ours, and this note was
carrying it without attribution until the literature was checked.
Murty and Vatwani **[MV17]** run it for the shift $n\mapsto n+h$:
$\Lambda=\mu*\log$, split at $y$, the factor $\mu^2(n+h)$ carried
throughout, and a reduction to sums they call $S_3$ and $S_4$ — where
$S_3$ is exactly the $\log$-weighted Möbius correlation sum in the
fixed class, which is $E_3$. Granting $EH_\Lambda$ and their
shifted-Möbius $EH_{\mu_h}$ they reach

$$
\sum_{n\le x}\Lambda(n)\Lambda(n+h)\;\sim\;
  \bigl(\SS(h)+o(1)\bigr)\Bigl(x-\sum_{n\le x}\Lambda(n)\mu(n+h)\Bigr),
$$

which is the right-hand side of Theorem [thm:C] with $E_3$ absorbed
into the $o(1)$; and their $A_h=\prod_{p\nmid h,\,p>2}(1-1/(p(p-1)))$
is [eq:AN]. [HL] transpose the chain to the Goldbach shift and their
(22) is its endpoint. Their own reference list carries **[MV17]**,
Vatwani's **[Vat19]**, and Pan's Goldbach-side attempt **[Pan]**;
this note carried **[MV17]** in its bibliography and cited it nowhere.

Three things are added here and the chain is not among them. The
residual is *named* instead of absorbed, so the statement is an
equality with an explicit error term rather than an asymptotic. It is
*unconditional* in the Corollary-1 regime — possible only because
Theorem [thm:A] discharges the other consumption of $EH_\mu$, which
**[MV17]** spend on $S_4$ by partial summation and [HL] spend on $E_4$
through their Lemma 4. And the three-way distinction of
[rem:threeway] is visible only once the $o(1)$ is a named object.

Theorem [thm:A] itself has no counterpart in either: both papers assume
their hypothesis on the flat branch, and the claim here is that the
branch does not need one. Nothing found in **[MV17]**, **[Vat19]**,
**[Pan]** or [HL]'s reference list bears on it. That is a negative
search over four items and not a survey.

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

What it does not open is the divisor switch. Theorem the companion no-go costs
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
  B(N) \;:=\; \sum_{\substack{k<K\\ (k,N)=1}} \mu^2(k)\,(\log k)\,
  \bigl|\Emu(N;k)\bigr| ,
\end{equation}

the restriction to squarefree $k$ being what the triangle inequality
on [eq:E3] leaves, since $\mu(k)=0$ elsewhere.
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
measured, which is what Remark the withdrawn level measurement recorded in v2 does.

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
factors, $k$ must avoid them, and by the dilate identity of the companion note so must $m$: for
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



#### Remark (the sign is locked, and the reason on record was the wrong one) {#rem:signlock}
<!-- evidence: audit_signlock.py -->

[rem:threshfam] records that at $N=2\cdot3^2\cdot5\cdot7\cdot11\cdot13
\cdot17$ all $513$ terms of $E_3$ are negative and $|E_3|=B(N)$
exactly; the arithmetic-family sweep recorded in v2 records
$\lvert\sum a\rvert/\ell^1=1.00000000$ on a primorial-like family at
every swept $\theta'$. Two measurements, taken for different reasons,
found the same thing, and neither asked whether it survives $N$
growing.

**It survives, across the range reachable.** On the family
$N=30030\,j$, $j=1,2,4,\dots,128$ — radical containing
$2,3,5,7,11,13$ at every point — the fraction of $k$ with $H(N;k)<0$ is
$1.000000$ at all eight $N$, over $60$ to $921$ nonzero terms, and
$\rho=\lvert\sum_k(\log k)H\rvert/\sum_k\lvert(\log k)H\rvert$ is
$1.0000000000$. The swept family $2^a5^b$ run through identical code
cancels, and cancels *more* as $N$ grows: $\rho=0.545340$, $0.554436$,
$0.453017$, $0.386453$, $0.358519$. The two families diverge rather
than converge.

**The reason on record does not survive.** [rem:threshfam] explains the
lock by saying the admissible $m$ "are almost all primes, so
$\mu(m)=-1$ dominates". The prime fraction among the $m$ that carry a
term falls monotonically — $0.898907$, $0.856926$, $0.809884$,
$0.767039$, $0.721740$, $0.683812$, $0.650484$, $0.617012$ — while the
lock does not weaken at all. At the largest $N$ nearly two in five
contributing $m$ are not prime and not one of $921$ sums changes sign.
The registered rule S3 asked whether the prime fraction falls below
$1/2$ while the lock holds; it does not, in this range, so **S3 holds
as registered and the prediction attached to it was wrong**. What is
wrong with the explanation is not its verdict but its content.

**The mechanism is forced, not statistical.** Let $q\mid\rad(N)$. If
$q\mid m$ then $q\mid N-mk$, so $N-mk$ is a prime power only if it is a
power of $q$ itself. Primality of $N-mk$ therefore *forces* $m$ to be
coprime to $\rad(N)$, with only the degenerate exceptions: at
$N=120120$, $31$ of the $9389$ contributing $m$ fail coprimality and
all of them are that case. So the inner sum is not over all $m$ but
over the $\rad(N)$-rough ones, and among those the numbers with one
prime factor outnumber those with two — both factors being at least
$17$ here, hence at least $289$ in product — by a wide margin at
accessible sizes. That is why $\mu(m)=-1$ wins term by term.

**And the mechanism was claimed here to have an expiry. It does not** —
see [rem:signlockmargin], which measured the margin and corrects this
paragraph. The argument made here was that the prime share of the
$\rad(N)$-rough numbers below $x$ is
$\sim(\log x)^{-1}\prod_{p\mid\rad(N)}(1-1/p)^{-1}$ and tends to zero,
so the bias producing the lock thins like $1/\log$. The share does
thin; the sign does not follow it. What decides the sign is the parity
of $\omega$, and the two are different quantities. The paragraph is
kept as written because the next remark is a correction of it.

**What this is not.** It is not a statement about $C(N)$, and it bears
on the Goldbach problem in one narrow way only: [prop:onesided]'s
threshold is $\SS(N)(1-\AAA(N))N$, and $1-\AAA(N)$ collapses on exactly
the family where the lock holds. Where the threshold is thinnest, the
sum it constrains has no cancellation to spend. That is a remark about
which of [prop:onesided]'s two inequalities is doing the work, not a
bound on anything.



#### Remark (the margin saturates, and the lock is about rough integers) {#rem:signlockmargin}
<!-- evidence: audit_signlock_margin.py -->
<!-- evidence: audit_signlock_why.py -->

[rem:signlock] closed by saying the mechanism "has an expiry" — that
the prime share among rough numbers thins like $1/\log$, so the lock
should end. **That reasoning is wrong, and the correction is the point
of this remark.** What thins like $1/\log$ is the share of the
contributing $m$ that are prime. What decides the sign is the parity of
$\omega$, and those are not the same quantity.

**The margin saturates rather than crosses.** Split $H(N;k)$ by the
sign of $\mu(m)$ and write $r_k$ for the positive $\Lambda$-mass over
the negative one, so that the lock is $\max_k r_k<1$. On
$N=30030\cdot2^{\,j}$, $j=0,\dots,9$,

$$
\max_k r_k = 0.191113,\ 0.289221,\ 0.428727,\ 0.523536,\ 0.662493,\
0.772145,\ 0.809604,\ 0.871314,\ 0.914678,\ 0.938479 ,
$$

over $60$ to $1999$ terms. The increments run
$0.098,\,0.140,\,0.095,\,0.139,\,0.110,\,0.037,\,0.062,\,0.043,\,0.024$:
the rise decelerates by a factor of five across the range. A straight
line in $\log N$ fitted to these ten points reaches $1$ at $N=10^{7.08}$
— and the lock still holds at $N=1.5\cdot10^7$, so that model is
refuted by the data it was fitted to. It is printed in the result file
for exactly that reason. The control family $2^a5^b$, run through the
same code, is not locked at any $N$: $\max_k r_k$ reads $2.279886$,
$2.406981$, $2.342519$, $2.698008$, $1.935344$, $1.882319$.

**Why it saturates.** Landau's count of squarefree $m\le x$ free of
prime factors below $z$ with $\omega(m)=j$ is
$\sim(x/\log x)\lambda^{j-1}/(j-1)!$ with
$\lambda=\log\log x-\log\log z$. The index starts at $j=1$, so the odd
$j$ — where $\mu=-1$ — carry $\cosh\lambda$ and the even $j$ carry
$\sinh\lambda$, and the ratio is $\tanh\lambda<1$ for every $\lambda$.
The bias is not a shortage of composites; it is that the first term of
the distribution is a prime. As $\lambda\to\infty$ the ratio approaches
$1$ and never reaches it, which is what the ten points do. Since
$\lambda$ moves like $\log\log$, the approach is as slow as the table
shows. This is a heuristic, not a proof: Landau's asymptotic is poor at
the sizes where $x=N/k$ is a few hundred, which is most of the range.

**And the primes contribute the roughness and little else.** Primality
of $N-mk$ enters twice — it forces $m$ coprime to $\rad(N)$, and it
selects which rough $m$ carry a term at all. Discarding the second and
counting the same index set by parity alone gives $q_k$, and $q_k$
tracks $r_k$ closely: across $k$ the correlation is $0.9135$, $0.9596$,
$0.9779$, $0.9785$, $0.9804$, $0.9781$, $0.9746$, $0.9720$ at
$N=60060$ upward, and the medians track it -- $0.282392$ against $0.257840$ at the top. The
unweighted ratio is locked too, at every $N$ ($\max_k q_k=0.899701$
against $\max_k r_k=0.914682$ at the top). **So the lock is, to within
a few per cent, the parity of $\omega$ on rough integers below $N/k$**,
and the selection by primality nudges it upward without unlocking it.

The registered rule W2 asked for correlation $\ge0.8$ at *every* $N$
and the median gap $\le0.15$. At the smallest point, $N=30030$ with
only $60$ terms, the correlation is $0.7423$. **W2 is refuted as
registered.** The refutation is confined to that point and the
threshold was too tight for a correlation over $60$ values; the
conclusion it was testing is supported at the other eight. That is
stated rather than repaired.

**What this settles.** The two earlier measurements found something
real and it is not about the Goldbach problem. It is an elementary
fact about short rough ranges, visible in $E_3$ because
$\rad(N)$-roughness is forced there. The one thing it leaves standing
is the coincidence of position recorded in [rem:signlock]:
[prop:onesided]'s threshold $\SS(N)(1-\AAA(N))N$ collapses on exactly
the family where the terms do not cancel. That remains a statement
about which of that proposition's two inequalities does the work, and
it is not a route to anything.


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


The published form of equation (18) drops a range. Two separate things
follow and they carry different credit. **The defect is not claimed as
new**: that the $n$-dependent restriction had been replaced by a fixed
outer range was reported to the authors by S. Zheleznov before this
note was written, and they have since prepared a corrected manuscript
— see [rem:movingswitch]. **The missing term and its treatment are what
this section sets down**, and nothing on record attributes them
elsewhere. The two are different kinds of work and neither replaces the
other. A correction that keeps the moving truncation repairs the
derivation *from the inside*: the range $m\le\alpha$ never appears, so
there is no term to bound. What follows is the *outside* version — the
published form is taken as it stands, the discrepancy is named, and it
is bounded. The second is not the first done again: the term **is** the
difference between the two formulations, and that is what lets a
theorem proved against the fixed cut transfer to the corrected one
([prop:movingcut]).

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

The omission recorded above had been reported to the authors by
S. Zheleznov before this note, and they have prepared a corrected
version of the
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

with $\rho_N$ as in Lemma the extraction lemma of the companion no-go. Splitting at $N-n=H$ with
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
with $\rho_N$ as in Lemma the extraction lemma of the companion no-go. Its hypothesis is met at both
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
threshold the companion no-go turns on, and it is not a coincidence: below it the
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

- **Theorem [thm:A].** The Möbius-weighted, fixed-class correlation sum
with weight $w_k=1$ is $\ll_A N(\log N)^{-A}$, unconditionally, for any
level $N^{\theta'}$ with $\theta'>1/2$, uniformly in the truncation
point. This is the result; everything below follows from it or bounds
what it leaves.

- **Corollary [cor:B].** Hence the $E_4$ consumption of $EH_\mu$ —
Huang–Li's Lemma 4, and the partial-summation step **[MV17]** makes at
the same place — is unnecessary, and the whole $EH_\mu$ demand of the
reduction collapses to the scalar $E_3$.

- **Corollary [thm:C].** That scalar is unconditionally equal to
$\tilde r(N)-\SS(N)(N-C(N))$ up to $O_A(N(\log N)^{-A})$; hence
$E_3\ll_A N(\log N)^{-A}$ is equivalent to **[HL]**'s equation (22),
already gives binary Goldbach for large even $N$, and gives the
asymptotic $\tilde r(N)\sim\SS(N)N$ exactly when $C(N)=o(N)$. The chain
is **[MV17]**'s; the unconditionality is Theorem [thm:A]'s.

- **Proposition [prop:onesided].** Goldbach needs only the one-sided
$E_3>-\SS(N)(1-\AAA(N))N(1+o(1))$, whose threshold is $\asymp N$ for
almost all even $N$ and never below $cN/(\log N\log\log N)$ — weaker
than the consumed bound by a factor $(\log N)^A$ for every $A$. Note
[rem:threshfam] records that the threshold is not a constant and that
the family every sweep here uses sits at its maximum.

- **Proposition [prop:nolog].** The demand needs no saving in $\log N$
at all, only a constant factor.

- **Section [sec:delta] and Proposition [prop:movingcut].** The
published form of **[HL]**'s equation (18) drops an $n$-dependent
constraint — reported to the authors independently, and before this
note, by S. Zheleznov. The missing term $\Delta$ is exhibited here and
closes under hypotheses already assumed. The authors' correction keeps
the constraint instead, so $\Delta$ does not arise there; the two
formulations differ by exactly a term of that shape
([eq:cutbridge]), and Theorem [thm:A] together with it gives the bound
against the corrected formulation.

- **What is not claimed.** Net progress toward Goldbach: zero.

## Relation to the literature {#sec:lit}

**The hypothesis is Murty–Vatwani's.** **[MV17]** posed the shifted
Möbius Elliott–Halberstam conjecture $EH_{\mu_h}$ and proved that, with
$EH_\Lambda$, it makes $\sum_{n\le x}\Lambda(n)\Lambda(n+h)\sim\SS(h)x$
and $\sum_{n\le x}\Lambda(n)\mu(n+h)=o(x)$ equivalent, and that the
twin prime conjecture follows. **[HL]** transposed the hypothesis and
the argument to the Goldbach shift. **[Vat19]** develops the
equidistribution side; the Goldbach-side precursor is Pan **[Pan]**,
which this note has not been able to consult and does not rely on.
Everything in Section [sec:mechanism] onward is written against
**[HL]**'s formulation because that is the one stated for $N-n$, but
the decomposition, the $\mu^2$ restriction and the split at $y$ are
**[MV17]**'s.

**What Theorem [thm:A] adds.** Both papers assume their hypothesis on
the flat branch. **[MV17]** obtain equidistribution for
$\Lambda(n)\mu(n+h)\log(n+h)$ from $EH_{\mu_h}$ by partial summation
(their Proposition 4.2); **[HL]** use their Lemma 4 at the same place.
Theorem [thm:A] is the assertion that the branch carries no hypothesis
at all. A search of **[MV17]**, **[Vat19]**, **[HL]**'s reference list
and the surrounding literature found nothing that states it; that is a
negative search over a short list and not a survey, and **[Pan]** in
particular is unchecked.

**The mechanism is standard.** A divisor switch that moves the work
onto a short variable, followed by Bombieri–Vinogradov, is ordinary
practice; the completion of the divisor sum is the step that makes the
cofactor short, and Section [sec:mechanism] says why it is not a
formality. Lemma [lem:mu] is the classical zero-free-region estimate
for $\sum\mu(a)/a$, and Lemma the extraction lemma of the companion no-go's bound on $\rho$ is
Huang–Li's Lemma 1, a form of Goldston–Yıldırım **[GY]**, carrying its
hypothesis $\log n\ll\log x$.

**The companion no-go.** That $E_3$ cannot be reached from $E_4$ by
reweighting — that no weight between $w_k=1$ and $w_k=\log k$ extracts
$C(N)$ by divisor switching and Bombieri–Vinogradov, even granting
$EH(N^{\theta_E})$ for every fixed $\theta_E<1$ — is proved separately
and is not reproduced here. In genre it is a precise form, for this
reduction, of the phenomenon Bombieri's asymptotic sieve **[Bom76]**
records, and a reader who finds it unsurprising is right to.

**Beyond the square-root barrier.** Lichtman **[Li23]** obtains a level
of distribution past $1/2$ for the primes, which postdates **[HL]**;
whether an analogue for $\mu$ on shifted primes would move
Corollary [thm:C]'s regime is open and nothing here bears on it.

## References

- **[HL]**  Jing-Jing Huang and Huixi Li,
*On the connection between the Goldbach conjecture and the
Elliott–Halberstam conjecture*, arXiv:2005.03811v2 [math.NT], 2022.

- **[GY]**  Daniel Goldston and Cem Y{ı}ld{ı}r{ı}m,
*Higher correlations of divisor sums related to primes I*, Integers
**3** (2003), A5.

- **[Pan]**  Cheng-Dong Pan, *A new attempt on Goldbach
conjecture*, Chinese Ann. Math. **3** (1982), 555–560.

- **[MV17]**  M. Ram Murty and Akshaa Vatwani, *Twin primes and
the parity problem*, J. Number Theory **180** (2017), 643-659.

- **[Vat19]**  Akshaa Vatwani, *Variants of equidistribution in
arithmetic progression and the twin prime conjecture*, Math. Z.
**293** (2019), 285-317.

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
