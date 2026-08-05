# Contact drafts (NOT sent — for the user to review, edit, and send)

*These are drafts only. Sending is the user's decision and action.
Both are framed as questions, lead with measurements, and disclose the
AI-assisted provenance.*

---

## Draft A — email to Prof. J. D. Lichtman

Subject: Numerical landscape of the EH_mu object in your Goldbach
framework — is this known?

Dear Professor Lichtman,

I am writing because a computational program I have been running lands
squarely on two of your papers, and I would value your judgment on
whether its measurements are known or of any interest.

Following Huang–Li (arXiv:2005.03811), Goldbach for large even N
reduces to the Mobius-twisted hypothesis EH_mu past the square-root
barrier, whose object is c(n) = Lambda(n) mu(N−n). We measured that
object's landscape directly: the fixed-residue discrepancy of c(n)
against the random-walk benchmark stays at 0.6–1.1 straight through
theta = 0.30 → 0.70 (five values of N, scales 1e8 and 1e9), with no
visible change at the barrier; the innermost sum Sum Lambda(n)mu(N−n)
follows a textbook half-normal law over 340 values of N; and in the
thin-progression regime (moduli beyond the square root of the range,
where as far as I can tell even GRH-strength bounds are void), 10,000
sampled progressions give Mobius sums at 0.74 of the half-normal
benchmark with an exactly Gaussian tail.

Everything regenerates from ~50 short numpy scripts in a public
repository (github.com/soke91/goldbach-ln2 — see MEASUREMENTS.md, with
one-shot reproduction via verify_all.py). The repository also contains
a proof-program sketch that attempted to route the k-face of the
Vaughan decomposition through the dispersion machinery of your
2309.08522; an internal adversarial review has since REFUTED its core
reductions (the verdict is published in the repo as REVIEW_VERDICT.md
— in short: our substitution violated Lemma 5.1's Q-premise, and our
pairs lack the congruence that powers your conductor collapse), so
please disregard it entirely — my questions are only about the
measurements:

1. Are direct numerics of the EH_mu object's discrepancy landscape
   (through the barrier, fixed residue) known in the literature?
2. Is the thin-progression behaviour (sub-half-normal at L > sqrt(y))
   expected/known folklore, or worth documenting?
3. Is there a known obstruction to the general shape "dilate-averaged
   Mobius over shifted primes" (the k-averaged analog of your
   2009.08969) that would make that direction hopeless?

Full disclosure: this program was AI-assisted (documented in the
repository, including 24 corrections and every failed route). I am an
amateur; the measurements are honest and reproducible, and I would be
grateful for even a one-line verdict.

With respect and thanks for your two papers, which this program kept
rediscovering as the frontier,
[user's name]

---

## Draft B — MathOverflow question

Title: Numerics for the Mobius-twisted Elliott–Halberstam object
through the square-root barrier — known?

Body:

By Huang–Li (arXiv:2005.03811), binary Goldbach for large even N
follows from BV plus EH_mu past x^{1/2} (fixed residue), whose object
is c(n) = Lambda(n) mu(N−n). I have direct numerics of this object's
landscape and could not find comparable computations in the
literature:

- fixed-residue discrepancy of c(n) vs the random-walk benchmark:
  mean ratio 0.6–1.1 across moduli x^{0.30} to x^{0.70} (through the
  barrier), five N's at scales 1e8/1e9, no behaviour change at 1/2;
- Sum_{n<N} Lambda(n) mu(N−n): half-normal over 340 values of N;
- thin-progression Mobius sums (moduli L > sqrt(range)): mean 0.74 of
  half-normal, Gaussian tail, 10,000 classes.

Everything is reproducible from short scripts (repo linked; one-shot
verify script included; AI-assisted provenance disclosed). The repo
also contains a proof-program sketch whose core reductions were
refuted by an internal adversarial review (verdict published in the
repo); this question concerns only the measurements.

Questions: (1) Are such numerics known/published? (2) Is the
thin-progression sub-half-normal behaviour folklore? (3) Is there a
known structural obstruction to dilate-averaged analogs of averaged
Mobius-on-shifted-primes results (Lichtman arXiv:2009.08969) that
would collapse this direction?

(I am aware numerics prove nothing; the question is about the
literature and known obstructions, not a claim.)

---

*Recommended order: send Draft A first; post Draft B only if no reply
in ~2 weeks.*

*Update (increment 143): the adversarial internal review returned its
verdict — the sketch's core reductions are refuted (REVIEW_VERDICT.md).
Both drafts have been updated to disclose this. The measurements are
unaffected (the review explicitly affirmed the numerics), so the
question-framed contact remains appropriate; the drafts now carry zero
implied claim about the sketch.*
