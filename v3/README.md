# v3 — the flat branch

One document, one result, and the evidence it cites.

```
paper/theorem_A.md   the note
code/                the scripts it cites
results/             their output, one file each
gate/gate.py         the checks; 0 is the condition for a commit
verify/              re-verification passes, when they run
```

The result is that the Möbius-weighted correlation sum in the fixed
class $n\equiv N\ (k)$, taken with the flat weight $w_k=1$, is
$\ll_A N(\log N)^{-A}$ unconditionally and uniformly in the truncation
point. Murty and Vatwani posed the hypothesis this note is about, for
the shift $n\mapsto n+h$; Huang and Li transposed it to $n\mapsto N-n$.
Both spend the hypothesis on the flat branch. The claim here is that
the branch does not need it.

---

## What changed from v2

v2 held three documents' worth of material in one file: the spine
around this theorem, the exact identities of the divisor-switch design
space, and the no-go built on them. The literature check that closed
the question of what in it was new made the ordering wrong. The chain
behind the identity for $E_3$ is Murty–Vatwani's, run on the Goldbach
shift by Huang–Li; what is not theirs is the flat branch. So the
headline is the flat branch, and the identity is a corollary of it
together with their chain.

| | v2 | v3 |
|---|---|---|
| headline | three results side by side | one theorem, then what it buys |
| the identity for $E_3$ | Theorem | Corollary, attributed |
| the design-space identities and the no-go | here | left in v2 |
| lines | 7,677 | 1,953 |

**v2 is frozen.** Its `paper/theorem_A.md` still contains the spine that
moved here, and that copy is superseded; what still projects from it is
the no-go, which becomes `deploy/papers/P2`. `deploy/PROJECTS_FROM`
records which generation each deployed paper is cut from, and the gate
reads it, so a correction that reaches one tree and not the other is
caught per paper.

---

## The evidence

Fifteen scripts, one result file each. Every printed decimal of three
places or more comes from one of them, and the gate checks that it
comes from the one the sentence names. Statements whose evidence is a
proof rather than a computation say `analytic` and cite nothing.

```bash
python code/<name>.py            # writes results/<name>.txt
python gate/gate.py              # 0 failures is the condition for a commit
```

## What is not claimed

Net progress toward the Goldbach conjecture: zero. The theorem removes
the part of the demand that carries no Goldbach content. The mechanism
— a divisor switch onto a short variable, then Bombieri–Vinogradov — is
ordinary practice, and the ingredients are classical.

Pan's 1982 Goldbach attempt is cited and has not been consulted; the
search behind "nothing states this" is a short list, not a survey, and
the note says so where it says it.
