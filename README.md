# A nonergodic flow rule with Galton–Watson stationary measure

Companion material for the paper of the same title (Y. Shmalo, October 2026). [Read the paper](paper/main.pdf).

A flow rule assigns to every infinite rooted tree T a unit flow Θ_T from the root to infinity, consistently with
passage to subtrees: if a child x of the root carries positive flow, the flow through T^x, normalized, is Θ_{T^x}.
It defines a Markov chain on the space of trees that moves from T to the subtree T^x of a child x of the root with
probability Θ_T(x). Lyons, Pemantle and Peres showed that the chain is ergodic when the flow is positive on every child
and a stationary measure absolutely continuous with respect to a Galton–Watson measure GW exists. Question 17.18 of
R. Lyons and Y. Peres, *Probability on Trees and Networks* (Question 3 of Section 10 of Lyons, Pemantle and Peres,
*Ergodic theory on Galton–Watson trees*, 1995), asks whether a flow rule with a stationary measure equivalent to GW
must give an ergodic chain. The paper shows that it need not.

**Theorem 1.** Let GW have offspring distribution p₂ = 3/5, p₃ = 2/5. There are a flow rule Θ and an
isomorphism-invariant Borel set B of trees such that

1. GW is Θ-stationary;
2. GW(B) = 1/2 and B is invariant, so the chain started from any stationary measure equivalent to GW, GW itself
   included, is not ergodic;
3. GW(·|B) and GW(·|Bᶜ) are stationary and ergodic, and GW = (GW(·|B) + GW(·|Bᶜ))/2 is the ergodic decomposition.

The rule follows a colouring of the vertices that is a function of the tree: a vertex with two children takes the
AND of the colours of its children, a vertex with three children their OR, and the root has colour 1 with probability
1/2. The colouring is the limit of the probability q_n that the root gets colour 1 when generation n is coloured at
random; q_n is a martingale, and the probability d_n that two independent random colourings of generation n give the
root different colours satisfies d_{n+1} = d_n (3 − d_n)²/10 ≤ (9/10) d_n. Every vertex has a child of its own colour,
and the flow is spread evenly over the children of the root that share its colour, so B, the event that the root has
colour 1, is invariant. With probability 3/5 some child of the root gets no flow; by the result of Lyons, Pemantle and
Peres, any example must leave some child without flow with positive probability. In the terminology of Johnson, Podder and Skerman the colouring is the
interpretation of the AND/OR tree automaton at the fixed point 1/2; its pivot tree has mean matrix diag(9/10, 9/10).

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/checks.py`: symbolic checks of the identities (SymPy); in exact rational arithmetic, the laws of q₀, …, q₃
  and the first two moments of q₄ against the recursion, the martingale identity on small trees, the mean matrix of
  the pivot tree, and one step of the chain from GW on the first two generations of the tree, computed from the
  definition of the flow rule; and a population-dynamics simulation of the law of q_n for n ≤ 80 (NumPy). It needs
  Python 3 with SymPy and NumPy and runs in under a minute. `code/checks_output.txt` is its output.
- `code/Bounds.lean`: the disagreement identities, the recursion for d_n, the bound d_n ≤ (1/2)(9/10)ⁿ and the
  entries of the pivot mean matrix, checked with Lean 4.33.1 and Mathlib (commit
  0df444a360eaa60ab8c11dca51a86af692955474). `code/lean_output.txt` is the output of the check.
- `SHA256SUMS.txt`: checksums of the files above.

## Running the checks

    python code/checks.py > code/checks_output.txt

For the Lean file, in a Lake project that depends on Mathlib at the commit above:

    lake env lean Bounds.lean

## Citation

Y. Shmalo, A nonergodic flow rule with Galton–Watson stationary measure, 2026,
https://github.com/yspennstate/nonergodic-galton-watson-flow-rule
