import Mathlib

/-!
The arithmetic behind Lemma 1(c) and Remark 3 of the paper.

The Galton-Watson tree has `p₂ = 3/5` and `p₃ = 2/5`; a vertex with two children takes the AND of their colours,
a vertex with three children the OR. Bernoulli(1/2) is a fixed point (`fixed_point`). Two independent colourings of
generation `n` give root colours that differ with probability `dₙ`; a pair of colours with uniform marginals that
differ with probability `d` has `P(1,1) = P(0,0) = (1 - d)/2`, and the AND of two independent such pairs differs
with probability `d - d²/2` (`and_disagreement`), the OR of three with probability `(3d - 3d² + d³)/4`
(`or_disagreement`). So `dₙ₊₁ = dₙ (3 - dₙ)² / 10` (`recursion`), which lies in `[0, 9dₙ/10]` when
`0 ≤ dₙ ≤ 1` (`step`); from `d₀ = 1/2` this gives `dₙ ≤ (1/2) (9/10)ⁿ` (`geometric`), a summable bound
(`summable_bound`). The pivot tree of the fixed point has mean `9/10` for either colour (`pivot_mean_one`,
`pivot_mean_zero`).
-/

theorem fixed_point : (3 / 5 : ℝ) * (1 / 2) ^ 2 + 2 / 5 * (1 - (1 - 1 / 2) ^ 3) = 1 / 2 := by
  norm_num

theorem and_disagreement (d : ℝ) : 2 * ((1 / 2) ^ 2 - ((1 - d) / 2) ^ 2) = d - d ^ 2 / 2 := by
  ring

theorem or_disagreement (d : ℝ) : 2 * ((1 / 2) ^ 3 - ((1 - d) / 2) ^ 3) = (3 * d - 3 * d ^ 2 + d ^ 3) / 4 := by
  ring

theorem recursion (d : ℝ) :
    3 / 5 * (d - d ^ 2 / 2) + 2 / 5 * ((3 * d - 3 * d ^ 2 + d ^ 3) / 4) = d * (3 - d) ^ 2 / 10 := by
  ring

theorem step {d : ℝ} (h0 : 0 ≤ d) (h1 : d ≤ 1) :
    0 ≤ d * (3 - d) ^ 2 / 10 ∧ d * (3 - d) ^ 2 / 10 ≤ 9 / 10 * d := by
  constructor
  · have := mul_nonneg h0 (sq_nonneg (3 - d))
    linarith
  · have hgap : 9 / 10 * d - d * (3 - d) ^ 2 / 10 = d ^ 2 * (6 - d) / 10 := by ring
    have hnn : 0 ≤ d ^ 2 * (6 - d) := mul_nonneg (sq_nonneg d) (by linarith)
    linarith

theorem geometric (d : ℕ → ℝ) (h0 : d 0 = 1 / 2) (hs : ∀ n, d (n + 1) = d n * (3 - d n) ^ 2 / 10) :
    ∀ n, 0 ≤ d n ∧ d n ≤ 1 / 2 * (9 / 10) ^ n := by
  intro n
  induction n with
  | zero => rw [h0]; norm_num
  | succ n ih =>
    obtain ⟨hlo, hhi⟩ := ih
    have hpow : (9 / 10 : ℝ) ^ n ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
    have hle1 : d n ≤ 1 := by linarith
    obtain ⟨s0, s1⟩ := step hlo hle1
    rw [hs n, pow_succ (9 / 10 : ℝ) n]
    constructor
    · exact s0
    · linarith

theorem summable_bound : Summable (fun n : ℕ => 1 / 2 * (9 / 10 : ℝ) ^ n) :=
  (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left (1 / 2)

/-- A pivotal vertex of colour 1: two children both of colour 1 (probability 3/10 given the colour), each pivotal;
or three children with exactly one of colour 1 (three cases of probability 1/10), the 1 pivotal. -/
theorem pivot_mean_one : (3 / 10 : ℝ) * 2 + 3 * (1 / 10 * 1) = 9 / 10 := by
  norm_num

/-- A pivotal vertex of colour 0: children (0,1) or (1,0) (probability 3/10 each), the 0 pivotal; or (0,0,0)
(probability 1/10), all three pivotal. The conditional probabilities: `(3/5)(1/4)/(1/2) = 3/10`,
`(2/5)(1/8)/(1/2) = 1/10`. -/
theorem pivot_mean_zero : (3 / 10 : ℝ) * 1 + 3 / 10 * 1 + 1 / 10 * 3 = 9 / 10 ∧
    (3 / 5 : ℝ) * (1 / 4) / (1 / 2) = 3 / 10 ∧ (2 / 5 : ℝ) * (1 / 8) / (1 / 2) = 1 / 10 := by
  norm_num

#print axioms fixed_point
#print axioms and_disagreement
#print axioms or_disagreement
#print axioms recursion
#print axioms step
#print axioms geometric
#print axioms summable_bound
#print axioms pivot_mean_one
#print axioms pivot_mean_zero
