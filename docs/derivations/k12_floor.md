# K.12: what the "10 % floor" is made of – an analytic verdict

Analysis of 24 Oct 2026. An analytic computation and a linear model of the K.12 scheme (scripts in
`docs/derivations/k12_diag/`, the author's runs of 24 Oct). No new test. The criteria and statuses of K.12 do not
change (R132): s16 and s08_R200 FAILED, s32 INSTRUMENT, the reading "converge" = floor 0.099 stays as it was computed.
Only the interpretation of what the floor is changes.

The author's question: is the floor (a) an irremovable property of the discrete operator (non-local screening at the
scale ℓ_k), or (b) a systematic error in the definition of the continuum source or in the mapping onto the grid?

**Answer: neither in pure form.** Non-local screening is analytically excluded for the operator. The floor is a sum
of four components with different behaviour in σ. Three of them are properties of the instrument, not of the causal
set. The fourth (a non-linear remainder) stays open.

## 1. The operator itself: no screening, only an enhancement at the source of order ℓ_k/σ

The expected flat BD operator in 2+1 (coefficients −1, 27/8, −9/4, ερ = 1) acts on a static field as a spatial
convolution with the kernel K(ρ) = ∫₀^∞ τ g(λτ³)/√(τ² + ρ²) dτ. Its Fourier symbol is exactly

    B(k) = pref · [1 + (2π/k) ∫₀^∞ τ g(λτ³) e^(−kτ) dτ].

Expansion in the moments M_m = ∫ τ^(1+m) g dτ:
- M₀ = 0 and M₂ = 0 exactly (numerically 7·10⁻¹⁷ and 1·10⁻¹⁶). Without this the symbol would have terms 1/k and k¹.
- M₁ = 1/(2π) cancels the self term exactly.
- What remains is B(k) = −k² · r(k), where **r(k) = 1 − 1.155 ℓ_k·k + O(k²)**; M₄ ≠ 0 gives an odd power k³,
  typical of an odd dimension.

Numerically r(k) = 0.989 / 0.966 / 0.892 / 0.798 / 0.581 for k = 0.01 / 0.03 / 0.1 / 0.2 / 0.5 (`kern.py`).
The rows of the instrument agree with the symbol to 0.3 % (`check.py`).

Consequences:
1. **r(0) = 1 exactly.** The far field sees exactly the whole charge, Gauss's law holds asymptotically. Screening
   would require r(k) → ∞ for k → 0 (of the type k² + m²). The operator has no such term. Option (a) is thereby
   excluded.
2. r(k) < 1: the operator is *weaker* at short wavelengths, so the field near the source comes out *steeper*. It is
   anti-screening, not screening. The Green function receives an additional term ~ℓ_k/r, so the excess flux falls as
   ℓ_k/r and in a self-similar family as ℓ_k/σ.
3. The exact solution without a grid and without a wall (`hankel.py`), the ratio of the flux to Gauss's law:

   | σ | r = 1.25σ | 3σ | 9σ |
   |---|---|---|---|
   | 8 | 1.159 | 1.064 | 1.016 |
   | 16 | 1.077 | 1.032 | 1.008 |
   | 32 | 1.038 | 1.016 | 1.004 |
   | 128 | 1.009 | 1.004 | 1.001 |

   The excess halves exactly at every doubling of σ and converges to zero. A dimensional argument confirms this in
   general: the operator has a single scale ℓ_k and its continuum limit is □_g − R/2 (K.11). In a self-similar family
   it therefore cannot give a difference independent of σ. Such a difference must come from something that scales with
   σ, i.e. from the instrument.

## 2. A linear model of the K.12 scheme

The curved operator on Π = 1 + u splits at first order into two parts:
- the operator acting on a field that is piecewise linear on the nodes (local limit 1 × ∇²u);
- the response of the operator to the geometry, which is described by a C² spline (local limit 2 × ∇²u, the term
  −R/2).

`lin2.py` / `lin3.py` solve 3∇²u = s in five variants. Collocation and least squares are as in K.12, two rows per
element, the wall is "u = 0 beyond R". In the flat geometry this is an exact consequence of the vacuum subtraction
E[B]f − E[B]1. **Caveat:** the non-local response of the geometry is represented by the flat operator on the spline,
not by the actual response of the curved operator.

Results (author, 24 Oct); grid K.12 = step σ/4 up to 4σ, beyond that ratio 1.12, 28 nodes; fine grid = step 1 ℓ_k:

| variant | σ | u(0)/cont. | slope 1.25σ | flux R/2 | flux 0.9R | Gauss deficit R/2 |
|---|---|---|---|---|---|---|
| everything on the spline, local (grid K.12) | 8–32 | 1.000 | 1.000 | 1.000 | 1.000 | −0.011 (= continuum) |
| field PL non-local + geometry local, fine | 8 / 16 / 32 | 1.022 / 1.009 / 1.001 | 1.049 / 1.025 / 1.013 | 1.005 / 0.999 / 0.995 | 1.002 / 0.990 / 0.974 | 0.004 / −0.001 / −0.002 |
| field PL non-local + geometry non-local, fine | 8 / 16 / 32 | 0.978 / 0.924 / 0.895 | 1.161 / 1.077 / 1.037 | 0.936 / 0.908 / 0.892 | 0.494 / 0.402 / 0.344 | 0.102 / 0.100 / 0.101 |
| the same, grid K.12 | 8 / 16 / 32 | 0.960 / 0.905 / 0.865 | 1.168 / 1.098 / 1.088 | 0.949 / 0.892 / 0.815 | 0.245 / 0.167 / 0.122 | 0.094 / 0.120 / 0.175 |
| everything on the spline, non-local (no PL), fine | 8 / 16 / 32 | 0.984 / 0.928 / 0.898 | 1.158 / 1.075 / 1.036 | 0.943 / 0.913 / 0.896 | 0.526 / 0.424 / 0.357 | 0.096 / 0.095 / 0.097 |
| field PL local + geometry local (limit ℓ_k → 0) | any | 0.845 | 1.147 | 0.647 | 0.483 | 0.294 |
| **K.12 (non-linear, Π(0) ≈ 0.71)** | 8 / 16 / 32 | depth 0.995 / 0.965 / 0.946 | 1.245 / 1.158 / 1.143 | 0.985 / 0.972 / 0.921 | 0.175 / 0.173 / 0.190 | 0.098 / 0.095 / 0.135 |

Check of the grid (σ = 8, step 1 → 0.5 ℓ_k): near the source the result is converged (slope 1.161 → 1.159, flux 3σ
1.054 → 1.056). The far field in the non-local variants is **not** converged (flux 0.9R 0.494 → 0.562, u(0)
0.978 → 0.991). It depends on how close to the wall the rows lie.

## 3. Decomposition of the floor

**(i) The operator, ∝ ℓ_k/σ, vanishes.** On the fine grid the slope at the source is 1.161 / 1.077 / 1.037. These are
the values of the exact solution without a wall (1.159 / 1.077 / 1.038) to 0.002. The wall does not reach the region
of motion.

**(ii) A self-similar grid, grows with σ.** A step at the source of σ/4 means h/ℓ_k = 2 / 4 / 8. The difference of the
slope between the grid of K.12 and the fine grid is +0.007 / +0.021 / +0.051. The response of the grid grows, the
response of the operator falls, and their sum stops. This is the "stopping" of the excess between σ = 16 and 32 that
the reading "converge" interpreted as a floor.

The mapping onto the grid has a deeper flaw still. The field on which the operator acts is piecewise linear.
Collocation of a second-order operator on such functions is inconsistent in the local limit: the second derivative
fits only at the nodes, and the rows inside the elements do not see it (the limit ℓ_k → 0 gives u(0) 0.845). At
h ≲ ℓ_k the non-locality smooths this out fully (PL and spline differ on the fine grid by less than 0.5 %). At h ≫ ℓ_k
it does not. The shape of the grid itself (28 nodes) is, however, not the problem: the local operator on the spline
gives exactly the continuum solution on the same grid.

**(iii) The boundary problem: a non-local operator with a wall against a local Dirichlet reference, independent of σ.**
- The Gauss deficit at R/2 ≈ 0.10 at all σ and on the fine grid. It is present also without the PL basis
  (0.096 / 0.095 / 0.097). It vanishes when the response of the geometry is local (≈ 0.00). K.12 has the same number
  (0.098 / 0.095).
- Mechanism: a logarithmic field has the slope u′(R) = flux/R at the wall. Its continuation by zero beyond the wall
  therefore has a kink, and by Gauss's law the kink carries the whole charge. A local operator does not see the kink
  at the boundary, a non-local one does (`check2.py`: B[r² − R²] = 4.02 inside, 3.20 at 20 ℓ_k from the wall, −7.2 at
  5 ℓ_k). Since the strength of the kink is proportional to the total flux and not to ℓ_k, the effect does not scale
  away in the family R = 12.5σ.
- The solution is by least squares over the whole domain, so it is carried inwards, into the flux at R/2 and into
  the depth of the core. The dependence on the step at the wall (point 2) confirms that it is a property of the
  implementation of the wall.
- In 3+1 the slope at the wall is ~flux/R², which is why the same instrument passed there (K.11c far_d4).
- In nature there is no wall. The correct reference would be the same non-local operator on an unbounded domain, or
  with an exterior continuation of the field.

This also explains the depth of the core, which moves *away* from the continuum value with σ. The wall gives a
σ-independent deficit far out (a shallower core), the operator a σ-decreasing excess at the source (a deeper core). At
σ = 8 they nearly cancel; with growing σ only the wall remains (model 0.978 → 0.924 → 0.895, K.12 0.995 → 0.965 →
0.946).

**(iv) A non-linear remainder, open.** The slope at 1.25σ in K.12 is 0.077 / 0.060 / 0.055 above the linear model on
the same grid. It falls slowly and the linear model does not explain it. It comes from the strong field (qσ² = 0.64,
Π(0) ≈ 0.71), i.e. from the non-linear self-consistent coupling or from the non-linear part of the response of the
geometry.

Summary for the slope at 1.25σ (the region of motion, which decides P_cs/P_a):

| σ | K.12 | operator (i) | grid (ii) | non-linear remainder (iv) |
|---|---|---|---|---|
| 8 | 1.245 | +0.159 | +0.007 | +0.077 |
| 16 | 1.158 | +0.077 | +0.021 | +0.060 |
| 32 | 1.143 | +0.038 | +0.051 | +0.055 |

## 4. Verdict

1. **Option (a), irremovable non-local screening: excluded.** The symbol is exact, r(0) = 1, the operator
   anti-screens at the source and its contribution vanishes as ℓ_k/σ (the fine grid confirms it to 0.002).
2. **Option (b), definition of the continuum source: excluded.** The source and the normalisation fit. The variant
   with a local geometry converges to u(0) 1.0006 and a flux of 0.995 at R/2 at σ = 32.
3. **Mapping onto the grid: yes, in part.** A self-similar grid (h = σ/4) gives an error growing with σ. This is the
   mechanism of the stopping that the reading "converge" interpreted as a floor. Aitken's estimate 0.099 is therefore
   not an asymptote, it is the crossover between a falling and a growing error.
4. **The main origin of the far deficit and of the shift of the core: the boundary problem.** A non-local operator
   with a wall is compared with a local Dirichlet problem. Shown in the model; in the real curved instrument only
   probable (the same numbers, but the model represents the response of the geometry).
5. **K.12 is not a falsification of the naive continuum of Path A.** The instrument K.12 could not measure the
   continuum limit of the operator, because with σ it also scaled the grid and the wall, and both give errors that do
   not vanish with ℓ_k/σ. Nor is it a confirmation: the non-linear remainder ~0.055 at σ = 32 stays unexplained.

## 5. What will decide (the author's runs)

Only the region of motion (for F2.7) and the non-linear remainder (iv) decide:
- **a strong field with the step at the source fixed in ℓ_k (≤ 1 ℓ_k instead of σ/4)** for σ = 16 and 32. The
  prediction of the linear decomposition for the slope at 1.25σ without the grid component: ≈ 1.137 and ≈ 1.093. If
  the excess keeps falling as ℓ_k/σ plus a small remainder, F2.7 can be repeated at a σ where |P_cs − P_a| ≤ P_a/16;
- **a weak-field family w08/w16/w32 with the same grid** (declared with step σ/4 it would inherit component (ii)).
  It separates the non-linear remainder: in a weak field the slope at 1.25σ should agree with the fine linear model;
- the deficit in the far field need not be corrected for F2.7. For any measurement of the far field in 2+1 (F1.9 in
  2+1, flux, Π(0)), however, another boundary problem is needed: an exterior continuation of the field beyond the
  wall instead of zero, or a comparison with a non-local reference.

**Result of 25 Oct (h16, step 1 ℓ_k):** P_cs/P_a = 1.111 (s16: 1.124), slope at 1.25σ 1.141 against the prediction
1.137, grid difference +0.018 against the model +0.021, the reading 'grid' supported. The decomposition thus held in
its first prediction.

Changing the parts of K.12 before a run is a new hypothesis after a failure (R132 allows it); the old criteria do not
change. No code yet.
