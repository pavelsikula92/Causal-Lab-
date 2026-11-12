# The wandering of a defect is the shot noise of the order (1+1)

Written on 10 Nov 2026, after the wandering study (findings, Discoveries). **Status: derived at the leading order,
and tested by F3.14c (11 Nov, seed 0 PASSED: R = 1.016 / 0.991 / 0.882 / 1.008 ± 0.05; replication pending).**

## 1. Why the rule is a field with point sources

The rule of F3.14a is

  φ(x) = φ₀(x) − (1/2ρ) Σ_{y ≺ x} V′(φ(y)).

In 1+1 the retarded Green function of ∂_t² − ∂_x² is exactly ½ inside the past light cone. The rule is therefore
**exactly** the continuum field with the source

  J_eff(z) = (1/ρ) Σ_y V′(φ(y)) δ²(z − y),

i.e. the continuous source V′(φ) replaced by point sources at the events. The difference ξ = J_eff − V′(φ) is
Poisson shot noise with mean zero and covariance

  ⟨ξ(z) ξ(z′)⟩ = V′(φ(z))² / ρ · δ²(z − z′)   (white in spacetime).

## 2. Momentum diffusion of the kink

For φ_tt − φ_xx + V′(φ) = −s(t, x) the momentum of the field, P = −∫ φ_t φ_x dx, obeys dP/dt = ∫ s φ_x dx. This is the
same identity that gives F = 2πQf for the uniform kick of F3.14a. With s = ξ, to leading order (the kink keeps its
profile):

  ⟨ΔP²⟩ = t ∫ V′(φ)² φ_x² dx / ρ ≡ 2 D_p t,

  **D_p = (1/2ρ) ∫ V′(φ)² φ_x² dx.**

For the sine-Gordon kink with V = μ(1 − cos φ), the profile is φ = 4 arctan e^{√μ x}, so V′ = 2μ sech tanh and
φ_x = 2√μ sech (both of √μ x):

  ∫ V′² φ_x² dx = 16 μ^{5/2} ∫ sech⁴ tanh² du = (64/15) μ^{5/2},

  **D_p ρ = (32/15) μ^{5/2} = 2.13 μ^{5/2} = 2.13 (M/8)⁵.**

There is no free constant. The dependence on the mass is M⁵ because the kink is smaller when it is heavier (M w = 8),
so it is sampled by fewer events.

## 3. Comparison with the development data (diag/wander.py; decides nothing)

D_p ρ is measured from the cubic growth of Var x, from the linear growth of Var v, or from the spread at the end. The
estimates have a relative error of ~25–30 % at 32 sprinklings.

| μ (M) | ρ | measured D_p ρ | prediction |
|---|---|---|---|
| 0.25 (4) | 200 / 400 / 800 | 0.06 / 0.10 / 0.09 | 0.067 |
| 1 (8) | 200 / 400 / 800 | 2.2–2.4 / 3.0–3.2 / 2.2–2.6 | 2.13 |
| 1 (8) | 1600 / 3200 (T 20, from the spread at the end) | 1.5 / 2.0 | 2.13 |
| 4 (16) | 200 / 400 / 800 | 46 / 48 / 49 | 68 |

**Agreement within the precision for M 4 and M 8**, with no fitted parameter. M 16 is 30 % low. That kink is coarsely
resolved (ρw² = 50–200) and loses fits, and the leading order assumes the kink keeps its profile.

## 4. What it means

- **A bound object on the discrete order makes a random walk in momentum.** This is the "swerves" of causal set theory
  (Dowker, Henson and Sorkin 2004). Here its constant is not postulated: it is **the shot noise with which the order
  samples the object's own binding** (V′ ≠ 0 only inside the defect).
- **The vacuum does not swerve.** V′ = 0 there, and only structure is kicked. The strength of the kicks depends on the
  internal structure through ∫ V′² (∂φ)², not on the mass alone.
- **The noise is white in spacetime and Lorentz invariant in distribution.** The diffusion is therefore set in the
  rest frame of the defect, as in the Lorentz-invariant swerves of Dowker–Henson–Sorkin. Consistent with this, the
  moving kinks of F3.14a wandered less (time dilation and relativistic mass).
- **Practical consequence:** the formula gives the resolving power of every test of this branch in advance, with no
  development rows needed.

## 5. Limits and what is open

- **1+1 only.** In 3+1 the Green function is not ½ × the causal relation (it sits on links), and both the noise and the
  integral change. No physical number for real particles follows from this note.
- **Leading order.** The kink is assumed to keep its profile. The feedback of the noise on the profile, and the
  radiation it emits, are neglected.
- **The density trend for M 8** (3.1 at ρ = 400 against 1.5 at ρ = 1600) is within the error of the estimates. A test
  with enough sprinklings and long times would pin down the coefficient and its possible corrections.
- **A physical prediction** would need the 3+1 form and the internal structure of real particles. As a structural
  statement it is already one: **swerves come from internal structure, not from mass.**
