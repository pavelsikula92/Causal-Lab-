# F6.8d preparation: where the mixed sector departs from covariance

**15 Nov 2026. Analysis of recorded data and a derivation; no new run, decides nothing.**

## 1. The question as recorded

F6.8b2 found two departures of the mixed sector from the conformal one, which a covariant quadratic action forbids:
- the local parts differ, Δ = K_loc,ab − K_loc ≈ 0.46 k²;
- the age coefficients differ, c_T,ab ≈ 0.195 against c_T ≈ 0.110 (ratio 1.78).

The paper and planned_tests ask: "from the age cut, which fixes the rest frame of the field, or from the action?"

## 2. Any covariant quadratic form, with any T-dependent coefficients, gives equal ratios

For the static plane waves ω = a cos kz, ν = b cos kz:
- R⁽¹⁾ = 2k²(3a + 2b) cos kz;
- C⁽¹⁾·C⁽¹⁾ = (4/3) b² k⁴ cos² kz: the Weyl tensor contains b only.

A covariant quadratic action R + R F₀ R + C F₂ C therefore contributes to (A, B, C) in the ratio (6, 8, 2) from R,
(9, 12, 4) φ₀ from R F₀ R, and only to C from C F₂ C. So B/8 = A/6 for any F₀, including an F₀ that depends on the age
T. The ratio of the age coefficients of the mixed and conformal sectors is therefore 1 for any covariant form, and the
recorded 1.78 means that the T-proportional part is not of covariant form. The same holds for Δ.

## 3. The explicit boundary term of the age cut does not carry either departure

In the mixed sector the cut s ≤ l₂(T² − T_xy²) depends on b and gives an explicit boundary term (F6.8b2,
_boundary_term). Refitting the recorded rows of F6.8b2 (results/F6.8b2, 28 Oct) without it:

| k | Δ with the boundary term | Δ without | c_T,ab with | c_T,ab without | c_T (conformal) |
|---|---|---|---|---|---|
| 0.05 | +0.00016 | +0.00040 | 0.155 | 0.160 | 0.110 |
| 0.07 | +0.00233 | +0.00276 | 0.197 | 0.201 | 0.110 |
| 0.1 | +0.00462 | +0.00550 | 0.195 | 0.200 | 0.110 |
| 0.14 | +0.00906 | +0.01077 | 0.195 | 0.199 | 0.110 |
| 0.2 | +0.01844 | +0.02193 | 0.195 | 0.199 | 0.109 |

The boundary term is at most 2.4 × 10⁻³ in K̂ units. Both departures are in the bulk of the past.

## 4. Why the bulk of an aged static field need not be covariant

The age term is the accumulation of long, nearly null pairs (F6.8a): a Lorentz-invariant retarded kernel acting on a
static field does not decay along the past light cone. The sum over the past of an eternal static field therefore
diverges linearly. Any finite age needs a cut, and every cut picks a frame: a covariant weight cannot suppress long
null pairs, because a long null pair is boost-equivalent to a short one. The T-proportional part is thus a property of
the field's history in its rest frame, not of a covariant quadratic action, and §2 does not constrain it. On this
reading the ratio 1.78 is physical (a static field that has existed for a time T in the rest frame of its source) and
not an artefact of how the cut is drawn.

This argument says nothing about Δ, the T-independent part. Δ can come either from short pairs (r of a few l_k), which
would be a non-covariance of the local action at O(k²), or from the long pairs whose accumulation also makes the age
term (the intercept of a non-covariant accumulation).

## 5. The question F6.8d can decide

**Where along the pair length does Δ accumulate?** Compute the angle- and s-integrated bulk densities of the mixed and
conformal sectors as functions of the spatial pair length r, and the cumulative difference D(R) in K̂ units.
- **Local origin:** D(R) reaches Δ by R ≈ a few l_k and stays there.
- **Age origin:** D(R) builds up over r ~ 1/k and beyond, together with the age term.

Development rows come first, then a card.

## 6. Development rows (15 Nov, after §1–5; decides nothing)

diag/f68d_dev.py, outputs diag/out/f68d_dev_k0.2.txt, f68d_dev_k0.1.txt, f68d_dev_k0.2b.txt, f68d_shell_far.txt.

**Mixed bulk density at long pair lengths** (k = 0.2, full s range, angle-integrated):

| r | 20 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|
| ρ_ab | −3.91e−4 | −8.47e−5 | −1.075e−4 | −3.05e−4 | −1.11e−3 | −4.3e−3 |
| levels 1 vs 2 | 5e−8 | 6e−6 | 7e−5 | 1.5e−4 | 1.8e−3 | 1.3e−2 |
| local exponent | | | | 1.50 | 1.87 | 1.96 |
| ρ_c | −3.27e−5 | −6.72e−6 | −4.26e−6 | −3.84e−6 | −3.75e−6 | −3.73e−6 |

- The conformal density tends to a constant, the linear age term.
- The mixed density grows as r², so the mixed age term grows as (kT)³ and not linearly.
- With F6.8e (optical density ∝ r⁴, age term ∝ (kT)⁵) this is a pattern: each factor of the optical index brings
  (kL)², the focusing of long intervals by ν.
- **Consequence.** The linear reading of F6.8b2 (c_T,ab ≈ 0.195; three ages at kT = 8–16) is pre-asymptotic. The
  verdict of F6.8b2 stands (R132).
- **A hypothesis to test, not a result.** A cubic age term of the form k² c₃ n³, read at fixed n by a linear fit
  over the ages, produces an intercept bias ∝ k². The recorded Δ ≈ 0.46 k² has exactly this form, so the local parts
  may be covariant (Δ_loc = 0), with Δ an artefact of the linear age model.

**Cumulative difference without an age cut** (cut in the pair length R instead of the time):

| | R 8 | 17 | 33 | 41 | 60 | 80 | 128 |
|---|---|---|---|---|---|---|---|
| D(R), k = 0.2 | +0.038 | +0.013 | +0.0045 | +0.0030 | −0.0026 | −0.014 | −0.084 |

At k = 0.1, D = +0.021, +0.012, +0.0036, +0.0011, −0.0011 at R ≈ 9, 17, 32, 65, 128. The local part settles to a
small value by kR ≈ 6–10, and then the long-pair accumulation takes over.
