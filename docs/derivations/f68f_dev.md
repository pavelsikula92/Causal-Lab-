# F6.8f-dev: the optical response of the expected BD action at finite amplitude (16–17 Nov 2026)

Development of the branch F6.8f (docs/derivations/substrate_and_action_audit.md, section 8; roadmap, working order 1).
It decides nothing. It does not change F6.8e3 or F6.8d. The cards built on it, K.18 (instrument) and F6.8f (test), are
frozen by docs/author_decisions.md, decision 34.

## 1. Question and quantities

Does the expected BD action give a usable static response at finite amplitude b and large age T?

- E_opt(k, T, b) = ∂E[S_BD]/∂b and K_eff = [E_opt(b) − E_opt(0)]/b, with the optical projection of F6.8e3.
- The medium is the optical sector of F6.8c: g = −dt² + n² dx², n = e^{b cos(k·x + θ)}, averaged over the phase θ.
- The second variation (F6.8c–e3) holds only for ε² = b²(kT)⁴ ≪ 1. The amplitude parameter of the branch is
  ε = b(kT)².

## 2. The exact first arrival (diag/f68f_dev.py; diag/out/f68f_step0.txt)

**Method.**
- The medium is stratified, so Snell's invariant q = n sin θ is conserved along a ray. The travel time and the
  transverse distance are 1D integrals in z: t = qX + ∫√(n² − q²) dz and X = ∫ q/√(n² − q²) dz.
- Ray families: monotone rays for q < e^{−|b|}; bouncing rays inside the allowed layer around a maximum of n, with the
  phase unfolded (start up or down, end up or down, number of half-periods).
- The first arrival is the smallest time over all families that reach the point. Caustics and multiple arrivals are
  included exactly, with no expansion in b.
- dt/db is taken along the arriving ray (Fermat): ∫ cos(kz + θ) n dl.
- Precision: turning-point integrals in u = m − a cos ψ (no end singularities); cos u − cos a as a product of sines;
  n² − q² through expm1; families parametrised by the half-width a of the layer and λ = −ln(q e^{|b|}). Double
  precision.

**Checks.**

| Check | What | Result |
|---|---|---|
| T1 | b = 10⁻⁹ against the flat distance | differences 1.3e-9, 4.0e-9, 1.7e-9, all first order in b |
| T2 | second-order coefficient, phase-averaged (16 phases, b = ±10⁻⁴), against the closed form t2_mean of F6.8c | 7.5e-10, 1.8e-7, 2.6e-8, 1.0e-7 (d = 3, 5, 2, 40), at the level of the finite difference |
| T3 | b = 0.3 and 0.6, against the polyline minimisation of F6.8c (n = 200 and 400, Richardson) | 6.8e-12, 3.2e-11 |
| T4 | the first caustic from the layer of largest index (b = 10⁻², k = 1) | one arrival (the axis ray) for X < π/(k√b); three from X = 1.01 π/(k√b), the bounced ones earlier than the axis ray (−1.3e-4, −1.1e-2, −0.18 at X/X_c = 1.01, 1.1, 1.5) |
| T5 | dt/db along the ray against a central difference in b | 3.1e-10, 3.0e-10 |

The first caustic for rays along the layer sits at b(kL)² = π², the paraxial value.

## 3. Variables, the BD moments and the age cut

### σ is the Lorentzian distance of the pair
- In the optical sector g = −dt² + h, h = n² dx² (ultrastatic, g_tt = −1).
- A timelike curve from x to y with spatial h-length ℓ in time Δt has proper time ∫√(1 − v²) dt ≤ Δt √(1 − (ℓ/Δt)²)
  (concavity; equality at constant speed).
- The supremum takes the shortest spatial path, ℓ = T_xy, the first arrival (global minimum, caustics included).
  Hence sup τ = √(Δt² − T_xy²) = σ, a coordinate-free quantity expressed through the Killing time t.
- With a conformal factor Ω ≠ 1, σ remains a valid integration variable but is no longer the proper time.
- F6.8c parametrises the past by s = l₂(dt² − T(x, y)²) ("fixed cone"). The γ of F6.8c/e/e3 is therefore at fixed σ,
  and the b² limit of this branch compares with it directly at any L. New at finite b: T_xy, dt and the measure are
  taken exactly.
- At fixed flat τ instead, the optical distance of a long pair shifts at first order by bT₁ with T₁ ~ L; τ² changes by
  about 2LbT₁ ≫ τ². Fixed τ is meaningful only for short intervals (section 5).

### The Jacobian
- s = l₂σ² at fixed spatial end points gives ds = 2l₂Δt dΔt. The pair measure dΔt becomes m ds,
  m = 1/(2l₂Δt), Δt = √(T_xy(b)² + s/l₂).
- m depends on b through T_xy (the terms (ln m)′, (ln m)″ of F6.8c). At finite b it is taken exactly.

### The BD moments
M_j[g] = ∫₀^∞ s^j g(s²) ds for the kernel functions of F6.8a (kernel_derivs; BD coefficients −1, 9, −16, 8; w = s²,
s = l₂τ²):

| j (power τ^2j) | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| g0 = g_k | −1.6e-16 | −3.2e-16 | −4.3e-16 | +0.500 | +2.659 |
| g1 = w g_k′ | +4.9e-16 | +8.9e-16 | +2.1e-15 | −1.000 | −6.647 |
| g2 = w² g_k″ | −2.0e-15 | −3.9e-15 | −9.9e-15 | +3.000 | +23.26 |

- Ṽ/V0 = A + Bσ² + Cσ⁴ + Dσ⁶ + …. The moments remove A, B and C (τ⁰, τ², τ⁴) to rounding; D (τ⁶) is the first
  survivor. At second order A − 1 ∝ ε² is the τ-independent part ∝ (kL)⁴ of F6.8e.
- For a long pair m ≈ (1/(2l₂L))(1 − s/(2l₂L²) + …): its s-dependence moves a coefficient to higher moments with
  factors 1/L².

### The age cut
The past of x is cut at Δt ≤ T. With the flat cut, s ≤ l₂(T² − r²), i.e. r ≤ R(s) = √(T² − s/l₂). Per direction:

    J(T) = ∫₀^∞ ds ∫₀^{R(s)} r² dr m F(s, r),     dJ/dT = ∫₀^{l₂T²} ds (T/R) R² m F(s, R(s)),

with F = Σ g_i(s²) P_i(s, r) the second-order integrand (exact).

- F6.8e/e3 compute the bulk term at r = T over the full s range, ∫ds T² m F(s, T).
- The rest is the layer of pairs with r near T: R − T = −s/(2l₂T) − s²/(8l₂²T³) − …. The layer term of order n
  carries (R − T)ⁿ ∂ⁿ_r[r m F]: n more powers of s and about T⁻²ⁿ.

| Question | Answer |
|---|---|
| Which terms cancel exactly? | In the bulk at full s range, A, B and C for all three kernel functions, to rounding. Truncation at s_max is exponentially small (g ∝ e^{−s²}) except in the layer r ≈ T. |
| Which cancel only as T → ∞? | In the layer: C first contributes at n = 1 (relative ~T⁻²), B at n = 2 (~T⁻⁴), A at n = 3 (~T⁻⁶), each times its own size. |
| Does the edge produce a contribution from A − 1? | Only at n = 3: about (A − 1)M₃s³/T⁶ against the bulk survivor DM₃. A − 1 ~ 10⁹ D at r = 400 (F6.8e card) and T⁶ = 4e15 at T = 400: ratio ~10⁻⁷, falling with T. |
| Is it the optical term of F6.8e3? | No. F6.8e3 is the n = 0 bulk density, the τ⁶ survivor D and higher. The edge terms are a separate, decreasing piece. |

The b-dependence of the cut (the explicit boundary term of F6.8c) is ≤ 2.4e-3 and falls as T⁻⁴ in the k = 0 null of
F6.8c; it is not computed at large n.

## 4. The layer of the flat age cut, numerically (diag/f68f_layer.py)

dJ/dT at second order (flat cut) with the instrument of F6.8e3 (lev1, its s and p nodes), against its bulk term alone.
k = 0.12; age slope −a₄ρ/k⁵ in units of k² per unit n (diag/out/f68f_layer.txt).

| n | bulk (r = T, as F6.8e3) | exact with layer | layer / bulk | (layer / bulk)·n² |
|---|---|---|---|---|
| 24 | 1.109022 | 1.489463 | +0.343 | 198 |
| 96 | 137.4478 | 142.4668 | +0.0365 | 337 |
| 384 | 33 320.06 | 33 399.55 | +0.00239 | 352 |

- The layer is a relative correction c/n², c ≈ 3.5e2: 0.24 % at n = 384, below 1e-4 from n ≈ 2 000. At n = 24 the
  higher layer terms are still large (198).
- Effect on the local exponent of F6.8e3 between n and 2n: log₂[(1 + c/4n²)/(1 + c/n²)] ≈ −1.08c/n² (−1.0e-2 at
  n = 192, −1.6e-4 at 1 536, −4e-5 at 3 072). The limit 5 is unchanged.
- At finite age the layer is not small: a finite-age response needs the full cut, not the bulk.

**Comparison with F6.8c2 at n ≈ 24 (k = 0.12).**
- Secant between n = 16 and 24 (p panels ending at p = n; diag/out/f68f_layer_secant.txt): bulk 0.459039, 0.707912,
  1.109022 and with the layer 0.664136, 0.985576, 1.489463 at n = 16, 20, 24; Simpson secant 0.733 (bulk), 1.016
  (with the layer).
- The F6.8c2 secant (K̂_bb(16) − K̂_bb(24))/(8k²) = (0.741976 − 0.598664)/(8·0.0144) = 1.244 contains the explicit
  b-dependent cut. Its stored boundary term (bnd = −1.3706e-5 at n = 16, −1.8363e-5 at n = 24) is a₄·4.656e-6/(8k⁴)
  = 0.0046 of it, so its interior is 1.2395.
- The F6.8c2 code as a derivative in the age (optical_row of F6.8c, mult 1, the F6.8c2 nodes, T = 197.5 and 199.9;
  diag/out/f68f_c2fd.txt): slope 2.6219 total at n = 23.844, 2.6134 interior, +0.0085 boundary. K̂_bb 0.598992 at
  T = 199.9 against the recorded 0.598664 at T = 200 (the expected trend).
- Bulk plus layer here: 1.489 at n = 24, about 0.57 of the F6.8c2 derivative.
- The integrands are identical: F6.8c `_pair_terms` (bb) and the integrand of diag/f68f_layer.py agree to all printed
  digits at four (CO, r) × four s points (values ±2e2 to ±4e2; the result is a small remainder after the moments).
- The quadratures of the F6.8e3 side are converged at n = 24: p, 48 against 16 nodes, 1.1090252 against 1.109022;
  s, N_S = 16, 24, 32, 48 agree to ~1e-5 relative at p = 0.5, 3, 10 (diag/out/f68f_sconv.txt).
- **Open.** The difference lies in the assembly or the quadrature of one of the two codes, not in the formulas. The
  F6.8c2 code at mult 1 has its only level check (mult 2) at n = 8. The deciding computation is
  `python -m diag.f68f_c2fd --k 0.12 --T 197.5,199.9 --mult 2` (float64). This is development and is not a verdict on
  F6.8c2 (R132).

## 5. Short intervals at finite b (diag/f68f_1a.py; diag/out/f68f_1a.txt)

**Instrument.** V(b, θ) = ∫d³w n³ [dt − T(x, w) − T(w, y)]₊ with T the exact first arrival:
- each segment minimised by Fermat (chord deviation u(1 − u)Σc_jT_j(2u − 1), 14 terms, Newton); against the Snell
  first arrival of section 2: 0 to 5e-16 (six random segments, b up to 0.16);
- prolate coordinates around the chord; the boundary root on every ray by a safeguarded secant on a bracket starting
  above σ;
- θ averaged on 8 phases; flat V0 reproduced to 1e-15.

Fixed σ (dt = √(T(x, y)² + τ²), as F6.8c/e3) and fixed τ (dt = √(L² + τ²)). R(b) = 2(⟨V⟩ − V0)/(b²V0) → γ; γ
extrapolated from b = 0.01 and 0.02 (R = γ + cb²).

| Row (k, CO, L, τ) | γ fixed σ: exact / F6.8e3 | γ fixed τ: exact / second order with dt′ = dt″ = 0 |
|---|---|---|
| 0.7, 0.4, 2, 1.2 (the test diamond of F6.8c) | −0.0159487 / −0.0159484 | 19.26881 / 19.26882 |
| 1, 0, 2, 1 | −0.089168 / −0.089202 | 48.98888 / 48.98885 |
| 1, 1, 2, 1 | 0.1213981 / 0.1213981 | 28.85486 / 28.85487 |
| 1, 0.7, 4, 1.5 | 0.19616 / 0.19600 | 74.22884 / 74.22868 |

- β (fixed σ and fixed τ) agrees to 6 digits in every row.
- The visible differences in rows 2 and 4 (3e-5, 1.6e-4) are the quadrature of the exact volume. The second-order
  instrument is constant to 9 digits over its levels (24,16,16), (32,24,24), (48,32,32); the exact volume moves towards
  it by a factor 3–4 per level (fixed τ; diag/out/f68f_1a_levels.txt): row 2 48.988877 → 48.988855 → 48.988849
  against 48.988846; row 4 74.228837 → 74.228725 → 74.228698 against 74.228682.

**Readings.**
- The second-order instrument of F6.8e3 is confirmed against a computation with no expansion in b and no closed forms
  of the optics, on four short intervals (k ⊥ the pair, k ∥ the pair, kL = 4).
- At fixed τ, γ is 100–1 200 times larger: the shift of the cone, which fixed σ removes (|δ(τ²)|/τ² is 5–10 % at
  b = 0.01). At b = 0.16 some phases give no interval at fixed τ (T(x, y) ≥ dt).
- At fixed σ, R(b)/γ − 1 is −0.2, −0.8, −3.3 % at b = 0.04, 0.08, 0.16 (row 2, ε up to 0.31), +0.3, +0.9, +3.3 %
  (row 4, ε up to 2.56), below 3e-4 (row 3). The departure grows like b², not like ε.

## 6. The scaling limit: needle volumes (diag/f68f_1b.py, f68f_1b_q.py, f68f_1b_exact.py, f68f_1b_cut.py)

### The limit and the ray
- κ = kL → ∞ at fixed ε = bκ² and p = k∥L (cos(k, pair) = p/κ → 0, k⊥ → k); q = kσ.
- With s = l/L ∈ [0, 1] and ψ = kδ (the transverse displacement along k⊥, in radians of phase):
  T = L + (L/κ²)S[ψ] + O(b²),

      S[ψ] = ∫₀¹ [ψ′²/2 + ε cos(ps + θ + ψ)] ds,   ψ(0) = ψ(1) = 0.

- Stationary rays solve ψ″ = −ε sin(ps + θ + ψ) (a driven pendulum). The first arrival is the global minimum of S.
  The ray is bent at finite ε (ψ = O(ε)), so the profile is taken along the bent ray. Every candidate has
  |ψ′(0)| ≤ ε + 2√ε.
- Convergence of the limit: the corrections are in q² and in b = ε/κ²; L/σ = κ/q. The limit needs κ → ∞ and q → 0
  independently.

### The plane wave (q → 0)
- The null Jacobi equation is the Jacobi equation of the optical metric along the ray: J″ = −ε cos(ps + θ + ψ(s))J in
  the plane of k, J″ = 0 across it.
- In Brinkmann coordinates around the ray (x = 0, y = (U, V, 0)), w = (u, v, X) is after x iff v ≥ ½XᵀM₁(u)X and
  before y iff v ≤ V − ½XᵀM₂(u)X; M₁ = E₁′E₁⁻¹ (E₁(0) = 0, E₁′(0) = 1), M₂ = −E₂′E₂⁻¹ (E₂(U) = 0, E₂′(U) = −1).
- ∫d²X[V − ½XᵀKX]₊ = πV²/√det K gives V = πV²∫du (det(M₁ + M₂))^{−1/2}. In the plane of k, M₁ + M₂ = W/(E₁E₂), W =
  E₁(1) constant; across it 1/s + 1/(1 − s). With σ² = 2UV, V0 = πσ⁴/24:

      A = V/V0 = 6 ∫₀¹ √( s(1 − s) E₁(s) E₂(s) / E₁(1) ) ds      (flat: A = 1).

- E₁(1) = ∂ψ(1)/∂ψ′(0) is the slope of the shooting map; E₁(1) = 0 is a conjugate point (a fold where two rays are
  born or die). A diverges only where the first-arrival ray itself has E₁(1) = 0; where the first arrival changes
  branch (the cut locus) A jumps.
- The branch count is computed from the roots of the shooting map, independently of the Jacobi fields. Positivity of
  E₁ on (0, 1] and E₂ on [0, 1) on the first-arrival ray is checked. At p = 0, θ = 0 the RK value equals the closed
  form (E₁ = sin(√ε s)/√ε) to 4e-15 (ε = 5: A = 1.28324).

### Working criteria of the protocol (16 Nov; a development protocol, the basis of K.18)

| Level | κ | b |
|---|---|---|
| 1 | 10√ε | 0.01 |
| 2 | 20√ε | 0.0025 |
| 3 | 40√ε | 0.000625 |

| Check | Criterion |
|---|---|
| volume, q = 1/4 → 1/8 | δ_A = \|Ā(1/8) − Ā(1/4)\| / max(1, \|Ā(1/8)\|, \|Ā(1/4)\|) ≤ 1e-3; Ā the mean over the same 128 phases θ_j = 2πj/128, compared pairwise |
| numerical error of the Jacobi fields and of the volume integral | δ_J ≤ 1e-4 (change of step and tolerance); the quadrature error reported separately |
| brute force at moderate ε | ≤ 1e-3 at ε = 2 and 4, outside caustics, between comparable configurations (finite-κ and q corrections included or extrapolated) |
| conjugate point | \|ε_c(1/8) − ε_c(1/4)\| ≤ 0.01 |

- Recorded with every mean: the distribution of the pairwise differences over θ, the switches of the minimising branch,
  the range of A.
- ε_c (conjugate point of a continued stationary branch, zero of E₁(1)) and ε_cut (the first arrival realised by
  another branch) are separate events. Near a cut locus the smooth criterion for A is not used.
- Grid: ε = 0, 1, 2, 4, 8, 9, 9.5, 9.8, 9.87, 9.95, 10, 10.5, 12, 16, 40; p = 0, 1, 2, 4, 8; q = 1/2, 1/4, 1/8.
  ε = 40 is exploratory if the q-convergence fails.

### Plane-wave results (diag/out/f68f_1b.txt)
128 phases. |A(N) − A(2N)| (N = 4 000) ≤ 1e-10 everywhere except ε = 9.87, p = 0 (1.9e-8). No phase has E₁ or
E₂ ≤ 0 inside the interval on the first-arrival ray.

Onset of a second ray (bisection in ε on the 128 phases):

| p | 0 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|
| ε of the onset | 9.8699 ± 0.0024 | 10.692 | 12.146 | 18.716 | none below 40 |

At p = 0 this is π² = 9.8696.

Ā (mean over θ):

| ε | p = 0 | p = 1 | p = 2 | p = 4 | p = 8 |
|---|---|---|---|---|---|
| 1 | 0.99950 | 0.99952 | 0.99958 | 0.99976 | 1.00002 |
| 4 | 0.99189 | 0.99225 | 0.99324 | 0.99622 | 1.00034 |
| 8 | 0.96405 | 0.96608 | 0.97132 | 0.98491 | 1.00133 |
| 9.87 | 1.57529 (max 85.5) | 0.94056 | 0.95288 | 0.97696 | 1.00202 |
| 12 | 0.86561 | 0.87276 | 0.89868 | 0.96558 | 1.00295 |
| 16 | 0.79911 | 0.80395 | 0.82034 | 0.93408 | 1.00514 |
| 40 | 0.64255 (0.627–0.684) | 0.64456 | 0.65120 | 0.69196 | 1.02802 (0.730–1.366) |

- Small ε: Ā − 1 ∝ ε² (p = 0: −5.0e-4 ε²); the sign changes between p = 4 and 8.
- p = 0, θ = 0, ε = π²: the straight ray is minimising and conjugate (a pitchfork, the cusp of the cut locus); A of
  that phase diverges (85 at ε = 9.87) and is finite on both sides. At p ≥ 1 the folds are generic and the minimising
  ray avoids them.
- ε = 40: Ā falls to 0.64 at small p with a narrow range although up to 5 rays exist; at p = 8 Ā exceeds 1 with a
  single ray. The plane-wave volume does not blow up with ε; the multi-ray structure changes.

### Finite q, expansion around each ray (diag/out/f68f_1b_q.txt)
A(q) = 32√2/(πq⁴) ∫ds √(s(1 − s)) ∫dψ [q²/2 − Φ]₊^{3/2}, Φ = G₁ + G₂ − S_min; the direction out of the plane of k is
integrated exactly. Φ is expanded around every ray to fourth order in the transverse offset (K, c₃, c₄ from the
variational equations of order 1–3), giving the factor 1 + c(5c₃²/(2K³) − c₄/K²), c = q²/2 − (S_b − S_min); every
ray with S_b − S_min < q²/2 adds its needle. K reproduces A_PW to all digits.

- Away from caustics δ_A ≤ 2e-5 at every p and ε, also at ε = 40 (1.3e-5 to 3.9e-5); the largest pairwise phase
  difference is ≤ 3e-4 at ε = 40.
- Near caustics the expansion is not valid (largest relative q² correction at q = 1/4 above 0.1; p = 0, ε = 9.5–10.5;
  weaker at p = 1 near ε = 10 and p = 4 near ε = 16).
- p = 0, ε ≥ 10.5: θ = 0 and π lie on the grid and are exact ties of the two branches born at the pitchfork; a tie on
  a grid point counts both needles at full weight (Ā about 0.011 above A_PW at every q, ε = 12: 0.8766 against 0.8656).
  In the continuum the double count holds only in a window of θ of width ∝ q² (the cut-locus term below).

### Finite q, exact needle (diag/f68f_1b_exact.py)
Φ from the true first-arrival actions: G₁ and G₂ are the lower envelopes over the whole families from x and from y
(3 001 slopes, RK4 2 000 steps), with no ray chosen in advance. The ψ window is ∝ √(s(1 − s)) and widened until
Φ > q²/2 at both ends; the slope range of the families covers the twice-widened window at every node (otherwise an
explicit coverage error); S(ψ) along a family is interpolated by quintic Hermite (d²S/dψ² = (dP/dv)/(dψ/dv)).

- Against the q² series (section 7) at three phases (ε = 4 and 2): ≤ 5e-9 (q = 1/2), ≤ 2e-8 (1/4), ≤ 1.2e-7 (1/8).
- The records of the caustic band below were made with the earlier version (constant window pad, cubic Hermite).
  Against the current version at ε = 9.5, p = 0 they differ by 1e-7 to 4e-7 (q = 1/2) and 2e-6 to 4e-6 (q = 1/8), far
  below δ_A.

Caustic band (Ā over 128 phases; p = 0: diag/out/f68f_1b_exact.txt, .npy; p = 1, 2, 4:
diag/out/f68f_1b_exact_20261116T212730Z, run on the author's machine, 8 processes, 1 354 s):

| ε | p | Ā(1/2) | Ā(1/4) | Ā(1/8) | δ_A | largest pairwise difference |
|---|---|---|---|---|---|---|
| 9.5 | 0 | 0.944737 | 0.947005 | 0.948510 | 1.5e-3 | 0.23 |
| 9.8 | 0 | 0.942161 | 0.950940 | 0.961661 | 1.1e-2 | 1.4 |
| 9.87 | 0 | 0.941976 | 0.953870 | 0.972224 | 1.8e-2 | 2.4 |
| 9.95 | 0 | 0.941720 | 0.956950 | 0.982186 | 2.5e-2 | 3.3 |
| 10 | 0 | 0.941389 | 0.957849 | 0.982021 | 2.4e-2 | 3.1 |
| 10.5 | 0 | 0.930312 | 0.927945 | 0.926008 | 1.9e-3 | 0.22 |
| 9.5 | 1 | 0.947645 | 0.947658 | 0.947677 | 1.9e-5 | 3.6e-2 |
| 9.8 | 1 | 0.942759 | 0.942771 | 0.942500 | 2.7e-4 | 1.4e-2 |
| 9.87 | 1 | 0.941514 | 0.941363 | 0.940838 | 5.3e-4 | 2.8e-2 |
| 9.95 | 1 | 0.940009 | 0.939486 | 0.938553 | 9.3e-4 | 8.1e-2 |
| 10 | 1 | 0.939013 | 0.938123 | 0.936897 | 1.2e-3 | 0.12 |
| 10.5 | 1 | 0.925201 | 0.917686 | 0.916659 | 1.0e-3 | 0.10 |
| 12 | 1 | 0.873560 | 0.872949 | 0.872798 | 1.5e-4 | 4.7e-3 |
| 10.5 | 2 | 0.944144 | 0.943353 | 0.943037 | 3.2e-4 | 1.2e-2 |
| 12 | 2 | 0.901029 | 0.899174 | 0.898793 | 3.8e-4 | 1.9e-2 |
| 16 | 4 | 0.934585 | 0.934257 | 0.934122 | 1.3e-4 | 2.9e-2 |

- p = 0 near the pitchfork (ε = 9.8–10): no convergence in q; Ā grows as q falls towards the divergence of the phase
  θ = 0 at π², finite at every finite q. The smooth criterion does not apply there. The p = 0 means contain the tie
  of θ = 0 on the grid.
- p = 1, 2, 4: δ_A ≤ 1e-3 except p = 1 at ε = 10 (1.2e-3) and 10.5 (1.0e-3), both next to the onset at p = 1 (10.69).
- ε = 16, p = 4: minΦ = −3.4e-7, beyond the self-check −1e-9 (cubic interpolation of the earlier version). The row is
  not used.

### Cut loci (diag/f68f_1b_cut.py)
At a cut locus θ_c the first arrival switches between rays a and b, S_b − S_a = α(θ − θ_c). At finite q both needles
count while |S_b − S_a| < c = q²/2, the second with weight (1 − |S_b − S_a|/c)². Integrated over θ, the extra volume
in the θ-mean per cut locus is dA(q) = (A_a + A_b)c/(6πα), α = |d(S_b − S_a)/dθ| (envelope theorem:
dS/dθ = −ε∫sin(ps + θ + ψ) ds along a fixed-end stationary ray).

p = 0 (diag/out/f68f_1b_cut_partial.txt):

| ε | θ_c | α | A_a = A_b | dA(1/2) | dA(1/4) | dA(1/8) |
|---|---|---|---|---|---|---|
| 9.95 | 0 | 3.203 | 6.0526 | 2.51e-2 | 6.27e-3 | 1.57e-3 |
| 10 | 0 | 4.075 | 4.7790 | 1.56e-2 | 3.89e-3 | 9.72e-4 |
| 10.5 | 0 | 8.880 | 2.2853 | 3.41e-3 | 8.53e-4 | 2.13e-4 |
| 12 | 0 | 15.93 | 1.3952 | 1.16e-3 | 2.90e-4 | 7.26e-5 |
| 16 | 0 | 25.68 | 0.9975 | 5.15e-4 | 1.29e-4 | 3.22e-5 |

- The only cut locus at p = 0 is the tie θ = 0 of the pitchfork. One phase counted twice on the 128-phase grid
  (A_b/128 = 1.395/128 = 0.0109 at ε = 12) is the 0.011 offset of the expansion table.
- Also computed: ε = 9.87, p = 0 has no cut locus; ε = 12 and 16 at p = 1 have one each at θ_c = 5.783185.
- Not computed: ε = 16, p = 2 and the points at ε = 40.

## 7. The age term at finite ε (diag/f68f_1c.py, f68f_1c_e3.py)

### The age derivative and the scaling limit
With the age cut at coordinate time T, dE[S]/dT is the integral over the events y of the past of x at time distance
exactly T with first arrival T_opt(x, y) ≤ T, with weight n(y)³ and the kernel of V(x, y):

    dE/dT ∝ P_x ∫_{T_opt ≤ T} d³y n(y)³ g(V(x, y)),   at Δt = T     (any b; in flat space the layer formula of §3).

In the scaling limit (κ = kT → ∞, b → 0 at fixed ε and p):
- only needles count (σ² = T² − T_opt² ≪ T², r ≈ T);
- n³, P_x and the measure tend to constants, removed by μ₀ = 0; their cross terms with the volume are O(b);
- per direction, with s = l₂σ² and q² = k²s/l₂:

      ρ(p) = (T/(2l₂)) ∫₀^∞ ds ⟨g₀(s² A(ε, p, θ; q))⟩_θ,

  A the needle volume of section 6. For ε below the first caustic there is one ray at every phase and A is analytic
  in q²: A = A₀ + A₁q² + A₂q⁴ + A₃q⁶ + …, A₀ = A_PW.

### The q² series
The coefficients follow from the Taylor series (order 8) of the first-arrival actions around the ray (truncated power
series in RK4, series reversion, Lagrange inversion):

    A_j = (32√2/π)(1/2)^{j+2} B_j (2j + 1) ∫ds √(s(1 − s)) x_{2j+1} (2/K)^{j+1/2},   B_j = Γ(j + ½)Γ(5/2)/Γ(j + 3).

Checks (diag/out/f68f_1c_check_20261116T220043Z.txt): A₀ equals the plane-wave A to 1e-12; A₁ equals the O(q²)
coefficient of the expansion of section 6 to ~3e-8 relative; the linear term of the action sum G₁ + G₂ at the ray
(stationarity) is ≤ 3e-10; against the exact needle as in section 6.

### The BD moments select the q⁶ order
With u = s(A(s)/A₀)^{1/2}, ∫ds g₀(s²A) = Σ_m d_m A₀^{−(m+1)/2} μ_m (d_m the coefficients of ds/du). The moments of
g₀ are μ₀ = μ₁ = μ₂ = 0, μ₃ = 1/2, μ₄ = 2.6587, μ₅ = 10. Every term below s³ vanishes, and

    ρ(p) − ρ_flat(p) = (T/(2l₂)) 4μ₃ (k²/l₂)³ ⟨ŝ₄ A₀^{−2}⟩_θ + O(k⁸),
    ŝ₄ = −e₄ + 5e₂e₃ − 5e₂³,  e₂ = a₁/2,  e₃ = a₂/2 − a₁²/8,  e₄ = a₃/2 − a₁a₂/4 + a₁³/16,  a_j = A_j/A₀.

This is the τ⁶ survivor D of section 3 at finite ε.

### The age function Φ
- With dΩ = 2π dp/κ and T = κ/k: dE/dn ∝ k⁵Φ(ε), ε = bn², Φ(ε) = ∫dp ⟨ŝ₄A₀^{−2}⟩_θ.
- Small ε: Φ = C₂ε², i.e. E ∝ b²n⁵, the (kT)⁵ law of F6.8e3.
- Finite b: E(n) ∝ ∫dn Φ(bn²) = (1/(2√b)) ∫₀^{bn²} ε^{−1/2}Φ(ε) dε. **The static limit of an old field exists if and
  only if ∫^∞ ε^{−1/2}Φ(ε) dε converges** (Φ decays faster than ε^{−1/2}).
- Range: the series needs one ray per phase (ε < π² at p = 0; up to 18.7 at p = 4). The behaviour of Φ at ε → ∞
  needs the multi-ray range, where every cut locus adds terms polynomial in c = q²/2 and the same moment selection
  applies; that extension is not built.

### The b → 0 limit against F6.8e3
At large n the F6.8e3 bulk density per direction must be

    bulk(p) = 4μ₃ k⁵ n⁵ C₂(p)/l₂⁴,   C₂(p) = lim_{ε→0} ⟨ŝ₄A₀^{−2}⟩_θ/ε²,

with no free constant (the F6.8e3 bulk is the full d²/db²; corrections O(1/n²), O(k²)). C₂ from ε = 0.1, 0.2 by
Richardson (diag/out/f68f_1c_20261116T225143Z.txt); k = 0.12.

| n | F6.8e3 / series, p-integrated | record |
|---|---|---|
| 384 | 1.002906 | diag/out/f68f_1c_e3_20261116T231219Z.txt |
| 1 536 | 0.999449 | diag/out/f68f_1c_e3_20261116T233417Z.txt |
| ∞ (c/n² per direction from 384 and 1 536) | 0.99922 | |

- The residual has a part c/n² (+3.7e-3 at n = 384) and a remainder of about −8e-4 that does not vanish with n,
  smooth in p (+7e-4 at p = 0.16 to −2.9e-3 at p = 12.4): O(k²) or the F6.8e3 quadrature, not separated (k = 0.06
  has lost digits; k = 0.085 not computed). The directions p > 12 carry little weight; their F6.8e3 values are
  quantised at ~1e-3 absolute at n = 1 536.
- Two independent codes agree in the b → 0 limit at the 1e-3 level with no free constant.

## 8. Brute force at finite κ and b in the real medium (diag/f68f_1d.py)

The chord-centred volume of section 5 does not reach the needle limit (the ray bends by ~ε/κ, beyond the needle width
q/(2κ)). This instrument works in the real 3+1 medium with no expansion in b, q or 1/κ.

- **Geometry.** The medium is stratified along k: the optical distance depends only on the separation Δz along k, the
  horizontal distance ρ and the start phase.
- **Families.** From x (phase θ) and from y (phase κ cos α + θ), rays as graphs z(ρ) with the exact Lagrangian
  n(z)√(1 + z′²) (RK4; the excess τ = T − ρ carried separately), with the slope derivatives.
- **Tables.** At each ρ node the first arrival τ(Δz) is the lower envelope over the family (quintic Hermite in z), on
  a window around the ray; interpolated by 7-point Lagrange in ρ and cubic Hermite in Δz.
- **Volume.** V = ∫dl dX dY n³h₊, h = Δt − T_x − T_y, Δt = √(T_xy² + σ²), normalised by πσ⁴/24. Per (l, X) the root
  Y* of h = 0 by Illinois regula falsi with a convergence check; Gauss in Y on [0, Y*], a fine grid in X, Gauss in l.
- **Verification.** The tables agree with exact distances to ≤ 1.2e-10 over all 4 297 evaluations of a flat run. Flat
  space (ε = 1e-9, κ = 20, q = 1/4): A = 1.0000031721 with the tables, 1.0000031764 with exact distances.
- **Extrapolation** in 1/κ² over levels doubling κ: two levels (4A₂ − A₁)/3; three levels (16R₂ − R₁)/15 with R₁, R₂
  the two-level values.

Against the q² series (diag/out/f68f_1d_20261116T233144Z, …233451Z, …233644Z, two levels; …20261117T072431Z, three
levels):

| ε | p | θ | q | A(10√ε) − series | A(20√ε) − series | A(40√ε) − series | ratios | extrapolated − series |
|---|---|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 1/4 | −5.66e-4 | −1.41e-4 | | 4.00 | +3.7e-8 |
| 4 | 2 | 2 | 1/4 | −1.19e-3 | −2.97e-4 | | 4.00 | −2.7e-7 |
| 2 | 0 | 0.5 | 1/8 | −1.72e-4 | −4.23e-5 | | 4.07 | +1.0e-6 |
| 4 | 0 | 3 | 1/8 | −1.48e-5 | −2.96e-6 | | 5.01 | +1.0e-6 |
| 4 | 4 | 1 | 1/8 | −3.74e-3 | −9.38e-4 | −2.35e-4 | 3.99, 3.99 | −4.0e-7 (three levels) |
| 4 | 1 | 5 | 1/8 | −4.82e-4 | −1.21e-4 | | 3.98 | −9.1e-7 |
| 4 | 2 | 2 | 1/8 | −1.20e-3 | −2.99e-4 | −7.49e-5 | 4.00, 3.99 | −1.5e-7 (three levels) |
| 2 | 0 | 0.5 | 1/2 | −2.17e-4 | −5.44e-5 | | 4.00 | −6.2e-8 |
| 4 | 0 | 3 | 1/2 | +3.02e-5 | +8.07e-6 | | 3.75 | +6.8e-7 |
| 4 | 4 | 1 | 1/2 | −3.71e-3 | −9.30e-4 | | 3.99 | −3.7e-6 |
| 4 | 1 | 5 | 1/2 | −5.10e-4 | −1.28e-4 | | 3.97 | −1.1e-6 |

- The finite-κ corrections are O(1/κ²): dividing b by 4 divides them by 3.97–4.07 wherever they exceed 1e-4. At
  ~1e-5 the instrument's own ~1e-6 error shows in the ratio (5.01, 3.75).
- After the extrapolation the brute force agrees with the series to ≤ 3.7e-6 at every point. This checks the chain
  rays → Jacobi fields → q² series → scaling limit against the real medium.
- The raw differences (3.7e-3 at 10√ε, 9.4e-4 at 20√ε, 2.4e-4 at 40√ε for ε = 4, p = 4) are physical O(b)
  corrections. K.18 applies the 1e-3 of the protocol to the extrapolated limit (decision 34).
- Cost 16–40 s per volume (about 35 s at 40√ε).

## 9. K_eff (diag/f68f_keff.py)

### From Φ to K_eff
In the scaling limit dE/dn = cΦ(bn²) (section 7); cC₂ is the b² coefficient of F6.8e3. At fixed k and age n = kT:
- E(b, n) = c∫₀ⁿ Φ(bm²) dm;
- E_opt = c∫₀ⁿ m²Φ′(bm²) dm = (c/2b^{3/2}) M(ε),  M(ε) = ∫₀^ε √u Φ′(u) du;
- E_opt(0) = 0 (Φ is even in ε: A(−ε, θ) = A(ε, θ + π));
- K_eff = (c/2b^{5/2}) M(ε);  K₂ = lim_{b→0} K_eff = (2/5)cC₂n⁵, the (kT)⁵ law of F6.8e3.

Hence

    R(ε) = K_eff/K₂ = M(ε) / ((4/5) C₂ ε^{5/2}).

- All K_eff curves of (b, k, T) collapse onto R(ε); the departures from the collapse are the O(b) corrections of
  section 8.
- R = 1 is the second-order extrapolation from F6.8e3.
- At fixed b: K_eff(T) = K₂(T)R(bk²T²); local exponent γ_T = d ln|K_eff|/d ln T = 2 d ln|M|/d ln ε (5 for the second
  order, 0 for a static limit).
- K_eff converges with T iff ∫^∞ √u Φ′ du converges (Φ approaches its limit faster than u^{−1/2}).

### Computation
- By parts: M(ε) = √εΦ(ε) − ½∫₀^ε u^{−1/2}Φ(u) du = √εΦ(ε) − ∫₀^{√ε} Φ(t²) dt.
- The t integral by Gauss–Legendre on panels in t = √u (Φ(t²) is smooth and even in t; in u the integrand is
  ~u^{3/2} at 0).
- Φ(ε) = Σ_p w_p⟨ŝ₄A₀^{−2}⟩_θ from diag/f68f_1c.py.
- Selftest: Φ = u² gives R = 1 and Φ = u³ gives M = (6/7)ε^{7/2}, both to 2e-16 (diag/out/f68f_keff_selftest_*).

### Development rows (diag/out/f68f_keff_20261117T001242Z.txt)
4 directions (p panels 0:12:3, 12:24:1), 16 phases, panels ending at ε = 0.5 and 3 (not judged ε of F6.8f), G = 4.

| ε | Φ(ε) | N = Φ/(C₂ε²) | R = K_eff/K₂ |
|---|---|---|---|
| 0.5 | −1.610e-10 | 1.0125 | 1.0139 |
| 3 | −9.231e-9 | 1.6126 | 1.6888 |

- γ_T(0.5 → 3) = 5.57 (second order 5); K_eff keeps the sign of C₂.
- Checks: N(0.2) = 1.0020; N at ε = 2.4e-3 is 1.000001; Φ = 0 exactly at ε = 0; one ray per phase, stationarity
  2.7e-10; R with G = 3 differs by 4.2e-5 (ε 0.5) and 4.5e-4 (ε 3); Φ with 32 against 16 phases differs by 3e-12
  (ε 0.5) and 2.3e-3 (ε 3), so 16 phases do not resolve ε = 3 (the card uses 128 against 64).
- Indication only: up to ε ≈ 3 the finite-amplitude optical age term grows faster than the second-order extrapolation,
  with no bend towards saturation at this resolution. Whether Φ turns over beyond the first caustic is outside this
  instrument.

## 10. The cards (docs/author_decisions.md, decision 34)

- **K.18** (docs/K18_card_draft.md; phases/K_instruments/K_18_needle_volume.py): the needle-volume instrument on
  ε = 0.1–4, p on the F6.8e3 panels, 128 phases, q = 1/2, 1/4, 1/8. Parts: single ray, Jacobi and A₀, series against
  the exact needle, q-convergence, b → 0 against F6.8e3 at n = 1 536, brute force at κ = 10, 20, 40√ε extrapolated in
  1/κ² (1e-3 on the limit).
- **F6.8f** (docs/F6_8f_card_draft.md; phases/F6_order_dynamics/F6_08f_bd_action_optical_finite_amplitude.py): R(ε)
  at the judged ε = 0.25, 1, 2, 4. Readings in the order B (FAILED: γ_T(2 → 4) ≥ 5, or a confirmed sign change of M,
  R < −1e-3 with both t quadratures), A (PASSED: sign kept, γ_T(2 → 4) ≤ 2.5, nonnegative and decreasing over the
  three intervals), C (UNDECIDED: young / weak fields, no decision on the static limit of an old field); separate null
  test |R(4) − 1| ≤ 0.05, not in the status. The criteria were written before any K_eff row of section 9 existed.
- Logic of the F6.8f readings on analytic Φ (the series replaced by a formula):

| Φ(ε) | γ_T over 0.25→1, 1→2, 2→4 | reading |
|---|---|---|
| −ε² | 5, 5, 5 | B |
| −ε²/(1 + ε²)^{0.8} | 4.18, 2.62, 1.82 | A |
| −ε²/(1 + ε²)^{1.2} | 3.76, 1.23, −0.56 | C (a peak inside ε ≤ 4 fails the nonnegative exponent) |
| −ε²/(1 + 0.3ε) | 4.71, 4.36, 4.02 | C |
| −ε² + 0.5ε³, −ε² + 0.27ε³ | sign change of M | B |

## 11. Open

1. The F6.8c2 code against bulk plus layer at n ≈ 24 (2.61 against 1.49; section 4): the mult 2 levels of the F6.8c2
   code at T = 197.5, 199.9.
2. The remainder −8e-4 of the b → 0 comparison (section 7): another k (0.085) separates O(k²) from the quadrature.
3. The multi-ray range (ε beyond the first caustic): the cut loci at p ≥ 1 and at ε = 40, and their terms in the
   moment selection. Needed for the static limit of an old field (a card F6.8g).
4. The exact-needle row ε = 16, p = 4 (section 6) with the current version.
5. The b-dependence of the cut at large n (section 3).

## 12. Scripts and records

All scripts follow the diag convention through diag/f68f_common.py: `python -m diag.<name>` from the repository root;
`--workers N` via core.stats.parallel.parallel_map; every job cached in diag/out/<name>_cache/ with the md5 of the code
in the key; output diag/out/<name>_<UTC time>.json and .txt with the md5 of the code; nothing overwritten. Scripts on
the 80-bit instrument of F6.8e3 (f68f_layer, f68f_1a) stop without an 80-bit long double (Linux / WSL);
f68f_c2fd runs in float64.

| Script | Computes | Verified against |
|---|---|---|
| diag/f68f_dev.py | section 2, checks T1–T5 | output identical to diag/out/f68f_step0.txt |
| diag/f68f_layer.py | section 4: bulk and layer; p and s quadrature options (`--p-list`, `--n-s`, `--p-panels`, `--bulk-only`) | one direction at n = 24 identical to the cached values; `--p-list 0.5 --n-s 16 --bulk-only` gives −6.0904258135e-06, equal to f68f_sconv.txt |
| diag/f68f_c2fd.py | the F6.8c2 code as a derivative in the age (`--mult`) | k = 0.12, T = 50: K̂_bb 0.8010366682 against 0.8010372305 in the F6.8c2 record of 29 Oct |
| diag/f68f_1a.py | section 5 (`--levels`) | main assembly reproduces f68f_1a.txt; `--levels` identical to f68f_1a_levels.txt (diag/out/f68f_1a_levels_20261117T004948Z.txt) |
| diag/f68f_1b.py | section 6, plane wave | ε = 4, 12 at p = 2 and the onset at p = 1 identical to f68f_1b.txt |
| diag/f68f_1b_q.py | section 6, finite q by expansion | ε = 4, 12 at p = 2 identical to f68f_1b_q.txt |
| diag/f68f_1b_exact.py | section 6, exact needle | the q² series (≤ 1.2e-7) |
| diag/f68f_1b_cut.py | section 6, cut loci | ε = 12, p = 0 identical to f68f_1b_cut_partial.txt |
| diag/f68f_1c.py | section 7, the q² series and Φ | plane wave (1e-12), O(q²) expansion (3e-8), exact needle |
| diag/f68f_1c_e3.py | section 7, b → 0 against F6.8e3 | records at n = 384 and 1 536 |
| diag/f68f_1d.py | section 8, brute force at finite κ | exact distances (4e-9 flat), the q² series |
| diag/f68f_keff.py | section 9, R(ε) | selftest exact |
| diag/run_resumable.py | a resumable start of a frozen test (disk cache of the pure worker functions; the frozen code unchanged) | two K.18 quick passes, data identical to the run without it |

**Reproducibility.** The F6.8c code on Linux (numpy 2.5.3) reproduces the F6.8c2 record of 29 Oct (Windows, numpy
2.2.6) at k = 0.12, T = 50 to 7e-7 relative in K̂_bb. The cause is not identified (numpy version or platform, or a
change of the code between 29 Oct and 10 Nov, before the git history). It is far below the difference of item 1 of
section 11.
