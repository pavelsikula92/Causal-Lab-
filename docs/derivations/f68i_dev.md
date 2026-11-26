# F6.8i-dev: the static response of an old field (22 Nov 2026; development, decides nothing)

The author chose this test as the next one (22 Nov, 09:16). F6.8i is the paper's planned test "the optical age term of
old fields, ε = b(kT)² ≫ π²: intervals with several rays; does an old field have a finite static response?"
(decisions 45, 47). It continues F6.8f (FAILED, reading B: no saturation up to ε = 4) and uses its machinery
(docs/derivations/f68f_dev.md).

## 1. The question in one formula

In the scaling limit (κ = kT → ∞ at fixed ε = bκ² and p = k∥T), the optical age term of the expected BD action is
(f68f_dev §7, §9):

    dE/dn ∝ k⁵ Φ(bn²),     Φ(ε) = ∫ dp ⟨ŝ₄ A₀⁻²⟩_θ,
    E(b, n) ∝ (1/(2√b)) ∫₀^{bn²} ε^{−1/2} Φ(ε) dε.

- **The static limit of an old field exists if and only if ∫^∞ ε^{−1/2} Φ(ε) dε converges**, i.e. if Φ falls faster
  than ε^{−1/2} at large ε.
- Equivalently for the stiffness: K_eff = K₂ R(ε) with R(ε) = M(ε)/((4/5)C₂ε^{5/2}), M(ε) = ∫₀^ε √u Φ′(u) du; a static
  limit needs M(ε) → const.
- F6.8f measured R = 1.003, 1.056, 1.248, 2.518 at ε = 0.25, 1, 2, 4: Φ grows faster than ε².

**F6.8i asks what Φ does beyond the first caustic (ε > π² at p = 0), where an interval contains several rays.**

## 2. What exists and what is missing

**Exists (f68f_dev §6–§9; K.18 PASSED):**
- the stationary rays of the driven pendulum ψ″ = −ε sin(ps + θ + ψ), all branches, the first arrival as the global
  minimum of S (caustics included exactly);
- the needle volume A(ε, p, θ; q) as a q² series around the first-arrival ray (diag/f68f_1c.py), checked against the
  exact needle and against brute force at finite κ (≤ 3.7e-6 after extrapolation);
- the plane-wave volume beyond the caustic up to ε = 40 (f68f_1b: Ā falls to 0.64 at p = 0, ε = 40; up to 5 rays);
- the cut loci at p = 0 and their O(q²) window term dA = (A_a + A_b)c/(6πα), c = q²/2.

**Missing:**
1. **The cut loci at the q⁶ order.** At a cut locus θ_c two rays tie; for finite q both needles count in a window
   |θ − θ_c| ≲ c/α. The BD moments remove every term below s³ for each phase, and they also remove the leading window
   term: integrated over θ the window gives (c/α)·∫dz[g₀(s²(A_a + (1 − z)²A_b)) − g₀(s²A_a)] with c ∝ s, and
   ∫ds s g₀(λs²) = λ⁻¹μ₁ = 0 for every λ. The window therefore contributes first at s³ (q⁶), the same order as the
   smooth part, through the q-dependence of A_a, A_b and the curvature of S_b − S_a in θ. These terms have to be derived
   (a series in c around each cut locus) and checked against the exact needle (f68f_1b_exact) at finite q.
2. **Near conjugate points** (K → 0 on the first-arrival ray) the series in q² is not valid. At p ≥ 1 the folds are
   generic and the first arrival avoids them (f68f_1b); at p = 0 the pitchfork at ε = π² is a cusp. The window around
   each such point needs a bound or a separate treatment.
3. **The p-integral at large ε.** The directions p ≫ √ε have one ray (p = 8: one ray up to ε = 40) and carry the
   tail of Φ; their weight has to be resolved.

## 3. Exploration E1: the smooth part beyond the caustic

The first-arrival ray per phase, without the cut-locus windows (diag/f68f_1c.py point, 16 phases), at
ε = 2–40 and p = 0.5, 2, 4, 8, 16 (Claude's container; scratch run, not a record). Results in §3.1 when it has run.
This is an indication of the trend of the smooth part only; the cut-locus terms are not included, and 16 phases do not
resolve ε ≳ 3 at the 1e-3 level (f68f_dev §9).

### 3.1 Results (22 Nov, 09:3x; Claude's container, 2 processes, ~20 s per point; scratch, not a record)

⟨ŝ₄A₀⁻²⟩_θ on 16 phases (the "smooth part": the first-arrival ray at every phase, no cut-locus windows). In brackets
the number of stationary rays at the worst phase, and min K on the first-arrival ray (K → 0 at a conjugate point,
where the q² series fails).

| ε | p = 0.5 | p = 2 | p = 4 | p = 8 | p = 16 |
|---|---|---|---|---|---|
| 2 | −6.5e-10 | −5.6e-10 | −3.6e-10 | −8.5e-11 | −2.1e-14 |
| 4 | −9.5e-9 | −5.9e-9 | −2.2e-9 | −3.4e-10 | −8.4e-14 |
| 8 | +3.0e-6 (±2.2e-6) | −4.1e-7 | +4.4e-8 | −1.4e-9 | −3.7e-13 |
| 12 | −2.0e-6 | −2.6e-6 | +2.5e-6 (min K 0.99) | −3.5e-9 | −9.6e-13 |
| 16 | −1.7e-7 [3 rays] | −5.3e-7 | **+8.2e-4 (min K 0.46)** | −6.7e-9 | −2.0e-12 |
| 24 | −1.2e-8 [3] | −4.2e-8 [3] | −1.4e-5 [3] | −1.8e-8 | −6.6e-12 |
| 40 | −1.9e-9 [3] | −6.5e-9 [3] | −4.0e-7 [3] | −7.7e-8 (min K 1.0) | −3.8e-11 |

Readings, with their limits (24 Nov: superseded as a trend by §5.2, the cut-locus windows are of the same order as
the smooth part):
- **Where the interval has several rays, the smooth part falls steeply with ε:** p = 0.5 from −1.7e-7 (ε 16) to
  −1.9e-9 (ε 40), p = 2 from −5.3e-7 to −6.5e-9, roughly ε^(−5). Up to the first caustic it grew faster than ε².
- **Where there is still one ray, it keeps growing,** roughly as ε^(2.3–2.5) (p = 8 and 16), until its own caustic.
  The single-ray range moves to larger p as ε grows (onset p = 1: 10.7, p = 2: 12.1, p = 4: 18.7, p = 8: none below 40),
  while the second-order weight C₂(p) falls very fast with p (p = 8 → 16: factor ~4 000 at ε = 2).
- So Φ(ε) could turn over once the directions that carry the weight have passed their caustic. **This is an
  indication, not a result:**
  - near conjugate points the series diverges (ε = 16, p = 4: +8.2e-4 from one phase with K = 0.46, four orders above
    its neighbours); these windows need D2 before any number of Φ can be trusted;
  - the cut-locus terms (D1) are missing, and they enter at the same order;
  - 16 phases: the phase spread is larger than the mean in several cells (e.g. ε = 8, p = 0.5: mean +3.0e-6, standard
    error 2.2e-6);
  - five directions only; no p-integral.

## 4. Plan

| Step | What | Check |
|---|---|---|
| E1 | the smooth part of ⟨ŝ₄A₀⁻²⟩ beyond the caustic, coarse grid | trend only |
| D1 | the cut-locus window at q⁶: series in c around each cut locus (envelope theorem for S_b − S_a; q-series of A_a, A_b); checks 1–3 done (§5.1, §5.2) | against the exact needle (f68f_1b_exact) at q = 1/2, 1/4, 1/8 on cut-locus phases |
| D2 | conjugate points: the window bound or the uniform (Airy/cusp) form | against the exact needle near the folds |
| D3 | Φ(ε) on ε up to ~40 (or until a trend is decided), p panels to p ~ 3√ε + 12 | phases 64 against 128, p panels, the b → 0 limit against F6.8f |
| Card | F6.8i: the reading of ∫ε^{−1/2}Φ (convergent: a static limit exists; divergent: none), with criteria fixed before the run | as K.18 / F6.8f |

**What F6.8i can decide.** Whether the optical static response of an old field converges with age (a static limit)
or keeps growing. If it converges, the old-field stiffness is finite and the far-field law and the short-range
prediction (tens of micrometres) can be computed. If it diverges, the static response of an old field is not a
property of the bare expected action, and the age has to be fixed by something else (decision 41: the balance of
action and entropy; the meaning of the age, F8.8).

## 5. D1: the cut-locus window at the q⁶ order (derivation, 22 Nov; not yet implemented)

**Setting.** Write the needle volume through the well function W(c) = ∫₀¹ds √(s(1−s)) ∫dψ [c − Φ(s, ψ)]₊^{3/2},
c = q²/2, Φ = G₁ + G₂ − S_min, so that A(q) = (8√2/π) W(c)/c² (f68f_dev §6). Near one ray, W(c) = Σ_j w_j c^{j+2} and
A = (8√2/π) Σ_j w_j c^j; the coefficients A_j of f68f_1c are these w_j up to the factor 2^j.

**Two wells.** At a cut locus θ_c the rays a and b tie; let d(θ) = |S_b(θ) − S_a(θ)| and let a be the first arrival.
At every s, Φ has a well of depth 0 around ψ_a(s) and a well of depth d around ψ_b(s), provided the half-segments of
ray b are themselves first arrivals from x and to y (to be checked on every cut locus; otherwise the second well is
shallower than d over part of s). For c > d both wells count:

    A(θ; q) = (8√2/π) [ Σ_j w_j^a c^j + Σ_j w_j^b (c − d)^{j+2} / c² ],

the second sum only for c > d. Away from the window (d > c) A is the one-ray series.

**The window integral.** With d(θ) = α|θ − θ_c| + β(θ − θ_c)² + γ₃|θ − θ_c|³ + … on each side (α, β, γ₃ one-sided),
the θ-integral over the window is, per side, an integral over d ∈ [0, c] with dθ = (1/α − 2βd/α³ + …) dd. Put d = cz:

    Δ(s) = Σ_sides ∫₀¹ c dz (1/α − 2βcz/α³ + …) [g₀(s²(A_a + X)) − g₀(s²A_a)],
    X = (8√2/π) Σ_j w_j^b c^j (1 − z)^{j+2},

with A_a, w^b and their θ-derivatives evaluated at θ_c + (side) cz/α + … (the variation of A_a and w^b across the
window enters at the next orders in c).

**Moment selection.** In the age term c = γs with γ = k²/(2l₂). A term of order γ^m in Δ has the form
∫ds s^m G(λs²) with G a combination of g₀, g₁ = w g₀′, g₂ = w² g₀″, … at scales λ that depend on z but not on s. The
BD moments of g₀, g₁, g₂ vanish for j ≤ 2 (f68f_dev §3), so **every window term of order γ¹ and γ² vanishes exactly**
(the leading one by μ₁[g₀] = 0, for every z), and the first survivor is γ³, the order of the smooth part.

**What the γ³ window term needs, per cut locus:** α, β, γ₃ on both sides; A_a and w^b up to second order in c (w₀, w₁,
w₂ of both rays); the first and second θ-derivatives of w₀^a, w₀^b, w₁^a, w₁^b at θ_c. All come from the existing
series code at θ_c and two or three neighbouring phases (the envelope theorem gives dS/dθ along each ray, f68f_dev §6).

**The smooth part with the jump.** ⟨ŝ₄A₀⁻²⟩_θ of the first-arrival ray jumps at θ_c (the first arrival changes
branch). Equally spaced phases converge only as 1/N_θ across a jump; the smooth part must be integrated piecewise
between the cut loci (Gauss–Legendre on each piece), with θ_c located by bisection on S_b − S_a.

**Checks for the implementation.**
1. On a cut-locus phase at finite q, the two-well A(θ; q) against the exact needle (f68f_1b_exact) at q = 1/2, 1/4, 1/8.
2. The window integral of A over θ at O(c) against dA = (A_a + A_b)c/(6πα) of f68f_dev §6 (the leading term).
3. The total ∫dθ (smooth + windows) of a quantity whose answer is known: at small ε (one ray, no cut loci) it must
   equal the present code; the γ¹ and γ² window terms must vanish numerically (exact moments).

### 5.1 Check 1 done (22 Nov, 11:3x; diag/f68i_d1.py, diag/out/f68i_d1_check1_20261122.txt)

ε = 12, p = 1, cut locus θ_c = 5.783185 (three stationary rays): the two-well series against the exact needle at five
phases across the window (d = 1.4e-1, 2.9e-2, 4.4e-6, 2.9e-2, 1.4e-1).

| d | q = 1/2 | q = 1/4 | q = 1/8 | one ray only (q = 1/4) |
|---|---|---|---|---|
| 0.144 (outside every window) | −8.9e-5 | −2.9e-7 | −8.3e-9 | the same |
| 0.029 (inside for q = 1/2, 1/4) | −1.8e-4 | −6.6e-7 | +1.0e-7 | −9.7e-3 |
| 4.4e-6 (the tie) | −3.6e-4 | +1.3e-5 | +2.6e-7 | −1.5 |

- **The two-well formula reproduces the exact needle** everywhere: the residual is the truncation of the series at
  q⁶ (it falls by ~300 from q = 1/2 to 1/4, as q⁸), largest at the tie.
- One ray alone misses the volume by up to 1.5 (a factor two) inside the window: the second needle is real and its
  weight (1 − d/c)² as derived.
- The windows are symmetric (θ_c ± 0.002 agree to 7e-7), as expected for a crossing.

### 5.2 Checks 2 and 3 (24 Nov; diag/f68i_d1_checks.py, diag/out/f68i_d1_checks_20261124T173517Z)

**Method.** Three cut loci, all symmetric ties θ_c = −p/2 (mod 2π), where the two lowest rays are mirror images:
(ε, p) = (12, 1), (16, 1), (12, 0). On each side of θ_c the rays, the gap d = S_second − S_first and the series
A₀…A₃ of both rays are computed at 9 Chebyshev nodes in δ = |θ − θ_c| ∈ (0, 1.25 c_max/α] (c_max = 1/8, q = 1/2) and
fitted by polynomials of degree 5 (d without a constant). Model checks, at every point: fit residual ≤ 3e-10 (d) and
≤ 7e-7 (A); a held-out node ≤ 4e-10 and ≤ 4e-7; continuity across θ_c (first on one side = second on the other)
≤ 5e-7; the third ray stays above c_max (gap ≥ 0.27).

**Check 2: the window integral over θ at finite q, with the exact needle.** (1/2π)∫dθ [A_exact − A_first-series] at
6 Gauss–Legendre nodes per side inside the window d < c, against the leading term dA = (A_a + A_b)c/(6πα) and against
the two-well model:

| (ε, p) | α | exact/dA − 1 at q = 1/2, 1/4, 1/8 | ratio per halving of q | exact − two-well (q = 1/2 … 1/8) |
|---|---|---|---|---|
| (12, 1) | 14.376 | 2.8e-2, 6.2e-3, 1.5e-3 | 4.49, 4.09 | 4.4e-7 … 1.0e-11 |
| (16, 1) | 24.426 | 2.2e-3, 5.4e-4, 1.3e-4 | 4.00, 4.04 | −1.5e-8 … −3.9e-11 |
| (12, 0) | 15.927 | 1.7e-2, 3.9e-3, 9.7e-4 | 4.29, 4.06 | 1.7e-8 … −3.6e-11 |

The leading window term of f68f_dev §6 holds, with corrections O(c) (factor 4 per halving of q); the two-well model
equals the exact needle to the truncation of the series. **Passed.**

**Check 3: the window in the age term.** With c = γs (γ = k²/(2l₂)):
- 3a, the moment selection of one ray with this file's quadrature: F(γ) = ∫ds g₀(s²A(γs))/(4μ₃(2γ)³) tends to
  ŝ₄A₀⁻² to 8e-6, 9e-5, 3e-6 (the second limited by rounding: ŝ₄A₀⁻² is 1.4e-5 there).
- 3b, the window part W(γ) = (1/2π)Σ_sides∫ds∫dδ [g₀(s²A_two-well) − g₀(s²A_first)] at γ = 0.0008 … 0.009, fitted by
  powers γ¹ … γ⁶. N is the same integral of |·|, the size of the γ¹ piece the moments remove.

| (ε, p) | W₁/N₁ | W₂/N₁ | W₃ (fit m 1–6; m 3–6 agrees to) | window term in ⟨ŝ₄A₀⁻²⟩_θ, W₃/(4μ₃·8) | smooth part, 64 equal phases (scale only) |
|---|---|---|---|---|---|
| (12, 1) | 1e-13 | −3e-10 | 1.8387e-4 (1e-5) | **+1.149e-5** | −4.93e-5 |
| (16, 1) | −7e-17 | 1e-13 | 1.9019e-6 (4e-7) | **+1.189e-7** | −4.18e-7 |
| (12, 0) | 2e-14 | −4e-11 | 8.3295e-5 (3e-6) | **+5.206e-6** | −3.16e-5 |

Doubling the s and z nodes changes W by ≤ 1.4e-6 (relative). **Replicated** on the author's PC (Windows, 6 workers,
994 s, diag/out/f68i_d1_checks_20261124T185616Z): check 2 identical in every printed digit, W(γ) to 7 digits, W₃ to
≤ 6e-6, the window terms identical in all five printed digits; only check 3a at (16, 1), the rounding-limited one,
moved (−4.4e-4 against −9.4e-5). **Passed:** the γ¹ and γ² window terms vanish (to 1e-13
and 3e-10 of the γ¹ size), and the first survivor is γ³, as derived.

**Readings (development, decides nothing).**
- **The window term is not negligible.** Per cut locus it is +1e-5 … +1e-7 in ⟨ŝ₄A₀⁻²⟩_θ, positive at all three points,
  against a negative smooth part; it is 16–28 % of the 64-phase smooth value. The size ratio is not yet known: the
  smooth part on equal phases does not converge across the jump at θ_c ((12, 1): −4.8e-6 on 16 phases (first run of
  this file), −4.9e-5 on 64).
- **E1 (§3) does not hold even as a trend in the multi-ray range:** it has no windows and 16 equal phases. Between
  ε = 12 and 16 at p = 1 both the smooth part and the window fall by about two orders of magnitude, which is consistent
  with a steep fall, but nothing is decided before the piecewise integration.
- **Limits:** symmetric ties only (mirror rays, the same α on both sides). Non-symmetric cut loci (expected at larger ε
  with more rays) are not tested; the model needs only "first" and "second" per side, so it applies to them unchanged,
  but the continuity check then compares different rays and has to be read per locus. The conjugate-point windows (D2)
  are untouched.

**Next.** The smooth part integrated piecewise between the cut loci (θ_c located by the tie, Gauss–Legendre on each
piece), plus W₃/(4μ₃·8) per cut locus: Φ-contributions at fixed (ε, p) with both parts. Then D2.

### 5.3 D3a: ⟨ŝ₄A₀⁻²⟩_θ at fixed (ε, p) with the cut loci (24 Nov; diag/f68i_d3a.py, diag/out/f68i_d3a_20261124T224000Z)

    ⟨ŝ₄A₀⁻²⟩_θ = (1/2π) Σ_pieces ∫_piece dθ ŝ₄A₀⁻² (first arrival)  +  Σ_loci W₃/(4μ₃·8)

The cut loci are those of diag/f68f_1b_cut (128-phase grid and bisection). Between neighbouring loci the first arrival
is one smooth branch.

**A finding about the integrand.** ŝ₄A₀⁻² of the first arrival peaks sharply at the cut loci and stays finite there.
At (12, 1) it is ≈ −1e-5 at |θ − θ_c| ≈ 0.2, −6e-4 at 0.04 and −3e-3 at 0.01 (A₀ rises from 1.05 to 1.51). Most of
⟨ŝ₄A₀⁻²⟩_θ comes from within ~0.1 rad of the cut locus: the same place where the windows act.
- One Gauss–Legendre rule per piece does not resolve this: 12 and 24 nodes differ by a factor 2 at (12, 1) (first run of
  this file, not kept).
- Each piece is therefore cut into 18 subintervals graded geometrically towards both cut loci (first 0.002 from θ_c,
  ratio 2.5), with Gauss–Legendre at 6 and 12 nodes on each.

| (ε, p) | smooth part (6 / 12 nodes differ by) | windows | **total** | 64 equal phases (E1 method) | min K |
|---|---|---|---|---|---|
| (4, 1), no cut locus | −6.378459e-9 (8e-17) | 0 | **−6.378459e-9** | −6.378459e-9 | 2.60 |
| (12, 1) | −3.176167e-5 (8e-13) | +1.149e-5 | **−2.0270e-5** | −4.93e-5 | 1.52 |
| (16, 1) | −3.916375e-7 (1e-14) | +1.189e-7 | **−2.7277e-7** | −4.18e-7 | 3.84 |
| (12, 0) | −1.461195e-5 (4e-13) | +5.206e-6 | **−9.406e-6** | −3.16e-5 | 1.88 |

**Readings (development, decides nothing).**
- **The method converges.**
  - The graded smooth part is stable to ≤ 3e-8 (relative) between 6 and 12 nodes.
  - The control without a cut locus agrees with equal phases to all digits.
  - No piece approaches a conjugate point: min K ≥ 1.5 against the flag at 0.5.
- **The windows remove 30–36 % of the smooth part**, with the opposite sign.
- **The E1 method is off by a factor 1.5–3.4** in the multi-ray range: equal phases across the jump and no windows.
- **At p = 1, ⟨ŝ₄A₀⁻²⟩_θ falls by a factor 74 between ε = 12 and 16** (a local exponent ≈ −15).
  - This is one direction p only.
  - Φ(ε) = ∫dp⟨ŝ₄A₀⁻²⟩_θ needs the whole p range, and E1 indicates that the single-ray directions at large p still grow.
  - The p-scan at fixed ε is the next run.
- **Limits:**
  - Cut loci closer than one grid step (0.049) can be missed; a missed pair would show as a large 6/12 difference on
    its piece.
  - The window model assumes a crossing with nonzero α. Near a cusp of the cut locus (where α → 0 and conjugate points
    appear) it does not apply. That is D2.

### 5.4 D3a p-scan at ε = 12 and 16 (25 Nov; run on the author's PC, 8 workers, 3684 s; diag/out/f68i_d3a_20261125T002324Z)

⟨ŝ₄A₀⁻²⟩_θ (smooth part + windows) over p at fixed ε. [n] = number of cut loci found on the 128-phase grid.

| p | ε = 12 | ε = 16 |
|---|---|---|
| 0 | −9.406e-6 [1] | −1.796e-7 [1] |
| 0.5 | −1.128e-5 [1] | −1.991e-7 [1] |
| 1 | −2.027e-5 [1] | −2.728e-7 [1] |
| 2 | −3.27e-4 [0], 6/12 differ 2.7e-4: **not converged** | −1.055e-6 [1] |
| 3 | +2.7e-4 [0], min K 0.37: **flagged (D2)** | −2.248e-5 [1] |
| 4 | −8.87e-7 [0], 6/12 differ 3e-7 | +2.8e-4 [0], min K 0.40: **flagged (D2)** |
| 6 | −2.060e-8 | −6.646e-8 |
| 8 | −3.467e-9 | −6.725e-9 |
| 12 | −1.703e-10 | −3.048e-10 |
| 16 | – | −2.012e-12 |

Readings (development, decides nothing):
- **Three regions in p.** Small p (a cut locus, three rays at the symmetric phase): small values, falling by a factor
  50–70 from ε = 12 to 16. Around p_c (where the symmetric 3-ray region ends): values of order ±3e-4, unconverged or
  flagged. Large p (one ray): converged, falling steeply with p, growing slowly with ε.
- **Φ(ε) = ∫dp⟨ŝ₄A₀⁻²⟩_θ is dominated by the region around p_c.** Its values are 10–1000 times larger than anything
  else in the scan.
- **ε = 12, p = 2 missed a cut locus.** The symmetric 3-ray region there is only |δ| < 0.02 wide, narrower than the
  128-phase grid step (0.049); D2a (§6.1) bisects for p_c directly. p_c = 2.44905 (ε 12), 3.72510 (ε 16), 4.38092 (ε 20).

## 6. D2: conjugate points

Near a phase where K → 0 on the first-arrival ray, the series coefficients grow like powers of 1/K and the radius of
convergence in c shrinks (E1: ε = 16, p = 4, K = 0.46 gave a value four orders above its neighbours). The window in θ
where the series fails shrinks too. Two options, to be decided after D1:
- **a bound:** show that the true θ-averaged contribution of the window (computed with the exact needle at finite q)
  is small against the rest, and exclude it;
- **a uniform form:** expand around the fold (the cubic normal form of Φ in ψ, an Airy-type needle) and integrate the
  window analytically.

### 6.1 D2a: the map of the cusp at p_c (25 Nov; diag/f68i_d2a.py; test run diag/out/f68i_d2a_20261125T082628Z)

p_c by bisection on the ray count at the symmetric phase θ_c = −p/2 (mod 2π). At p = p_c ∓ Δ: the 3-ray half-width w,
the window slope α, min K, and the smooth part on a mesh graded towards θ_c (first step 1e-4, ratio 2, 32 subintervals,
Gauss–Legendre 6 and 12). The window term only where the window fits inside the 3-ray region.

Test run, ε = 12 (Claude's container, 2 workers, 4356 s):

| | Δ = 0.2 | Δ = 0.05 | slope (2 points) |
|---|---|---|---|
| smooth part, below p_c (6/12 differ) | −1.786e-2 (2e-11) | −2.2046 (2e-9) | Δ^−3.47 |
| smooth part, above p_c | −3.808e-3 (9e-11) | −0.5762 (3e-9) | Δ^−3.62 |
| min K, below / above | 0.261 / 0.134 | 0.067 / 0.033 | Δ^0.99 / Δ^1.00 |
| w (below) | 4.556e-3 | 5.945e-4 | Δ^1.47 |
| α (below) | 5.950 | 2.997 | Δ^0.49 |
| third-ray gap (below) | 8.806e-3 | 5.789e-4 | Δ^1.96 |
| window term (below) | +8.56e-3 | not uniform (skipped) | |

The full run (ε = 12, 16, 20; Δ = 0.4, 0.2, 0.1, 0.05, 0.025) is on the author's PC.

Readings: the graded mesh converges at the cusp (6/12 to ≤ 3e-9). **The q⁶ coefficient is not integrable over p at
p_c**: ∫dp |p − p_c|^−3.5 diverges. The series in q² is not uniform near the cusp, and Φ(ε) as defined in §1 does not
exist. The exponents are those of the cusp catastrophe (§6.2).

### 6.2 The cusp normal form (25 Nov; diag/f68i_d2b_nf.py)

**The form.** At a cusp of the first arrival three rays merge (pitchfork at the symmetric phase). Let a be the
amplitude of the soft Jacobi mode. If the path through (s, ψ) is the ray plus aJ(s), then Φ(s, ψ) = V(ψ/J(s)), with

    V(a) = λa⁴/4 + σa²/2 + τa,   σ ∝ (p − p_c),   τ ∝ (θ − θ_c).

With λ = 1 (this fixes the unit of a), the needle factorises:

    W(c) = [∫ds √(s(1−s)) J(s)] · w(c; σ, τ),   w = ∫da [c − (V − V_min)]₊^{3/2},   A = N′ w/c².

Predictions of the form, all against D2a:
- min K ∝ σ ∝ Δ;
- third-ray gap σ²/4 ∝ Δ²;
- tie-gap slope α = 2τ′√|σ| ∝ Δ^½;
- 3-ray half-width w = 2(|σ|/3)^{3/2}/τ′ ∝ Δ^{3/2};
- q⁶ coefficient: A₃A₀⁻³ ∝ σ^{−5} per phase, over a θ-width ∝ σ^{3/2}, so ∝ Δ^{−3.5}.

**Quantitative check against the D2a geometry (ε = 12).** σ comes from the third-ray gap: −0.1877 (Δ 0.2) and
−0.04812 (Δ 0.05), i.e. σ′ = 0.94 and 0.96 per unit p. Then τ′ = α/(2√|σ|) = 6.867 and 6.830 rad⁻¹. The predicted
half-width w = 2(|σ|/3)^{3/2}/τ′ is 4.559e-3 and 5.95e-4 against the measured 4.556e-3 and 5.945e-4.

**Check against the exact needle (ε = 12, θ = θ_c; diag/f68f_1b_exact, scratch probe, 25 Nov).**
- At p = p_c the local slope d ln A/d ln q is −0.478, −0.483, …, −0.498 for q from 0.5 down to 0.031. The form
  predicts A = N′ 4^{1/4} B c^{−1/4} ∝ q^{−1/2}, with B = ∫(1 − x⁴)^{3/2}dx; this gives N′ ≈ 1.194.
- With N′ from the cusp and σ from the ray gaps (no fit to these curves), the form gives A(q) at Δ = −0.05, +0.05,
  +0.2 to ≤ 0.3 % for q ≤ 0.09, and to 2–3 % at q = 0.5. This includes the non-monotone A(q) of the double well below
  p_c. (σ above p_c is taken as σ′Δ with σ′ from below.)

**The age term.** A = c^{−1/4} F(σc^{−1/2}, τc^{−3/4}) holds exactly for the form. So the age term over the unfolding
plane,

    I(γ) = ∫dσ dτ J(σ, τ; γ),   J = ∫₀^∞ds g₀(s² A(γs; σ, τ)),

is homogeneous: I(γ) = C γ^{11/7} (s ~ γ^{1/7}, c ~ γ^{8/7}, σ ~ c^{1/2}, τ ~ c^{3/4}).
- Convergence: at large σ, J ∝ γ³σ^{−5}; over τ this gives σ^{−3.5}.
- The BD moments μ₀..μ₂ remove integer powers below γ³, not γ^{11/7}. In the exchanged order,
  C = (4/7) G(9/7) FP∫F^{−9/7}, with G(z) = ∫v^{z−1}g₀(v)dv and G(9/7) = −0.05771 (G(½) = G(1) = G(3/2) = 0, G(2) = 1).
- **If C ≠ 0, every cusp of the first arrival adds a term of order γ^{11/7} to the age term, which dominates the smooth
  γ³ for γ = k²/(2l₂) → 0.** Per unit p and θ the term is C (2π σ′τ′)⁻¹ N′^{−1/2} (γ/√N′)^{11/7}. Counted like §1
  (T = κ/k, dΩ = 2π dp/κ), it gives dE/dn ∝ k^{15/7} at fixed ε, instead of k⁵.
- The singularity index 1/4 of A is Berry's index of the cusp in catastrophe optics. The first arrival is a global
  minimum, so it cannot sit on a fold: the cusps (the end points of the cut-locus segments) are its only generic
  degenerations in the (p, θ) plane.

**C (diag/out/f68i_d2b_nf_20261125T101212Z; Claude's container, 2 workers, 435 s).**
- Checks of the form:
  - w against the closed form at the cusp: 1.4e-15;
  - paraxial limit: 1e-9;
  - adaptive quadrature of tilted double wells: 1.5e-11;
  - J at 48 against 96 s-nodes: ≤ 5e-6.
- I(γ = 1) = −8.82941e-3. The σ-nodes 12 → 20 and τ-nodes 20 → 32 change it by 2e-6 relative.
- Tails:
  - the σ-tails beyond |σ| = 16 (power fits p = 3.65 and 3.38, against the predicted 3.5) add −1.9e-5;
  - the τ-range beyond t = τ/|σ|^{3/2} = 10 is bounded by 1.1e-5 (1.2e-3 relative).
- I(γ = 1/4) = −9.99968e-4 gives C = −8.8324e-3. The exponent from the two γ is 1.57118 (11/7 = 1.57143); the grids
  are fixed while the structures move, so this checks the quadrature.

    **C = −8.83e-3 (±0.1 %): the cusp term does not vanish.**

- J changes sign over the plane and its window part cancels to 1e-5 of its amplitude, so C is a small net number.
  It is still 10⁻²-sized against a γ³ smooth part.
- An earlier run of this file integrated τ up to 50|σ|^{3/2}. There J is at its rounding floor (~1e-11) and the
  rows for |σ| > 15 were noise; that run was discarded, not kept.

**What this means (development, decides nothing).**
- In the age term of an old field (ε > π²), each cusp of the first arrival contributes
  C (2πσ′τ′)⁻¹ N′^{−1/2} (γ/√N′)^{11/7} per (p, θ).
  - This is the leading term for small γ = k²/(2l₂), ahead of the smooth γ³ (q⁶) term of §1.
  - At fixed ε it gives dE/dn ∝ k^{15/7} instead of k⁵.
- The quantity F6.8i asks about (Φ(ε) and ∫ε^{−1/2}Φ dε) must therefore be redefined around the cusp amplitude
  Ψ(ε) = Σ_cusps (σ′τ′)⁻¹ N′^{−1/2−11/14}. **This is the author's decision.**
- Open before any reading:
  1. D2a full run: the exponents over five Δ at ε 12, 16, 20.
  2. σ′, τ′, N′ along ε, and the number of cusps per ε (the end points of every cut-locus segment).
  3. The direct check: the age term from the exact needle on a (p, θ) patch around p_c at several γ, against
     C γ^{11/7} plus the smooth part.
  4. The young-field side: below the first caustic there is no cusp, and F6.8f stands.

### 6.3 D2a full run and the cusp constants (25 Nov; author's PC, 6 workers, 18520 s; diag/out/f68i_d2a_20261125T145244Z)

ε = 12, 16, 20; Δ = 0.4, 0.2, 0.1, 0.05, 0.025 on both sides of p_c. The smooth part converges at every point (6/12
differ ≤ 8e-8 relative). p_c = 2.449046 (ε 12), 3.725109 (ε 16), 4.380924 (ε 20).

**Local slopes between neighbouring Δ (from 0.4 down to 0.025).** The cusp predicts −3.5 for the smooth part and 1.5
for w.

| ε | smooth part, below | smooth part, above | w (below) |
|---|---|---|---|
| 12 | −3.415, −3.464, −3.484, −3.492 | −3.811, −3.661, −3.581, −3.540 | 1.414, 1.459, 1.480, 1.490 |
| 16 | −3.536, −3.515, −3.502, −3.489 | −3.729, −3.627, −3.572, −3.548 | 1.446, 1.472, 1.483, 1.487 |
| 20 | −3.553, −3.516, −3.492, −3.466 | −3.663, −3.601, −3.569, −3.566 | 1.448, 1.470, 1.478, 1.476 |

Both sides approach −3.5, from opposite directions. min K ∝ Δ^{0.98–1.00} and α ∝ Δ^{0.49–0.51} on every ε.

**The normal form against the geometry, all 14 points with a window.** σ = −2√(third gap) and τ′ = α/(2√|σ|). The
predicted half-width 2(|σ|/3)^{3/2}/τ′ agrees with the measured w to 1.5e-4 relative or better (to 1e-3 at
Δ = 0.025, where w is 2e-4).

**Cusp constants** (λ = 1). σ′ and τ′ are extrapolated linearly to Δ → 0 from the two smallest Δ. N′ is from the
exact needle at p_c and θ_c (scratch probe, 25 Nov), taken at q = 0.0625. N′ still drifts by ≤ 0.6 % down to q = 0.031,
and at ε 16 and 20 the local slope overshoots −0.5 slightly (−0.503, −0.511). That is the resolution of the exact needle
at the smallest q, and is open.

| ε | σ′ | τ′ | N′ | Ψ = (σ′τ′)⁻¹ N′^{−9/7} | per cusp C Ψ/(2π) |
|---|---|---|---|---|---|
| 12 | 0.970 | 6.815 | 1.196 | 0.1202 | −1.69e-4 |
| 16 | 1.513 | 6.383 | 1.267 | 0.0764 | −1.07e-4 |
| 20 | 1.811 | 6.026 | 1.337 | 0.0631 | −8.9e-5 |

So per cusp, per unit p and θ/2π, the age term is C Ψ (2π)⁻¹ γ^{11/7}. Over ε 12–20, Ψ ∝ ε^{−1.26}.

Open for the static question:
- Ψ near ε → π⁺², where the cut locus is born at p = 0 and its two cusps (±p_c) merge.
- The number of cusps at larger ε (new cut loci after the next caustics).
- The direct check of γ^{11/7} from the exact needle.

### 6.4 D2c: the cusp constants along ε (25 Nov; diag/f68i_d2c.py; Claude's container, 2 workers)

Records (three batches; one run of the whole grid exceeds the container's 2 h limit):
diag/out/f68i_d2c_20261125T211315Z (ε 9.9–12, 2879 s), f68i_d2c_20261125T214143Z (ε 14–25, 1702 s),
f68i_d2c_20261125T225751Z (ε 30–39, 4563 s). Test record f68i_d2c_20261125T175202Z.

For each ε, ray counts on both symmetric lines over p ∈ [0.02, 12]; p_c by bisection; σ′, τ′ from the geometry at
p_c − Δ (Δ = D, D/2, D = min(0.05, p_c/20)), extrapolated to Δ → 0; N′ from the exact needle at q = 0.0625 (§6.3).

**Number of cusps.** On every ε in 9.9–39 the line θ = −p/2 has 3 rays for p < p_c and 1 above; the line −p/2 + π
has 1 ray everywhere. So there is one cut locus, ending in one cusp pair ±p_c, up to ε = 39. The cut-locus search at
p_c/2 (f68f_1b_cut, 128 phases) reports extra loci off the symmetric line at ε 35 (18) and ε 39 (39). They are not
real: all but one have α = 0, zero gap and equal A on both sides (one ray counted twice at a grid point), the remaining
one at ε 39 has a gap of 48.8 (no tie), and a direct ray count at five of the phases (ε 39) gives 1 ray. They sit on
the 1-ray side, θ − θ_c ∈ [1.8, 3.1] in absolute value, and appear only as ε approaches 4π²; the search code is shared
with F6.8f and is left unchanged.

**The normal form.** The predicted half-width 2(|σ|/3)^{3/2}/τ′ agrees with the measured w to ≤ 2e-4 at every ε from
10 to 39; at ε 9.9 (w ≈ 2–5e-6) to 0.4 % and 0.9 %, the resolution of the bisection there.

| ε | ε − π² | p_c | σ′ | τ′ | N′ | slopes ln A / ln q | Ψ |
|---|---|---|---|---|---|---|---|
| 9.9 | 0.030 | 0.31332 | 0.1293 | 7.0857 | 1.1581 | −0.497, −0.498 | 0.9036 |
| 10 | 0.130 | 0.64718 | 0.2547 | 7.0724 | 1.1597 | −0.497, −0.498 | 0.4589 |
| 10.25 | 0.380 | 1.09567 | 0.4325 | 7.0393 | 1.1647 | −0.497, −0.498 | 0.2700 |
| 10.5 | 0.630 | 1.39846 | 0.5507 | 7.0068 | 1.1690 | −0.496, −0.497 | 0.2120 |
| 11 | 1.130 | 1.84139 | 0.7326 | 6.9433 | 1.1796 | −0.498, −0.499 | 0.1590 |
| 12 | 2.130 | 2.44905 | 0.9703 | 6.8218 | 1.1963 | −0.497, −0.498 | 0.1200 |
| 14 | 4.130 | 3.21902 | 1.2863 | 6.5979 | 1.2310 | −0.497, −0.498 | 0.0902 |
| 16 | 6.130 | 3.72511 | 1.5134 | 6.3943 | 1.2670 | −0.501, −0.503 | 0.0762 |
| 20 | 10.13 | 4.38092 | 1.8328 | 6.0321 | 1.3365 | −0.507, −0.511 | 0.0623 |
| 25 | 15.13 | 4.87923 | 2.1011 | 5.6401 | 1.4186 | −0.514, −0.520 | 0.0538 |
| 30 | 20.13 | 5.20186 | 2.2662 | 5.2961 | 1.4899 | −0.513, −0.518 | 0.0499 |
| 35 | 25.13 | 5.42694 | 2.3666 | 4.9884 | 1.5527 | −0.507, −0.509 | 0.0481 |
| 39 | 29.13 | 5.56039 | 2.7047 | 4.7664 | 1.6772 | −0.551, −0.566 | 0.0399 |

At ε 12 and 16 this agrees with §6.3 to 0.3 %; at ε 20 σ′ is 1.2 % larger (1.833 against 1.811), Ψ 1.3 % smaller.

**The birth of the cusp (ε → π²⁺).** p_c²/(ε − π²) = 3.23, 3.21, 3.16, 3.10, 3.00 at ε 9.9–11 (fit p_c² = 3.234 e
− 0.207 e², e = ε − π², on e < 1.2): the cut locus and its cusps are born as a pitchfork at p = 0. σ′/p_c stays at
0.394–0.398 (0.413 at ε 9.9), τ′ → 7.09 and N′ → 1.156 are finite. So Ψ ∝ 1/σ′ ∝ (ε − π²)^{−1/2}:

Ψ √(ε − π²) = 0.1657, 0.1665, 0.1683, 0.1690, 0.1751 at ε 10–12 (0.1575 at ε 9.9),

and Ψ ≈ 0.165 (ε − π²)^{−1/2} at the birth. This singularity is integrable in ε. The ε 9.9 point is 5 % lower; it is
the point with the smallest window (above) and is not used for the limit. Whether the limit is exactly the ε 10–11
value is not checked below ε 9.9.

**Large ε.** Ψ falls from 0.120 (ε 12) to 0.048 (ε 35); Ψ √(ε − π²) rises slowly (0.175 → 0.241), so Ψ falls more
slowly than (ε − π²)^{−1/2} there. The local slope of ln A against ln q at the cusp drifts away from −1/2 with ε
(−0.51 to −0.52 at ε 20–30), so N′ at q = 0.0625 is uncertain by 1–2 % there. At ε 39 the drift is large (−0.55,
−0.57; N′ 1.677 → 1.748 from q 0.0625 to 0.031) and σ/Δ changes by 3.7 % between the two Δ: near the next caustic
(4π² = 39.48) the cusp regime is narrower than q = 0.03, and Ψ(39) = 0.040 is not a cusp constant. It is not used.

**What this gives for the static question.** Per cusp the age term is C Ψ(ε) (2π)⁻¹ γ^{11/7} per unit p and θ/2π,
with one cusp pair for π² < ε < 4π² (up to ε 35 checked; ε 39 is not resolved at these q), and Ψ finite on
(π², 4π²) with an integrable (ε − π²)^{−1/2} singularity at the birth. Open:
- the direct check of γ^{11/7} from the exact needle on a (p, θ) patch around a cusp;
- ε ≥ 4π² (new caustics, new cut loci);
- how F6.8i is reformulated with this term (the author's decision).

### 6.5 Correction: p_c was too low in §6.3 and §6.4 (26 Nov; diag/f68i_d2c.py revisions 2–2.2)

Records: test diag/out/f68i_d2c_20261126T074303Z (ε 9.9, 20; reproduced on the author's PC,
f68i_d2c_20261126T092239Z, identical in every printed digit), f68i_d2c_20261126T085437Z (ε 10–14),
f68i_d2c_20261126T111759Z (ε 16, 25, 30), f68i_d2c_20261126T103445Z (ε 35, 39; reproduced on the author's PC,
f68i_d2c_20261126T122048Z, identical in every printed digit); window check
f68i_d2c_bracket_20261126T100536Z.

**The error.** The ray count of f68f_1b_cut.rays uses a grid of B.N_V = 400 slopes v. Just below the cusp the three
rays are closer than one grid cell, so the count falls from 3 to 1 before p_c, and the bisection placed p_c too low.
σ, τ′ and the needle were then measured at points shifted from the cusp. The shift is the reason for
- the drift of σ/Δ with Δ in the D2a run (§6.3);
- the 1–2 % spread between σ′ of §6.3 and §6.4 at ε 20;
- the overshoot of the needle slopes past −1/2 (−0.51 to −0.57 at ε ≥ 16), which §6.3 put down to the needle's
  resolution and §6.4 to a narrowing cusp regime;
- the 5 % low Ψ √(ε − π²) at ε 9.9.

The evidence, before revision 2: a fit σ = σ′(Δ + d) + b(Δ + d)² with a free origin d to the D2a rows below the cusp
gives d = 7.9e-7 (ε 12), 1.58e-4 (ε 16), 4.0e-4 (ε 20). With the origin moved by d, the local slopes of the smooth part
(scratch, from the D2a record) approach −3.5 cleanly. At ε 20 they are −3.558, −3.526, −3.513, −3.506 below the cusp
and −3.658, −3.590, −3.548, −3.525 above it; at ε 16, −3.538, −3.519, −3.510, −3.505 below and −3.726, −3.623, −3.563,
−3.532 above. min K approaches 1 on both sides. So the D2a conclusion of §6.3 (the cusp normal form, the −3.5
exponent) stands, and its p_c values and constants are replaced by the table below.

**Revision 2.** B.N_V = 6400 in every job. σ is measured at p_c − Δ for Δ = 4D, 2D, D, D/2. The same fit with a free
origin is recorded alongside the two-point extrapolation; its d is now ≤ 6.2e-6 on every ε (it was up to 4e-4). The
predicted half-width agrees with the measured w to ≤ 1.4e-3 at every point with ε ≥ 10, and to 0.5 % at ε 9.9 (w ≈ 2e-6).

**Revision 2.2.** At ε 35 and 39 the exact envelope (f68f_1b_exact.phase_A, shared with F6.8f, unchanged) refuses
the needle window at the true p_c. The cause is not another valley. In the middle nodes (s ≈ 0.5) the valley of Φ at the
cusp is flat, with Φ < q²/2 over |ψ − ψ_ray| < 0.57, while the window can be widened only to 4W = 0.575, where W is the
paraxial width. A refused needle is retried with guide rays ψ_ray ± 3W(s), which only move the window edges outward.
Where both windows work (ε 20, 30) they agree to 2e-9 in N′.

**Corrected constants** (σ′ two-point; free-origin fit and its d in brackets; N′ at q = 0.0625):

| ε | ε − π² | p_c | σ′ [fit, d] | τ′ | N′ | slopes ln A / ln q | Ψ |
|---|---|---|---|---|---|---|---|
| 9.9 | 0.030 | 0.3135853 | 0.1232 [0.1232, +9.5e-7] | 7.0857 | 1.1577 | −0.496, −0.497 | 0.94902 |
| 10 | 0.130 | 0.647195 | 0.2543 [0.2543, −1.8e-6] | 7.0724 | 1.1596 | | 0.45963 |
| 10.25 | 0.380 | 1.095731 | 0.4310 [0.4311, −6.2e-6] | 7.0393 | 1.1644 | | 0.27104 |
| 10.5 | 0.630 | 1.398455 | 0.5507 [0.5508, −4.6e-6] | 7.0068 | 1.1691 | | 0.21201 |
| 11 | 1.130 | 1.841530 | 0.7266 [0.7268, −4.9e-6] | 6.9432 | 1.1783 | −0.496, −0.497 | 0.16051 |
| 12 | 2.130 | 2.449046 | 0.9703 [0.9705, −4.2e-6] | 6.8218 | 1.1963 | | 0.11998 |
| 14 | 4.130 | 3.219028 | 1.2856 [1.2858, −2.9e-6] | 6.5979 | 1.2308 | | 0.09027 |
| 16 | 6.130 | 3.7252667 | 1.4993 [1.4995, −2.5e-6] | 6.3942 | 1.2638 | −0.498, −0.498 | 0.07720 |
| 20 | 10.13 | 4.3813227 | 1.7901 [1.7901, −2.5e-7] | 6.0318 | 1.3264 | −0.498, −0.499 | 0.06441 |
| 25 | 15.13 | 4.8798131 | 2.0299 [2.0299, −9.0e-7] | 5.6396 | 1.4008 | −0.499, −0.499 | 0.05663 |
| 30 | 20.13 | 5.2023443 | 2.2021 [2.2020, +4.1e-7] | 5.2956 | 1.4731 | −0.500, −0.500 | 0.05211 |
| 35 | 25.13 | 5.4271580 | 2.3364 [2.3360, +2.9e-6] | 4.9881 | 1.5444 | −0.501, −0.500 | 0.04907 |
| 39 | 29.13 | 5.5623074 | 2.4262 [2.4255, +4.6e-6] | 4.7633 | 1.6012 | −0.501, −0.501 | 0.04724 |

(Slopes not listed are in the records; they lie in −0.496 to −0.498.) Compared with §6.4: p_c is higher by up to
6e-4. The largest changes are σ′ −10 % and Ψ +18 % at ε 39, and σ′ −2.4 % and Ψ +3.4 % at ε 20. At ε 12 the change
is below 0.1 %. N′ no longer drifts with q at any ε: the slopes are −1/2 within 0.004.

Per cusp, C Ψ/(2π) = −1.69e-4 (ε 12), −1.08e-4 (ε 16), −9.05e-5 (ε 20). Over ε 12–20, Ψ ∝ ε^{−1.22}.

**The birth (ε → π²⁺), corrected.** On e = ε − π² < 1.2:
- p_c²/e = 3.240 − 0.214 e;
- σ′/p_c = 0.3928 + 0.0016 e;
- τ′ = 7.089 − 0.130 e;
- N′ = 1.157 + 0.019 e;
- Ψ √e = 0.1653 + 0.0047 e.

ε 9.9 now lies on the same line (Ψ √e = 0.1655). So Ψ ≈ 0.165 (ε − π²)^{−1/2} at the birth, as §6.4 said, and the
ε 9.9 anomaly of §6.4 was the p_c error.

**Large ε, corrected.** Ψ √(ε − π²) rises smoothly: 0.191 (ε 16), 0.205, 0.220, 0.234, 0.246, 0.255 (ε 39). The cusp
regime holds at q ≥ 0.031 up to ε 39. The statement in §6.4 that Ψ(39) is not a cusp constant is withdrawn: its slopes
are −0.501, and σ/Δ changes by 0.04 % between Δ = D and D/2 (3.7 % in §6.4) and by 0.3 % from 4D to D/2.

**Superseded:**
- §6.3: the values of p_c at ε 16 and 20, and the constants table. Its slope tables of the smooth part and of w,
  and the normal-form check, are not affected.
- §6.3: the reading of the N′ drift as needle resolution.
- §6.4: the constants table.
- §6.4: the ε 9.9 remark.
- §6.4: the ε ≥ 20 drift paragraph and the exclusion of ε 39.
- §6.4: the parenthesis "ε 39 is not resolved at these q".

The cusp count and the artefact loci of §6.4 are not affected.

**What this gives for the static question.** Per cusp, the age term is C Ψ(ε)(2π)⁻¹ γ^{11/7} per unit p and θ/2π.
There is one cusp pair on π² < ε < 4π², checked at 13 values of ε from 9.9 to 39. Ψ is finite on (π², 4π²), with an
integrable (ε − π²)^{−1/2} singularity at the birth, and Ψ √(ε − π²) rises from 0.165 to 0.255. Open:
- the direct check of γ^{11/7} from the exact needle on a (p, θ) patch;
- ε ≥ 4π²;
- the reformulation of F6.8i (the author's decision).

## 7. The author's decision on the cusp term (26 Nov, 13:47)

- **F6.8i stays as it is (option C).** It keeps its definition of Φ(ε) (the γ³ coefficient), its criteria and its
  weight ε^{−1/2}. It measures the smooth part. Because the cusp term may dominate the asymptotics, a result of F6.8i
  alone is not a verdict on whether an old field has a finite static response.
- **A new card with two parts (option B), frozen only after the direct check below:**
  E(γ, ε) = C Ψ(ε)(2π)⁻¹ γ^{11/7} + Φ(ε) γ³ + …, where the first term is the cusp term and the second the smooth part.
  - The card reports Ψ(ε), Φ(ε) and the sum separately.
  - It allows further terms if the direct computation finds them.
  - The ε-weight of the γ^{11/7} term must be derived from the definition of the total response. It must not be taken
    over from ε^{−1/2} by substituting Ψ for Φ.
- **The direct check comes first (development, approved):**
  1. The exact needle around a real cusp, without the normal form as a substitute.
  2. The scaling in γ, against γ^{11/7} and against γ³. C is not fitted; the exponent and the predicted coefficient
     must both agree.
  3. A control without a cusp: the same computation must give γ³.
  4. The birth of the cusp: points closer to ε = π², so that the integrable singularity of Ψ does not rest on one
     numerical point.
  5. The cusp count: the cusp pair is counted exactly once, with its multiplicity in the convention of Φ.
- **Not decided:** whether the static response of an old field converges. π² < ε < 4π² covers only the first cusp
  pair. A verdict needs:
  - the correct weight of the cusp term;
  - Ψ(ε) on this interval;
  - the same computation beyond the next caustics;
  - a check of whether the number of cusps or their contributions change with ε in a way that affects convergence.

  No extrapolation is made to ε ≥ 4π².
