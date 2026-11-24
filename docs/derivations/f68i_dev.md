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

## 6. D2: conjugate points (open)

Near a phase where K → 0 on the first-arrival ray, the series coefficients grow like powers of 1/K and the radius of
convergence in c shrinks (E1: ε = 16, p = 4, K = 0.46 gave a value four orders above its neighbours). The window in θ
where the series fails shrinks too. Two options, to be decided after D1:
- **a bound:** show that the true θ-averaged contribution of the window (computed with the exact needle at finite q)
  is small against the rest, and exclude it;
- **a uniform form:** expand around the fold (the cubic normal form of Φ in ψ, an Airy-type needle) and integrate the
  window analytically.
