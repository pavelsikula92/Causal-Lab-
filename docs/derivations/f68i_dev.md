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

### 3.1 Results

(pending)

## 4. Plan

| Step | What | Check |
|---|---|---|
| E1 | the smooth part of ⟨ŝ₄A₀⁻²⟩ beyond the caustic, coarse grid | trend only |
| D1 | the cut-locus window at q⁶: series in c around each cut locus (envelope theorem for S_b − S_a; q-series of A_a, A_b) | against the exact needle (f68f_1b_exact) at q = 1/2, 1/4, 1/8 on cut-locus phases |
| D2 | conjugate points: the window bound or the uniform (Airy/cusp) form | against the exact needle near the folds |
| D3 | Φ(ε) on ε up to ~40 (or until a trend is decided), p panels to p ~ 3√ε + 12 | phases 64 against 128, p panels, the b → 0 limit against F6.8f |
| Card | F6.8i: the reading of ∫ε^{−1/2}Φ (convergent: a static limit exists; divergent: none), with criteria fixed before the run | as K.18 / F6.8f |

**What F6.8i can decide.** Whether the optical static response of an old field converges with age (a static limit)
or keeps growing. If it converges, the old-field stiffness is finite and the far-field law and the short-range
prediction (tens of micrometres) can be computed. If it diverges, the static response of an old field is not a
property of the bare expected action, and the age has to be fixed by something else (decision 41: the balance of
action and entropy; the meaning of the age, F8.8).
