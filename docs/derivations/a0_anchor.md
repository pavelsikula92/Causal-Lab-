# F6.7: the anchor a₀ = cH₀/2π – what the text of the theory claims, what a derivation would mean and where the risk is

An analysis on paper, 23 Oct 2026 (the author's assignment: "F6.7 and the Jacobson route"). No simulation. Source:
*Physics of Causality* (sections 1.6, 3.4, 7.1–7.5, 30, 33.1–33.2, 38.4–38.6), results of the laboratory and four
external data (at the end). None of it is a result of a test.

## 1. What the text claims

- a₀ ≡ cH₀/2π is the **anchor of the bridge to SI**: "inputs, not derived quantities" (7.2), "converts update rates into
  units of acceleration; the thresholds of the kernel do not depend on this scale" (7.3), "does not affect the
  behaviour of the kernel" (7.1). Motivation: H₀⁻¹ is the causal depth of the observable universe, c/H₀ the
  corresponding horizon, 2π "a complete phase closure" (7.3).
- Units: t₀ = H₀⁻¹, R₀ = c/H₀ (7.5), so accelerations are reported in units of cH₀.
- At the same time: in galaxies the substrate passes into **compression** (D_eff < 3), and this gives flat rotation
  curves and the BTFR without dark matter (3.4, 33). "The reorganisation occurs at a scale given by the global causal
  normalisation … because a₀ is cosmologically anchored, the transition is universal" (33.1). Classification
  χ_obs = ln(g/a₀): χ < 0 is compression (38.4).
- The plateau D_eff ≈ 3 arises from the saturation of the degree at a finite capacity K₀ (30.2).

## 2. A contradiction in the status of a₀

The text gives a₀ two roles that exclude each other:

1. **a reporting convention** (7.1–7.3): it does not change the dynamics, the thresholds of the kernel do not depend
   on it;
2. **a physical threshold of the transition** (33.1, 38.4): compression begins where g falls below a₀, and only this
   gives the flat curve and the BTFR.

If 1 holds, a₀ cannot fix where a physical transition occurs, and the claim about galaxies has no content. If 2 holds,
a₀ is a physical quantity of the kernel and "cH₀/2π" is **a claim about the dynamics**, not a convention. Then it can be
tested.

The laboratory can work meaningfully only with option 2. In it F6.7 translates into a single question:

> In units of the kernel (acceleration in units of cH₀) the transition 3 → 2 must occur at a*_kernel = 1/(2π).
> To derive F6.7 means to compute a*_kernel from the dynamics of the order, not to choose it.

The factor 2π is so far a motive ("phase closure"), not a computation.

## 3. Why the transition must occur at a fixed acceleration (and what follows)

Conservation of the flow in an effective dimension D gives g(r) ∝ M / r^(D−1). If the flow from a mass M spreads
three-dimensionally up to a radius r_t and two-dimensionally beyond (D_eff = 2, e.g. in the plane of a disc):

    r ≤ r_t:  g = GM/r²,        r > r_t:  g = (GM/r_t²) · (r_t/r),   v² = g r = GM/r_t   (a flat curve).

The circular velocity is flat, but its height depends on r_t. The BTFR (v⁴ = GMa₀, observed with a small scatter)
follows only if r_t is the radius where the Newtonian g falls to a fixed value:

    GM/r_t² = a₀  ⇒  r_t = √(GM/a₀)  ⇒  v⁴ = GMa₀.

If the transition occurred at a fixed length or density, v² ∝ M (or another slope) would follow, not v⁴ ∝ M. A sharp
transition gives a coefficient of exactly 1, a smooth transition another numerical factor. **The content of the theory
of galaxies is therefore one sentence: the dimension of the substrate changes at a fixed acceleration a₀.** The
dynamics must supply it (F6.6 with F6.4), not the bridge to SI.

## 4. The basis of the compression mechanism after the tests of the laboratory

The plateau D_eff ≈ 3 is defined in the text by saturation at a finite degree K₀ (30.2). The tests in fact chose
substrate A: a Lorentz-invariant causal set (F0.3 refuted K₀; axiom audit, Bombelli–Henson–Sorkin: a finite degree and
Lorentz invariance exclude each other). In a causal set embedded in a manifold the dimension (Myrheim–Meyer, growth of
balls) is the dimension of the manifold, hence constant. Compression D_eff < 3 therefore needs **a new basis in
substrate A**:
- either a dynamics of the order in which the effective dimension changes (F6.4: an ensemble of orders with a global
  action, F6.5),
- or an effective dimension of the flow (how many directions the flow actually uses), not the dimension of the
  geometry. That can be measured in the laboratory: F3.6 showed that a bundle with a carried direction merges into a
  smaller number of world lines.

Without one of them F6.6 has no mechanism and F6.7 has nothing to derive.

## 5. Time dependence: H₀, H(z), or Λ?

| Variant | a₀ today | a₀ at z = 2 | Note |
|---|---|---|---|
| a constant cH₀/2π (H₀ as a number) | 1.04–1.13·10⁻¹⁰ m s⁻² (H₀ = 67.4–73) | the same | only reporting; physically unjustified why today's H₀ |
| causal depth, a₀(z) = cH(z)/2π | the same | ≈ 3× (E(2) ≈ 3.0 for Ω_m = 0.315) | directly from the motivation in 7.3 (the depth was smaller earlier) |
| Λ: a₀ = c·c√(Λ/3)/2π | 0.86–0.94·10⁻¹⁰ (√Ω_Λ ≈ 0.83) | the same | constant; the theory would have to explain why Λ |

The observed RAR scale: g† = 1.20 ± 0.02 (stat.) ± 0.24 (syst.) ·10⁻¹⁰ m s⁻² (McGaugh, Lelli, Schombert 2016).
All three variants are within the systematic error today; today's value cannot separate them. The **evolution**
separates them.

The motivation of 7.3 ("causal depth") leads to the variant H(z). Rotation curves at z ≈ 1–2 (Genzel et al. 2017),
according to Milgrom's analysis, do not tolerate a much larger a₀ (at 4a₀ the fraction of "phantom" mass would come out
around 0.5, the data show baryon-dominated discs).

**The author's objection (23 Oct):** Milgrom assumes the same geometry everywhere, Causal Physics does not (D_eff
changes), and D_eff = 2 is not confirmed in the laboratory. The objection holds for the shape of the interpolation:
Milgrom's function is not our law and a direct comparison of the numbers does not hold. The direction of the effect,
however, does not depend on the mechanism. If compression begins at g = a₀(z) and a₀(z) is three times larger at z = 2,
compression begins at a larger radius and in a higher part of the galaxy, and the enhancement of the rotation above the
baryons is larger. The data show a smaller one. The tension therefore persists for every mechanism in which the
transition occurs at a₀(z) ∝ H(z). It disappears if the scale of the transition does not grow with the epoch. The
analysis of the data is, moreover, disputed (pressure support of the discs, the split of the mass between disc and
halo).

## 6. Conclusion for F6.7

1. The text of the theory must decide whether a₀ is a convention or a physical threshold. The laboratory can test only
   a threshold.
2. F6.7 is then a computation of the number a*_kernel = 1/(2π) in units of cH₀ from the dynamics, not a choice. The
   factor 2π is so far a motive.
3. The compression mechanism is built on K₀, which the tests excluded. It needs a new basis in substrate A (F6.4, or
   the effective dimension of the flow).
4. The risky prediction of the theory is the evolution of a₀ with the epoch. The variant H(z) follows from the
   motivation of causal depth and is in tension with the data at z ≈ 2 for every mechanism with a transition at a₀(z).
   The variant Λ agrees with the data, but needs its own justification. The decision is the author's.

## External data

- S. McGaugh, F. Lelli, J. Schombert, *Radial Acceleration Relation in Rotationally Supported Galaxies*,
  PRL 117, 201101 (2016), arXiv:1609.05917: g† = 1.20 ± 0.02 ± 0.24 ·10⁻¹⁰ m s⁻².
- M. Milgrom, *High-redshift rotation curves and MOND* (2017), arXiv:1703.06110.
- R. Genzel et al., Nature 543, 397 (2017): rotation curves at z ≈ 0.6–2.6.
- Planck 2018: H₀ = 67.4 km s⁻¹ Mpc⁻¹, Ω_m = 0.315; local measurements H₀ ≈ 73.
