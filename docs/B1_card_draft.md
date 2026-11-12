# Card B1: inertia of a domain wall of a field of the causal order in 2+1 (intrinsic instrument)

**Revision 2, frozen 12 Nov 2026 (decision 21), implemented as F3.15**
(phases/F3_structure/F3_15_wall_inertia_2p1.py; seed 0 decides, seed 1 replicates). Revision 2 takes the proposals of
revision 1 unchanged, with the details the implementation needs (below and in the docstring of F3.15) and one change
after the smoke run: telegraph Γ = 1, 2 (decision 21). Development rows: docs/findings.md, Undecided, "B1 development
rows" (12 Nov). Structure taken from F3.14a/b (docs/H5_card_draft.md, section 8), transferred to 2+1. Decisions 16, 18,
20, 21.

## 1. Hypothesis

The retarded Green function of the 2+1 order, with the proper time read from the order alone (the longest chain and
the weight of decision 17), gives a domain wall of a field with inertia: it keeps its identity, moves at the launched
velocity, contracts by γ, carries the energy 8γ per unit length, and its inertial mass from a kick equals its rest
energy per unit length. In 1+1 this is F3.14b (replicated); B1 asks whether it survives the transverse dimension and
the intrinsic weight.

## 2. Model

- φ(v) = φ₀(v) − (1/ρ) Σ_{u ≺ v} w(L_uv) [μ sin φ(u) + f χ(t_u)], μ = 1; w(L) of decision 17 (variant B), L from the
  forward dynamic programme over the links (`diag/b0_chain_field.py`, `diag/b1_dev.py`).
- φ₀: the free (massless) wave of the initial data of the sine-Gordon kink extended along y, launched at v (declared
  coordinate launch, as in F3.14a).
- Domain: the pyramid of dependence |x| ≤ T + 8 − t, |y| ≤ T + 2.5 − t; no wall enters the dynamics.
- Arms on the same sprinkling: rest; launched v = 0.3 and 0.6; kick (at rest, f = 0.05 for t ∈ [1, 3)).

## 3. Readouts (coordinates declared, as in F3.14a)

On the band |y| < 2, per slab of width 1:
- the kink fit (centre x_c, width w) and Q from the plateaus 5–7 from the centre;
- the energy per unit length from local regressions in cells 0.5 × 1 (field1d.energy on the band), **with the noise
  bias of the squared estimates removed** (E[ĉ²] = c² + Var ĉ, E[cos(ĉ₀)] = cos c₀ e^{−Var/2}; new against F3.14a,
  because the per-event noise of the intrinsic field is ~0.3; the plain estimator is biased upward by tens of percent
  at low ρ). The same readout on the exact solution at the same events gives E_ex, w_ex, Q_ex.
- velocity: slope of x_c; p = γv; Δ(γv) of the kick paired with the rest arm, slopes over t ∈ [4, T − 1).
- **slab rule (revision 1):** the slabs t ∈ [0, T − 1) are read, the last slab is not, for every arm and every
  readout. Reason (development rows, not a result): at v = 0.6 the plateau window x_c + 7 leaves the pyramid
  |x| ≤ T + 8 − t in the last slab; the kink fit failed there in 7 of 8 sprinklings and the energy read 12.3 instead of
  10. The rule is geometric and is fixed for all arms, so it does not choose among results.

## 4. Criteria (frozen, revision 2)

Tolerance rule (10 Nov): expected 2 SE at the card's sample size + the largest development offset + a margin. Sample size:
**32 sprinklings** (section 5.2). In the column "derivation" the first number is 2 SE at n = 32 from the
development sd per sprinkling, the second the largest development offset (|mean| of the rows, either readout window).

| level | condition | proposed tolerance | derivation (2 SE₃₂ + offset) |
|---|---|---|---|
| instrument I1 | vacuum (B0 wave, 4 sprinklings): \|λ_eff\| + 2 SE | < 0.02 | B0: 0.0068 at ρ 25, 0.0038 at ρ 50 |
| instrument I2 | dispersion, k = (1,0), (0,1), (1,1), (2,1) separately, 4 sprinklings: \|ω/ω_c − 1\| + 2 SE | < 0.015 | B0: 1.004–1.007 at ρ 25, ≤ 1.002 at ρ 50 |
| instrument I3 | the rest wall at a second density ρ 60 (8 sprinklings): \|w/w_ex − 1\| + 2 SE | < 0.05 | ρ 30: +0.016 ± 0.002 |
| instrument I4 | forward against backward kernel; linear extension t + 0.3x; pyramid X against X − 2 | 1e-12 | B0: 3e-15 |
| rest | drift x_c(last read) − x_c(first) | 0.15 | 0.059 + 0.055 |
| rest | \|E − E_ex\|/8 (debiased, the same readout on the exact solution) | 0.05 | 0.003 + 0.019 |
| rest | \|w/w_ex − 1\| | 0.05 | 0.002 + 0.016 |
| rest | \|Q − Q_ex\| | 0.01 | 0.0002 + 0.0001 |
| inertia | \|v_meas − v\|, v = 0.3 and 0.6 separately | 0.03 (as 1+1) | 0.006 + 0.0065 |
| contraction | \|C − 1\|, each v | 0.05 (1+1: 0.1) | 0.0063 + 0.0062 |
| contraction | \|w/w_ex − 1\|, each v | 0.05 | 0.0064 + 0.014 |
| energy | \|e_B\| (debiased), each v | 0.05 (as 1+1) | 0.0087 + 0.030 (0.030 with the last slab; 0.013 under the slab rule) |
| masses | \|r_M − 1\| | 0.10 (as 1+1) | 0.030 + 0.052 |
| instrument | telegraph Γ = 1, 2, launched v 0.3 and kick, 8 sprinklings: (p_cons − p_Γ) − 2 SE > 0 and p_Γ ≤ 0.5 p_cons, p over t ∈ [4, 7) | as F3.14b | — |
| instrument | V = 0 (v 0.3, 8 sprinklings): D_L = L_cons − L_V0, L = E_loc(last read)/E_loc(first) | D_L − 2 SE > 0.2 | — |
| PASSED only | causality-blind (v 0.3, 8 sprinklings): the weights of each source kept, the receiving events drawn uniformly from all later events | D_L − 2 SE > 0.2 | — |

Each condition is evaluated as |·| + 2 SE (PASSED side) and |·| − 2 SE (FAILED side) against its tolerance.
**The binding condition is r_M:** the development value 0.948 ± 0.028 (per-sprinkling sd 0.085); if the offset −0.05
is real, at n = 32 |r_M − 1| + 2 SE ≈ 0.08 against 0.10. The offset is not resolved (1.8 SE); it comes from
Δ(γv) = 0.0816 ± 0.0025 against the continuum 0.0785 (+4 %, 1.2 SE) and from the rest energy +1.2 %. The readout window
of the kick slope moves r_M by 0.013 ([4, 8): 0.961; [4, 7): 0.948).

Verdict rule as F3.14a: INSTRUMENT_FAILED → FAILED (any inertia/mass condition with |·| − 2 SE above its tolerance)
→ PASSED (all hold) → UNDECIDED. Every condition on its own, nothing averaged, nothing chosen after the run.

## 5. Decisions (made: decision 21, 12 Nov)

1. ρ 30, T 8; I3 at ρ 60 (rest arm, 8 sprinklings).
2. 32 sprinklings × 4 arms.
3. Controls on 8 sprinklings: telegraph Γ = 1, 2 (changed from 2, 5 after the smoke run, decision 21); V = 0; the
   causality-blind sum with the weights kept and the receiving events randomised.
4. Band |y| < 2, Y = T + 2.5; the y-box dependence is not part of B1.
5. Slab rule: the last slab is not read.
6. I1 from the B0 vacuum wave run inside F3.15 (4 sprinklings).

Verdict order: INSTRUMENT_FAILED (any instrument row, including the rest rows, the telegraph and V = 0) → FAILED (any
of velocity, C, w/w_ex, e_B for either v, or r_M with |·| − 2 SE above the tolerance) → PASSED (all PASSED rows,
the blind control separated, Q of the moving arms) → UNDECIDED.

## 6. Reading declared now (draft)

**PASSED:** the inertia of a defect of a field of the order survives the transverse dimension and the intrinsic Green
weight: with the chain weight, the 2+1 order gives a wall with the inertial mass equal to its energy per unit length.
**FAILED:** the intrinsic weight or the transverse dimension breaks inertia at the tested scale.

**Limits:** a wall has no 1/r and no velocity along itself; the potential is chosen; the launch uses coordinates; the
noise of single events under the chain weight is ~1.5 × that of the coordinate reference (B0).

## 7. Development rows (12 Nov; decide nothing)

ρ 30, T 8, 8 sprinklings, variant B, slab rule applied; mean ± SE (sd per sprinkling):

| quantity | rest | v 0.3 | v 0.6 |
|---|---|---|---|
| v_meas − v | drift +0.016 ± 0.059 (sd 0.167) | +0.0019 ± 0.0058 (sd 0.0165) | −0.0019 ± 0.0050 (sd 0.0141) |
| C − 1 | — | −0.0062 ± 0.0051 | −0.0019 ± 0.0063 |
| w/w_ex − 1 | +0.0159 ± 0.0019 | +0.0096 ± 0.0040 | +0.0140 ± 0.0064 |
| E (debiased): E/8 − 1 or e_B | +0.0122 ± 0.0034 | +0.0132 ± 0.0037 | +0.0022 ± 0.0087 |
| plain (E − E_ex)/(8γ) | +0.074 | +0.067 | +0.039 |
| Q − Q_ex | −0.0001 ± 0.0002 | −0.0000 ± 0.0002 | −0.0001 ± 0.0002 |

Kick (f 0.05, τ 2): Δ(γv) = 0.0816 ± 0.0025 (sd 0.0072), continuum 0.0785; r_M = 0.948 ± 0.028.
The debiased energy read per slab rises with time (rest: 8.01 in the first slab, 8.2–8.5 in the last ones; its sd
grows from 0.04 to ~0.4): the energy of the field fluctuations in the window ±5 adds to the wall; the plain
estimator is still +7 % (the debiasing removes most of the bias but not all).
