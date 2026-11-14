# Card C1 (F5.3): the field of the C0 instrument on a curved 3+1 order

**14 Nov 2026, 19:00. Frozen by decision 29**, under the author's instruction of 18:43 ("as soon as the test is
ready, you can start it straight away").

**Test F5.3**: phases/F5_gravity/F5_03_field_curved_order_3p1.py, md5 4ff14c8a75c4dba678b999a669d287bd.

**Helpers:**
- diag/c1_conf.py 7f05b01eefdbef0ce2899846fe8ae4eb
- diag/c1_apex.py e12a8e24b2983f522abf4ae3b53d31e6
- diag/c0_fast.py 1dc98997a7e414a7d7eeb0b7a8b6a5da
- those of F3.16 / F3.16b (decisions 27, 28) for diag/c0_apex.py, diag/c0_minvar.py and diag/c0_links4.py.

The md5 sums are of the LF files.

## 1. Hypothesis

On a curved 3+1 order, the field of the C0 rule (MV8 + counterterm, order-only) follows the exact mean of the rule.
- That mean is Mecke's formula for the inhomogeneous sprinkling, with the curved volume of every interval.
- With rising density it converges to a Green function in which a near-null pair weighs ρ Ω_u⁴ ⟨Ω⁴⟩_I^(−1/2): the
  curvature along the whole null segment enters.
- This is not the local conformally coupled field (nor any local ξR; development).

The test claims both:
- the agreement with the exact mean (P1);
- the departure from the local model with the same finite-density smearing, in the predicted direction (P2).

## 2. Model

- **Metric**: g = Ω(x)² η, Ω = 1 − 0.15 e^(−x²/(2·1.5²)).
- **Domain**: the cones are flat, so the pyramid |x|, |y|, |z| ≤ 7 − t (T 6) is past-complete.
- **Pairs**: a flat sprinkling of proper density ρ and its curved child, thinned with probability Ω⁴. The seeds are
  derive_seed(base, "F5.3", ρ, k) for the sprinkling and derive_seed(base, "F5.3-thin", ρ, k) for the thinning.
- **Rule**: the C0 rule unchanged (MV8, c and δ at the proper ρ), linear, μ 1, data φ0 = A.
- **Runs**: ρ 30 and 60, 8 pairs each; base seed 0 decides.

## 3. Readouts and predictions

- **D** = ⟨φ⟩_curved − ⟨φ⟩_flat in bins of t and x (width 0.5) over the inner events, as mean and SE over the pairs.
- **P_exact**: diag/c1_apex.mean_field, mode 0.
  - It integrates the exact mean with the diamond-averaged ⟨Ω⁴⟩ (a Gauss product rule over the unit diamond).
  - Its flat case equals the 1D (A3) of C0 to 3e-4.
  - It is converged in the grid (dt 0.05 against 0.025: identical to 1e-3).
- **P_local**: mode 1, the conformally coupled Green function Ω_u⁻¹ Ω_v⁻¹ G_η smeared by the same flat kernel. At
  ρ 1 000 it follows the PDE with ξ = 1/6 within 0.016.
- **The flat parents** against the flat (A3), per time bin.

## 4. Criteria (the form of decision 28)

All tests use Student t with n − 1 degrees of freedom and Bonferroni family-wise α 1 % over all items of both
densities.

- **INSTRUMENT**: the flat parents per time bin (t ≥ 0.5, |x| ≤ 2.5): |d| / √(SE² + 0.002²) < t_crit.
- **P1, the exact mean**: every bin (t ≥ 1, |x| ≤ 2.5, both densities) has |D − P_exact| / √(SE² + 0.005²) < t_crit.
  - The floor 0.005 covers the quadrature of the diamond average and the grid.
  - Any bin above it: **FAILED**.
- **P2, not the local field**: in the centre window (|x| < 0.75, t ≥ 4.5) at each density, Δ = the mean over the pairs
  of (D − P_local). PASSED side Δ − 2 SE > 0; FAILED side Δ + 2 SE < 0.
  - Predicted Δ (P_exact − P_local): +0.141 at ρ 30, +0.148 at ρ 60.
- **Verdict order**: INSTRUMENT_FAILED → FAILED → PASSED → UNDECIDED.

**On the development rows** (seeds 21000+ and 31000+, 8 pairs per density; the same readouts):
- P1: max z 3.48 (ρ 30) and 2.51 (ρ 60), against t_crit 8.70; χ²/bin 0.58 / 0.64.
- The local model scores χ²/bin 33.6 / 32.6, max z 17.7 / 15.5.
- P2: Δ = +0.102 ± 0.020 and +0.146 ± 0.010.

## 5. Reading declared now

**PASSED.** On a curved order the order-only field follows its exact mean, a Green function weighted by the curved
volume along each near-null pair. It departs, as predicted, from the local conformally coupled field with the same
smearing. The C0 instrument therefore does not carry a local covariant field equation on curved orders; whether it
does in the limit of a vanishing pair length is a separate question (C1 continues).

**FAILED.** Either some bin departs from the exact mean (the derivation, or the model of the substrate, is wrong), or
the field is closer to the local model than the exact mean says.

**Limits:**
- one conformally flat, static, planar well (a 0.15, L 1.5);
- linear field;
- T ≤ 6;
- two densities.

The theory's metric −Π²dt² + Π⁻²dx², with tilted cones, is not tested here (C1b).

## 6. Development rows and an open item

- **Development**: docs/findings.md, "C1 development rows"; diag/c1_conf.py, diag/c1_apex.py.
- **Open (not a criterion here)**: the flat parents of the development show −0.023 ± 0.008 (ρ 30) and −0.029 ± 0.009
  (ρ 60) against (A3) at t 4.75–5.25. Pooled with C0 (docs/findings.md, F3.16b) the last flat slabs sit at about −0.01,
  2–3 SE below the exact (A3). The derivation says the mean there is exact. The instrument level here allows this
  size: z ≈ 3.3 against t_crit ≈ 5.8.
