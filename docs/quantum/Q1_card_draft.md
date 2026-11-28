# Q.1 — the sum-over-chains propagator on a 1+1 sprinkling (card, draft 0, 28 Nov 2026)

**Status: draft. Nothing here is frozen.** The author freezes the card. The first test of the quantum section
(decision 49; docs/quantum/README.md section 7). Development: diag/q1_dev.py, record
diag/out/q1_dev_20261128T093204Z.

## Question

A particle on a causal set propagates from x to y as a sum over chains x = z₀ ≺ z₁ ≺ … ≺ z_n = y. Each relation
carries an amplitude a and each intermediate element an amplitude b (S. Johnston, arXiv:0806.3083, section 3; 1+1:
the causal matrix C):

    K = aC + a²b C² + a³b² C³ + … = aC (1 − abC)⁻¹,   a = 1/2,  b = −m²/ρ.

For two fixed points at proper time τ and a Poisson sprinkling of density ρ, the expected number of chains with k
intermediate elements is (ρV)^k/(k!)² with V = τ²/2. So E K(x, y) = ½ J₀(mτ) **exactly at every density**: the
continuum retarded propagator of a free scalar of mass m in 1+1 (Johnston eq. 3.23).

The mean is therefore a code check. The question of the card is about **one order**:

> Does K(x, y) of a single sprinkling approach ½ J₀(mτ) as the density grows, at the Poisson rate N^{−1/2}
> (self-averaging)?

A quantum instrument that carries the continuum propagator only on average, but not on one realised order, could not
be used for interference or measurement on the order that actually exists.

## Model

- **Region:** the light-cone square u, v ∈ [0, 1] (u = t + x, v = t − x; dt dx = du dv/2), with x = (0, 0) and
  y = (1, 1) added as fixed elements (τ = 1, V = 1/2).
- **Sprinkling:** N ~ Poisson(N_exp) interior elements, ρ = 2 N_exp. Relation: z ≺ w iff u_z < u_w and v_z < v_w.
- **Propagator:** the forward recursion k_j = a(C_xj + b Σ_{i<j} k_i C_ij), in the order of u. This is K(x, z) for
  every z; K(x, y) is the last entry. Memory O(N), time O(N²).
- **Judged points:** m ∈ {2, 5, 10} (mτ; ½J₀ = +0.11195, −0.08880, −0.12297), N_exp ∈ {1000, 4000, 16000}
  (m²/ρ ≤ 0.05).
- **Sprinklings:** 2000 per point. Seeds derive_seed(base, "Q.1", m, N_exp, k), independent of the development seeds.

## Criterion (proposed; frozen only by the author)

For each m: s(m) = the least-squares slope of ln sd[K(x, y)] against ln N_exp over the three densities.

- **PASSED:** −0.65 ≤ s(m) ≤ −0.40 at all three m. One order approaches the continuum propagator at the Poisson rate.
- **FAILED:** s(m) > −0.25 at some m. No self-averaging: the instrument holds only on average.
- **UNDECIDED:** otherwise, including s < −0.65 (faster than Poisson, unexpected; explained first).

The slope's statistical error with 2000 sprinklings is ≈ 0.01 (relative error of an sd 1/√(2·2000), over a factor 16
in N). Development, 400 sprinklings: s = −0.50 / −0.53 / −0.50 over N_exp 1000–16000 (−0.47 … −0.68 between adjacent
densities from N_exp = 250).

## Instruments (INSTRUMENT_FAILED if any fails)

- **I1 (the recursion):** K(x, ·) from the recursion against the dense solve aC(1 − abC)⁻¹ on 20 orders of N = 50–300,
  max |difference| ≤ 1e-12.
- **I2 (the exact mean):** at every judged point, |mean K(x, y) − ½J₀(m)| ≤ 3.5 SE. Development: |z| ≤ 1.43 at all
  12 points.
- **I3 (the null, non-negative weights):** the same sum with b = +m²/ρ, at m = 5, N_exp = 4000. Its mean is ½I₀(mτ),
  positive and without zeros; |mean − ½I₀(5)| ≤ 3.5 SE.

## Reported, not judged

- sd[K(x, y)] and sd/|½J₀| per point (development at N_exp 16000: 2.1 % / 6.0 % / 7.2 % at m 2 / 5 / 10).
- The rms of K(x, z) − ½J₀(mτ_xz) over all interior z.
- The null's sd.

## Reading declared now

- **PASSED:** the sum-over-chains propagator is a usable quantum instrument on a single 1+1 order. Its error falls as
  N^{−1/2} at fixed mτ (m ≪ √ρ). This is a calibration; it has no physical conclusion of its own.
- It says nothing about:
  - 3+1 (Q.1b, links; the mean is not exact there);
  - interference (Q.2) or the Born rule (Q.3);
  - masses near the discreteness scale (m²/ρ → 1).
- **FAILED:** a single order does not carry the propagator. Q.2–Q.4 would then need an ensemble or another
  instrument, and that is decided before they are designed.
