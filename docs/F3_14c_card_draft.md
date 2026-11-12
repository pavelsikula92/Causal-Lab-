# Card F3.14c: swerves of a field defect, the momentum diffusion of the shot noise of the order (1+1)

**FROZEN** by the author on 10 Nov 2026, 22:56 (revision 1). Test file `phases/F3_structure/F3_14c_defect_swerves.py`;
smoke run OK. No number or condition below changes after this point (R132).

**Result:** PASSED on 11 Nov 2026 (run by the author, seed 0, 11 018 s): R = 1.016 / 0.991 / 0.882 / 1.008 (± 0.05).
Replication pending.

## Hypothesis

A kink of F3.14a at rest, with no force, makes a random walk in momentum with

  D_p = (1/2ρ) ∫ V′(φ)² φ_x² dx = (32/15) μ^{5/2} / ρ.

This is the shot noise with which the order samples the defect's own binding (docs/derivations/swerves_of_a_defect.md;
derived at the leading order, with no free constant). The test checks the magnitude and the scaling with ρ and μ.

## Model and readout

**Model:**
- the rule of F3.14a, unchanged;
- flat 1+1 trapezoid |x| ≤ T + 10 − t, T = 30;
- kink at rest;
- cells μ ∈ {0.25, 1} × ρ ∈ {400, 800};
- reported: μ = 4 at ρ = 800 (coarsely resolved).

**Readout (frozen):**
1. The kink fit per slab of width 1 gives the centre x_c.
2. The velocity v_k is the slope of x_c over windows of 5 slabs (6 windows).
3. The increments dv_k = v_k − v_0 remove the initial spread of velocities.
4. The slope of Var(dv_k) against dt_k equals 2D_v; a constant readout floor drops out.
5. R = M² D_v ρ / ((32/15) μ^{5/2}), with the SE by bootstrap.

## Instrument check

On every sprinkling an exact kink is moved along a **synthetic random walk in velocity with the predicted D_v** and
sampled at the same events. The estimator applied to its readout must return the D_v of the true synthetic path. The
ratio readout / truth is paired on the same paths. In the smoke run it was 0.99–1.01.

## Criteria (frozen 10 Nov 2026, 22:56)

The sample is set from the development data. At 32 sprinklings the SE of R was 0.17–0.6, so 1024 sprinklings per cell
give SE ≈ 0.03–0.11.

- **Sample:** 1024 new sprinklings per cell.
- **INSTRUMENT_FAILED:**
  - more than 2 % of the fits are lost in any cell; or
  - |readout/truth − 1| + 2 SE ≥ 0.1 in any cell.
- **PASSED:**
  - |R − 1| + 2 SE < 0.3 in each of the 4 cells; and
  - |R(800) − R(400)| − 2 SE < 0.3 for each μ (no density dependence beyond 1/ρ).
- **FAILED:** |R − 1| − 2 SE > 0.3 in any cell.
- **UNDECIDED** otherwise.

**Why 0.3:** the prediction is of leading order (it ignores the feedback on the profile and the radiation), and the
expected 2 SE is up to ~0.2.

**Reported:** μ = 4; R from the cubic growth of Var x (a second estimator); the exponent; the mean drift.

**Cost:** about 30–40 min with 8 processes.

## Reading declared now

**PASSED:** the wandering of a bound defect is the shot noise of the order sampling its own binding, with the predicted
constant and scaling (1/ρ, μ^{5/2}, i.e. M⁵ at M w = 8). The swerves of causal set theory are then derived for a bound
object, set by its internal structure, and zero in the vacuum.

**FAILED:** the leading-order shot-noise picture is wrong or incomplete in magnitude.

**Limits:**
- 1+1 only (in 3+1 the Green function sits on links);
- the potential is chosen;
- no number for real particles.
