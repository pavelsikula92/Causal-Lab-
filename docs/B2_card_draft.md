# Card B2: walls of a field in a prescribed curved 2+1 order of the theory (intrinsic instrument)

**Revision 2, frozen 12 Nov 2026 (decision 22), implemented as F5.2c**
(phases/F5_gravity/F5_02c_walls_curved_2p1.py; seed 0 decides, seed 1 replicates). Revision 1 (development, the same
day) set the design and the readouts; revision 2 fixes the tolerances from the development rows by the rule of 10 Nov
(section 7). Decisions 12, 16, 17, 20, 22. Derivation: docs/derivations/walls_in_curved_order.md. Order and continuum:
diag/b2_curved.py. Development: diag/b2_dev.py.

## 1. Hypothesis

The field rule of B1 (F3.15, PASSED: the Green function read from the order alone, the longest-chain weight), applied
unchanged on an exact sprinkling of the 2+1 metric of the theory (narrower cones and more events), moves walls of a
field as the continuum field theory in that metric does. The continuum makes three sharp predictions:

- **the kink does not fall.** A wall uniform along y has the passive weight E − σ (energy minus tension); the kink has
  σ = E and is an exact static solution at every position. A test particle would fall 0.44 in the same time.
- **the breather wall falls.** ω = 0.7, σ/E = 0.22: it falls as the continuum (about 0.78 of a test particle).
- **a launched kink follows the cones of the order.** Its proper speed stays constant, so its coordinate speed is
  Π² v (0.313 against 0.300 on a flat order, launched at 0.3 from Π = 1 towards larger Π).

B2 is the prescribed-order stage (decision 16): the order is given, not grown; the active mass is B3. The light bending
γ is built into an exact sprinkling (the relation is that of the metric) and is not tested here; it belongs to the
stage where a rule produces the relation (F7.1b-type, B3).

## 2. Model

- **Order:** Poisson sprinkling of density ρ per proper volume of ds² = −Π² dt² + Π⁻² (dx² + dy²) with
  Π² = (x + ℓ)/ℓ, ℓ = 24, in the domain of dependence of the rectangle −13 ≤ x ≤ 22, |y| ≤ 12.5, T 8 (it contains the
  band |y| < 2 and the readout windows of every arm up to t = 7). The relation in closed form (the optical metric is the
  hyperbolic plane of radius ℓ); links by the shell method with the exact relation, checked against the brute-force
  transitive reduction. Π = 1 at x = 0; across the domain the coordinate density of events varies 1.7×, the
  coordinate speed of light from 0.46 to 1.9.
- **Field rule:** exactly as F3.15 (variant B, the frozen constants a, c, β; μ = 1).
- **φ₀:** the free wave (μ = 0) of the same initial data in the same metric, from the y-independent continuum solver
  (leapfrog, dx 0.005; convergence 1e-5). As in B1, the initial data are the only place where coordinates enter.
- **Arms** (initial data at x₀ = 0, s = proper distance along x):
  kink at rest 4 atan e^s; kink launched at proper velocity 0.3 towards +x; breather wall ω = 0.7 at maximal amplitude
  (φ_t = 0). The three arms share each curved sprinkling. **Flat control:** the breather on a flat order (ℓ = ∞, the
  same rectangle), on its own sprinklings.

## 3. Readouts (coordinates declared, as in B1)

Band |y| < 2, slabs of width 1, centres 0.5 … 6.5 (the slab rule of B1: the last slab is not read). Every readout is
repeated on the exact continuum field at the same events (suffix _ex; paired).
- kink: the fit of field1d.fit_kink (centre x_c, width w, Q);
- breather: the centre x_b of a fit of the profile 4 atan(A(t)/cosh(k (x − x_b))), A quadratic in t over the slab,
  over |x| < 5 (reported also: the centroid of the noise-debiased energy density over |x| < 7, the readout of B1);
- D_K = [x_c(last) − x_c(first)] − [x_c,ex(last) − x_c,ex(first)] (kink at rest);
- dv_L = slope(x_c) − slope(x_c,ex) over the read slabs (launched kink);
- r_B = Σ_t d(t) d_ex(t) / Σ_t d_ex(t)², d = x_b − x_b(first slab) (the breather fall projected on the continuum fall);
- D₀ = [x_b(last) − x_b(first)] on the flat order, divided by the mean curved continuum fall d_ex(last);
- cones = slope(x_c) − 0.3 for the launched kink (flat cones: 0; the continuum: +0.013).

## 4. Criteria (frozen, revision 2)

Tolerance rule (10 Nov): 2 SE at the card's sample size from the development sd per sprinkling + the largest
development offset + a margin. Sample size **32 curved sprinklings** (three arms each), the flat control on 8.

| level | condition | tolerance | derivation (2 SE at n + offset) |
|---|---|---|---|
| instrument I1 | the links equal the brute-force transitive reduction (curved and flat, ~1100 events) | identical | identical in development |
| instrument I4 | linear extension t + 0.2x; the smaller rectangle (xr 6 against 8) | 1e-12 | as B1 |
| instrument I0 | \|N/(ρ V_proper) − 1\| + 2 SE | < 0.01 | −0.002 ± 0.0015 (n 4) |
| instrument | kink at rest \|Q − Q_ex\| + 2 SE | < 0.01 | 0.0002 + 0.0004 |
| instrument | kink at rest \|w/w_ex − 1\| + 2 SE | < 0.05 (as B1) | 0.0023 + 0.016 |
| instrument | flat breather \|D₀\| + 2 SE (in units of the curved continuum fall) | < 0.25 | n 8: 0.11 + 0.03 |
| the kink stays | \|D_K\| | 0.10 (23 % of the fall of a test particle, 0.44) | sd 0.057 (B1 rest drift: 0.12): 0.02–0.04 + 0.02 |
| the launched kink moves as the continuum | \|dv_L\| | 0.02 | 0.0047 + 0.0086 |
| the launched kink follows the curved cones (PASSED only) | (slope − 0.3) − 2 SE > 0 | continuum +0.013 | dev +0.022 ± 0.007 |
| the breather falls as the continuum | \|r_B − 1\| | 0.15 | 0.047 + 0.039 |

Each condition on its own, |·| + 2 SE on the PASSED side and |·| − 2 SE on the FAILED side. Verdict order as B1:
INSTRUMENT_FAILED → FAILED (D_K, dv_L or r_B with |·| − 2 SE above the tolerance) → PASSED (all, including the cones) →
UNDECIDED. A kink that falls like a test particle gives D_K ≈ −0.44 (FAILED); a breather that does not fall gives
r_B ≈ 0 (FAILED); a kink that moves along flat cones gives dv_L ≈ −0.013 (the cones condition not met: UNDECIDED).
The binding conditions are r_B (the breather oscillates about once in the run) and the cones (the effect is 0.013).

## 5. Decisions (made: decision 22, 12 Nov)

1. The order: the exact sprinkling of the metric of the theory, Π² = (x + ℓ)/ℓ, ℓ 24, ρ 30, T 8, the rectangle
   −13 ≤ x ≤ 22, |y| ≤ 12.5 (flat control −16 ≤ x ≤ 22). Not the F7.1b rule (coordinate proxy, undecided γ); that is B3.
2. Arms: kink at rest, kink launched at 0.3, breather ω 0.7, on the same 32 curved sprinklings; the breather on 8 flat
   sprinklings.
3. The breather readout: the profile fit (section 3); the energy centroid is reported only.
4. γ is not part of B2.
5. The slab rule of B1.

## 6. Reading declared now

**PASSED:** the order of the theory, prescribed and read only through its relation, moves walls of a field as the
continuum metric does: no spurious force on a wall whose tension equals its energy, the right fall of a wall with little
tension, and motion along the narrower and wider cones. The passive response of extended bodies in 2+1 comes from the
order.
**FAILED:** at the tested scale the intrinsic rule on a curved order gives a force the continuum does not have (the
kink moves), or the wrong fall of the breather, or the wrong motion along the cones.

**Limits:** one profile (Π² linear in x, ℓ 24) and one density; the order is an exact sprinkling of a chosen metric,
not produced by matter (B3); walls only (no 1/r, B4); the breather readout is coarse (the period, 9 time units, is
longer than the run).

## 7. Development rows (12 Nov; decide nothing)

ρ 30, ℓ 24, T 8; diag/b2_dev.py, diag/b2_breather_readout.py; diag/out/b2_dev_*.json, b2_dev_cache/, b2_dev_fields/.
Mean ± SE (sd per sprinkling, n):

| quantity | development | continuum |
|---|---|---|
| kink at rest D_K | +0.020 ± 0.029 (sd 0.057, n 4) | 0 (a test particle −0.44) |
| kink at rest w/w_ex − 1; Q − 1 | +0.016 ± 0.003; −0.004 | |
| launched kink dv_L | +0.0086 ± 0.0067 (sd 0.013, n 4) | slope 0.3132 |
| launched kink slope − 0.3 | +0.022 ± 0.007 | +0.013 (flat cones: 0) |
| breather r_B (profile fit) | 0.961 ± 0.054 (sd 0.132, n 6) | 1; the fitted centre moves −0.78 by t 6.5 |
| breather r_B (energy centroid, \|x\| < 7) | 0.82 ± 0.16 (sd 0.40, n 6) | the centroid moves −0.30 |
| flat breather D₀ (profile fit) | −0.14, +0.08 (n 2) | 0 |
| density N/(ρV) − 1 | −0.002 ± 0.0015 (n 4) | |

Cost: ~900 s per curved arm and ~1150 s per flat arm on one core of the container (N 109 000 and 127 000).
The first development round read the breather by the energy centroid of B1 (|x| < 7); its per-sprinkling sd 0.40 would
need ~200 sprinklings for a tolerance of 0.15, so the readout was replaced by the profile fit before the card was
frozen (development decides nothing; the energy centroid stays as a reported line).
