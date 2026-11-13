# Card B3: active mass from the causal growth of a 2+1 order (B3a, the active side)

**Revision 2, frozen 13 Nov 2026 (decision 23), implemented as F5.2d**
(phases/F5_gravity/F5_02d_growth_active_2p1.py; seed 0 decides). Revision 1 (the same morning) was the design; revision 2
fixes the tolerances from the development rows by the rule of 10 Nov. Decisions 15, 16, 17, 19, 22, 23; the author's
choices of 13 Nov (B3a: the cones cut in the rest frame of the static source; B3 = the active side, the fall of a test
body moved to B4). Helpers: diag/b3_growth.py (growth), diag/b3_continuum.py (continuum).

## 1. Hypothesis

A 2+1 order that grows event by event, with its coordinate density ρ Π⁻¹ and its cones Π² set by a potential L = ln Π
that matter produces from its realised causal past through the Green weight of the order (the chain weight of B0, the
load = trace of T), becomes around a static source the order of the theory, −Π² dt² + Π⁻² (dx² + dy²), as much as the
source's load says (active mass): the potential, the slowed clocks and the narrowed cones of the continuum response
box_g L = κ q.

## 2. Model

As the docstring of F5.2d. In short: candidates in time order; parents by the cone rule |Δr| < Π_m² Δt with Π_m from
the binned potential at the midpoint of the pair (B3a, a declared coordinate proxy: the cone at the parent alone made
the cones of long near-null links too wide where the potential deepens in time, and the potential ran away, r_L 1.6);
field and potential by the chain weight; the load by the identity X = (π − φ) box φ − box f, f = φ(2π − φ)/2 (zero on
both vacua; f = φ² of decision 19 is 4π² on a plateau and its retarded cancellation amplified the weight's calibration
error); gravity switched on from t 3 to 4 (then the flat free wave φ₀ is exact where the order is curved; with gravity
from t 0 the rule differed from the theory by 1.6–2.5× in the continuum); realised iff mark < exp(Lmin − L).
Source: one kink wall at rest, μ 1 (load 16) and μ 0.25 (load 8); κ 0.01; ρ 30; T 9; each sprinkling also at κ = 0.

## 3. Readouts

r_L (the potential of the order projected on the continuum potential, band |y| < 2, |x| < 6, t ∈ [5, 8)); u = r_L(μ 1) −
r_L(μ 0.25); clocks = the near/far ratio of longest-chain lengths over Δt ≈ 2.75 (κ against κ = 0) against the continuum
∫Π dt; cones = the same ratio of the reach in y of the causal future against Π² at the midpoint; density = realised /
expected.

## 4. Criteria (frozen, revision 2)

| level | condition | tolerance | derivation (2 SE₃₂ + development offset) |
|---|---|---|---|
| instrument | A0: new marks for the latest 10 % of candidates leave every earlier event unchanged (2 sprinklings) | exact | exact |
| instrument | no event with L < Lmin (−0.7) | 0 | lowest −0.36 |
| instrument | \|density − 1\| + 2 SE, each μ | < 0.01 | 0.001 + 0.0007 |
| instrument | \|cones/cont − 1\| + 2 SE, each μ (the rule realises the cones it was given) | < 0.03 | 0.005 + 0.003 |
| potential | \|r_L − 1\|, μ 1 | 0.06 | 0.012 + 0.006 (per slab offsets to 0.08) |
| potential | \|r_L − 1\|, μ 0.25 | 0.12 | 0.036 + 0.040 |
| active mass | \|u\| | 0.12 | 0.037 + 0.046 |
| clocks | \|clocks/cont − 1\|, each μ | 0.06 | 0.015 + 0.024 |
| clocks (PASSED only) | (1 − clocks) − 2 SE > 0, each μ | — | continuum 0.122 and 0.046 |

Every condition on its own; |·| + 2 SE on the PASSED side, |·| − 2 SE on the FAILED side. Verdict order:
INSTRUMENT_FAILED → FAILED → PASSED → UNDECIDED. 32 sprinklings.

## 5. Reading declared now

**PASSED:** a causal growth rule that reads only the realised past (counts, chains, independent marks) and lets the
density and the cones follow the potential that matter produces through the order's own Green weight makes, around a
static wall, the order of the theory in 2+1, in proportion to the wall's load: active mass from the order, with the
relation changing (not only the density, as in 1+1).
**FAILED:** the grown order does not carry the potential, the clocks or the active mass of the continuum response.
**Limits:** the cones are cut in coordinates in the rest frame of the static source (B3a; the intrinsic reading of the
frame is B3b); one source geometry (a planar wall), one κ, one density; the fall of a test body is not part (B4); the
potential is chosen; κ is chosen.

## 6. Development rows (13 Nov; decide nothing)

ρ 30, 4 sprinklings, diag/b3_growth.py, diag/out/b3_dev_20261113T071448.json; mean ± SE (sd):

| quantity | μ 1 (load 16) | μ 0.25 (load 8) |
|---|---|---|
| r_L (t 5–8) | 0.994 ± 0.017 (sd 0.035) | 1.040 ± 0.050 (sd 0.100) |
| r_L per slab 5.5 / 6.5 / 7.5 | 1.076 / 1.049 / 0.961 | 1.157 / 1.121 / 0.999 |
| L at x = 0 (last slab) | −0.216 (continuum −0.221) | −0.094 (−0.096) |
| clocks near/far | 0.899 ± 0.022 (cont. 0.878) | 0.956 ± 0.012 (0.954) |
| cones near/far | 0.759 ± 0.008 (0.762) | 0.910 ± 0.006 (0.908) |
| density − 1 | −0.0006 ± 0.0016 | −0.0007 ± 0.0017 |
| kink drift | +0.04 ± 0.11 (κ 0: +0.04 ± 0.14) | −0.01 ± 0.04 |
| width w/w_ex − 1, last slab | −0.26 ± 0.01 (the discrete kink breathes after the switch) | 0.00 ± 0.01 |

Cost ~200 s per growth (one core of the container). Before these rows: the cone at the parent alone (r_L 1.6 at ρ 10)
and the load identity with φ² (asymmetric potential, the kink drifted) were replaced (section 2).
