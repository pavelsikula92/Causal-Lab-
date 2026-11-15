# Astrophysical layer: points and observations for future development

Opened 15 Nov 2026 at the author's request. **Status: notes only. Nothing here is a result, a test or a part of the
theory yet.**

Astrophysical tests come after the classical base. In particular they wait for the asymptotic static propagator of
the action of the order, which fixes the far field of a point source (F8.1, P1), and for the age terms of the sectors
(F6.8d, F6.8e2). An astrophysical test before that would have no derived prediction to test.

## 1. Motivation (the author, 15 Nov)

- The author expects the theory to bear on dark matter and dark energy. Astrophysical tests are to decide it.
- Gravity is still a mystery, as fire and electricity once were. Understanding it may carry technological potential.
  If it exists, it lies where the theory differs from general relativity: gravity as the response of the causal
  order, which might be changed by other means than mass. This is a question, not a claim.

## 2. Where the present theory could depart from Newton and Einstein

**Features of the present work that might lead to departures from GR.** None is a prediction yet.

- **Non-locality of the field response.** On a curved 3+1 order the order-only field follows its exact mean, a Green
  function weighted by the curved volume along each near-null pair, and not the local conformally coupled field (F5.3,
  C1, replicated). This is a statement about the response of a field, not about the fall of a body: that a localised
  structure accelerates towards a source has not been shown in 3+1 (F5.1e; open question).
- **Age dependence of the static response of the action.** The conformal age term is linear in T (F6.8a). In
  development, the mixed and optical sectors grow faster, ∝ (kT)³ and ∝ (kT)⁵. Three steps stand between this and a
  prediction:
  1. F6.8d and F6.8e2 must decide the exponents.
  2. An independent test must decide whether the age dependence is physical (the history of a static field in its rest
     frame) or an artefact of the age cut.
  3. The asymptotic static propagator must be derived from the action.

  Only then can a departure from GR, its scale and its dependence on the field strength be stated. Development numbers
  (the sign change of the optical stiffness near kT ≈ 43; the breakdown of the expansion in the field strength b near
  b²(kT)⁴ ~ 1) are not to be read as scales of a prediction.

## 3. Old GTS scripts (archive; reviewed 15 Nov, not to be revived as they are)

### 3.1 CROSS-SCALE VALIDATOR (GTS, Phase 4 kernel)

**What it does.** It maps the baryonic acceleration g_bar to χ = ln(g_bar / a₀), with a₀ = cH₀/2π. From χ, a tanh
gives an effective dimension D_eff:
- compression towards 2.05 below a₀ (meant as dark matter);
- a plateau D = 3 in between;
- expansion towards 3.6 above a "vacuum rigidity" χ = 36.5015 (white dwarfs, neutron stars).

From D_eff it builds a ratio Φ_GTS/Φ_Newton = (3/D_eff)^p. It prints the ratio for point objects and two galaxies.

**Why it would not pass today:**
- **It inserts the law it looks for.** a₀ = cH₀/2π is the MOND coincidence put in by hand. The steepness 3.0 is chosen
  "to align with the MOND interpolation". VACUUM_RIGIDITY = 36.5015 has six significant digits: it is placed so that
  the expansion starts at white dwarfs. The laboratory must not put the law it looks for into the result (roadmap,
  Aim).
- **The ratio is not derived.** (3/D_eff)^p with a hand-chosen p does not follow from "flux ∝ 1/r^(D−1)", which would
  change the r-dependence, not a factor at one radius.
- **It is not a test.** It has no criterion and no comparison with data: the rotation curves are not compared with
  observed velocities.
- **Physical obstacles:**
  - the local acceleration is not an invariant (a freely falling observer feels none), so "dimension by g" needs a
    defined frame;
  - the expansion regime would change the gravitational redshift of white dwarfs, which is measured (Sirius B,
    about 80 km/s) and agrees with GR;
  - neutron stars are constrained by binary pulsars and NICER mass–radius measurements.

**What remains:** the intuition that the vacuum responds differently to weak and strong fields, with a scale where
gravity departs from Newton. Today such a scale must come out of the theory (the propagator), not out of a parameter.
Whether the age dependence of the action provides one is open (section 2).

### 3.2 REDSHIFT WD VALIDATION, test 8 (catalog mapping)

**What it does.** It reads white-dwarf catalogues (Jiménez-Esteban 2023, MWDD, Gaia samples) and converts log g
(cgs) to g. It estimates the radius from r = √(GM/g) and computes the GR gravitational redshift
z_GR = (1 − 2GM/rc²)^(−1/2) − 1. It then maps χ_obs = ln(g/a₀) through a tanh into "I_v, Π, Z*, D_eff, c_eff". It
allows a departure z_GTS = z_GR · Z*^β only where D_eff > 3 (the "rigidity barrier"). If a catalogue has an observed
redshift or radial-velocity column, it reports the residuals obs − GR and obs − GTS. The catalogues used had none.

**Assessment:**
- **The mapping is inserted, not derived:**
  - the pivot χ = 36.5 is the same tuned number;
  - the gain is 30;
  - α_g = 0.007297 is the fine-structure constant without a reason;
  - the amplitude 0.15 and β are free.
- **The barrier rule cannot fail.** With a tuned pivot, "deviations only where D_eff > 3" can be placed wherever the
  data allow.
- **Without observed redshifts it tests nothing.** The headers it was run on contain none.
- **Worth keeping: the data layer.** Robust catalogue ingestion (headers starting with '#', column normalisation),
  log g → g, r from M and g, and z_GR are the right infrastructure for a later test.

**Why white dwarfs are a good test bed.** Gravitational redshifts of white dwarfs are measured:
- Sirius B and 40 Eri B;
- wide binaries and common-proper-motion pairs, where the space velocity is known from the companion;
- large spectroscopic samples of DA white dwarfs.

They test the mass–radius relation together with z_GR. A theory that departs from GR in compact objects must pass
them.

## 4. What an astrophysical test must be (methodology of the laboratory)

1. **A derived prediction.** The far-field law and any scale come from the propagator of the action of the order
   (F8.1). The parameters are fixed by the laboratory, never fitted on the astrophysical data used for the test (F8:
   derived, calibrated and predicted kept apart).
2. **A frozen card.** The data set, the selection, the criterion and the null model (GR / Newton) are fixed before
   the evaluation.
3. **Several independent witnesses together:**
   - rotation curves (for example SPARC) and the radial-acceleration relation;
   - gravitational lensing (galaxies, clusters, the Bullet cluster);
   - cluster dynamics;
   - the cosmic microwave background and baryon acoustic oscillations;
   - the solar system (perihelia, Cassini γ);
   - binary pulsars;
   - white-dwarf redshifts.

   MOND-like alternatives pass the first and fail others. A replacement of dark matter has to pass all at once.
4. **Dark energy.** Whether the history of a field matters on cosmological scales needs, after the steps of section 2,
   the dynamic (non-static) response of the order: the fourth item of the author's order of work (dynamic and nonlinear
   gravity).

## 5. Order (when the classical base stands)

1. Asymptotic static propagator → far field of a point source → P1 (deviation from 1/r and from γ = 1, with the scale).
2. A first confrontation where the theory is already constrained: the solar system and white-dwarf redshifts (the
   theory must not fail where GR is tested precisely).
3. Galaxies: rotation curves and lensing with the same prediction, no parameter refit.
4. Cosmology: after the dynamic response of the order.
