# Astrophysical layer: points and observations for future development

Opened 15 Nov 2026 at the author's request. **Status: notes only. Nothing here is a result, a test or a part of the
theory yet.**

Astrophysical tests come after the classical base, in the author's order of work of 15 Nov 2026 (roadmap):
- the age terms of the sectors (F6.8e2, F6.8d);
- the asymptotic static propagator;
- the covariance of the sectors;
- dynamic and nonlinear gravity.

Without the propagator the theory has no derived far field of a point source. An astrophysical test before that would
test parameters put in by hand, not the theory.

## 1. The motive (the author, 15 Nov)

Gravity is still a mystery, as fire and electricity once were. Understanding came first and technology followed, by
decades, through exact and testable equations.

The theory may explain dark matter and dark energy, but that is for astrophysical tests to decide. There may also be
a technological potential: if gravity is the response of the causal order, the order might be influenced otherwise
than by matter. This is speculation.

## 2. What the present theory could bring to these scales (open, not results)

- **Non-locality.** On a curved order the order-only field follows a Green function weighted by the curved volume
  along whole near-null pairs, and not a local field (F5.3, C1, replicated). This is a finite-density result on one
  conformally flat well.
- **Memory of age.** A static field responds according to how long it has existed (the age terms of the action, F6.8).
  - Development: in the optical sector the formal second-order coefficient grows as (kT)⁵; in the mixed sector as
    (kT)³ (F6.8e2 and F6.8d decide).
  - The age ties a departure from Newton to a time, in the end to the age of the source or of the universe. This is
    reminiscent of a₀ = cH₀/2π in MOND, but the analogy is not derived.
  - The development numbers point rather to large scales, and the dependence on the field strength (the breakdown near
    b²(kT)⁴ ~ 1) goes the opposite way from MOND, where weak fields depart.
- **What must be derived first:** the static propagator from the action of the order (F8.1), giving the far field of a
  point source and the scale (if any) at which it departs from Newton, with no parameter fitted to galaxies.

## 3. The old GTS validators (archive; ideas, not tests)

Two scripts from the GTS period were shown on 15 Nov.

**"Cross-scale validator".** A tanh map from the baryonic acceleration g_bar, χ = ln(g_bar/a₀), to an effective
dimension D_eff: 2.05 below a₀ ("compression", dark-matter-like), 3 on a plateau, up to 3.6 above a threshold near
white dwarfs ("expansion"). From D_eff it builds a ratio of potentials Φ/Φ_N = γ^p.

**"Redshift WD validation" (test 8).** An observational layer that does not evolve a causal kernel.
- It maps catalogue values (logg, mass) into the phase space χ_obs = ln(g/a₀).
- It computes the GR gravitational redshift and a mapped D_eff, c_eff, Z*.
- Its rule: deviations from GR are allowed only where D_eff > 3 (the "rigidity barrier").
- If a catalogue has an observed redshift column, it computes the residuals obs − GR and obs − GTS.

**Why neither is a test by today's rules:**
- **The law sought is put in.**
  - a₀ = cH₀/2π (the MOND coincidence) is inserted by hand.
  - The steepness 3.0 is chosen to "align with MOND".
  - VACUUM_RIGIDITY = 36.5015 (six significant digits) and CHI_PIVOT = 36.5 are tuned so that the expansion regime
    switches on at white dwarfs.
  - GAIN, POWER, IV_AMP and β_excess are free.
  - The roadmap rule: the laboratory must not put the law it looks for into the result.
- **The ratio is not derived.** γ^p with a chosen p does not follow from "flux ~ 1/r^(D−1)". A change of dimension
  would change the r-dependence, not a factor at one point.
- **No criterion, no comparison.** Rotation curves are not compared with observed velocities; the redshift residuals
  have no criterion fixed before the evaluation.
- **Not covariant.** A local acceleration is not an invariant: a freely falling observer does not feel it. A dimension
  set by g needs a defined frame.
- **Constrained regimes.**
  - The gravitational redshift of white dwarfs is measured (Sirius B about 80 km/s) and agrees with GR within its
    errors.
  - Neutron stars are constrained by binary pulsars and NICER mass–radius measurements.
  - Any expansion regime there has to pass these first.

**What stays useful:**
- **The intuition.** The vacuum may respond differently to weak and strong fields, with a scale at which gravity
  departs from Newton. This is legitimate, but today the scale must come from the theory.
- **The data pipeline of test 8.** It reads WD catalogues: Jiménez-Esteban 2023, MWDD, Gaia samples, logg → g,
  r = √(GM/g), and the robust CSV header handling. It can be reused for a real test once a derived prediction exists.
- **White dwarfs as a test bed.** Their gravitational redshift is an independent observable, with g ~ 10⁶ m/s², where
  any strong-field departure would show.

## 4. How an astrophysical test must be built (when the time comes)

1. **The prediction from the theory:** the static propagator (F8.1), with every parameter fixed in the laboratory and
   nothing fitted to the astrophysical data. Derived, calibrated and predicted are kept apart (gate F8).
2. **A frozen card** with the data set, the selection, the observable, the criterion and the null model (GR, and MOND
   as an alternative), fixed before the data are read.
3. **The data:**
   - galaxy rotation curves (e.g. the SPARC catalogue) for the weak field;
   - white-dwarf gravitational redshifts for the strong field.
4. **Checks where MOND-like ideas break:** weak lensing, the Bullet Cluster (offset of mass and gas), the acoustic peaks
   of the CMB, the solar system and binary pulsars. A theory that replaces dark matter must pass all of them, not only
   the rotation curves.
5. **Replication on an independent catalogue,** and negative results kept.

## 5. Open points to come back to

- At what scale, if any, does the far field of the order's action depart from Newton, and does it depend on the field
  strength, the acceleration or the age?
- Does the memory of age give a dependence on the history of a galaxy (old versus young systems)?
- Dark energy: is there a route from the growth of the order (the passage of time) to an accelerated expansion, and
  what would it predict for w(z)?
- Technology: what would it take to change the order other than by matter, and is any such handle in the theory?
