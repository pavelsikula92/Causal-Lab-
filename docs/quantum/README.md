# Quantum layer — points and observations for future development

Opened 9 Nov 2026 at the author's request. **Status since 28 Nov 2026: the quantum layer is a section of the
laboratory of its own (decision 49): its tests (IDs Q.n), cards, planned tests, findings and results are kept apart
from the classical ones (section 7). Nothing here is a result yet.** The notes below (sections 1-6) are the notes of
9 Nov, kept as they were. The quantum layer is deliberately postponed until the classical base stands (the author's decision,
9 Nov): motion of matter from conservation (F5.1g, F5.1h), the field equations (the Jacobson route), waves and frame
dragging from the dynamics of the order, and the strong field. A quantum test before that would not tell whether a
failure comes from the phase or from the classical rule.

## 1. The stance of the theory (what is already written)

- **Ontology (paper, part II):** the passage of time is the growth of the order; the future is open because it has not
  been realised. Superposition is the coexistence of admissible causal continuations; measurement is a causal
  interaction that biases the routing towards one of them.
- **First version (February 2026):** "the substrate does not select a future trajectory; histories with higher
  accumulated resistance carry lower transition weight" — a statistical, non-teleological selection.
- **Methodology:** M2 (readouts never act on the dynamics), M4 (the state of an event depends only on its causal past;
  the same result in every frame and every linear extension), M5 (finite without a wall).

## 2. Observations from the discussion of 9 Nov

- **Realism.** The world beyond a horizon exists; a horizon only says that it never sends the observer a signal. The
  order is absolute (M4); a horizon is a relation between the world line of an observer and the order (which events fall
  into its past cone). "What I do not perceive does not exist" is one interpretation of quantum mechanics, not quantum
  mechanics.
- **Observation is a causal event.** Not perceiving has no effect. Measuring is an interaction: it acts on the future
  of the observed (the routing of what is not yet realised), never on its realised past or internal order.
- **Indeterminacy as the unwritten part of the world** (the author): quantum behaviour is a corner/edge case of the
  world where geometry, and with it time, is not yet written. Observation, a causal event, determines a position at
  that level of the world — a place in the causal order — not a point in the SI / Euclidean picture, which is only a
  coarse image of the order and has no meaning below the discreteness scale.
- Related views in the literature: the growing (evolving) block universe (G. Ellis); "becoming" in causal set theory
  (R. Sorkin); relational quantum mechanics (C. Rovelli: facts are relations between systems) is probably closest to the
  causal intuition. Realist alternatives to keep in mind: Bohmian mechanics (non-local, positions always defined),
  objective collapse (GRW; Penrose–Diósi gravitational collapse), many worlds.

## 3. What a quantum layer must deliver (the hard tests)

1. **Interference: amplitudes, not only weights.** Sharing a flow among continuations with non-negative weights (F5.1f,
   F5.1g) is probability, not quantum mechanics. Quantum mechanics needs contributions that can cancel (a phase).
2. **The Born rule:** why probabilities are |amplitude|².
3. **Bell correlations without signalling.** Correlations between distant measurements stronger than any local theory
   with predetermined values allows, while no signal travels faster than light. A realist account must give up
   something (locality of correlations, a single outcome, or measurement independence). Rideout and Sorkin formulated
   "Bell causality" for classical sequential growth; classical growth obeys it and cannot violate Bell inequalities.
   A quantum growth must. This is where the realist reading of the theory will be decided.
4. **The classical limit / decoherence:** why the realised order looks classical at large scales.
5. **ħ.** In the Jacobson route the temperature of a local horizon must be proportional to the acceleration (T = α κ)
   for the field equations to have their form; the value α = ħ/2π (Unruh) is quantum and fixes G = 1/(4ħη). This is
   the single point where the classical base touches the quantum one. Open idea (speculation): if the substrate has
   its own quantum of action (energy × time per event), ħ would be a property of the substrate, not an external input.

## 4. Known formal tools on causal sets (to read before designing anything)

- **Quantum measure / decoherence functional** on the space of histories (Sorkin): the sum-over-histories form of
  quantum mechanics on a causal set; quantum sequential growth is open.
- **Propagators as sums over chains** (S. Johnston, 2008–2010): a particle propagating on a causal set as a sum over
  chains with amplitudes per link and per element reproduces the retarded and Feynman propagators of a free scalar
  field in the continuum limit. This is the closest known step to the flows of the laboratory: replace the
  non-negative shares by amplitudes.
- **Sorkin–Johnston vacuum state** of a free field on a causal set.
- **Entanglement entropy** on causal sets (Sorkin, Yazdi): needs a truncation of the spectrum; an alternative to the
  number of horizon molecules as the entropy of the Jacobson route.
- **Unruh effect** on a causal set (a detector coupled to the SJ state) — the quantum input of the Jacobson route.

## 5. Ideas for future laboratory tests (not designed, not scheduled)

- **Double slit on a causal set:** a flow with amplitudes through two openings of a region; interference fringes as a
  function of the phase rule; the classical (non-negative) flow as the null model.
- **Bell-type correlations in growth:** two separated regions grown from a common past; can a growth rule produce
  correlations above the classical bound without signalling?
- **Unruh detector:** response of an accelerated detector on a sprinkled order with the SJ state; T ∝ κ?
- **Gravitationally induced collapse / entanglement (Penrose–Diósi, BMV):** external experiments; the theory should
  say what it predicts for them once the quantum layer exists.

## 6. Rules for this layer (9 Nov; the first rule was replaced by decision 49, section 7)

- Nothing here starts before the classical pillars have their verdicts.
- Every test gets a card, frozen criteria and a null model (the classical, non-negative flow), as in the rest of the
  laboratory; R132 holds.
- The interpretation (section 2) is a guide for building rules, not a result; it is confirmed only through tests of
  section 3.

## 7. The quantum section of the laboratory (28 Nov 2026, decision 49)

**Why now.** The author opened the layer on 28 Nov, before all classical pillars have verdicts. The reason given in
the discussion: the classical age term does not close by itself (hypothesis H3: no continuum limit for old fields),
and records and measurement need a quantum account (H4). The danger named on 9 Nov stays: a quantum test built on an
open classical rule could not tell which of the two failed. So the first quantum tests are **calibrations of quantum
instruments on the bare substrate** (a Poisson sprinkling of flat spacetime), which depend on no open classical
pillar. Tests that use the classical rules (fields, defects, the age term) come only after these instruments pass.

**Layout (kept apart from the classical section).**

| What | Where |
|---|---|
| notes and the stance (this file) | `docs/quantum/README.md` |
| planned tests, ID;Question;Level;Depends on;Note | `docs/quantum/planned_tests.csv` |
| cards | `docs/quantum/Q<n>_card_draft.md` |
| findings (verdicts and their reading) | `docs/quantum/findings.md` |
| test files | `phases/Q_quantum/Q_<nn>_<name>.py` (registered when the first card is frozen) |
| records of full runs | `results/Q_quantum/<ID>/` |
| development scripts and records | `diag/q<n>_*.py`, `diag/out/q<n>_*` |

**IDs.** Q.1, Q.2, ... with revision letters as elsewhere (Q.1b). The phase Q is added to the registry of the
laboratory (`core/registry.py`: the phase list and the file-name pattern) when the first quantum test file is
written. That is a change of the core, so records from then on carry a new core hash, as with decision 48.

**Rules (replacing the first rule of section 6).**
- Every quantum test has a card with a frozen question, instruments, criteria and a **null model**: the classical,
  non-negative flow, or the same construction with non-negative weights. A quantum effect counts only where it
  differs from the null.
- A quantum test may use only instruments that passed their own calibration (Q.1 first).
- Interpretation (section 2) guides the construction of rules. It is not a result.
- The standing rules of the laboratory hold unchanged (M1-M7, R132, smoke test before every run).

**The first line of tests (docs/quantum/planned_tests.csv).**
1. Q.1: the sum-over-chains propagator (Johnston) on a 1+1 sprinkling: exact mean, self-averaging of one order.
2. Q.1c: the vacuum of one order: the imaginary part of Johnston's Feynman propagator (the SJ state) against the
   continuum SJ state of the same diamond (the author chose it before interference, 28 Nov).
2b. Q.1b: the same in 3+1 with the link matrix; the mean is not exact there.
3. Q.2: two-route interference on a sprinkled order, against the non-negative null.
4. Q.3: the quantum measure of a two-route experiment (Sorkin's decoherence functional): the second-order interference
   term I₂ and the vanishing of the third-order term I₃ (the Born-rule test).
5. Q.4: the first arrival as an amplitude sum. Does the cusp term γ^{11/7} of the classical age term survive when the
   minimum over rays is replaced by a sum with phases (H3, test 3)?
6. Q.5: Bell correlations in a growth of two separated regions from a common past.
7. Q.6: records: redundancy of a fact in sub-orders (H4, test 3).
