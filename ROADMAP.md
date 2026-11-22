# Roadmap

How to work with this laboratory, and what is open. State of 22 Nov 2026.

- Paper: *Physics of Causality*, v2, doi:[10.5281/zenodo.24427611](https://doi.org/10.5281/zenodo.24427611)
  (v1: 10.5281/zenodo.23853689).
- Code and records: this repository (`run.py`, `core/`, `phases/`, `results/`).
- Interactive tour: `Demo/causal_physics.py` (Windows build in the GitHub releases).

## The idea in five lines

The world is a causal order of events; nothing else is primitive. Geometry is order plus number (a Lorentz-invariant
Poisson sprinkling; a substrate with a rest frame is not supported). Gravity is the response of the order to matter:
clocks slow, cones narrow, events crowd, and conserving causal throughput gives γ = 1 in 3+1. Matter is a persistent
pattern of the order: flows that carry their four-velocity, and defects of a field on all events whose inertial mass
equals their rest energy. The expected Benincasa–Dowker action of the order and an order-only matter action together
give the linearised Einstein equation with a non-local spin-0 form factor counted from intervals.

## The rules of the laboratory (keep them)

1. **Pre-registration.** Every test has a card (`docs/*_card_draft.md`) with its question, instruments, criteria and
   readings, written and frozen before the run. Criteria are never changed after a result is known.
2. **Four verdicts:** PASSED, FAILED, UNDECIDED, INSTRUMENT_FAILED. Negative and undecided results are kept and
   reported (`docs/findings.md`).
3. **Exercise before freezing.** Every judged point is run at full resolution in development before the card is
   frozen; numbers seen before freezing are disclosed in the card.
4. **Provenance.** `python run.py <ID>` saves every full run with the hashes of the test file and of the core, the seed
   and the library versions (`results/<phase>/<ID>/`). Nothing is overwritten. Frozen test files are never edited; a
   change is a new revision with a new ID.
5. **Decisions are written down** (`docs/author_decisions.md`): every choice of design, criterion or interpretation,
   with its reason. New decisions continue this log, each with its author's name.
6. **Levels are declared:** mechanism / chain / prediction; ensemble level (expected values) versus single
   sprinklings; development (`diag/`, decides nothing) versus frozen tests (`phases/`).

## Where things are

| What | Where |
|---|---|
| verdicts and their reading | `docs/findings.md` |
| list of tests and their latest state | `docs/test_list.md` (`python run.py report`) |
| planned tests | `docs/planned_tests.csv` |
| the author's decisions 1–48 | `docs/author_decisions.md` |
| derivations | `docs/derivations/` (start with `path_to_predictions.md`, `order_variation_OV1-3.md`) |
| frozen tests | `phases/<phase>/` |
| development scripts | `diag/` (run as `python -m diag.<name>` from the repository root) |
| the paper's LaTeX and figures | `paper/` (`latexmk -pdf main`; figures from the records: `paper/figures/make_figures.py`) |

## What is open, in order of priority (state of 22 Nov 2026)

1. **F6.8i: the static response of old fields** (`docs/derivations/f68i_dev.md`). The optical age term grows as
   (kT)⁵ at second order (F6.8e3) and faster at finite amplitude ε = b(kT)² ≤ 4 (F6.8f, FAILED: no saturation). The
   question is whether ∫ε^(−1/2)Φ(ε)dε converges beyond the first caustic (ε > π², several rays per interval).
   - Done: E1 (the smooth part beyond the caustic falls steeply where there are several rays; indication only), D1
     derivation and check 1 (the two-well needle volume near a cut locus reproduces the exact needle).
   - Next: D1 checks 2–3 (the window integral; the q² and q⁴ window terms must vanish by the BD moments), the q⁶ window
     term, the smooth part integrated piecewise between cut loci, D2 (conjugate points), then the card.
   - Why it matters: it decides whether old fields have a finite stiffness, which the far-field law, β of the Sun and
     the short-range prediction (tens of micrometres) depend on.
2. **F6.12: β of the young field** (`docs/derivations/order_variation_OV3.md` §10, `docs/F6_12_card_draft.md`, draft).
   β − 1 = −κ₃/8, κ₃ the cubic static action along the exponential metric (ω, ν) = (−U, 2U). Stage β1 done (conformal
   cubic → 0.9995, `diag/ov3_beta1.py`); next β2: interval volumes to third order in the cone perturbation (a Fermat
   term T₃ for a stratified medium), then the test.
3. **F6.8h:** the spin-2 form factor φ₂ measured independently (static transverse-traceless shear).
4. **T3 in 3+1:** the Noether energy of a field defect conserved under the order's dynamics, as a frozen test
   (2D development: `docs/derivations/order_variation_OV2.md` §13–14, `diag/ov2_t3_*`).
5. **Matter in 3+1:** a field defect with a dynamics in 3+1 (the 3+1 Green function exists: F3.16b, F5.3) and a fall
   rule that works in 3+1 (longest chains fall in 2+1, F5.1d, but not in 3+1, F5.1e).
6. **The selection principle beyond linear order:** the balance of action and entropy, Γ_λ = −(1/λ) ln Σ e^(−λS)
   (decision 41); the entropic correction to the stiffness at quadratic order is not measured.
7. **The value of G** and the meaning of the age T of a static field (F8.8).

## How to start

    python -m pip install numpy scipy matplotlib     (Python 3.10+)
    python run.py core                  the checks of the core
    python run.py list                  all tests
    python run.py results F6.8g         a record
    python Demo/causal_physics.py       the tour

Read the paper's Part V and `docs/findings.md` (section "Holds" first), then `docs/derivations/f68i_dev.md` for the
work in progress. Keep the rules above; they are what makes the results worth something.
