# Test inventory (25 Nov 2026)

From `python run.py list` on the author's PC (25 Nov, 11:49) and this repository (branch claude/c0-minvar). The
inventory is a reading of the records. It decides nothing.

## 1. The original laboratory

The tests marked "[older version of the test]" that predate this repository (first commit 10 Nov 2026) come from the
original laboratory. They are the tests that rejected the earlier published theory of causality, and they are kept
as its history:
- K.12 (s16, s32, s08_R200, h16), K.13a, K.14, K.15;
- F0.6 (fixed), F1.5, F1.7, F1.9a;
- F3.6, F3.6c, F3.7, F3.8, F3.9, F3.10, F3.10b, F3.11;
- F6.3, F6.3b, F6.4a, F6.8b;
- F8.1.

Their test files are not in the git history, so the records cannot be re-run from the current code.

## 2. Marks that are not history: line endings

On the author's PC, K.18b, F3.16b and F5.2c2 show "older version". Their files are CRLF there, and with the CR removed
they are byte-identical to this repository. The provenance hash (`core/provenance.sha256_file`) hashes raw bytes, so a
change of line endings alone breaks the link. F6.9a shows the same effect the other way round (current on the PC,
"older" in this repository). Normalising line endings in the hash would remove this; it is a change to core, for the
author to decide.

## 3. F5.1e: two configurations

The records on the PC show:

| record | memory / horizon | status |
|---|---|---|
| 30 Oct, seed 0 | 8 / 8 | FAILED |
| 1 Nov, seed 0 | 20 / 20 (the version in git) | UNDECIDED |
| 18 Nov, seed 1 | 8 / 8 | FAILED |

The PC file has 8 / 8, git has 20 / 20. The author (25 Nov): the horizon was enlarged for the measurement. So the
"FAILED (replicated)" in the list comes from the old horizon, and the measurement configuration (20) has one record:
UNDECIDED.

## 4. What stands (PASSED)

- **Instruments:** K.2, K.4–K.9, K.11, K.16b, K.16c, K.18b.
- **Field:** F1.2b (d3), F1.3, F1.4.
- **Motion and inertia:** F2.5, F4.1, F3.14b, F3.14c, F3.15, F3.13a (static, flow d4), F3.13b2.
- **Gravity:** F5.1d, F5.2a3, F5.2b, F5.2c2, F5.2d, F5.3.
- **The BD action:** F6.8a, F6.8b2, F6.8c2, F6.8e3, F6.8g, F6.9a, F6.10a, F6.10b.
- **Tensor coupling:** F7.4a.

## 5. Negative results of the current programme (FAILED)

- **The F3.3 line (b–h):** a localised persistent structure does not arise from load, memory, threshold or clock
  feedback.
- F3.13d (spontaneous defects).
- F5.1, F5.1c, F5.1g (free fall without the apex construction); F5.1e (see §3).
- F7.1a (capacity does not give γ = 1).
- F6.8f (no static saturation of the optical sector up to ε = 4).
- **Instrument failures with a successor:** K.16 → K.16b, K.18 → K.18b, F6.8c → F6.8c2, F6.8e/e2 → F6.8e3,
  F1.2 → F1.2b, F3.12a → F3.12a2, F3.13b → F3.13b2, F3.14a → F3.14b, F3.16 → F3.16b, F5.2c → F5.2c2.

## 6. Open

- **UNDECIDED:**
  - F2.7, the closing chain source → field → motion;
  - F1.2b (d4), F3.13c, F6.4c;
  - F7.1b and F7.1c (replicated);
  - F5.1e in the measurement configuration.
- **INSTRUMENT_FAILED with no successor:** F6.8d (the age term of the mixed sector), F5.1f, K.11c (far_d4_g2),
  F3.13a (flow_d3).
- **Never run:** K.1, K.12 (new variants w, h, hw), F1.4b, F1.4c, F3.3h, F5.1b, F7.4, F0.6 (dynamic),
  F3.12b (confirm), F6.3 (A).

## 7. The cusps of F6.8i and the chain

f68i_dev §6.2 found that each cusp of the first arrival adds a term of order γ^{11/7} to the age term of an old
field (ε > π²); C = −8.83e-3 for the normal form. Consequences for the records:
- **F6.8f and F6.8g stand.** They are young fields (below the first caustic), where there are no cusps.
- **F6.8d** (the mixed age term) meets the same cusps for old fields. It is a candidate for revision once F6.8i is
  settled.
- **F8.1** (far-field law, old record) rests on the age term. It has to be recomputed once it is known how the cusp
  term enters.

## 8. Order of work (approved by the author, 25 Nov)

1. Finish F6.8i:
   - D2a full run;
   - the cusp constants σ′, τ′, N′ along ε and the number of cusps;
   - the direct check of the age term from the exact needle.
2. F6.8d.
3. F2.7.
4. F8.1.
