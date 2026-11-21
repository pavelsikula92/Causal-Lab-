# Results of the Causality Lab

Every full run of a test is saved here with what is needed to reproduce it: the parameters, the seed, the hash of the
test file and of the core, the versions of the libraries, the printed lines and the data. Nothing is overwritten.

## Layout

    results/<phase folder>/<ID>/<UTC time>_<part>_seed<seed>_<test hash 8>.json   (and the same lines as .txt)
    results/<phase folder>/<ID>/quick/...                                          reduced runs; they never count

The phase folders are those of `phases/`. Each has a README.md with the objective of the phase and the question of
every test in it.

| Folder | Phase | Objective (short) |
|---|---|---|
| `K_instruments/` | K | instruments and their calibration; no physical conclusion of their own |
| `F0_substrate/` | F0 | what the substrate is: a Lorentz-invariant causal set or a finite-degree graph |
| `F1_field/` | F1 | transport of influence and the field of a source |
| `F2_motion/` | F2 | motion along maximal chains in a given order |
| `F3_structure/` | F3 | what matter is: flows, bound bundles, locked defects, defects of a field of the order |
| `F4_inertia/` | F4 | inertia and contact (merging with four-momentum conserved) |
| `F5_gravity/` | F5 | whether matter falls in an order shaped by matter; equivalence |
| `F6_order_dynamics/` | F6 | the dynamics and the action of the order; the field equation of the order |
| `F7_light/` | F7 | how the response of the order fixes gamma; the tensor coupling |
| `F8_data/` | F8 | confrontation with data |

Other entries:
- `history_3x/`: results of the predecessor GTS Lab 2.x and 3.x (`history_3x/overview.csv`).
- `sha_aliases.json`: maps hashes of test files whose code is unchanged and whose text was edited only (the English
  translation of 30 Oct 2026), so that their records stay current.

## Reading

- `python run.py results <ID>`: all records of a test.
- `python run.py report`: regenerates `docs/test_list.md` with the latest state of every test.
- `docs/findings.md`: the verdicts and their reading.

## The earlier flat layout

Until 21 Nov 2026 the records were in `results/<ID>/`. They were moved unchanged (file names and contents) into the
phase folders. `core.provenance` writes to the new layout and still reads the old one, so nothing is lost if an old
folder reappears. Documents written before that date (the decisions, cards and derivations) cite paths as
`results/<ID>/...`; the same file is now at `results/<phase folder>/<ID>/...`.

One frozen test keeps its own cache path: F1.9a writes `results/F1.9a/cache/` (a cache of intermediate values, not a
record). A run of F1.9a will recreate that folder.
