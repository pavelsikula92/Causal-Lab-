# Findings of the quantum section

The verdicts of the quantum tests (IDs Q.n) and their reading, kept apart from the classical findings
(`docs/findings.md`). The notes and the rules of the section are in `docs/quantum/README.md`, the planned tests in
`docs/quantum/planned_tests.csv`.

State of 28 Nov 2026: the section is open (decision 49); one verdict (Q.1).

## Holds

- **Q.1: one sprinkled 1+1 order carries Johnston's sum-over-chains propagator; its error falls as ρ^{−1/2} at fixed
  mτ, both at a fixed pair and over the order's own elements (PASSED, calibration)**
  - **Setup:** card docs/quantum/Q1_card_draft.md draft 1, frozen (decision 50). Light-cone square, τ_xy = 1, m = 2, 5,
    10, N_exp = 1000, 4000, 16000, 2000 realisations per point. Run by the author 28 Nov on the frozen file (sha256
    d58fd228, core 6eac15c4), 870 s. Record results/Q_quantum/Q.1/20261128T102351Z_main_seed0_d58fd228.json.
  - **Instruments (all OK):**
    - I1: the recursion against the dense solve, 1.0e-15;
    - I2: the exact mean ½J₀(m) at all nine points, max |z| 2.40 (tol 3.5);
    - I3: the non-negative null (b = +m²/ρ) has mean 13.6208 against ½I₀(5) = 13.6199 (+0.07 SE). It is positive,
      without the zeros of J₀.
  - **Criterion:** all six slopes in band [−0.65, −0.40] (bootstrap 95 %).

    | m | M1 (spread at the pair) | M2 (error of one order over its elements) |
    |---|---|---|
    | 2 | −0.501 [−0.517, −0.484] | −0.490 [−0.506, −0.478] |
    | 5 | −0.502 [−0.518, −0.487] | −0.494 [−0.503, −0.485] |
    | 10 | −0.510 [−0.526, −0.493] | −0.504 [−0.510, −0.498] |

  - **Reported:**
    - The slopes of median |e| and q90 |e| are −0.494 … −0.512. The tails narrow at the same rate as the bulk; there
      are no rare large errors that stay.
    - The relative error of one order at N_exp = 16000 is 2.1 % / 6.1 % / 7.3 % of |½J₀| at m = 2 / 5 / 10, and 8.3 %
      / 24 % / 30 % at N_exp = 1000. The precision is set by m²/ρ: at a given density, heavier fields are carried less
      precisely by one order.
  - **Reading (declared on the card):**
    - The first instrument of the quantum section is calibrated. The propagator does not appear only in the average
      over orders: it is carried by the single realised order, with Poisson self-averaging.
    - This is a calibration; it has no physical conclusion of its own.
    - It says nothing about 3+1 (Q.1b), interference against the non-negative flow (Q.2), the Born rule (Q.3), or
      m²/ρ → 1.
  - **Replication:** not yet (a run with another seed, --seed K).

## Does not hold

(none yet)

## Open

(none yet)
