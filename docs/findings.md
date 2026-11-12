# Findings and discoveries

State of 21 Oct 2026: results of GTS Lab 2.x and 3.x and the first runs of Causality Lab (calibrations K.2–K.9, K.11, K.11b, F1.2, F1.2b, F1.3, F1.4, F1.5 preliminarily, F2.5 d4 and d3, F2.7, K.12 s16, s08_R200 and s32, F6.3 V1 and V2, F6.3b V1, V2 and A, K.14 all 5 parts, F3.6 rest and move, F3.10 unit, attraction, momentum, universality and law, F3.10b law, K.13a all three parts, F1.7 V, W and N, F0.6 fixed, F3.7 search, F3.8 search, F3.9 search, K.11c, F3.3 all three channels, F3.3b q90, F3.3c main_q99, F3.3d the whole scan, F3.3e search and confirm, F3.3f capacity and source, F3.3g A); the development diagnosis of K.11b and the derived conflict of (a) with H1 in 3+1. Entries added later carry their own dates.
The full list of tests and their latest states is in `docs/test_list.md`.

The strongest result: paths in the field Π follow the geodesic, the flux of a source satisfies Gauss's law and the
field of a source has a Newtonian shape far out with a strength according to Gauss's law (F1.3); one order reproduces
the field in 2+1 (F1.2b).
The largest gap: no matter-like structure has yet arisen from the dynamics by itself.

## Holds

- **F5.2b: under a causal growth rule of the order (candidates realised with probability Pi^2 computed from the realised
  causal past, ln Pi = -(kappa/2 rho) sum of the load 2V), matter thins its future order, and a test defect falls towards
  it as the continuum response predicts, in proportion to the source's energy (r_B = 1.029 / 1.017 / 1.009 for M_A 4 /
  8 / 16) and independently of its own mass (u = -0.0007 +- 0.0081), also for a moving system (r_sep 1.014 +- 0.008)
  (seed 0 PASSED; replication pending)** (mechanism, 1+1, rho 800, kappa 0.003, two fields interacting only through the
  order, kinks at -3 and +3, T 20, 256 sprinklings each paired with kappa = 0; card docs/F5_2b_card_draft.md rev. 1,
  run by the author 11 Nov without changes, 5065 s; results/F5.2b/).
  - **Instrument:**
    - A0 and A1 hold;
    - density: realised / Pi^2 = 0.9998 +- 0.0001;
    - intervals near A: count / proper area 0.9972 +- 0.0023;
    - null: the B displacement at kappa = 0 is zero within 2 SE for every combination.
  - **Level 1 (growth):** count / coordinate area 0.6211 +- 0.0027. The counted proper times near matter are 38 % shorter.
  - **Level 2 (active response), r_B = mean d_B / d_cont:** 1.0288 +- 0.0170 (M_A 4), 1.0165 +- 0.0094 (8), 1.0085 +-
    0.0059 (16). B moved -0.559 / -1.113 / -2.160 against the continuum -0.543 / -1.095 / -2.142, i.e. 1 : 1.99 : 3.87.
  - **Level 3 (universality):** u = r_B(M_B 4) - r_B(M_B 8) = -0.0007 +- 0.0081.
  - **A2 (boosted launch, v 0.3):** separation change -1.772 against -1.748, r_sep 1.0139 +- 0.0079.
  - **Reported:**
    - r_A (A towards B) 1.010 / 0.992 / 0.965 / 0.977;
    - third law M_A d_A / (M_B |d_B|) 0.965-0.978 for the static combinations (the 1.835 of the launched arm is not a
      meaningful ratio there);
    - L at the end 2-10 % below the continuum in magnitude (partly the slab average of a growing field);
    - the excess of the falls over the continuum is now 1-3 %, smaller for the heavier source.

  Reading (declared before the run). A causal growth rule of the order, which uses only the realised past, the counts and
  independent marks, makes matter thin its future order (slower counted clocks). A test defect then falls towards the
  source in proportion to the source's energy and independently of its own mass, as box_g L = -kappa q predicts, and
  this holds also for a moving system. With F3.14b and F5.2a3, the inertial, passive and active mass come from one field
  and one rule (1+1). All three levels of the author's list are met, including level 3, which closes the test of the
  active mass. The rule is a candidate for the Axiom of the dynamics of the order.

  Limits:
  - 1+1: a scalar response (R2), no 1/r, and the potential of a static source grows;
  - the causal relation does not change (conformal flatness), as it will have to in 2+1;
  - kappa is chosen (G is not derived);
  - the potential is chosen;
  - the card was run without a separate freezing note; its text was not changed between the smoke run and the run.

- **F3.14c: a field defect at rest on the 1+1 order makes a random walk in momentum with the diffusion constant
  predicted by the shot noise of the order, D_p = (32/15) mu^(5/2) / rho, with no free constant: R = 1.016 / 0.991 /
  0.882 / 1.008 (+- 0.05) in the four cells (seed 0 PASSED; replication pending)** (mechanism, 1+1, flat, kink at rest,
  T 30, 1024 sprinklings per cell, mu 0.25 / 1 x rho 400 / 800; card docs/F3_14c_card_draft.md frozen 10 Nov 22:56;
  run by the author 11 Nov, 11 018 s; results/F3.14c/). Derivation: docs/derivations/swerves_of_a_defect.md.
  - **Instrument:** no fit lost in the decided cells. The readout check (an exact kink moved along a synthetic random
    walk of the predicted D_v, read like the data) returns its truth: 1.000 +- 0.000 (paired).
  - **R = D_p,meas rho / ((32/15) mu^(5/2))**, from the variance of the windowed velocity increments:
    - M 4: 1.016 +- 0.053 (rho 400), 0.991 +- 0.052 (800);
    - M 8: 0.882 +- 0.048 (400), 1.008 +- 0.055 (800);
    - density scaling: R(800) - R(400) = -0.025 +- 0.074 (M 4), +0.126 +- 0.072 (M 8).
  - **Reported:**
    - M 16 at rho 800: R = 0.898 +- 0.045 (1.7 % of fits lost).
    - The second estimator (cubic growth of Var x): 1.077 / 1.007 / 0.880 / 0.990 / 0.883.
    - Exponent of the spread 1.44-1.48 (random walk in momentum: 1.5).
    - Mean drift zero within 1 SE.

  Reading (declared before the run). The wandering of a bound defect on the order is the shot noise with which the
  order samples its own binding, with the predicted constant and the predicted scaling: 1/rho, and mu^(5/2), i.e. M^5
  at M w = 8, over a factor 32 in D_p between M 4 and M 8, without any fitted parameter. The swerves of causal set
  theory are thereby derived for a bound object, set by its internal structure (int V'^2 (d phi)^2), and zero in the
  vacuum.

  Not resolved at this precision: the corrections beyond the leading order. M 8 at rho 400 lies 2.5 SE below 1, and M 16
  lies 10 % below.

  Limits:
  - 1+1 only (in 3+1 the Green function sits on links);
  - the potential is chosen;
  - no number for real particles follows.

- **F3.14b (revision of F3.14a for resolving power): a defect of a field of the 1+1 order keeps its identity, its
  velocity and its Lorentz contraction, and its inertial mass equals its rest energy to 0.8 %: r_M = 0.9988 +- 0.0077
  (PASSED, REPLICATED: seed 0 and seed 1)** (mechanism, 1+1, deciding density 800 with 128 sprinklings, controls on 32,
  rho 100 / 200 / 400 with 32 each; criteria of card docs/H5_card_draft.md rev. 4 unchanged; run by the author 10 Nov,
  2241 s; results/F3.14b/).
  - **Replication (seed 1, run by the author 10 Nov, 1079 s): PASSED.**
    - Instrument OK; I3 eps 0.0255 / 0.0197 / 0.0120 / 0.0086; drift at 800 +0.004 +- 0.013.
    - H5-B at 800: v 0.3 +0.0013 +- 0.0014, v 0.6 -0.0007 +- 0.0008; C - 1 within 0.0007.
    - **H5-C: r_M = 0.9895 +- 0.0079** (seed 0: 0.9988 +- 0.0077).
    - Telegraph stops; blind L 0.375.
    - **The reported deceleration after the kick did not recur:** gamma v [4, 7) - [7, 10) = -0.0003 +- 0.0011 (seed 0:
      +0.0046 +- 0.0011). The open point below is closed as a fluctuation of seed 0 plus the readout effect of +0.0011
      found in the continuum.
  - **Instrument (all OK):**
    - I3 eps 0.0249 / 0.0192 / 0.0110 / 0.0088 (monotonic);
    - H5-A at 800: drift +0.015 +- 0.015 (2 SE 0.03 against the tolerance 0.1), (E - E_ex)/8 +0.0013, M_rest 8.021;
    - H5-E: the telegraph kink stops (p 0.002-0.005 against 0.319 and 0.081);
    - V = 0: L 0.0005.
  - **H5-B at 800:** v 0.3: v_meas - v +0.0015 +- 0.0014, C - 1 -0.0008 +- 0.0005; v 0.6: +0.0006 +- 0.0009, -0.0007 +-
    0.0009; e_B +0.0016.
  - **H5-C at 800:** Delta p 0.0781 +- 0.0006 against 0.0785; M_in 8.011 against M_rest 8.021; r_M 0.9988 +- 0.0077.
    At rho 100 / 200 / 400: 0.994 / 1.016 / 1.014.
  - **Blind control:** L 0.373 against 1.003.
  - **Reported:**
    - The field momentum after the kick is 0.6304 against 2 pi f tau = 0.6283.
    - gamma v over [4, 7) exceeds gamma v over [7, 10) by +0.0046 +- 0.0011, about 6 % of the kick. The momentum of the
      field is conserved, so this is probably the readout of the centre being shifted by the uniform vacuum oscillation
      left after the switch-off (period 2 pi).
      - **Checked the same evening (decides nothing).** The continuum kicked kink, read with the same fit on 600 000
        random points, gives +0.0011 (P conserved at 0.6283; the vacuum oscillates with an amplitude ~0.08). The
        readout explains a quarter of the effect.
      - The discrete excess +0.0035 +- 0.0011 (3 SE) is open. Earlier values: +0.003 to +0.007 at rho 100-400 in F3.14a
        and F3.14b, each within 1-2 SE.
      - Candidate: after the switch-off the vacuum oscillates, so V' != 0 everywhere and the shot noise of the order
        (docs/derivations/swerves_of_a_defect.md) acts on the whole field. The kink may then pass momentum to that
        radiation while the total momentum stays conserved. Needs a dedicated diagnostic (the kink's local momentum
        against time, several densities).
    - At rho 100 and 200, v 0.3 runs +0.015 / +0.013 (2 SE) fast: a discreteness effect at low density.

- **F5.2a3 (revision of F5.2a2 for resolving power): in a 1+1 order whose only change is the number of events, kinks of
  mass 4, 8 and 16 fall as their continuum bodies to 0.1-0.8 %, and alike within ~1.5 % (PASSED, REPLICATED: seed 0
  and seed 1)** (mechanism, 1+1,
  512 sprinklings at density 800, 32 at 400 reported; criteria of card docs/F5_2a_card_draft.md rev. 2 unchanged; run
  by the author 10 Nov, 2316 s; results/F5.2a3/).
  - **Instrument OK:** I4 <= 1.5e-15; flat |D| <= 0.048; counting null D +0.0003 +- 0.0067 against -2.145 with counting.
  - **Falls:** a = 0.9991 +- 0.0010 / 1.0001 +- 0.0028 / 0.9976 +- 0.0079.
  - **Equivalence (frozen statistic, mean of the paired ratios a_i / a_j - 1):** (4, 8) +0.0018 +- 0.0021, (8, 16)
    +0.0291 +- 0.0072, (4, 16) +0.0359 +- 0.0090. All within the tolerance.
  - **The launched arm:** 1.003 +- 0.015 of its continuum. The implementation error of F5.2a is gone.
  - **The non-zero eta of the pairs with M 16 is a bias of the statistic, not a violation (derived after the run;
    decides nothing).**
    - The mean of a ratio is biased by the variance of its denominator: E[1/a] ~ (1/a)(1 + CV^2).
    - The fall of M 16 has CV = 0.387 / 2.15 = 0.18, so the bias is CV^2 = 0.032, which matches the measured +0.029 and
      +0.036.
    - The ratio of the means (SE by bootstrap) gives (4, 8) -0.0010 +- 0.0021, (8, 16) +0.0025 +- 0.0062 and (4, 16)
      +0.0015 +- 0.0075. The paired differences a_i - a_j give the same. No violation is resolved, at the level of
      ~1.5 % (2 SE).
    - The same bias explains the replication of F5.2a2 (+0.040, +0.046).
  - **Replication (seed 1, run by the author 10 Nov, 3535 s; results/F5.2a3/20261110T202629Z_all_seed1_81a79075.json):
    PASSED.**
    - Instrument OK; counting null +0.012 +- 0.007.
    - a = 0.9965 +- 0.0010 / 0.9949 +- 0.0028 / 0.9928 +- 0.0078.
    - Frozen eta: +0.0046 / +0.0277 / +0.0373, with the same bias of the mean of ratios.
    - Ratio of means, or paired differences: (4, 8) +0.0015 +- 0.0022, (8, 16) +0.0021 +- 0.0059, (4, 16) +0.0037 +-
      0.0075.
    - Launched arm 1.028.
    - **Observed in both runs (decides nothing):** at rho 800 all masses fall 0.1-0.7 % less than their continuum
      bodies (M 4: 0.9991 and 0.9965). A small discreteness correction, common to all masses, so it does not touch the
      equivalence.
  - **Lesson for future cards:** compare means (a ratio of means, or paired differences of D), not means of ratios,
    when the denominator is noisy.

- **F5.2a2: in a 1+1 order whose only change is the number of events per coordinate area, the kink of F3.14a under
  the unchanged rule falls as its continuum body, and kinks of mass 4, 8 and 16 fall alike: the weak equivalence
  principle from counting, within 2 % (seed 0 PASSED; replication seed 1 UNDECIDED on two WEP pairs with M 16: not
  replicated)** (mechanism, 1+1, Phi = g L tanh(x/L), g 0.05,
  L 10, trapezoid |x| <= 26 - t, T 10, 128 sprinklings at the deciding density 800, 32 at 400 reported; PASSED 10 Nov by
  the frozen criteria of card docs/F5_2a_card_draft.md rev. 2, unchanged from F5.2a; corrected sprinkling; run by the
  author, 463 s; results/F5.2a2/).
  - **Replication (seed 1, run by the author 10 Nov, 466 s): UNDECIDED.**
    - Instrument OK.
    - Falls: a = 0.9995 +- 0.0021 / 1.0026 +- 0.0060 / 0.9931 +- 0.0165.
    - eta (4, 8) +0.0004 +- 0.0048 (passes); eta (8, 16) +0.0397 +- 0.0156 and (4, 16) +0.0461 +- 0.0197. Each is above
      0.06 with 2 SE and not FAILED.
    - The launched arm: ratio 0.989 +- 0.035.
    - Not replicated. The runs are not pooled (R132).
  - **Design error, the same as in F3.14a:** the resolving power was set from 16 development sprinklings, which
    underestimated the wandering of the narrow kink M 16.
    - Expected SE of eta (4, 16): 0.010. Actual: 0.015 (seed 0) and 0.020 (seed 1).
    - sd of D for flat M 16 at rho 800: 0.28 in development, 0.41-0.44 here.
    - The tolerance 0.06 then sat within 2 SE of an offset of 0.02-0.04.
  - **Rule:** phi(x) = phi_0(x) - (1/2 rho) sum over the causal past of mu sin phi. It is the rule of F3.14a,
    unchanged; it knows nothing about the metric. The order is that of flat spacetime with event density rho Omega^2.
  - **Instrument:** I4 <= 1.8e-15; flat control |D| <= 0.023 for every mass; the fall is resolved; counting null
    (every term divided by Omega^2) D -0.022 +- 0.013 against -2.154 with counting.
  - **Fall** (a = D / D_cont, the continuum body with the same initial data and readout):
    - M 4: 0.9999 +- 0.0020;
    - M 8: 1.0042 +- 0.0054;
    - M 16: 1.0189 +- 0.0149.
  - **Equivalence:** eta (4, 8) -0.0017 +- 0.0041, (8, 16) +0.0055 +- 0.0115, (4, 16) +0.0077 +- 0.0146 (tolerance
    0.06 with 2 SE).
  - **Reported:**
    - The raw D_4 / D_8 - 1 = -0.024 +- 0.004 is the continuum finite-size effect (the wide kink falls 2.2 % less than a
      point in the continuum too); the order reproduces it.
    - Launched against the field: ratio to the continuum 0.93 +- 0.03.
    - At rho 400: M 16 a = 0.961 +- 0.042, eta (8, 16) +0.11 +- 0.05. The narrow kink (w = 0.5) is coarsely resolved at
      the lower density.

  Reading (declared before the run). In an order whose only change is the number of events, a bound defect falls as
  its continuum body in the geometry the number defines, and defects of mass 4, 8 and 16 fall alike within the
  resolving power of 6 %. The measured violation is |eta| <~ 1-3 % (2 SE); this is not a claim at the level of
  MICROSCOPE. With the counting removed it does not fall: the fall comes from the number of events. With F3.14a,
  inertia and passive gravitational mass come from the same rule.

  Limits:
  - 1+1, where density and clock rate go together (the kink falls towards the sparser order), unlike 3+1;
  - the geometry is prescribed (no back-reaction: F5.2b);
  - the potential is chosen (H5-D);
  - the launch, the readout and the counting null use coordinates.

- **F5.2a2: a defect of a field of the order falls in a 1+1 order whose only change is the number of events, as its
  continuum body in the geometry that the number defines, and kinks of mass 4, 8 and 16 fall alike within 1.5 %; with the
  counting removed it does not fall** (mechanism, 1+1, Phi = g L tanh(x/L), g 0.05, L 10, sprinkling of density
  rho exp(2 Phi) in |x| <= 26 - t, T 10, 128 new sprinklings at the deciding density 800, 32 at 400 reported; PASSED
  10 Nov by the criteria of card docs/F5_2a_card_draft.md rev. 2, frozen 13:53, unchanged from F5.2a; run by the author,
  463 s; results/F5.2a2/). Revision of F5.2a (implementation error in the sprinkling, no conclusion). Not yet
  replicated with another seed.
  - **Rule:** that of F3.14a, unchanged; it knows nothing about the metric. In 1+1 every metric is conformally flat,
    so the order is that of flat spacetime and only the number of events per coordinate area changes.
  - **Instrument (all OK):**
    - I4 <= 1.8e-15;
    - flat control |D| <= 0.023 for every mass;
    - fall resolved;
    - counting null D -0.022 +- 0.013 against -2.154 curved (the same order with every term divided by Omega^2).
  - **Fall (a = D / D_cont):** 0.9999 +- 0.0020 (M 4), 1.0042 +- 0.0054 (M 8), 1.0189 +- 0.0149 (M 16).
  - **WEP (eta of a, paired):** (4, 8) -0.0017 +- 0.0041, (8, 16) +0.0055 +- 0.0115, (4, 16) +0.0077 +- 0.0146.
  - **The finite size of the body is reproduced.** The raw ratio of falls (4, 8) is -0.0237 +- 0.0040: the wide kink
    (w = 2) falls 2.4 % less than the narrow one, as the continuum extended body does (-0.022). The order reproduces
    the tidal and non-equilibrium effect of an extended body, not only the geodesic.
  - **Reported, rho 400 (32 sprinklings):** M 4 and M 8 as at 800. The narrowest kink (M 16, w = 0.5) falls 0.961 +-
    0.042 and eta (8, 16) +0.11 +- 0.05: a resolution effect of the narrow kink at the lower density (about 400 events
    per unit width squared), gone at 800. The launched arm against the field: 0.93 +- 0.03 (rho 800), 1.06 +- 0.09
    (rho 400) of its continuum.

  Reading (declared before the run). In an order whose only change is the number of events, a bound defect falls as
  its continuum body in the geometry the number defines, and defects of different mass fall alike. This is the weak
  equivalence principle from counting, at the resolving power of this experiment (6 %; achieved about 1.5 %; not a claim
  at the level of MICROSCOPE). The fall comes from the counting: with the counting removed the kink does not fall. With
  F3.14a, inertia and passive gravitational mass come from the same rule.

  Limits:
  - 1+1, where the density and the clock rate are the same function, so the kink falls towards the sparser order,
    opposite to 3+1;
  - the geometry is prescribed (no back-reaction: F5.2b);
  - the potential is chosen;
  - the launch, the readout and the counting null use coordinates.

- **F3.14a (seed 0 PASSED, replication seed 1 INSTRUMENT_FAILED on the drift alone; not replicated): a topological
  defect of a field on all events of a 1+1 causal order, evolved by the retarded Green function
  of the order (1/2 x the causal relation), has inertia: it keeps its identity at rest, keeps its velocity and Lorentz
  contraction when launched, and answers a force with an inertial mass equal to its rest energy (r_M = 0.994 +- 0.019),
  while an overdamped dynamics stops it and a causality-blind sum dissolves it** (mechanism, 1+1, sine-Gordon kink,
  trapezoid |x| <= 26 - t, T 10, ~170 000 events per arm at the deciding density 400, 32 new sprinklings per density
  100 / 200 / 400; PASSED 10 Nov by the criteria of card docs/H5_card_draft.md rev. 4, frozen by the author at 12:05;
  run by the author, 176 s; results/F3.14a/20261110T111802Z_all_seed0_e37cf249.json). Hypothesis H5 of
  docs/cosmogenesis.md; decisions 14 (D1-D4, D2').
  - **Replication (seed 1, run by the author 10 Nov, 193 s; results/F3.14a/20261110T112302Z_all_seed1_e37cf249.json):
    INSTRUMENT_FAILED by one condition, the drift of the rest kink.**
    - The drift was +0.0575 +- 0.0452, so |mean| + 2 SE = 0.148 > 0.1.
    - Every other condition was met as in seed 0: I1-I4, Q, width, energy, H5-E, V = 0; H5-B at both velocities; H5-C
      r_M = 0.9953 +- 0.0217; blind control D_L 0.618 +- 0.004.
    - The PASSED of seed 0 stands under its own criterion but is **not replicated**. The two runs are not pooled after
      the fact (R132).
  - **Design error (found after the replication):** the drift tolerance 0.1 lay at the resolving power of the test. At
    32 sprinklings the spread of the end position of the rest kink is sd 0.26-0.27 (both seeds), so 2 SE = 0.09 alone
    nearly fills the tolerance. Seed 0 passed with 0.098. The card estimated SE ~ 0.025 from the 16 development
    sprinklings, where the spread was sd 0.144; why the development sample was narrower is open (an F-test puts it at
    p ~ 0.01). Pooled descriptively, the two seeds give a mean drift of +0.027 +- 0.032: no systematic drift is seen.
  - **Observed (descriptive, decides nothing): the kink wanders.** The spread of its displacement grows with time faster
    than linearly.
    - Rest kink at rho 400: sd 0.03 / 0.08 / 0.14 / 0.26 at dt 2 / 4 / 6 / 9; the local exponent is ~1.3-1.5, between
      ballistic (1) and a random walk in momentum (1.5).
    - It falls slowly with the density: 0.43 / 0.31 / 0.26 at dt 9 for rho 100 / 200 / 400.
    - It is smaller for a moving kink: 0.11 at v = 0.6, roughly as time dilation and the relativistic mass 1/gamma^3
      would make it.
    - A candidate for a momentum diffusion of a bound object from discreteness, the "swerves" of Dowker, Henson and
      Sorkin (2004), which have not been computed for a bound object. To be separated from readout noise and from the
      noisy radiation background of the field before anything is claimed.
  - **Rule:** phi(x) = phi_0(x) - (1/2 rho) sum over y < x (the order), t_y >= 0, of [sin phi(y) + f chi(t_y)]; phi_0 the
    free wave of the initial data (the launch, the only place where coordinates enter the dynamics).
  - **Instrument (all OK at rho 400):**
    - I4: order, engine and region agree to <= 1.05e-15;
    - I1: delta1 +0.0022 +- 0.0019, error last slab 0.075;
    - I2: |omega / sqrt(k^2 + 1) - 1| <= 0.00013 at every k;
    - I3: eps 0.0254 / 0.0178 / 0.0135 (monotonic);
    - H5-A: Q - Q_ex 0.0000, drift -0.004 +- 0.047, w / w_ex - 1 +0.0005, (E - E_ex)/8 +0.0027;
    - H5-E: the telegraph kink stops (p_Gamma 0.002-0.007 against 0.312 launched and 0.077 kicked; separation 10-57 SE);
    - V = 0: L 0.0005 against 1.001.
  - **H5-B (rho 400):**
    - v 0.3: v_meas - v -0.0000 +- 0.0035, C - 1 -0.0004 +- 0.0013, e_B +0.0021 +- 0.0010;
    - v 0.6: v_meas - v +0.0006 +- 0.0021, C - 1 -0.0019 +- 0.0022, e_B +0.0024 +- 0.0018;
    - Q - Q_ex 0.0000.
  - **H5-C:** Delta p 0.0784 +- 0.0015 against 0.0785 for M = 8; M_in 7.985 against M_rest = E 8.032; r_M = 0.9941 +-
    0.0186. At rho 100 / 200: r_M 1.015 +- 0.049 / 0.988 +- 0.034. The apparent rho^(-1) trend of the development rows
    (1.065 / 1.035 / 1.013, 16 sprinklings) did not recur: no discreteness correction of the mass ratio is resolved at
    the 2 % level.
  - **Blind control:** L 0.376 against 1.001 (D_L 0.625 +- 0.003); width 54 at the end.
  - **Reported:** the free motion after the kick is kept (gamma v [4, 7) - [7, 10) +0.0006 +- 0.0038); the field
    momentum after the kick 0.619 against 2 pi f tau = 0.628.

  Reading (declared before the run). The order has a conservative dynamics of its own: a sum over the causal past,
  stable where the BD operator in retarded form is not. In it, a defect of a field on all events has inertia with
  M_in = M_rest within 2 %, without a carrier flow and without any source after its launch, while a relaxational
  dynamics with the same potential does not give it. The causal structure carries it: a sum over a random set of the
  same size dissolves the kink.

  Limits, stated before the run:
  - in a Lorentz-invariant sprinkling with a covariant rule, a moving kink is in distribution a boosted resting kink,
    so H5-B is largely fixed by symmetry once H5-A holds; the content is the stability of the discrete sum, the absence
    of resolved discreteness corrections, and the separation from relaxation and from a causality-blind sum;
  - 1+1 (no transverse directions);
  - the potential is chosen (H5-D derives it from the finite capacity of the flow);
  - the internal target space is an input;
  - the launch uses coordinates;
  - nothing about gravity (F5.1i).

- **F3.13b2: a hedgehog of the carried orientation is locked: a dynamics that relaxes every unprotected distortion
  does not remove it, randomising its core does not erase its charge, and its core keeps a cost that a twist without
  charge loses** (mechanism, 3+1, revision of F3.13b, box t [0, 20], half 7, ~54 900 events, ~1 045 sources,
  16 new sprinklings, PASSED 31 Oct). Rules approved by the author: the internal orientation n (F3.13a) and its
  alignment n ← N[(1 − κ) n + κ m̂] towards the flow-weighted mean of the causal past with rest-frame lag ≤ 2 (the
  gradient flow of the configuration cost E = Σ(1 − n·n′)).
  - **Instrument OK:** vacuum 0 on 336 of 336 in every κ; ≥ 50 readable shells in every deciding arm; the alignment
    works: the twist at κ = 0.3 relaxes, mean cost 0.0270 → 0.0094.
  - **Charge:** the hedgehog reads +1 on 139 of 139 (κ = 0.1) and 205 of 205 (κ = 0.3) readable shells; after the
    core (r < 2) is randomised at 6 ≤ t < 7, the shells R ≥ 3 read +1 on 56 of 56 and 97 of 97. No core event in any
    arm. κ = 0 (reported): 90 of 90, after the erase 31 of 31.
  - **Cost:** in the last slab at κ = 0.3 the core cost of the hedgehog is 0.275 ± 0.075 against 0.021 ± 0.005 of the
    twist. The hedgehog's core cost rises from 0.22 to 0.27–0.38 while the twist's mean cost falls threefold.
  - **Reported:** after the erase the core cost jumps (0.46–0.59) and relaxes back; at κ = 0.3 it ends at 0.25 against
    0.27 of the undisturbed hedgehog (F3.13b: 0.31 against 0.34). Readability of the hedgehog 0.27 / 0.41 / 0.61 for
    κ = 0 / 0.1 / 0.3.

  Reading (declared before the run): the charge is locked. The alignment removes what nothing protects (the twist
  relaxes towards the vacuum), but it cannot remove the hedgehog, and a complete randomisation of its core is followed by
  a core with the same charge and about the same cost. The core cost is the candidate rest mass in the author's sense:
  the cost of a configuration, not the number of events, and it seems to belong to the charge rather than to the history
  of the core. The first locked structure of the laboratory. Limits: the orientation and the alignment are new
  assumptions; the defect is put in by hand; cold flow at rest, no field; readout to t = 17; the erase is a
  randomisation, not an optimised attempt to unwind; an S² hedgehog's gradient cost grows with its size in 3D (Derrick:
  its cost concentrates in a shrinking core, as seen), so a finite mass independent of the box would need a full
  orientation (a Skyrmion) or a core held at the discreteness scale. Next: the cost of a moving hedgehog (does it grow
  as γ, like a rest energy?), a hedgehog and an anti-hedgehog (annihilation, release of the cost), size and κ.

- **F3.13a (static_d3, static_d4, flow_d4): a defect of the carried orientation has an integer charge that the
  instrument reads exactly, and in 3+1 the cold flow of the laboratory carries a hedgehog with its degree +1 intact**
  (mechanism, 31 Oct; the author's principle of a particle as a topologically locked defect; new carried state
  approved by the author: an internal orientation n ⟂ u, carried unchanged, merged as the flow-weighted sum projected
  on the rest space of the merged u).
  - **Instrument (static, PASSED in 2+1 and 3+1):** the degree of the map shell → sphere of orientations reads vacuum
    0, vortex/hedgehog +1, anti −1, double +2 and the twisted control 0 on all 96 shells of 32 sprinklings. Under noise
    (random rotation of each n up to θ) it never reads a wrong integer: a shell whose field is not resolved
    (neighbours differ by ≥ π/2) is unreadable instead. Readability of the defect: 3+1 1.00 / 1.00 / 0.08 at
    θ = 0 / π/8 / π/4; 2+1 1.00 / 0.98 / 0.69.
  - **flow_d4 (3+1, box t [0, 12], half 7, ~33 000 events, ~1 050 sources, 16 sprinklings, PASSED):** the hedgehog
    reads +1 on 78 of 78 readable shells (readable 0.54), anti −1 on 78 of 78, vacuum 0 on 144 of 144, twist 0 on 54
    of 54, at t = 5, 7, 9; no core event (the orientations never cancelled at a merge).
  - **Reading (declared before the run):** the transport of the laboratory, merging included, does not destroy the
    charge: the defect is a configuration of the flow that the flow carries. This is the necessary condition only: the
    cold flow moves almost straight up and merges only neighbouring streams, so the charge was carried, not challenged.
    Whether it is locked against a perturbation that tries to unwind it, and what erasing it costs, is F3.13b. Limits:
    the orientation is a new assumption; the defect is put in by hand; 3+1 readout over t = 5–9 only; readability
    0.54, just above the limit.

- **F5.1d: a history along longest chains that starts at the apex of a thrown ball falls in an order curved by matter
  with the Newtonian acceleration; no bending rule is put in** (mechanism, 2+1, order of F7.1b σ0 = 0.45, box T = 46,
  32 + 32 sprinklings of ~68 600 events, 40 chains, rule of F5.1c m = H = 8, PASSED 30 Oct; replication with seed 1
  PASSED 31 Oct).
  - **Instruments OK:** flat order R = +0.146 ± 0.098 (|R| + 2 SE = 0.34 ≤ 1/2), longest chains without sag
    R_max = +0.104 ± 0.120; loaded longest chains read the geometry, R_max = 0.742 ± 0.075.
  - **Loaded order:** R = a / a_pred = **1.069 ± 0.051** (criterion 1/2 < R ∓ 2 SE < 2: 0.97–1.17), a_pred = −0.01109
    (φ = −σ/3, fall from rest at r_apex = 8.31). 39.1 of 40 chains reach the end. Displacement from the apex,
    t = 0, 2, …, 22: 0.00 −0.10 −0.33 −0.62 −0.92 −1.26 −1.68 −2.15 −2.66 −3.20 −3.75 −4.30.
  - **Reported, not evaluated:**
    - The fall is stronger than what the longest chains read of the geometry: R / R_max = 1.44.
    - Velocity at the start: the rising half of the history moves outwards at +0.026, the fit from the apex gives
      v0 = −0.069 ± 0.011. The history thus gets an inward velocity jump of ~0.09 at the apex; in the flat order the jump
      is +0.02 (rising half +0.016, v0 +0.037). The criterion reads the acceleration, not the jump.
    - The flat order drifts slightly outwards (+0.43 over 22), within the tolerance.
    - ⟨v²⟩ of the links of the history: 0.369 loaded, 0.569 flat.
    - Development rows (8–10 chains, decided nothing) gave R 1.67 / 2.23; the upper bound R − 2 SE ≥ 2 was added
      before the full run.
  - **Replication with seed 1 (31 Oct, new sprinklings, the same criteria; PASSED, instruments OK):**
    R = **1.101 ± 0.064** (0.97–1.23); flat R = +0.113 ± 0.105; R_max = 0.879 ± 0.073 loaded, +0.007 ± 0.109 flat;
    39.0 of 40 chains to the end. Displacement from the apex (t = 0, 2, …, 22): 0.00 −0.12 −0.35 −0.66 −1.03 −1.42
    −1.86 −2.36 −2.89 −3.45 −4.05 −4.64. The reported quantities repeat too: the jump at the apex (rising half
    +0.036, v0 = −0.080 ± 0.012; flat +0.020 → +0.026), R / R_max = 1.25, ⟨v²⟩ 0.368 / 0.568; the flat order drifts
    outwards and back (+0.34 at t = 16, +0.18 at 22). Both runs together: R = 1.07 ± 0.05 and 1.10 ± 0.06.
  - **Reading:** for the first time in the laboratory a flow falls in an order shaped by matter only through inertia
    (it continues its history) and the longest chains of the order, without a bending rule and without a field put
    into the dynamics. The start condition matters: a history begun at rest (F5.1c) does not fall cleanly, a history
    read from a free state of the order (the apex of a thrown ball) does. The rule "start from the apex" arose after
    F5.1c FAILED; its criteria were fixed before the run (R132). Replicated on new sprinklings. Open: the velocity jump
    at the apex and the excess over R_max (1.44 and 1.25). The same rule in 3+1 does not fall (F5.1e, FAILED 30 Oct,
    below): the result holds in 2+1 only.
  - **Note 8 Nov (F5.1-diag2, decides nothing):** in 3+1 the rule takes the velocity of the region in which it reads
    the order (here a cylinder at rest in the frame of the source); covariance of the rule in 2+1 is not tested.

- **F4.1 (d3 and d4): flows that meet stick together: contact transfers momentum, the merged flow stops in the
  centre-of-momentum frame and the kinetic energy turns into rest mass; the rule has no bounce** (mechanism, 2+1,
  32 sprinklings, ~19 200 events, PASSED 29 Oct).
  - **Head-on collision ±0.6:** mixed fraction of the flow f_mix = 0.821 ± 0.013 (range 0.62–0.94), velocity of the
    merged flow −0.008 ± 0.013 (sticking: 0), the pure flows keep ±0.600, bounce 0. Rest mass of the merged events /
    sum of the parts 1.216 (γ = 1.25 for equal parts; unequal parts give less).
  - **Instruments OK:** without contact (shifted by 2 r + 2) f_mix = 0.003 ± 0.002 (in 2 of 32 sprinklings 0.07 and
    0.03, chains stray across the gap); the recomputation in a frame β' = 0.6 chooses the same; four-momentum on the
    mixed events to 2e−16.

  Reading: a collision of flows is perfectly inelastic at the level of chains. Momentum is shared, kinetic energy
  becomes rest mass (E = mc² directly on the substrate), and the inertia of the merged flow is its energy. Matter in
  this model cannot bounce by contact: an elastic channel would have to come from fields or from a bound structure,
  which the model does not have (open).
  - **d4 (3+1, 16 sprinklings, ~81 500 events, sources 227 + 221, PASSED 29 Oct):** f_mix = 0.774 ± 0.008
    (0.72–0.85), v_mix = −0.003 ± 0.011, pure flows ±0.600, bounce 0, rest mass/parts 1.216; without contact exactly
    0; frame and four-momentum (2e−16) OK. The same behaviour as in 2+1, so not an artefact of the dimension.

- **F6.8c2: the third coefficient also arises by itself: 𝒞(0) = 2.002 (EH 2); all three coefficients of the static EH
  form (6, 8, 2) come out of counting intervals, and γ(0) = 0.997 to 1.000 follows from the action of the order**
  (chain, 3+1, PASSED 29 Oct; revision of F6.8c, reading at fixed n = kT on new k = 0.06 / 0.085 / 0.12 / 0.17,
  n = 6 / 8 / 10, T ≤ 167).
  - **Instruments I1–I5 OK:** ⟨T2⟩ ≤ 7.9e−5; V'' 6.6e−4; uniform null 3.6e−15; conformal K̂ against F6.8a 1.6e−10; nulls
    k = 0 ≤ 2.9e−4 in units of K̂_bb; level 2 ≤ 6.7e−4 (k = 0.06, T = 133); the method on the conformal sector
    K0 = 0.99976 (F6.8a 1.00005).
  - **K0_bb:** 0.99971 / **1.00118** / 0.97108 at n = 6 / 8 / 10 (criterion at n = 8). Conformal by the same method
    0.99976 / 0.99976 / 0.99977. At n = 10 the result is pulled by the oldest row (k = 0.06, T = 167: 0.9188 against
    0.9366 at n = 8), where level 2 also changes most. Together with F6.8c (a break of the scaling at T ≳ 170) this
    points to a numerical limit at large age; undecided.
  - **Reported:** the ageing of the optical sector grows faster than linearly. Increments D(n) = (K0 − K̂)/k² at
    k = 0.12: 13.90 / 14.28 / 15.02 / 18.00 / 27.95 for n = 6 / 8 / 10 / 16 / 24 (the increments grow). The data do
    not fix the form, however: a local fit s n + q n² gives different coefficients depending on the range of n
    (n = 6–10: q ≈ 0.044, s ≈ −0.4; n = 8–16 at k = 0.14–0.2: q ≈ 0.019, s ≈ 0). The conformal sector does the same at
    n = 6–10 (s ≈ −0.28, q ≈ 0.023), and becomes linear only from n ≈ 8 (c_T 0.110; confirmed analytically). The
    earlier description "≈ 0.02 k⁴T² with no linear term" was overstated and is corrected in the paper too. The row
    n = 24 has T = 200, i.e. in the region of suspect numerics.

  Reading: γ = 1 at linear order follows from the discrete BD action of the order in all three static sectors, nothing
  of GR put in. The prediction of the age term (the far-field law, the window of ℓ_k) now depends on the asymptotics of
  the optical sector for kT → ∞ (F6.8e, analytically like c_T).

- **F6.8b2: the mixed coefficient also arises by itself: c_mix = 7.983 (EH 8); in the local part the mixed sector
  stands in the covariant ratio to the conformal one, beyond the leading order the sectors differ** (chain, 3+1,
  deterministic expectation value, PASSED 28–29 Oct; ερ = 1, ℓ_k = 1; revision of F6.8b, whose instruments failed).
  - **What is computed:** ∂²/∂a∂b of the mean BD action in g = e^{2af}(−dt² + e^{2bh}dx²), f = h = cos kz averaged over
    the phase. The optical index changes the order: to first order in b through the optical length of a straight
    segment (Fermat). The volume of an interval through the spatial projection of the diamond (a prolate spheroid),
    the difference of the delays formed pointwise, the boundary term of the age cut on its own quadrature. The
    conformal K̂ on the same nodes and ages. The number 8 is only in the denominator of K̂_ab.
  - **Instruments I1–I4 OK:** geometry against a spherical quadrature 3.5e−5 / 5.8e−5, uniform null 8.6e−14; conformal
    K̂ against F6.8a 1.6e−10; null k = 0 at T = 40 / 160 / 320: +6.9e−9 / −1.2e−7 / −8.0e−7, in units of K̂ at k = 0.05
    2.6e−4 (tolerance 2e−3); level 2 against 1 ≤ 3.2e−5.
  - **Criterion:** max_{k≤0.1}|Δ| = 4.6e−3 (PASSED ≤ 0.01). Ages 8 / 12 / 16 × 1/k.
    K_loc,ab = 0.9683 / 0.9464 / 0.9045 / 0.8387 / 0.7261 against the conformal 0.9681 / 0.9441 / 0.8999 / 0.8297 /
    0.7076 for k = 0.05 / 0.07 / 0.1 / 0.14 / 0.2. Fit over k ≤ 0.14: **K0_ab = 0.9979, c_mix = 8 K0_ab = 7.983**,
    c1 = 3.55, c0 = 1.15, rms 2.2e−4.
  - **Reported, not evaluated:**
    - Δ ≈ 0.46 k² purely analytically for k = 0.07–0.2 (Δ/k² = 0.475 / 0.462 / 0.462 / 0.461).
    - The age term in the mixed sector c_T,ab ≈ 0.195 against the conformal 0.110 (ratio 1.78).
    - k = 0.05 (the longest age of the grid) deviates from both regularities: Δ = 1.6e−4, c_T,ab = 0.155. It has no
      effect on the verdict (k = 0.1 decides), not read.
    - k = 0.3: K̂_ab changes sign and the linear model of the age does not hold (residual 1.7e−2).

  Reading: two of the three coefficients of the EH quadratic form (6 and 8) arise from counting intervals and stand in
  the covariant ratio at leading order. A covariant quadratic action (local or non-local) gives B/8 = A/6 for every k,
  so Δ ∝ k² and the differing age terms are not covariant. The mean BD action is exactly diffeomorphism invariant, but
  a static field cut at the age T is not a covariant configuration (the cut fixes the rest frame of the field). The
  differences may therefore come from the cut, not from the action; which of these holds is open (F6.8d). The field
  of a point mass is therefore composed directly from A, B, C: γ = B/(2C) − 1, Φ/Φ_GR = −2/(A − B²/4C). A and B are
  computed, C is missing (F6.8c).

- **F6.8a: the discrete BD action itself carries the Einstein–Hilbert conformal stiffness: c_conf = 6.0003 (EH 6),
  nothing of EH enters the computation** (chain, 3+1, deterministic expectation value, PASSED 28 Oct; ερ = 1,
  ℓ_k = 1).
  - **What is computed:** the second derivative of the mean BD action with respect to a in the metric e^{2af}η,
    analytically from the kernel gk and the proper volumes of intervals; f = cos(kz) averaged over the phase; the exact
    Fourier transform of the diamond. The number 6 is only in the denominator of the comparison K̂.
  - **Instruments I1–I5 OK:** transform against quadrature of the kernel 6.5e−8; derivative against finite differences
    of the kernel operator 3e−4; phase average 2.5e−11; null k = 0: |J0/J| ≤ 5.4e−5; level 2 against 1 ≤ 3.3e−5.
  - **Local part:** K_loc = 0.9904 / 0.9826 / 0.9680 / 0.9438 / 0.8994 / 0.8287 / 0.7056 / 0.4687 / −0.0152 for
    k = 0.025 / 0.035 / 0.05 / 0.07 / 0.1 / 0.14 / 0.2 / 0.3 / 0.5. Fit K_loc = K0 − k²(c1 ln(1/k) + c0) over
    k ≤ 0.1: **K0 = 1.00005**, c1 = 3.98, c0 = 0.91, rms 2.3e−5.
  - **The correction for a finite ℓ_k is non-analytic:** 1 − K_loc ≈ (kℓ_k)²(4.0 ln(1/kℓ_k) + 0.9). The stiffness
    changes sign at k ≈ 0.49/ℓ_k (wavelength 12.7 ℓ_k).
  - **The age term (non-local, new):** K̂(k, T) = K_loc(k) − c_T k³T, where T is the time for which the static field
    has existed. Linear in T (residual ≤ 1e−3), c_T = 0.098–0.107 for all nine k, for k ≤ 0.2 c_T = 0.1035 ± 0.0003.
    Predicted analytically before the computation (long intervals perpendicular to k in phase with the field; the BD
    moments kill it up to τ⁴, not at τ⁶); the asymptotic integral 0.1012–0.1022 for k = 0.05–0.2 (agreement to 2 %,
    the rest is the finite T). For k ≤ 0.035 the asymptotic integral is numerically unreliable (0.147 and 0.088: the
    integrand is ~k⁶ and is lost in round-off); a descriptive comparison, it does not affect the verdict.

  Reading: the local part of the discrete BD action is exactly EH in the static conformal sector. With the
  diffeomorphism invariance of the mean action and locality in the other sectors, γ = 1 follows (the rule of F7.4a
  derived, not put in); locality in the mixed sector was confirmed by F6.8b2 (at leading order). The relative size of
  the age term is 0.1 k³ℓ_k²T: unobservable for a Planckian ℓ_k; for a smeared action with a large ℓ_k a bound
  ℓ_k ≲ √(0.1 L³/T) (deviation below 10⁻²). What T is in a real situation is open.

- **F7.4a: when every flow narrows the cone in its own rest frame and the capacity deficits add, two counter-streaming
  flows gravitate according to T_ab: the clock field grows as Γ_v²(1 + v²) and the cones are anisotropic**
  (mechanism, 3+1, 16 sprinklings per configuration, ~55 700 events, PASSED 28 Oct; decision 13, path A).
  - **Instrument:** flat control OK, static amplitude −0.017 ± 0.020 against the prediction, isotropic cones read
    isotropic (K = 1.002 ± 0.004).
  - **Amplitude at ±0.6c:** R = A(0.6)/A(0) = 2.401 ± 0.053. The exact prediction of the rule is 2.361 (linearly
    2.125). The alternatives are far away: energy in the mean frame 1.5625 (16 SE), laboratory density 1.25 and rest
    density 1. The pressure of the streams (T_00 + T_ii) therefore gravitates too, not only the energy.
  - **Anisotropy:** K = along/across = 1.221 ± 0.012 against K_T 1.192. Every scalar rule of the mean frame gives 1.
    The slight excess (+0.030 ± 0.012) is within the tolerance.
  - **Reported (one stream at 0.6c):** R = 2.538 (prediction 2.519). The cone is tilted, the deficit against the
    stream is 17.6× larger than along the stream (prediction 12, scalar rules 1): the substrate drags the cones with
    the flow. The deficit along the stream fits exactly (−0.065 against −0.065). The deficit against the stream is,
    at a strong narrowing (c ≈ 0.46), read beyond the measure (−1.24 against −0.79). Most likely a limit of the cone
    reading at strong narrowing; only reported.
  - Reading: the rule that path A needs works on a causal set. The tensor structure arises from the response of the
    order without an inserted field equation, Green function or trace reversal. The test does not say which rule nature
    uses (that pressure and motion gravitate is known from observation). It says that a substrate with this rule gives
    the coupling of general relativity at linear order, including frame dragging. The rule itself remains an input
    (own frame of the flow and addition of the capacity deficits, book 5.3, F3.6).
- **K.16c: the blind instrument of γ works also in 3+1, in a box that 3+1 can carry** (calibration, 16 sprinklings
  per configuration, ~52 000 events, box t [0, 9], |x| <= 9, PASSED; null control OK). The 3+1 dictionary:
  density → φ − 3ψ, cones → φ + ψ, diamond πτ⁴/24 → φ. Field φ = −0.1·exp(−r²/32): γ = −1 reads −1.001 ± 0.001 /
  −1.001 ± 0.001; γ = 1/3 (the alternative "narrower cones only") 0.325 ± 0.037 / 0.249 ± 0.041; γ = +1
  0.994 ± 0.049 / 1.013 ± 0.061. The instrument for F7.1c.

- **K.16b: the blind instrument of γ reads clocks, cones and density from the order and returns the γ = −1, 0, +1 put
  in** (calibration, 2+1, 32 sprinklings per configuration, ~18 000 events, PASSED; null control OK). The instrument
  receives only the events, the relation of the substrate and the proper density 1, not the metric. Dictionary: clocks
  (number of elements of a diamond at rest) → φ, boundary of the relation (speed of light) → φ + ψ, density per
  coordinate volume → φ − 2ψ (in 3+1 φ − 3ψ). Field φ = −0.1·exp(−r²/128), ψ = γφ:
  - γ = −1: A (cones + density) −1.000 ± 0.001, B (clocks + cones) −1.000 ± 0.001. Unchanged cones by themselves mean
    γ = −1, hence exact;
  - γ = 0: −0.039 ± 0.019 and −0.068 ± 0.021;
  - γ = +1: +0.921 ± 0.026 and +0.891 ± 0.044 (a residual shift towards zero ~0.1).

  Range: at a = 0.1 the instrument is within the tolerance 0.25. In a weaker field a = 0.03 (only reported) it reads
  γ = +1 as 0.56 ± 0.07 and 0.27 ± 0.08, so for |φ| ≲ 0.03 it is unreliable. Candidate couplings Π → relation must
  therefore be measured in a field |φ| ~ 0.1. The instrument uses coordinates and the known density (a declared
  readout); it serves substrates built on an embedding. A coordinate-free reading (radar echo) is the next step.
  The previous version K.16 (16 sprinklings) narrowly FAILED: B for γ = +1 1.092 ± 0.080, |deviation| + 2 SE = 0.253
  against 0.25. A diagnosis on the same sprinklings (it decided nothing) showed that this is not regression dilution
  but two shifts of the reading: the midpoint of the gap at each anchor reads the speed of light ~0.005 lower even
  without a field, and a diamond of height 5 averages the field over a width of 3.75. K.16b is a revision of the
  instrument with new sprinklings and the same criteria: a pooled boundary of the relation and diamonds of height 3.

- **F3.10b law: the rule transmits the partner's field as a force with the right dependence on distance** (mechanism,
  2+1, d = 8 / 12 / 16 / 20, 64 sprinklings per d, window of linear response [1.25 d; 2.5 d], PASSED 23 Oct, seed 2;
  replication missing). The momentum received by the pair (B − A) against the first-order prediction from branch A
  (field of the real partner, bundles in isolation, no data from B): exponent k_meas = 0.845, k_pred = 0.946,
  Δk = −0.101 ± 0.051, |Δk| + 2 SE = 0.203 ≤ 1/4. A continuous static segment in the same windows gives 0.927 (the
  kernel alone 1; the difference comes from the retarded build-up of the field in the window). Ratios measured /
  predicted: 0.931 ± 0.018 / 0.908 ± 0.030 / 0.964 ± 0.038 / 1.021 ± 0.051; shift of the separation in the window
  0.4–2.1, far below d/4. Instrument on all 256 sprinklings: C = A, mirror and frame bitwise, 0 ties. The reading is
  limited:
  - the kernel 1/(2πτ) and the sign of the coupling are inputs; the test shows that the discrete dynamics of the field
    converts it into a force faithfully (to ~10 %), not that the 1/r law follows from the order;
  - at d = 8 and 12 the force is 7–9 % smaller than the prediction (2–4 SE), at d = 16 and 20 it fits; cause unknown
    (estimate of the gradient from the parents against a central difference, second order of the shift);
  - it explains the trend of the old part F3.10 law (0.61 → 0.98): that came mainly from the faulty quasi-static
    prediction.
  The part force (ρ_d = 1 for every d) has its own verdict and with the same seed computes the same data.
- **F3.10 attraction: two separate bundles attract each other through their retarded fields** (mechanism, 2+1,
  32 sprinklings ~41 500 events, d = 12, PASSED, seed 0; replication with seed 8 PASSED). Relative acceleration after
  subtracting the drift of isolated bundles a = −0.0080 ± 0.0011. The prediction from the retarded field of a static
  segment (ε = 0.1 locked from F3.8/F3.9, λ and μ from the calibration of a single bundle) is −0.0111, ratio 0.72.
  Instrument: C = A choice by choice (coupling only through the field), the recomputation in a frame β' = 0.6 and the
  mirror x → −x choose identically. No new mechanism: the dynamics is the field from F3.8. The reading is limited:
  - the attraction follows from the rule bending the carried direction towards slower clocks while the partner's
    field slows the clocks; the kernel 1/(2πτ) and the sign of the coupling are inputs (F3.8), not outputs;
  - the attraction alone says nothing yet about the force law, momentum or universality (separate parts below and in
    the plan).
  Replication with seed 8 (23 Oct, PASSED; independent sprinklings and calibration): a = −0.0092 ± 0.0014
  (a + 2 SE = −0.0064), accelerations from momentum +0.0038 and −0.0041 (seed 0: +0.0037 and −0.0036),
  Ṗ_L + Ṗ_R = −0.051 ± 0.042; instrument OK, 0 ties. The independent calibration gave φ_ref 26.9 instead of 24.7,
  i.e. λ 9 % smaller (0.00371): the calibration noise of λ between seeds is ~9 %. It does not matter for the sign of
  the attraction nor for F3.10b (the same λ in the simulation and the prediction); ratios to the quasi-static
  approximation (here 0.92) are therefore only indicative for this reason too.

  Part momentum (d = 12, flows 1 and 2, window of linear response [15, 30], 32 sprinklings, PASSED, seed 0): the bundles
  exchange opposite and equal momenta. Ṗ_L = +1.552 ± 0.091, Ṗ_R = −1.545 ± 0.120, sum +0.007 ± 0.158 against a
  tolerance of 0.387; the resolution of the balance is thus ~20 % of a single change. Accelerations from momentum
  +0.0094 and −0.0047: a bundle twice as heavy accelerates at half the rate (this follows from Ṗ_L = −Ṗ_R and M ∝ flow,
  nothing more). The reading is limited:
  - the coupling is symmetric (the same λ, the same kernel, each bundle is a source and a receiver in proportion to its
    flow), so action and reaction are built into the rule for static sources; the test shows that retardation,
    discreteness and losses at merging do not break it in linear response. The field carries no momentum; at faster
    motion retardation could break the balance; this cannot be resolved here;
  - the separation shortens in the window by 1.49 ± 0.18 (linear response holds); the quasi-static approximation gives
    the ratios 1.06 (centroid) and 0.84 (momentum L) here, without a verdict.
  Replication with seed 8 (23 Oct): UNDECIDED, hence **not confirmed**. Ṗ_L = +1.467 ± 0.056, Ṗ_R = −1.543 ± 0.152,
  sum −0.076 ± 0.179; |sum| + 2 SE = 0.434 above the tolerance 0.376, |sum| − 2 SE below it. The balance is compatible
  with zero in both runs (+0.007 ± 0.158 and −0.076 ± 0.179), but the test lacks the power for it: the error of Ṗ_R
  alone is almost half the tolerance. Seed 0 passed narrowly (0.323 against 0.387). A decision would need a new run with
  a larger number of sprinklings fixed in advance; the criteria do not change (R132).

  Part universality (d = 12, an external source R with flow 2 that responds only to its own field, L with flow 1 and 2
  on the same sprinklings, window of linear response, 32 sprinklings, PASSED, seed 0): acceleration from the carried
  momentum a_P,L = +0.00919 ± 0.00046 (flow 1) and +0.00904 ± 0.00062 (flow 2), paired difference
  +0.00014 ± 0.00066 against a tolerance of 0.00152. The centroids agree (+0.0092 ± 0.0026 and +0.0105 ± 0.0022). The
  received momentum doubles with the flow (1.55 → 3.02), the acceleration stays. The external field is bitwise the same
  in both configurations. The reading is limited: the bending does not depend on the flow of the carrier, so the
  universality of the response is built into the rule; the test shows that the own field (stronger for flow 2),
  merging and discreteness do not break it. It is not an independent confirmation of the equivalence principle. The
  quasi-static approximation (it assumes a moving R, here R stands still, so the ratios are not comparable) gives 0.82
  and 0.78, without a verdict; a_P,L agrees within the error with part momentum (R responds), 0.0094.
  Replication with seed 8 (23 Oct): UNDECIDED, hence **not confirmed**. a_P,L = +0.00908 ± 0.00029 (flow 1) and
  +0.00980 ± 0.00052 (flow 2), paired difference −0.00072 ± 0.00052 against a tolerance of 0.00157;
  |difference| + 2 SE = 0.00176. In seed 0 the difference was +0.00014 ± 0.00066 (the opposite sign), so together there
  is no systematic difference, but the replication did not meet the criterion. Flow 1 gave the same acceleration in
  both runs (0.00919 and 0.00908).

  Part law (d = 8 / 12 / 16 / 20, 24 sprinklings per d) ended UNDECIDED. All d attract; the centroid exponent is
  k = 0.77 ± 0.57 against k_pred = 1.11 (the same fit of the prediction), the centroid noise at d = 16 and 20 is
  36–45 %. Diagnosis without a verdict: the acceleration from momentum Ṗ/M is roughly 5× more precise. Its ratio to the
  point prediction grows with the distance (0.61 / 0.76 / 0.91 / 0.98, mean of L and R) and the exponent is ~0.58
  instead of 1.11. Momentum balance of a symmetric pair: Ṗ_L + Ṗ_R is within ~10 % of a single change at all d.
  Run with seed 8 (23 Oct, superseded part, for the record only): UNDECIDED, k = 1.51 ± 0.49 against the approximation
  1.12; the centroid at d = 16 is unresolved (−0.0011 ± 0.0010). Acceleration from momentum 0.0052 / 0.0042 / 0.0032 /
  0.0029 (seed 0: 0.0058 / 0.0042 / 0.0037 / 0.0034) – the dependence reproduces, flatter than the approximation;
  F3.10b decides.

  The author's QA review (23 Oct, after the runs of attraction and law) found a fault of the prediction, not of the
  dynamics. The quasi-static prediction inserts the current distance into the formula of a static source and treats
  the bundle as a segment. With the calibrated λ and μ, point segments at d = 8 pass through each other in the window
  (smallest r = −4.7) and at d = 12 the separation shortens from 12 to 5.5. The ratios a/a_pred and the trend of the
  diagnosis are therefore distorted by the prediction itself; the explanation by an "extended source" is only one of
  the possibilities. The verdict of law stays UNDECIDED (R132), law is not replicated and is replaced by F3.10b: a
  first-order prediction from the real discrete field of the real partner in the window of linear response.
  Further findings of the QA:
  - the model conserves the scalar flow, not the four-momentum: when different directions merge, |P| > F and the excess
    is discarded; measured ~1–3·10⁻⁴ of the flow, ~0.03–0.05 % of the x-momentum per event, negligible against the
    effects. The field carries no momentum;
  - the processing order is a linear extension of the causal order (now verified in the code); a boost preserves the
    order, so the frame check compares the same computation;
  - the checks of C, mirror and frame now compare also the flow and directions, a difference starting with a numerical
    tie is reported separately; a deterministic part of unit was added (λ = 0, flow balance, sign of the kick). unit
    PASSED 23 Oct (4 sprinklings, all OK, 0 ties; 93–96 % of the kicked events point towards the partner; excess of
    four-momentum at merging 1·10⁻⁵ to 3·10⁻⁴ of the flow); again on independent sprinklings (seed 8) PASSED 23 Oct
    (all OK, 0 ties, 87–97 % of the kicked events towards the partner, excess 8·10⁻⁵ to 1.7·10⁻⁴).
- **F3.6 (rest and move): a flow that carries its timelike direction as a state stays at the world line of the source
  and moves with it (inertia); a direction taken from the last step does not** (mechanism, 32 sprinklings 2+1,
  ~32 800 events, source of 204 events, PASSED). Expansion speed R_rms in the window [5, 11] in multiples of c:
  uniform splitting (F3 series) 0.770 ± 0.003, direction from the last step 0.393 ± 0.007, carried direction
  **0.054 ± 0.006** (localisation tolerance 1/4). The carried bundle holds the whole flow (~300 in each layer), its
  width grows slowly from 3.0 to 4.0 over 24 ℓ (roughly diffusively). The carriers merge from 73 to 27 in one layer:
  the bundle contracts to a smaller number of world lines. The recomputation in a reference frame β' = 0.6 chose the
  same at every event (the rule does not pick a frame).

  The kick of one link is ~0.45–0.52 rapidity. Momentum therefore cannot be a displacement between neighbouring events
  (such a memory decays and the flow spreads at a speed of ~0.4), it must be an internal state of the flow. The reading
  is limited:
  - the rule uses geometry at the scale of a link, not only order and counts (criterion (a) of the audit), and the
    choice of the descendant needs the future (b), like the whole F3 series;
  - localisation does not mean a bound object: the bundle holds because the flow is conserved and carries a direction.
    Nothing binds its width.

  Part move (source β = 0.5, ~49 200 events, PASSED):
  - the carried direction is localised (v_s 0.050 ± 0.006) and stands still in the frame of the source: centroid
    velocity +0.0075 ± 0.0043 against the tolerance β/4 = 0.125 and the alternative "stays at rest in the laboratory"
    −0.5;
  - the width evolves the same as with a source at rest (2.94 → 3.74), as the principle of relativity requires;
  - the direction from the last step spreads (0.349 ± 0.012), uniform splitting too (0.782).

  The hypothesis "localisation in an invariant substrate requires a memory of direction (momentum)" thus has support,
  and the carried direction also gives inertia. Replication of rest with seed 1 (23 Oct, PASSED): carried direction
  0.0534 ± 0.0073, step 0.367 ± 0.012, uniform splitting 0.774 ± 0.004; width 2.94 → 3.92, carriers 73 → 26 in one
  layer, the frame β' = 0.6 chooses identically – agrees with seed 0 in all numbers. Replication of move with seed 1
  (23 Oct, PASSED): carried direction localised (v_s 0.0528 ± 0.0056), stands still in the frame of the source
  (v_c +0.0045 ± 0.0061 against the tolerance 0.125 and the alternative −0.5), width 2.89 → 3.77; step 0.369, uniform
  splitting 0.777. F3.6 is thereby closed: the first result of the laboratory that was not put into the rule,
  confirmed on independent data.
- **F1.7 V and W: the solver of the equation B_g f(Π) = κρ_s (3+1, H1, continuum) is verified, and in a weak field all
  candidates with f(1) = 0 give a Newtonian response, equation (a) none** (mechanism, PASSED 23 Oct). V: the symbolic
  forms reproduce C5 and λ₁(γ) = 2 − γ(d_s − 1); solver against the analytic linear solution 3.3·10⁻⁹; ln Π at σ = 4
  and 16 gives the same profile Π(r/σ) to 1.1·10⁻¹² (Π(0) = 0.6253 at u_lin(0) = −0.3). W (u_lin(0) = −10⁻⁴): ln Π,
  Π − 1, (Π⁻¹ − 1)/(−1), (Π² − 1)/2, (Π³ − 1)/3 differ from the linear solution by 2·10⁻⁵ to 1.8·10⁻⁴ (limit 10⁻²);
  the negative control f = Π has a linear coefficient of exactly 0. As expected this does not distinguish the
  candidates; N decides. γ = 1 is an input of the ansatz here, not a result.
- **F1.7 N: strong field – every candidate has a regular static solution only up to a finite source strength; ln Π is
  the second best of five** (mechanism, RECORDED 23 Oct; 3+1, H1, continuum, Gaussian source, strength as u_lin(0)).
  Last strength with a regular solution / first without one:
  - (Π⁻¹ − 1)/(−1): −0.5 / −0.7; ln Π: −0.3 / −0.5; Π − 1 and (Π² − 1)/2: −0.2 / −0.3; (Π³ − 1)/3: −0.1 / −0.2;
  - the order is monotonic in k of the family (Π^k − 1)/k (ln Π is k → 0): smaller k, wider range;
  - beside the weak-field branch a second, strong-field branch appears (Π(0) ≈ 0.09–0.40); both approach each other
    towards the limiting strength (for k = −1 at −0.5: 0.41 and 0.28), which looks like a fold (two solutions merge
    and vanish);
  - the external charge grows above the linear value (ln Π: 1.02 / 1.08 / 1.22 / 1.51 at −0.03 / −0.1 / −0.2 / −0.3):
    the non-linear field itself contributes to the mass seen from outside.
  Reading: a finite limit is not a fault in itself (general relativity also has a compactness limit of static bodies),
  but the source has no pressure here and the conversion of u_lin(0) to the compactness GM/(Rc²) has not been done, so
  a comparison with neutron stars is not yet possible. The candidates are distinguished only by the second order of the
  external field (Derived consequences: β). Limitation: γ = 1 is an input.
- **K.15 (all four parts): the tensor of the carried flow T^ab = (1/ΔV₄) Σ F u^a u^b τ_link measures dust, counter-
  streaming beams and an isotropic gas correctly, and with the merge rule F = |P| the 4-momentum is conserved exactly**
  (beam: calibration, 3+1, 64 sprinklings ~28 700 events, PASSED 23 Oct). In the rest frame of the beam (Landau frame)
  p = 0 and π = 0, in the laboratory E/ρ = γ² = 1.5625, |q|/E = v = 0.6, S_xx/E = v² = 0.36 and p/E = 0.12 in each cell
  2 and 3 without scatter; frame check (tensor from boosted coordinates) to 1e-9. Both merge rules give the same,
  because equal directions merge without loss. It verifies the algebra of the decomposition (E, q, p, π).
  counter (two counter-streaming populations ±v, PASSED 23 Oct): q_x/E = −0.001 ± 0.009 (0), S_yy = 0 exactly,
  p/E = 0.1127 ± 0.0006, S_xx/E = 0.338, π_xx/E = 0.225 – all ~6 % below the value of the emitted directions
  (0.12 / 0.36 / 0.24), within the tolerance. The deficit is explained by mixing: 9.1 % of the flow in the window no
  longer has the velocity v (rule F = |P|; summation 7.5 %), because counter-streaming chains merged. Inelastic merging
  turns directed motion into rest mass (with F = |P| energy is conserved, pressure drops): the population cools. Cells 2
  and 3 give p/E 0.109 and 0.113 with a scatter of 0.025 and 0.017.
  iso (isotropic gas, PASSED 23 Oct): q_x, q_y, q_z / E = +0.003 / +0.000 / +0.009 (± 0.005), π_xx, π_yy, π_xy / E
  below 0.0012 ± 0.0017 (isotropy), p/E = 0.1119 ± 0.0005 against 0.12 (7 % lower; mixing 18.7 % of the flow, merging
  7.9 %).
  conserve (isotropic gas, PASSED 23 Oct): with the merge rule F = |P| the 4-momentum on the cuts t = 2.5–4.5 is
  conserved to 6·10⁻¹⁷ (exactly); today's rule (flow = sum of flows) loses on average 0.67 % of energy (at most 1.9 %)
  in 2 units of time at 7 % merged flow. No flow in the window ends at a wall. Consequence: today's transfer rule (F3
  series) violates the conservation of 4-momentum whenever different directions merge; the derived rule F = |P|
  satisfies it exactly. Whether to adopt it as a rule of the laboratory is the author's decision.
  Replication with seed 2 (23 Oct, all four parts PASSED): beam exactly; counter p/E 0.1132 ± 0.0006 (seed 0: 0.1127),
  q_x/E −0.015 ± 0.009, mixing 8.2 %; iso p/E 0.1114 ± 0.0005 (0.1119), components of q and π to 0.004 ± 0.005,
  mixing 19.8 %; conserve F = |P| 7.6·10⁻¹⁷, the summation rule loses 0.81 % of energy (at most 2.0 %).
  K.15 is closed.
- **K.13a: the discrete BD operator in the H1 geometry (γ = 1) converges at first order to the continuum limit
  □_g − R/2, the density and the cone part each separately** (calibration, 3+1, all three parts PASSED 23 Oct;
  deterministic quadrature, replication with a seed is meaningless). lin, ratio to the limit at σ = 4 / 8 / 16 / 32 ℓ_k:
  - w = 1, r = 0: 0.743 / 0.894 / 0.961 / 0.987; r = σ: 0.806 → 0.992; r = 2σ: 0.537 → 0.977;
  - Gaussian w: r = 0: 0.596 → 0.977; r = σ: 0.689 → 0.986; r = 2σ: 0.700 → 0.986;
  - the deviation falls sharply in all six rows, local exponent log2 1.56–1.65 (slower than (ℓ_k/σ)²); at σ = 32 it is
    0.8–2.4 % against the limit d_alt/4 (12.5–57.5 %); main quadrature against control ≤ 0.02 (limit d_alt/16);
  - the density part P and the cone part Q (response to a change of the speed of the cones, i.e. to the spatial part
    of the metric) converge each separately; the coefficient of ∇²u for w = 1 goes to −1 (the term −R/2 = −∇²u is
    reproduced).
  The substrate therefore carries also the cone part, on which γ depends: a geometry with bending of light can be
  expressed on a causal set. For the consequence for equation (a) see Derived consequences. geom and valid (σ = 4):
  - optical time against the optimised path converges quadratically (log2 ratios 1.97–2.00 at γ = 1 and 0, limit
    1.75);
  - volumes of intervals with curved cones against Monte Carlo (20 diamonds, 400 000 samples): max |z| 2.73;
  - thin intervals against mpmath: 40/40, worst relative error 2.1·10⁻¹⁰.
  valid:
  - Π = 1 against K.4 at a relative level of ~10⁻¹² (rows r = 0, 4, 8);
  - γ = −1 against the central difference of the instrument K.11: differences 2·10⁻⁷ to 2.4·10⁻⁵ against tolerances
    ~10⁻³–10⁻².
  At σ = 4 the discrete operator gives 0.74 (w = 1) and 0.66 (Gaussian w) of the continuum values on the axis: this is
  the systematics ℓ_k/σ, whose decrease with σ is tested by lin.
- **K.14: the number of molecules of a horizon measures the area of the horizon, in a curved geometry the proper area**
  (calibration, all 5 parts PASSED; the entropy instrument for the Jacobson route; the construction is taken from
  Barton et al. 2019, not a claim of the theory itself). Flat space against the constant a^(d) from an independent
  quadrature, per unit area in units of ℓ^(d−2):
  - 1+1: 0.3325 ± 0.0043 and 0.3241 ± 0.0043 at ρ = 1 and 4 (a = 1/3; 20 000 sprinklings);
  - 2+1: 0.2182 ± 0.0089 and 0.3612 ± 0.0127 (a ρ^(1/3) = 0.2189 and 0.3474);
  - 3+1: 0.1716 ± 0.0014 and 0.3407 ± 0.0027 (a ρ^(1/2) = 0.1732 and 0.3464; a^(4) = √3/10).

  All |z| <= 2.2. The two z = −2.17 at ρ = 4 were not confirmed as systematic by a control run in the container (1+1,
  ρ = 4, 60 000 sprinklings: 0.3329 ± 0.0025, z = −0.17). Conformally flat space Π(y) = 0.75 + 0.25 cos(2πy/40)
  (Π 0.5–1): exponent of the proper area Π^k
  - 2+1: k = 1.060 ± 0.065 (target 1), C/a = 0.992, model without free parameters χ² 18.7 / 16 bands, coordinate area
    1151, volume 778;
  - 3+1: k = 1.964 ± 0.055 (target 2), C/a = 0.959, χ² 26.0 / 16, coordinate area 11 746, volume 1 261.

  The first assumption of Jacobson's argument (entropy proportional to the proper area with a universal constant) is
  thus supplied by the causal set at least where the density changes and the cones stay flat. Curved cones, matter with
  momentum and the Clausius relation remain (docs/derivations/jacobson_route.md). In 3+1 the amplitude is 4 % below the
  leading order; it may be the finite scale ℓ/λ, another ratio λ would decide.
- **F6.3b V1 and V2: the exact construction of the dynamic order is verified in vacuum and in the limits** (validation,
  both PASSED). V1: exact vacuum = links of the sprinkling in 4 sprinklings N = 36 000 (no link at a wall occurred this
  time), brute force with variable Π without a difference in 6 sprinklings. V2: at Π_sat = 1 the solution ends after
  the 1st iteration and is element by element equal to the fixed order in all 4 sprinklings; windows 1.5 and 3.0 give
  identical links on the solution with matter (~1.95–1.97 million, vacuum ~1.52 million: matter adds ~28 % of links).
  The fraction of events whose Π changes between iterations falls from ~0.5 to 0.26 in 9 iterations (first
  sprinkling). Part A may run.
- **F6.3 V1: the order from cones dependent on the clocks is exactly a causal set in vacuum** (validation, PASSED). At
  Π = 1 the exact construction (transitive closure of the local relation, Hasse diagram) gives links element by element
  identical to the links of the sprinkling in 4 sprinklings N = 36 000 (~1.52 million links each); with variable Π (a
  clock well down to 0.5) it agrees with brute force (full closure) in 6 sprinklings of 1 200 events without a single
  difference (~20 000 links each); all links go forward in time (acyclicity).
- **The field of a localised source is far away the Green function with a strength according to Gauss's law** (F1.3,
  mechanism, expected non-local operator, window 10–40 l_k): 3+1 g = −1.02 (1/r), strength 1.0045 of the prediction
  Q/4π; 2+1 g = −0.08 (ln r), strength 1.055 of the prediction Q/2π. Two resolutions agree. The relative deviation from
  the continuum solution grows towards the edge (3+1: 2 %, 7 %, 20 % at r = 10, 20, 40), because the continuum u → 0 at
  the wall and the field carries a constant shift from the edge; neither the shape nor the strength is affected.
- **Integral A5′: the source creates a field according to a pre-computed operator** (F1.2b, mechanism, N = 32 000, one
  sprinkling, smooth Gaussian source σ = 0.05). It holds in 2+1; in 3+1 the result is UNDECIDED (see below). The
  averaged operator applied to the field returns on average the source at all six averaging scales (max |z| 1.32 in
  2+1, 1.16 in 3+1), and the error falls with the same exponent as for a field without a source (p_test / p_cal =
  0.63 / 0.48 in 2+1, 0.76 / 0.81 in 3+1). The null model with permuted kernel weights failed (z up to 19.7 in 2+1). In
  3+1 seed 0 passed only narrowly (the null model failed at a single scale, z = 3.23); the replication with seed 1 ended
  UNDECIDED: fields A and B both satisfied it (max |z| 1.51; p_test 0.51 ± 0.21 against p_cal 0.44 ± 0.15), but the
  null model satisfied A too (max |z| 2.64). In 3+1 at N = 32 000 (L/ℓ = 1.6–7.1) one order therefore does not reliably
  distinguish the field from permuted weights; this is a lack of resolution, not a contradiction with A5′. The previous
  version with a sharp source (F1.2) ended as an instrument failure, see below.
- **Gauss's law for the flux of a persistent source** (F1.1, E45a): 2+1 and 3+1, N up to 256 000,
  J/S_w = 1.003–1.004; a model without conservation of the flux 0.77–0.88.
- **Paths follow the geodesic of the metric Π²η** (F2.2, F2.4; E46, E46b, replicated). In the conformal representation
  every maximal chain follows it; in the weighted representation the form of the action decides
  (Σ Π² → Π⁴η: 0.0770 against 0.0776).
- **In 3+1 the form of the action decides the path just as in 2+1** (F2.5 d4, mechanism, weighted representation,
  inserted field of a soft point source, 256 seeds per size, PASSED). Rule E46b part B at N0 = 64 000: SPECIFIC.
  H1 (Σ Π) = 0.0524 ± 0.0026 against P1 = 0.0432 (geodesic of Π²η, band ± 0.0108); H2 (Σ Π²) = 0.0799 ± 0.0025
  against P2 = 0.0776 (geodesic of Π⁴η); H2 − H1 = 0.0275 ± 0.0023, so the functionals clearly differ. At 32 000 also
  SPECIFIC: H1 0.0466 ± 0.0026, H2 0.0760 ± 0.0026. The instruments passed (mirror identity exactly 0, weights exactly,
  1+1 control E44b 0.0678 ± 0.0102 against 0.0584, no invalid seed). Three caveats: (1) the condition "the error of H1
  does not grow" passed only thanks to the statistical margin: pointwise the error grew from 0.0035 to 0.0093 (allowed
  up to 0.0109); (2) at 64 000 H1 lies only 0.0016 below the upper edge of the band; (3) for a large fraction of the
  seeds H1 = H2 or H1 = 0, because in 3+1 the chains are short and the heaviest chain is often the same for both
  weights, or straight. The geodesic law Π²η in 3+1 thus holds according to the rule given in advance, but the precision
  of the reading of H1 in 3+1 is on the edge and did not improve with N.
- **Π can be read back from the density of chains** (K.3, E45c, replicated); the geometry is readable from the order in
  3+1 from N = 32 000 (K.1, E54).
- **The causal cone as an instrument** (K.2, calibration, 2+1 and 3+1): delay 0.982 and 0.986 of the GR value, bending
  0.874 and 0.844 (in E47 0.913 and 0.919). The delay is exact; the bending reads systematically 10–15 % lower, because
  the local linear fit of the slope smooths the steep part of the 1/y profile. The instrument distinguishes γ = 1 from
  γ = 0.5, smaller differences not.
- **The calibration instruments of Causality Lab passed**: expected operator (K.4), unbiasedness of the operator (K.5),
  flux identities (K.6), mirror identity of chains (K.7), probe (K.8), null models of topology (K.9).
- **A causal set in the emergent geometry Π²η gives □_g − ½R in the continuum limit** (K.11 local, calibration): ratio
  of the operator to the local form of equation (a) 0.77 → 0.89 → 0.95 (2+1) and 0.87 → 0.954 → 0.9998 (3+1) for field
  widths 4, 8, 16 l_k; □_g without −½R would give 0.41 and 0.29. The proper volume of intervals (mean of Π^D from the
  exact projection of the diamond) agrees with Monte Carlo (K.11 diamond, max |z| 3.22 and 2.41; the coarsest level of
  the quadrature did not pass). Convergence by rows passed (K.11 conv). In 2+1 the curved operator is calibrated also
  for a distant field (K.11b far_d3: volume-weighted remainder Q_far/Q_src = +0.015, −0.0001, +0.002 at three levels).
- **The curved operator is calibrated in 3+1 also for a distant field, with a smooth (C²) geometry** (K.11c,
  calibration, 20 Oct). The criterion of K.11b unchanged (|Q_far(L3) − Q_far(L2)| <= 0.01 Q_src), only the geometry
  between nodes is a C² spline instead of a broken line. 3+1, grid step 1 l_k / ratio 1.08 (far_d4_g1), PASSED:
  Q_far/Q_src = −0.0107, −0.0132, −0.0101 (L1, L2, L3), |L3 − L2| = 0.0031; smooth reference at L2 −0.0106 (deviation of
  the spline 0.0026); L1 already suffices. The same field with the broken-line geometry gives +2.26, so the cause of the
  failure of K.11b is confirmed. 2+1 passed on both grids (far_d3_g1, far_d3_g2). 3+1 on the coarser grid step 2 l_k /
  ratio 1.12 (far_d4_g2): INSTRUMENT FAILED according to the criterion (|L3 − L2| = 0.0122 and the deviation of the
  spline from the reference at L2 0.0121, both above 0.01). Descriptively: L1 +0.0322, L3 +0.0301 and the reference at
  L2 +0.0300 agree; only the spline at L2 deviates (+0.0179), almost entirely in a single row r = 43 l_k (3.8e-5 against
  6.0e-5 in the others). The spline itself is therefore fine, the L2 level on the sparse grid is unreliable; even on the
  g1 grid L2 is farthest from the others. Rule: compute fields in 3+1 with the spline on the g1 grid (calibrated); in
  2+1 both grids suffice. This unblocks 3+1 for F1.4 and F1.5 (their 3+1 parts do not exist yet).
- **Superposition in equation (a) fails at second order, as (a) predicts** (F1.4, 2+1, mechanism): a core
  (σ = 8 l_k) and a shell (r0 = 30, width 4 l_k); the deficit Δ = u12 − u1 − u2 is positive everywhere (the common well
  is shallower than the sum: a source in the well of the other weakens), in the centre +0.046 against +0.057 of the
  continuum prediction of (a); RMS(Δ_cs − Δ_a) = 0.0101 at a tolerance of 0.0111, i.e. narrowly. Π(0): core 0.7099
  (continuum 0.7099), shell 0.915 (0.901), both 0.671 (0.668). Two caveats: (1) the quadrature of F1.4
  [160, 8, 48, 6, 24, 1] lies, according to the later K.11b, 1.3 % of the source outside its tolerance; (2) a shell of
  width 4 l_k is sensitive to non-locality (K.11 local at width 4: 0.77), its field came out about 14 % weaker than the
  continuum one and the deficit 18 % smaller; the comparison is with the local continuum equation, the criterion does
  not change. A replication with the quadrature according to K.11b would remove caveat (1).
- **PRELIMINARY: the field in equation (a) has the coupling κ/D and the non-linearity according to (a)** (F1.5, 2+1,
  mechanism). A run of an older version of the F1.5 file (quadrature [160, 8, 48, 6, 24, 1], node step 1 l_k,
  48 nodes, 94 rows), not the current version with the quadrature according to K.11b and step 2 l_k. P1: k = 0.3363
  against 1/D = 0.3333 (□_g without −½R would give 0.99). P2: strength of the distant field per unit source
  s(q) = 0.826, 0.709, 0.560, 0.402 against (a) 0.825, 0.709, 0.562, 0.408 for Π(0) = 0.89 to 0.58; RMS of the
  difference 0.0034 at a tolerance of 0.10. Π(0) fits the continuum (a) to 0.003. Preliminary, because the quadrature
  lies, according to K.11b, 1.3 % of the source outside its tolerance and a confirming run of the current version has
  not taken place yet; it cannot be presented as an independent confirmation. The distant field of the flat operator is
  3 % stronger than the continuum one (A_flat 65.9 against 63.9; F1.3 at σ = 1 l_k: 5.5 %).
- **Inertial and passive mass are the same** when the body factorises (F5.0, E44c, 1+1, conditionally).
- **The overlap of pasts is a signature of geometry**: 1.42–1.66 against 1.00 (F0.1, E43).
- **Shielding by saturation** (F3.1, E55a): a source at capacity limits its own outflow, the external flux drops
  exactly in the ratio C/F0_tube (k/k_pred = 1.0054), nothing saturates outside. Without feedback the region of
  resistance does not spread.
- **Without feedback no structure arises** (F3.2, E55b): Δ_pers = +0.006 ± 0.006, E_loc = +0.004 ± 0.004, equivalent
  to zero within a limit of 0.05; the ground state for F3.3.
- Capacity limits the flow, not the number of links (F1.8, E40, E41-A), but the verdict is not robust to the seed
  (K.10, E49).

## Does not hold

- **F5.1g: a flow whose memory is the common flow-past of its parents does not hold together and depends on the
  region; shared equally among its admissible continuations it spreads at ~0.65 c and fills the box, like radiation,
  not like matter** (mechanism, flat 3+1, density 1, box t [0, 24], |x| <= 9, |y|, |z| <= 7, 8 sprinklings of
  ~84 700 events, 2 sources per velocity (v 0 and 0.3), launch tube r_s 1.5 for 5 units, then free for 10; branches
  G0 (memory = parents), G1 (carried direction of F3.6), G2 (the hypothesis, regions B and A), G3 (borrowed memory);
  FAILED: M5 on 10 Nov by the criteria frozen in the card; run by the author, 55 231 s with 8 processes).
  - **Instrument (G1) OK:** Delta −0.016 ± 0.008 (v 0) and −0.010 ± 0.011 (v 0.3), width at the end 1.84 / 1.74,
    growth 0.080 / 0.052; the carried direction keeps its velocity and holds together (as F3.6).
  - **M5 fails:** M5a: of ~24 000 flow events only 579 are untainted in both regions (0.00 of the flow, the test needs
    0.10), and 129 of them differ; M5b: the width in region A is 29 % (v 0) and 23 % (v 0.3) smaller than in B
    (± 0.024), and the centre differs by 0.180 ± 0.046 at v 0.3.
  - **Descriptive (decides nothing):** G2 spreads ballistically: the width across grows from 1.85 to 4.5 in the first
    4 units (~0.65 per unit) and then stops at ~4.8–5.2 against the walls (half width 7); G0 and G3 the same (5.8 and
    4.8). The flow involves ~24 000 events (G0, G2) against 240 for G1. The centre at v 0.3: G1 moves 1.54 → 4.41 in
    10 units (0.29); G2 moves only 1.22 → 2.00 in the first 4 units (~0.2), then jumps to 6.9 once the back of the flow
    reaches the near wall, so Delta(G2, 0.3) = +0.350 ± 0.032 is a product of the box. Memory diagnostics of G2: |M|
    1.3–1.4, age 3.5–3.6 links, 2.6–2.7 flow parents, eta_common 0.57–0.58, 13.3–13.4 admissible children, D_min > 0
    in 12 % of the events (G0: 47 %). Growth G0 − G2 +0.073 ± 0.035 (v 0) and +0.042 ± 0.030 (v 0.3): the common memory
    slows the spread only slightly; G3 − G2 −0.017 ± 0.051 and +0.014 ± 0.049: whether the memory is the event's own
    or borrowed makes no difference.
  - **Design errors found after the run:** the growth of the width was read over the second half of the walk, when the
    flow had already reached the walls, so the frozen growth criterion measured the box (0.08 instead of ~0.65); the
    sources at v 0.3 started 6.75 from the near wall and 11.25 from the far one, so the walls cut the flow
    asymmetrically. Neither changes the verdict (M5 is checked first and fails in both of its parts).
  - **Reading:** equal sharing among the admissible continuations branches the flow at every event (13 admissible
    children), so the flow spreads at a large fraction of c whatever its memory; the zero-deficit condition and the
    common flow-past do not confine it. A free flow of this kind behaves like radiation. Coherence was shown only by a
    flow that does not branch (G1, the carried direction, which uses coordinates, M3). In the theory matter is a
    locked defect, i.e. a bound structure; the inertia and the fall of matter should be tested on a bound object, not
    on a free flow (a new hypothesis, R132).

- **F3.13d: no defect arises by itself: the alignment of disordered orientations leaves no locked defect at the scale
  R = 2** (mechanism, 3+1, box t [0, 30], half 7, ~82 200 events, 1 043 sources, 16 sprinklings, FAILED 1 Nov by the
  criteria approved before the run; deciding kappa = 0.9, changed from 0.3 with the author's approval before the full
  run because kappa 0.3 was unreadable; instrument OK: vacuum 0 on 21 566 of 21 566 readable shells, random_0.9 5 916
  readable shells).
  - **Sites:** S = 0.000 ± 0.000 (fail < 0.1); not one shell of 5 916 readable at kappa 0.9 (1 109 at kappa 0.3) read a
    non-zero degree; persistence undefined (no site); net charge on R = 4 read twice, 0 both times.
  - **Smoothing:** mean cost per carrying event 0.215 -> 0.042 (kappa 0.9), 0.488 -> 0.091 (0.3), 0.667 -> 0.311 (no
    alignment); readability rises to 0.37-0.40 at kappa 0.9 from t ~ 17; the vacuum loses readability late (0.68 at
    t = 27) as merging thins the flow.

  Reading. With an orientation carried by the flow and its alignment (the two assumptions added in F3.13a/b), disorder
  at the sources does not leave behind a locked defect that can be read at the scale R = 2 inside a region of radius 5:
  whatever knots exist at the scale of single events are untied, or leave through the boundary, before the field is
  resolved. The hedgehog of F3.13b2 is therefore a configuration that is stable once made but is not made by this
  dynamics: in this branch the defect remains put in by hand. Limits: one size of the region (sources r < 5), one scale of
  reading (R = 2), a cold flow at rest without a field, and a field resolved only from t ~ 13 on; a larger region could
  keep defects that here escape. The branch of a carried orientation does not, as tested, create particles by itself.

- **F3.12b search: capacity as an elastic wall does not hold a bundle in its own field; without merging a warm bundle
  (η = 0.06) comes closest to a stationary width of all F3 runs, but even there the classification gives no bound arm**
  (mechanism, exploration, 2+1, F3.11 set-up, 32 sprinklings of ~49 000 events, 201 source events, RECORDED 31 Oct;
  candidate by the rule in the header: none; instrument OK: frame check identical in all 17 arms on 8 sprinklings, rest
  mass of the stream arms conserved to 3e−16, balances ≤ 2e−15, field_0 contracts against free_0).
  Tolerance of boundness s(free_0)/4 = 0.0103 (W1 [12, 28]) and 0.0063 (W2 [28, 44]). Slopes of the width W1 / W2:
  - free_0 +0.041 / +0.025; field_0 −0.010 / −0.038 (collapse, R 3.0 → 2.5); free_cl_0.06 +0.055 / +0.056;
  - streams without merging and without capacity: cl_0.03 −0.002 (bound) / −0.017 (collapse); **cl_0.06 +0.011 ±
    0.003 / −0.004 ± 0.003 (both unresolved): R holds 3.39–3.50 from t ≈ 15 to 45**; cl_0.1 +0.033 / +0.032 (grows);
  - capacity C = 2 / 4 / 8 (reflections per stream 0.39–0.41 / 0.35–0.37 / 0.20–0.31): at η = 0 collapses like
    field_0, slower (W2 −0.020 to −0.023); at η = 0.03 bound or unresolved in W1, collapse in W2; at η = 0.06 grows in
    W1 (+0.020 to +0.036) and unresolved in W2 (−0.006 to +0.001), plateau R ≈ 3.5–4.0; at η = 0.1 grows.
  - The internal spread grows in every arm with a field (field_0 from 0 to 0.14, cl_0.06 from 0.06 to 0.14): the
    discrete bending heats the bundle, as in F3.7.

  Reading (declared before the run): reflection off overloaded events does not balance the field in this form. The
  capacity pushes the bundle out early (a larger plateau) but neither stops the collapse of a cold bundle nor holds a
  warm one better than streams without capacity. The new fact is the role of merging: without it, internal motion alone
  nearly balances the field at η ≈ 0.06 (F3.11 with merging only gave plateaus followed by contraction). What prevents
  a stationary state now is not cooling but heating: the spread rises steadily, so a plateau is a slow drift between
  two regimes, not an equilibrium. Limits: 2+1 (logarithmic potential), conformal bending (γ = −1), rest-mass source.

- **F3.12a2: the topological cost of rerouting is not the primary resistance of the substrate; capacity is**
  (calibration, revision of F3.12a, 2+1, box t [0, 10], |x| ≤ 10, ~4 000 events, ~85 000 links, 32 new sprinklings,
  lens 3 ≤ t ≤ 7, w = 2.5, RECORDED 31 Oct with the classification "capacity dominates"; instrument OK: unperturbed
  I_flow 2.4e−8, balance 3.3e−10, all programs solved). The author's hypothesis (31 Oct): the impedance is the
  topological cost of rerouting the flow, I_v = C(H → H′ | Φ, Q, K0); his definitions: I_flow = the smallest L1
  change of the flow that the new network admits (linear program), I_topo = created and removed links separately.
  Split I_flow = I_forced + I_cascade; I_cascade at σ0 = 0.1 / 0.2 / 0.3:
  - T (narrowed cones, fewer paths, 3528 / 5513 / 6848 links removed): 0.0045 / **0.0134 ± 0.0003** / 0.0196;
  - K (capacity cut to the same forced change, the same paths): 0.0383 / **0.0525 ± 0.0011** / 0.0598;
  - P (throughput: narrowed cones and added events, 659 / 1170 / 1641 links created): 0.0042 / 0.0125 / 0.0170.
  Criterion at σ0 = 0.2: K − 2 SE = 0.050 > 2 (T + 2 SE) = 0.028, capacity dominates; the same at all strengths and
  for P against K. F3.12a (instrument failed on its tolerance) gave the same numbers to the second decimal.

  Reading (declared before the run): resistance is the capacity effect; the topological candidate is not supported as
  the primary one. When paths vanish the network reroutes almost freely: beyond the flow that has to leave the lost
  links, only ~1 % of the link flux changes (it grows with σ0, so topology has a small cost of its own). A cut in
  capacity at the same forced change costs 3–9× more. Adding events (the throughput principle) adds no resistance of
  its own. Limits: one instrument (L1 distance to the uniform flow), one geometry (a lens), 2+1. For a bound structure
  (F3.12b) the candidate is therefore capacity: a flow that meets a full capacity must be deflected outwards with
  exact four-momentum, i.e. elastically, not merged (F3.9 had capacity without an outward momentum and gave only a
  plateau).

- **F5.1e: in 3+1 the same history from the apex does not fall: it takes an inward velocity at the apex and then keeps
  it, without acceleration; the longest chains between fixed events read the geometry well** (mechanism, 3+1, order of
  F7.1c σ0 = 0.3, w = 4, box t [0, 24], |x| ≤ 8, 16 + 16 sprinklings of ~105 000 events, 60 chains, rule of F5.1c/F5.1d
  m = H = 8, FAILED 30 Oct).
  - **Instruments OK:** flat R = +0.139 ± 0.138, longest chains R_max = −0.011 ± 0.071; loaded longest chains read the
    geometry, R_max = 0.850 ± 0.058 (better than in 2+1).
  - **Loaded order:** R = a / a_pred = **−0.134 ± 0.121** (R + 2 SE = 0.11 < 1/2: does not fall), a_pred = −0.02187
    (φ = −σ/2, fall from rest at r_apex = 4.15). 59.5 of 60 chains reach the end. Displacement from the apex,
    τ = 0, 1, …, 12: 0.00 −0.07 −0.18 −0.28 −0.39 −0.48 −0.57 −0.65 −0.72 −0.80 −0.88 −0.94 −1.01 (the prediction
    −1.62).
  - **Descriptively (not evaluated):** the displacement is almost linear. The history moves inwards at ~0.1 from the
    apex on and slows slightly (−0.10 at the start, −0.06–0.08 at the end), hence R slightly negative. At the apex the
    velocity jumps from the rising half (+0.067) to v0 = −0.102 ± 0.021, a jump of −0.17; in the flat order the jump is
    −0.045 (+0.070 → +0.025). The order is felt once, at the start, as a kick, not as a continuous acceleration. The
    same kick appeared in 2+1 (F5.1d, −0.09), where an acceleration followed it.
  - **Reading (declared before the run):** the rule that falls in 2+1 loses the gravity of the order in 3+1; a property
    of the rule at m = H = 8, not of the order, whose longest chains between fixed events read the geometry (R_max 0.85).
    A possible cause (a hypothesis, not tested): in 3+1 the effect of the field on the lengths L(p, q) over the horizon
    is below one link, so the choice of q* is decided by the tie-breakers rather than by the geometry. A new rule or a
    rescaled one would be a new hypothesis with its own criteria (R132). F5.1d stays PASSED for 2+1.
  - **Diagnosis 7 Nov (decides nothing):** the tie-breaker hypothesis above does not hold; see "Fall diagnostic"
    under Undecided (descriptive).

- **F5.1c: a flow that continues its history as the longest chain moves inwards in an order curved by matter, but does
  not fall with the Newtonian acceleration** (mechanism, 2+1, order of F7.1b σ0 = 0.45, 32 + 32 sprinklings of
  ~47 700 events, 40 chains, memory m = 8, horizon H = 8 links, FAILED 29 Oct).
  - **Instruments OK:** flat order R = +0.148 ± 0.081, shift +0.17 over 22, R_max = −0.31 ± 0.09; loaded
    R_max = 0.863 ± 0.064 (this confirms that the narrow failure of the control in F5.1 carried, 0.630, was scatter).
  - **Loaded order:** R = a / a_pred = 0.337 ± 0.049 (R + 2 SE = 0.435 < 1/2, does not fall), v0 = −0.158 ± 0.010; 39
    of 40 reach the end. Inward shift −4.0 over 22, more than the geometry (−2.7), but with a different shape:
    t = 0, 2, …, 22: −0.01 −0.17 −0.46 −0.82 −1.25 −1.72 −2.22 −2.74 −3.28 −3.68 −3.86 −4.01.
  - **Descriptively (post-hoc analysis, decides nothing):** three phases. Up to t ≈ 8 a strong inward acceleration,
    about 2× a_pred; this is the time when a virtual history at rest serves as p (the same "early kick" as in the
    development of F5.1b). For t = 8–16 the acceleration is about 0.5× a_pred. From t ≈ 16 a slowdown; there the
    horizon H = 8 links exceeds the upper edge of the box T = 32 and the rule passes to a smaller depth. A fit over the
    whole 22 mixes these phases.
  - **Reading:** the rule carries inertia (the velocity holds) and feels the gravity of the order, but a clean Newtonian
    acceleration does not come out of the whole course. The initial condition (a virtual history at rest is not a
    geodesic) and the edge of the box are two recognised sources of distortion; by R132 this does not change the
    verdict. The revision F5.1d follows.

- **F5.1 intrinsic: a chain that carries its origin and keeps to the layer of the longest chains from it does not fall
  in an order curved by matter, it moves away instead; the longest chains between fixed points in the same order read
  gravity correctly** (mechanism, 2+1, order of F7.1b σ0 = 0.45, 32 + 32 sprinklings of ~47 700 events, 40 chains,
  FAILED 29 Oct).
  - **Instruments OK:** flat order without drift (R = 0.150 ± 0.069) and without bending of the longest chains
    (R_max = −0.04 ± 0.13); in the loaded order the longest chains bend outwards, R_max = 0.833 ± 0.076 (prediction of
    the geometry φ = −σ/3: 1).
  - **Loaded order:** R = a / a_pred = −0.299 ± 0.068, outward shift +0.95 over 22 (a fall would be −2.7). The chain
    reaches the end only in 24.7 of 40 (events without a kept descendant).
  - **Probable cause (descriptively):** the layer L(origin, k) = min selects descendants to which the shortest longest
    chain leads, and those are the sparser regions farther from the source. The choice "stay on the longest chain from
    the origin" step by step thus pulls outwards, not along the geodesic.

  Reading: the order carries gravity (the longest chains between fixed points read it), but the local rule "memory of
  the origin" does not pass it on to matter. Free fall needs a memory of velocity, i.e. a discrete parallel transport
  of direction, and that has not yet been built from order and counts. Three attempts at a rule from order alone failed
  (two in development: speed 0.4c and jumps along long links, the third with this run). An open problem of coupling
  matter to the gravity of the order.

- **F6.8c: the optical sector ν² FAILED according to the registered reading; the faulty part is the reading model, not
  necessarily the action** (chain, 3+1, 29 Oct, FAILED; the verdict and the criteria stand). Instruments I1–I4 OK: ⟨T2⟩
  against a numerical Fermat ray ≤ 7.9e−5, V'' against an independent computation 6.6e−4, uniform null 3.6e−15,
  conformal K̂ against F6.8a 1.6e−10, null k = 0 ≤ 2.1e−4 in units of K̂_bb, level 2 ≤ 3.2e−4.
  - The registered reading (linear fit in T per k, extrapolation T → 0) gave K0_bb = 1.559, c1 = 98.8, c0 = −156,
    rms 4.2e−2. The fit carries no information: K̂_bb falls with T faster than linearly.
  - K̂_bb at k = 0.05 / 0.07 / 0.1 / 0.14 / 0.2 and kT = 8 / 12 / 16: 0.9527 / 0.9035 / 0.6886; 0.9167 / 0.9029 /
    0.8646; 0.8485 / 0.8327 / 0.8084; 0.7373 / 0.7076 / 0.6659; 0.5372 / 0.4769 / 0.3927.
  - **Descriptively (not evaluated):**
    - At k = 0.14 and 0.2, K̂_bb ≈ K_loc − q k⁴T² with q ≈ 0.019 and no linear term. At k = 0.1 it already deviates
      slightly, at k = 0.07 and 0.05 (T ≥ 171) the scaling breaks down. Suspected numerical loss in the bending terms at
      T ≳ 200 (the null k = 0 does not see the bending), undecided. (This form was later withdrawn as overstated, see
      F6.8c2.)
    - The reading at fixed n = kT (K̂ = K0 − k²[c1 ln(1/k) + c0(n)]) does not assume a form of the age dependence.
      At n = 8: conformal K0 = 0.99993 (F6.8a: 1.00005, the method fits), mixed 0.99988, **optical 0.99914**
      (k ≤ 0.14, residual 1.4e−4; k = 0.07–0.2: 0.99966). Hence C(0) ≈ 1.998 and γ(0) = 2K0_ab/K0_bb − 1 ≈ 1.0015.
      The data had already been seen, so this is decided in the revision F6.8c2 on new k.
    - If a term ~(kT)² in the ν² sector is physical and persists asymptotically, it dominates the long-term behaviour
      at a cosmological T. A pure k⁴ term gives only a contact term after Fourier inversion, but the decay scale
      (qℓ_k²T²)^{1/4} would require ℓ_k ≲ 10⁻³⁴ m for 52 μm. To be verified independently (an asymptotic integral as for
      c_T).

  Reading: F6.8c did not decide the question C(0) = 2. The linear model of the age, taken over from the conformal and
  mixed sectors, does not hold for the optical sector. The prediction F8.1 (linear age term) is therefore not launched
  yet.

- **F7.1a: the candidate from capacity (a loaded event keeps only slower links) gives γ ≈ 1/2, not 1, independently of
  the strength of the coupling, as the derivation before the run predicted** (mechanism, 2+1, 16 sprinklings of
  ~18 000 events, FAILED; instrument K.16b, flat control OK: without load it reads γ −0.998 / −0.999 and all readings
  zero). Load σ = σ0·exp(−r²/128), a link kept when its speed <= exp(−σ) at the parent event:
  - σ0 = 0.2 (42 % of links kept): γ_A (cones + density) +0.465 ± 0.030, γ_B (clocks + cones) +0.433 ± 0.033;
    difference of the routes 0.03 ± 0.05, the order is metric; from the prediction 1/2 −0.04 and −0.07, from GR −0.54
    and −0.57;
  - σ0 = 0.3 (34 %): +0.490 ± 0.013 and +0.531 ± 0.029;
  - the read speed of light fits the rule to thousandths (ln c −0.198 against −0.198 in the centre): chains through the
    kept links reach up to the narrowed cone; the density stays flat (|d| <= 0.03).

  Reading: removing links is not enough. A rule that only narrows the cones and keeps the events gives ψ = φ/2 (in 3+1
  φ/3), however strong the narrowing, and the geometry from it is a consistent metric (clocks and cones give the same
  γ). For γ = 1 the substrate must at the same time add events near matter: a density per coordinate volume e^(−φ) in
  2+1 (e^(−2φ) in 3+1) together with a narrowing of the cones to e^(2φ). Capacity (axiom 4) as the only coupling
  Π → relation is thereby excluded for γ.

- **F3.11 search: a hot bundle in its own field gives no bound state; at η = 0.15 it holds a plateau for ~30 units of
  time and then slowly contracts, like a cooling self-gravitating gas** (mechanism, exploration, 2+1, 32 sprinklings,
  ~49 100 events, 202 source events, RECORDED; candidate according to the rule in the header: none; instrument OK:
  frame check in all branches, balance of the free branches ≤ 5·10⁻¹⁵, field_0 contracts relative to free_0).
  The field and bending of F3.8 (ε = 0.1 locked), merging B and the hot source sym from F3.6c. Tolerance of boundness
  s(free_0)/4 = 0.0085 (W1) and 0.0103 (W2). Slope of the width in W1 [12, 28] / W2 [28, 44]:
  - free_0 +0.034 / +0.041 (grows);
  - field_0 −0.016 ± 0.006 (unresolved) / −0.037 ± 0.007 (collapse); agrees with F3.8 ε = 0.1 (−0.004 / −0.036);
  - **field_0.15 +0.009 ± 0.007 / −0.019 ± 0.006 (both unresolved)**: R holds 3.4–3.6 from t ≈ 13 to 37, then falls
    to 3.3–3.4; without the field (free_0.15) it grows from 3.5 to 5.3;
  - field_0.3 +0.062 / +0.034 and field_0.45 +0.146 / +0.090: they grow, the field only brakes them (against free
    0.111 / 0.112 and 0.179 / 0.135).

  The internal spread of rapidities in the fields grows: field_0 from 0 to 0.136, field_0.15 from 0.05 to 0.13 (without
  the field it falls to 0.035). Contraction therefore heats the bundle (the potential energy of the field passes into
  internal motion) and merging cools it. The plateau at η = 0.15 is a quasi-equilibrium that cooling slowly shifts
  towards contraction. This corresponds to F3.9 (capacity: a plateau, then slow contraction) and to the contraction of
  a star without an energy source.

  Reading: a resistance of the substrate with range (a field) and internal motion can balance each other, but only
  temporarily. Merging with F = |P| does not erase energy, it converts internal motion into rest mass, which then no
  longer pushes. A stationary bound state would need a non-dissipative pressure or an energy source, i.e. a new
  mechanism. A finer grid of η would not remove the cooling. Caveats: 2+1 (logarithmic potential), the bending of F3
  corresponds to γ = −1, the field is driven only by rest mass.

- **F3.6c (eta03 and eta06): merging with exact 4-momentum (F = |P|) does not hold a bundle from a hot source; it only
  slows the spreading roughly by half** (mechanism, 2+1, 32 sprinklings, ~61 000 events, 202 source events, both parts
  FAILED; instrument fine in all conditions). Part eta03 (η = 0.3). Expansion speed in the window [6, 18] (centroid
  frame, R_rms around the centroid):
  - sym: B 0.112 ± 0.005, A 0.114 ± 0.005, free motion C 0.220 ± 0.002 (boundness threshold C/4 = 0.055);
  - asym: B 0.117 ± 0.007, A 0.118 ± 0.007, C 0.214 ± 0.004;
  - cold (control): A = B = C bit by bit, 0.049 ± 0.003 (diffusive jitter of the position, the same as F3.6).

  Balance of B 4·10⁻¹⁵ and 1·10⁻¹⁵, B − A −0.0011 ± 0.0007 and −0.0010 ± 0.0016 (no measurable difference). Merging
  (~136 merges, mean equivalent rapidity 0.24) cools the streams: the internal spread falls from 0.30 to 0.13 at t = 5
  and further slowly to 0.088 at t = 29. The number of streams falls from 78 to 35. The increase of rest mass in B is
  3.9 %, the old rule A loses 3.7 % of energy. The width of B grows linearly (3.05 → 5.56 over 24), i.e. ballistically,
  not diffusively.

  Reading according to the author's rules: the rule F = |P| works as bookkeeping (4-momentum exactly), does not address
  the width and in this regime does not differ measurably from the old rule. Merging is only cooling: streams that once
  separate never meet again and nothing returns them. For the Jacobson route F = |P| remains necessary (A loses
  energy); it changes nothing about localisation.

  Part eta06 (η = 0.6, FAILED):
  - sym: B 0.237 ± 0.007, A 0.249 ± 0.007, C 0.465 ± 0.003 (threshold C/4 = 0.116);
  - asym: B 0.243 ± 0.007, A 0.251 ± 0.006, C 0.440 ± 0.003.

  Balance of B 4·10⁻¹⁵ and 1·10⁻¹⁵. Here B differs measurably from A: B − A = −0.0115 ± 0.0023 and −0.0079 ± 0.0021,
  i.e. B is 3–5 % narrower. A loses 12 % of energy, B converts it into rest mass (+14 %). The internal spread falls
  from 0.60 to 0.27 at t = 5 and then almost stops (0.22 at t = 29). The number of streams settles around 50: merging
  is quickly exhausted.

  Summary of both parts: the ratio v_s(B)/v_s(C) is 0.51–0.55 at both η. Merging cools only the streams that meet right
  at the source; afterwards the diverging streams never meet. F = |P| is the right bookkeeping and at higher η narrows
  the bundle slightly (more weight to heavier, already merged streams), but it gives no binding. What holds the width
  stays open.

- **F3.9 search: a field with range together with a finite channel capacity gives no bound state; capacity only damps
  the rebound and holds a plateau briefly, then the bundle slowly collapses** (mechanism, exploration, 32 sprinklings
  2+1, RECORDED; candidate according to the rule in the header: none). F_ref = 7.0, φ_ref = 43.5. Slope of the width in
  W1 / W2:
  - free: +0.052 / +0.031;
  - replication of F3.8: ε 0.1 −0.010 / −0.031 (collapse), ε 0.5 −0.030 / +0.235 (collapse and rebound) — agrees with
    F3.8;
  - **cap ε 0.1 κ 1: +0.0018 ± 0.0030 (bound in W1) / −0.027 ± 0.004 (collapse in W2)**: the width holds ~3.4 from
    t ≈ 13 to 27, then falls to 2.95;
  - cap ε 0.1 κ ¼: +0.018 / −0.019 (rises to 3.9 and returns);
  - cap ε 0.5: collapse and rebound, capacity damps the rebound (W2 +0.205 / +0.102 / +0.030 for κ 4 / 1 / ¼ against
    +0.235 without capacity);
  - cap ε 0.01: like free.

  The frame check passed in all branches. Reading (descriptive): capacity in this form sends the excess along less
  parallel paths, but with the same carried direction. It spreads the positions but does not carry momentum outwards,
  so it is not pressure in the proper sense. The merging of carriers moreover averages directions, i.e. loses relative
  motion (inelastic, "cools"). The slow contraction after the plateau corresponds to this, similar to the contraction
  of cooling self-gravitating matter.
- **F3.7 search: bending towards slower clocks with clocks from the local load does not bind the bundle but heats it;
  the covariant narrowing of the cones does not bind it either** (mechanism, exploration, 32 sprinklings 2+1,
  ~49 100 events, RECORDED; candidate for confirm according to the rule in the header: none). Clocks
  Π = 1 − 0.5 min(F/C, 1)², C = κ F_ref, F_ref = 6.0. Slope of the width in the windows W1 [12, 28] and W2 [28, 44]:
  - free (F3.6): 0.037 ± 0.007 and 0.038 ± 0.007;
  - lens κ = 4 / 1 / ¼: 0.12 / 0.47 / 0.59 in W1, 0.43 / 0.43 / 0.33 in W2. All grow faster than the free bundle. At
    κ = ¼ and 1 the bundle first contracts (R 3.0 → 2.5 until t ≈ 9), then explodes all the way to the wall. The
    bending thus points inwards, but without damping the paths fly through the centre and the noise of the gradient
    estimate from the parents heats the bundle (mean bending 0.05–0.30 rapidity per step, Π falls down to 0.5);
  - narrow κ = 4 / 1: spreads at 0.17 / 0.37 up to the wall. At κ = ¼ the cone in the core is so narrow that the
    fallback choice applies and the core behaves like free (0.044 ± 0.007), while the edge spills (~2 000 carriers).

  The frame check (β′ = 0.6) passed in all branches. Reading: with clocks directly from the local load the field is
  rough and short-ranged (only inside the bundle), its gradient is mainly the noise of the discrete load. Binding would
  need a smooth field with range, i.e. a field driven by the flow through a field equation (A5′, the retarded Green
  function of the causal set), or damping (exchange of momentum). Both are a new hypothesis after a failure, with their
  own criteria.
- **F0.6 fixed: narrowed cones drawn in coordinates (F3.3h) select a preferred frame** (mechanism, relativity gate,
  128 sprinklings N = 72 000, FAILED; all 6 deciding comparisons EXCEEDED). The same events, the same order, driving and
  readout in the rest frame of the source; only the rule draws the cones in the frame β:
  - s_bg = 0: in the rest frame 86 % of the sprinklings survive (N late 8.6), at β = 0.3 and 0.6 all vanish. The frame
    effect −8.56 ± 0.48 is as large as the whole effect of the mechanism (+8.56; tolerance 2.14).
  - s_bg = 0.2: N late 116 → 78.9 → 18.7 for β = 0 / 0.3 / 0.6 (frame effect −37 ± 1 and −97 ± 1 against a tolerance
    of 29). The formation drifts along the axis of motion: centroid +0.224 ± 0.002 and +0.339 ± 0.002 against the
    tolerance r0/4 = 0.0375; R_late grows by 0.08 and 0.19.
  - The change of frame changes the number of allowed descendants for 93–96 % of the events in the saturated region
    (5 % at s_bg = 0).

  The positive control (Π_sat = 1, cone = light cone) was bitwise identical in both frames and the balance fits: the
  difference comes from the rule, not from the transformation. Consequence: results built on this rule (F3.3h, the
  compact domain of F6.3b A, whose local relation draws the cones the same way) are not physics. The idea of narrowing
  the cones does not fall with it, only this coordinate form does. A covariant version would narrow the cone relative
  to the carried direction of the flow (the rest frame of matter, F3.6).
- **Feedback through capacity (C_v = C·Π_v) does not create structure** (F3.3 capacity, mechanism, 2+1,
  128 sprinklings, Π_sat = 0.5, m = 2, FAILED). The feedback saturates 54 % of the events (without it 0.9 %), but
  uniformly: paired differences against disabled feedback on the same graphs d_pers = −0.007 ± 0.008,
  d_loc = −0.010 ± 0.004 (both EQUIVALENT to zero, limit 0.05); with it E_loc = −0.004 ± 0.001 (more uniform than
  chance), Δ_pers = +0.001 ± 0.002. After the source is switched off the saturation vanishes within one layer
  (0.54 → 0.006 → 0). The branch without feedback repeated F3.2 with new seeds (Δ_pers +0.008 ± 0.008,
  E_loc +0.006 ± 0.004). Interpretation (descriptive, not a test): in one causal pass the state of an event is given
  only by its own load, so channel A only shifts the saturation threshold (for m = 2 to F/C ≈ 0.54) and the saturation
  copies the envelope of the flow; a saturated event moreover sends on less, which lowers the load of its successors.
  An event does not pass its state to its descendants (the store ΔQ is not transferred), so a persistent structure has
  no carrier.
- **Feedback through the source (s_v·Π_v) removes the saturation** (F3.3 source, mechanism, 2+1, 128 sprinklings,
  Π_sat = 0.5, m = 2, FAILED according to the rule "saturation vanishes"): saturated fraction 0.9 % without feedback,
  with feedback 1 × 10⁻⁶ (saturated events only in 2 of 128 sprinklings). Negative feedback: near capacity the clocks
  slow down and the source drops to as little as half, so the load does not exceed the capacity. Without saturated
  events there is nothing to localise.
- **Routing towards low impedance (weights Π_k/I_k, I = 1) removes the saturation** (F3.3 routing, mechanism, 2+1,
  128 sprinklings, Π_sat = 0.5, m = 2, stationary state, FAILED according to the rule "saturation vanishes").
  Saturated fraction 0.93 % without feedback, with feedback 0 in all 128 sprinklings; smallest Π 0.546 > Π_sat. The
  stationary state converged (21 passes), the instruments passed. Interpretation (descriptive): an event sends flow
  preferentially to descendants with faster clocks, i.e. less loaded ones, so the load evens out below capacity
  (negative feedback, as the header of the test assumed). All three channels of F3.3 therefore FAILED: A floods, B and C
  remove the saturation; there is no winner for step D (combination).
- **Memory through a transferred store does not create structure (at C = Q_0.9)** (F3.3b q90, mechanism, 2+1,
  128 sprinklings, FAILED). The store Q is transferred along links to the descendants uniformly and offered again;
  nothing else changes. The memory is measurable: saturated fraction 8.0 % → 9.9 %, mean held store +32 %. The
  structure is not: paired d_pers = +0.0006 ± 0.0016, d_loc = −0.0015 ± 0.0014 (both EQUIVALENT to zero with a large
  margin, limit 0.05); absolutely Δ_pers +0.001 ± 0.002, E_loc −0.008 ± 0.002 (saturated events are slightly more
  uniform than chance, −0.007 even without memory). After the source is switched off saturation and store vanish within
  one layer of 0.25 in both branches (10.5 % → 0.1 %). The instruments passed, including bitwise agreement when
  evaluated in the order of the boosted frame. Interpretation (descriptive): uniform splitting among descendants
  spreads the store into the future cone just like the flow, so the memory dilutes and does not create a carrier that
  would stay concentrated along a timelike direction.
- **A copied pulse (a rule from the main parent) floods the region but creates no localised structure** (F3.3c
  main_q99, mechanism, 2+1, 128 sprinklings, FAILED). The pulse ω = 1 does not split, each descendant receives the
  pulse of the parent from which most flow arrives, it lives only in saturated events and acts as a source. Saturation
  0.85 % → 30.5 %, but uniformly: paired d_pers = 0.000 ± 0.008, d_loc = −0.004 ± 0.004 (both EQUIVALENT to zero);
  with the pulse Δ_pers +0.005 ± 0.002, E_loc +0.003 ± 0.002. New against F3.3–F3.3b: the structure survives switching
  off the source. Saturated fraction after switching off 44 % → 8.1 % → 1.6 % → 0.05 % → 0 in layers of 0.25 (without
  the pulse it vanishes at once), and the pattern after switching off resembles the pattern before it (r_s +0.66,
  +0.41, +0.12). The pulse is therefore a memory that does not dilute and carries the position too, but under driving it
  forms an extended flooding, not a localised node. The instruments passed, including bitwise agreement in the order of
  the boosted frame.
- **A materialisation threshold ΔQ_crit = 0.25 gives a weak but real localisation, below the limit 0.05** (F3.3d
  crit025, mechanism, 2+1, 128 sprinklings, FAILED). The model F3.3c main_q99 with a single change: the pulse lives only
  at ΔQ ≥ 0.25, otherwise it vanishes immediately. The flooding disappears: saturation 0.92 % → 4.4 % (F3.3c 30.5 %),
  the pulse carries 2.3 % of the events. Paired d_pers = +0.010 ± 0.008 and d_loc = +0.019 ± 0.004 (both EQUIVALENT,
  limit 0.05). The localisation is at the same time 4.7 SE above zero: the first non-zero d_loc in the F3 series (with
  the pulse E_loc +0.020 ± 0.003, F3.2 +0.004 ± 0.004), roughly 2.5× below the limit. The price is memory: after the
  source is switched off saturation 6.5 % → 0.27 % → 0.01 % → 0 and the pulse 3.6 % → 0.11 % → 0 (F3.3c 44 % → 8.1 % →
  1.6 %); r_s with the last segment before switching off +0.13 in the first layer (few cells, only reported).
  Interpretation (descriptive, not a test): the threshold trades persistence for concentration; the pulse survives only
  in the most heavily loaded places, which are clustered, but without driving the load does not stay above the
  threshold. The instruments passed, including bitwise agreement.
- **A stronger threshold ΔQ_crit = 0.5 does not increase the localisation but lowers it** (F3.3d crit050, mechanism,
  2+1, 128 sprinklings, FAILED). The pulse carries 0.16 % of the events, saturation 0.92 % → 1.2 %. Paired
  d_pers = +0.002 ± 0.004, d_loc = +0.007 ± 0.002 (both EQUIVALENT, limit 0.05; d_loc 3.0 SE above zero, about a third of
  the value at 0.25). The memory after switching off vanishes within one layer (saturation 1.3 % → 0.01 % → 0, pulse
  0.19 % → 0). Together with crit025: above 0.25 both localisation and memory weaken with the threshold, there is no
  window above 0.25 in the scan given in advance.
- **A weaker threshold ΔQ_crit = 0.1 brings back flooding and memory, the localisation drops** (F3.3d crit010,
  mechanism, 2+1, 128 sprinklings, FAILED). Saturation 0.92 % → 17.3 %, the pulse carries 14.8 % of the events. Paired
  d_pers = +0.007 ± 0.007, d_loc = +0.008 ± 0.004 (both EQUIVALENT, limit 0.05). After switching off saturation
  27.7 % → 3.6 % → 0.41 % → 0, pulse 24.3 % → 2.8 % → 0.28 %, pattern similar to the pattern before switching off
  (r_s +0.58, +0.25, +0.16): between F3.3c and crit025.
- **An extreme threshold ΔQ_crit = 1.0 practically switches the pulse off** (F3.3d crit100, FAILED): the pulse carries
  < 0.01 % of the events, saturation 0.92 % → 0.93 %, d_pers = +0.001 ± 0.001, d_loc = +0.000 ± 0.001 (EQUIVALENT);
  model = F3.2.
- **The scan of the materialisation threshold given in advance FAILED over the whole range** (F3.3d, ΔQ_crit 0.1 /
  0.25 / 0.5 / 1.0, with F3.3c as ΔQ_crit = 0). Paired d_loc: −0.004 / +0.008 / +0.019 / +0.007 / +0.000; saturation
  30.5 / 17.3 / 4.4 / 1.2 / 0.9 %; saturation in the first layer after switching off 8.1 / 3.6 / 0.27 / 0.01 / 0 %. The
  localisation has a single maximum around 0.25 (agrees with the author's estimate 0.15–0.25), in size about 2.5× below
  the limit 0.05; memory falls monotonically with the threshold. A node that would concentrate and at the same time
  persist without driving does not arise within the range of the scan. The maximum is searched for and verified by
  F3.3e.
- **Not even the optimum of the threshold reaches the limit** (F3.3e search, exploration, RECORDED, 2+1, 64 new
  sprinklings, grid ΔQ_crit 0.05–0.40 in steps of 0.025 on the same causal sets). The quantity searched for,
  min(d − 2 SE), has a broad, shallow peak around 0.15–0.2 (agrees with the author's estimate): best ΔQ_crit = 0.2,
  d_pers = +0.034 ± 0.009, d_loc = +0.015 ± 0.006, min(d − 2 SE) = +0.002; passing would need > 0.05. d_pers is positive
  by 3–3.6 SE over the whole range 0.1–0.275 (in F3.3d at 0.1 and 0.25 zero within the error; difference ~1.3 σ). The
  memory after switching off falls monotonically with the threshold (saturation in the 1st layer 5.3 % at 0.05 → 0.62 %
  at 0.2 → 0.03 % at 0.4). The instruments passed. The record was not saved at the first run (outdated
  core\provenance.py on the author's machine); the numbers are from the printout, the run is repeated.
- **The materialisation threshold does not pass even at its optimum on independent data** (F3.3e confirm, mechanism,
  2+1, 128 new sprinklings, ΔQ_crit = 0.2 chosen from the data, FAILED). Paired d_pers = −0.009 ± 0.007,
  d_loc = +0.010 ± 0.004 (both EQUIVALENT, limit 0.05). Saturation 0.92 % → 6.1 %, the pulse carries 3.9 % of the
  events; after switching off saturation 10.4 % → 0.57 % → 0.01 % → 0. The positive d_pers of the search
  (+0.034 ± 0.009) did not repeat; d_loc stays weakly positive (2.2 SE). The threshold model (a copied pulse from the
  main parent, C = Q_0.99) on a fixed causal order is thereby closed over the whole range and at the optimum:
  localisation exists, but is ~5× below the limit, and persistence above the level without a pulse does not arise.
- **Threshold + feedback through capacity creates a self-sustaining but extended saturated phase, not a node** (F3.3f
  capacity, mechanism, 2+1, 128 sprinklings, ΔQ_crit = 0.2 from the data of F3.3e, C·Π, Π_sat = 0.5, m = 2, FAILED).
  Saturation: without the pulse 1.0 %, threshold only 7.6 %, threshold + channel 99.1 %; all saturated events have
  Π = Π_sat = 0.5. Paired against the branch without the pulse d_pers = +0.004 ± 0.009, d_loc = +0.036 ± 0.005 (both
  EQUIVALENT, limit 0.05). d_loc is the largest in the F3 series, but with 99 % saturated events E_loc (Moran's I is
  symmetric under χ ↔ 1 − χ) measures the clustering of the remaining 1 % of unsaturated places, i.e. "holes", not a
  node. New: the state survives switching off the source – saturation 99 % → 51 % → 34 % → 17 % → 6.8 % → 3.4 % in
  layers of 0.25 (the threshold alone 12.7 % → 0.9 % → 0.04 % → 0), the pulse carries the same, the pattern weakly
  similar to the pattern before switching off (r_s ≈ +0.1). Interpretation (descriptive): positive feedback C·Π with
  the pulse as a source gives a bistable medium – the whole region passes into the saturated phase with slow clocks,
  which sustains itself and drains only slowly (probably from the absorbing edge). Under homogeneous driving it floods
  the whole region; whether a bounded "droplet" of this phase can exist this test does not measure. The instruments
  passed, including bitwise agreement.
- **Threshold + feedback through the source removes the saturation** (F3.3f source, mechanism, 2+1, 128 sprinklings,
  s_eff = (s + P_in)·Π, FAILED according to the rule "saturation vanishes"). In the branch with the channel 0 saturated
  events in all 128 sprinklings (smallest Π 0.519 > Π_sat); the branches without the pulse (1.0 %) and with the
  threshold only (7.6 %) identical with part capacity (the same sprinklings). Slow clocks damp the source and the pulse
  below capacity before the threshold ignites the pulse – negative feedback as in F3.3 source. Step D is thereby
  closed: the positive channel (C·Π) gives a self-sustaining extended phase, the negative one (s·Π) nothing.
- **The bistable saturated phase has no boundary: a droplet either fades or floods the region** (F3.3g A, mechanism,
  2+1, 128 sprinklings, model F3.3f capacity without a field and without feedback on routing, ignition of a disc
  r0 = 0.15 (s_ign = 5, t ∈ [3; 4)) on a background s_bg = 0–0.8, FAILED for all values). At s_bg ≤ 0.2 the droplet
  vanishes (at 0.2 in about 2.5 time units); at 0.25 it fades over about 4.5 time units (93 % vanished by the late
  window), at 0.3 and more it fills the whole region (R_rms → 0.35–0.40, ring at the edge 18–94 % saturated). The
  transition lies between 0.25 and 0.3 and no standing droplet arises at it. What matters is how this happens: right
  after the ignition is switched off the droplet dissolves into a sparse cloud – R_rms jumps from 0.15 to ~0.31 within
  0.25–0.5 time units (roughly at the speed of light), the compactness B drops from ~0.8 to ~0.45 and the saturation
  inside C_drop from ~0.5 to 0.13–0.25 independently of the background. The saturated "phase" therefore does not form a
  sharp domain but a sparse cloud, whose density then falls (fades) or grows (floods) according to the background.
  Interpretation (descriptive): an event splits the flow uniformly among the descendants in the whole future cone, so
  the influence of every saturated region spreads at the speed of light and nothing holds a boundary; only a narrowing
  of the cones can limit this. Part B (perturbations) is not run: its premise (a droplet from part A) did not occur.
  The instruments passed.
- **Note on the instrument Δ_pers (open, not a verdict):** the branch without the pulse, i.e. the same model F3.2 on
  different independent sets of sprinklings, gives Δ_pers +0.006 ± 0.006 (F3.2), +0.008 ± 0.008 (F3.3),
  +0.000 ± 0.007 (F3.3d), −0.009 ± 0.008 (F3.3e search) and +0.021 ± 0.006 (F3.3e confirm). The scatter is larger than
  the quoted SE correspond to (χ² ≈ 11 around the weighted mean +0.008 at 4 degrees of freedom, p ≈ 0.03); the
  bootstrap SE may underestimate Δ_pers roughly 1.5×. This does not change the verdicts of F3 (all are EQUIVALENT with a
  margin), but claims of the type "x SE above zero" for d_pers and d_loc must be taken with caution until a calibration
  of the scatter on independent sets decides it.
- Local growth does not create spacetime: 0 of 7 rules (F6.1).
- Three-dimensional space did not arise from the dynamics (F6.2).
- The limit on the number of links (degree ≤ K₀) does not hold.
- The original local law Π = f(F/C) gives an acceleration ~1/r³; it was replaced by A5′.
- The pointwise A5′ is unprovable at an achievable N, noise decides (F1.6, E45b).
- A local greedy step does not follow the geodesic (F2.9).
- Saturation does not spread according to the uncorrected profile of the flow (F0(r_sat) = C).
- **Reading Π from the density of chains in 3+1 at the lowest density** (K.3, part d4): at ρ = 4·10⁶ an error of 0.114
  against a limit of 0.10. In 3+1 the reading of Π is usable only from ρ = 1.6·10⁷; in 2+1 it passed over the whole
  range.

## Instrument failures (no conclusion about the theory)

- **F5.2a: the free fall of a field defect in a 1+1 order of nonuniform density came out UNDECIDED, but the run had an
  implementation error, found after the run: the sprinkling had 0.808 of the density the rule assumes. No conclusion**
  (mechanism, 1+1, 128 sprinklings at rho 800, run by the author 10 Nov, 326 s; results/F5.2a/).
  - **Verdict by the frozen criteria:** UNDECIDED.
    - Instrument OK: I4 <= 1.4e-15; flat control |D| <= 0.018; counting null D -0.020 +- 0.013 against -2.14 curved.
    - a = D / D_cont 1.0174 +- 0.0023 (M 4), 0.9997 +- 0.0054 (M 8), 1.0304 +- 0.0143 (M 16: 0.059 > 0.05 with 2 SE,
      not FAILED).
    - eta (4, 8) +0.020 +- 0.004, (8, 16) -0.011 +- 0.012, (4, 16) +0.013 +- 0.015.
    - The launched arm (reported) gave D +0.717 +- 0.014 against +0.4235 in the continuum.
  - **The error:** field1d_curved.sprinkle_curved drew the number of points for the area of the trapezoid but placed
    them in the whole box [0, T] x [-X, X] and dropped those outside. The density was therefore rho (2XT - T^2)/(2XT)
    exp(2 Phi) = 0.808 rho exp(2 Phi) (measured 0.795 rho), while the rule normalises with rho. Every arm ran with an
    effective coupling 0.808 mu: kinks 11 % wider and lighter than the continuum references they were compared with.
    The development rows of the card had the same error.
  - **How it was found:** the launched arm stayed at ratio 1.7 against its continuum from rho 400 to 6400 (a bias that
    does not shrink with density). A deterministic quadrature of the same rule on a null grid agreed with the
    continuum (and so did two continuum readouts, to <= 0.2 %); the field of the Poisson version differed from the
    quadrature by up to 0.35 even in a flat order. The point count settled it: 318 per unit area at rho 400.
  - **Corrected (sprinkle_curved2: points placed in the trapezoid, then thinned):** in a development check (16
    sprinklings, rho 800, decides nothing), D / D_cont = 0.996 +- 0.005 / 0.994 +- 0.016 / 1.022 +- 0.040 and the
    launched arm 0.95 +- 0.09.
  - **Follow-up:** the revision F5.2a2 with the same frozen criteria and new seeds. F3.14a is not affected: its
    sprinkling (field1d.sprinkle_trapezoid) places the points in the trapezoid.
  - **Lesson for the instrument:** a test of a new sprinkling must count its points against the intended density
    before any dynamics; added to the instrument checks of future cards.

- **F5.1f-a: INSTRUMENT_FAILED (the proxy braked too little on these histories); descriptively the shared flow is
  not finite in 3+1: the same walker takes the same path in two regions in only 10 % of the histories** (mechanism,
  flat 3+1, density 1, box t [0, 27], |x| <= 10, |y|, |z| <= 8.5, 8 sprinklings of ~156 000 events, 128 histories of
  ~9 links (v 0 and 0.3), 8 walkers each, memory m = 8; rule: the flow shared equally among the children c of the
  last event with min L(p, c); run by the author 9 Nov, 10 662 s with 8 processes).
  - **Instrument:** the proxy of F5.1d on the same histories braked by Delta_h = −0.092 ± 0.034 at v 0.3 (+2 SE =
    −0.024, the criterion needs < −0.05); in F5.1-diag3c on shorter histories (~4 links, history velocity 0.278) it
    braked by −0.243 ± 0.036. The threshold was calibrated on different histories; this is a design error. By the
    criteria the run gives no conclusion about the rule.
  - **Descriptive (decides nothing):** M5 identity (walker 0 in region A = in region B): 0.102 of 128 histories (v 0
    0.078, v 0.3 0.125), against the required 0.95. Per step: 316 children of v in the window, of which 120–127 have
    their diamond inside the box; the admissible set has 6.3–6.4 members; the smallest deficit is > 0 on average 0.41–
    0.42 (in 26 % of the walkers at v 0 its mean exceeds 0.5), while for the launch histories of F5.1-diag4 it was 0 in
    every case. The flow disperses fast: per-walker SD of the velocity 0.37–0.38, spread at the end of 10 units 3.2–
    3.5 along the throw and 3.7–3.8 across. Centre: Delta_h −0.041 ± 0.011 (v 0) and −0.029 ± 0.017 (v 0.3); at v 0 it
    depends on the sign of the throw direction (−0.007 for +x, −0.075 for −x), i.e. a drift fixed in the laboratory.
  - **Reading (a hypothesis):** equal sharing among the children with the smallest coarse L lets the walkers step to
    strongly boosted children (rapidities up to ~3 in F5.1-diag4); after a few steps their own history p -> v is no
    longer a timelike longest chain, the deficit rises above 0, the admissible set is then decided by which children
    the box admits, and the paths depend on the region. One link carries no proper time (a link has no events inside),
    so L(p, c) does not distinguish an aligned child from a boosted one well enough to confine a flow of single
    events. A rule for a free flow needs either weights that suppress boosted continuations by an invariant of the
    order, or a flow that is a bundle of events rather than a single chain (F3.6); both are new hypotheses (R132).

- **F3.13b: the readability requirement failed; every readable shell shows a locked charge and the hedgehog keeps a
  core cost that the twist loses** (mechanism, 3+1, box t [0, 20], half 7, ~54 800 events, ~1 050 sources,
  16 sprinklings, INSTRUMENT FAILED 31 Oct; the verdict and the criteria stand). New rule approved by the author:
  alignment of orientations, n ← N[(1 − κ) n + κ m̂] towards the flow-weighted mean of the causal past (the gradient
  flow of the configuration cost E = Σ(1 − n·n′)).
  - **What failed:** the hedgehog arms (all κ together) were readable on 0.445 of their triples against 0.5 (κ = 0:
    0.26; 0.1: 0.46; 0.3: 0.62). The requirement was a fraction; the taller box added late slabs in which the cold
    bundle has thinned. The other instrument conditions passed: vacuum 0 on 335 of 335 in every κ; the twist at
    κ = 0.3 relaxed, mean cost 0.0277 → 0.0097. Revision F3.13b2: a count instead of a fraction (≥ 50 readable shells
    per arm), decision on κ > 0, new sprinklings.
  - **Descriptively (not evaluated):**
    - The hedgehog read +1 on 89 / 153 / 207 of 89 / 153 / 207 readable shells (κ = 0 / 0.1 / 0.3); the erase arms
      +1 on 81 / 123 / 186 of all of theirs, on shells R ≥ 3 after the erase 25 / 59 / 115 of 25 / 59 / 115; no core
      event in any arm.
    - Core cost at κ = 0.3 in the last slab: hedgehog 0.343 ± 0.049, twist 0.011 ± 0.003 (30×). Over time the
      hedgehog's core cost grows (κ = 0: 0.27 → 0.63; κ = 0.3: 0.23 → 0.34–0.41) while the twist's falls: the cost
      concentrates in a shrinking core, as Derrick's theorem expects for a directional field.
    - After the erase the core cost jumps (0.54–0.59) and relaxes back; at κ = 0.3 it ends at 0.31 against 0.34 of the
      undisturbed hedgehog: the core re-forms with about the same cost, as if the cost belonged to the charge rather
      than to the history of the core.

- **F3.13a flow_d3: in 2+1 too few rings were readable; every readable ring read the right charge** (mechanism,
  2+1, box t [0, 24], half 10, ~9 600 events, ~156 sources, 32 sprinklings, INSTRUMENT FAILED 31 Oct; the verdict and
  the criteria stand). The defect arm was readable on 30 % of its (sprinkling, slab, ring) triples against a limit of
  50 %; the limit came from a single development sprinkling (59 %), a design error. Descriptively: defect +1 on 258 of
  258 readable rings, anti −1 on 258 of 258, vacuum 0 on 606 of 606, twist 0 on 459 of 459, at every slab t = 5–21;
  double +2 on 14 of 15; no core event. After merging a cold bundle keeps few world lines, so rings of width 1.5 catch
  few carriers. A revision would read one wide annulus around the centre (any closed loop around the core carries the
  winding); not run, since flow_d4 (3+1) passed.

- **F3.12a: the instrument failed on its own tolerance; the numbers point clearly to capacity, not topology, as the
  source of the resistance to rerouting** (calibration, 2+1, box t [0, 10], |x| ≤ 10, ~4 000 events, ~85 000 links,
  32 sprinklings, lens 3 ≤ t ≤ 7, INSTRUMENT FAILED 31 Oct; the verdict and the criteria stand).
  - **What failed:** the unperturbed program returned the original flow to 1.3e−7 of the total link flux and the
    balance held to 4.6e−9, against limits of 1e−9. The limit was a design error: the linear program (HiGHS) has a
    feasibility tolerance of 1e−7. All programs were solved. Revision F3.12a2: only the tolerances (1e−6), new
    sprinklings.
  - **Descriptively (not evaluated):** I_cascade, the change of the flow beyond the flow that is forced to move, at
    σ0 = 0.1 / 0.2 / 0.3:
    - T (narrowed cones, fewer paths): 0.0043 / 0.0129 / 0.0194 (± 0.0002–0.0005);
    - K (capacity cut to the same forced change, same paths): 0.0391 / 0.0539 / 0.0618 (± 0.0009–0.0011);
    - P (throughput: narrowed cones and added events, 648–1594 new links): 0.0041 / 0.0115 / 0.0167, like T.
    By the rule declared for the verdict this would read "capacity dominates" at every strength (K about 3–9× T).
    I_forced at σ0 = 0.2: 0.105 (T, K), 0.117 (P).
  - **Reading:** confirmed by the revision F3.12a2 on new sprinklings (Does not hold, above).

- **F5.1 carried: the instrument failed narrowly; a flow with a carried direction alone does not fall** (mechanism, 2+1,
  order of F7.1b σ0 = 0.45, 32 + 32 sprinklings of ~47 650 events, 40 chains, INSTRUMENT FAILED 29 Oct; the verdict and
  the criteria stand).
  - **What failed:** the control of the longest chains in the loaded order R_max = 0.630 ± 0.074, the lower bound 0.482
    is just below 1/2 (criterion 1/2 < R_max ∓ 2 SE < 2). In part intrinsic (other sprinklings) R_max = 0.833 ± 0.076
    came out; the difference 0.20 is about 1.9 σ. Flat order OK: R = +0.003 ± 0.042, R_max = −0.13 ± 0.12.
  - **Descriptively (not evaluated):** the flow with a carried direction in the loaded order R = −0.106 ± 0.038, shift
    +0.18 over 22 (a fall would be −2.7), v0 = −0.005 ± 0.005; 25.1 of 40 chains reached the end. The flow keeps its
    direction and does not perceive the gravity of the order.
  - **Reading:** a carried direction in flat coordinates is inertia without geometry; the author himself expected this
    in advance ("lost in math"). The control R_max evidently lies between 0.6 and 0.85: the longest chains read the
    geometry, but not completely (discreteness of the chains). The band of the control stays for F5.1b and F5.1c as it
    was set.

- **F6.8b: the null k = 0 grew with the age, the mixed sector undecided** (chain, 28 Oct, INSTRUMENT FAILED; the verdict
  and the criteria stand). Null J_ab = +5.9e−7 / +9.3e−6 / −4.2e−4 at T = 40 / 160 / 640. Causes: the boundary term of
  the age cut on 8 nodes of the panel (error ~T²) and the rounding of the difference of delays in thin intervals (~T⁴).
  The criterion of the null divided by the smallest |J_ab| (a design error). The fit K0_ab = 1.78 is an artefact. Fixed
  in the revision F6.8b2.

- **K.17: the quadratic form of the BD action cannot be measured directly on a 3+1 sprinkling in the box of the
  laboratory; the noise is 10× to 10⁴× above the limit of feasibility** (calibration, 28 Oct, box t ∈ [0, 8], |x| ≤ 8,
  ρ = 1, a = 0.2, w = 2.5, rows t ∈ [6.5; 8], r ≤ 6, ~1360 rows, 48 sprinklings of ~33 000 elements, FAILED; instrument
  OK: the configuration a = 1e−9 gives S exactly equal to S(0) for all ε).
  - Q_pred = 2.73 (Einstein–Hilbert). Number of sprinklings needed for 10 % on the variance of
    Q = (S₊ + S₋)/2 − S₀: ε = 0.3: sd(Q) 1045, M = 1.5·10⁷; ε = 0.1: sd 140, M = 2.7·10⁵; ε = 0.03: sd 28,
    M = 1.1·10⁴ (limit 1000).
  - Correlation of S₀ and S₊ only 0.55 / 0.67 / 0.85: the thinning e^{3af} changes the elements, so correlated sampling
    cancels little noise.
  - The mean Q is biased by orders of magnitude: 2623 / 473 / 76 against 2.73. According to the criterion it is not
    evaluated (relation through a straight segment, non-locality, edges of the box); it falls with growing smearing.

  Reading: a direct measurement of c_ΦΦ, c_ΦΨ, c_ΨΨ from a sprinkling is out of reach of the laboratory. Nothing follows
  for the theory. F6.8 therefore computes the mean value of the action deterministically (F6.8a), so it does not depend
  on K.17.

- **F6.4a: the calibration of 2D orders with the BD action reproduces the manifold-like phase and the transition, but
  not the height of the crystalline phase** (calibration, N = 50, ε = 0.12, 20 values of β, 4 chains of 500 + 3000
  sweeps, FAILED; checks of the instrument OK: exact equilibrium N = 4 TV 0.0078, β = 0 against 20 000 directly sampled
  orders within the tolerance). Against Surya 2012:
  - phase I (β = 0.1): χ 0.4981 (0.499), h 10.237 (10.232), S 3.835 (3.846), d_MM 2.005: agreement to hundredths;
  - transition: the largest fluctuations var(S) and the longest autocorrelation at β = 2.4 (published 2.2–2.5);
  - phase II (β = 3.5): S −41.07 (−41.37) and χ 0.593 (0.579) within the tolerance, **height h 3.70 against 4.18**
    (tolerance 0.3).

  In phase II 2 % of the moves are accepted and χ is not monotonic in β (0.593 / 0.576 / 0.605 for β = 3.5 / 4.0 / 4.5),
  which points to frozen, preparation-dependent states. The difference in height may therefore be a property of the
  preparation (thermalisation, initial states), not an error of the action. By R132 the verdict holds, the tolerance is
  not changed retroactively. The phase that F6.4b needs (manifold-like) and the position of the transition fit very well
  descriptively; formally the instrument is not calibrated.
- **F6.4c: the height of the crystalline phase does not depend on the length of thermalisation nor on the start; no
  preparation heads towards the published 4.18** (calibration, β = 3.5, 5 preparations of 8 chains, 3000 measured
  sweeps, UNDECIDED according to the criteria: the largest difference of the long preparations 0.179 above the limit of
  independence 0.15, but far below the limit of dependence (difference − 2 SE = +0.010 against 0.3); instrument OK:
  repetition of F6.4a h 3.742, S −41.17, exact equilibrium TV 0.0088). Height h:
  - P1 repetition 3.742 ± 0.038; P2 long thermalisation 3.784 ± 0.032; P3 cooling from β = 2 3.887 ± 0.099;
    P4 start from a chain 3.750 ± 0.025; P5 start from an antichain 3.708 ± 0.008; shift between the halves of the
    measurement ≤ 0.014.
  - The distribution is bimodal: most chains 3.68–3.92, two cooled chains (P3) 4.31 and 4.37 and one from P1 4.00.
    Phase II therefore has at least two long-lived states and cooling sometimes ends in the higher one.
  - The mean of the long preparations 3.78, 0.40 below 4.18. The most stable preparations (P4, P5 from opposite
    extremes) give 3.71–3.75 with a small error. No preparation and no single chain (at most 4.37) reproduces the
    published value as a mean.

  Reading: the difference of F6.4a was not caused by a short thermalisation (P2 with a tenfold one gives the same as
  P1). The equilibrium height of phase II in the code of the laboratory is ~3.7–3.8. The value 4.18 in the paper
  corresponds rather to a preparation that ended in the higher state, or to another definition of the height; this
  cannot be decided from the available text. Phase I (height included) and the transition fit to hundredths, so the
  action, the moves and the measurement of the height are fine. The discrepancy is limited to the metastable
  crystalline phase.

- **F1.9a: β of the discrete operator outside the source, INSTRUMENT FAILED only in the check of the quadrature at
  σ = 4** (mechanism, 25 Oct, 252 jobs, 5.8 h). The only unmet condition: |c_main − c_control| in the row 6σ at σ = 4
  is 0.090 (H1) and 0.076 (γ = −1) against a tolerance of 1/16. At σ = 8 and 16 it is 0.023 and 0.004. The wall
  (difference ≤ 0.002), the control γ = −1 (deviation 0.014 and 0.012 at σ = 16) and the gate of the instrument (K.13a
  lin PASSED with the same core) are fine. According to the protocol β is not decided. Descriptively, without a
  verdict:
  - H1, decision rows: c_disc(σ = 4, 8, 16) at 4σ −0.114, −0.029, −0.002 against the continuum value of the mixture
    +0.008; at 6σ +0.016, +0.047, +0.055 against +0.058. The deviation falls by a factor of 3.3–3.9 per doubling of σ
    (roughly (ℓ_k/σ)²), so it would otherwise meet PASSED (≤ 1/4 at σ = 16, does not grow). β_disc at σ = 16: 1.001 and
    0.973. The deviation from 1 at 6σ is a property of the mixture (continuum value 0.971, |∇L| falls with the exponent
    −2.09 instead of −2), not of discreteness.
  - The control γ = −1 goes to 2 in the same way: 4σ +2.081, +2.045, +2.022 (continuum 2.008).
  - Row 3σ (only reported, tail of the source): H1 −0.419, −0.281, −0.224 against −0.195, the deviation falls.

  Reading: the coarsest point of the ladder failed (a field of width 4 ℓ_k with the control quadrature
  [80, 12, 8, 3, 8, 4]). By R132 the criterion is not changed retroactively. A formal verdict on β would require a new
  test with a protocol given in advance (e.g. a ladder σ = 8, 16, 32 with the same bands). The value of such a test is
  confirmatory, the band Δc = 1 is far from the measured numbers.

- **F6.3b A: in the order from narrowed cones a compact, slowly growing domain of the saturated phase arises, but only
  with an external flow; a standalone droplet does not arise** (mechanism, 32 sprinklings, N = 36 000, INSTRUMENT
  FAILED: the ensemble of solutions is not stationary at s_bg = 0.1, C_late differs between the halves of the averaging
  window by −0.117 ± 0.014 against a tolerance of 0.05). Verdicts by background (criterion of F3.3h A):
  - s_bg = 0: the dynamic order VANISHED (100 %), but later than the fixed one: the domain holds until t ≈ 6–7 (the
    fixed one until ≈ 4.5).
  - s_bg = 0.1 and 0.2: UNDECIDED. The domain survived on all sprinklings and is a domain of the phase (C_drop late
    0.69 ± 0.01 and 0.84 ± 0.004, compactness B 0.87–0.95), but it keeps growing: the growth of R 0.035 ± 0.002 and
    0.039 ± 0.001 is right around r0/4 = 0.0375 and the radius keeps creeping up at the end of the window.
  - s_bg = 0.3 and 0.4: FAILED (growth above r0/4 and ignition of the background, ring 0.08 and 0.31).

  The fixed order with the same driving dissolves into a sparse cloud (as F3.3g); the dynamic order has 7–43 % more
  links than vacuum. Descriptively: the narrowing of the cones turns the dissolution at the speed of light into a compact
  domain that grows at a speed of ~0.01. The reading is limited (axiom audit 22 Oct): the cones are drawn in coordinates
  (F0.6, quick runs: frame-dependent), the pulse acts as a source of flow (in contradiction with decision 4) and the
  solution is self-consistent over the whole spacetime. Only the negative reading is therefore valid: even with these
  advantages no droplet arises without an external flow (s_bg = 0). The domain at s_bg > 0 is, according to decision 3,
  a structure maintained by the external flow. A repetition with a fixed stationarity instrument is not recommended
  until F0.6 passes.
- **F6.3 V2: the limit Π_sat → 1 differed by one link because of a fault of the old instrument at the wall**
  (validation, FAILED; revision F6.3b). At Π_sat = 1 the dynamic solution (the vacuum order of the exact construction)
  differed in 1 of 4 sprinklings from the fixed order of core.order.links by 1 link out of 1 966 510. Analysis: the pair
  (4725, 5857) lies at the wall of the cylinder (r = 0.492 and 0.498), half of its Alexandrov interval is outside the
  sprinkled region, effective ρV = 15.4 (nominally 30.2), direct count 0 – it is a real link. core.order.links selects
  candidates by the nominal volume of the interval in an unbounded region (λ = 25), and therefore omits such links at
  the wall; the exact construction finds them. F6.3 V1 passed only because its 4 sprinklings had no such link. The fault
  concerns only the outer wall (beyond the absorption radius 0.4 in the F3 tests, where the flow leaves), so it
  practically does not affect earlier results; core.order.links stays unchanged. Part (b) passed: windows 1.5 and 3.0
  give identical links on the solution with matter (~1.96 million) in all 4 sprinklings. The revision F6.3b uses the
  exact construction in all branches and V1 admits only empty links at the wall.
- **K.12 s16: the source → field link in 2+1 is not calibrated at σ = 16 ℓ_k for a test of the F2.7 type, but the
  deviation falls with σ** (calibration, self-similar family σ = 16, q = 0.0025, R = 200, FAILED). The geodesic in the
  field of the causal set P_cs = 3.771 against the continuum P_a = 3.354, ratio 1.124 (tolerance P_a/16, i.e. 1.0625);
  at σ = 8 (field of F2.7) 1.206. The excess thus fell from 0.206 to 0.124 when σ doubled (formally ~σ^(−0.73); whether
  it falls to zero, s32 decides). The quadrature converged (control and main level differ by 0.0001 against
  P_a/64 = 0.052), Π(0) 0.7202 against 0.7099. Ratio of the slopes cs/continuum at r = 1.25σ / 2.5σ / 3.75σ: 1.158 /
  1.097 / 1.049 (at σ = 8: 1.245 / 1.129 / 1.056). The Gauss deficit of the flux 1 − flux(R/2)/flux(3σ) stays ~0.10 at
  both σ (continuum −0.014): in the self-similar family R = 12.5σ, so an edge effect would scale with σ and not
  diminish; s08_R200 distinguishes it.
- **K.12 s08_R200 and s32: the excess of the field of the causal set in the region of motion is not a wall effect and
  according to the rule given in advance has a floor of ~10 % (s32 not yet certified)** (calibration, 2+1, equation
  (a), 24 Oct). s08_R200 FAILED: P_cs/P_a = 1.227 at the wall R = 200 against 1.206 at R = 100. The reading "wall" = not
  supported: in the same region r = 48 the Gauss deficit 0.088 → 0.074, excess 0.206 → 0.227 (both changes below a
  quarter). Closer to the wall the deficit depends on the wall (r = 64: 0.190 → 0.110; r = 80: 0.414 → 0.143), but the
  flux of the causal set falls with r even at R = 200 (0.116 at r = 24 → 0.099 at r = 80, continuum 0.106 constant):
  the field outside the source is not harmonic.
  s32: P_cs/P_a = 1.105 (quadrature agreeing to 0.0003 against a tolerance of 0.105), but INSTRUMENT FAILED, because the
  control solution from vacuum did not converge within 30 iterations. Reading "convergence" (preliminary): excesses
  0.206 / 0.124 / 0.105 for σ = 8 / 16 / 32, Aitken floor 0.099 > 1/16 → **floor**. The depth of the core Π(0) moves
  away from the continuum value 0.7099 with growing σ (0.7113 / 0.7202 / 0.7256). Gauss deficit at R/2: 0.098 / 0.095 /
  0.135. Reading according to the formulation given in advance: if the floor is confirmed, the discrete operator in
  this regime (2+1, equation (a), self-consistently, strong field Π(0) ≈ 0.71) does not converge to the assumed
  continuum one; the wall in the region of motion is excluded; the cause is open (boundary conditions, normalisation of
  the source, scaling of the family, or the operator). A difference independent of σ in a self-similar family points
  rather to another continuum equation than to a local correction of order ℓ_k/σ.
  Next (declared 24 Oct before the runs): s32 again with max_iter = 60 (certification) and the weak-field family
  w08/w16/w32 (q ten times smaller, Π(0) ≈ 0.97): is the floor a property of the linear operator, or of the non-linear
  self-consistent coupling?
- **Decomposition of the K.12 "floor": non-local screening by the operator is excluded, the floor is a sum of
  components of the instrument and an open non-linear remainder** (analysis and linear model, 24 Oct,
  `docs/derivations/k12_floor.md`, scripts `docs/derivations/k12_diag/`; statuses and reading of K.12 unchanged, R132).
  Exact symbol of the expected BD operator 2+1 on a static field: B(k) = −k²(1 − 1.155 ℓ_k k + …), M₀ = M₂ = 0,
  r(0) = 1. The distant field sees the whole charge and the operator anti-screens at the source ∝ ℓ_k/σ (without a grid
  and a wall: flux at 1.25σ 1.159 / 1.077 / 1.038 for σ = 8 / 16 / 32). Linear model of the K.12 scheme (the author's
  runs, fine grid 1 ℓ_k): slope at the source 1.161 / 1.077 / 1.037, i.e. the operator to 0.002, and the wall does not
  reach the region of motion. On the self-similar grid of K.12 (h = σ/4) +0.007 / +0.021 / +0.051 are added: the grid
  error grows with σ, the response of the operator falls, hence the excess stops (Aitken 0.099 is a transition, not an
  asymptote). The Gauss deficit at R/2 ≈ 0.10 independently of σ and on the fine grid, even without the PL basis; it
  vanishes with a local response of the geometry. Origin: the non-local operator sees the kink of the field continued
  by zero beyond the wall (in 2+1 it carries the whole flux), while the reference is a local Dirichlet problem. The far
  field of the non-local variants is not converged in the grid step at the wall (σ = 8: flux 0.9R 0.494 → 0.562 at step
  1 → 0.5). The same mechanism explains the depth of the core, which moves away with σ (model 0.978 / 0.924 / 0.895,
  K.12 0.995 / 0.965 / 0.946). The definition of the source is fine (the variant with a local geometry converges: u(0)
  1.0006 at σ = 32). Open remainder: the slope at 1.25σ in K.12 is 0.077 / 0.060 / 0.055 above the linear model on the
  same grid (strong field, non-linearity). Caveat: the model represents the non-local response of the geometry by a flat
  operator on a spline. Conclusion: K.12 is neither a falsification of the continuum limit of Path A nor its
  confirmation. A strong and a weak field with a step at the source ≤ 1 ℓ_k decide.
- **K.12 h16: a step of 1 ℓ_k at the source removed the grid component, as the decomposition predicted; the field still
  did not pass** (calibration, 2+1, σ = 16, q = 0.0025, R = 200, 76 nodes, 25 Oct, FAILED). P_cs/P_a = 1.111 (s16 with
  the grid σ/4: 1.124; tolerance 1.0625). Slope at 1.25σ 1.141 against the prediction 1.137 (operator +0.077 and the
  non-linear remainder +0.060 from the analysis). Reading 'grid' (declared 24 Oct) = **supported**: difference of the
  grids +0.018 against the linear model +0.021. Quadrature agreeing (|L2 − L1| = 0.0000), the solution converged
  (6 iterations). Gauss deficit at R/2 0.077 (s16 0.095), Π(0) 0.7193. After subtracting the grid an excess above the
  tolerance remains at σ = 16. It consists of the operator (∝ ℓ_k/σ) and the non-linear remainder, which h32 and the
  weak-field hw16/hw32 decide.
- **K.12 closed on 25 Oct without h32 and the weak-field family** (a decision after the analysis; the run of h32 was
  interrupted by a computer crash). For the theory it gives: the discrete operator converges to the continuum one and
  does not screen (analytically, r(0) = 1), the floor was largely an error of the instrument (grid, wall) and the
  prediction of the decomposition was confirmed on h16. What stays open is the non-linear remainder ~0.06 in the slope
  at 1.25σ in a strong field of 2+1 with the old equation (a). Statuses of the parts unchanged. Lesson for further
  fields: a step at the source fixed in ℓ_k, a distant wall or an outer continuation of the field, a reference with the
  symbol of the operator.
- **F1.2, integral A5′ with a sharp source**: residual of the expected field 0.0875 q against a tolerance of 0.05. A
  sharp edge of the source cannot be represented at the scale of the non-locality ℓ_k. I set the tolerance without a
  measurement on the full grid (my error). The revision F1.2b changed only the shape of the source to a Gaussian, the
  criteria stayed the same. Descriptively: A and B were met, the null model failed (z up to 21), but without a valid
  instrument this does not count.

- **The curved operator in 3+1 does not converge in the distant field** (K.11b far_d4): volume-weighted remainder of
  the distant rows Q_far/Q_src = 3.5, 2.3, 1.9 at three levels of the quadrature (tolerance 0.01). The spurious source
  is thus larger than the real one and falls only slowly with resolution; F1.4 and F1.5 in 3+1 cannot be computed with
  this instrument. Cause found during the development of K.11c (19 Oct, development run, not a test): the remainder
  does not arise from the long intervals themselves, but from the piecewise linear interpolation of the geometry Π
  between nodes. The operator responds to the curvature of the geometry (the term −R/2, i.e. ∇²Π^D) and a linear
  interpolation has ∇²Π at the nodes as delta functions; long nearly light-like intervals are thin needles that collect
  them, in 3+1 with a weight ~r. Row r = 40 l_k: linear geometry +3.4e-3 and on doubling any parameter of the
  quadrature 1.3e-3 to 4.8e-3; smooth geometry −4.15e-5, stable to 0.2 %; C² spline through the same nodal values
  −4.152e-5; C¹ spline (PCHIP) +4.4e-4. The fix of the instrument is in K.11c (geometry as a C² spline); K.11b stays
  FAILED according to its criterion. In 2+1 the linear and the smooth geometry differ in the distant rows by about 1e-4
  (r = 30, grid of F1.5); K.11c far_d3 measures how much this makes in the volume-weighted remainder, i.e. how large a
  caveat it brings to F1.4 and F1.5.

## Discoveries

- **Growth probe (direction A, 1+1): with the causal growth rule "a candidate event is realised with probability Pi^2,
  ln Pi = -(kappa/2 rho) sum over the realised causal past of the load 2V", matter thins its future order, the counted
  proper times shrink near it, and a second defect falls towards it in proportion to its mass (1 : 1.97 : 3.86 for
  M_A 4 : 8 : 16), independently of its own mass, and the two attract each other mutually (M_A dx_A = M_B dx_B)**
  (development, DESCRIPTIVE, decides nothing; approved by the author 11 Nov 11:59 with his corrections; diag/growth_probe.py;
  rho 800, kappa 0.003, kinks A at -3 and B at +3 as two separate fields interacting only through the order, T 20,
  32 sprinklings per combination, each paired with kappa = 0 on the same candidates and marks; run in the container
  11 Nov, 925 s; diag/out/growth_probe_20261111T111825.json).
  - **Rule and normalisation (after the author's correction):**
    - every candidate carries an independent mark u fixed at the sprinkling, and is realised iff u < e^{2L};
    - L = -(kappa / 2 rho) sum over the realised y < x of q, with q = 2 V_A + 2 V_B and rho the candidate (= proper)
      density;
    - the realised density rho Pi^2 carries the factor Pi^2, so the continuum law is (d_t^2 - d_x^2) L = -kappa Pi^2 q,
      i.e. box_g L = -kappa q in ds^2 = Pi^2 (dt^2 - dx^2);
    - the fields use the unchanged rule of F3.14a on the realised events.
  - **A0 / A1:**
    - a change of later marks leaves every earlier event unchanged;
    - two linear extensions give the same realised set and fields to <= 3.7e-13.
  - **A3, density and intervals:**
    - realised fraction / mean Pi^2 per x bin = 1.000 +- 0.002;
    - for realised timelike pairs near A, the realised count inside the interval / (rho x proper area) = 1.004-1.006,
      against 0.49-0.72 for the coordinate area;
    - the counting measures the proper times of the metric Pi^2 (dt^2 - dx^2), so clocks near matter slow down.
  - **A4, continuum** (three coupled wave equations): L at the end within 3-7 % of the continuum (partly the slab average
    of a growing field). The displacements below are within 2-5 % of the continuum.
  - **A5, active source** (M_B = 8): B moves towards A by
    - -0.571 +- 0.020 (M_A 4; continuum -0.543);
    - -1.126 +- 0.025 (M_A 8; -1.095);
    - -2.201 +- 0.030 (M_A 16; -2.142).
  - **A6, universality** (M_A = 8): B moves by -1.086 +- 0.008 (M_B 4) and -1.126 +- 0.025 (M_B 8). The continuum gives
    -1.068 and -1.095 (finite size); the ratios to it are 1.017 and 1.029.
  - **Mutual attraction:** A moves towards B by +0.61 (M_B 4) and +1.10 to +1.17 (M_B 8). M_A dx_A = M_B dx_B within
    the noise (4.4 / 4.6; 17.6 / 17.6), as Newton's third law requires.
  - **A7, null:** every displacement is relative to the paired kappa = 0 run.
  - **A2, covariance with a boosted launch:** not yet run.

  Reading (development).
  - All three levels of the author's list are met descriptively:
    1. growth: the source changes the probability of realising future events;
    2. active response: the change of the order moves a test defect;
    3. universality: the active mass is the energy of the source, and the acceleration does not depend on the test
       mass.
  - With F3.14b and F5.2a3, inertial, passive and active mass come from one field and one causal growth rule, in 1+1.
  - The rule is causal, independent of the order of computation, and uses only the order, counts and marks. It is a
    candidate for the Axiom of the dynamics of the order.

  Limits:
  - 1+1: the response is a scalar (R2), with no 1/r; the potential of a static source grows linearly in space and time;
  - in 1+1 the cones need not change (conformal flatness); in 2+1 the relation itself will have to change;
  - the discrete falls exceed the continuum by 2-5 %, not yet explained;
  - kappa is a chosen coupling (G is not derived).

- **Wandering study (1+1): a field defect at rest on the causal order makes a random walk in momentum, with a diffusion
  constant D_p proportional to 1/rho (= l^2) and rising steeply with the mass, D_p ~ 2.7 l^2 / w^5 = 2.7 l^2 (M/8)^5**
  (development, DESCRIPTIVE, decides nothing; diag/wander.py; flat 1+1, rule of F3.14a unchanged, kink at rest, T 30,
  trapezoid |x| <= 40 - t, mu 0.25 / 1 / 4 (M 4 / 8 / 16, w 2 / 1 / 0.5), rho 200 / 400 / 800, 32 sprinklings each (up to
  1.2 million events per run); run in the container 10 Nov, 596 s; diag/out/wander_20261110T171302.json).
  - **Time law:**
    - the spread of the position grows with an exponent 1.45-1.55 between dt 9 and 29 (1.32 for M 16 at rho 200, where 10
      of 32 fits are lost), i.e. Var x ~ t^3, a random walk in momentum;
    - the earlier exponent 1.4 (to dt 9) was a mixture with early terms;
    - the variance of the fitted local velocity grows linearly in t.
  - **Two independent readouts agree on the diffusion of the velocity, D_v** (from Var x = (2/3) D_v t^3, and from
    Var v = 2 D_v t):
    - mu 1: 2.2e-4 / 1.8e-4 (rho 200), 1.06e-4 / 1.19e-4 (400), 5.0e-5 / 4.5e-5 (800);
    - mu 4 at 800: 2.38e-4 / 2.45e-4.
  - **Density:** D_p rho is constant within each mass: mu 1: 2.84 / 2.71 / 2.57; mu 4: 45.9 / 48.1 / 48.7; mu 0.25:
    0.06-0.10. So D_p ~ 1/rho: Poisson noise of the order.
  - **Mass:** D_p rho = 0.09 / 2.7 / 48 for M 4 / 8 / 16, an exponent 4.2-4.9 per doubling. Dimensional analysis with
    the width as the only scale (M w = 8) gives M^5.
    - One dimensionless law, sd/w = A (rho w^2)^(-1/2) (t/w)^(3/2), collapses all nine cases with A = 0.16-0.18 (0.12-0.13
      for the poorly resolved M 16).
    - The heavier, narrower kink sees the graininess more.
  - **No mean drift** (all within 1 SE of zero).

  Reading (development). A bound object on a Lorentz-invariant discrete order receives random kicks of momentum from
  the discreteness. Their rate grows as the object becomes small compared with the discreteness (rho w^2). This is the
  "swerves" of causal set theory (Dowker, Henson and Sorkin 2004), here for a bound object, with a measured law, not a
  postulated diffusion constant.

  Limits:
  - 1+1 only; the law is dimension-dependent, and no extrapolation to 3+1 is justified;
  - the potential is chosen;
  - whether the kicks come from the dynamics of the defect or from the noisy radiation of its own field is not
    separated.

  **Practical consequence:** the wandering sets the resolving power of every test of this branch. The replications of
  F3.14a and F5.2a2 failed because 16 development sprinklings underestimated it. The revisions must size their samples
  from this law with a margin.

- **F5.1-diag2: the longest-chain rule of F5.1c/d/e has no inertia of its own; the velocity of the history goes to the
  velocity of the region in which the rule reads the order** (diagnostic, DESCRIPTIVE, decides no registered test;
  flat 3+1, density 1, box t [0, 28], |x| ≤ 11, |y|, |z| ≤ 7.5, 8 sprinklings of ~138 500 events, 16 starts each, the
  registered rule m = H = 8 unchanged, only the cylinder of LocalLinks moves with u along the throw; the cylinder at
  rest reproduces the registered rule in 8 of 8 checks; run by the author 8 Nov, 10 019 s with 8 processes).
  - **Late velocity (s ∈ [6, 12] after the start of the walk, along the throw; SE over sprinklings), arms (v, u):**
    R = 4: (0, 0) +0.025 ± 0.025, (0.3, 0) −0.027 ± 0.027, (0.3, 0.3) +0.321 ± 0.029, (0, 0.3) +0.305 ± 0.032;
    R = 6: −0.014 ± 0.026, +0.049 ± 0.026, +0.351 ± 0.025, +0.272 ± 0.033. The thrown history had +0.237 (measured
    velocity of the history with the nominal 0.3), the history at rest +0.052.
  - **Paired differences (same start, same history):** co-moving − at rest +0.348 ± 0.042 (R 4) and +0.302 ± 0.025
    (R 6); pulled − at rest +0.280 ± 0.047 and +0.286 ± 0.043; thrown − at rest in the cylinder at rest −0.052 ± 0.025
    and +0.063 ± 0.031. The history takes over ~95 % of the velocity of the cylinder and keeps none (R 4) or at most a
    third (R 6) of its own.
  - **At the wall (within 0.5 of the cylinder surface, mean per step):** horizon set Q 32 % (R 4) / 21 % (R 6), tie set
    of min L(p, q) 34 % / 26 %, the chosen target 71–76 % / 53–58 %, mean distance of the target 0.90 R / 0.87 R.
  - **Reading against the outcomes declared before the run:** the cylinder is the medium (both paired differences
    > 0 by 7–12 SE). Not as declared: the effect is not smaller at R = 6; the velocity locks to u at both radii. The
    laboratory itself is not the medium, the time top of F5.1e was not the cut (the thrown history loses its velocity
    in the cylinder at rest also here).
  - **Mechanism (inferred, not tested separately):** both criteria of the target are minima, min L(p, q) and then the
    fewest events n(p, q). A diamond between p and q that leaves the region is cut by it, which lowers both L and n.
    The rule therefore prefers targets whose diamonds stick out of the region (the target is at the wall in 3/4 of the
    steps, while only 1/3 of the candidates are), and the velocity follows the region. In the continuum the volume of a
    diamond does not depend on its boost; the preference comes only from the cut. The horizon set (all q with
    L(v, q) = H) is unbounded along the light cone, so any finite region cuts it.
  - **Consequences:** the "medium at rest" of the 3+1 fall (F5.1-diag) is the cylinder fixed at the apex; together
    with the turn at the apex it gives the terminal velocity of F5.1e. In 2+1 (F5.1d) the cylinder was also at rest in
    the frame of the source; the fall there was measured from rest in that frame, so F5.1d PASSED does not show that the
    rule is covariant (R132: the verdict stays). A rule of motion that does not depend on a region (or reads it in the
    frame of the history itself) is a new hypothesis with its own criteria (F5.1f), and inertia in the flat order must
    be its first instrument.

- **F3.8 search: a retarded field with range contracts the bundle for the first time, but no bound state arose:
  ε = 0.1 slowly collapses, ε = 0.5 collapses and rebounds** (mechanism, exploration, 32 sprinklings 2+1,
  ~49 200 events, RECORDED; candidate according to the rule in the header: none). Field φ = Σ F/(2πτ) through the flow
  in the past (τ >= 1), ln Π = −λφ, λ = ε/φ_ref, φ_ref = 44.15. Slope of the width in W1 [12, 28] and W2 [28, 44]
  (tolerance of boundness s_free/4 ≈ 0.010):
  - free: +0.042 ± 0.005 and +0.040 ± 0.006 (grows);
  - ε = 0.01: +0.037 and +0.027 (grows, a little slower);
  - **ε = 0.1: −0.004 ± 0.005 (unresolved, |s| + 2 SE = 0.014 just above the tolerance) and −0.036 ± 0.005
    (collapse)**; the width stops at ~3.35 and then falls to 2.7;
  - ε = 0.5: −0.042 ± 0.009 (collapse, R 3.0 → 2.1 until t ≈ 23), then +0.278 ± 0.017 (rebound and expansion to 7.6).

  The frame check passed in all branches. It is the first mechanism in the laboratory in which the bundle's own field
  overcomes the diffusion of discreteness and pulls the bundle together. The course corresponds to cold dust without
  pressure: under gravity it collapses (ε = 0.1) and without damping it flies apart after the collapse (ε = 0.5). An
  equilibrium (bound, stationary) state needs something that resists compression: pressure, a velocity dispersion, or
  the finite capacity (saturation) of the theory. That is a new hypothesis after a failure, with its own criteria.
  Caveat in advance: in 2+1 the potential is logarithmic, a physical reading requires 3+1.
- **The bending of a path is a property of the substrate**: the shortest chains bend too, and so does the average over
  all maximal chains, which is moreover about 5× less noisy.
- **Shielding by saturation** (above): without feedback the saturation does not spread beyond the driven region.
- **The heavy tail of the flow is a property of the topology of links, not of degrees**: a graph with the same degrees
  and swapped links has a narrower distribution (99th percentile 4.45 against 5.06) and hardly saturates.
- **The link graph of a causal set is rigid**: swaps of links that preserve a valid order succeed in only 0–0.8 % of
  attempts.
- **The response of the flow to an impulse decays exponentially** at a rate of about 1 per unit of time, almost
  independently of N.
- **Full bending of light requires c_eff = c·Π² near matter**: the weight Π alone gives half the bending, the conformal
  density zero.
- **The sum of Π lies above the continuum prediction, the sum of Π² on it, in 2+1 and in 3+1** – not yet explained. In
  2+1 by about 10 %; in 3+1 (F2.5 d4) by 8 % at N0 = 32 000 and by 22 % at 64 000 (difference −0.0058 ± 0.0037, not
  significant), while H2 fits to 3 %. The same pattern in two dimensions points rather to a property of the reading of
  Σ Π (small differences of the weights, coincidences of chains) than to chance. F2.5 d3 (2+1, 768 seeds per size,
  RECORDED) showed that it does not vanish with N: excess H1 − P1 = +0.0044 ± 0.0016, +0.0029 ± 0.0015,
  +0.0044 ± 0.0013 at N0 = 8 000, 16 000, 32 000 (+10 %, +7 %, +10 %; trend e(8 000) − e(32 000) = 0.0000 ± 0.0021 →
  PERSISTENT), while H2 − P2 = +0.0004 to +0.0011 (below 1.5 %). Rule B is SPECIFIC at every size. The excess is
  therefore systematic and independent of N in this range. The theory claims that maximal chains of Σ Π converge to the
  geodesic of Π²η (the continuum limit of the action); for Σ Π² the data confirm it, for Σ Π a shift of about 10 %
  remains that does not fall with N. Either the convergence of Σ Π is very slow, or the limit of Σ Π is not exactly the
  geodesic of Π²η; undecided.
- **Motion follows the geodesic of the field that the causal set created itself; that field is, however, steeper than
  the continuum (a)** (F2.7, chain, 2+1, 512 seeds, UNDECIDED). Source σ = 8, q = 0.01; the self-consistent field (a)
  on the causal set converged (13 iterations, residual 3.7e−5), Π(0) = 0.7113 against the continuum 0.7099. Motion
  (b = 10, T = 40): H1 = +2.027 ± 0.071 against P_a = +1.677 (geodesic of the continuum solution of (a)), ratio 1.21;
  |H1 − P_a| ± 2 SE = 0.207 to 0.493 against the band d_alt/4 = 0.419 → undecided according to the criterion (but it
  is neither the null field nor □_g without −R/2 with P = 3.49). Reported outside the criterion: the geodesic in the
  field of the causal set P_cs = +2.023, so H1 − P_cs = +0.004 ± 0.071 – the motion lies exactly on the geodesic of the
  field in which the chains were weighted. The whole difference against P_a arises in the source → field step: Π_cs
  along the path r = 10…30 is 0.757 → 0.883, the continuum (a) 0.746 → 0.856, i.e. the rise is ~15 % steeper
  (P_cs/P_a = 1.21). H2 = +3.740 ± 0.068 against P2_a = 3.322 (ratio 1.13, reported; P2 in the field of the causal set
  not computed). Tension with F2.5 d3: there Σ Π lay ~10 % above the geodesic of the given field, here on it (a
  different setting; undecided whether it is the same phenomenon). K.11c calibrated in 2+1 the distant flux balance of
  the operator, not the shape of the self-consistent profile at r = 10–30.
  Analysis of the saved profile (22 Oct, without a new run): in 2+1 (a) has the form ∇²Π = (κ/3)Π²ρ_s, so outside the
  source the flux r·dΠ/dr must be constant (Gauss's law; the continuum solution 0.1197 from r ≈ 30 up to the wall). The
  field of the causal set does not have this property: r·dΠ/dr peaks at 0.129 near r ≈ 22 and then falls towards the
  wall (0.124 at r = 36, 0.111 at 56, 0.073 at 79, 0.006 at 99); Π_cs reaches 0.9993 already at r = 89 (continuum
  0.986). The ratio of the slopes cs/a falls smoothly from 1.28 (r = 4) through 1.24 (10), 1.13 (20), 1.06 (30) to 1.0
  near r ≈ 43 and below 1 further on. The depth of the well agrees at the same time (1 − Π(0): 0.2887 against 0.2901).
  It is therefore not a different constant coupling (that would scale the slope the same everywhere), but a profile
  that is not harmonic outside the source and "feels" the wall R = 100: the remaining candidates are an edge effect of
  the rows at the wall (a cut-off past in a curved geometry, which the vacuum subtraction E[B_flat]1 does not cancel
  exactly) or screening by the non-local operator; a scan of R distinguishes them.
- **A bistable saturated phase** (F3.3f capacity, an observation to a FAILED result): a copied pulse with a threshold
  and with the capacity feedback C·Π creates the first state in the F3 series that persists on its own for several
  layers after the source is switched off (half the saturation after ~1 layer of 0.25, a trace still after 1.25). It is
  extended, not localised.
- **The noise of the operator is spatially correlated** (F1.2b, an observation, not a test): the error of the averaged
  field falls with the exponent p ≈ 0.5–0.8, far below d/2 for uncorrelated noise (1.5 in 2+1, 2 in 3+1). An exact
  reading of the field therefore needs a large L/ℓ, i.e. a denser sprinkling, not more rows of the operator.

## Derived consequences (no test)

- **The swerves of a defect are the shot noise of the order sampling its own binding: D_p = (1/2 rho) int V'(phi)^2
  phi_x^2 dx = (32/15) mu^(5/2) / rho for the sine-Gordon kink, with no free constant** (derived 10 Nov at the leading
  order; docs/derivations/swerves_of_a_defect.md).
  - **Why:** in 1+1 the retarded Green function is exactly 1/2 inside the cone, so the rule of F3.14a is exactly the
    continuum field with point sources (1/rho) sum V'(phi_y) delta(z - y). The difference from the continuous source
    is Poisson shot noise, white in spacetime, of variance V'^2 / rho. The momentum balance dP/dt = int s phi_x dx then
    gives the diffusion.
  - **Compared with the wandering study (development):** D_p rho measured
    - 0.06-0.10 for M 4 (prediction 0.067);
    - 1.5-3.2 for M 8 across rho 200-3200 (2.13);
    - 46-49 for M 16 (68; the narrow kink is coarsely resolved).
    - Agreement within the precision (~25-30 %) for M 4 and 8.
  - **Content:**
    - the vacuum does not swerve (V' = 0);
    - kicks come from internal structure, int V'^2 (d phi)^2, not from mass alone;
    - the noise is Lorentz invariant in distribution, so the diffusion is set in the rest frame (as in the swerves of
      Dowker, Henson and Sorkin), consistent with the smaller wandering of the moving kinks of F3.14a;
    - the formula gives the resolving power of future tests in advance.
  - **Limits:** 1+1 only (in 3+1 the Green function sits on links); leading order (no feedback on the profile); no
    number for real particles.

- **The fall rule of F5.1c/d/e is a proxy by the criteria of the axiom audit (9 Nov, the author's agreement)**
  (classification of a rule, no test). The rule chooses the next step from the events H = 8 links ahead of the
  history, i.e. from its future, against M4 ("the state of an event depends only on its causal past") and criterion
  (b) of the audit. Its result depends on the region in which the order is read, against M5 ("finite and well defined
  without a wall of the region or a regulator") and criterion (d): the horizon set is unbounded along the light cone,
  and the diagnostics F5.1-diag2 and F5.1-diag3 showed that the history takes the velocity of the cylinder and follows
  the shape of the admissible region. Its last tie-breaker and the step of chain_between choose the earliest event,
  i.e. by a linear extension (M4). The rule was introduced after the audit (F5.1c, 29 Oct) without being checked
  against criteria (a)–(d); this is a design error of the laboratory. The verdicts stay: F5.1d PASSED and F5.1e FAILED
  by their criteria (R132); what changes is what may be read from them as physics. The fall in 2+1 shows that a rule
  of this family can carry the gravity of the order into a moving history; it does not show that matter falls by a
  rule of the theory. A rule that satisfies M3–M5 uses only the past of the history (Axiom of the action: Δτ between
  successive events; Axiom of momentum: the memory of the flow) and finite intervals; F5.1-diag3c, already running,
  is evaluated by its frozen criteria as a property of the proxy.

- **Carrying by the flow alone does not give Newton; it leads back to the geodesic with inertia** (on paper, 29 Oct, on
  the author's picture of tilted cones). In the river picture of a static field (Painlevé–Gullstrand, space flows
  towards the mass at a speed v(r)) the rule "keep the velocity relative to the local flow" gives the acceleration
  (ẋ·∇)v, i.e. zero for a body released from rest; Newton ẍ = ∇(v²/2) is the geodesic equation. The frame of a
  "balanced surrounding" is, by the equivalence principle, every freely falling frame, whatever its velocity: the
  surrounding does not select a velocity by itself, it needs a memory (inertia). For a static source the tilting and
  the narrowing of the cones are the same order in other coordinates; physically new is the tilting only with rotation
  (Lense–Thirring), hence F5.6. From this follows F5.1c: a history continuing as the longest chain of the order.
- **The substrate rule of F7.4a is not a separate postulate: at linear order it is the response of an order with the
  Einstein–Hilbert (BD) stiffness to the action of flows A6** (derived 28 Oct, sympy;
  `docs/derivations/rule_from_action.md`).
  - The own frame and the addition follow from A6 and from linearity: δS_m/δh_ab = ½T_ab, because the proper time of a
    flow depends on h only through u_a u_b.
  - Lorentz invariance leaves a single number γ = −c_ΦΨ/(2c_ΨΨ). Einstein–Hilbert gives c_ΦΦ = 0, c_ΦΨ = −2(d − 1) and
    c_ΨΨ = (d − 1)(d − 2), hence γ = 1/(d − 2).
  - The throughput principle is a consequence in 3+1: dust has no stress, so h̄_ij = 0.
  - It is new only in the continuum (a textbook chain). A discrete proof through a sprinkling is out of reach (K.17
    FAILED, noise); F6.8 therefore computes the mean value of the BD action deterministically (F6.8a conformal sector,
    F6.8b mixed).
- **In 3+1 the throughput principle at linear order is the same as the tensor coupling to T_ab (linearised general
  relativity). A scalar equation as a complete dynamical law is excluded** (derived 27 Oct, sympy;
  `docs/derivations/tensor_coupling.md`).
  - h̄_ab = U[(1 − γ)η_ab + 2(1 + γ)u_a u_b]: γ = −1 is a pure trace coupling (Nordström), γ = 1 a pure coupling to
    T_ab.
  - Throughput gives γ = 2/(d_s − 1), the tensor coupling γ = 1/(d_s − 2). They agree only for d_s = 3.
  - A scalar Π radiates quadrupolarly 1/6 of the power of general relativity (circular binary) and has no g_0i. Double
    pulsars and Gravity Probe B exclude it.
  - Fork for the non-linear order:
    - A: throughput is a plateau law and the completion is general relativity (bootstrap, Jacobson). The throat and the
      2PN predictions of the scalar equation then drop out.
    - B: throughput holds exactly and the theory deviates from 2PN onwards.
  - The author chose A (decision 13, 28 Oct).
- **The vacuum equation in 3+1 (H1, ln Π) has an exact first integral r²L′ = M e^L (1 − L) and the closed solution
  Ei(1 − L) = Ei(1) + eM/r. A point mass thus has no horizon, it has a photon sphere at 3.070 M, b_crit 5.2435 M
  (+0.91 % against general relativity) and a throat at 2.341 M** (derived 27 Oct, sympy and a numerical ODE agreeing to
  5·10⁻⁸; `docs/derivations/gravity_equations.md`).
  - 2PN: −g_tt = 1 − 2x + 2x² − x³ and g_rr = 1 + 2x + 2x² + x³ (x = M/r). General relativity has −3/2 x³ and 3/2 x²
    and 1/2 x³, so the first deviation is in g_rr at second order.
  - Below the throat the area of the spheres grows again and g_tt → 0 only as r → 0, at an infinite proper distance.
  - Limitation: a continuum limit without a discrete test. The source has no pressure, so it is not clear whether a
    body can expose the throat at all.
- **The throughput principle is contained in the ontology of *Physics of Causality*, it just had not been joined into
  a formula. Section 26.6, however, contradicts it** (analysis 27 Oct, `docs/derivations/why_gamma1.md`, section 7).
  The book: c is the maximal stable throughput without congestion (5.3), c_eff falls with impedance (5.3, 7.7), the
  Einstein limit has a uniform maximal throughput (26.1), gravity is a redistribution of the density of updates in a
  substrate with a conserved flow (23), rigid structures increase branching in their surroundings (27.5), Π falls with
  saturation (6.5), matter is Π → 0 (19.1). The author's formulation ("throughput according to saturation, the more
  saturated a level of the cone, the more materialised") as a formula: q = S · c_eff = const., S = saturation of the
  layer. In 3+1 with clocks Π: S = Π⁻², c_eff = Π², hence H1 and γ = 1. Contradiction: 26.6 claims that gravitational
  dilation does not change the causal order. For γ = 1 the order must change (a fixed order gives γ = −1). The
  adjustment of the text must be decided by the author.

- **A candidate for the principle of γ = 1: conservation of the causal throughput ρ_c · c = const. gives
  γ = 2/(d_s − 1), i.e. γ = 1 exactly in 3+1, and satisfies the H1 metric in 3+1 exactly at all orders** (derived
  27 Oct, sympy, no test; `docs/derivations/why_gamma1.md`). Intuition from traffic and continuity (q = ρv): where the
  substrate narrows the cones (c), the events densify (ρ_c per coordinate volume) so that the number of events that the
  light front passes per unit time and area stays the same. ln(ρ_c c) = 2φ − (d_s − 1)ψ. In 2+1 it gives γ = 2, where
  general relativity has no Newtonian limit. In 3+1 the clocks Π and ρ_c c = 1 together give exactly H1 (A·B = 1,
  c = Π², ρ_c = Π⁻²), Schwarzschild only at first order ((1 − M²/4r²)²). One physical condition thus replaces the two
  choices of the H1 ansatz. Caveats: a postulate, not a derivation; referred to the Killing time and isotropic
  coordinates (a covariant formulation is missing, a candidate is counting the area of the front in events, K.14, which
  would meet the focusing principle); it requires the creation of events near matter (axiom 7). Proposed test: F7.1b in
  2+1 with the prediction γ = 2.

- **In the known ensembles of orders with the BD action the causal structure is flat, matter can change only the
  density in them, not the cones; γ = 1 cannot come out of them** (literature search and derivation 27 Oct, full text
  of Cunningham and Surya, arXiv:1908.11647). In 3D they do not use 3-orders (their "cones" are cubes, not Lorentz
  cones), but the occupation of sites of a fixed lattice on a flat R × T² (n = 300, ε = 0.1, transition β_c ≈ 1.98,
  first order, poor thermalisation near the transition). The relation there is determined by the flat lattice. 2D
  orders, in turn, are exactly causal sets embeddable in a 2D Minkowski interval. In both ensembles the action or matter
  can therefore change only the number of elements in a region, i.e. the volume (conformal) factor Ω²η. That is the same
  class of geometry as A6 on a fixed order (γ = −1, derived 27 Oct). F6.4b in these ensembles would not measure γ, it
  would only show the conformal response again. For the cones to be able to respond to matter, the sample space must
  contain orders with a non-flat causal structure. General posets, however, are dominated by Kleitman–Rothschild
  structures, and a manifold-like ensemble with free cones in 3D has not yet been constructed by anyone (an open problem
  of the field).

- **The action of axiom 6 is compatible with H1 only if the order itself narrows the cones to c_eff = cΠ²; on a fixed
  order it gives γ = −1** (derived 27 Oct, sympy, no test). The action S = −ħω_B Σ Π Δτ gives in the limit the
  Lagrangian Π dτ_u, where dτ_u is the proper time of the metric g_u in which Δτ is measured. The effective metric of
  motion is therefore Π² g_u:
  - on a fixed sprinkling (g_u = η): Π√(dt² − dx²) = √(Π²dt² − Π²dx²), i.e. the conformal Π²η, PPN γ = −1, light does
    not bend. This agrees with what was measured: maximal chains of Σ Π near the geodesic of Π²η (F2.5), the bending
    rule of F3.8–F3.10;
  - H1 (−Π²dt² + Π⁻²dx²) comes out exactly for g_u = −dt² + Π⁻⁴dx², i.e. for an order whose cones have the speed
    c_eff = cΠ²: Π√(dt² − Π⁻⁴dx²) = √(Π²dt² − Π⁻²dx²).

  Axiom 6 is therefore not wrong but underdetermined. Its physical content depends on what Δτ is (the open decision A1,
  `docs/derivations/variational_audit.md`) and on the dynamics of the order (axiom 7). On a fixed order A6 and H1 give
  two different geometries, and only an order whose cones matter narrows can reconcile them. It is the same requirement
  as F7.1 and the focusing principle (`docs/derivations/spatial_part.md`). Consequence for the laboratory: γ cannot be
  obtained by adjusting the axioms on a fixed order, only in an ensemble or growth of orders that respond to matter
  (F6.4b).

- **The horizon route to a₀ (Deser–Levin, Jacobson, Milgrom) does not give a direct derivation of the observed a₀**
  (derived 24 Oct, revised by the author, `docs/derivations/a0_horizon.md`). Consistently following Jacobson (κ of an
  auxiliary horizon) the temperature cancels and general relativity arises; a modification requires additionally
  identifying κ with the physical g. The full Deser–Levin temperature as G_eff makes the deep MOND regime non-existent
  (no solution for g_N < cH). The difference of temperatures (Milgrom 1999) as G_eff gives
  g_N = cH_dS(√(1 + (g/cH_dS)²) − 1): deep MOND with a₀ = 2cH_dS (9–11× more than observed) and a constant radial term
  ~cH_dS ≈ 5–7·10⁻¹⁰ m/s² even at high acceleration, in strong tension with planetary residuals (Iorio: 10⁻¹⁵ Mars,
  10⁻¹⁴ m/s² Saturn; Carloni et al.: perihelia exclude a universal acceleration of the size of the Pioneer anomaly,
  radar limit without caveats 3·10⁻⁹). The factor 1/(2π) does not arise from the route. This does not refute Milgrom's
  modified inertia nor all horizon mechanisms. K.14-dS only as a separate geometric test.
- **The focusing principle (every light-like direction k focuses in proportion to the energy flux of matter T_kk)
  enforces γ = 1 and Newton's law in a weak static field, without a scalar field equation and without the H1 ansatz**
  (derived 23 Oct, sympy, `docs/derivations/focusing_gamma.py`). With independent potentials φ (clocks) and ψ (space),
  R_kk = ∇²φ + ∇²ψ + nⁱnʲ∂ᵢ∂ⱼ(ψ − φ); equal focusing in all directions for matter at rest gives ψ = φ and
  2∇²φ = 8πGρ. For moving matter it requires the tensor field Σ F u u (gravitomagnetism). The spatial part of the metric
  is in a causal set the shape of the cones, i.e. a deformation of the order around matter; on a fixed sprinkling it is
  missing. Consistency of the principle requires that matter conserves 4-momentum: the derived correction of the merge
  rule F = |P|. The momentum of the field is carried by a deformation of the order that propagates. Plan:
  `docs/derivations/spatial_part.md`.
- **The Jacobson condition (Clausius ⇔ R_kk = 8πG T_kk) at linear order requires that the field carries the
  4-velocity of the sources: H1 as a scalar Π in a fixed frame violates it for moving matter, H1 written with the own
  4-velocity of the source satisfies it exactly** (derived 23 Oct, sympy, `docs/derivations/jacobson_linear_check.py`;
  weak field h_ab = −2U(η_ab + 2 n_a n_b), dust at a velocity v = 0.6, four light-like directions: differences −8.2 /
  +32.9 / +2.4 / −6.8 against exactly 0). The covariant form is linearised general relativity for dust. For many
  sources the field is a tensor, h̄_ab ∝ Σ F u_a u_b G_ret, with gravitational waves with polarisations + and × and with
  an asymmetric deformation of the cones. The momentum of matter supplied by F3.6 is therefore a necessary condition of
  the Jacobson route (step 3) and of gravitational waves. It holds at linear order and in 3+1; in 2+1 general
  relativity has no local gravity. Details: `docs/derivations/jacobson_route.md`, section 6.
- **The status of the anchor a₀ = cH₀/2π is contradictory in the text of the theory and the compression mechanism
  stands on the excluded K₀** (analysis of F6.7 on paper, 23 Oct, `docs/derivations/a0_anchor.md`). The text gives a₀
  the role of a reporting convention (it does not change the dynamics) and of a physical threshold of the transition
  D_eff 3 → 2. Only the threshold can be tested: F6.7 then means computing a*_kernel = 1/(2π) in units of cH₀ from the
  dynamics. The BTFR v⁴ = GMa₀ comes out only if the dimension changes at a fixed acceleration (at a fixed length
  v² ∝ M comes out). The plateau D_eff ≈ 3 is given in the text by saturation at K₀, which the tests excluded;
  compression needs a new basis (F6.4, or the effective dimension of the flow). The variant a₀(z) ∝ H(z), to which the
  motivation of "causal depth" leads, is in tension with rotation curves at z ≈ 2 for every mechanism with a transition
  at a₀(z); the constant variant (H₀ or Λ) agrees with the data. The decision is the author's.
- **The perihelion selects f = ln Π: in the equation B_g f(Π) = κρ_s with H1 (γ = 1) the PPN parameter is
  β = (1 − f″(1))/2, so for the family (Π^k − 1)/k β = 1 − k/2** (derived 23 Oct after F1.7 N, sympy: the outer vacuum
  solution Π = 1 − U + cU², U = M/r, inserted into the radial equation to second order). Results: ln Π: c = 1/2, β = 1
  (Π = e^(−U) to second order, g₀₀ = 1 − 2U + 2U², as in general relativity); Π − 1: β = 1/2; (Π⁻¹ − 1)/(−1): β = 3/2;
  (Π² − 1)/2: β = 0; (Π³ − 1)/3: β = −1/2. The perihelion precession is proportional to (2 + 2γ − β)/3; measurements
  (planetary ephemerides) give |β − 1| of order 10⁻⁴, hence |k| ≲ 2·10⁻⁴ and f″(1) = −1. Consequences:
  - the linear control Π − 1 (and with it the withdrawn proposal B_g Π − Π·B_g 1, which agrees with it at linear order)
    would give a perihelion precession 1/6 larger than general relativity (β = 1/2) and is excluded by observation;
  - ln Π is no longer just a motive (composition of clock ratios): it is the only member of the family compatible with
    the perihelion. Observation, however, fixes only f″(1); another f with f(1) = 0, f′(1) = 1, f″(1) = −1 gives the
    same at the first two orders;
  - ln Π in N does not have the widest range of solutions (k = −1 has a wider one), but k = −1 has β = 3/2.
  It holds under the assumption of H1 (γ = 1 as an input) and of the continuum limit of the operator; a numerical check
  from the outer profiles of N (fit 1 − U + cU²) is cheap and not yet done.
  On 23 Oct the author adopted ln Π as the canonical scalarisation (docs/author_decisions.md, item 1, fourth round).
- **The bending rule of the F3 series (F3.8–F3.10) is motion along geodesics of the conformal geometry Π²η (PPN
  γ = −1), so it does not bend light** (derived 23 Oct after F3.10 universality). A chain bends its carried direction
  towards −∇ ln Π at a rate |∇ ln Π| per unit proper time. For a chain with velocity v perpendicular to the gradient of
  a static field: the kick of one link Δu_x = |∇ ln Π| τ, the laboratory time of the link γτ, so the acceleration is
  |∇ ln Π| / γ² = |∇ ln Π| (1 − v²). This is exactly what the Euler–Lagrange equation gives for L = −Π √(1 − v²), i.e.
  the geodesics of the metric Π²η. General relativity gives (1 + v²), for light twice the Newtonian bending; Π²η gives
  zero for light (the null geodesics of a conformal geometry are flat). Consequences:
  - the universality of free fall (F3.10 universality) is automatic for a geodesic rule; an independent test of the
    equivalence principle requires bodies with a different composition of mass (binding energy, the Nordtvedt effect),
    i.e. bound structures, which we do not have yet (F3.7–F3.9);
  - the dynamics of F3 realises exactly the geometry that the consequence of equation (a) given above excludes by
    observation (Cassini γ = 1). Bending of light requires also the spatial part of the metric (H1: optical time
    L − (1 + γ)εU), which the rule of F3.8 does not use. K.13a does not decide it: it verifies that the discrete BD
    operator on the H1 geometry (γ = 1) converges to the continuum limit; only with a verified instrument will it be
    decided which field equation is compatible with γ = 1 (the contradiction of H1 and (a) in 3+1 above). K.13a PASSED
    23 Oct: the obstacle to the bending of light is not the substrate but equation (a).
  Cheaply verifiable: a fast probe perpendicular to the field, bending as a function of v (rule 1 − v², GR 1 + v²).
- **The field equation (a) and the hypothesis H1 cannot both hold in 3+1** (derived 19 Oct from the local limit
  □_g − R/2, which K.11 confirmed numerically only for the conformal geometry Π²η; for other geometries it relies on
  Belenchia et al. 2016). For the static metric −Π²dt² + Π^(−2γ)dx² (γ is the PPN parameter; Π²η has γ = −1, H1 has
  γ = 1) the symbolic computation gives the local form (□_g − R/2)Π:
  2+1: γ = −1: 3Π⁻²∇²Π; γ = 0: 2∇²Π + Π⁻¹|∇Π|²; γ = 1: Π²∇²Π + 2Π|∇Π|²;
  3+1: γ = −1: 4Π⁻²∇²Π + 2Π⁻³|∇Π|²; γ = 0: 2∇²Π + Π⁻¹|∇Π|²; γ = 1: 2Π|∇Π|² (no term with ∇²Π).
  The linear coupling is κ/c(γ), c(γ) = 2 − γ(d_s − 1), d_s the number of spatial dimensions (D = d_s + 1). In 3+1
  therefore (a) in the geometry Π²η (coupling κ/D, in 2+1 preliminarily F1.5) gives γ = −1, i.e. no bending of light and
  no Shapiro delay (Cassini: γ = 1); in the geometry H1 it gives the equation 2Π|∇Π|² = κρ_s, so outside the source Π
  is constant and no Newtonian field arises. In 2+1 H1 and (a) do not exclude each other (coupling κ). Consequence: F7.1
  in the sense "H1 with equation (a)" has no solution in 3+1; either the form of the field equation changes (what the
  operator acts on, or whether it carries −R/2), or H1. The decision is the author's. The derivation of the coupling
  κ/D → κ = 16π (chapter 4.7) holds only for the conformal geometry.
  Confirmed by the laboratory's own instrument (K.13a, 23 Oct): in the H1 geometry in 3+1 the coefficient of ∇²u in
  B_g 1 goes to −1, so in B_g Π = B_0[εu] + ε δB_u[1] the linear terms ε∇²u − ε∇²u cancel (c(1) = 0). The contradiction
  therefore no longer rests only on Belenchia et al. 2016; it holds at first order and at σ ≥ 4 ℓ_k with a deviation
  below 2.4 % at σ = 32.
- **The response of the BD operator does not distinguish the argument f of the equation B_g f(Π) = ρ_s** (derived
  22 Oct, leads to the cancellation of the planned K.13b). The operator is linear in the field. For f(1) = 0,
  E_g[f(Π)] = ε f′(1) B_0[u] + O(ε²), so the geometry (H1 or another) does not show at linear order. The difference
  ln Π − (Π − 1) at second order is −½ ε² B_0[u²], i.e. the flat operator on u² (∇²(u²) only in the continuum); the
  geometry enters only at ε³. The candidates f are distinguished only by a principle
  (docs/derivations/ln_pi_potential.md: additivity of clock ratios, c² ln N is the exact acceleration potential of
  static observers) and by the whole non-linear equation (F1.7 N).
- **The first order of the operator in the family −Π²dt² + Π^(−2γ)dx² is linear in γ** (sympy C6, C7, 22 Oct):
  δ(□_g − R/2)w = 2γ u∇²w + (1 − γ(d_s − 2)) ∇u·∇w + (1 − γ(d_s − 1)) w∇²u = (1 − γd_s) P + (1 + γ) Q.
  P = (−2u∇²w + (d_s − 1)∇u·∇w + d_s w∇²u)/D is the response to a change of the density √−g with flat cones (the
  conformal instrument K.11). Q = (2u∇²w + 2∇u·∇w + w∇²u)/D is the response to a change of the speed of the cones
  Π^(1+γ). The discrete instrument K.13a computes both parts separately.

- **The text of the theory and the laboratory stand on two different substrates; the tests have in effect made the
  choice** (axiom audit 22 Oct, `docs/axiom_audit.md`, corrected after the author's remark). *Physics of Causality*
  assumes a finite degree K₀, locality and dimension from graph distance, time as a number of updates and the Lorentz
  factor from clocks at rest relative to the substrate (a rest frame). The laboratory uses a Lorentz-invariant
  sprinkling with infinitely many links per event. By the theorem of Bombelli, Henson and Sorkin a finite degree and the
  invariance of the substrate exclude each other. The tests have already decided in favour of the invariant substrate:
  K₀ does not hold (F0.3), capacity limits the flow, not the links (F1.8), local growth of the graph type did not create
  a manifold (F6.1, F6.2). What remains unspoken: the text of the theory carries remnants of a substrate with a rest
  frame, and so do the rules of the laboratory (narrowed cones in coordinates, splitting of the flow according to the
  number of descendants; cf. F0.6).
- **The resistance functional R = Σ (I/Π) Δτ is not yet a functional of the field and the variational principle alone
  selects nothing** (audit 24 Oct, `docs/derivations/variational_audit.md`). For a fixed graph and I independent of Π,
  R is strictly monotonic in u = ln Π (δR = −Σ I e^(−u) η Δτ ≠ 0), so it has no interior stationary point. If
  Δτ = Π Δt, Π drops out of R entirely (R = Σ I Δt). Open decisions: A1 (what Δτ is), A2 (an operational definition of
  I). Every one-component equation A∇²L + B|∇L|² = f is variational after the integrating factor μ′/μ = (2B − A′)/A, so
  the existence of an action selects neither β nor ln Π. The static second variation of equation (a) is positive
  (∫ |∇η|² + (2κ/3) ρ Π η²) and says nothing about the stability of a structure. A test of the hypothesis D_eff ≈ 3
  belongs to phase C (stability of the causal set against a discrete action of the Benincasa–Dowker type, a protocol
  given in advance). The literature connections (Kleitman–Rothschild, phase behaviour) are an unverified search.

- **β in the discrete field: at second order it is determined only by two already calibrated parts and the source
  convention does not affect it** (derived 25 Oct, no test; the author left the decision on F1.9 to the laboratory). For
  L = ln Π = εL₁ + ε²L₂ the equation B_g L = κρ_s gives B₀[L₁] = κρ_s and B₀[L₂] = −dB_{L₁}[L₁] + (a source term
  ∝ ρ_s). β is read outside the source, where ρ_s = 0. Any convention of the source at second order (proper or
  coordinate density, the factor √−g = e^{(1−3γ)L}) therefore changes only the active mass (a harmonic response 1/r),
  not β. Outside the source ∇²L₁ = 0 holds; of the first-order response δ(□_g − R/2)w there remains
  (1 − γ(d_s − 2))∇u·∇w = (1 − γ)|∇L₁|² (d_s = 3). For S = c|∇L₁|², L₂ = −(c/2)L₁² and β = 1 − c/2, so for H1 (γ = 1)
  c = 0 and β = 1. The candidates f = (Π^k − 1)/k add k|∇L₁|² (β = 1 − k/2, in agreement with F1.7). The discrete β at
  second order thus depends only on the first-order response in H1 (K.13a: converges to the continuum form, 0.8–2.4 % at
  σ = 32) and on the flat operator (K.4, F1.3). **Not verified:** K.13a measured with Gaussian fields at the source
  (r ≤ 2σ), not with a 1/r field outside the source, where β is read. The non-local tails fall quickly there in 3+1, but
  a direct measurement of c(r) is missing. A possible test (F1.9a): c(r) = dB_{L₁}[L₁]/|∇L₁|² in rows outside the
  source, prediction 0, control with the same quadrature γ = −1 with prediction 2. Source convention: the laboratory
  keeps using ρ_s as a fixed coordinate function without √−g (F1.7, K.12). For the definition of mass (κ ↔ G) it must be
  fixed separately.

## Undecided (descriptive)

- **B1 development rows (2+1, intrinsic instrument, variant B; 12 Nov 2026; decide nothing): the domain wall keeps
  its velocity, contraction and charge; the inertial mass from a kick is 0.948 ± 0.028 of the rest energy (not resolved
  from 1).** `diag/b1_dev.py`, ρ 30, T 8, 8 sprinklings, arms rest / v 0.3 / v 0.6 / kick (f 0.05, t ∈ [1, 3)) on the
  same sprinkling; `diag/out/b1_dev_20261112T023039.json`. Slab rule t ∈ [0, T − 1) (the last slab at v 0.6 leaves the
  pyramid; set geometrically for all arms, recorded in the card draft, revision 1). Mean ± SE (sd per sprinkling):
  - velocity: v 0.3 +0.0019 ± 0.0058 (sd 0.0165); v 0.6 −0.0019 ± 0.0050 (sd 0.0141); rest drift +0.016 ± 0.059;
  - contraction C − 1: −0.006 ± 0.005 and −0.002 ± 0.006; width against the exact readout +1.6 % (rest), +1.0 %, +1.4 %;
  - charge Q − Q_ex: within 2e-4 in all arms;
  - energy per unit length (noise-debiased): rest +1.2 % ± 0.3 %, v 0.3 +1.3 % ± 0.4 %, v 0.6 +0.2 % ± 0.9 %; the plain
    estimator +7 / +7 / +4 %; the debiased energy per slab rises with time (rest 8.01 → 8.2–8.5, sd 0.04 → 0.4):
    fluctuation energy in the readout window;
  - kick: Δ(γv) 0.0816 ± 0.0025 against the continuum 0.0785 (+4 %, 1.2 SE); r_M = 0.948 ± 0.028, per-sprinkling sd
    0.085 (with the slope window [4, 8) instead of [4, 7): 0.961);
  - cost: N ≈ 80 000, 455 s per run per core.
  Reading: in development the transverse dimension and the intrinsic weight leave the inertial quantities of the wall
  within a few percent of the continuum; whether r_M deviates is not resolved. They fix the proposed tolerances and the
  sample size of B1 (32 sprinklings; binding condition r_M) in docs/B1_card_draft.md, revision 1. Nothing is decided.

- **D2 development rows (2+1): the gradient of the field from a compact past neighbourhood is accurate to 1-8 %, but an
  intrinsic neighbourhood (the links) is not compact and fails; and the pointwise T_ab of one scalar field does not fix
  a unique rest frame where the gradient is spacelike** (development, DESCRIPTIVE, decides nothing; decision 16 D2;
  diag/d2_gradient.py; run in the container 11 Nov evening).
  - **Convention (development, to be fixed before freezing):** signature (+, -, -), L = 1/2 (d phi)^2 - V,
    T_ab = d_a phi d_b phi - g_ab L, trace T = -1/2 (d phi)^2 + 3V in 2+1 (the same convention gives T = 2V in 1+1, the
    load of F5.2b).
  - **Set-up:** exact fields evaluated on the events of a sprinkling (so that only the estimator is tested): a sine-Gordon
    wall moving at 0.3, a Klein-Gordon plane wave k = (1, 0.5) of amplitude 0.5, a moving Gaussian lump; 300 targets;
    fit of phi(u) - phi(v) on the displacement with linear and quadratic terms, from the causal past only.
  - **E1, coordinate neighbourhood 0 < t_v - t_u <= h (development proxy):** at h = 0.5 the rms error of the gradient is
    0.4-1.4 % (wall), 4-8 % (wave), 2.5-6 % (lump) of |grad phi|; the error of (d phi)^2 2-25 %, of the trace 0-9 %. It
    grows with h (truncation of the one-sided quadratic fit) and does not shrink with rho at fixed h (the fields are
    exact; there is no noise to average). With an evolved field the shot noise will favour a larger h: a bias-variance
    choice still to be measured.
  - **E2, the past links of v and of its parents (an intrinsic neighbour set):** errors 25-125 % of |grad phi|, bias up
    to 1.8 in T_xx. The links run along the light cone to the edge of the box (median 1 000 / 3 400 neighbours at rho
    50 / 200): the set is not compact, so a local fit fails.
  - **Structural observations (on paper, checked against the formulas):**
    - A compact neighbourhood of an event needs a timelike direction (a frame): Lorentz invariant neighbour sets of a
      causal set are not compact. The gradient needs a frame, and the frame of D3 is to come from T_ab, which needs the
      gradient: a circle to be broken.
    - For one scalar field, T^a_b = d^a phi d_b phi - delta^a_b L has the eigenvector grad phi and a degenerate
      eigenspace orthogonal to it. Where grad phi is timelike the rest frame is unique (u ~ grad phi); where it is
      spacelike (every wall and kink, moving or not) every timelike u orthogonal to grad phi is an eigenvector: in 2+1
      the boost along the wall is free. This is physical for a wall (it has no velocity along itself), but it means that
      the rest frame of a defect cannot be read pointwise; it needs T_ab summed over a region (the momentum of the
      defect), or another source of the frame.
    - Where u ~ grad phi (timelike gradient), the rapidity of a link e relative to u is relational:
      cosh(eta) = |delta phi along the link| / (tau_link |grad phi|), with |grad phi|^2 = (d phi)^2 a scalar. The scalar
      (d phi)^2 = 1/2 box(phi^2) + phi V'(phi) needs only box(phi^2), because box phi = -V'(phi) is the field rule itself.
  - **The load without a gradient (decision 19; diag/d2_load.py):** in the growth rule the load enters only through its
    Green integral, and box(phi^2) = -2 phi J + K (K = 2 (d phi)^2) gives
    **G*K = phi^2 - chi_0 + 2 G*(phi J)** (chi_0 the free wave of the data of phi^2). So G*T = -G*K/4 + 3 G*V needs phi
    at the event, the free wave of the data and two more Green sums on the same chains: no gradient, no frame.
    - Checked on the intrinsic field (weight B), T 4, 4 sprinklings, 400 targets near the centre line; ratio of means
      (SE over sprinklings), rms per target relative to the continuum:
      - wall at rest: identity / continuum 0.997 +- 0.007 (rho 25), 0.998 +- 0.002 (rho 50); rms 0.053 / 0.037;
      - wall at 0.3: 1.001 +- 0.010 / 0.996 +- 0.004; rms 0.057 / 0.040;
      - the same discrete sum of the EXACT K (same weights): 0.998-1.009, but rms 0.08-0.10: the identity is about half
        as noisy as summing the kinetic term itself;
      - G*T by the identity route / continuum: walls 1.004-1.005 (rho 25), 1.001-1.009 (rho 50);
      - plane wave k = (1, 0), amplitude 0.1: identity / continuum 1.070 +- 0.013 (rho 25), 1.031 +- 0.008 (rho 50);
        G*V of the evolved field is off by the same amount (1.059 / 1.029) while the direct sum of the exact K is
        0.997 / 0.991. The excess belongs to the evolved wave (its amplitude), not to the identity, and halves when rho
        doubles.
    - Reading (development): the kinetic part of the 2+1 load is available from the order without derivatives, to
      better than 1 % on walls at rho 50. The remaining error of waves is the error of the field itself.
  - **The rest frame of matter (decision 19; diag/d3_frame.py; exact fields on pyramid sprinklings, rho 50 / 150,
    25 targets per cell; error = relative rapidity eta between the estimated and the true four-velocity; window of
    chain lengths [k, k + 2) with k = 1.3 rho^(1/3) tau, tau = 1, 2, 3):**
    - **F1a, world tube seen from the event** (core events, V above a threshold, in the past at chain distance k; mean
      direction): on the core 0.02-0.10; **off the core it points from the source to the event**: 0.25-0.8 at a
      rest-frame distance of 1-2. Unusable off the core.
    - **F1b, world tube from pairs of core events** (fully intrinsic in the selection): 0.04-0.37 (typically 0.1-0.2),
      no improvement with rho or tau (the chords inside a tube of finite width tilt).
    - **F2r, summed T_ab with the gradient of the E1 fit** (h = 3 rho^(-1/3)) over the loaded events of the past with
      L <= k: lumps 0.002-0.04 on and off the core (moving at rapidity 0, 0.31, 0.55); a second lump at rest 4.5 away
      leaves the frame local (0.005-0.027); walls 0.05-0.10 (the boost along the wall is free: an expected, physical
      degeneracy).
    - **Caveats:** (1) all test fields move rigidly. Then u is an exact eigenvector of the pointwise T_ab, so any sum
      keeps it: with the exact gradient F2 is exact and tests nothing, and with the E1 gradient only the gradient noise
      is tested. Non-rigid configurations (an evolved field, radiation, two bodies with different velocities that
      overlap) are still needed. (2) F2r uses the coordinate gradient of the development stage: the circle
      gradient <-> frame is not broken by it.
    - Reading (development): coarse-grained T_ab is the best of the three routes on rigid configurations; the fully
      intrinsic route through the world tube (F1b) reaches only ~0.1-0.2 in rapidity.
- **B0 development rows (2+1): with the Green function 1/(2 pi tau) (coordinate tau, a declared development proxy) the
  field rule carries a domain wall as the 1+1 kink and reproduces the 2+1 dispersion including k_y != 0; the longest
  chain estimates proper time with alpha_3 ~ 1.2, the same at a given realised density whatever the thinning, but it
  underestimates the Green weight near the cone by 20-40 %** (development, DESCRIPTIVE, decides nothing; decision 16;
  diag/b0_green2p1.py, diag/b0_chain_calib.py; run in the container 11 Nov).
  - **Field (pyramid |x| <= T + 8 - t, |y| <= T + 3 - t, T 6, 8 sprinklings, rho 50 / 100, ~83 000 / 166 000 events,
    cap tau >= 0.68 rho^(-1/3)):**
    - wall at rest: velocity -0.003 / -0.001 (+- 0.004);
    - wall launched at 0.3: 0.293 +- 0.004 / 0.302 +- 0.003;
    - rms against the 1+1 kink at the end 0.117 / 0.095 (falls with rho);
    - plane waves omega / sqrt(k^2 + 1) = 0.993-0.996 (rho 50) and 0.996-0.997 (rho 100) for k = (1, 0), (0, 1), (1, 1),
      (2, 1), approaching 1;
    - vacuum: rms ratio 1.13-1.15 against the exact 1.00-1.04, relative error 0.26 / 0.20 at amplitude 0.01. Not
      clearly growing; a longer run is needed.
    - Cost: ~9 s / ~33 s per run (brute force over pairs).
  - **Longest chain** (L links between y < x, intervals sprinkled in their rest frame; 400 pairs per cell):
    - alpha = L / ((rho_real)^(1/3) tau), for s = rho_real^(1/3) tau >= 4: 1.19 / 1.20 / 1.23 at rho_real 25 / 50 / 100.
      It grows slowly with s (finite size); for s < 1 it is 2.0-2.5 (L = 1 for short intervals).
    - **Thinning:** at the same realised density, alpha is the same whatever the thinning: 1.202 / 1.211 / 1.207 (rho_real
      50 from rho 50 p 1, rho 100 p 0.5, rho 200 p 0.25); 1.188 / 1.187 (rho_real 25). The longest chain reads the
      realised density, not the candidate density.
    - The spread of tau_hat / tau for s >= 4 is 0.10-0.12.
    - **Near the cone (s < 2):** mean 1/tau_hat / mean 1/tau = 0.60-0.81. The Green weight from a single alpha is biased
      low by 20-40 %.
    - Interval count / (pi/12 rho_real tau^3) = 0.90-1.04 (diagnostic).

  Reading (development).
  - The 2+1 Green function works as a field evolution: walls behave as kinks, and the dispersion converges, also for
    waves with k_y != 0.
  - The longest chain is a sound intrinsic clock that follows the realised density.
  - A single alpha is not enough near the cone. Proposal: a weight table w(L) = E[1/(2 pi tau) | L], calibrated under
    the volume measure of the past (tau^2 d tau). [Superseded 11 Nov evening, decision 17: the calibration condition is
    E[w(L) | tau] = G(tau) with two moment conditions; the calibration is not "exact by construction", the accuracy of
    the Green integral is a consequence to be tested.]
  - Open:
    - the cost: L for all pairs in the past is out of reach at full size, and needs an algorithm;
    - the stability of the vacuum over longer times.
  - **Addendum, 11 Nov evening (diag/b0_chain_alg.py, container):**
    - **Algorithm for L (the author's plan):** dynamic programming from each source over the past links, in time order,
      L_s(s) = 0, L_s(v) = 1 + max over links u -> v of L_s(u). The links (a numba version of the shell method, lam 25)
      are identical to core.order.links.links_csr (7 832 events, 163 607 links). L_s(v) equals an independent longest
      chain inside the interval on 373 / 373 random pairs.
    - **Costs (rho 100, pyramid T, X = T + 8, Y = T + 3):** links 33-65 per event (they grow with the box: the shell
      reaches far along the cone); building the links 3.4 / 10 / 22 / 44 s at N 52 000 / 82 000 / 119 000 / 165 000; the
      DP 2.7 / 8 / 13 / 25 ms per source, so 2 / 11 / 26 / 68 min for all N sources on one core.
    - **The table depends on the source (confirmed):** the posterior mean E[1/(2 pi tau) | L] differs between a
      near-cone pair measure (~ tau d tau, which dominates in a finite box) and the bulk measure (~ tau^2 d tau): w(1) =
      0.195 / 0.138 rho^(1/3), w(2) = 0.084 / 0.077; from L = 4 within 3 %. So the posterior table is not valid for all
      sources.
    - **Calibration under decision 17 (diag/b0_weight_calib.py; s_max = 6, beta from the rms relative residual over s in
      [2, 20], both declared in the script before the run; 40 000 + 8 000 intervals, an independent sample of the same
      size):** a = 0.2733 +- 0.0038, c = 0.2131 +- 0.0005, beta = 1.203 +- 0.039 (bootstrap 200); equivalent tail alpha
      2 pi c = 1.339; rms relative tail residual 0.021.
      - Moments of r over [0, 6] / moments of G on the independent sample: first +0.0018 +- 0.0025, second -0.0001 +-
        0.0031 (zero within the calibration uncertainty).
      - **Links alone (L = 1):** their share of the error is +0.082 +- 0.002 (first moment) and +0.059 +- 0.001
        (second), compensated by L >= 2.
      - **Pointwise E[w~ | s] / G (independent sample):** 1.135 at s 0.5-1, **1.458 at 1-1.5**, 1.255 at 1.5-2, 0.99 at
        2-2.5, **0.93 at 2.5-3.5**, 0.95-0.99 at 3.5-6, 0.995-1.004 beyond 6. The local error is large; only the integral
        is unbiased.
    - **Validation on pyramid sprinklings (diag/b0_weight_valid.py; |x|, |y| <= T + 2 - t, T 4; 8 sprinklings per
      density; 150 targets each at t in [2.5, 3.5], |x|, |y| < 1; L by backward DP over the links; none of these
      configurations was used in the fit):**
      - **Density scaling and boost independence:** E[w~] / E[G] per bin of s is the same at rho 50 / 100 / 200 and in
        the three rapidity classes (eta < 0.5, 0.5-1.5, > 1.5) to about 0.01, and equals the calibration curve
        (e.g. 1.469 / 1.482 / 1.478 at s 1-1.5; 0.940 / 0.945 / 0.944 at s 3-4).
      - **Green integral, ratio of means sum I_hat / sum I** (coordinate tau with the cap | intrinsic A, alpha_3 1.2 |
        intrinsic B, the table):
        - rho 50: 0.92-0.93 | 0.99-1.01 | 0.978-0.995 (SE 0.005-0.01);
        - rho 100: 0.953-0.957 | 1.019-1.022 | 1.006-1.013 (SE 0.003-0.008);
        - rho 200: 0.952-0.965 | 0.999-1.018 | 0.996-1.006 (SE 0.002-0.004).
        The coordinate variant is biased by the cap (-7 / -5 / -4 %); A is about 1 % above B at every density.
      - **Dependence on the source (paired with the homogeneous source in the same sprinkling, B):** at rho 200 the walls
        of width 1 and 0.5 at 0 and 35 degrees -0.3 to -0.4 % (+- 0.13-0.31 %), the blob of width 1 -1.0 +- 0.3 %;
        at rho 100 all within 1.3 SE of zero; at rho 50 -0.8 to -1.7 %. The coordinate variant shows the same pattern
        (blob -1.3 +- 0.5 % at rho 200), A slightly more (blob -1.9 +- 0.3 %).
      - **Outside the declared class (diagnostic):** a wall of width 0.1, -1.9 +- 0.9 % at rho 200; a step, -0.3 +- 0.2 %.
      - rms per-target relative error (shot noise of the sum): 4-6 % at rho 200, similar for all three variants.
      - Cost: links 1.2 / 1.6 / 5.5 s, backward DP 1-5 ms per target (N 14 000 / 28 000 / 56 000).
    - **Long vacuum (diag/b0_vacuum_long.py; coordinate tau with the cap, rho 50, T 12, 8 random waves, 4 sprinklings
      per amplitude):** A / A_exact grows steadily to 1.36 / 1.39 / 1.36 at t = 12 for amplitudes 0.003 / 0.01 / 0.03;
      lambda_eff = +0.024 / +0.024 / +0.025 (first half) and +0.028 / +0.031 / +0.025 (second half) per unit time,
      SE 0.001-0.002. Exponential growth at a constant rate, independent of the amplitude (linear), not an initial
      transient and not a constant bias. Relative error grows about linearly to 0.50.
      - **Predictions written before the runs at rho 25 and 100 (amplitude 0.01):** multiplicative shot noise of the
        sum (heating) gives lambda ~ 1 / rho (0.05 at 25, 0.013 at 100); the near-cone cap bias gives lambda ~
        rho^(-1/3) (0.033 / 0.021); an instability of the rule itself gives lambda independent of rho.
      - **Result (amplitude 0.01, 4 sprinklings):** lambda_eff (first / second half) = 0.035 / 0.043 at rho 25, 0.024 /
        0.031 at rho 50, 0.020 / 0.022 at rho 100 (SE 0.0005-0.003). 1/rho is excluded (0.012 expected at 100). A power
        law gives exponent about 0.4-0.5; lambda_inf + b/rho (lambda_inf ~ 0.013) fits too. Three points do not
        separate them. lambda is larger in the second half at every density.
      - **Mechanism under test:** the cap removes weight in a thin shell at the cone (the coordinate variant's integral
        is -7 / -5 / -4 % at rho 50 / 100 / 200, scaling like the cap). A kernel concentrated at the cone has a
        Fourier transform that is complex for timelike (omega > k) modes, so its deficit can give the massive modes a
        growth rate. **Predictions written before the run (rho 50, same sprinklings and modes, cap factor 0.5 and 2):**
        if the cap drives the growth, lambda follows the cap (~0.012 at 0.5, ~0.05 at 2); if shot noise drives it, a
        smaller cap gives more noise and a larger lambda (the opposite direction).
      - **Result:** lambda_eff (first / second half) = 0.0137 / 0.0184 at cap factor 0.5, 0.024 / 0.031 at 1, 0.046 /
        0.057 at 2 (SE 0.002-0.003). lambda follows the cap almost proportionally, as predicted for the cap mechanism and
        opposite to shot noise. A linear extrapolation to cap 0 leaves 0.003 (first half) / 0.006 (second half): small,
        not separated from zero by three points.
      - **Reading (development):** the growth of the vacuum seen with coordinate tau comes mainly from the cap of the
        development proxy, not from the rule. It does not close the vacuum question: the intrinsic weight B is unbiased
        only in the integral; near the cone it has its own deficit (E[w~] / E[G] ~ 0.42 at s < 0.5) and an excess at
        s 0.5-2. The vacuum must be measured with the intrinsic variants (stage 3) before B0 is closed.
    - **Intrinsic field, stage 3 (diag/b0_chain_field.py; forward DP from every source over the links of its future,
      general source J, threshold eps_J = 0):** instrument checks: the forward kernel equals the backward DP of the
      validation to 3e-15 (30 targets); the field under the linear extension t + 0.3 x equals the time-ordered field to
      3e-15.
      - **rho 25, T 6, 4 sprinklings, the same sprinklings and initial data for the three variants** (coordinate with the
        cap | A, alpha_3 1.2 | B, the table; ~100 s per intrinsic run against 2 s for the coordinate one):
        - wall at rest: velocity -0.013 / -0.009 / -0.008 (+- 0.01); launched at 0.3: 0.300 / 0.293 / 0.298 (+- 0.011);
          rms against the 1+1 kink at the end 0.16-0.17 / 0.24-0.25 / 0.26 (the intrinsic variants distort the wall
          more);
        - plane waves, omega / omega_c: 0.993-0.996 / 1.009-1.035 / 1.004-1.007; amplitude last / first quarter (energy
          proxy): 1.10-1.19 / 1.12-1.16 / 0.96-0.98;
        - vacuum (8 waves, |x|, |y| <= T + 4 - t): lambda_eff +0.030 +- 0.002 / +0.036 +- 0.002 / +0.007 +- 0.002.
      - Reading (development): the table B removes most of the growth of the vacuum and of the waves, and gives the best
        dispersion; A (one alpha) is biased in omega by 1-3.5 % and grows like the coordinate variant. B leaves a small
        decay of single waves (2-4 % over T 6) and a small residual growth of the vacuum, and distorts the wall more than
        the coordinate variant.
      - **rho 50, T 6, 4 sprinklings (run by the author, 11 Nov evening, diag/out/b0_chain_field_20261111T191554.json;
        ~410 s per intrinsic run):**
        - wall at rest: velocity -0.001 / +0.003 / +0.003 (+- 0.002-0.004); launched at 0.3: 0.297 / 0.300 / 0.302
          (+- 0.004-0.006); rms against the kink at the end 0.107-0.113 / 0.193-0.195 / 0.203;
        - plane waves, omega / omega_c: 0.993-0.996 / 1.008-1.030 / **0.9996-1.0020**; amplitude last / first quarter:
          1.07-1.15 / 1.05-1.07 / **0.970-0.988**;
        - vacuum: lambda_eff +0.025 +- 0.003 / +0.020 +- 0.001 / **+0.0038 +- 0.0006**.
      - **From rho 25 to rho 50 (variant B):** omega / omega_c 1.004-1.007 -> 1.000-1.002 (converges to 1); the decay of
        single waves 2.4-3.8 % -> 1.2-3.0 %; lambda of the vacuum 0.0068 -> 0.0038 (about 1/rho^0.8, towards zero);
        rms of the wall 0.26 -> 0.20 (factor 1.29; the coordinate variant 0.165 -> 0.11, factor 1.5). A keeps a bias of
        omega of 1-3 % that does not shrink (one alpha is not enough). The velocity of the wall is right in all three.
      - Reading (development): B is the intrinsic instrument that converges in dispersion, amplitude and vacuum. Open:
        whether the larger rms of the wall under B is noise (the scatter of the chain weight) or a distortion of the
        profile; to be separated by the profile averaged over sprinklings against the scatter between them.
      - **Noise or distortion (diag/b0_wall_profile.py; wall at rest, rho 25, T 6, 4 sprinklings, last slab, |y| < 2,
        bins of 0.25 in x):** per-event rms in the core |x| < 2: coordinate 0.20, B 0.30; in the tails 0.13 / 0.22. The
        bias of the profile (mean over sprinklings) is 0.09 / 0.07 in the core and 0.03 / 0.04 in the tails, with an SE
        of ~0.03-0.08: not resolved for either variant. The scatter of the bin means between sprinklings (it includes the
        wandering of the wall) is 0.13 / 0.16. So the larger rms of the wall under B is mostly the noise of single events
        (the scatter of the chain weight), not a resolved distortion of the profile. The coordinate variant shows a bias
        of up to +0.15 at the centre (its cap).

- **Load candidates (1+1): 2V(phi_v)/rho, the trace of the energy-momentum tensor per event, summed over the world tube
  of a kink gives M times its proper time in every frame and in a nonuniform order; the BD-based kinetic estimators are
  unusable per event** (development, DESCRIPTIVE, decides nothing; docs/derivations/field_and_throughput.md, section 5;
  decision 15; diag/load_candidates.py; 1+1, rho 400, eps 0.1, 8 sprinklings per arm, slab t in [4, 6), |x - x_kink| <
  6; run in the container 10 Nov, 557 s).
  - **Target:** the sum over the slab of the per-event contribution equals -M dtau (M = 8 sqrt(mu)).
  - **V2 = -2V(phi_v)/rho** (sum / (-dtau), discrete / the same estimator on the exact field):
    - rest: 8.02 +- 0.06 / 8.00;
    - v 0.3: 8.09 +- 0.04 / 8.03;
    - v 0.6: 7.86 +- 0.06 / 7.95;
    - mu 4 (M 16): 16.21 +- 0.18 / 15.96;
    - curved order of F5.2a2 (target M with Omega at the centre): 7.84 +- 0.07.
    - Zero in the vacuum (|.| <= 7e-5 per unit area). The centroid follows the kink (-0.01, 1.497, 3.021 against 0,
      1.5, 3.0).
    - For a plane wave (amplitude 0.3, k 1): -0.0446 per area, equal to the continuum trace m^2 A^2 / 2 on the exact
      field. V2 measures the trace of T, not the energy and not the Lagrangian.
  - **The BD-based kinetic estimators are unusable at this density.**
    - grad (difference form): 4.4 +- 4.7 (rest), 3.3 +- 4.3 (v 0.3), -0.9 +- 7.9 (v 0.6), the same on the exact field.
    - phiB (not shift invariant): wrong and noisy (vacuum per area -9.9 +- 17.6).
    - The per-event fluctuation of the BD sum (probe P-a) dominates.

  - **Independent check by the author** (new sprinklings, --base 1, 16 per arm, 10 Nov). V2, sum / (-dtau):
    - rest 7.96 +- 0.04; v 0.3 8.06 +- 0.07; v 0.6 7.88 +- 0.09; M 16: 16.03 +- 0.25; curved 7.69 +- 0.11;
    - vacuum <= 1e-4 per area; plane wave -0.0443 (exact -0.0446).
    - grad and phiB are again unusable (rest: 1.0 +- 4.3; vacuum per area of phiB -6.5 +- 13.6).
  - **Note on the curved target:** the target in the curved arm used Omega at the centre but not the velocity of the
    falling kink. At t ~ 5 the kink moves at v ~ g t ~ 0.25, so dtau carries an extra sqrt(1 - v^2) ~ 0.97 and the
    target is ~7.75. Both runs (7.84, 7.69) agree with it.
  Reading (development). In 1+1 the kinetic term of a scalar field is conformal, so the trace of its energy-momentum
  tensor is exactly 2V(phi). The local, frame-free, phase-periodic quantity 2V(phi_v)/rho is therefore the trace per
  event. On a soliton (dust: no net stress) it is the rest-mass density, and its sum over the world tube is M dtau, the
  action of a flow of the Axiom of the action of a flow (S = -sum m dtau).

  It is consistent with F5.2a2: in 1+1 the counting couples to matter only through Omega^2 V'.

  Limits:
  - the identity is special to 1+1; in 3+1 the trace has a kinetic part, and a discrete gradient estimator better than
    the BD sums is needed;
  - the trace is the scalar source (enough for matter at rest); moving matter and pressure need the tensor (the
    throughput principle in 3+1 is the tensor coupling).

- **H5-0 development rows (1+1): the retarded Green function of the order (1/2 x the causal relation) carries a
  sine-Gordon kink with its charge, velocity, contraction and energy, and its inertial mass approaches its rest energy
  as the density grows: M_in = 8.62 / 8.34 / 8.13 against E = 8.09 / 8.06 / 8.03 at rho = 100 / 200 / 400; the BD
  operator in retarded form diverges** (development, DESCRIPTIVE, decides nothing; card docs/H5_card_draft.md rev. 3;
  1+1 Poisson sprinkling in the domain of dependence |x| <= X - t of the initial slice, X 26 (kinks) / 20 (waves), T 10,
  16 sprinklings per density, `diag/h5_instrument.py`, Fenwick engine; run in the container 10 Nov, 66 s).
  - **Rule:** phi(x) = phi_0(x) - (1/2 rho) sum over y < x (the order), t_y >= 0, of [sin phi(y) + f chi(t_y)];
    phi_0 the free wave of the initial data (d'Alembert, the launch); V = 1 - cos phi.
  - **BD in retarded form (probe P-a):** the sum over the past has mean -1.0 (normalisation right) but per-event sd
    0.5-1.1 with sum |K| ~ 7-11; a constant field diverges within ~2 units of time. Recorded as the failure of the
    differential form as an evolution; the BD operator stays an instrument for other tests (decision 14).
  - **I4:** two linear extensions (t, t + 0.5x), two engines (Fenwick, brute-force blocks) and two regions (the same
    events) agree to <= 1.4e-15 relative.
  - **I1 vacuum (amplitude 0.01, k in (0.25, 4]):** rms last / first 1.041 / 0.988 / 0.971 against the exact 1.033 /
    0.982 / 0.967; relative error in the last slab 0.150 / 0.107 / 0.074 (falls as rho^(-1/2)). No growth.
  - **I2 dispersion:** omega / sqrt(k^2 + 1) within 0.001 at every k (0.5, 1, 2, 4) and density; at rho = 400 the
    mass shift omega - k is reproduced to 1.000 / 1.000 / 0.999 / 0.999.
  - **I3 and H5-A (rest):** rms error against the exact kink 0.121 / 0.116 / 0.046; drift of the centre over 9 units
    +0.11 +- 0.09 / +0.07 +- 0.09 / -0.03 +- 0.04 (the common drift of probe P-b was noise); Q 0.996 (readout on the
    plateaus 5-7 from the centre: 0.997 on the exact kink); w 0.999; E 8.093 / 8.058 / 8.026 (the same readout on the
    exact field 8.008-8.011).
  - **H5-B (free motion):** v = 0.3 gives 0.318 +- 0.014 / 0.302 +- 0.009 / 0.302 +- 0.003; v = 0.6 gives 0.597 +-
    0.007 / 0.597 +- 0.006 / 0.599 +- 0.003; w gamma 0.995-1.004; E / gamma 7.98-8.13 (readout on exact 7.98-8.01).
  - **H5-C (kick f = 0.05 for 1 <= t < 3, paired with the rest kink):** Delta(gamma v) 0.0726 +- 0.0042 / 0.0750 +-
    0.0036 / 0.0770 +- 0.0027 against 0.0785 for M = 8, i.e. M_in = 2 pi Q f tau / Delta(gamma v) = 8.62 / 8.34 / 8.13;
    M_in / E = 1.065 / 1.035 / 1.013 (SE ~ 0.035-0.06): the excess roughly halves with each doubling of the density.

  Reading (development): the order's own retarded Green function is a working conservative dynamics in 1+1. It is
  stable, it reproduces the dispersion, it converges, and it carries a topological defect with the continuum
  properties. Two limits apply. In a Lorentz-invariant sprinkling with a covariant rule, a moving kink is in
  distribution a boosted resting kink (the lesson of F3.13c): once H5-A holds, H5-B is largely fixed by symmetry. M_in
  = E follows from the momentum balance of the field, which holds on average when the rule is invariant in
  distribution. What is not fixed by symmetry: that the discrete sum is stable at all (BD is not), the size of the
  discreteness corrections (M_in / E - 1 ~ rho^(-1)?), and the separation from a relaxational dynamics (H5-E, not yet
  run). The frozen numbers of H5-0 and H5-A to H5-C are proposed from these rows (card, section 3).

- **F5.1-diag3c: the proxy rule without the cylinder keeps most of the velocity of its history, but slows down
  gradually even in the whole box; UNDECIDED by the frozen criteria** (diagnostic of the proxy of F5.1c/d/e,
  DESCRIPTIVE; flat 3+1, box t [0, 26], |x| <= 11, |y|, |z| <= 9.5, 8 new sprinklings of ~206 500 events, 16 distinct
  starts per sprinkling and velocity (256 of 256 distinct), starts centred so that the walk is symmetric in the box;
  engine 16 of 16; run by the author 9 Nov, 11 721 s with 8 processes).
  - **Frozen primary Delta_h = v_walk − v_hist** (slope over the whole walk minus the velocity of the history handed
    to the rule; SE over sprinklings, over walks in brackets): A0 (cylinder R 4) v 0 +0.009 ± 0.020, **v 0.3 −0.243 ±
    0.036** (instrument OK); Binf v 0 −0.001 ± 0.019 (0.025), **v 0.3 −0.022 ± 0.027** (0.020). KEPT needs |Delta_h| +
    2 SE < 0.05 for both v: v 0 0.039 yes, v 0.3 0.076 no; LOST needs Delta_h + 2 SE < −0.05: +0.032 no. Verdict
    UNDECIDED. Paired Binf − A0 −0.010 ± 0.023 (v 0) and **+0.221 ± 0.036** (v 0.3).
  - **Reported, not evaluated:** in Binf at v 0.3 the walk decelerates steadily: quadratic fit a = −0.0119 ± 0.0033,
    first three steps +0.262, late half +0.216 against the history +0.278 (late − history −0.061 ± 0.030); at v 0 a =
    +0.0025 ± 0.0030. The launch chains were nearly straight this time (first half − whole −0.031 ± 0.025 at v 0.3,
    against −0.091 in F5.1-diag3b, whose starts lay off the middle of the box); the explanation of the bend in
    F5.1-diag3b as a property of chain_between alone is therefore not confirmed. Shared events between walks 0.30–0.57.

  Reading. Without the cylinder the proxy keeps the velocity to within its precision on average over the walk, but it
  loses ~0.06 by the late half even in the whole box: what remains of the boundaries of the laboratory (the box and the
  time cap t_v + 9, both at rest) still brakes it weakly. This agrees with its classification as a proxy (M5); the
  rule for F5.1f is built from the past only (F5.1-diag4).

- **F5.1-diag4 (2+1): the set of admissible continuations of a history is finite without a wall: about 8 children of
  the last event, within ~3 units of time and space, always with no deficit against the memory** (diagnostic,
  DESCRIPTIVE, decides nothing; flat 2+1, density 1, 4 sprinklings of ~35 000 events (box t [0, 34], |x|, |y| <= 16;
  repeated with t [0, 40] and larger windows), 320 histories (first half of a longest chain, v 0 and 0.3), memory
  m = 8; run in the container 9 Nov).
  - **The question (the author's reading of superposition):** only p (m steps back) and v (the last event) are fixed;
    the next layer is the set of admissible continuations, the children c of v with L(p, c) = L(p, v) + 1 (the
    history through v stays a longest chain from p), or the smallest deficit if none reaches it. A rule that shares
    the flow among them is finite (M5) only if this set is finite.
  - **Result:** the children of v in the window grow with the window without limit (3.2, 8.4, 17.9, 30.8, 43.4, 56,
    75 for dt_w = 1.5, 2.5, 4, 6, 8, 10, 13), the admissible set does not: 3.2, 6.0, 7.7, 7.9, 8.0 and then 8.7 / 7.6
    (v 0 / 0.3) unchanged from dt_w = 8 to 13. The set equals the set of the largest window in 0.89–0.97 of the
    histories at dt_w = 4–6 and in 1.00 from dt_w = 8; it does not depend on the box (12 or 16: 1.00). The smallest
    deficit is 0 in every history. The set reaches ~3.1 in time and ~3.0 in space from v, and contains children up to a
    rapidity ~2.4 relative to the history.
  - **Reading:** min L(p, c) over the children tames the light cone: the set is finite without a wall, as the
    continuum suggests (deficit ~ distance along the cone), but its size (~8) and spread (rapidity up to 2.4) are
    large. A rule that shares the flow equally among them is finite; its velocity noise per step will be large, and
    whether the centre of the shared flow keeps the velocity is the question of F5.1f-a. 3+1 is to be run by the
    author (diag/band.py --dim 3).

- **F5.1-diag3b: without the cylinder the walk keeps the velocity of the history it is given; the frozen reference
  (the whole launch chain) was itself bent by the cylinder, so the verdict is UNDECIDED** (diagnostic, DESCRIPTIVE,
  UNDECIDED by the criteria frozen before the run; flat 3+1, the box, order, rule and engine of F5.1-diag3, 8 new
  sprinklings of ~206 400 events, 16 distinct starts each (128 of 128 distinct), v 0 and 0.3; engine 16 of 16; run by
  the author 9 Nov, 10 029 s with 8 processes).
  - **Frozen primary Delta = v_walk − v_ref** (least-squares slope over s ∈ [0, 12] minus the velocity of the whole
    launch chain; SE over sprinklings, over walks in brackets): A0 v 0 −0.016 ± 0.017, **A0 v 0.3 −0.285 ± 0.011**
    (instrument OK); Binf v 0 +0.034 ± 0.022 (0.023), **Binf v 0.3 −0.099 ± 0.017** (0.025). KEPT needs
    |Delta| + 2 SE < 0.1 for both v: v 0 0.078 yes, v 0.3 0.133 no; LOST needs Delta + 2 SE < −0.1: −0.065 no.
    Verdict UNDECIDED. Paired Binf − A0 +0.050 ± 0.031 (v 0) and **+0.186 ± 0.020** (v 0.3).
  - **Reported, not evaluated (found after the run):** the history handed to the rule (the first half of the launch
    chain) moves at v_hist = +0.204 ± 0.017, the whole chain at v_ref = +0.295 ± 0.003. Relative to the history it was
    given, the walk in Binf keeps its velocity: v_walk − v_hist = −0.009 ± 0.016 (v 0.3), acceleration −0.003 ±
    0.005; in the cylinder −0.194 ± 0.021. The bend is not the cut: on a small flat order (40 chains, 9 Nov) the first
    half minus the whole chain was −0.222 ± 0.037 both with the cylinder of chain_between at radius 4 and at radius 7;
    it belongs to chain_between itself (it steps along the tube of maximal chains to the earliest child). The
    history at rest in Binf accelerates along the throw, a = +0.018 ± 0.005 (v_walk − v_hist +0.052 ± 0.019), towards
    the middle of the box (its start lay at x = −+2.5), as in C2 of F5.1-diag3. Per-walk SD of the slope ~0.25 (late
    slope ~0.35). Shared events between walks of one arm 0.21–0.30.

  Reading. The braking of F5.1e / F5.1-diag2 is the cylinder: without it the walk continues with the velocity of its
  history to within ±0.03 (2 SE). This comparison was not frozen and decides nothing (R132); the frozen reference was
  wrong because the launch chain built by chain_between is bent. Repeated relative to the history handed to the rule,
  with centred starts and a tolerance of 0.05, as F5.1-diag3c.

- **F5.1-diag3: without the cylinder the thrown history keeps most of its velocity (Binf − A0 = +0.24 ± 0.08), but
  the frozen inertia criterion is not met because of the statistics; where the admissible region is asymmetric about
  the history, the target follows it** (diagnostic, DESCRIPTIVE, UNDECIDED by the criteria frozen before the run; flat
  3+1, density 1, box t [0, 26], |x| ≤ 11, |y|, |z| ≤ 9.5, 8 sprinklings of ~206 400 events, 8 starts each, v 0 and
  0.3, the registered rule m = H = 8 computed by a new engine without bitsets (future cone of p, one bit per tied
  target), which reproduces walk_history exactly on the cylinder R = 4 in 16 of 16 checks; time cap t_v + 9; run by
  the author 8 Nov, 14 205 s with 8 processes).
  - **Late velocity (s ∈ [6, 12]), late − history, SE over sprinklings (in brackets SE over the 64 walks):**
    A0 (cylinder R 4 at rest): v 0 +0.011, −0.014 ± 0.096; v 0.3 −0.021, **−0.288 ± 0.043** (0.033).
    Binf (whole box, only diamonds inside): v 0 +0.044, +0.018 ± 0.084 (0.047); v 0.3 +0.219, **−0.048 ± 0.087**
    (0.044). Bplain (whole box, all candidates): v 0 +0.061, +0.035 ± 0.064; v 0.3 +0.209, −0.058 ± 0.070.
    C2 (diamonds ≥ 2 from the faces): v 0 +0.226, +0.175 ± 0.068; v 0.3 +0.295, +0.029 ± 0.036.
  - **Paired differences of v_late:** Binf − A0 +0.033 ± 0.067 (v 0) and **+0.240 ± 0.077** (v 0.3; walk SE 0.050);
    Binf − Bplain −0.017 ± 0.046 and +0.010 ± 0.082; C2 − Binf **+0.144 ± 0.060** (v 0) and +0.044 ± 0.096 (v 0.3).
  - **Readings frozen before the run:** inertia restored in Binf (|late − hist| + 2 SE < 0.1): not met for v 0
    (0.19) nor v 0.3 (0.22), with the walk SE 0.11 and 0.14; the point estimates are compatible with keeping the
    velocity and with losing up to ~0.2. Admissibility matters (|Binf − Bplain| > 2 SE): no. The boundary still reaches
    the decision (|C2 − Binf| > 2 SE): yes at v 0, no at v 0.3. The alternative reading "Binf does not keep the
    velocity, so the loss is in the rule" cannot be drawn either.
  - **Reported, not evaluated:** the events between p and the tied targets, n(p, q), are ~1 210 in the cylinder and
    ~2 110 in the box: the cylinder removed ~43 % of every diamond (the cut of F5.1-diag2 seen directly). The ties of
    min L(p, q) have 67–73 candidates in A0, 133–137 in Binf, 200–213 in Bplain, 50–54 in C2; their excess
    L(p, q) − L(p, v) − L(v, q) is 0.02–0.11 (aligned candidates exist), the fewest events leave one candidate in
    99–100 % of steps. Binf removes 43–50 % of the horizon set, C2 84–85 %. In C2 the history at rest drifts along the
    throw direction for both signs of the throw (+0.33 for the starts at x = −2.5, +0.13 for x = +2.5), i.e. towards
    the middle of the box, where the admissible region is larger. Per-walk SD of the late velocity ~0.35.

  Reading. The braking of F5.1e / F5.1-diag2 is mostly the cylinder: in the box the thrown history loses 0.05 ± 0.09
  of 0.27 instead of 0.29 ± 0.04. Whether the rule keeps the velocity to within 0.1 is undecided at this precision.
  The same sensitivity as in F5.1-diag2 appears in another form: the choice follows the shape of the readable region
  wherever it is asymmetric about the history (C2). A computation of the rule therefore needs a region much larger
  than the diamonds and a flat arm with the same geometry, so that the drift of the region is subtracted.
  - **Design error found after the run:** the start was the event nearest to a point with a jitter of 0.5, so starts
    of the same sprinkling collapsed onto one event: in Binf v 0.3, 22 of 64 walks repeat another walk exactly. The
    effective number of walks is smaller and the SE over walks is too optimistic; the SE over sprinklings stands.
    Repeated with distinct starts as F5.1-diag3b.

- **Fall diagnostic (diag/fall_layers, 7 Nov): the tie-breakers do not decide the fall in either dimension; in 3+1
  the history behaves like a body in a medium at rest in the laboratory (a kick at the apex, then a terminal velocity),
  and the suspect is the computation cylinder fixed at the apex** (diagnostic, DESCRIPTIVE, decides nothing; 2+1 with
  the parameters of F5.1d, 3 + 3 sprinklings, 12 histories each, run in the container; 3+1 with the parameters of
  F5.1e, 8 + 8 sprinklings of ~105 000 / ~98 000 events, 60 histories each, run by the author; the variant "earliest"
  reproduces the registered rule exactly in both).
  - **Layer 1, geometry (tube of longest chains between fixed events):** 3+1 loaded sag +0.432 ± 0.084 (rule scale,
    dT 14.8) and +0.386 ± 0.083 (control, dT 16), flat −0.033 ± 0.085 and −0.117 ± 0.087; loaded − flat 0.47 ± 0.12
    and 0.50 ± 0.12 against the geometry 0.61 and 0.71. 2+1 loaded +0.59 ± 0.30 and +0.89 ± 0.28 against 0.53 and
    0.68, flat −0.25 ± 0.27 and −0.05 ± 0.35. The order carries the gravity in both dimensions.
  - **Layer 2, the decision in a step:** level-1 ties (minimal L(p, q)) have 63 / 71 candidates in 3+1 (loaded /
    flat) and 21 / 43 in 2+1, in 95–100 % of steps; the fewest events n(p, q) leave more than one candidate in 1 % of
    steps in 3+1 and in 3–27 % in 2+1. The last tie-breaker shifts the radial velocity of the target by +0.056 ± 0.015
    (3+1 loaded), +0.020 ± 0.007 (3+1 flat), +0.008 ± 0.005 (2+1 loaded).
  - **Layer 3, the fall:** R = a / a_pred for earliest / latest / random: 3+1 loaded −0.14 / −0.21 / −0.19 (± 0.13),
    flat +0.32 / +0.28 / +0.29 (± 0.20); 2+1 loaded 1.10 / 1.07 / 1.14 (± 0.38), flat 0.59 / 0.46 / 0.60 (± 0.36).
    The tie-breaker changes R by less than its error in every sprinkling. The hypothesis written in F5.1e (the choice
    is decided by the tie-breakers rather than by the geometry) does not hold.
  - **Shape in 3+1 (earliest, loaded − flat displacement from the apex at τ = 0, 2.2, 4.4, 6.0, 8.2, 10.4, 12):**
    0, −0.27, −0.53, −0.62, −0.73, −0.85, −0.90 against the prediction ½ a_pred τ² = 0, −0.05, −0.21, −0.40, −0.74,
    −1.18, −1.58. The velocity jump at the apex is −0.168 ± 0.019 loaded and −0.018 ± 0.035 flat. The flat history
    moves out by 0.21 in the first ~3–4 units of time and then stays (0.21–0.25). Flat throws (48 walks): v 0 history
    −0.073 → walk −0.004 ± 0.051; v 0.3 history +0.176 → walk +0.061 ± 0.051, acceleration along −0.011 ± 0.006. In
    2+1 the flat throws instead gained velocity (+0.22 → +0.41 ± 0.10, then −0.021 ± 0.005).
  - **The "wide" arm is not readable:** the cylinder ×1.5 (3+1) and ×2 (2+1) reaches beyond the box (|x| ≤ 8 / 18),
    so it measures the cylinder and the edge of the box together (3+1 loaded +0.14 ± 0.18, flat +0.50 ± 0.21).

  Reading (a hypothesis, not tested). In 3+1 the geometry turns the history inwards at the apex, and then its velocity
  relaxes in ~4 units of time: in the flat order to rest, in the loaded order to a terminal inward velocity, hence the
  linear displacement of F5.1e. This is how a body moves in a medium at rest in the laboratory. A candidate for that
  frame: the cylinder of LocalLinks is fixed at the apex, with radius 4 in 3+1, narrower than the reach of the horizon
  (H = 8 links ≈ 7.4 units of time); it removes more of the candidates that continue the velocity and so pulls the
  target back towards its axis. In 2+1 the cylinder (12) is wider than the reach (~9.6) and the rule fell. Not tested:
  the velocity kept in a flat order as a function of the cylinder radius, with a box large enough. F5.1e stays FAILED
  (R132).

- **F3.13c: the moving hedgehog keeps its charge and moves with its source; its rest-frame cost is 0.99 of the cost at
  rest at β = 0.3 and 0.87 at β = 0.6, but the difference measures the finite box, not the defect** (mechanism, 3+1,
  box t [0, 20], x [−7, 19], y, z [−7, 7], ~101 900 events, 16 sprinklings, UNDECIDED 31 Oct by the criteria fixed
  before the run; instrument OK: vacuum 0 on 240 of 240, ≥ 50 readable shells in every defect arm, frame check 4 of 4).
  - **Charge and motion:** +1 on 162 / 147 / 70 of 240 shells (all readable ones) at β = 0 / 0.3 / 0.6; centroid
    velocity in the source frame +0.0036 ± 0.0029 / −0.0011 ± 0.0034 / +0.0010 ± 0.0030 (limit β/4); no core event.
  - **Cost:** r(0.3) = 0.986 ± 0.029 (meets the PASSED condition), r(0.6) = 0.866 ± 0.025 (|r − 1| + 2 SE = 0.18 > 0.15,
    |r − 1| − 2 SE = 0.08 < 0.15: undecided). Reported: E_lab / (γ E(0)) = 1.036 ± 0.041 and 0.949 ± 0.031.
  - **Diagnosis after the run (decides nothing):** the arms had 1 065 / 906 / 768 source events. The source region
    (0 ≤ t_s < 2, r_s < 5 in the source frame) lies partly below the laboratory time t = 0 when it moves: the fraction
    inside the box is 1.000 / 0.859 / 0.732 (Monte Carlo), i.e. 1 047 / 899 / 766 expected sources, as observed. The
    moving hedgehogs started with 14 % and 27 % of their source ball cut off on the −x side, and the substrate below
    t = 0 was missing from their causal past. The slab 12 ≤ t_s < 14 at β = 0.6 was also cut by the top of the box
    (0.96 of the ball r_s < 6 inside).

  Reading. The criteria fixed before the run give UNDECIDED and stay. The design, however, could not answer its
  question: the sprinkling is Lorentz invariant, the sources and orientations were given in the source frame, and the
  dynamics is covariant (frame check: identical flows and costs after a boost of 0.6). In such a set-up an arm moving
  with β differs from the arm at rest only by the shape of the region in its own frame, here the laboratory box; any
  deviation of r from 1 is therefore a boundary effect, and here it follows the size of the cut. A defect that moves
  relative to nothing is a resting defect seen from another frame: that its rest-frame cost is the same is a property of
  the covariant rules, and the frame check is the test of it. This is my design error: the test asked a question whose
  answer was fixed by symmetry, and the only free element was the box. A revision with congruent regions would give
  r = 1 by construction and is not proposed. What the run does show: a moving flow carries the locked defect with its
  charge and keeps it with its source. The physical content of a rest mass (inertia, release of the cost) needs motion
  relative to something: another defect (annihilation of a hedgehog with an anti-hedgehog), a medium, or a push.

- **F8.1 (run 29 Oct): the numbers carry no information, because the input was the invalid linear fit of the optical
  sector from F6.8c** (prediction, RECORDED; instruments I1–I3 OK, but I1 accepted an input with the status FAILED,
  which is a design error). The coefficients c of the optical sector came from the linear fits of F6.8c (13.2 / 1.33 /
  0.50 / 0.46 / 0.45 for k = 0.05–0.2), so ρ = −1.50 ± 1.75 (sign undetermined); β₄ and the bounds on ℓ_k (1.4e−20 m
  for a cosmological T) as well as γ_loc = 0.28 (from K0_bb = 1.559) are not read and do not go into the paper. F8.1
  will be revised only after F6.8e (asymptotics of the optical age term); the revision will require inputs with the
  status PASSED.

- **F7.1b: a discrete order with conserved causal throughput (narrower cones and more events) gives γ ≈ 1.8–1.9 in
  2+1, close to the prediction 2 and far from 1 and from 1/2; formally UNDECIDED because of statistics** (mechanism,
  2+1, 16 sprinklings of ~20 000 events; instrument OK: flat control zero, the exact metrics γ = 2 and γ = 1 clearly
  distinguished: 1.97 / 1.84 against 0.86 / 0.83). Event density exp(σ), cones exp(−σ) by the capacity rule, σ0 = 0.3:
  - γ_A (cones + density) 1.813 ± 0.073, γ_B (clocks + cones) 1.944 ± 0.077; the routes agree (difference
    0.13 ± 0.11), the order is metric;
  - against the exact metric γ = 2 read with the same instrument: A −0.156 ± 0.084, B +0.108 ± 0.086. The criterion
    PASSED (|deviation| + 2 SE <= 0.25) is not met (0.32 and 0.28), nor is FAILED (|deviation| − 2 SE > 0.25);
  - the reading by bins follows the continuum expectations: density 0.28 / 0.30 / 0.24 (0.30 / 0.28 / 0.25), cones
    −0.30 / −0.29 / −0.26, clocks −0.106 / −0.085 / −0.084 (−0.099 / −0.093 / −0.082).

  σ0 = 0.45 (only reported): near the centre the order stops behaving like a metric. The cones read −1.06 instead of
  −0.45 and clocks + cones give γ 3.8 ± 0.4, cones + density 1.1. Reachability through the kept links in the dense,
  strongly narrowed core lags behind the nominal cone. The book describes this as a loss of fidelity of the geometric
  reconstruction outside the plateau (5.2, 26.6).

  Replication (seed 1, new sprinklings, the same criteria; UNDECIDED, instrument OK): γ_A 1.768 ± 0.029,
  γ_B 1.960 ± 0.080, the controls γ = 2 read 1.879 / 1.763 and γ = 1 read 0.929 / 0.822. Against the control γ = 2:
  A −0.111 ± 0.059 (|deviation| + 2 SE = 0.229 <= 0.25, this route passes), B +0.197 ± 0.103 (0.403, does not pass). In
  both runs the shift has the same sign: route A reads 0.11–0.16 lower than the exact metric, route B 0.11–0.20 higher.
  This points to a small systematic deviation of the discrete order (reachability through the kept links and the
  density near the centre), not to chance. Both routes together give γ ≈ 1.85 ± 0.1.

  Reading: the throughput principle gives in a discrete order a metric geometry with γ around 2 in 2+1, as
  γ = 2/(d_s − 1) predicts. Capacity without adding events gives 1/2 (F7.1a), a fixed order −1. A formal verdict needs a
  smaller error (replication), physics needs the 3+1 version with the prediction γ = 1.

- **F7.1c: a discrete order with conserved causal throughput gives γ = 1.000 ± 0.045 and 0.972 ± 0.074 in 3+1,
  exactly as predicted. Formally UNDECIDED, because the control it is compared with undercounted by chance in one
  route** (mechanism, 3+1, 16 sprinklings of ~54 000 events; instrument OK: flat control zero, the exact γ = 1 and 1/3
  distinguished). σ = σ0·exp(−r²/32), density exp(σ), cones exp(−σ):
  - σ0 = 0.2 (55 % of links kept): γ_A (cones + density) 1.000 ± 0.045, γ_B (clocks + cones) 0.972 ± 0.074; the routes
    agree (0.027 ± 0.087). The reading by shells fits the continuum expectations to thousandths: density 0.197 / 0.163 /
    0.133 (0.197 / 0.171 / 0.129), cones −0.193 / −0.169 / −0.126, clocks −0.095 / −0.085 / −0.061 (−0.098 / −0.085 /
    −0.064);
  - **narrower cones only** (σ0 = 0.2, without added events; reported): 0.349 ± 0.017 and 0.328 ± 0.032, prediction
    1/3;
  - σ0 = 0.3 (reported): 0.924 ± 0.027 and 0.901 ± 0.045, near the centre the density lags (0.258 against 0.295);
  - controls in the same run: the exact γ = 1 reads 0.984 ± 0.071 and **0.781 ± 0.078**, the exact γ = 1/3 reads 0.309
    and 0.266. The calibration K.16c read the same metric as 0.994 and 1.013, so route B of the control is ~2.4 SE low
    here.

  Why UNDECIDED: against the control γ = 1 route A is +0.016 ± 0.085 (|deviation| + 2 SE = 0.19, passes), route B
  +0.192 ± 0.107 (0.41, does not pass). By R132 the verdict holds. The comparison with the calibration K.16c
  (A 0.006 ± 0.067, B −0.041 ± 0.096) is a computation after the run, it does not change the verdict.

  Replication (seed 1, new sprinklings, the same criteria; UNDECIDED, instrument OK): γ_A 0.948 ± 0.047,
  γ_B 0.979 ± 0.058 (the routes agree, 0.031 ± 0.074); σ0 = 0.3 0.975 / 1.043; narrower cones only 0.318 / 0.310
  (prediction 1/3); controls: γ = 1 reads 0.965 / 0.900, γ = 1/3 0.290 / 0.217. Against the control γ = 1:
  A −0.017 ± 0.070 (0.157, passes), B +0.079 ± 0.107 (0.293, narrowly does not pass). In both runs route B therefore did
  not pass, and the reason is the error of the comparison (SE ~0.1, of which ~0.08–0.09 from the control), not a shift
  of the candidate. The criterion |deviation| + 2 SE <= 0.25 leaves only ~0.04 for a real deviation at 16 sprinklings:
  the test was statistically weak for a verdict. That is a fault of the design, not of the result.

  Reading: all three predictions of the principle in 3+1 fit the raw values. A fixed order gives −1 (control), the
  narrowing of the cones alone 1/3, conserved throughput 1, and the order is metric. A caveat remains: the test verifies
  that a discrete order built according to the principle (density = 1/c) has a geometry with γ = 1. Why the substrate
  conserves throughput it does not derive.

## Open questions

- Will the author explicitly adopt the Lorentz-invariant substrate to which the tests lead (F0.3, F1.8, F6.1, F6.2),
  and rewrite the axioms accordingly (locality, dimension, time, derivation of the Lorentz factor)? (axiom audit,
  section 0)

- Does one order in 3+1 distinguish the field from the null model at a higher N (larger L/ℓ)?
- Is the positive bias of A5′ at large scales (in both dimensions, not significant) chance or a systematic effect?
- Does a persistent localised structure arise with feedback without tuning (F3.3)?
- Does such a structure have inertia and a gravitational field (F4, F5)?
- Why do maximal chains of Σ Π lie about 10 % above the geodesic of Π²η, while Σ Π² lies exactly on the geodesic of
  Π⁴η, and why does this not change with N (F2.5 d3)? Is the continuum limit of the action W = Σ Π really the geodesic
  of Π²η?
- Why does the self-consistent field (a) on the causal set in 2+1 not have a constant flux r·dΠ/dr outside the source
  and flatten towards the wall R, while near the source it is steeper (F2.7)? Is it an edge effect of the instrument (it
  moves with the wall at a larger R), or non-local screening (it stays at the same r)?
- Is non-locality a property of nature or of the instrument (F0.4)?
- Will the narrowing of the future cones at the place of slow clocks (c_eff = cΠ² from F7) hold the saturated region
  together, when uniform splitting of the flow into the whole cone dissolves it at the speed of light (F3.3g A)?
- What is the dynamics of the order: a local rule, or a weighted sum over orders (F6)?
