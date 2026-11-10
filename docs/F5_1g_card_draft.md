# F5.1g — Inertial coherence of a distributed causal flow (card, 9 Nov 2026)

Status: **revision 2, approved by the author 9 Nov and implemented as phases/F5_gravity/F5_01g_coherence.py; the numbers of the launch and of the regions are frozen in its params (the author chose to run without 2+1 development rows).** The marks *(to be fixed)* below refer to those params: T_L 5, r_s 1.5,
S 10, box |x| <= 9, |y|, |z| <= 7, dt_w 6 (region B), box |x| <= 7.5, |y|, |z| <= 5.5, dt_A 4.5 (region A).

## Why

- The fall rule of F5.1c/d/e is a proxy (axiom audit, 9 Nov): it reads the future and depends on the wall (M4, M5).
- F5.1-diag4: the admissible continuations of a single history (children c of v with L(p, c) = L(p, v) + 1) form a
  finite set for a launch history (2+1: ~8, 3+1: ~21–25).
- F5.1f-a (INSTRUMENT_FAILED; descriptive): a single-event flow shared equally among them disperses (spread 3.2–3.8 in
  10 units, per-walker velocity SD 0.37) and its admissible set then depends on the region (M5 identity 0.10). One
  link carries no proper time; a single chain cannot hold a direction.
- F3.6 (PASSED): a flow that carries its timelike direction stays localized and moves with its source, but the
  direction is carried in coordinates (embedding), and such a direction does not fall (F5.1 FAILED).
- The author's decisions (9 Nov): (1) the free flow first, gravity later (F5.1h homogeneous gradient, F5.1i
  equivalence of different bundles, F5.1j Newtonian limit); (2) memory is **distributed**: every event reconstructs it
  from its local causal past, no global state of the bundle; (3) the memory of an event is the **common past of all
  the parents that bring it flow**, not the memory of a dominant parent; (4) that common past is taken **along the
  flow** (the past in the whole order has an infinite maximal boundary along the past light cones and would violate
  M5 by construction).

## Question

Does a flow whose events take their memory from the maximal common flow-past of their flow parents, and pass their
flow equally to the children with no deficit against every memory event, keep (a) the velocity of its centre and
(b) its coherence in a flat 3+1 order, (c) without depending on the region in which the order is read?

## Definitions (the rule; only the order, counts and the flow)

- Flow amounts Φ_pv ≥ 0 on links p → v; inflow of v: Φ_v = Σ_p Φ_pv + s_v (s_v: source, only during the launch).
- Flow parents P_Φ(v) = {p : Φ_pv > 0}; flow-past J⁻_Φ(v) = the events from which flow reached v along flow links
  (v included).
- Common flow-past C_Φ(v) = ⋂_{p ∈ P_Φ(v)} J⁻_Φ(p); memory M_v = ∂C_Φ(v) = its maximal elements (no y in C_Φ(v) with
  x ≺ y). With a single parent p, M_v = {p}. A source event without parents carries the memory of the source (below).
- Deficit of a child: D(c) = Σ_{x ∈ M_v} [L(x, c) − L(x, v) − 1] ≥ 0 (0: the history through v stays a longest chain
  from every memory event).
- Admissible children (review point 1): A(v) = {c : D(c) = 0}. **Zero deficit is a condition of continuing, not the
  best of the bad options:** if no child has D = 0, A(v) = ∅ and the flow stays in v (ΔQ_v = Φ_v, recorded).
- Passing: the inflow of v is shared equally among A(v); streams that meet in an event add (Axiom of continuity:
  Σ Φ_in + s_v = Σ Φ_out + ΔQ_v).
- The rule is deterministic (equal sharing, no random choice): the same sprinkling and source give the same flow.
- No coordinates, no time order, no orientation vector and no weights enter the rule.

## Launch (the "throw")

A source world tube of radius r_s *(to be fixed)* moving with velocity v along ±x for a launch time T_L *(to be
fixed)*; the tube is the initial condition (defined in the laboratory, as every launch), not part of the rule.
Initialisation (review point 4):
- every event of the tube receives s_v = 1 as inflow; the source term creates no parent;
- the memory of a tube event is built like any memory, but its "flow-past" includes the tube: C(v) = the common past
  of its real flow parents together with the tube events in its causal past; an event of the first layer of the tube
  (no tube event and no flow parent in its past) has an empty memory and keeps its inflow (ΔQ, it does not pass);
- the events above the first layer therefore start with the tube as their history, exactly as a thrown body starts
  with the history of the hand;
- after T_L the source terms are zero and the flow is free.
Velocities v = 0 and 0.3; flat order of density 1; computing region a box with a time window above the current
events *(sizes to be fixed after the development rows)*.

## Branches (same sprinklings, same sources)

| branch | memory | role |
|---|---|---|
| G0 | M_v = P_Φ(v) (the parents only, depth 1) | negative control: reproduces the dispersal of F5.1f |
| G1 | F3.6 carried direction (coordinates; one child of smallest rapidity) | positive control of the estimators (known inertia and coherence); not a candidate rule (M3) |
| G2 | M_v = ∂C_Φ(v) (this hypothesis) | the test |
| G3 | G2 with the memories permuted among the flow events of the same time slab | negative control: does causal memory matter or only the presence of some memory |

## Readouts

Flow-weighted, per time slab after the release, in the laboratory (only readouts, M2):
- centre: velocity of the flow-weighted centre along the throw, s ∈ [0, S] after the release; Δ = v_centre − v.
- width: flow-weighted RMS across and along the throw; growth rate over the second half of the walk.
- orientation coherence: spread of the step velocities of the flow (flow-weighted SD).
- diagnostics of the memory: age of the common flow-past (L from M_v to v), |M_v|, the number of distinct flow parents,
  η_common(v) = |C_Φ(v)| / |⋃_p J⁻_Φ(p)| (the share of the parents' flow-past that is common; a readout, not a
  parameter), share of the parents whose own last events lie within depth d of M_v, number of admissible children per
  event, **share of the events with D_min > 0** (the flow held), flow held and flow lost at the boundary (ΔQ).

## Criteria (to be frozen after the development rows; SE by bootstrap over sprinklings)

- G1 is reproduced first (review point 2): development rows with the same sources, readout times, centre and width
  as the test; if G1 does not pass its instrument condition in the full run, the card ends INSTRUMENT_FAILED, and G1
  is not adjusted to pass.
- INSTRUMENT_FAILED if G1 does not show inertia and coherence: |Δ(G1, v)| + 2 SE ≥ 0.05 for either v, or its width
  growth rate + 2 SE ≥ 0.2 (F3.6 carried spread at 0.16–0.20).
- finite (M5), in two parts (review point 3); G2 computed in two regions A ⊂ B on the same sprinklings and sources:
  - **M5a, interior identity:** an event is boundary-affected in a region if a child of it was dropped because its
    diamond left the region; the interior of the comparison is the set of events of both runs whose flow-past contains
    no boundary-affected event of either run. Because the rule is deterministic, the flow on the interior must be
    identical in A and B (to rounding, relative 1e-9) in every sprinkling;
  - **M5b, global statistics:** paired |Δ_A − Δ_B| + 2 SE < 0.05 and |w_A − w_B| / w_B + 2 SE < 0.1 for both v, with
    the flow lost at the boundary reported separately for A and B.
  The identity of single walkers of F5.1f is not a criterion here.
- inertia (G2): |Δ(G2, v)| + 2 SE < 0.05 for v = 0 and 0.3.
- coherence (G2): width growth rate of G2 + 2 SE < 0.2 and smaller than G0 by more than 2 SE.
- memory matters: width growth of G3 exceeds that of G2 by more than 2 SE.
- Verdicts (review point 7), checked in this order:
  - INSTRUMENT_FAILED: G1 fails its condition;
  - FAILED: M5 — M5a not identical, or M5b violated;
  - FAILED: inertia — Δ(G2, 0.3) + 2 SE < −0.05, or |Δ(G2, 0)| − 2 SE > 0.05;
  - FAILED: coherence — width growth of G2 − 2 SE ≥ 0.2, or G2 not more coherent than G0 by 2 SE;
  - FAILED: memory content — G3 not less coherent than G2 by 2 SE (a failure of this hypothesis about the causal
    content of the memory, not a proof that memory as a mechanism does not exist);
  - PASSED: all of finite, inertia, coherence and memory content;
  - UNDECIDED: the intervals do not decide.

## Reading declared now

- PASSED: a flow whose memory is its own common past, read by the order alone, keeps its velocity and holds together;
  it is the candidate for F5.1h (the same rule in a homogeneous gradient of the order: σ linear, density e^σ, cones
  e^−σ).
- FAILED, not finite: the common flow-past still lets the admissible set depend on the region.
- FAILED, inertia or coherence: the common flow-past is too old or too narrow to hold a direction (the diagnostics of
  the memory say which).
- G3 ≈ G2: what holds the flow together is the presence of a memory, not its causal content.

## Known risks (stated in advance)

- The common flow-past can be old when branches split early; then the memory is weak (measured).
- With a single parent the memory is one link, as in G0; coherence must come from the merging of streams.
- Cost: one forward relaxation per distinct memory event; the number of flow events grows with the spreading.
- With zero deficit as a condition, a large share of held flow (D_min > 0, F5.1f: ~40 % of the steps of single
  walkers) would stop the flow; this is measured and is itself an answer.
- Even if F5.1g passes, it gives a candidate mechanism of inertia, not gravity; gravity starts in F5.1h with the same
  rule unchanged.
- The launch tube sets the initial memory; its radius and duration are parameters of the initial condition, not of
  the rule, and are reported.
