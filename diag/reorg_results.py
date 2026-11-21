"""Reorganise results/ into phase folders (results/<phase folder>/<ID>/) with a README.md per phase.
Usage: python reorg_results.py <repository root> <docs/test_list.md> [--apply]
Without --apply it only prints the plan. Moves directories (rename); deletes nothing; never overwrites."""
import re
import shutil
import sys
from pathlib import Path

PHASES = [
    ("K", "K_instruments", "Instruments and calibration",
     "Instruments that read geometry and physics from a causal order, and their calibration against known answers: "
     "dimension, proper time from chains, cones and density (the blind instrument), the Benincasa-Dowker operator, "
     "horizon molecules, the energy-momentum tensor of flows, and interval volumes. A calibration has no physical "
     "conclusion of its own; it makes the physics tests possible, and a physics test that depends on a failed "
     "calibration cannot conclude."),
    ("F0", "F0_substrate", "Substrate",
     "What the substrate is: events and their order. Whether a substrate with a finite number of neighbours per event "
     "(and therefore a rest frame) is viable, or whether a Lorentz-invariant causal set (a Poisson sprinkling) is "
     "required; geometry read from order and number."),
    ("F1", "F1_field", "Transport and field",
     "How influence is transported along the order and conserved at every event, and what field a source produces: "
     "Gauss's law, the far field, superposition and its failure at second order, and the role of capacity."),
    ("F2", "F2_motion", "Motion in a field",
     "How a history moves in a given order: whether maximal (longest) chains follow the geodesics of the effective "
     "geometry, in 2+1 and 3+1."),
    ("F3", "F3_structure", "Formation of structure and matter",
     "What matter is in a causal order: flows that carry their four-velocity, attempts at a bundle bound by its own "
     "field, locked defects of a carried orientation, and defects (kinks, walls) of a field on all events evolved by "
     "the retarded Green function of the order - their identity, inertia, swerves - and the Green function itself."),
    ("F4", "F4_inertia", "Inertia and contact",
     "How matter resists and exchanges motion: flows that meet merge with their four-momentum conserved."),
    ("F5", "F5_gravity", "Gravity and equivalence",
     "Whether matter falls in an order shaped by other matter, without a rule for falling added: longest chains in a "
     "loaded order, field defects in orders of non-uniform density, under a causal growth rule and on prescribed "
     "curved orders; the active side (how matter shapes the order); the equivalence of gravitational and inertial "
     "mass."),
    ("F6", "F6_order_dynamics", "Dynamics of the order",
     "How the order itself responds and what action weights it: local growth rules, the phases of the "
     "Benincasa-Dowker action, the expected BD action to second order in a static field (the Einstein-Hilbert "
     "stiffness of its conformal, mixed and optical sectors, the non-local age terms, the finite-amplitude response), "
     "the vacuum control, and the source of the order from a matter action that depends on the order alone - together "
     "the field equation of the order."),
    ("F7", "F7_light", "Light and the response of the order",
     "How the response of the order fixes the post-Newtonian gamma: fixed order, narrowing cones, conservation of causal "
     "throughput; and the tensor coupling (pressure and motion gravitate, cones tilt)."),
    ("F8", "F8_data", "Confrontation with data",
     "Predictions of the theory compared with observation: the static field of a point mass near the non-locality "
     "scale, the solar system, and the astrophysical consequences."),
]
ID_RE = re.compile(r"^(K|F\d)\.\d+[a-z]?\d?$")


def questions(test_list):
    q = {}
    for line in Path(test_list).read_text(encoding="utf-8").replace("\r", "").splitlines():
        if line.startswith("| ") and not line.startswith("| ID") and not line.startswith("|--"):
            c = [x.strip() for x in line.strip().strip("|").split("|")]
            if len(c) >= 4:
                q[re.sub(r"\s*\(.*\)$", "", c[0]).split()[0]] = (c[2], c[3])
    return q


def readme(phase, folder, name, objective, ids, q):
    lines = [f"# {phase}: {name}", "", "## Objective", "", objective, "",
             "## Contents", "",
             f"One folder per test, `{folder}/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_"
             "<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines "
             "and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is "
             "overwritten.", "",
             "The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in "
             "`docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.", "",
             "| Test | Question | Level |", "|---|---|---|"]
    for i in ids:
        qq, lev = q.get(i, ("(not in docs/test_list.md)", ""))
        lines.append(f"| {i} | {qq} | {lev} |")
    return "\n".join(lines) + "\n"


def main():
    root, tl = Path(sys.argv[1]), sys.argv[2]
    apply = "--apply" in sys.argv
    res = root / "results"
    q = questions(tl)
    order = {p: i for i, (p, *_) in enumerate(PHASES)}
    def key(t):
        ph, rest = t.split(".", 1)
        m = re.match(r"(\d+)(.*)", rest)
        return (order.get(ph, 99), int(m.group(1)), m.group(2))
    moves = []
    for d in sorted(res.iterdir()):
        if d.is_dir() and ID_RE.match(d.name):
            ph = d.name.split(".")[0]
            folder = next(f for p, f, *_ in PHASES if p == ph)
            moves.append((d, res / folder / d.name))
    for src, dst in moves:
        if dst.exists():
            raise SystemExit(f"target exists, nothing done: {dst}")
    print(f"{len(moves)} test folders to move")
    for p, folder, name, obj in PHASES:
        ids = sorted([s.name for s, t in moves if t.parent.name == folder] +
                     [x.name for x in (res / folder).glob("*") if x.is_dir() and ID_RE.match(x.name)], key=key)
        print(f"  {folder}: {len(ids)} tests")
        if apply:
            (res / folder).mkdir(exist_ok=True)
            (res / folder / "README.md").write_text(readme(p, folder, name, obj, ids, q), encoding="utf-8")
    if apply:
        for src, dst in moves:
            shutil.move(str(src), str(dst))
        print("moved")


if __name__ == "__main__":
    main()
