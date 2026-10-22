"""
core.cli
========

Commands behind run.py.

    python run.py list                     phases, tests, latest status
    python run.py <ID>                     run a test (all its parts)
          --part P                         only one part (e.g. --part d3)
          --quick                          small smoke run: checks that the test runs, decides nothing
          --seed K                         independent replication with base seed K (default 0)
          --workers N                      worker processes (default: CAUSALITY_WORKERS or all cores)
          --set key=value                  change a parameter: any with --quick; in a full run only the
                                           parameters the part declares settable (python run.py show <ID>)
    python run.py show <ID>                the test's header: question, prediction, criterion, ...
    python run.py results <ID>             the saved results of a test
    python run.py core                     fast checks of the core and the instruments
    python run.py report                   regenerate docs/seznam_testu.md
"""

from __future__ import annotations

import ast
import sys
import time
from typing import Any, Dict, List, Optional

from . import provenance, registry
from .stats.parallel import default_workers
from .testing import Context, Result, STATUS_CS, LEVEL_CS, check_overrides, merged_params, parts_of, settable_of


def _out(msg: str) -> None:
    print(msg, flush=True)


def _parse(argv: List[str]) -> Dict[str, Any]:
    opts: Dict[str, Any] = {"part": None, "quick": False, "seed": 0, "workers": None, "set": {}}
    rest: List[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--part":
            opts["part"] = argv[i + 1]
            i += 2
        elif a == "--quick":
            opts["quick"] = True
            i += 1
        elif a == "--seed":
            opts["seed"] = int(argv[i + 1])
            i += 2
        elif a == "--workers":
            opts["workers"] = int(argv[i + 1])
            i += 2
        elif a == "--set":
            k, v = argv[i + 1].split("=", 1)
            try:
                opts["set"][k] = ast.literal_eval(v)
            except (ValueError, SyntaxError):
                opts["set"][k] = v
            i += 2
        else:
            rest.append(a)
            i += 1
    opts["rest"] = rest
    return opts


def latest(tid: str, part: Optional[str]) -> Optional[Dict[str, Any]]:
    """Newest non-quick record of a test part."""
    recs = [r for r in provenance.load_all(tid) if r.get("part") == part]
    return recs[-1] if recs else None


def status_line(tid: str, entry: Dict[str, Any]) -> str:
    """Latest status per part, with replication and 'older code' marks."""
    test, path = entry["test"], entry["path"]
    cur = provenance.sha256_file(path)
    bits = []
    for part in parts_of(test):
        recs = [r for r in provenance.load_all(tid) if r.get("part") == part]
        if not recs:
            bits.append(f"{part or ''}: -" if part else "-")
            continue
        r = recs[-1]
        seeds = {x["seed"] for x in recs if x["status"] == r["status"]}
        mark = STATUS_CS[r["status"]]
        if len(seeds) > 1:
            mark += " (replikováno)"
        if r["test_sha256"] != cur:
            mark += " [starší verze testu]"
        bits.append(f"{part}: {mark}" if part else mark)
    return "; ".join(bits)


def cmd_list() -> int:
    tests = registry.discover()
    for ph, folder, name, part in registry.PHASES:
        ids = sorted([t for t in tests if registry.phase_of(t) == ph], key=registry.sort_key)
        _out(f"{ph:3s} {name} – {part}")
        if not ids:
            _out("      (zatím bez testů; plán v docs/roadmapa.md)")
        for tid in ids:
            t = tests[tid]["test"]
            old = f" ({t['old_id']})" if t.get("old_id") else ""
            _out(f"      {tid:7s}{old:9s} {t['title']}  [{LEVEL_CS[t['level']]}]  {status_line(tid, tests[tid])}")
    return 0


def cmd_show(tid: str) -> int:
    tests = registry.discover()
    if tid not in tests:
        raise SystemExit(f"unknown test {tid}")
    t = tests[tid]["test"]
    for k in ("id", "title", "phase", "theory_part", "question", "method", "prediction", "criterion",
              "null_model", "level", "meters", "depends_on", "old_id", "estimate"):
        if k in t:
            _out(f"{k:12s} {t[k]}")
    _out(f"{'params':12s} {t['params']}")
    for part, pp in (t.get("parts") or {}).items():
        _out(f"{'part ' + part:12s} {pp}")
    _out(f"{'quick':12s} {t['quick']}")
    for part, keys in (t.get("settable") or {}).items():
        _out(f"{'settable ' + part:12s} {keys}  (--set in a full run)")
    return 0


def cmd_results(tid: str) -> int:
    recs = provenance.load_all(tid, include_quick=True)
    if not recs:
        _out(f"{tid}: no results yet")
        return 0
    for r in recs:
        q = " QUICK" if r.get("quick") else ""
        _out(f"{r['time_utc']}  part {r.get('part') or '-'}  seed {r['seed']}{q}  {r['status']:17s} "
             f"test {r['test_sha256'][:8]}  {r['duration_s']} s  {r['_file']}")
    return 0


def cmd_run(tid: str, opts: Dict[str, Any]) -> int:
    tests = registry.discover()
    if tid not in tests:
        raise SystemExit(f"unknown test {tid} (python run.py list)")
    entry = tests[tid]
    test, mod, path = entry["test"], entry["module"], entry["path"]
    parts = parts_of(test)
    if opts["part"]:
        if opts["part"] not in parts:
            raise SystemExit(f"{tid} has parts {parts}")
        parts = [opts["part"]]
    if opts["set"] and not opts["quick"]:
        # a full run takes from --set only the parameters its part declares settable; checked for every part
        # BEFORE anything runs, so a long chain never stops halfway on a bad option
        if len(parts) > 1:
            raise SystemExit("--set in a full run needs --part (settable parameters are declared per part)")
        why = check_overrides(test, parts[0], False, opts["set"])
        if why:
            raise SystemExit(why)
    workers = default_workers(opts["workers"])
    # fingerprints taken BEFORE the run: a file replaced while a long run is going must not be recorded as the
    # version that ran (found 19 Oct 2026: F1.5 ran for ~2.7 h while a newer F1.5 file was delivered)
    hashes = {"test_sha256": provenance.sha256_file(path), "core_sha256": provenance.core_sha256()}
    for part in parts:
        p = merged_params(test, part, opts["quick"], opts["set"])
        ctx = Context(part=part or "", seed=opts["seed"], quick=opts["quick"], workers=workers, progress=_out)
        head = f"{tid} {test['title']}" + (f"  part={part}" if part else "")
        mode = "QUICK (decides nothing)" if opts["quick"] else "FULL"
        _out(head)
        _out(f"  {mode}; seed {opts['seed']}; workers {workers}; level {test['level']}")
        _out(f"  params {p}")
        if opts["set"] and not opts["quick"]:
            _out(f"  set at the command line (settable {settable_of(test, part)}): {opts['set']}")
        t0 = time.time()
        res: Result = mod.run(p, ctx)
        for line in res.lines:
            _out(line)
        f = provenance.save(tid, path, part, opts["seed"], opts["quick"], p, res.status, res.lines, res.data, t0,
                            hashes=hashes)
        _out(f"\nRESULT {tid}{(' [' + part + ']') if part else ''}: {res.status} ({STATUS_CS[res.status]})"
             f"{'  (quick run)' if opts['quick'] else ''}")
        _out(f"saved: {f.relative_to(provenance.ROOT)}  ({time.time() - t0:.0f} s)\n")
    return 0


def main(argv: List[str]) -> int:
    if not argv or argv[0] in ("-h", "--help", "help"):
        _out(__doc__)
        return 0
    cmd = argv[0]
    opts = _parse(argv[1:])
    if cmd == "list":
        return cmd_list()
    if cmd == "show":
        return cmd_show(opts["rest"][0])
    if cmd == "results":
        return cmd_results(opts["rest"][0])
    if cmd == "core":
        from core_tests.checks import main as checks
        return checks(opts["rest"])
    if cmd == "report":
        from .report import write_report
        p = write_report()
        _out(f"written: {p.relative_to(provenance.ROOT)}")
        return 0
    return cmd_run(cmd, opts)
