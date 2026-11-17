"""
F6.8f development: common plumbing of the diag scripts (decides nothing)
=======================================================================

The conventions of the other diag scripts (diag/b2_dev.py, diag/c0_minvar.py):
  * run from the repository root as a module, python -m diag.<name> [options], with argparse options;
  * --workers N (0 = CAUSALITY_WORKERS or all cores) through core.stats.parallel.parallel_map, results in the order of
    the jobs, every job a pure function of its arguments;
  * every job cached in diag/out/<name>_cache/<key>.json (atomic write), so an interrupted run resumes; the key
    contains the md5 of the code the job depends on, so a changed code never reuses an old value;
  * the output of a run is diag/out/<name>_<UTC time>.json (arguments, md5 of the code, data) and .txt (the printed
    summary); nothing is overwritten.
Scripts that use the 80-bit instrument of F6.8e3 refuse to run without an 80-bit long double (numpy on Windows has
none), as F6.8e3 itself does: run them under Linux or WSL.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from typing import Any, Callable, Dict, Iterable, List

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "diag", "out")


def md5_file(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def code_md5(paths: Iterable[str]) -> Dict[str, str]:
    """md5 of the files a computation depends on, keyed by the path relative to the repository root."""
    return {os.path.relpath(p, ROOT).replace(os.sep, "/"): md5_file(p) for p in paths}


def ld_guard() -> None:
    """Stop if numpy has no 80-bit long double (the F6.8e3 instrument needs it)."""
    if np.finfo(np.longdouble).nmant < 63:
        raise SystemExit("this script uses the 80-bit instrument of F6.8e3; numpy on this platform has no 80-bit long "
                         "double (Windows): run it under Linux or WSL")


def _plain(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return _plain(x.tolist())
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.complexfloating, complex)):
        return [float(np.real(x)), float(np.imag(x))]
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def cached(name: str, code: Dict[str, str], fn: Callable, args: tuple) -> Any:
    """fn(*args), cached per job in diag/out/<name>_cache/ under a key of (name, code md5, args)."""
    key = hashlib.sha256(json.dumps([name, code, _plain(list(args))], sort_keys=True).encode()).hexdigest()[:24]
    path = os.path.join(OUT, f"{name}_cache", key + ".json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)["value"]
    value = _plain(fn(*args))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump({"args": _plain(list(args)), "value": value}, fh)
    os.replace(path + ".tmp", path)
    return value


def save(name: str, args: Dict[str, Any], code: Dict[str, str], data: Any, lines: List[str], t0: float) -> str:
    """Write diag/out/<name>_<UTC time>.json and .txt; returns the stem."""
    os.makedirs(OUT, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    stem = os.path.join(OUT, f"{name}_{stamp}")
    i = 1
    while os.path.exists(stem + ".json"):           # two runs in the same second: never overwrite
        stem = os.path.join(OUT, f"{name}_{stamp}_{i}")
        i += 1
    rec = {"script": name, "time_utc": stamp, "duration_s": round(time.time() - t0, 1), "args": _plain(args),
           "code_md5": code, "data": _plain(data)}
    with open(stem + ".json", "w") as fh:
        json.dump(rec, fh, indent=1)
    with open(stem + ".txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return os.path.relpath(stem, ROOT).replace(os.sep, "/")


def floats(s: str) -> List[float]:
    return [float(v) for v in s.split(",") if v.strip()]


def pairs(s: str) -> List[tuple]:
    """'9.87:0,16:4' -> [(9.87, 0.0), (16.0, 4.0)]"""
    out = []
    for item in s.split(","):
        a, b = item.split(":")
        out.append((float(a), float(b)))
    return out
