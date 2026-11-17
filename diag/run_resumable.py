"""
Resumable start of a frozen test in Claude's workspace (an execution aid; decides nothing, changes no frozen code)
==================================================================================================================

The workspace stops a background job after 2 h; K.18 needs more. This wrapper imports the test module under the name
the registry uses (core/registry.py: phases.<folder>.<stem>), replaces its module-level worker functions (the names
given with --workers-fn, the pure functions the test passes to parallel_map) by a disk-cached version, and then starts
the ordinary CLI (core.cli, the same as run.py). parallel_map forks on Linux, so the workers see the cached version.

Cache: diag/out/run_resumable_cache/<test id>/<key>.pkl, key = sha256(function name, md5 of the test file and of the
files given with --code, repr(args)); pickle keeps the exact return values and their types. An interrupted run is
started again with the same command; the jobs already done are read back, the rest is computed. The final record is
written by the CLI as usual; the code md5 in the key guarantees that a changed file never reuses a value.

    python -m diag.run_resumable K.18 --module phases.K_instruments.K_18_needle_volume \
        --workers-fn _jacobi,_series,_exact,_brute,_e3 --code diag/f68f_1b.py,... -- K.18 --workers 2
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import importlib
import os
import pickle
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def md5(path):
    with open(os.path.join(ROOT, path), "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def main():
    argv = sys.argv[1:]
    cut = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("test_id")
    ap.add_argument("--module", required=True)
    ap.add_argument("--workers-fn", required=True)
    ap.add_argument("--code", default="")
    a = ap.parse_args(argv[:cut])
    mod = importlib.import_module(a.module)
    files = [os.path.relpath(mod.__file__, ROOT)] + [c for c in a.code.split(",") if c]
    code = {f.replace(os.sep, "/"): md5(f) for f in files}
    cdir = os.path.join(ROOT, "diag", "out", "run_resumable_cache", a.test_id)
    os.makedirs(cdir, exist_ok=True)
    print(f"[run_resumable] {a.test_id}: cache {os.path.relpath(cdir, ROOT)}; code {code}", flush=True)
    for name in a.workers_fn.split(","):
        fn = getattr(mod, name)

        @functools.wraps(fn)
        def wrapped(*args, _fn=fn, _name=name):
            key = hashlib.sha256(repr((_name, sorted(code.items()), args)).encode()).hexdigest()[:32]
            path = os.path.join(cdir, key + ".pkl")
            if os.path.exists(path):
                with open(path, "rb") as fh:
                    return pickle.load(fh)
            value = _fn(*args)
            with open(path + ".tmp", "wb") as fh:
                pickle.dump(value, fh)
            os.replace(path + ".tmp", path)
            return value

        wrapped.__module__ = mod.__name__
        wrapped.__qualname__ = name
        setattr(mod, name, wrapped)
    from core.cli import main as cli
    return cli(argv[cut + 1:])


if __name__ == "__main__":
    raise SystemExit(main())
