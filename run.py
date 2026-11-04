"""
Causality Lab -- start a test from a terminal (no installation needed).

    python run.py list
    python run.py F1.2 --part d3
    python run.py F1.2 --quick
    python run.py core
    python run.py help           all commands
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main() -> int:
    # a console that cannot encode a character (e.g. Greek letters on a Windows code page) prints "?" instead of
    # stopping the run with UnicodeEncodeError
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    from core.cli import main as cli
    return cli(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
