"""Dump the complete, unfiltered QuadoaCore API to a file.

    python tools/dump_api.py [--grep TERM ...]

The API is documented only in these docstrings, so this is the reference for
questions like "can a configuration be created programmatically?". Note Quadoa
calls a configuration a *realisation* in its docs, so search for both.

Needs a free license seat.
"""

from __future__ import annotations

import argparse
import inspect
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "quadoa_api_reference.txt"))
    ap.add_argument("--grep", nargs="*", default=None,
                    help="also print methods whose name or doc matches any term")
    args = ap.parse_args()

    from beamdown.config import load_config
    from beamdown.session import _import_quadoa

    cfg = load_config()
    quadoa = _import_quadoa(cfg.trace.quadoa_folder)
    core = quadoa.QuadoaCore()
    try:
        methods = [
            (name, member)
            for name, member in inspect.getmembers(core, predicate=inspect.ismethod)
            if not name.startswith("_")
        ]
        lines = []
        for name, member in methods:
            doc = inspect.getdoc(member) or "(no docstring)"
            lines.append(f"{name}\n{'-' * len(name)}\n{doc}\n")
        text = "\n".join(lines)
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"{len(methods)} public methods -> {args.out}")

        if args.grep:
            terms = [t.lower() for t in args.grep]
            print(f"\nmethods matching {args.grep}:")
            hits = 0
            for name, member in methods:
                doc = (inspect.getdoc(member) or "").lower()
                if any(t in name.lower() or t in doc for t in terms):
                    hits += 1
                    first = (inspect.getdoc(member) or "").splitlines()
                    sig = first[0] if first else ""
                    print(f"  {name}\n      {sig}")
            if not hits:
                print("  (none)")
        else:
            print("\nall method names:")
            for i in range(0, len(methods), 3):
                print("  " + "".join(f"{n:38s}" for n, _ in methods[i:i + 3]))
    finally:
        import gc

        del core
        gc.collect()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
