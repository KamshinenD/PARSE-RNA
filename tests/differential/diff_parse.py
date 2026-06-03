#!/usr/bin/env python3
"""Differential check: C++ pdb parser vs the Python oracle (pdb_parser.parse_pdb).

For each PDB, both sides emit "<res_id>\\t<atom>\\t%.3f\\t%.3f\\t%.3f" sorted, and
we diff. Exit non-zero if any structure diverges.

Usage:
    diff_parse.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY_SRC = HERE.parents[2] / "prototyped-pair-finder-main" / "src"
sys.path.insert(0, str(PY_SRC))

from pair_finder.io.pdb_parser import parse_pdb  # noqa: E402


def python_rows(pdb_path):
    residues = parse_pdb(Path(pdb_path), use_cache=False)
    rows = []
    for res_id, res in residues.items():
        for name, atom in res.atoms.items():
            c = atom.coords
            rows.append(f"{res_id}\t{name}\t{c[0]:.3f}\t{c[1]:.3f}\t{c[2]:.3f}")
    return sorted(rows)


def cpp_rows(binary, pdb_path):
    out = subprocess.run([binary, "dump-parse", str(pdb_path)],
                         capture_output=True, text=True, check=True)
    return sorted(l for l in out.stdout.splitlines() if l)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        py = python_rows(pdb)
        cpp = cpp_rows(binary, pdb)
        pset, cset = set(py), set(cpp)
        only_py = pset - cset
        only_cpp = cset - pset
        name = Path(pdb).stem
        if not only_py and not only_cpp and len(py) == len(cpp):
            print(f"  OK   {name:8s}  {len(py):>6} atoms identical")
        else:
            overall_ok = False
            print(f"  DIFF {name:8s}  py={len(py)} cpp={len(cpp)} "
                  f"only_py={len(only_py)} only_cpp={len(only_cpp)}")
            for r in list(only_py)[:4]:
                print(f"        py-only : {r}")
            for r in list(only_cpp)[:4]:
                print(f"        cpp-only: {r}")
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
