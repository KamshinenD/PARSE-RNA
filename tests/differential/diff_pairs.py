#!/usr/bin/env python3
"""Differential: C++ full pipeline SELECTED pairs vs the real Python
find_pairs().pairs. Run under the scipy venv.

Compares the final selected pair set (keyed by sorted res_id pair) and lw_class.

Usage: VENVPY diff_pairs.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import os
import subprocess
import sys
from pathlib import Path

PY_SRC = os.environ.get(
    "PAIR_FINDER_SRC",
    "/Users/kdewan2/Desktop/Projects/find-pair-score-rna/prototyped-pair-finder-main/src")
sys.path.insert(0, PY_SRC)
JSON_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/json")
PDB_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/pdb")

from pair_finder.finder import FinderConfig, PairFinder  # noqa: E402


def python_pairs(pdb):
    fr = PairFinder(JSON_DIR, FinderConfig(min_score=0.0, selection_strategy="helix_priority")).find_pairs(pdb)
    out = {}
    for c in fr.pairs:
        out[tuple(sorted([c.res_id1, c.res_id2]))] = c.lw_class or "None"
    return out


def cpp_pairs(binary, pdb):
    res = subprocess.run([binary, "dump-pairs", str(PDB_DIR / f"{pdb}.pdb")],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        out[(f[0], f[1])] = f[2]
    return out


def main():
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall = True
    tot_py = tot_cpp = tot_match = tot_lwbad = 0
    for pdb in pdbs:
        o = python_pairs(pdb)
        c = cpp_pairs(binary, pdb)
        only_py, only_cpp = set(o) - set(c), set(c) - set(o)
        lw_bad = [(k, o[k], c[k]) for k in (set(o) & set(c)) if o[k] != c[k]]
        ok = not only_py and not only_cpp and not lw_bad
        overall &= ok
        tot_py += len(o); tot_cpp += len(c)
        tot_match += len(set(o) & set(c)) - len(lw_bad); tot_lwbad += len(lw_bad)
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {pdb:8s} py={len(o)} cpp={len(c)} only_py={len(only_py)} "
              f"only_cpp={len(only_cpp)} lw_bad={len(lw_bad)}")
        if not ok:
            for k in list(only_py)[:4]:
                print(f"      only_py : {k} {o[k]}")
            for k in list(only_cpp)[:4]:
                print(f"      only_cpp: {k} {c[k]}")
            for k, ov, cv in lw_bad[:4]:
                print(f"      lw {k}: py={ov} cpp={cv}")
    print(f"\nTOTAL py={tot_py} cpp={tot_cpp} exact_match={tot_match} lw_mismatch={tot_lwbad}")
    print("RESULT:", "ALL MATCH" if overall else "DIVERGENCE")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
