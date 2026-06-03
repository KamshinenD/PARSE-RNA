#!/usr/bin/env python3
"""Differential: C++ empirical pair scores vs the real Python Scorer.

Runs find_pairs(include_context=True) + Scorer.score_finder_result and compares
per-pair scores (keyed by sorted res_id pair). Run under the scipy venv.

Usage: VENVPY diff_scores.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import sys
from pathlib import Path

PY_SRC = "/Users/kdewan2/Desktop/Projects/find-pair-score-rna/prototyped-pair-finder-main/src"
sys.path.insert(0, PY_SRC)
JSON_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/json")
PDB_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/pdb")

import subprocess  # noqa: E402

from pair_finder.finder import FinderConfig, PairFinder  # noqa: E402
from pair_finder.scoring import Scorer  # noqa: E402

TOL = 0.05


def python_scores(pdb):
    fr = PairFinder(JSON_DIR, FinderConfig(min_score=0.0,
                                           selection_strategy="helix_priority")).find_pairs(
        pdb, include_context=True, score=False)
    # Fresh Scorer per PDB: its cached HBondFinder keys by res_id, so reusing it
    # across structures would serve stale atoms for shared res_ids (production
    # uses one scorer per structure).
    ss = Scorer.load_default().score_finder_result(fr, pdb_id=pdb)
    out = {}
    for p in ss.pair_scores:
        out[tuple(sorted([p.res_id1, p.res_id2]))] = p.score
    return out, ss


def cpp_scores(binary, pdb):
    res = subprocess.run([binary, "dump-scores", str(PDB_DIR / f"{pdb}.pdb")],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        out[(f[0], f[1])] = float(f[2])
    return out


def main():
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall = True
    for pdb in pdbs:
        o, ss = python_scores(pdb)
        c = cpp_scores(binary, pdb)
        only_py, only_cpp = set(o) - set(c), set(c) - set(o)
        bad = [(k, o[k], c[k]) for k in (set(o) & set(c)) if abs(o[k] - c[k]) > TOL]
        ok = not only_py and not only_cpp and not bad
        overall &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {pdb:8s} py={len(o)} cpp={len(c)} only_py={len(only_py)} "
              f"only_cpp={len(only_cpp)} score_bad={len(bad)}  (py overall={ss.overall})")
        if not ok:
            for k in list(only_py)[:3]:
                print(f"      only_py : {k} {o[k]}")
            for k in list(only_cpp)[:3]:
                print(f"      only_cpp: {k} {c[k]}")
            for k, ov, cv in bad[:6]:
                print(f"      score {k}: py={ov:.3f} cpp={cv:.3f}")
    print("RESULT:", "ALL MATCH" if overall else "DIVERGENCE")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
