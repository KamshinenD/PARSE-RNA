#!/usr/bin/env python3
"""Differential: C++ whole-structure score vs the real Python Scorer.

Compares overall / pairs_score / residues_score and per-residue backbone scores.
Run under the scipy venv.

Usage: VENVPY diff_structure.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import subprocess
import sys
from pathlib import Path

PY_SRC = "/Users/kdewan2/Desktop/Projects/find-pair-score-rna/prototyped-pair-finder-main/src"
sys.path.insert(0, PY_SRC)
JSON_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/json")
PDB_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/pdb")

from pair_finder.finder import FinderConfig, PairFinder  # noqa: E402
from pair_finder.scoring import Scorer  # noqa: E402

TOL = 0.05


def python_struct(pdb):
    fr = PairFinder(JSON_DIR, FinderConfig(min_score=0.0,
                                           selection_strategy="helix_priority")).find_pairs(
        pdb, include_context=True, score=False)
    ss = Scorer.load_default().score_finder_result(fr, pdb_id=pdb)
    res = {r.res_id: round(r.score, 2) for r in ss.residue_scores}
    return ss.overall, ss.pairs_score, ss.residues_score, res


def cpp_struct(binary, pdb):
    out = subprocess.run([binary, "dump-structure", str(PDB_DIR / f"{pdb}.pdb")],
                         capture_output=True, text=True, check=True)
    overall = pairs = residues = None
    res = {}
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if f[0] == "OVERALL":
            overall, pairs, residues = float(f[1]), float(f[2]), float(f[3])
        elif f[0] == "RES":
            res[f[1]] = float(f[2])
    return overall, pairs, residues, res


def main():
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        po, pp, pr, pres = python_struct(pdb)
        co, cp, cr, cres = cpp_struct(binary, pdb)
        only_py, only_cpp = set(pres) - set(cres), set(cres) - set(pres)
        res_bad = [(k, pres[k], cres[k]) for k in (set(pres) & set(cres))
                   if abs(pres[k] - cres[k]) > TOL]
        agg_bad = abs(po - co) > TOL or abs(pp - cp) > TOL or abs(pr - cr) > TOL
        ok = not only_py and not only_cpp and not res_bad and not agg_bad
        overall_ok &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {pdb:8s} overall py={po}/cpp={co} pairs={pp}/{cp} res={pr}/{cr} "
              f"res_bad={len(res_bad)} only_py={len(only_py)} only_cpp={len(only_cpp)}")
        if not ok:
            for k, pv, cv in res_bad[:6]:
                print(f"      res {k}: py={pv} cpp={cv}")
            for k in list(only_py)[:3]:
                print(f"      only_py: {k} {pres[k]}")
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
