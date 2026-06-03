#!/usr/bin/env python3
"""Differential check: C++ slot-based H-bond finder vs the Python HBondFinder.

For each PDB, enumerates every residue pair (canonical res_id order) and runs
find_between on both sides, then compares the resulting H-bond sets:
(donor_res, donor_atom, acceptor_res, acceptor_atom) identity, distance,
chosen slot indices, and alignment score.

Usage: diff_hbonds.py <pairfinder_binary> <pdb1> [pdb2 ...]
Needs the Python package importable (PAIR_FINDER_SRC or default sibling path).
"""
import os
import subprocess
import sys
from pathlib import Path

PY_SRC = os.environ.get(
    "PAIR_FINDER_SRC",
    str(Path(__file__).resolve().parents[3] / "prototyped-pair-finder-main" / "src"),
)
sys.path.insert(0, PY_SRC)

from pair_finder.hbond.finder import HBondFinder  # noqa: E402
from pair_finder.io.pdb_parser import parse_pdb  # noqa: E402

DIST_TOL = 1e-3
ALIGN_TOL = 1e-3

# The C++ loader now follows find_pair_2's GEMMI conventions (drops waters,
# renames legacy atom spellings). Apply the same to the parse_pdb oracle so this
# stays an apples-to-apples test of the FINDER logic on an identical residue set.
_WATERS = {"HOH", "WAT", "H2O", "OH2", "SOL"}
_ATOM_RENAME = {"OP1": "O1P", "OP2": "O2P", "OP3": "O3P", "OL": "O1P", "OR": "O2P",
                "O1'": "O4'", "C5A": "C5M", "O5T": "O5'", "O3T": "O3'"}


def _norm_atom(name):
    name = name.replace("*", "'")
    return _ATOM_RENAME.get(name, name)


def _is_water(res):
    parts = res.res_id.split("-")
    return len(parts) >= 2 and parts[1].upper() in _WATERS


def python_hbonds(pdb_path):
    residues = parse_pdb(Path(pdb_path), use_cache=False)
    for res in residues.values():
        renamed = {_norm_atom(n): a for n, a in res.atoms.items()}
        res.atoms = renamed
    ordered = [residues[k] for k in sorted(residues) if not _is_water(residues[k])]
    finder = HBondFinder()
    out = {}
    for i in range(len(ordered)):
        for j in range(i + 1, len(ordered)):
            for hb in finder.find_between(ordered[i], ordered[j]):
                a, b = ordered[i].res_id, ordered[j].res_id
                ki, kj = (a, b) if a < b else (b, a)
                key = (ki, kj, hb.donor_res_id, hb.donor_atom,
                       hb.acceptor_res_id, hb.acceptor_atom)
                out[key] = (hb.distance, hb.h_slot_idx, hb.lp_slot_idx,
                            hb.alignment_score)
    return out


def cpp_hbonds(binary, pdb_path):
    res = subprocess.run([binary, "dump-hbonds", str(pdb_path)],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        key = (f[0], f[1], f[2], f[3], f[4], f[5])
        out[key] = (float(f[6]), int(f[7]), int(f[8]), float(f[9]))
    return out


def compare(o, c):
    ok = set(o) & set(c)
    only_py, only_cpp = set(o) - set(c), set(c) - set(o)
    bad = 0
    sample = None
    for k in ok:
        do, ho, lo, ao = o[k]
        dc, hc, lc, ac = c[k]
        if abs(do - dc) > DIST_TOL or ho != hc or lo != lc or abs(ao - ac) > ALIGN_TOL:
            bad += 1
            if sample is None:
                sample = (k, o[k], c[k])
    return only_py, only_cpp, bad, sample


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        stem = Path(pdb).stem
        o = python_hbonds(pdb)
        c = cpp_hbonds(binary, pdb)
        only_py, only_cpp, bad, sample = compare(o, c)
        ok = (not only_py and not only_cpp and bad == 0)
        overall_ok &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {stem:8s}  py={len(o)} cpp={len(c)} "
              f"only_py={len(only_py)} only_cpp={len(only_cpp)} bad={bad}")
        if not ok:
            for k in list(only_py)[:3]:
                print(f"        only_py : {k}  {o[k]}")
            for k in list(only_cpp)[:3]:
                print(f"        only_cpp: {k}  {c[k]}")
            if sample:
                print(f"        mismatch: {sample}")
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
