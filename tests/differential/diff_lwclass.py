#!/usr/bin/env python3
"""Differential check: C++ LW edge classifier vs Python classify_lw_class.

For every valid candidate (canonical res_id order), builds face-voter H-bonds and
runs classify_lw_class on both sides, comparing (lw_class, swapped). C++ computes
frames from the PDB; the Python oracle uses ls_fitting + pdb_atoms JSON with the
production (PairCache) base types.

Usage: diff_lwclass.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import json
import os
import subprocess
import sys
from pathlib import Path

PY_SRC = os.environ.get(
    "PAIR_FINDER_SRC",
    str(Path(__file__).resolve().parents[3] / "prototyped-pair-finder-main" / "src"),
)
sys.path.insert(0, PY_SRC)
JSON_DIR = Path(os.environ.get("PAIR_FINDER_JSON_DIR",
                               "/Users/kdewan2/Desktop/Projects/find_pair_2/data/json"))

import numpy as np  # noqa: E402

from pair_finder.classification.edge_classifier import classify_lw_class  # noqa: E402
from pair_finder.core import Residue  # noqa: E402
from pair_finder.hbond.finder import HBondFinder  # noqa: E402
from pair_finder.io.json_loaders import load_reference_frames  # noqa: E402
from pair_finder.validation import GeometricValidator  # noqa: E402


def _extract_res_name(res_id):
    name = res_id.split("-")[1]
    return name[1] if (name.startswith("D") and len(name) == 2) else name


def python_lwclass(stem):
    frames = load_reference_frames(JSON_DIR / "ls_fitting" / f"{stem}.json")
    data = json.load(open(JSON_DIR / "pdb_atoms" / f"{stem}.json"))
    res_atoms = {}
    for atom in data[0]["atoms"]:
        res_atoms.setdefault(atom["res_id"], []).append(atom)
    residues = {}
    for res_id, atoms in res_atoms.items():
        r = Residue(res_id=res_id, base_type=_extract_res_name(res_id))
        for a in atoms:
            r.add_atom(a["atom_name"], np.array(a["xyz"], dtype=np.float64))
        residues[res_id] = r

    validator = GeometricValidator()
    finder = HBondFinder()
    res_ids = list(frames.keys())
    origins = {r: np.asarray(frames[r].origin, dtype=np.float64) for r in res_ids}
    out = {}
    for i in range(len(res_ids)):
        for j in range(i + 1, len(res_ids)):
            ri, rj = res_ids[i], res_ids[j]
            if float(np.linalg.norm(origins[ri] - origins[rj])) > 15.0:
                continue
            r1, r2 = residues.get(ri), residues.get(rj)
            if r1 is None or r2 is None:
                continue
            n1, n2 = r1.get_glycosidic_n(), r2.get_glycosidic_n()
            if n1 is None or n2 is None:
                continue
            v = validator.validate(frames[ri], frames[rj], n1.coords, n2.coords)
            if not v.is_valid:
                continue
            ra, rb = (ri, rj) if ri < rj else (rj, ri)
            voters = finder.find_face_voters(residues[ra], residues[rb])
            lw, swapped = classify_lw_class(
                frames[ra], frames[rb], residues[ra], residues[rb], hbonds=voters,
            )
            out[(ra, rb)] = (lw, bool(swapped))
    return out


def cpp_lwclass(binary, pdb):
    res = subprocess.run([binary, "dump-lwclass", str(pdb)],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        out[(f[0], f[1])] = (f[2], bool(int(f[3])))
    return out


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        stem = Path(pdb).stem
        o = python_lwclass(stem)
        c = cpp_lwclass(binary, pdb)
        only_py, only_cpp = set(o) - set(c), set(c) - set(o)
        bad = [(k, o[k], c[k]) for k in (set(o) & set(c)) if o[k] != c[k]]
        ok = (not only_py and not only_cpp and not bad)
        overall_ok &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {stem:8s}  py={len(o)} cpp={len(c)} only_py={len(only_py)} "
              f"only_cpp={len(only_cpp)} class_mismatch={len(bad)}")
        if not ok:
            for k in list(only_py)[:3]:
                print(f"        only_py : {k} {o[k]}")
            for k in list(only_cpp)[:3]:
                print(f"        only_cpp: {k} {c[k]}")
            for k, ov, cv in bad[:6]:
                print(f"        mismatch {k}: py={ov} cpp={cv}")
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
