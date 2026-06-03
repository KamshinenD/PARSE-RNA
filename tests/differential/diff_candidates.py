#!/usr/bin/env python3
"""Differential check: C++ candidate finder + validator vs Python PairCache.

C++ computes frames from the PDB then finds/validates candidates; the Python
oracle builds PairCache from the JSON (ls_fitting + pdb_atoms). Compares the
candidate set (keyed by the sorted res_id pair) and each pair's validation
geometry (dorg, d_v, plane_angle, dNN, direction cosines, quality, is_valid).

Usage: diff_candidates.py <pairfinder_binary> <pdb1> [pdb2 ...]
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
JSON_DIR = Path(os.environ.get("PAIR_FINDER_JSON_DIR",
                               "/Users/kdewan2/Desktop/Projects/find_pair_2/data/json"))

# Avoid importing pair_finder.finder.cache (pulls scipy); replicate PairCache's
# logic directly from the same pieces (frames JSON + pdb_atoms JSON + validator).
# KDTree query_ball_point(r=15) == brute-force origin distance <= 15 (same metric).
import json  # noqa: E402

import numpy as np  # noqa: E402

from pair_finder.core import Residue  # noqa: E402
from pair_finder.io.json_loaders import load_reference_frames  # noqa: E402
from pair_finder.validation import GeometricValidator  # noqa: E402

DTOL = 1e-3      # distances (dorg, d_v, dNN)
COSTOL = 2e-4    # direction cosines
QTOL = 3e-3      # quality score (carries the amplified plane_angle/180 term)
# plane_angle = degrees(acos(|dir_z|)); near |dir_z|=1 (parallel bases) acos is
# hyper-sensitive, so the ~1e-6 frame agreement (computed C++ vs JSON-loaded
# Python) amplifies to ~1e-2 deg here. This region is far from the 70-deg
# validity threshold, so is_valid is unaffected (checked exactly below).
PLANE_TOL = 0.05


def _extract_res_name(res_id):
    name = res_id.split("-")[1]
    return name[1] if (name.startswith("D") and len(name) == 2) else name


def python_candidates(stem):
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
            key = (ri, rj) if ri < rj else (rj, ri)
            out[key] = (v.dorg, v.d_v, v.plane_angle, v.dNN, v.dir_x, v.dir_y,
                        v.dir_z, v.quality_score, bool(v.is_valid))
    return out


def cpp_candidates(binary, pdb):
    res = subprocess.run([binary, "dump-candidates", str(pdb)],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        key = (f[0], f[1])
        out[key] = (float(f[2]), float(f[3]), float(f[4]), float(f[5]),
                    float(f[6]), float(f[7]), float(f[8]), float(f[9]),
                    bool(int(f[10])))
    return out


def bad_value(o, c):
    # o, c: (dorg,d_v,plane,dNN,dx,dy,dz,qscore,valid)
    for i in (0, 1, 3):  # dorg, d_v, dNN
        if abs(o[i] - c[i]) > DTOL:
            return True
    if abs(o[2] - c[2]) > PLANE_TOL:  # plane_angle (arccos-amplified)
        return True
    for i in (4, 5, 6):
        if abs(o[i] - c[i]) > COSTOL:
            return True
    if abs(o[7] - c[7]) > QTOL:
        return True
    return o[8] != c[8]  # is_valid must match exactly


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        stem = Path(pdb).stem
        try:
            o = python_candidates(stem)
        except Exception as e:  # noqa: BLE001
            print(f"  SKIP {stem:8s}  (python: {e})")
            continue
        c = cpp_candidates(binary, pdb)
        only_py, only_cpp = set(o) - set(c), set(c) - set(o)
        bad = sum(1 for k in (set(o) & set(c)) if bad_value(o[k], c[k]))
        nvalid_py = sum(1 for v in o.values() if v[8])
        nvalid_cpp = sum(1 for v in c.values() if v[8])
        ok = (not only_py and not only_cpp and bad == 0)
        overall_ok &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {stem:8s}  py={len(o)} cpp={len(c)} valid(py/cpp)="
              f"{nvalid_py}/{nvalid_cpp} only_py={len(only_py)} "
              f"only_cpp={len(only_cpp)} bad={bad}")
        if not ok:
            for k in list(only_py)[:3]:
                print(f"        only_py : {k}  {o[k]}")
            for k in list(only_cpp)[:3]:
                print(f"        only_cpp: {k}  {c[k]}")
            for k in (set(o) & set(c)):
                if bad_value(o[k], c[k]):
                    print(f"        mismatch {k}\n          py ={o[k]}\n          cpp={c[k]}")
                    break
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
