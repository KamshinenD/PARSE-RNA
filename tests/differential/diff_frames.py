#!/usr/bin/env python3
"""Differential check: C++ base frames vs the ls_fitting oracle (generate_modern_json).

For each PDB, compares per-residue rotation_matrix / translation / rms_fit from
the C++ `dump-frames` output against json_legacy ls_fitting/<pdb>.json.

Usage: diff_frames.py <pairfinder_binary> <pdb1> [pdb2 ...]
The ls_fitting JSON is looked up by stem under LS_FITTING_DIR.
"""
import json
import subprocess
import sys
from pathlib import Path

LS_FITTING_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/json/ls_fitting")
ROT_TOL = 2e-3
TRANS_TOL = 2e-2
RMS_TOL = 5e-3


def oracle(pdb_stem):
    data = json.load(open(LS_FITTING_DIR / f"{pdb_stem}.json"))
    out = {}
    for e in data:
        out[e["res_id"]] = (e["rotation_matrix"], e["translation"], e["rms_fit"])
    return out


def cpp(binary, pdb_path):
    res = subprocess.run([binary, "dump-frames", str(pdb_path)],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        res_id, rms = f[0], float(f[1])
        vals = [float(x) for x in f[2:11]]
        rot = [vals[0:3], vals[3:6], vals[6:9]]
        trans = [float(x) for x in f[11:14]]
        out[res_id] = (rot, trans, rms)
    return out


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        stem = Path(pdb).stem
        o = oracle(stem)
        c = cpp(binary, pdb)
        common = set(o) & set(c)
        only_o, only_c = set(o) - set(c), set(c) - set(o)
        max_rot = max_tr = max_rms = 0.0
        worst = None
        bad = 0
        for rid in common:
            ro, to, rmso = o[rid]
            rc, tc, rmsc = c[rid]
            dr = max(abs(ro[i][j] - rc[i][j]) for i in range(3) for j in range(3))
            dt = max(abs(to[k] - tc[k]) for k in range(3))
            drms = abs(rmso - rmsc)
            if dr > max_rot:
                max_rot, worst = dr, rid
            max_tr = max(max_tr, dt)
            max_rms = max(max_rms, drms)
            if dr > ROT_TOL or dt > TRANS_TOL or drms > RMS_TOL:
                bad += 1
        ok = (bad == 0 and not only_o and not only_c)
        overall_ok &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {stem:8s}  res common={len(common)} only_o={len(only_o)} "
              f"only_c={len(only_c)} | maxΔrot={max_rot:.2e} maxΔtr={max_tr:.2e} "
              f"maxΔrms={max_rms:.2e} bad={bad}")
        if not ok and worst:
            print(f"        worst rot diff at {worst}")
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
