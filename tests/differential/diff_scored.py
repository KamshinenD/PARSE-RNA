#!/usr/bin/env python3
"""Differential: C++ full scoring pass (classify + score + O2 correction) vs the
REAL Python find_pairs all_candidates. Run under the scipy venv.

Compares per valid candidate (canonical res_id key): pair_category, lw_class,
template_rmsd, quality_score, num_hbonds, has_strong_base_hbond, display label,
confidences, is_ambiguous, and _precorr_lw (the O2-correction marker).

Usage: VENVPY diff_scored.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import math
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

DTOL = 1.5e-3


def python_scored(pdb):
    fr = PairFinder(JSON_DIR, FinderConfig(min_score=0.0,
                                           selection_strategy="helix_priority")).find_pairs(pdb)
    out = {}
    for c in fr.all_candidates:
        if not c.validation.is_valid:
            continue
        key = tuple(sorted([c.res_id1, c.res_id2]))
        rmsd = c.template_rmsd
        out[key] = (
            c.pair_category,
            c.lw_class or "None",
            rmsd,
            round(c.quality_score, 3),
            c.num_hbonds,
            bool(c.has_strong_base_hbond),
            c.lw_class_display or "None",
            tuple(round(x, 3) for x in (c.lw_class_confidence or [])),
            bool(c.is_ambiguous),
            getattr(c, "_precorr_lw", None) or "None",
        )
    return out


def cpp_scored(binary, pdb):
    res = subprocess.run([binary, "dump-scored", str(PDB_DIR / f"{pdb}.pdb")],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        rmsd = None if f[4] == "None" else (math.inf if f[4] == "inf" else float(f[4]))
        confs = tuple(float(x) for x in f[9].split(",")) if f[9] != "None" else tuple()
        out[(f[0], f[1])] = (f[2], f[3], rmsd, float(f[5]), int(f[6]), bool(int(f[7])),
                             f[8], confs, bool(int(f[10])), f[11])
    return out


def rmsd_eq(a, b):
    if a is None or b is None:
        return a is None and b is None
    if math.isinf(a) or math.isinf(b):
        return math.isinf(a) and math.isinf(b)
    return abs(a - b) <= DTOL


def cmp(o, c):
    # o, c: (cat, lw, rmsd, qual, nhb, strong, disp, conf, amb, precorr)
    if o[0] != c[0] or o[1] != c[1] or o[4] != c[4] or o[5] != c[5]:
        return "class/cat/hb"
    if not rmsd_eq(o[2], c[2]):
        return "rmsd"
    if abs(o[3] - c[3]) > DTOL:
        return "quality"
    if o[6] != c[6] or o[8] != c[8] or o[9] != c[9]:
        return "display/amb/precorr"
    if len(o[7]) != len(c[7]) or any(abs(a - b) > DTOL for a, b in zip(o[7], c[7])):
        return "confidence"
    return None


def main():
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall = True
    for pdb in pdbs:
        o = python_scored(pdb)
        c = cpp_scored(binary, pdb)
        only_py, only_cpp = set(o) - set(c), set(c) - set(o)
        bad = [(k, cmp(o[k], c[k])) for k in (set(o) & set(c)) if cmp(o[k], c[k])]
        ok = not only_py and not only_cpp and not bad
        overall &= ok
        tag = "OK  " if ok else "DIFF"
        corr = sum(1 for v in o.values() if v[9] != "None")
        print(f"  {tag} {pdb:8s} py={len(o)} cpp={len(c)} only_py={len(only_py)} "
              f"only_cpp={len(only_cpp)} bad={len(bad)} o2corr={corr}")
        if not ok:
            for k in list(only_py)[:3]:
                print(f"      only_py: {k} {o[k]}")
            for k in list(only_cpp)[:3]:
                print(f"      only_cpp: {k} {c[k]}")
            for k, why in bad[:6]:
                print(f"      [{why}] {k}\n        py ={o[k]}\n        cpp={c[k]}")
    print("RESULT:", "ALL MATCH" if overall else "DIVERGENCE")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
