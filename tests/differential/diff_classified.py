#!/usr/bin/env python3
"""Differential check: C++ full LW classification (geometry + template RMSD +
ambiguity resolution) vs the Python finder classify block.

Replicates finder._score_candidates' classification for each valid candidate
(canonical res_id order) and compares (final lw_class, template_rmsd, swapped).
The finder helpers are copied here verbatim because importing finder.py pulls
scipy (PEP668 system Python).

Usage: diff_classified.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import json
import math
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
FP2 = Path("/Users/kdewan2/Desktop/Projects/find_pair_2")
IDEALIZED_DIR = FP2 / "basepair-idealized"
EXEMPLAR_DIR = FP2 / "basepair-catalog-exemplars"

import numpy as np  # noqa: E402

from pair_finder.classification.aligner import TemplateAligner  # noqa: E402
from pair_finder.classification.edge_classifier import (  # noqa: E402
    canonicalize_lw, classify_lw_class, match_saenger,
)
from pair_finder.core import Residue  # noqa: E402
from pair_finder.core.residue import normalize_base_type  # noqa: E402
from pair_finder.hbond.finder import HBondFinder  # noqa: E402
from pair_finder.io.json_loaders import load_reference_frames  # noqa: E402
from pair_finder.validation import GeometricValidator  # noqa: E402

# Load quality_intrinsic.py directly (importing pair_finder.finder pulls scipy).
import importlib.util as _ilu  # noqa: E402
_qi_spec = _ilu.spec_from_file_location(
    "quality_intrinsic", Path(PY_SRC) / "pair_finder" / "finder" / "quality_intrinsic.py")
_qi = _ilu.module_from_spec(_qi_spec)
_qi_spec.loader.exec_module(_qi)
intrinsic_confidence = _qi.intrinsic_confidence

_qs_spec = _ilu.spec_from_file_location(
    "quality_scorer", Path(PY_SRC) / "pair_finder" / "finder" / "quality_scorer.py")
_qs = _ilu.module_from_spec(_qs_spec)
_qs_spec.loader.exec_module(_qs)
QUALITY_SCORER = _qs.QualityScorer()

_QS_BASE_ATOMS = frozenset(
    {"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"})


def _sugar_fb(r):
    return 0.0 if r >= 3.0 else round(0.30 * max(0.0, 1.0 - r / 3.0), 3)


def _hh_fb(r):
    return 0.0 if r >= 4.0 else round(0.35 * max(0.0, 1.0 - r / 4.0), 3)


def selection_quality(validation, sequence, strict, lw, template_rmsd, amb_classes, amb_rmsds):
    hbl = [{"donor_atom": hb.donor_atom, "acceptor_atom": hb.acceptor_atom,
            "distance": hb.distance, "alignment": hb.alignment_score,
            "context": ("base_base" if hb.donor_atom in _QS_BASE_ATOMS
                        and hb.acceptor_atom in _QS_BASE_ATOMS else "other")}
           for hb in strict]
    if amb_classes:
        scores = []
        for cls, rmsd in zip(amb_classes, amb_rmsds):
            s = QUALITY_SCORER.compute_score(validation, sequence, hbonds=hbl,
                                             rmsd=rmsd, lw_class=cls)
            if s == 0.0 and rmsd is not None:
                if len(cls) == 3 and "S" in cls:
                    s = _sugar_fb(rmsd)
                elif len(cls) == 3 and cls[1] == "H" and cls[2] == "H":
                    s = _hh_fb(rmsd)
            scores.append(s)
        return max(scores)
    s = QUALITY_SCORER.compute_score(validation, sequence, hbonds=hbl,
                                     rmsd=template_rmsd, lw_class=lw)
    if s == 0.0 and template_rmsd is not None and lw is not None:
        if len(lw) == 3 and "S" in lw:
            s = _sugar_fb(template_rmsd)
        elif len(lw) == 3 and lw[1] == "H" and lw[2] == "H":
            s = _hh_fb(template_rmsd)
    return s

_ALL_LW_CLASSES = ("cWW", "cWH", "cWS", "cHH", "cHS", "cSS",
                   "tWW", "tWH", "tWS", "tHH", "tHS", "tSS")

_FACE_PRIORITY = {"W": 0, "H": 1, "S": 2}
_SYMMETRIC = {"cWW", "tWW", "cHH", "tHH", "cSS", "tSS"}
_RMSD_RESOLUTION_THRESHOLD = 0.15
RMSD_TOL = 1e-3


def compute_template_rmsd(aligner, res1, res2, lw_class):
    if len(lw_class) == 3:
        o, a, b = lw_class[0], lw_class[1], lw_class[2]
        if _FACE_PRIORITY[a] > _FACE_PRIORITY[b]:
            lw_class = f"{o}{b}{a}"
            res1, res2 = res2, res1
    sequence = normalize_base_type(res1.base_type) + normalize_base_type(res2.base_type)
    path = aligner.loader.find_template(lw_class, sequence)
    if path:
        rmsd, _ = aligner.align_to_template(res1, res2, path)
        return rmsd
    if lw_class in _SYMMETRIC:
        rev = sequence[1] + sequence[0]
        path = aligner.loader.find_template(lw_class, rev)
        if path:
            rmsd, _ = aligner.align_to_template(res2, res1, path)
            return rmsd
    return None


def resolve_by_o2prime(classes, res1, res2, hbonds):
    if not hbonds:
        return None
    r1 = r2 = False
    for hb in hbonds:
        if hb.donor_atom == "O2'" and hb.donor_res_id == res1.res_id:
            r1 = True
        elif hb.acceptor_atom == "O2'" and hb.acceptor_res_id == res1.res_id:
            r1 = True
        if hb.donor_atom == "O2'" and hb.donor_res_id == res2.res_id:
            r2 = True
        elif hb.acceptor_atom == "O2'" and hb.acceptor_res_id == res2.res_id:
            r2 = True
    if not r1 and not r2:
        return None
    edges0 = (classes[0][1], classes[0][2])
    edges1 = (classes[1][1], classes[1][2])
    for pos in (0, 1):
        if edges0[pos] == edges1[pos]:
            continue
        s0, s1 = edges0[pos] == "S", edges1[pos] == "S"
        if not (s0 ^ s1):
            continue
        if r1 or r2:
            return classes[0] if s0 else classes[1]
    return None


def classification_confidence(aligner, r1, r2, sequence, primary_class, hbonds):
    orientation = primary_class[0]
    class_rmsds = []
    rmsd_by_class = {}
    for cls in (c for c in _ALL_LW_CLASSES if c[0] == orientation):
        r = compute_template_rmsd(aligner, r1, r2, cls)
        if r is not None:
            class_rmsds.append((cls, r))
            rmsd_by_class[cls] = r
    geo_uncertain = "|" in primary_class
    is_ambig = geo_uncertain and len(class_rmsds) >= 2
    if is_ambig:
        parts = primary_class.split("|")
        parts.sort(key=lambda x: rmsd_by_class.get(x, float("inf")))
        top_two = parts[:2]
        display = f"ambiguous ({top_two[0]}|{top_two[1]})"
        named = top_two
    else:
        display = primary_class
        named = [primary_class]
    confs = []
    for cls in named:
        score = intrinsic_confidence(
            hbonds=hbonds, plane_angle=0.0, d_v=0.0,
            template_rmsd=rmsd_by_class.get(cls), lw_class=cls,
            sequence=sequence, rmsd_by_class=rmsd_by_class)
        confs.append(round(score, 3))
    return display, confs, is_ambig


def resolve_ambiguity(classes, rmsds, res1, res2, hbonds):
    if len(classes) != 2:
        return None
    if hbonds:
        sr = match_saenger(res1, res2, hbonds)
        if sr is not None:
            f1, f2, _, _ = sr
            sc, _ = canonicalize_lw(classes[0][0], f1, f2)
            if sc in classes:
                return sc
    o = resolve_by_o2prime(classes, res1, res2, hbonds)
    if o is not None:
        return o
    r0, r1 = rmsds[0], rmsds[1]
    if r0 is not None and r1 is not None and abs(r0 - r1) > _RMSD_RESOLUTION_THRESHOLD:
        return classes[0] if r0 < r1 else classes[1]
    return None


def python_classified(stem):
    frames = load_reference_frames(JSON_DIR / "ls_fitting" / f"{stem}.json")
    data = json.load(open(JSON_DIR / "pdb_atoms" / f"{stem}.json"))
    res_atoms = {}
    for atom in data[0]["atoms"]:
        res_atoms.setdefault(atom["res_id"], []).append(atom)
    residues = {}
    for res_id, atoms in res_atoms.items():
        name = res_id.split("-")[1]
        bt = name[1] if (name.startswith("D") and len(name) == 2) else name
        r = Residue(res_id=res_id, base_type=bt)
        for a in atoms:
            r.add_atom(a["atom_name"], np.array(a["xyz"], dtype=np.float64))
        residues[res_id] = r

    validator = GeometricValidator()
    finder = HBondFinder()
    aligner = TemplateAligner(IDEALIZED_DIR, EXEMPLAR_DIR)
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
            strict = finder.find_between(residues[ra], residues[rb])
            lw, swapped = classify_lw_class(
                frames[ra], frames[rb], residues[ra], residues[rb], hbonds=voters)
            a1, a2 = (residues[rb], residues[ra]) if swapped else (residues[ra], residues[rb])

            template_rmsd = None
            amb_classes, amb_rmsds = None, None
            if "|" in lw:
                classes = lw.split("|")
                rmsds = [compute_template_rmsd(aligner, a1, a2, cls) for cls in classes]
                resolved = resolve_ambiguity(classes, rmsds, a1, a2, strict)
                if resolved is not None:
                    lw = resolved
                    template_rmsd = rmsds[classes.index(resolved)]
                elif "cWW" in classes and any(r is None for r in rmsds):
                    lw = "cWW"
                    template_rmsd = rmsds[classes.index("cWW")]
                else:
                    template_rmsd = rmsds[0]
                    amb_classes, amb_rmsds = classes, rmsds
            else:
                template_rmsd = compute_template_rmsd(aligner, a1, a2, lw)
            sequence = a1.base_type + a2.base_type
            display, confs, is_ambig = classification_confidence(
                aligner, a1, a2, sequence, lw, strict)
            quality = selection_quality(v, sequence, strict, lw, template_rmsd,
                                        amb_classes, amb_rmsds)
            out[(ra, rb)] = (lw, template_rmsd, bool(swapped), display, tuple(confs),
                             bool(is_ambig), round(quality, 3))
    return out


def cpp_classified(binary, pdb):
    res = subprocess.run([binary, "dump-classified", str(pdb)],
                         capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        f = line.split("\t")
        rmsd = None if f[3] == "None" else (math.inf if f[3] == "inf" else float(f[3]))
        confs = tuple(float(x) for x in f[6].split(",")) if f[6] else tuple()
        out[(f[0], f[1])] = (f[2], rmsd, bool(int(f[4])), f[5], confs, bool(int(f[7])),
                             float(f[8]))
    return out


def rmsd_eq(a, b):
    if a is None or b is None:
        return a is None and b is None
    if math.isinf(a) or math.isinf(b):
        return math.isinf(a) and math.isinf(b)
    return abs(a - b) <= RMSD_TOL


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    binary, pdbs = sys.argv[1], sys.argv[2:]
    overall_ok = True
    for pdb in pdbs:
        stem = Path(pdb).stem
        o = python_classified(stem)
        c = cpp_classified(binary, pdb)
        only_py, only_cpp = set(o) - set(c), set(c) - set(o)
        cls_bad, rmsd_bad, conf_bad, qual_bad = [], [], [], []
        for k in (set(o) & set(c)):
            if o[k][0] != c[k][0] or o[k][2] != c[k][2]:
                cls_bad.append((k, o[k], c[k]))
            elif not rmsd_eq(o[k][1], c[k][1]):
                rmsd_bad.append((k, o[k], c[k]))
            else:
                conf_ok = (o[k][3] == c[k][3] and o[k][5] == c[k][5]
                           and len(o[k][4]) == len(c[k][4])
                           and all(abs(a - b) <= 1.5e-3 for a, b in zip(o[k][4], c[k][4])))
                if not conf_ok:
                    conf_bad.append((k, o[k], c[k]))
                elif abs(o[k][6] - c[k][6]) > 1.5e-3:
                    qual_bad.append((k, o[k], c[k]))
        ok = (not only_py and not only_cpp and not cls_bad and not rmsd_bad
              and not conf_bad and not qual_bad)
        overall_ok &= ok
        tag = "OK  " if ok else "DIFF"
        print(f"  {tag} {stem:8s}  py={len(o)} cpp={len(c)} only_py={len(only_py)} "
              f"only_cpp={len(only_cpp)} class/swap_bad={len(cls_bad)} "
              f"rmsd_bad={len(rmsd_bad)} conf_bad={len(conf_bad)} qual_bad={len(qual_bad)}")
        if not ok:
            for k in list(only_py)[:2]:
                print(f"        only_py : {k} {o[k]}")
            for k, ov, cv in (cls_bad + rmsd_bad + conf_bad + qual_bad)[:6]:
                print(f"        mismatch {k}: py={ov} cpp={cv}")
    print("RESULT:", "ALL MATCH" if overall_ok else "DIVERGENCE")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
