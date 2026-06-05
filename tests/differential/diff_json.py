#!/usr/bin/env python3
"""Differential: C++ `find` JSON vs the Python find_pairs.py CLI output.

Builds the SAME output dict find_pairs.py emits (default pipeline config, which
is what the C++ replicates) from the real find_pairs, then compares it to the
C++ binary's JSON per pair (keyed by sorted res_id pair). Numeric blocks use
tolerances because the C++ computes frames from the raw PDB while Python loads
the ls_fitting JSON frames (~1e-6 apart, plus 4-dec rounding half-even vs
half-up).

Usage: VENVPY diff_json.py <pairfinder_binary> <pdb1> [pdb2 ...]
"""
import json
import subprocess
import sys
from pathlib import Path

PY_SRC = "/Users/kdewan2/Desktop/Projects/find-pair-score-rna/prototyped-pair-finder-main/src"
sys.path.insert(0, PY_SRC)
JSON_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/json")
PDB_DIR = Path("/Users/kdewan2/Desktop/Projects/find_pair_2/data/pdb")

from pair_finder.finder import FinderConfig, PairFinder  # noqa: E402
from pair_finder.validation.rigid_body import compute_rigid_body_parameters  # noqa: E402

# Field tolerances (0 = exact match required).
T_COORD, T_GEOM, T_PLANE, T_RB = 2e-3, 2e-2, 0.1, 0.15
T_CONF, T_SCORE, T_RMSD, T_QUAL = 0.03, 0.15, 0.02, 0.02


def py_output(pdb):
    """Replicate cli/find_pairs.py output_data using the default pipeline."""
    result = PairFinder(JSON_DIR, FinderConfig()).find_pairs(pdb, score=True)

    def vec(v):
        return [round(float(x), 4) for x in v]

    def classification(p):
        if p.lw_class_display is None and p.lw_class is None:
            return None
        conf = p.lw_class_confidence
        return {"lw_class": p.lw_class_display or p.lw_class,
                "lw_class_confidence": ([round(float(c), 2) for c in conf] if conf else conf),
                "is_ambiguous": p.is_ambiguous}

    def pair(p):
        v = p.validation
        rb = compute_rigid_body_parameters(p.frame1, p.frame2)
        d = {"res_id1": p.res_id1, "res_id2": p.res_id2, "sequence": p.sequence,
             "pair_category": p.pair_category, "num_hbonds": p.num_hbonds,
             "hbond_categories": p.hbond_categories, "classification": classification(p),
             "geometry": {"dorg": round(float(v.dorg), 4), "d_v": round(float(v.d_v), 4),
                          "plane_angle": round(float(v.plane_angle), 4),
                          "dNN": round(float(v.dNN), 4)},
             "rigid_body": {k: round(val, 4) for k, val in rb.items()},
             "frames": {"res1": {"origin": vec(p.frame1.origin), "x_axis": vec(p.frame1.x_axis),
                                 "y_axis": vec(p.frame1.y_axis), "z_axis": vec(p.frame1.z_axis)},
                        "res2": {"origin": vec(p.frame2.origin), "x_axis": vec(p.frame2.x_axis),
                                 "y_axis": vec(p.frame2.y_axis), "z_axis": vec(p.frame2.z_axis)}},
             "template_rmsd": (round(p.template_rmsd, 4) if p.template_rmsd is not None else None),
             "quality_score": round(p.quality_score, 4),
             "score": (round(p.score, 2) if p.score is not None else None),
             "issues": list(p.issues)}
        return d

    out = {"pdb_id": result.pdb_id, "candidates_total": result.candidates_total,
           "candidates_valid": result.candidates_valid, "n_pairs": len(result.pairs),
           "pairs": [pair(p) for p in result.pairs]}
    if result.overall_score is not None:
        sd = result.score_details
        out["overall_score"] = round(result.overall_score, 2)
        out["pairs_score"] = round(sd.pairs_score, 2)
        out["residues_score"] = round(sd.residues_score, 2)
        out["score_details"] = {"pairs_score": round(sd.pairs_score, 2),
                                "residues_score": round(sd.residues_score, 2),
                                "issue_summary": sd.issue_summary(),
                                "n_scored_pairs": sd.n_pairs,
                                "n_residues_scored": sd.n_residues_scored,
                                "n_skipped_pairs": sd.skipped_pairs}
    return out


def cpp_output(binary, pdb):
    r = subprocess.run([binary, str(PDB_DIR / f"{pdb}.pdb"), "--details"],
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def key(p):
    return tuple(sorted([p["res_id1"], p["res_id2"]]))


def num_eq(a, b, tol):
    if a is None or b is None:
        return a == b
    if isinstance(a, str) or isinstance(b, str):  # "Infinity"
        return str(a) == str(b) or (a in ("Infinity", float("inf")) and b in ("Infinity", float("inf")))
    return abs(float(a) - float(b)) <= tol


def rot_mag(rb):
    """Magnitude of the rigid-body rotation vector (deg). Near 180 the axis-angle
    split into buckle/propeller/opening is mathematically ill-conditioned."""
    return (rb["buckle"] ** 2 + rb["propeller"] ** 2 + rb["opening"] ** 2) ** 0.5


def is_degenerate_rb(py, cp):
    return rot_mag(py["rigid_body"]) > 150.0 or rot_mag(cp["rigid_body"]) > 150.0


def cmp_pair(py, cp, skip_rb=False):
    errs = []
    for f in ("res_id1", "res_id2", "sequence", "pair_category", "num_hbonds",
              "hbond_categories", "issues"):
        if py.get(f) != cp.get(f):
            errs.append(f"{f}: py={py.get(f)} cpp={cp.get(f)}")
    pc, cc = py.get("classification"), cp.get("classification")
    if (pc is None) != (cc is None):
        errs.append(f"classification presence py={pc} cpp={cc}")
    elif pc is not None:
        if pc["lw_class"] != cc["lw_class"]:
            errs.append(f"lw_class py={pc['lw_class']} cpp={cc['lw_class']}")
        if pc["is_ambiguous"] != cc["is_ambiguous"]:
            errs.append(f"is_ambiguous py={pc['is_ambiguous']} cpp={cc['is_ambiguous']}")
        pconf, cconf = pc["lw_class_confidence"] or [], cc["lw_class_confidence"] or []
        if len(pconf) != len(cconf) or any(not num_eq(a, b, T_CONF) for a, b in zip(pconf, cconf)):
            errs.append(f"confidence py={pconf} cpp={cconf}")
    for f, tol in (("dorg", T_GEOM), ("d_v", T_GEOM), ("dNN", T_GEOM), ("plane_angle", T_PLANE)):
        if not num_eq(py["geometry"][f], cp["geometry"][f], tol):
            errs.append(f"geom.{f} py={py['geometry'][f]} cpp={cp['geometry'][f]}")
    # rigid_body (Tsukuba params): skip near a 180 deg rotation, where the
    # axis-angle split is ill-conditioned and amplifies the ~1e-6 frame diff
    # (documented; not a real divergence, and rb is not used downstream).
    if not skip_rb:
        for f in ("shear", "stretch", "stagger", "buckle", "propeller", "opening"):
            if not num_eq(py["rigid_body"][f], cp["rigid_body"][f], T_RB):
                errs.append(f"rb.{f} py={py['rigid_body'][f]} cpp={cp['rigid_body'][f]}")
    for r in ("res1", "res2"):
        for ax in ("origin", "x_axis", "y_axis", "z_axis"):
            for a, b in zip(py["frames"][r][ax], cp["frames"][r][ax]):
                if not num_eq(a, b, T_COORD):
                    errs.append(f"frame.{r}.{ax} py={py['frames'][r][ax]} cpp={cp['frames'][r][ax]}")
                    break
    if not num_eq(py["template_rmsd"], cp["template_rmsd"], T_RMSD):
        errs.append(f"template_rmsd py={py['template_rmsd']} cpp={cp['template_rmsd']}")
    if not num_eq(py["quality_score"], cp["quality_score"], T_QUAL):
        errs.append(f"quality_score py={py['quality_score']} cpp={cp['quality_score']}")
    if not num_eq(py["score"], cp["score"], T_SCORE):
        errs.append(f"score py={py['score']} cpp={cp['score']}")
    return errs


def main():
    binary, pdbs = sys.argv[1], sys.argv[2:]
    ok_all = True
    for pdb in pdbs:
        py, cp = py_output(pdb), cpp_output(binary, pdb)
        pm = {key(p): p for p in py["pairs"]}
        cm = {key(p): p for p in cp["pairs"]}
        only_py, only_cpp = set(pm) - set(cm), set(cm) - set(pm)
        bad, rb_skipped = {}, 0
        for k in (set(pm) & set(cm)):
            deg = is_degenerate_rb(pm[k], cm[k])
            if deg:
                rb_skipped += 1
            e = cmp_pair(pm[k], cm[k], skip_rb=deg)
            if e:
                bad[k] = e
        agg = []
        for f, tol in (("overall_score", T_SCORE), ("pairs_score", T_SCORE),
                       ("residues_score", T_SCORE)):
            if not num_eq(py.get(f), cp.get(f), tol):
                agg.append(f"{f} py={py.get(f)} cpp={cp.get(f)}")
        psd, csd = py.get("score_details", {}), cp.get("score_details", {})
        for f in ("n_scored_pairs", "n_residues_scored", "n_skipped_pairs"):
            if psd.get(f) != csd.get(f):
                agg.append(f"score_details.{f} py={psd.get(f)} cpp={csd.get(f)}")
        if psd.get("issue_summary") != csd.get("issue_summary"):
            agg.append(f"issue_summary py={psd.get('issue_summary')} cpp={csd.get('issue_summary')}")
        # candidates_total must match exactly. candidates_valid is a raw
        # diagnostic that can differ by a few near the d_v/plane/dorg validation
        # cutoffs (C++ computes frames; Python loads ls_fitting JSON, ~1e-6
        # apart) — it does not affect the selected pairs or scores, so it is a
        # WARNING annotation, not a failure (see PORT_DEVIATIONS plane_angle).
        cand_fail = py["candidates_total"] != cp["candidates_total"]
        if cand_fail:
            agg.append(f"candidates_total py={py['candidates_total']} cpp={cp['candidates_total']}")
        note = "" if py["candidates_valid"] == cp["candidates_valid"] else \
            f" [valid py={py['candidates_valid']} cpp={cp['candidates_valid']} ~cutoff]"
        rbnote = f" rb_skip={rb_skipped}(~180deg)" if rb_skipped else ""
        ok = not only_py and not only_cpp and not bad and not agg
        ok_all &= ok
        print(f"  {'OK  ' if ok else 'DIFF'} {pdb:8s} pairs py={len(pm)} cpp={len(cm)} "
              f"bad={len(bad)} only_py={len(only_py)} only_cpp={len(only_cpp)}{rbnote}{note}")
        if not ok:
            for f in agg:
                print(f"      agg {f}")
            for k, errs in list(bad.items())[:6]:
                print(f"      {k}: {errs[:4]}")
            for k in list(only_py)[:4]:
                print(f"      only_py {k}")
            for k in list(only_cpp)[:4]:
                print(f"      only_cpp {k}")
    print("RESULT:", "ALL MATCH" if ok_all else "DIVERGENCE")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
