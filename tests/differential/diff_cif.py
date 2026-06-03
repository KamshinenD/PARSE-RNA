#!/usr/bin/env python3
"""CIF<->PDB parity: the GEMMI loader must give identical results for the same
structure whether read from .pdb or .cif (proves CIF support matches the
validated PDB path).

Usage: diff_cif.py <pairfinder_binary> <pdb_file> <cif_file> [<pdb> <cif> ...]
Args are alternating (pdb, cif) pairs of the SAME structure.
"""
import json
import subprocess
import sys


def run(binary, path):
    r = subprocess.run([binary, path, "--details"], capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def main():
    binary = sys.argv[1]
    pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
    ok_all = True
    for pdb_path, cif_path in pairs:
        d_pdb, d_cif = run(binary, pdb_path), run(binary, cif_path)
        # pdb_id differs (file stem); compare everything else.
        d_pdb.pop("pdb_id", None)
        d_cif.pop("pdb_id", None)
        ok = d_pdb == d_cif
        ok_all &= ok
        name = cif_path.split("/")[-1]
        print(f"  {'OK  ' if ok else 'DIFF'} {name:14s} "
              f"pdb n_pairs={d_pdb['n_pairs']} cif n_pairs={d_cif['n_pairs']} "
              f"overall pdb={d_pdb.get('overall_score')} cif={d_cif.get('overall_score')}")
        if not ok:
            for k in sorted(set(d_pdb) | set(d_cif)):
                if d_pdb.get(k) != d_cif.get(k):
                    pv, cv = d_pdb.get(k), d_cif.get(k)
                    if k == "pairs":
                        print(f"      pairs differ: pdb has {len(pv)}, cif has {len(cv)}")
                    else:
                        print(f"      {k}: pdb={pv} cif={cv}")
    print("RESULT:", "ALL MATCH" if ok_all else "DIVERGENCE")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
