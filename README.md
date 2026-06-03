# pairfinder (C++)

A fast, faithful C++ port of the Python RNA base-pair finder + scorer
(`prototyped-pair-finder-main`). Goal: same results as the Python tool, far
faster, open-sourceable. Follows the architecture/pattern of the reference
`x3dna/` codebase (typed atoms, value types, dependency-injected algorithms,
observer/strategy patterns), but ports its **logic from the current Python**
(x3dna is behind).

## Build
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/pairfinder <input.pdb>
```
Requires a C++20 compiler and CMake ≥ 3.20.

## Layout (mirrors x3dna)
```
include/pairfinder/   public headers
  geometry/           vector3d
  core/               atom_type, atom, residue, chain, structure
  io/                 (pdb/cif parser, json reader/writer)        [todo]
  algorithms/         frames, hbonds (slot), validation,
                      edge classification, pair finder, selection [todo]
  scoring/            issues, prosco, penalty, scorer             [todo]
  protocols/          find_pair + score protocol                  [todo]
src/                  implementations mirroring include/
apps/main.cpp         CLI
tests/                unit tests + differential harness vs Python
```

## Porting principle: faithful first, differentially tested
The Python implementation is the **oracle**. Each layer is ported bottom-up and
its output diffed against Python on a fixed PDB set; divergence is driven to
zero before moving up. Refactors come only after parity.

### Order
1. [x] skeleton + CMake + core value types
2. [x] io: structure parse -> Structure. GEMMI v0.6.5 loader (PDB / mmCIF / .gz),
       replicating find_pair_2's conventions (first model, alt-loc ' '/'A'/'1',
       modified-nt HETATM filter, OP1->O1P atom renames). diff_cif: .cif==.pdb
       byte-identical (1GID, 1EHZ). Reads full ribosomes (8ife 80S, 219K atoms).
3. [x] core: reference frames (per-base templates + type-check gate + altloc)
       (diff vs ls_fitting oracle: 400 PDBs exact)
4. [x] algorithms: slot-based H-bond detection
       (diff vs Python HBondFinder: 25 PDBs exact, incl protein-RNA)
5. [x] candidate finding + geometric validation
       (diff vs Python PairCache: 120 PDBs exact)
6. [x] edge classification (geometry + template RMSD + ambiguity + confidence)
       6a geometric classifier (220 PDBs exact)
       6b template RMSD  6c ambiguity + Saenger  6d confidence (v3 two-term)
       (full classify diff vs Python: lw_class + rmsd + swap + display + conf)
7. [x] selection (helix-priority phases)
       carbonyl-O2 correction + full scoring pass + RNAChains + _CWWHelixPhase
       + GlobalOptimal MWIS + HelixPriority
       (diff SELECTED pairs vs REAL find_pairs().pairs: 150 PDBs / 3904 pairs exact)
8. [x] scoring: load reference JSON, 6 issues, ProSco, penalty,
       per-pair / per-residue / per-RNA
       (diff_scores 70 PDBs exact; diff_structure overall+per-residue exact,
        incl. negative-seq / insertion-code residues now scored both ends)
9. [x] CLI: pdb -> JSON of pairs + scores
       `pairfinder <in.pdb> [--out FILE] [--details] [--no-score]` emits the
       find_pairs.py-parity JSON (pairs: sequence, pair_category, num_hbonds,
       hbond_categories, classification, geometry, rigid_body, frames, score,
       issues; structure: overall/pairs/residues + score_details).
       (diff_json vs the real find_pairs CLI: per-pair content + scores exact
        within frame-precision tol; candidates_valid may differ ±1-2 at the
        validation cutoff — does not affect selected pairs/scores)

### Scope (lean)
Only the runtime tool: parse -> find -> classify -> score. Reference-data
generation (KDE/penalty fitting) stays in Python; the C++ scorer **loads** the
Python-produced `parameter_distributions.json`, `prosco_distributions.json`,
`penalty_weights.json` and applies them.
```
