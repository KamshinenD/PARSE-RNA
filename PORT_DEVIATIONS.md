# Port deviations & subtleties to revisit

Running log of every place the C++ port deviates from the Python source, makes a
deliberate engineering trade-off, or relies on a non-obvious faithfulness fix.
Each is differential-verified unless noted. Review before open-sourcing.

## Deliberate engineering trade-offs

1. **KD-tree → brute-force neighbor search** (Step 5, `candidate_finder.cpp`).
   Python `PairCache` uses `scipy.spatial.KDTree.query_ball_point(r=15)`. C++ uses
   an O(n²) brute-force `origin distance ≤ 15`. Identical metric/radius → identical
   candidate set (120 PDBs exact). NOTE: **the reference x3dna `base_pair_finder.cpp`
   also uses brute-force O(n²)** (squared-distance early-out) — so this matches the
   reference, not just the Python. PROFILED on the 80S ribosome (8ife.cif, 17K
   residues): `score_candidates` (find+classify) is only **0.20s** — the O(n²)
   candidate finder is NOT the bottleneck. So this stays brute-force; no spatial
   index needed here. (The real bottleneck was `rna_chains`; see below.)

   1b. **`rna_chains` spatial grid [2026-06-03, FIXED].** `RNAChains::from_structure`
   built backbone chains by scanning ALL residues to find each O3'↔P neighbor,
   repeated per residue added → O(n²) (~290M `connected()` calls on the ribosome).
   Profiling showed it at **2.13s of 2.94s (72%)** — the true ribosome bottleneck,
   not candidate finding. Replaced the linear scans with a uniform spatial hash
   (cell = 2.75Å backbone cutoff; `grid_p` on P/PA, `grid_o3` on O3'); each
   extension queries the 27 neighbor cells and distance-filters, picking the
   lowest residue index (identical tie-break to the old first-in-order scan).
   Result: **2.13s → 0.007s** (~300×), chains byte-identical → diff_pairs /
   diff_structure / diff_json still exact. **80S ribosome end-to-end: 2.94s →
   ~0.7s — under the <1s goal.** (Reference x3dna find_pair_app: 12s.)

2. **Bifurcation code intentionally omitted** (Step 4, `hbond/finder.cpp`).
   Python's `_try_select_bond`/`_find_alternative_slots`/`can_add_bond` form an
   orphan island — unreachable from `find_between`/`select_optimal` (which use the
   MIS selector). This is the PI's design (verified against the PI's original
   backup; not a local edit). **If `select_optimal` is ever rewired to allow a slot
   to serve two angularly-separated bonds, the bifurcation island must be ported.**
   x3dna's own `SlotOptimizer::select_optimal` DOES bifurcate — so x3dna and the
   current Python differ here; we follow Python.

## Faithfulness fixes (subtle, easy to get wrong)

3. **Pipeline base types are stripped (`DC→C`, `DT→T`, `DA→A`)** (Step 6).
   Production residues come from `PairCache._extract_residue_name`, NOT the parser's
   kept name ("DC"). This flips `is_rna_pyrimidine` (DC normalizes to "C" → treated
   as an RNA pyrimidine!) and changes DNA H-bond capacity keys. C++ strips base
   types via `core::pipeline_base_type(res_id)` **after** frame computation (frames
   need the original DA/PSU for template selection). Applied in `dump-lwclass` /
   `dump-classified`; the real pipeline (Step 9) must do the same.

4. **Three "normalize" notions — do not conflate** (Steps 3/5/6).
   - parser `base_type` = identity-ish upper/keep ("DC" stays "DC").
   - `core::NucleotideRegistry` (modified_nucleotides.json) → frame template selection.
   - `core::BaseTyping` (ligand_hbond_db.json parent + curated + DNA-prefix) →
     the Python *core* `normalize_base_type` used by classification/glycosidic.
   These come from different data sources and can disagree (e.g. PSU: registry
   canonical='P' for Atomic_P template, but BaseTyping normalizes PSU→U).

## Numerical / tolerance notes

5. **Computed frames (full precision) vs JSON-loaded frames (6-decimal) [2026-06-05].**
   C++ computes frames; Python production *loads* the ls_fitting JSON. As of
   2026-06-05 the C++ uses a **faithful port of x3dna's exact `LeastSquaresFitter`**
   (`geometry/superposition.cpp` — covariance → 4x4 quaternion key matrix → Jacobi
   largest-eigenvector → rotation), i.e. the *same method* that generated the JSON,
   at full double precision. So the residual difference is NOT the algorithm — it
   is that **the ls_fitting JSON stores frames rounded to 6 decimals** (e.g.
   `-0.567567`), while the C++ keeps full precision. PROOF: swapping the C++ from
   the old Horn impl to x3dna's exact fitter left `candidates_valid` unchanged
   (6V3A: 5303 both), still 44 above Python's 5259 — the gap is entirely the
   6-decimal JSON rounding flipping candidates that sit exactly on a validation
   cutoff. Benign: zero effect on selected pairs or scores. To close it one would
   bump the JSON precision in find_pair_2's `generate_modern_json` (not a C++
   change). Only `plane_angle` amplifies the per-value diff (see #6).

6. **`plane_angle` differential tolerance = 0.05°** (Step 5, `diff_candidates.py`).
   `plane_angle = deg(acos(|dir_z|))`; near `|dir_z|=1` (parallel bases) acos is
   hyper-sensitive, so the ~1e-6 frame difference amplifies to ~1e-3°. Provably
   benign: only near `plane_angle≈0`, far from the 70° validity cutoff, so
   `is_valid` always matches exactly. Carried into `quality_score` (the `/180`
   term) → `diff_candidates` uses `QTOL=3e-3`.

7. **x3dna quaternion fitter used for Kabsch RMSD** (Step 6b, `template_aligner.cpp`).
   Template RMSD uses `geometry::superpose().rms` (now x3dna's `LeastSquaresFitter`
   port, see #5) instead of Python's SVD Kabsch.
   Both return the optimal-alignment RMSD (global minimum is unique) → match to
   ~1e-9. Template RMSD is frame-independent (template + target coords only).

## Test-scope notes (not code deviations)

8. **`diff_hbonds` uses non-production base types.** The Step-4 H-bond unit test
   feeds Python `parse_pdb` residues (base_type "DC"), matching C++ parse_pdb — it
   validates the finder *logic*. Production H-bonds use stripped base types; that
   path is exercised by the Step-6 `diff_classified` face voters (stripped).

9. **Exemplar template dir effectively unused** (Step 6b). Files live in
   `basepair-catalog-exemplars/<lw>/...` but the loader does a *flat* lookup →
   returns None (matches Python exactly). If exemplars are reorganized flat, both
   implementations' behavior changes.

17. **altloc-B-only residues: C++ (raw PDB) scores them, Python (precomputed
    `pdb_atoms` JSON) doesn't.** Python's runtime scorer loads atoms from
    `json/pdb_atoms/<id>.json` (`cache._load_residues`), NOT the raw PDB. That
    preprocessing dropped residues whose every atom is alternate-conformer 'B'
    with no blank/'A' copy (e.g. 6C8M `A-G-15`, a chain-A terminal G that is
    altloc-B-only → absent from the JSON). C++ `dump-structure` parses the raw
    PDB and keeps it, so it scores one residue Python can't. This is a
    data-source divergence (same class as the ls_fitting-frame differentials),
    independent of #16 — a predecessor/parse change cannot add residue
    membership. Surfaced as the single broad-sweep diff (6C8M overall
    60.14/59.95). DECISION PENDING: leave as differential noise (C++ is the more
    complete tool), or have the scorer apply the same altloc-B-only drop for
    JSON parity.

18. **Step 9 `find` JSON — three faithful notes.**
    (a) **`rigid_body` block uses a THIRD parameter algorithm.** The CLI's
    per-pair `rigid_body` is `validation/rigid_body.py compute_rigid_body_parameters`
    (flip diag(1,-1,-1) on frame2, halve the relative rotation, express in the
    mid-frame), ported to `algorithms/rigid_body.{hpp,cpp}`. This is NOT the
    scorer's `scoring::compute_pair_parameters` (scorer_parameters/parameters.py,
    hinge/half-angle) — they give different numbers and are not interchangeable.
    (b) **`candidates_valid` may differ ±1-2 vs the Python CLI.** C++ computes
    frames from the raw PDB; the Python CLI loads ls_fitting JSON frames (~1e-6
    apart). A candidate sitting on a validation cutoff (d_v≤2.5 / plane≤70 /
    dorg≤15) can flip `is_valid` (e.g. 8GXB plane=69.667, d_v=2.4927). Same root
    as the diff_candidates 0.05° plane tol; it does NOT change the selected pairs
    or any score, so diff_json treats `candidates_valid` as a warning annotation
    (`candidates_total` is required exact).
    (c) **`template_rmsd: Infinity`.** Python `json.dump` emits a non-standard
    `Infinity` token for inf template_rmsd (--details, template-missing pairs);
    nlohmann would serialize inf as `null`, so the C++ emits the string
    `"Infinity"` to preserve the signal. diff_json treats both as equal.

19. **Parser swapped to GEMMI, replicating find_pair_2 conventions [2026-06-03].**
    The hand-rolled fixed-column PDB parser was replaced with a **GEMMI v0.6.5**
    loader (`io/pdb_parser.cpp`) — the same library find_pair_2 uses — so the port
    reads **PDB / mmCIF / .gz** and (critically) full ribosomes. It replicates
    find_pair_2's `convert_gemmi_structure` exactly: first model only; chain =
    gemmi chain name; seq = auth `seqid.num`; insertion code appended; alt-loc
    filter keep `' '/'A'/'1'`; HETATM kept only for modified nucleotides (via
    `NucleotideRegistry`), waters/ligands dropped; atom-name normalization
    `OP1→O1P, OP2→O2P, OP3→O3P, O1'→O4', OL→O1P, OR→O2P, C5A→C5M, O5T→O5', O3T→O3'`
    and `'*'→'\''`. CMake builds GEMMI **in-tree via FetchContent**, cloning
    v0.6.5 from upstream git (`EXCLUDE_FROM_ALL`/`SYSTEM`, static `gemmi_cpp`
    under our own `build/_deps`) — **zero dependency on find_pair_2** (no build
    artifact, no source seed). An `add_test()` shim drops GEMMI's unconditional
    `cpptest` from our CTest suite. Requires CMake ≥ 3.28.
    **Consequences vs the old parser:** (a) now matches find_pair_2/`pdb_atoms`
    exactly rather than Python `parse_pdb` (raw); (b) **resolves #17** — altloc-B-only
    residues are now dropped, matching the Python pipeline; (c) waters/ligands no
    longer appear as residues (harmless: they never got frames). Re-validated: all
    pipeline differentials (frames/candidates/lwclass/classified/scored/pairs/scores/
    structure/json) still exact. `diff_hbonds` oracle updated to drop waters + apply
    the same atom renames (it had used raw `parse_pdb`); now exact again (1EHZ 250
    H-bonds). New `diff_cif` proves `.cif` and `.pdb` of the same structure (1GID,
    1EHZ) yield byte-identical JSON. **Ribosome:** 8ife.cif (human 80S, 219,097
    atoms) runs end-to-end in **2.94s** (vs reference x3dna `find_pair_app` 12s),
    2157 pairs, overall 79.02.

## Step 6/7 gaps discovered while starting selection

**STATUS: #11-13 RESOLVED + verified (diff_scored vs real find_pairs all_candidates,
incl. all O2-affected PDBs exact). #14 (pre-selection filter) is part of Step 7.**

15. **Phosphate atom naming: raw PDB `OP1/OP2` vs pdb_atoms-JSON `O1P/O2P`
    [FIXED IN BOTH PYTHON + C++ — user-approved real fix].**
    generate_modern_json normalizes phosphate oxygens to O1P/O2P. The modified-base
    H-bond DB lists only OP-spelling, so a modified base's phosphate capacity
    lookup missed under O2P → production silently dropped modified-base phosphate
    H-bonds (e.g. 2MG-10 in 1EHZ: num_hbonds 1 should be 2).
    REAL FIX (both ends, identical): alias OP1<->O1P / OP2<->O2P / OP3<->O3P when
    extending capacities from the ligand DB — Python `hbond/geometry.py`
    `_PHOSPHATE_ALIAS` in `_extend_capacities_from_database`, C++
    `chemistry.cpp` `phos_alias` in the DB loop. (Standard bases already had both
    spellings; this mirrors that for DB-derived modified residues.) Now both FIND
    the bond and agree (verified via diff_scored vs real find_pairs).
    NOTE: this CHANGES production output (find_pairs now reports more modified-base
    phosphate H-bonds) — reference data / DSSR-comparison numbers / notebooks
    generated before this may need regeneration.
    Atom-name parity is kept separately by `core::apply_pipeline_view` renaming
    OP1->O1P etc. so C++ reports the same atom names as the pdb_atoms JSON.

16. **Backbone scoring dropped residues with negative seq numbers / insertion
    codes [FIXED IN BOTH PYTHON + C++ — user-approved real fix].**
    The suite predecessor lookup parsed the sequence number with
    `int(res_id.split("-")[2])`, which throws for an empty 3rd field (negative
    seqs: `"B-A--10"` → `parts[2]=""`) and for insertion / alt-conf codes
    (`"A-G-20A"` → `int("20A")`). Those residues silently got NO suiteness score
    (build_suite_for → None), and `compute_all_torsions`' sort_key lumped them at
    seq 0 (wrong prev/next → wrong torsions even when they were scored).
    REAL FIX (both ends, identical semantics): new `parse_res_seq` reconstructs
    the full sequence field as `"-".join(parts[2:])` and splits a signed-integer
    prefix from trailing insertion-code letters → `(chain, num:int, icode)`.
    Predecessor is now the residue immediately preceding within the chain ordered
    by `(num, icode)`, accepted when `prev.num == cur.num` or `prev.num == cur.num-1`
    (so `20 → 20A → 21` chains correctly; gaps ≥2 still break the suite, preserving
    chain-break behavior). Python: `identifiers.parse_res_seq`,
    `features.build_predecessor_index` (now maps `res_id → predecessor res_id`),
    `features.build_suite_for`, `torsions.sort_key`. C++: `parse_res_seq` +
    `ResSeq` in `scoring/torsions.{hpp,cpp}` (used by the torsion sort) and
    `build_predecessor_index` in `scoring/scorer.cpp`.
    NOTE: this CHANGES production output — structures with negative-seq /
    insertion-code residues now contribute those residues to the backbone
    (residues) score (e.g. 2VNU overall 0 → 27.29). Verified equal across both
    ends on the negative-seq set (2BX2, 2C0B, 2H0S, 2JEA, 2VNU) and many
    insertion-code structures (1HMH, 1MME, 299D, 2FMT, 1SER, 1HNW, 2CSX, …) plus
    a broad small-structure sweep (38/39 exact; the lone diff is #17, unrelated).

11. **Carbonyl-O2 correction is a production post-pass NOT yet in the C++ Step 6. [RESOLVED]**
    `finder._apply_o2_correction` (carbonyl_margin.py `committed_w_face_residues`
    + `correct_cis_o2`) runs after classify/score and can change a cis pyrimidine
    pair's class S→W (e.g. cWS→cWW), set `_precorr_lw`, and re-derive
    template_rmsd / quality_score / confidence. The C++ classifier and the
    `diff_classified` oracle BOTH omit it, so they agree with each other but
    NOT with the true production `lw_class` on O2-correctable pairs. Must port
    carbonyl_margin (2-pass over all candidates: build `committed` from pass-1
    classes, correct in pass-2) and feed corrected classes to selection. The
    selection `_filter_valid` grandfather branch reads `_precorr_lw`.

12. **`pair_category` (base-base vs other) gates classification in production.**
    `_score_candidates` only runs the LW pipeline for `pair_category=="base-base"`
    (from `_classify_hbond_category`/`_classify_pair_category` over the strict
    H-bonds; adjacent same-chain base-base → "adjacent-base-base"). The C++
    dump-lwclass/dump-classified classify ALL valid candidates. For the
    classification differentials that's fine (both sides classify all), but
    SELECTION needs non-base-base valid candidates to carry `lw_class=None`
    (they enter selection unclassified). Must port pair_category before selection.

13. **`has_strong_base_hbond` / `num_hbonds`** (set in `_score_candidates`):
    needed by selection `_class_valid` and the O2 `committed` set.
    has_strong_base_hbond = any base-base strict bond with distance ≤
    `max_hbond_distance` (default 3.5). num_hbonds = len(hbonds_list).

14. **Pre-selection filter**: valid candidates → `num_hbonds >= min_hbonds`
    (default 0) → `not is_stacking_signature` (per-class dNN/d_v thresholds,
    `validation/lw_class_filters.py`) → strategy.select_with_details.

Remaining Step 7 modules to port: carbonyl_margin, pair_category + stacking
filter, RNAChains (chains.py, spatial chain connectivity), `_CWWHelixPhase`
(Phase 1), `GlobalOptimalStrategy` (Phase 2 MWIS), `HelixPriorityStrategy`.
Done + validated so far: QualityScorer (selection quality_score, 90 PDBs exact
via diff_classified).

## Open-source TODO

10. **Hardcoded absolute resource paths [RESOLVED 2026-06-03].** All runtime data
    is now **vendored under `pair-finder-main-cpp/resources/`** (~3.6M: templates/,
    config/{modified_nucleotides,ligand_hbond_db}.json, basepair-idealized/,
    basepair-catalog-exemplars/, reference/{parameter_distributions,prosco_
    distributions,penalty_weights}.json, richardson_suites.json). A new
    `core::resource_locator` (`pairfinder::resources`, x3dna ResourceLocator
    pattern) resolves the root: `$PAIRFINDER_RESOURCES` → exe-relative
    (`<exe>/resources`, `<exe>/../resources`, `<exe>/../share/pairfinder/resources`)
    → cwd-relative → compile-time `PAIRFINDER_SOURCE_RESOURCES` (dev fallback).
    main.cpp path helpers + `io/pdb_parser`'s registry use it; per-resource
    `PAIRFINDER_*` env overrides retained. CMake stages `resources/` next to the
    binary (POST_BUILD) and `install()`s to `share/pairfinder/resources`. VERIFIED:
    `cmake --install` to an isolated prefix, run from elsewhere with ALL resource
    env vars unset → correct output (1GID 150 pairs / 72.78). No find_pair_2 /
    prototyped-pair-finder-main paths remain in the C++ (only the differential
    *test oracles* still reference them, which is correct — they ARE the oracle).
