# Lightweight tests with no external framework (CTest-driven).
# Each test is a small executable that returns non-zero on failure.

set(_pf_test_inc ${CMAKE_CURRENT_SOURCE_DIR})

macro(pf_add_unit_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE pairfinder_core)
    target_include_directories(${name} PRIVATE ${_pf_test_inc})
    add_test(NAME ${name} COMMAND ${name})
endmacro()

pf_add_unit_test(core_smoke_test)
pf_add_unit_test(suiteness_test)
pf_add_unit_test(superposition_test)
pf_add_unit_test(validate_pair_test)
pf_add_unit_test(parse_res_seq_test)
pf_add_unit_test(base_typing_test)
pf_add_unit_test(nucleotide_registry_test)
pf_add_unit_test(hbond_chemistry_test)
pf_add_unit_test(saenger_test)
pf_add_unit_test(selection_planarity_test)

# Differential tests vs the Python / x3dna oracles. These need the find_pair_2
# reference data; registered only when present so a bare checkout still builds.
find_package(Python3 COMPONENTS Interpreter QUIET)
set(_fp2_pdb "/Users/kdewan2/Desktop/Projects/find_pair_2/data/pdb")
if(Python3_FOUND AND EXISTS "${_fp2_pdb}")
    # Frame fit vs ls_fitting: covers modified bases (1EHZ), the RMSD type-check
    # gate (4E8M/EPE), and altloc/occupancy handling (1M8W).
    add_test(NAME diff_frames
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_frames.py
            $<TARGET_FILE:pairfinder>
            ${_fp2_pdb}/1EHZ.pdb ${_fp2_pdb}/100D.pdb ${_fp2_pdb}/2CKY.pdb
            ${_fp2_pdb}/4E8M.pdb ${_fp2_pdb}/1M8W.pdb)

    # CIF<->PDB parity: same structure from .cif and .pdb must give identical
    # JSON (proves the GEMMI loader's CIF path matches the validated PDB path).
    set(_fp2_root "/Users/kdewan2/Desktop/Projects/find_pair_2")
    if(EXISTS "${_fp2_root}/1gid.cif" AND EXISTS "${_fp2_root}/1ehz.cif")
        add_test(NAME diff_cif
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_cif.py
                $<TARGET_FILE:pairfinder>
                ${_fp2_pdb}/1GID.pdb ${_fp2_root}/1gid.cif
                ${_fp2_pdb}/1EHZ.pdb ${_fp2_root}/1ehz.cif)
    endif()

    # Slot-based H-bond finder vs the Python HBondFinder: DNA (100D) + a tRNA
    # with many modified bases (1EHZ). Needs the Python package importable.
    add_test(NAME diff_hbonds
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_hbonds.py
            $<TARGET_FILE:pairfinder>
            ${_fp2_pdb}/100D.pdb ${_fp2_pdb}/1EHZ.pdb)

    # Candidate finding + geometric validation vs the Python PairCache (frames
    # computed in C++; oracle uses the ls_fitting + pdb_atoms JSON).
    add_test(NAME diff_candidates
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_candidates.py
            $<TARGET_FILE:pairfinder>
            ${_fp2_pdb}/100D.pdb ${_fp2_pdb}/1EHZ.pdb)

    # LW edge classification vs Python classify_lw_class (geometry + H-bond face
    # voting, aligner=None production path).
    add_test(NAME diff_lwclass
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_lwclass.py
            $<TARGET_FILE:pairfinder>
            ${_fp2_pdb}/100D.pdb ${_fp2_pdb}/1EHZ.pdb)

    # Full classification (geometry + template RMSD + ambiguity resolution) vs
    # the Python finder classify block. Needs the basepair template dirs.
    if(EXISTS "/Users/kdewan2/Desktop/Projects/find_pair_2/basepair-idealized")
        add_test(NAME diff_classified
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_classified.py
                $<TARGET_FILE:pairfinder>
                ${_fp2_pdb}/100D.pdb ${_fp2_pdb}/1EHZ.pdb)
    endif()

    # End-to-end tests vs the REAL find_pairs (need a scipy-capable Python).
    # PAIRFINDER_SCIPY_PYTHON overrides; defaults to the o2_ablation venv.
    set(_scipy_py "$ENV{PAIRFINDER_SCIPY_PYTHON}")
    if(NOT _scipy_py)
        set(_scipy_py "/Users/kdewan2/Desktop/Projects/find-pair-score-rna/tmp/o2_ablation/venv/bin/python")
    endif()
    if(EXISTS "${_scipy_py}")
        # Full scoring pass (classify + score + O2 correction) vs find_pairs all_candidates.
        add_test(NAME diff_scored
            COMMAND ${_scipy_py}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_scored.py
                $<TARGET_FILE:pairfinder> 1EHZ 100D 8GXB 6DN1 8CAZ)
        # Final selected pairs vs find_pairs().pairs (whole pipeline).
        add_test(NAME diff_pairs
            COMMAND ${_scipy_py}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_pairs.py
                $<TARGET_FILE:pairfinder> 1EHZ 100D 8GXB 6DN1 8CAZ 1ASY)
        # Empirical pair scores vs the real Scorer.
        add_test(NAME diff_scores
            COMMAND ${_scipy_py}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_scores.py
                $<TARGET_FILE:pairfinder> 1EHZ 100D 8GXB 8CAZ 1ASY)
        # Whole-structure score (pairs + residue suiteness + overall) vs Scorer.
        add_test(NAME diff_structure
            COMMAND ${_scipy_py}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_structure.py
                $<TARGET_FILE:pairfinder> 1EHZ 100D 8GXB 6DN1 8CAZ 1ASY)
        # Step 9: full `find` JSON (pairs + classification + geometry + rigid_body
        # + frames + scores + issues + score_details) vs the find_pairs.py CLI.
        # The broad set: DNA, tRNA (modified bases), group-I intron, ambiguous /
        # O2-correction / negative-seq / insertion-code hard cases, mid-size RNAs,
        # and three rRNA-rich ribosomes (~1300 pairs each). ~60s under the venv.
        add_test(NAME diff_json
            COMMAND ${_scipy_py}
                ${CMAKE_CURRENT_SOURCE_DIR}/differential/diff_json.py
                $<TARGET_FILE:pairfinder>
                100D 1EHZ 8GXB 6DN1 8CAZ 1ASY 2VNU 1HMH 299D 1GID
                3SUY 7T3K 4MEH 6VMY 9BZ1 7V9X 7C7A
                8VRL 6DZP 7AQD)
        set_tests_properties(diff_json PROPERTIES TIMEOUT 600)
    endif()
endif()
