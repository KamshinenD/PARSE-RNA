/**
 * @file selection.hpp
 * @brief Two-phase pair selection (port of finder/selection.py HelixPriorityStrategy).
 *
 * Phase 1 (_CWWHelixPhase): cWW canonical/wobble pairs via helix-priority
 * (helix continuity > raw score; isolated pairs stricter). Phase 2
 * (GlobalOptimalStrategy): slot-aware maximum-weighted independent set over the
 * remaining candidates. Face-aware conflicts (a residue may pair on W and H
 * faces simultaneously). Pre-filtered by the stacking-signature filter.
 */
#ifndef PAIRFINDER_ALGORITHMS_SELECTION_HPP
#define PAIRFINDER_ALGORITHMS_SELECTION_HPP

#include <vector>

#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>

namespace pairfinder::algorithms::selection {

/// Final selected pairs (subset of the scored candidates, with selection-time
/// ambiguity resolution applied to lw_class). min_score from FinderConfig.
std::vector<classification::ScoredCandidate> select_pairs(
    std::vector<classification::ScoredCandidate> candidates, const RNAChains& chains,
    double min_score = 0.0);

}  // namespace pairfinder::algorithms::selection

#endif  // PAIRFINDER_ALGORITHMS_SELECTION_HPP
