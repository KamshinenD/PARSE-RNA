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
/// Per-pair helix membership: 1 if the pair sits in a detected helix segment
/// (>=2 stacked consecutive pairs), else 0. Same index order as `pairs`. Mirror
/// of HelixDetector.detect + the length>=2 test used for the sugar-pucker flag.
std::vector<char> helix_membership(
    const std::vector<classification::ScoredCandidate>& pairs, const RNAChains& chains);

std::vector<classification::ScoredCandidate> select_pairs(
    std::vector<classification::ScoredCandidate> candidates, const RNAChains& chains,
    double min_score = 0.0);

}  // namespace pairfinder::algorithms::selection

#endif  // PAIRFINDER_ALGORITHMS_SELECTION_HPP
