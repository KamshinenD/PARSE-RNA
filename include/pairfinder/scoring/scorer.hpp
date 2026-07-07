/**
 * @file scorer.hpp
 * @brief Empirical 0-100 base-pair quality scorer (port of scoring/scorer.py +
 *        issues.py + features.py + canonical.py). Pair scoring (Step 8a):
 *        score = 100 - sum(pair_weight[issue] * severity) over the 6 issues,
 *        severities driven by ProSco density tables + graded H-bond count.
 *
 * Loads reference data (parameter_distributions.json + prosco_distributions.json
 * + penalty_weights.json) — Python generates them; the C++ only consumes them.
 */
#ifndef PAIRFINDER_SCORING_SCORER_HPP
#define PAIRFINDER_SCORING_SCORER_HPP

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/core/structure.hpp>

namespace pairfinder::scoring {

namespace classification = pairfinder::algorithms::classification;
namespace hbond = pairfinder::algorithms::hbond;

/// One triggered issue and its weighted penalty contribution (weight × severity).
/// The Černý three-tier provenance (ProSco / Z' / tier) is carried for reporting
/// and downstream highlighting (e.g. PyMOL); it does not affect the score.
struct IssuePenalty {
    std::string issue;
    double weight = 0.0;            ///< pair_weight[issue] × severity
    double severity = 0.0;          ///< Černý severity in [0,1]
    std::optional<double> prosco;   ///< ProSco of the worst contributing value (none if N/A)
    std::optional<double> zprime;   ///< Z' of the worst contributing value (none if N/A)
    std::string tier;               ///< "Preferred"/"Allowed"/"Of Concern"; "" if not ProSco-scored
};

/// Per-pair scoring result.
struct PairScore {
    std::string res_id1, res_id2;
    std::string bp_type, lw_class;
    double score = 0.0;
    double penalty = 0.0;
    std::vector<IssuePenalty> issues;  ///< triggered issues, sorted by descending penalty
};

/// One suite torsion flagged as genuinely off its target conformer. Reported
/// only when it *fires*: ProSco below the Preferred cutoff (5) against the
/// empirical per-conformer distribution. No Z'/penalty — suiteness is the score.
struct AngleDeviation {
    std::string angle;              ///< suite dihedral name (alpha .. delta_prev)
    double value = 0.0;             ///< measured angle, degrees
    double target = 0.0;            ///< conformer-center angle, degrees
    double gap = 0.0;               ///< signed measured − target, degrees
    std::optional<double> prosco;   ///< density percentile (lower = more anomalous)
};

/// Refinement-facing backbone recommendation for one residue (mirror of Python
/// BackboneRecommendation). tier: "fixable" (≥1 firing torsion, correct these) /
/// "review" (low suiteness but no firing torsion — jointly distorted, rebuild) /
/// "flag_only" (outlier, nearest conformer is only a hint).
struct BackboneRecommendation {
    std::string tier;
    std::string target_conformer;
    std::vector<AngleDeviation> deviations;
};

/// Per-residue backbone score (suiteness * 100).
struct ResidueScore {
    std::string res_id;
    double suiteness = 0.0;
    double score = 0.0;
    std::optional<BackboneRecommendation> recommendation;  ///< none for ok residues
};

/// Whole-structure score.
struct StructureScore {
    double overall = 0.0;
    double pairs_score = 0.0;
    double residues_score = 0.0;
    int skipped_pairs = 0;  ///< selected pairs that were not scorable
    std::vector<PairScore> pair_scores;
    std::vector<ResidueScore> residue_scores;
};

class Scorer {
public:
    Scorer(const std::filesystem::path& distributions_path,
           const std::filesystem::path& prosco_path,
           const std::filesystem::path& weights_path,
           const std::filesystem::path& richardson_suites_path);
    ~Scorer();
    Scorer(Scorer&&) noexcept;

    /// _is_scorable: base-base, not ambiguous, single (no '|') lw_class.
    static bool is_scorable(const classification::ScoredCandidate& pair);

    /// Score one scorable pair (Step 8a). res1/res2 in the pair's residue order.
    PairScore score_pair(const classification::ScoredCandidate& pair,
                         const core::Residue& res1, const core::Residue& res2,
                         hbond::HBondFinder& finder, const hbond::HBondChemistry& chem,
                         const core::BaseTyping& typing) const;

    /// Whole-structure score: pair scores + residue suiteness + overall.
    StructureScore score_structure(const std::vector<classification::ScoredCandidate>& selected,
                                   const core::Structure& structure, hbond::HBondFinder& finder,
                                   const hbond::HBondChemistry& chem,
                                   const core::BaseTyping& typing) const;

    double w_pairs() const;
    double w_residues() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_SCORER_HPP
