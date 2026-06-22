/**
 * @file selection_planarity_test.cpp
 * @brief Shared-residue planarity selection and H-bond count fallback (1GID regression).
 */
#include <string>
#include <unordered_set>
#include <vector>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/quality_scorer.hpp>
#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/selection.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/core/structure.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::algorithms::RNAChains;
using pairfinder::algorithms::classification::HBondPatterns;
using pairfinder::algorithms::classification::QualityScorer;
using pairfinder::algorithms::classification::ScoredCandidate;
using pairfinder::algorithms::classification::ScoringHBond;
using pairfinder::algorithms::ValidationResult;
using pairfinder::algorithms::hbond::HBondChemistry;
using pairfinder::algorithms::selection::select_pairs;
using pairfinder::core::BaseTyping;
using pairfinder::core::Structure;
using pairfinder::tests::check;
using pairfinder::tests::check_near;
using pairfinder::tests::finish;

ScoredCandidate make_candidate(const std::string& res_id1, const std::string& res_id2,
                               const std::string& lw_class, double quality_score,
                               double plane_angle, double d_v, char base1 = 'A',
                               char base2 = 'A') {
    ScoredCandidate c;
    c.res_id1 = res_id1;
    c.res_id2 = res_id2;
    c.res_name1 = std::string(1, base1);
    c.res_name2 = std::string(1, base2);
    c.lw_class = lw_class;
    c.quality_score = quality_score;
    c.has_strong_base_hbond = true;
    c.template_rmsd = 1.5;
    c.validation.is_valid = true;
    c.validation.plane_angle = plane_angle;
    c.validation.d_v = d_v;
    c.validation.dNN = 9.0;
    c.validation.quality_score = quality_score * 10.0;
    return c;
}

bool pair_set_has(const std::vector<ScoredCandidate>& selected, const std::string& a,
                  const std::string& b) {
    for (const auto& p : selected) {
        if ((p.res_id1 == a && p.res_id2 == b) || (p.res_id1 == b && p.res_id2 == a))
            return true;
    }
    return false;
}

}  // namespace

int main() {
    const RNAChains chains = RNAChains::from_structure(Structure{});

    // 1GID-style: A114·A206 beats A115·A206 on planarity despite lower quality.
    {
        std::vector<ScoredCandidate> cands;
        cands.push_back(make_candidate("B-A-114", "B-A-206", "tHS", 0.666, 7.748, 0.739));
        cands.push_back(make_candidate("B-A-115", "B-A-206", "cWS", 0.754, 12.649, 1.749));
        const auto selected = select_pairs(std::move(cands), chains, 0.0);
        check(selected.size() == 1, "shared-residue clique picks exactly one pair");
        check(pair_set_has(selected, "B-A-114", "B-A-206"),
              "more coplanar A114·A206 selected over A115·A206");
        check(!pair_set_has(selected, "B-A-115", "B-A-206"),
              "A115·A206 not selected when A114·A206 is more coplanar");
    }

    // H-bond count fallback: unknown non-cWW class expects 1 bond, not 2.
    {
        const auto db = pairfinder::resources::config("ligand_hbond_db.json");
        const BaseTyping typing(db);
        const HBondChemistry chem(db);
        const HBondPatterns patterns(pairfinder::resources::basepair_idealized(), chem);
        const QualityScorer scorer;
        ValidationResult vr;
        vr.is_valid = true;

        ScoringHBond hb;
        hb.donor_atom = "N6";
        hb.acceptor_atom = "N3";
        hb.distance = 2.97;
        hb.alignment = 1.0;
        hb.base_base = true;

        const double one_bond =
            scorer.compute_score(vr, "AC", {hb}, 2.0, "tZZ", patterns, typing);
        ScoringHBond hb2 = hb;
        hb2.donor_atom = "N1";
        const double two_bonds =
            scorer.compute_score(vr, "AC", {hb, hb2}, 2.0, "tZZ", patterns, typing);
        check_near(one_bond, two_bonds, 0.02,
                   "unknown non-cWW class: one bond already satisfies expected count");
    }

    return finish("selection_planarity_test");
}
