/**
 * @file saenger_test.cpp
 * @brief Unit tests for classification::match_saenger.
 */
#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/hbond.hpp>
#include <pairfinder/algorithms/saenger.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/core/resource_locator.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::algorithms::classification::match_saenger;
using pairfinder::algorithms::hbond::HBond;
using pairfinder::core::BaseTyping;
using pairfinder::core::Residue;
using pairfinder::tests::check;
using pairfinder::tests::finish;

std::vector<HBond> au_cww_hbonds() {
    return {
        {"A-A-1", "A-U-2", "N6", "O4", 2.9, 0, 0, 0.8},
        {"A-U-2", "A-A-1", "N3", "N1", 2.8, 0, 0, 0.8},
    };
}

}  // namespace

int main() {
    const auto db = pairfinder::resources::config("ligand_hbond_db.json");
    BaseTyping typing(db);

    Residue a("A-A-1", "A");
    Residue u("A-U-2", "U");

    {
        const auto m = match_saenger(a, u, au_cww_hbonds(), typing);
        check(m.found, "A-U cWW: found");
        check(m.face1 == 'W', "A-U cWW: face1 W");
        check(m.face2 == 'W', "A-U cWW: face2 W");
    }

    {
        const auto m = match_saenger(a, u, {}, typing);
        check(!m.found, "empty hbonds: not found");
    }

    {
        Residue g("A-G-1", "G");
        Residue c("A-C-2", "C");
        const std::vector<HBond> wrong = {
            {"A-G-1", "A-C-2", "N2", "O2", 2.9, 0, 0, 0.8},
        };
        const auto m = match_saenger(g, c, wrong, typing);
        check(!m.found, "partial wrong pattern: not found");
    }

    {
        // Same atom pairs with donor/acceptor on res2 still match via reverse set.
        const std::vector<HBond> swapped = {
            {"A-U-2", "A-A-1", "O4", "N6", 2.9, 0, 0, 0.8},
            {"A-A-1", "A-U-2", "N1", "N3", 2.8, 0, 0, 0.8},
        };
        const auto m = match_saenger(a, u, swapped, typing);
        check(m.found, "swapped direction: found");
        check(m.face1 == 'W' && m.face2 == 'W', "swapped direction: cWW faces");
    }

    return finish("saenger_test");
}
