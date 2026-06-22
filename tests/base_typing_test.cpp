/**
 * @file base_typing_test.cpp
 * @brief Unit tests for core::BaseTyping normalization and glycosidic-N lookup.
 */
#include <string>

#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/resource_locator.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::core::BaseTyping;
using pairfinder::tests::check;
using pairfinder::tests::finish;

}  // namespace

int main() {
    const auto db = pairfinder::resources::config("ligand_hbond_db.json");
    BaseTyping typing(db);

    check(typing.normalize("PSU") == "U", "PSU -> U");
    check(typing.normalize("DA") == "A", "DA -> A");
    check(typing.normalize("G") == "G", "G -> G");
    check(typing.normalize("XYZNOTBASE") == "XYZNOTBASE", "unknown kept");

    check(typing.is_purine("G"), "G purine");
    check(typing.is_purine("2MG"), "2MG purine");
    check(!typing.is_purine("U"), "U not purine");
    check(typing.is_pyrimidine("C"), "C pyrimidine");
    check(typing.is_pyrimidine("5MC"), "5MC pyrimidine");

    check(typing.glycosidic_n_name("G") == "N9", "G glycosidic N9");
    check(typing.glycosidic_n_name("C") == "N1", "C glycosidic N1");
    check(typing.glycosidic_n_name("PSU") == "N1", "PSU glycosidic N1");
    check(typing.glycosidic_n_name("UNKNOWN") == "", "unknown glycosidic empty");

    // Curated map works even with empty DB path.
    BaseTyping empty_typing("");
    check(empty_typing.normalize("PSU") == "U", "empty path: PSU -> U");
    check(empty_typing.glycosidic_n_name("A") == "N9", "empty path: A N9");

    return finish("base_typing_test");
}
