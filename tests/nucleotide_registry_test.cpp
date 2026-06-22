/**
 * @file nucleotide_registry_test.cpp
 * @brief Unit tests for core::NucleotideRegistry (modified_nucleotides.json).
 */
#include <iostream>
#include <string>

#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/resource_locator.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::core::NucleotideRegistry;
using pairfinder::tests::check;
using pairfinder::tests::finish;

}  // namespace

int main() {
    const auto cfg = pairfinder::resources::config("modified_nucleotides.json");
    if (cfg.empty()) {
        std::cout << "SKIP: modified_nucleotides.json not found\n";
        return 0;
    }

    NucleotideRegistry reg(cfg);
    check(reg.size() > 0, "registry non-empty");

    {
        const auto info = reg.lookup("G");
        check(info.has_value(), "G lookup");
        if (info) {
            check(info->canonical == 'G', "G canonical");
            check(info->is_purine, "G purine");
            check(!info->is_modified, "G not modified");
        }
        check(reg.template_filename("G") == "Atomic_G.pdb", "G template");
    }

    {
        const auto info = reg.lookup("PSU");
        check(info.has_value(), "PSU lookup");
        if (info) {
            check(info->canonical == 'P', "PSU canonical P");
            check(!info->is_purine, "PSU not purine");
        }
        check(reg.template_filename("PSU") == "Atomic_P.pdb", "PSU template");
    }

    {
        const auto info = reg.lookup("5MC");
        check(info.has_value(), "5MC lookup");
        if (info) {
            check(info->canonical == 'C', "5MC canonical C");
            check(info->is_modified, "5MC modified");
        }
        check(reg.template_filename("5MC") == "Atomic.c.pdb", "5MC template");
    }

    check(!reg.lookup("NOTARESIDUE").has_value(), "unknown lookup nullopt");
    check(reg.template_filename("NOTARESIDUE").empty(), "unknown template empty");

    return finish("nucleotide_registry_test");
}
