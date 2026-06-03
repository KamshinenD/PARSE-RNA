/**
 * @file core_smoke_test.cpp
 * @brief Sanity checks for the core value types (no framework; returns 1 on fail).
 */
#include <cmath>
#include <iostream>

#include <pairfinder/core/residue.hpp>
#include <pairfinder/core/structure.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace {

int failures = 0;

void check(bool cond, const char* what) {
    if (!cond) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

}  // namespace

int main() {
    using namespace pairfinder;

    // geometry
    const geometry::Vector3d a{1, 0, 0};
    const geometry::Vector3d b{0, 2, 0};
    check(std::abs(a.distance_to(b) - std::sqrt(5.0)) < 1e-9, "distance_to");
    check(std::abs(a.cross(b).z - 2.0) < 1e-9, "cross");
    check(std::abs(a.dot(b)) < 1e-12, "dot orthogonal");

    // atom typing (prime + altname normalization)
    check(core::atom_type_from_name("O2'") == core::AtomType::O2_PRIME, "O2' type");
    check(core::atom_type_from_name("O2*") == core::AtomType::O2_PRIME, "O2* type");
    check(core::atom_type_from_name("O2P") == core::AtomType::OP2, "O2P type");
    check(core::atom_type_from_name("XYZ") == core::AtomType::UNKNOWN, "unknown type");

    // base-type helpers
    check(core::is_purine("G") && core::is_purine("PSU") == false, "is_purine G/PSU");
    check(core::is_pyrimidine("PSU"), "PSU is pyrimidine (normalized U)");
    check(core::normalize_base_type("DA") == "A", "DA -> A");
    check(core::normalize_base_type("XYZ") == "XYZ", "unknown kept");

    // residue
    core::Residue g("A-G-1", "G");
    g.add_atom("N9", {0, 0, 0});
    g.add_atom("C1'", {1.5, 0, 0});
    check(g.res_id() == "A-G-1", "res_id format");
    check(g.chain() == "A" && g.seq_num() == "1", "chain/seq from res_id");
    check(g.has_atom("N9") && !g.has_atom("N1"), "has_atom");
    g.add_atom("N9", {9, 9, 9});  // overwrite, not duplicate
    check(g.atoms().size() == 2 && g.get_atom("N9")->coords.x == 9.0, "add_atom overwrite");
    check(g.has_atom_type(core::AtomType::C1_PRIME), "has_atom_type C1'");
    check(g.glycosidic_n() != nullptr && g.glycosidic_n()->name == "N9", "glycosidic N9");

    // structure
    core::Structure s("test");
    core::Chain ch("A");
    ch.add_residue(std::move(g));
    s.add_chain(std::move(ch));
    check(s.residue_count() == 1, "residue_count");

    if (failures == 0) std::cout << "core_smoke_test: all passed\n";
    return failures == 0 ? 0 : 1;
}
