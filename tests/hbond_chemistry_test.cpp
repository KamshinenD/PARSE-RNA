/**
 * @file hbond_chemistry_test.cpp
 * @brief Unit tests for hbond geometry helpers and capacity tables.
 */
#include <string>
#include <unordered_map>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/geometry/vector3d.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::algorithms::hbond::HBondChemistry;
using pairfinder::algorithms::hbond::angle_between;
using pairfinder::algorithms::hbond::compute_base_normal;
using pairfinder::algorithms::hbond::hb_normalize;
using pairfinder::algorithms::hbond::rotate_vector;
using pairfinder::geometry::Vector3d;
using pairfinder::tests::check;
using pairfinder::tests::check_near;
using pairfinder::tests::finish;

}  // namespace

int main() {
    check_near(angle_between({1, 0, 0}, {0, 1, 0}), 90.0, 1e-9, "angle_between 90");

    {
        const Vector3d tiny{1e-12, 0, 0};
        const Vector3d n = hb_normalize(tiny);
        check_near(n.x, tiny.x, 1e-15, "hb_normalize tiny unchanged");
    }

    {
        const Vector3d v{1, 0, 0};
        const Vector3d out = rotate_vector(v, {0, 0, 1}, 90.0);
        check_near(out.x, 0.0, 1e-9, "rotate +X to +Y x");
        check_near(out.y, 1.0, 1e-9, "rotate +X to +Y y");
    }

    {
        std::unordered_map<std::string, Vector3d> ring;
        ring["C2"] = {0, 0, 0};
        ring["C4"] = {1, 0, 0};
        ring["C6"] = {0, 1, 0};
        const Vector3d normal = compute_base_normal(ring);
        check_near(normal.norm(), 1.0, 1e-9, "base_normal unit");
        check(normal.z > 0, "base_normal +z handedness");
    }

    const auto db = pairfinder::resources::config("ligand_hbond_db.json");
    HBondChemistry chem(db);

    {
        const auto cap = chem.donor_capacity("G", "N1");
        check(cap.has_value() && *cap > 0, "G N1 donor");
    }
    {
        const auto cap = chem.acceptor_capacity("U", "O4");
        check(cap.has_value() && *cap > 0, "U O4 acceptor");
    }
    {
        check(!chem.donor_capacity("G", "N99").has_value(), "G unknown atom donor");
        check(!chem.acceptor_capacity("U", "O99").has_value(), "U unknown atom acceptor");
    }

    return finish("hbond_chemistry_test");
}
