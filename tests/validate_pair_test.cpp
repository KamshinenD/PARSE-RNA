/**
 * @file validate_pair_test.cpp
 * @brief Unit tests for algorithms::validate_pair geometry checks.
 */
#include <array>
#include <cmath>

#include <pairfinder/algorithms/candidate_finder.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/geometry/vector3d.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::algorithms::ValidationResult;
using pairfinder::algorithms::ValidationThresholds;
using pairfinder::algorithms::validate_pair;
using pairfinder::core::ReferenceFrame;
using pairfinder::geometry::Vector3d;
using pairfinder::tests::check;
using pairfinder::tests::check_near;
using pairfinder::tests::finish;

ReferenceFrame identity_frame(Vector3d origin) {
    ReferenceFrame f;
    f.origin = origin;
    f.rotation = std::array<std::array<double, 3>, 3>{
        std::array<double, 3>{1, 0, 0},
        std::array<double, 3>{0, 1, 0},
        std::array<double, 3>{0, 0, 1},
    };
    f.valid = true;
    return f;
}

ReferenceFrame tilted_z_frame(Vector3d origin, double tilt_deg) {
    const double rad = tilt_deg * 3.14159265358979323846 / 180.0;
    const double c = std::cos(rad);
    const double s = std::sin(rad);
    ReferenceFrame f;
    f.origin = origin;
    // Rotate +Z toward +Y by tilt_deg about the x-axis.
    f.rotation = std::array<std::array<double, 3>, 3>{
        std::array<double, 3>{1, 0, 0},
        std::array<double, 3>{0, c, -s},
        std::array<double, 3>{0, s, c},
    };
    f.valid = true;
    return f;
}

}  // namespace

int main() {
    const ValidationThresholds thr;

    // Canonical parallel geometry: all checks pass.
    {
        const auto f1 = identity_frame({0, 0, 0});
        const auto f2 = identity_frame({10, 0, 0});
        const Vector3d n1{0, 0, 0};
        const Vector3d n2{10, 0, 0};
        const ValidationResult r = validate_pair(f1, f2, n1, n2, thr);
        check_near(r.dorg, 10.0, 1e-9, "valid: dorg");
        check_near(r.d_v, 0.0, 1e-9, "valid: d_v");
        check_near(r.plane_angle, 0.0, 1e-9, "valid: plane_angle");
        check_near(r.dNN, 10.0, 1e-9, "valid: dNN");
        check(r.distance_check, "valid: distance_check");
        check(r.d_v_check, "valid: d_v_check");
        check(r.plane_angle_check, "valid: plane_angle_check");
        check(r.dNN_check, "valid: dNN_check");
        check(r.is_valid, "valid: is_valid");
        const double expected_q = r.dorg + thr.d_v_weight * r.d_v +
                                  r.plane_angle / thr.plane_angle_divisor;
        check_near(r.quality_score, expected_q, 1e-9, "valid: quality_score");
    }

    // Origin distance too large.
    {
        const auto f1 = identity_frame({0, 0, 0});
        const auto f2 = identity_frame({20, 0, 0});
        const ValidationResult r = validate_pair(f1, f2, {0, 0, 0}, {20, 0, 0}, thr);
        check(!r.distance_check, "dorg_fail: distance_check");
        check(!r.is_valid, "dorg_fail: is_valid");
    }

    // Large vertical separation along mean helix axis.
    {
        const auto f1 = identity_frame({0, 0, 0});
        const auto f2 = identity_frame({0, 0, 5});
        const ValidationResult r = validate_pair(f1, f2, {0, 0, 0}, {0, 0, 5}, thr);
        check_near(r.d_v, 5.0, 1e-9, "dv_fail: d_v");
        check(!r.d_v_check, "dv_fail: d_v_check");
        check(!r.is_valid, "dv_fail: is_valid");
    }

    // Glycosidic N too close.
    {
        const auto f1 = identity_frame({0, 0, 0});
        const auto f2 = identity_frame({10, 0, 0});
        const ValidationResult r = validate_pair(f1, f2, {0, 0, 0}, {2, 0, 0}, thr);
        check_near(r.dNN, 2.0, 1e-9, "dnn_fail: dNN");
        check(!r.dNN_check, "dnn_fail: dNN_check");
        check(!r.is_valid, "dnn_fail: is_valid");
    }

    // Tilted z-axes (45°) -> plane_angle ~ 45°, still within default threshold.
    {
        const auto f1 = identity_frame({0, 0, 0});
        const auto f2 = tilted_z_frame({10, 0, 0}, 45.0);
        const ValidationResult r = validate_pair(f1, f2, {0, 0, 0}, {10, 0, 0}, thr);
        check_near(r.plane_angle, 45.0, 1e-6, "tilt45: plane_angle");
        check(r.plane_angle_check, "tilt45: plane_angle_check");
        check(r.is_valid, "tilt45: is_valid");
    }

    // Large plane tilt fails validation (> 70° default).
    {
        const auto f1 = identity_frame({0, 0, 0});
        const auto f2 = tilted_z_frame({10, 0, 0}, 80.0);
        const ValidationResult r = validate_pair(f1, f2, {0, 0, 0}, {10, 0, 0}, thr);
        check(r.plane_angle > thr.max_plane_angle, "tilt80: plane_angle over threshold");
        check(!r.plane_angle_check, "tilt80: plane_angle_check");
        check(!r.is_valid, "tilt80: is_valid");
    }

    return finish("validate_pair_test");
}
