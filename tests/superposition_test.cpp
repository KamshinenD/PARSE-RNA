/**
 * @file superposition_test.cpp
 * @brief Unit tests for geometry::superpose (Kabsch/Horn).
 */
#include <array>
#include <cmath>
#include <vector>

#include <pairfinder/geometry/superposition.hpp>
#include <pairfinder/geometry/vector3d.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::geometry::SuperpositionResult;
using pairfinder::geometry::Vector3d;
using pairfinder::tests::check;
using pairfinder::tests::check_near;
using pairfinder::tests::finish;

SuperpositionResult run(const std::vector<Vector3d>& from, const std::vector<Vector3d>& to) {
    return pairfinder::geometry::superpose(from, to);
}

bool near_matrix(const std::array<std::array<double, 3>, 3>& got,
                 const std::array<std::array<double, 3>, 3>& expected, double tol) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (std::abs(got[i][j] - expected[i][j]) > tol) return false;
    return true;
}

}  // namespace

int main() {
    // Identity fit.
    {
        const std::vector<Vector3d> pts = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const auto r = run(pts, pts);
        check(r.valid, "identity: valid");
        check_near(r.rms, 0.0, 1e-9, "identity: rms");
        check_near(r.translation.x, 0.0, 1e-9, "identity: tx");
        check_near(r.translation.y, 0.0, 1e-9, "identity: ty");
        check_near(r.translation.z, 0.0, 1e-9, "identity: tz");
        const std::array<std::array<double, 3>, 3> eye = {{
            {{1, 0, 0}},
            {{0, 1, 0}},
            {{0, 0, 1}},
        }};
        check(near_matrix(r.rotation, eye, 1e-9), "identity: rotation");
    }

    // Pure translation.
    {
        const std::vector<Vector3d> from = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        const std::vector<Vector3d> to = {{5, 0, 0}, {6, 0, 0}, {5, 1, 0}};
        const auto r = run(from, to);
        check(r.valid, "translation: valid");
        check_near(r.rms, 0.0, 1e-9, "translation: rms");
        check_near(r.translation.x, 5.0, 1e-9, "translation: tx");
        check_near(r.translation.y, 0.0, 1e-9, "translation: ty");
        check_near(r.translation.z, 0.0, 1e-9, "translation: tz");
    }

    // 90° rotation about Z maps +X onto +Y.
    {
        const std::vector<Vector3d> from = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const std::vector<Vector3d> to = {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}};
        const auto r = run(from, to);
        check(r.valid, "rot90z: valid");
        check_near(r.rms, 0.0, 1e-8, "rot90z: rms");
        // Apply R to from[0]: should land on to[0].
        const Vector3d mapped{
            r.rotation[0][0] * from[0].x + r.rotation[0][1] * from[0].y +
                r.rotation[0][2] * from[0].z + r.translation.x,
            r.rotation[1][0] * from[0].x + r.rotation[1][1] * from[0].y +
                r.rotation[1][2] * from[0].z + r.translation.y,
            r.rotation[2][0] * from[0].x + r.rotation[2][1] * from[0].y +
                r.rotation[2][2] * from[0].z + r.translation.z,
        };
        check_near(mapped.x, to[0].x, 1e-8, "rot90z: mapped x");
        check_near(mapped.y, to[0].y, 1e-8, "rot90z: mapped y");
        check_near(mapped.z, to[0].z, 1e-8, "rot90z: mapped z");
    }

    // Too few points.
    {
        const std::vector<Vector3d> from = {{0, 0, 0}, {1, 0, 0}};
        const std::vector<Vector3d> to = {{0, 0, 0}, {1, 0, 0}};
        const auto r = run(from, to);
        check(!r.valid, "too_few: invalid");
    }

    // Size mismatch.
    {
        const std::vector<Vector3d> from = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        const std::vector<Vector3d> to = {{0, 0, 0}, {1, 0, 0}};
        const auto r = run(from, to);
        check(!r.valid, "size_mismatch: invalid");
    }

    return finish("superposition_test");
}
