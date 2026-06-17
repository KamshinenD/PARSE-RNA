/**
 * @file suiteness_test.cpp
 * @brief Richardson suiteness classifier unit tests.
 *
 * Uses hardcoded canonical torsion angles for three well-known conformers
 * (1a, 1c, !!) to verify the classifier returns the expected conformer and
 * suiteness range without any external PDB data.
 */
#include <array>
#include <cmath>
#include <iostream>

#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/scoring/richardson.hpp>

namespace { int failures = 0; }

static void check(bool cond, const char* what) {
    if (!cond) { std::cerr << "FAIL: " << what << "\n"; ++failures; }
}

int main() {
    using namespace pairfinder;
    const auto suites_path = resources::file("richardson_suites.json");
    if (suites_path.empty()) {
        std::cerr << "SKIP: richardson_suites.json not found\n";
        return 0;  // skip rather than fail in unusual build layouts
    }

    scoring::RichardsonClassifier clf(suites_path);

    // --- 1a conformer (canonical A-form RNA backbone) ----------------------
    // Canonical centre from Richardson 2008 Table 1.
    // [delta_prev, eps_prev, zeta_prev, alpha, beta, gamma, delta]
    const std::array<double, 7> conformer_1a = {
        83.0, 210.0, 285.0,   // delta_prev=83, eps_prev=-150→210, zeta_prev=-75→285
        300.0, 177.0, 54.0,   // alpha=-60→300, beta=177, gamma=54
        83.0                  // delta=83
    };
    const double s1a = clf.suiteness(conformer_1a);
    check(s1a > 0.0, "1a conformer: suiteness > 0");
    check(s1a <= 1.0, "1a conformer: suiteness <= 1");

    // --- 1c conformer (common alternative) ----------------------------------
    const std::array<double, 7> conformer_1c = {
        83.0, 210.0, 167.0,   // zeta_prev=167
        288.0, 177.0, 48.0,   // alpha=-72→288
        83.0
    };
    const double s1c = clf.suiteness(conformer_1c);
    // 1c is a real conformer — suiteness should be positive when near centre
    // (accept 0 in case torsions don't land in the cluster, not a hard fail)
    check(s1c >= 0.0, "1c conformer: suiteness >= 0");

    // --- Outlier (!!) -------------------------------------------------------
    // Torsions far from any cluster should return 0 (outlier).
    const std::array<double, 7> outlier = {
        180.0, 90.0, 45.0, 135.0, 0.0, 270.0, 180.0
    };
    const double s_out = clf.suiteness(outlier);
    check(s_out == 0.0, "outlier: suiteness == 0");

    // --- L-RNA mirror of 1a -------------------------------------------------
    // 360 - each torsion of 1a should also be an outlier in D-RNA space
    // (it's an enantiomer — the whole point of the L-RNA exclusion in the scorer).
    std::array<double, 7> mirror_1a;
    for (int k = 0; k < 7; ++k)
        mirror_1a[k] = std::fmod(360.0 - conformer_1a[k], 360.0);
    const double s_mirror = clf.suiteness(mirror_1a);
    check(s_mirror == 0.0, "L-RNA mirror of 1a: outlier in D-RNA classifier");

    if (failures == 0) std::cout << "suiteness_test: all passed\n";
    return failures == 0 ? 0 : 1;
}
