/**
 * @file test_helpers.hpp
 * @brief Minimal helpers for lightweight CTest executables (no framework).
 */
#ifndef PAIRFINDER_TESTS_TEST_HELPERS_HPP
#define PAIRFINDER_TESTS_TEST_HELPERS_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace pairfinder::tests {

inline int& failures() {
    static int n = 0;
    return n;
}

inline void check(bool cond, const char* what) {
    if (!cond) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures();
    }
}

inline void check_near(double got, double expected, double tol, const char* what) {
    if (std::abs(got - expected) > tol) {
        std::cerr << "FAIL: " << what << " (got " << got << ", expected " << expected
                  << ", tol " << tol << ")\n";
        ++failures();
    }
}

inline int finish(const char* name) {
    if (failures() == 0) {
        std::cout << name << ": all passed\n";
        return EXIT_SUCCESS;
    }
    std::cerr << name << ": " << failures() << " failure(s)\n";
    return EXIT_FAILURE;
}

}  // namespace pairfinder::tests

#endif  // PAIRFINDER_TESTS_TEST_HELPERS_HPP
