/**
 * @file richardson.hpp
 * @brief Richardson RNA backbone suite classifier (port of scoring/richardson.py).
 *
 * Loads richardson_suites.json (54 conformers + widths + triage/sieve +
 * satellite overrides) and classifies a 7D suite
 * [delta_prev, eps_prev, zeta_prev, alpha, beta, gamma, delta] into a suiteness
 * score in [0,1] (0 for triage/sieve/4D/7D outliers).
 */
#ifndef PAIRFINDER_SCORING_RICHARDSON_HPP
#define PAIRFINDER_SCORING_RICHARDSON_HPP

#include <array>
#include <filesystem>
#include <memory>

namespace pairfinder::scoring {

class RichardsonClassifier {
public:
    explicit RichardsonClassifier(const std::filesystem::path& suites_path);
    ~RichardsonClassifier();
    RichardsonClassifier(RichardsonClassifier&&) noexcept;

    /// Suiteness in [0,1] (0 for outliers). suite = 7 angles, each already %360.
    double suiteness(const std::array<double, 7>& suite) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_RICHARDSON_HPP
