/**
 * @file scorer.cpp
 * @brief Empirical pair scorer (port of scorer.py + issues.py + features.py + canonical.py).
 *
 * Severity uses the Cerny method (Cerny et al. NAR 2026 gkaf1335):
 *   severity = 0              if ProSco >= 5  (Preferred tier)
 *   severity = min(1, |Z'|/5) if ProSco <  5
 * Z' is an asymmetric non-parametric standard score loaded from z_tables.json.
 * The 9 individual parameters are scored independently (shear, stretch, stagger,
 * buckle, propeller, opening, distance, hbond_angles, incorrect_hbond_count).
 */
#include <pairfinder/scoring/scorer.hpp>

#include <pairfinder/scoring/clashes.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <nlohmann/json.hpp>

#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/selection.hpp>
#include <pairfinder/scoring/hbond_angles.hpp>
#include <pairfinder/scoring/pair_parameters.hpp>
#include <pairfinder/scoring/richardson.hpp>
#include <pairfinder/scoring/torsions.hpp>

namespace pairfinder::scoring {

namespace {

using json = nlohmann::json;
using geometry::Vector3d;

// ---- Sugar-pucker checks (Richardson/Kapral): the base-phosphate perpendicular
//      (Pperp) and the two report-only categories built on it. ----
constexpr double kPperpCutoff = 2.9;   // Å: >= C3'-endo, < C2'-endo
const std::unordered_set<std::string> kHelicalWcBp = {"G-C", "C-G", "A-U", "U-A"};

/// Base-phosphate perpendicular distance (Å): perpendicular from the 3' P (the
/// next residue's phosphorus) to the extended C1'->N glycosidic vector. Long
/// (>=2.9) = C3'-endo, short (<2.9) = C2'-endo. nullopt if any atom is missing.
std::optional<double> compute_pperp(const core::Residue* residue,
                                    const core::Residue* next_res) {
    if (!residue || !next_res) return std::nullopt;
    const core::Atom* c1 = residue->get_atom("C1'");
    const core::Atom* n  = residue->glycosidic_n();
    const core::Atom* p  = next_res->get_atom("P");
    if (!c1 || !n || !p) return std::nullopt;
    const Vector3d a = c1->coords;
    const Vector3d axis = n->coords - a;
    const double norm = axis.norm();
    if (norm < 1e-6) return std::nullopt;
    const Vector3d d = axis * (1.0 / norm);
    const Vector3d v = p->coords - a;
    const Vector3d perp = v - d * v.dot(d);
    return perp.norm();
}

int pucker_from_pperp(double pperp) { return pperp >= kPperpCutoff ? 3 : 2; }
const char* pucker_label(int p) { return p == 3 ? "C3'-endo" : "C2'-endo"; }

// Cerny severity thresholds.
constexpr double kProscoPreferred  = 5.0;   // ProSco >= 5 -> Preferred, severity = 0
constexpr double kZprimeThreshold  = 5.0;   // |Z'| = 5 -> severity = 1

// Backbone recommendation thresholds (mirror backbone_scorer.py).
constexpr double kLowSuitenessMax     = 0.5;  // suiteness below -> residue flagged
constexpr double kFallbackWidthFactor = 1.0;  // no ProSco cell -> fires beyond 1 width

const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
const std::unordered_set<std::string> kSelfReciprocal = {
    "cWW", "tWW", "cHH", "tHH", "cSS", "tSS"};
const std::unordered_set<std::string> kSignFlip = {"shear", "stagger", "buckle"};

// H-bond distance policy (hbond_distance_policy.py). Watson-Crick cells judge
// every base-base bond; all other cells judge only the strongest (primary) bond.
const std::unordered_set<std::string> kWcDistanceCells = {"C-G/cWW", "A-U/cWW"};
constexpr double kMinPhysicalDistance = 2.2;  // below this: not a real H-bond
constexpr double kMaxStrongDistance   = 3.5;  // strong (primary) bond cutoff

// 10 individual parameters — mirrors Python ISSUE_KEYS.
const std::array<const char*, 10> kIssueKeys = {
    "shear", "stretch", "stagger", "buckle", "propeller", "opening",
    "distance", "hbond_angles", "incorrect_hbond_count", "clash"};

// canonicalize_bp -> (canonical_bp, swapped).
std::pair<std::string, bool> canonicalize_bp(const std::string& bp, const std::string& lw) {
    if (!kSelfReciprocal.count(lw)) return {bp, false};
    const auto dash = bp.find('-');
    if (dash == std::string::npos) return {bp, false};
    const std::string left = bp.substr(0, dash), right = bp.substr(dash + 1);
    if (left <= right) return {bp, false};
    return {right + "-" + left, true};
}

double maybe_flip(double v, const std::string& source, const std::string& param, bool swapped) {
    if (!swapped) return v;
    if (source == "pairs" && kSignFlip.count(param)) return -v;
    return v;
}

// True if (canonical bp_type, lw_class) is a Watson-Crick distance cell.
bool is_wc_distance_cell(const std::string& bp, const std::string& lw) {
    const auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
    (void)swapped;
    return kWcDistanceCells.count(canon_bp + "/" + lw) > 0;
}

// distances_to_score (hbond_distance_policy.py): the base-base distance(s)
// whose severity should be taken (max downstream). WC cells return every
// physical bond; other cells return only the single strongest bond within
// [kMinPhysicalDistance, kMaxStrongDistance], or empty when none qualifies.
std::vector<double> distances_to_score(const std::vector<double>& distances,
                                       const std::string& bp, const std::string& lw) {
    std::vector<double> phys;
    for (const double d : distances)
        if (d >= kMinPhysicalDistance) phys.push_back(d);
    if (phys.empty()) return {};
    if (is_wc_distance_cell(bp, lw)) return phys;
    std::vector<double> strong;
    for (const double d : phys)
        if (d <= kMaxStrongDistance) strong.push_back(d);
    if (strong.empty()) return {};
    return {*std::min_element(strong.begin(), strong.end())};
}

/// Cerny severity, ramped continuously off the Preferred boundary:
///
///     severity = min(1, |Z'|/5) * clamp((5 - ProSco) / 5, 0, 1)
///
/// The plain Cerny form (0 if ProSco >= 5, else min(1,|Z'|/5)) is a STEP: the
/// two gates are redundant, so by the time ProSco falls below 5 a typical
/// cell's |Z'| is already ~2.4-3.0 and severity jumps 0 -> ~0.5 with nothing
/// between. A pair was therefore either exactly 100 or had already paid 3-5
/// points, leaving 96-100 unreachable (0.05% of canonical cWW pairs).
///
/// The (5 - ProSco)/5 factor is 0 exactly at the boundary and 1 once ProSco
/// reaches 0, so severity eases in from zero and reduces to the plain Cerny
/// value in the far tail -- the region the method actually calibrates.
///
/// None ProSco -> 0 (no penalty when cell unavailable).
/// None Z' with ProSco < 5 -> the ramp alone (was: maximum severity).
/// Mirrors _severity_from_prosco_zprime in scoring/issues.py.
double severity_from_prosco_zprime(std::optional<double> prosco,
                                   std::optional<double> zprime) {
    if (!prosco || *prosco >= kProscoPreferred) return 0.0;
    const double ramp =
        std::min(1.0, std::max(0.0, (kProscoPreferred - *prosco) / kProscoPreferred));
    if (!zprime) return ramp;
    return std::min(1.0, std::abs(*zprime) / kZprimeThreshold) * ramp;
}

/// Severity plus the ProSco/Z' that produced it (kept for per-issue reporting).
struct SevDetail {
    double severity = 0.0;
    std::optional<double> prosco;
    std::optional<double> zprime;
};

/// Černý three-tier label (Černý et al. NAR 2026 gkaf1335):
///   Preferred   = ProSco >= 5
///   Allowed     = ProSco < 5 and |Z'| <= 5
///   Of Concern  = ProSco < 5 and |Z'| > 5  (or Z' unavailable)
/// Returns "" when the parameter is not ProSco/Z'-scored (e.g. hbond count,
/// where prosco is none) — triggered issues of that kind carry an empty tier.
std::string prosco_tier(std::optional<double> prosco, std::optional<double> zprime) {
    if (!prosco) return "";
    if (*prosco >= kProscoPreferred) return "Preferred";
    if (!zprime) return "Of Concern";
    return std::abs(*zprime) <= kZprimeThreshold ? "Allowed" : "Of Concern";
}

// H-bond angle geometry (angle_at_atom / hbond_angle) lives in the shared
// header scoring/hbond_angles.hpp so the CLI per-pair emitter uses it too.

struct HBF {
    double distance;
    std::optional<double> donor_angle, acceptor_angle;
};

}  // namespace

struct Scorer::Impl {
    json prosco;       // prosco_distributions
    json z_tables;     // z_tables (pairs + hbonds, from z_tables.json)
    json categorical;  // categorical_distributions
    json bb_prosco;    // backbone prosco_distributions["backbone"] (per-conformer torsion)
    std::unordered_map<std::string, double> pair_weights;
    std::unique_ptr<RichardsonClassifier> richardson;

    /// ProSco bin lookup with bp/lw -> _ANY/lw -> _ANY/_ANY fallback.
    std::optional<double> prosco_lookup(std::optional<double> value, const std::string& source,
                                        const std::string& param, const std::string& bp,
                                        const std::string& lw) const {
        if (!value || std::isnan(*value)) return std::nullopt;
        auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
        const double v = maybe_flip(*value, source, param, swapped);
        const auto src = prosco.find(source);
        if (src == prosco.end()) return std::nullopt;
        const auto par = src->find(param);
        if (par == src->end()) return std::nullopt;
        for (const std::string& key : {canon_bp + "/" + lw, "_ANY/" + lw, std::string("_ANY/_ANY")}) {
            const auto it = par->find(key);
            if (it == par->end()) continue;
            if (!it->contains("prosco_per_bin") || (*it)["prosco_per_bin"].empty()) continue;
            const double lo = (*it)["support_lo"].get<double>();
            const double hi = (*it)["support_hi"].get<double>();
            const int n_bins = (*it)["n_bins"].get<int>();
            int idx = static_cast<int>((v - lo) / (hi - lo) * n_bins);
            idx = std::max(0, std::min(n_bins - 1, idx));
            return (*it)["prosco_per_bin"][idx].get<double>();
        }
        return std::nullopt;
    }

    /// Z' asymmetric standard score with bp/lw -> _ANY/lw -> _ANY/_ANY fallback.
    /// Z' = (value - M) / sigma, where sigma = sl (lower half) or su (upper half).
    std::optional<double> zprime_lookup(std::optional<double> value, const std::string& source,
                                        const std::string& param, const std::string& bp,
                                        const std::string& lw) const {
        if (!value || std::isnan(*value)) return std::nullopt;
        auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
        const double v = maybe_flip(*value, source, param, swapped);
        const auto src = z_tables.find(source);
        if (src == z_tables.end()) return std::nullopt;
        const auto par = src->find(param);
        if (par == src->end()) return std::nullopt;
        for (const std::string& key : {canon_bp + "/" + lw, "_ANY/" + lw, std::string("_ANY/_ANY")}) {
            const auto it = par->find(key);
            if (it == par->end() || !it->contains("M")) continue;
            const double M  = (*it)["M"].get<double>();
            const double sl = (*it)["sl"].get<double>();
            const double su = (*it)["su"].get<double>();
            if (v >= M) return su > 1e-9 ? std::optional<double>((v - M) / su) : std::nullopt;
            return sl > 1e-9 ? std::optional<double>((v - M) / sl) : std::nullopt;
        }
        return std::nullopt;
    }

    int expected_hbond_count(const std::string& bp, const std::string& lw) const {
        auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
        (void)swapped;
        const auto src = categorical.find("pairs");
        if (src == categorical.end()) return 1;
        const auto par = src->find("num_base_hbonds");
        if (par == src->end()) return 1;
        for (const std::string& key : {canon_bp + "/" + lw, "_ANY/" + lw, std::string("_ANY/_ANY")}) {
            const auto it = par->find(key);
            if (it == par->end() || !it->contains("frequencies")) continue;
            const auto& freqs = (*it)["frequencies"];
            if (freqs.empty()) continue;
            std::string best_k;
            double best_v = -1.0;
            for (auto fit = freqs.begin(); fit != freqs.end(); ++fit)
                if (fit.value().get<double>() > best_v) { best_v = fit.value().get<double>(); best_k = fit.key(); }
            try { return std::stoi(best_k); } catch (...) { return 1; }
        }
        return 1;
    }

    double hbond_count_severity(const std::string& bp, const std::string& lw, int actual) const {
        const int canonical = expected_hbond_count(bp, lw);
        if (canonical <= 0 || actual >= canonical) return 0.0;
        return static_cast<double>(canonical - actual) / canonical;
    }

    // --- backbone recommendation (per-conformer ProSco; mirror backbone_scorer.py) ---

    /// ProSco of a suite torsion against its empirical per-conformer cell, with
    /// circular unwrap onto the branch used at build time (unwrap_center).
    std::optional<double> bb_prosco_lookup(const std::string& angle,
                                           const std::string& conformer,
                                           double value) const {
        if (conformer.empty() || bb_prosco.is_null()) return std::nullopt;
        const auto ait = bb_prosco.find(angle);
        if (ait == bb_prosco.end()) return std::nullopt;
        const auto cit = ait->find(conformer);
        if (cit == ait->end()) return std::nullopt;
        const json& cell = *cit;
        if (!cell.contains("prosco_per_bin") || cell["prosco_per_bin"].empty())
            return std::nullopt;
        const double uc = cell["unwrap_center"].get<double>();
        const double vu = uc + signed_angle_gap(value, uc);
        const double lo = cell["support_lo"].get<double>();
        const double hi = cell["support_hi"].get<double>();
        const int n = cell["n_bins"].get<int>();
        int idx = static_cast<int>((vu - lo) / (hi - lo) * n);
        idx = std::max(0, std::min(n - 1, idx));
        return cell["prosco_per_bin"][idx].get<double>();
    }

    /// Per-residue backbone recommendation. deviations = only torsions that fire
    /// (ProSco < 5), ranked most-anomalous first. tier fixable / review / flag_only.
    std::optional<BackboneRecommendation> build_recommendation(
            const std::array<double, 7>& suite, const SuiteResult& sr) const {
        std::string tier, target;
        if (sr.is_outlier) {
            tier = "flag_only";
            target = richardson->nearest_conformer(suite);
        } else if (sr.suiteness < kLowSuitenessMax) {
            tier = "fixable";
            target = sr.conformer;
        } else {
            return std::nullopt;
        }
        BackboneRecommendation rec;
        rec.target_conformer = target;
        const std::array<double, 7>* center =
            target.empty() ? nullptr : richardson->center_for(target);
        if (center) {
            const auto& widths = richardson->normal_widths();
            for (int k = 0; k < 7; ++k) {
                const double val = suite[k];
                const double gap = signed_angle_gap(val, (*center)[k]);
                const std::string angle = kSuiteAngleNames[k];
                const std::optional<double> ps = bb_prosco_lookup(angle, target, val);
                const bool fires = ps ? (*ps < kProscoPreferred)
                                      : (std::abs(gap) > kFallbackWidthFactor * widths[k]);
                if (fires) {
                    AngleDeviation d;
                    d.angle = angle; d.value = val; d.target = (*center)[k];
                    d.gap = gap; d.prosco = ps;
                    rec.deviations.push_back(std::move(d));
                }
            }
            std::sort(rec.deviations.begin(), rec.deviations.end(),
                      [](const AngleDeviation& a, const AngleDeviation& b) {
                          const double ka = a.prosco ? *a.prosco : 1000.0 - std::abs(a.gap);
                          const double kb = b.prosco ? *b.prosco : 1000.0 - std::abs(b.gap);
                          return ka < kb;
                      });
        }
        if (tier == "fixable" && rec.deviations.empty()) tier = "review";
        rec.tier = tier;
        return rec;
    }
};

Scorer::Scorer(const std::filesystem::path& distributions_path,
               const std::filesystem::path& prosco_path,
               const std::filesystem::path& weights_path,
               const std::filesystem::path& richardson_suites_path)
    : impl_(std::make_unique<Impl>()) {
    impl_->richardson = std::make_unique<RichardsonClassifier>(richardson_suites_path);
    {
        std::ifstream in(distributions_path);
        json d = json::parse(in);
        if (d.contains("categorical_distributions")) impl_->categorical = d["categorical_distributions"];
    }
    {
        std::ifstream in(prosco_path);
        json d = json::parse(in);
        if (d.contains("prosco_distributions")) impl_->prosco = d["prosco_distributions"];
    }
    // Auto-discover z_tables.json alongside prosco_distributions.json (mirrors Python).
    {
        const auto z_path = prosco_path.parent_path() / "z_tables.json";
        if (std::filesystem::exists(z_path)) {
            std::ifstream in(z_path);
            impl_->z_tables = json::parse(in);
        }
    }
    // Auto-discover backbone_prosco_distributions.json (per-conformer torsion
    // ProSco used for backbone recommendations; mirrors Python scorer).
    {
        const auto bb_path = prosco_path.parent_path() / "backbone_prosco_distributions.json";
        if (std::filesystem::exists(bb_path)) {
            std::ifstream in(bb_path);
            json d = json::parse(in);
            if (d.contains("prosco_distributions") &&
                d["prosco_distributions"].contains("backbone"))
                impl_->bb_prosco = d["prosco_distributions"]["backbone"];
        }
    }
    {
        std::ifstream in(weights_path);
        json d = json::parse(in);
        const json& raw = d["penalty_weights"];
        // Weights already sum to 100 from the Python generator, but renormalise
        // defensively in case the file is from a different generation.
        double total = 0.0;
        for (const char* k : kIssueKeys) total += raw.value(k, 0.0);
        for (const char* k : kIssueKeys)
            impl_->pair_weights[k] = total <= 0 ? 100.0 / kIssueKeys.size()
                                                : raw.value(k, 0.0) * 100.0 / total;
    }
}

Scorer::~Scorer() = default;
Scorer::Scorer(Scorer&&) noexcept = default;

namespace {
/// Residue name out of a "<chain>-<name>-<seq>" res_id.
std::string res_name_of(const std::string& res_id) {
    const auto a = res_id.find('-');
    if (a == std::string::npos) return {};
    const auto b = res_id.find('-', a + 1);
    return res_id.substr(a + 1, (b == std::string::npos ? res_id.size() : b) - a - 1);
}
}  // namespace

bool Scorer::is_scorable(const classification::ScoredCandidate& p,
                         const core::BaseTyping& typing) {
    if (p.pair_category != "base-base") return false;
    if (p.is_ambiguous) return false;
    if (p.lw_class.empty() || p.lw_class.find('|') != std::string::npos) return false;
    // DNA is not scored: the ProSco/Z' cells, the canonical H-bond counts and
    // the penalty weights are all calibrated on RNA. A pair is excluded when
    // EITHER side is DNA, so RNA-DNA hybrid pairs are excluded too. Matches
    // the backbone side, which already skips DNA (see no_rna_suite).
    // Excluded here rather than at candidate generation so the pair is still
    // found, classified and emitted with a null score — same convention as an
    // ambiguous pair — and counted in skipped_pairs.
    if (typing.is_dna(res_name_of(p.res_id1))) return false;
    if (typing.is_dna(res_name_of(p.res_id2))) return false;
    return true;
}

PairScore Scorer::score_pair(const classification::ScoredCandidate& pair, const core::Residue& res1,
                             const core::Residue& res2, hbond::HBondFinder& finder,
                             const hbond::HBondChemistry& chem, const core::BaseTyping& typing) const {
    const PairParameters params = compute_pair_parameters(*pair.frame1, *pair.frame2);
    const std::string bp = typing.normalize(pair.res_name1) + "-" + typing.normalize(pair.res_name2);
    const std::string& lw = pair.lw_class;

    // H-bonds (base-base) with angles.
    std::vector<HBF> hbonds;
    int num_base = 0;
    for (const auto& hb : finder.find_between(res1, res2)) {
        if (!kBaseAtoms.count(hb.donor_atom) || !kBaseAtoms.count(hb.acceptor_atom)) continue;
        const core::Residue* dres = (hb.donor_res_id == pair.res_id1) ? &res1 : &res2;
        const core::Residue* ares = (hb.donor_res_id == pair.res_id1) ? &res2 : &res1;
        const core::Atom* da = dres->get_atom(hb.donor_atom);
        const core::Atom* aa = ares->get_atom(hb.acceptor_atom);
        HBF f{hb.distance, std::nullopt, std::nullopt};
        if (da && aa) {
            f.donor_angle    = hbond_angle(*dres, hb.donor_atom,    da->coords, aa->coords, chem);
            f.acceptor_angle = hbond_angle(*ares, hb.acceptor_atom, aa->coords, da->coords, chem);
        }
        hbonds.push_back(f);
        ++num_base;
    }

    // Helper: Cerny severity for a continuous parameter, plus the ProSco/Z' that
    // produced it (kept so each triggered issue can report its tier).
    auto cerny_sev = [&](std::optional<double> val, const std::string& source,
                         const char* param) -> SevDetail {
        const auto ps = impl_->prosco_lookup(val, source, param, bp, lw);
        const auto zp = impl_->zprime_lookup(val, source, param, bp, lw);
        return {severity_from_prosco_zprime(ps, zp), ps, zp};
    };

    std::unordered_map<std::string, SevDetail> issues;

    // ---- 6 pair geometry parameters, each individually ----
    const std::array<std::pair<const char*, double>, 6> pair_params = {{
        {"shear",     params.shear},
        {"stretch",   params.stretch},
        {"stagger",   params.stagger},
        {"buckle",    params.buckle},
        {"propeller", params.propeller},
        {"opening",   params.opening},
    }};
    for (const auto& [param, val] : pair_params) {
        const auto d = cerny_sev(val, "pairs", param);
        if (d.severity > 0.0) issues[param] = d;
    }

    // ---- H-bond distance: WC cells judge all base-base bonds; other cells
    //      judge only the strong (primary) bond (hbond_distance_policy.py).
    //      Keep the ProSco/Z' of the worst (max-severity) scored bond. ----
    std::vector<double> bond_distances;
    bond_distances.reserve(hbonds.size());
    for (const auto& hb : hbonds) bond_distances.push_back(hb.distance);
    SevDetail dist;
    for (const double d : distances_to_score(bond_distances, bp, lw)) {
        const auto sd = cerny_sev(d, "hbonds", "distance");
        if (sd.severity > dist.severity) dist = sd;
    }
    if (dist.severity > 0.0) issues["distance"] = dist;

    // ---- H-bond angles: max over donor + acceptor angles across all bonds;
    //      keep the ProSco/Z' of the worst (argmax) angle. ----
    SevDetail ang;
    for (const auto& hb : hbonds) {
        const auto sd_d = cerny_sev(hb.donor_angle,    "hbonds", "donor_angle");
        if (sd_d.severity > ang.severity) ang = sd_d;
        const auto sd_a = cerny_sev(hb.acceptor_angle, "hbonds", "acceptor_angle");
        if (sd_a.severity > ang.severity) ang = sd_a;
    }
    if (ang.severity > 0.0) issues["hbond_angles"] = ang;

    // ---- H-bond count: graded (canonical - actual) / canonical (no ProSco/Z') ----
    {
        const double s = impl_->hbond_count_severity(bp, lw, num_base);
        if (s > 0.0) issues["incorrect_hbond_count"] = SevDetail{s, std::nullopt, std::nullopt};
    }

    // ---- Steric clash: graded on heavy-atom overlap (no ProSco/Z') ----
    // Only a contact past kClashThreshold is a clash. Scoring the UNGATED
    // overlap charges score for ordinary van der Waals proximity and dilutes
    // the parameter, whose 2.7x good/poor enrichment is measured at
    // overlap > kClashThreshold.
    //
    // worst_overlap() is the same worst SCORED contact pair_clash() gates on,
    // so applying the threshold here is equivalent and computes the contact
    // once. The raw value is kept on the PairScore as the feature behind the
    // severity. Mirrors features.py (pair_clash) + issues.py.
    const std::optional<double> clash_overlap = worst_overlap(res1, res2);
    if (clash_overlap && *clash_overlap > kClashThreshold) {
        const double s = clash_severity(*clash_overlap);
        if (s > 0.0) issues["clash"] = SevDetail{s, std::nullopt, std::nullopt};
    }

    // Build weighted penalty list in kIssueKeys order, then sort by descending penalty.
    std::vector<IssuePenalty> weighted;
    double penalty = 0.0;
    for (const char* k : kIssueKeys) {
        const auto it = issues.find(k);
        if (it == issues.end()) continue;
        const SevDetail& d = it->second;
        // A weights file generated before `clash` existed simply has no entry:
        // treat it as weight 0 (detected and reported, but costs no score)
        // rather than throwing.
        const auto wit = impl_->pair_weights.find(k);
        const double w = (wit == impl_->pair_weights.end() ? 0.0 : wit->second) * d.severity;
        IssuePenalty ip;
        ip.issue = k;
        ip.weight = w;
        ip.severity = d.severity;
        ip.prosco = d.prosco;
        ip.zprime = d.zprime;
        ip.tier = prosco_tier(d.prosco, d.zprime);
        weighted.push_back(std::move(ip));
        penalty += w;
    }
    std::stable_sort(weighted.begin(), weighted.end(),
                     [](const IssuePenalty& a, const IssuePenalty& b) {
                         return a.weight > b.weight;
                     });

    PairScore out;
    out.res_id1 = pair.res_id1;
    out.res_id2 = pair.res_id2;
    out.bp_type = bp;
    out.lw_class = lw;
    out.score   = std::max(0.0, std::min(100.0, 100.0 - penalty));
    out.penalty = penalty;
    out.issues  = std::move(weighted);
    out.clash_overlap = clash_overlap;
    return out;
}

namespace {

// Longest O3'(i-1) - P(i) accepted as the phosphodiester bond joining two
// residues into one suite. Real links are extremely tight: over a 300-structure
// sample the median is 1.60 A and the 99th percentile 1.62 A, with only 0.06%
// above 2.0 A. Those stragglers are not stretched bonds, they are residues that
// are not joined at all.
constexpr double kMaxPhosphodiester = 2.0;

// Map each res_id to its backbone predecessor's res_id: the residue
// immediately preceding it within the same chain, ordered by (seq, icode),
// when their integer sequence numbers differ by at most one (equal when only
// the insertion code advances, e.g. 20 -> 20A -> 21) AND the two are actually
// bonded. Numbering alone is not enough: 1HMH labels separate strands with the
// same numbers plus an 'L' suffix, so 48 of its 87 suites were being built
// between residues 26 to 40 A apart, which classify as outliers and dragged its
// suiteness from the deposited 0.66 down to 0.283.
std::unordered_map<std::string, std::string> build_predecessor_index(
    const std::unordered_map<std::string, ResidueTorsions>& tors,
    const std::unordered_map<std::string, const core::Residue*>& res_by_id) {
    std::map<std::string, std::vector<std::tuple<long, std::string, std::string>>> by_chain;
    for (const auto& [rid, t] : tors) {
        (void)t;
        const auto p = parse_res_seq(rid);
        if (!p) continue;
        by_chain[p->chain].emplace_back(p->num, p->icode, rid);
    }
    std::unordered_map<std::string, std::string> predecessor;
    for (auto& [chain, residues] : by_chain) {
        (void)chain;
        std::sort(residues.begin(), residues.end());
        for (std::size_t i = 1; i < residues.size(); ++i) {
            const long prev_num = std::get<0>(residues[i - 1]);
            const long cur_num  = std::get<0>(residues[i]);
            if (prev_num != cur_num && prev_num != cur_num - 1) continue;
            const std::string& prev_id = std::get<2>(residues[i - 1]);
            const std::string& cur_id  = std::get<2>(residues[i]);
            const auto pit = res_by_id.find(prev_id);
            const auto cit = res_by_id.find(cur_id);
            if (pit == res_by_id.end() || cit == res_by_id.end()) continue;
            const auto* o3 = pit->second->get_atom("O3'");
            const auto* ph = cit->second->get_atom("P");
            if (!o3 || !ph) continue;
            if ((o3->coords - ph->coords).norm() > kMaxPhosphodiester) continue;
            predecessor[cur_id] = prev_id;
        }
    }
    return predecessor;
}

double mod360(double v) {
    double a = std::fmod(v, 360.0);
    if (a < 0) a += 360.0;
    return a;
}

// build_suite_for: [delta_prev, eps_prev, zeta_prev, alpha, beta, gamma, delta].
std::optional<std::array<double, 7>> build_suite(
    const std::string& res_id, const std::unordered_map<std::string, ResidueTorsions>& tors,
    const std::unordered_map<std::string, std::string>& pred_index) {
    const auto pit = pred_index.find(res_id);
    if (pit == pred_index.end()) return std::nullopt;
    const auto cit  = tors.find(res_id);
    const auto prit = tors.find(pit->second);
    if (cit == tors.end() || prit == tors.end()) return std::nullopt;
    const ResidueTorsions& cur  = cit->second;
    const ResidueTorsions& pred = prit->second;
    const std::array<std::optional<double>, 7> n = {pred.delta, pred.epsilon, pred.zeta,
                                                     cur.alpha, cur.beta, cur.gamma, cur.delta};
    std::array<double, 7> suite{};
    for (int k = 0; k < 7; ++k) {
        if (!n[k]) return std::nullopt;
        suite[k] = mod360(*n[k]);
    }
    return suite;
}

double round2(double v) { return std::round(v * 100.0) / 100.0; }

// DNA residues have no RNA backbone suite — exclude them so they don't deflate
// residues_score for DNA/RNA complexes (matches the wwPDB rnasuiteness RNA-only
// convention and the Python backbone scorer). res_id is chain-resname-num
// (e.g. N-DA-102).
//
// The test is data-driven rather than name-based: modified DNA carries no ``D``
// prefix (2'-fluoro AF2/GF2/CFZ/UFT, BRU, ...), and matching only D<base> let
// those through to be scored as RNA. 3P4A is built entirely from 2'-fluoro
// residues, so it scored 0.811 while the wwPDB — correctly, for an RNA-only
// metric — reports 0.
// A residue carries an RNA backbone suite only if it is a polymer nucleotide.
// DNA has no RNA suite, and a free nucleotide or cofactor (ATP, GTP, ... typed
// LIGAND in the ligand DB) is not a polymer residue at all: it may still pair,
// which is why it is parsed, but it has no place in an RNA backbone score.
bool no_rna_suite(const std::string& res_id, const core::BaseTyping& typing) {
    const auto a = res_id.find('-');
    if (a == std::string::npos) return false;
    const auto b = res_id.find('-', a + 1);
    const std::string name =
        res_id.substr(a + 1, (b == std::string::npos ? res_id.size() : b) - a - 1);
    return typing.is_dna(name) || typing.is_ligand(name);
}

}  // namespace

StructureScore Scorer::score_structure(
    const std::vector<classification::ScoredCandidate>& selected, const core::Structure& structure,
    hbond::HBondFinder& finder, const hbond::HBondChemistry& chem,
    const core::BaseTyping& typing, bool collect_backbone) const {
    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    const auto tors       = compute_all_torsions(structure);
    const auto pred_index = build_predecessor_index(tors, res_by_id);
    // 3'-successor (invert predecessor) for the Pperp phosphate, and helix
    // membership over the selected pairs — both drive the report-only sugar
    // pucker flags (do not affect the score).
    std::unordered_map<std::string, std::string> succ_index;
    succ_index.reserve(pred_index.size());
    for (const auto& [res, pred] : pred_index) succ_index[pred] = res;
    const auto chains = pairfinder::algorithms::RNAChains::from_structure(structure);
    const auto in_helix = pairfinder::algorithms::selection::helix_membership(selected, chains);

    StructureScore ss;
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const auto& sc = selected[i];
        if (!is_scorable(sc, typing)) { ++ss.skipped_pairs; continue; }
        const auto r1 = res_by_id.find(sc.res_id1);
        const auto r2 = res_by_id.find(sc.res_id2);
        if (r1 == res_by_id.end() || r2 == res_by_id.end()) continue;
        PairScore ps = score_pair(sc, *r1->second, *r2->second, finder, chem, typing);
        // "Unusual sugar pucker": C2'-endo (by δ) in a helical cWW G-C/A-U pair.
        if (in_helix[i] && ps.lw_class == "cWW" && kHelicalWcBp.count(ps.bp_type)) {
            for (const core::Residue* res : {r1->second, r2->second}) {
                if (no_rna_suite(res->res_id(), typing)) continue;  // DNA is natively C2'-endo
                const auto tit = tors.find(res->res_id());
                if (tit == tors.end() || !tit->second.delta) continue;
                const double delta = *tit->second.delta;
                if (impl_->richardson->pucker(delta) != 2) continue;  // only C2'-endo
                const core::Residue* next_res = nullptr;
                const auto sit = succ_index.find(res->res_id());
                if (sit != succ_index.end()) {
                    const auto nit = res_by_id.find(sit->second);
                    if (nit != res_by_id.end()) next_res = nit->second;
                }
                const auto pp = compute_pperp(res, next_res);
                UnusualPucker up;
                up.res_id = res->res_id();
                up.delta = delta;
                up.pperp = pp;
                up.pperp_confirms = (pp && *pp < kPperpCutoff);
                ps.unusual_pucker.push_back(std::move(up));
            }
        }
        ss.pair_scores.push_back(std::move(ps));
    }
    for (const auto& [rid, t] : tors) {
        const auto suite = build_suite(rid, tors, pred_index);
        if (!suite) continue;
        const SuiteResult sr = impl_->richardson->classify(*suite);
        // Full torsion row for every classifiable residue (incl. DNA, outliers,
        // L-RNA) — captured before the RNA-only / L-RNA scoring filters below so
        // it mirrors extract_backbone_torsions.py (torsions_all.csv).
        if (collect_backbone) {
            BackboneTorsion bt;
            bt.res_id = rid;
            if (const auto rit = res_by_id.find(rid); rit != res_by_id.end())
                bt.base_type = typing.normalize(rit->second->base_type());
            bt.alpha = t.alpha; bt.beta = t.beta; bt.gamma = t.gamma;
            bt.delta = t.delta; bt.epsilon = t.epsilon; bt.zeta = t.zeta;
            bt.conformer = sr.conformer;
            bt.suiteness = sr.suiteness;
            bt.suite_bin = sr.bin;
            bt.is_outlier = sr.is_outlier;
            ss.backbone_torsions.push_back(std::move(bt));
        }
        if (no_rna_suite(rid, typing)) continue;  // RNA-only backbone score (see helper)
        // ...and the suite's delta/epsilon/zeta come from the PREDECESSOR, so a
        // residue whose 5' neighbour is not RNA polymer has no RNA suite either.
        // 3P4C/3P4D alternate RNA with 2'-fluoro residues, so every suite there
        // was scored off a DNA sugar; 6E0O/4HOS/4HOR begin with an initiating
        // ATP/UTP/CTP, which is covalently linked but is a free nucleotide, not
        // a polymer residue. The wwPDB reports ~0 for both cases.
        if (const auto pit = pred_index.find(rid);
            pit != pred_index.end() && no_rna_suite(pit->second, typing))
            continue;
        if (sr.suiteness <= 0.0) {
            // L-RNA enantiomer check: a rejected suite may be mirror-image RNA
            // (e.g. a Spiegelmer aptamer). Negate every torsion (360 - t) and
            // re-classify; if the mirror fits a D-RNA cluster it is L-RNA —
            // exclude it (D-RNA suiteness doesn't apply; the backbone is
            // genuinely high quality). A real D-RNA outlier won't fit even when
            // mirrored, so it still scores 0. Mirrors the Python backbone scorer.
            std::array<double, 7> mirror;
            for (int k = 0; k < 7; ++k)
                mirror[k] = std::fmod(360.0 - (*suite)[k], 360.0);
            if (impl_->richardson->suiteness(mirror) > 0.0) continue;
        }
        ResidueScore rs;
        rs.res_id   = rid;
        rs.suiteness = sr.suiteness;
        rs.score    = std::max(0.0, std::min(100.0, sr.suiteness * 100.0));
        rs.recommendation = impl_->build_recommendation(*suite, sr);
        // Sugar pucker: δ, Pperp, and the report-only OUTLIER (δ/Pperp disagree).
        if (t.delta) {
            const double delta = *t.delta;
            rs.delta = delta;
            const int pd = impl_->richardson->pucker(delta);
            const core::Residue* residue = nullptr;
            if (const auto rit = res_by_id.find(rid); rit != res_by_id.end())
                residue = rit->second;
            const core::Residue* next_res = nullptr;
            if (const auto sit = succ_index.find(rid); sit != succ_index.end())
                if (const auto nit = res_by_id.find(sit->second); nit != res_by_id.end())
                    next_res = nit->second;
            if (const auto pp = compute_pperp(residue, next_res)) {
                rs.pperp = pp;
                if (pd == 2 || pd == 3) {
                    const int pp_pucker = pucker_from_pperp(*pp);
                    if (pp_pucker != pd) {
                        PuckerOutlier po;
                        po.pperp = *pp;
                        po.delta = delta;
                        po.pucker_by_pperp = pucker_label(pp_pucker);
                        po.pucker_by_delta = pucker_label(pd);
                        rs.pucker_outlier = std::move(po);
                    }
                }
            }
        }
        ss.residue_scores.push_back(std::move(rs));
    }

    auto mean = [](const auto& v, auto getter) {
        if (v.empty()) return 0.0;
        double s = 0.0;
        for (const auto& x : v) s += getter(x);
        return s / v.size();
    };
    const double pm = mean(ss.pair_scores,    [](const PairScore& p)    { return p.score; });
    const double rm = mean(ss.residue_scores, [](const ResidueScore& r) { return r.score; });
    ss.pairs_score    = round2(pm);
    ss.residues_score = round2(rm);

    // Whole-structure RNA-RNA clashes (report-only; includes residues that are
    // not a base pair). structure_clashes takes an owning map, so build one
    // from the parsed chains.
    {
        std::unordered_map<std::string, core::Residue> owned;
        owned.reserve(res_by_id.size());
        for (const auto& [rid, res] : res_by_id) owned.emplace(rid, *res);
        ss.clashes = structure_clashes(owned, typing);
    }
    return ss;
}

}  // namespace pairfinder::scoring
