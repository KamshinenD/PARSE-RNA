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

#include <pairfinder/scoring/pair_parameters.hpp>
#include <pairfinder/scoring/richardson.hpp>
#include <pairfinder/scoring/torsions.hpp>

namespace pairfinder::scoring {

namespace {

using json = nlohmann::json;
using geometry::Vector3d;

// Cerny severity thresholds.
constexpr double kProscoPreferred  = 5.0;   // ProSco >= 5 -> Preferred, severity = 0
constexpr double kZprimeThreshold  = 5.0;   // |Z'| = 5 -> severity = 1
constexpr double kPi = 3.14159265358979323846;

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

// 9 individual parameters — mirrors Python ISSUE_KEYS.
const std::array<const char*, 9> kIssueKeys = {
    "shear", "stretch", "stagger", "buckle", "propeller", "opening",
    "distance", "hbond_angles", "incorrect_hbond_count"};

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

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

/// Cerny severity: 0 if ProSco >= 5 (Preferred); else min(1, |Z'|/5).
/// None ProSco -> 0 (no penalty when cell unavailable).
/// None Z' with ProSco < 5 -> 1 (maximum severity: outside distribution, no Z' data).
double severity_from_prosco_zprime(std::optional<double> prosco,
                                   std::optional<double> zprime) {
    if (!prosco) return 0.0;
    if (*prosco >= kProscoPreferred) return 0.0;
    if (!zprime) return 1.0;
    return std::min(1.0, std::abs(*zprime) / kZprimeThreshold);
}

// ---- H-bond angle geometry (hbond_angles.py) ----
double angle_at_atom(const Vector3d& neighbor, const Vector3d& atom, const Vector3d& partner) {
    const Vector3d v1 = neighbor - atom, v2 = partner - atom;
    const double n1 = v1.norm(), n2 = v2.norm();
    if (n1 < 1e-10 || n2 < 1e-10) return 0.0;
    double c = v1.dot(v2) / (n1 * n2);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * 180.0 / kPi;
}

std::optional<double> hbond_angle(const core::Residue& res, const std::string& atom_name,
                                  const Vector3d& atom_pos, const Vector3d& partner_pos,
                                  const hbond::HBondChemistry& chem) {
    const std::vector<std::string>* nbrs = &chem.connectivity(res.base_type(), atom_name);
    if (nbrs->empty() && !res.base_type().empty())
        nbrs = &chem.connectivity(res.base_type().substr(0, 1), atom_name);
    if (nbrs->empty()) return std::nullopt;
    std::optional<double> best;
    for (const auto& nn : *nbrs) {
        const core::Atom* na = res.get_atom(nn);
        if (!na) continue;
        const double a = angle_at_atom(na->coords, atom_pos, partner_pos);
        if (!best || a > *best) best = a;
    }
    return best;
}

struct HBF {
    double distance;
    std::optional<double> donor_angle, acceptor_angle;
};

}  // namespace

struct Scorer::Impl {
    json prosco;       // prosco_distributions
    json z_tables;     // z_tables (pairs + hbonds, from z_tables.json)
    json categorical;  // categorical_distributions
    std::unordered_map<std::string, double> pair_weights;
    double w_pairs = 0.5, w_residues = 0.5;
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
        if (d.contains("bucket_weights")) {
            impl_->w_pairs    = d["bucket_weights"].value("W_PAIRS",    0.5);
            impl_->w_residues = d["bucket_weights"].value("W_RESIDUES", 0.5);
        }
    }
}

Scorer::~Scorer() = default;
Scorer::Scorer(Scorer&&) noexcept = default;
double Scorer::w_pairs() const { return impl_->w_pairs; }
double Scorer::w_residues() const { return impl_->w_residues; }

bool Scorer::is_scorable(const classification::ScoredCandidate& p) {
    if (p.pair_category != "base-base") return false;
    if (p.is_ambiguous) return false;
    if (p.lw_class.empty() || p.lw_class.find('|') != std::string::npos) return false;
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

    // Helper: Cerny severity for a continuous parameter.
    auto cerny_sev = [&](std::optional<double> val, const std::string& source,
                         const char* param) -> double {
        const auto ps = impl_->prosco_lookup(val, source, param, bp, lw);
        const auto zp = impl_->zprime_lookup(val, source, param, bp, lw);
        return severity_from_prosco_zprime(ps, zp);
    };

    std::unordered_map<std::string, double> issues;

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
        const double s = cerny_sev(val, "pairs", param);
        if (s > 0.0) issues[param] = s;
    }

    // ---- H-bond distance: WC cells judge all base-base bonds; other cells
    //      judge only the strong (primary) bond (hbond_distance_policy.py) ----
    std::vector<double> bond_distances;
    bond_distances.reserve(hbonds.size());
    for (const auto& hb : hbonds) bond_distances.push_back(hb.distance);
    double dist_sev = 0.0;
    for (const double d : distances_to_score(bond_distances, bp, lw))
        dist_sev = std::max(dist_sev, cerny_sev(d, "hbonds", "distance"));
    if (dist_sev > 0.0) issues["distance"] = dist_sev;

    // ---- H-bond angles: max over donor + acceptor angles across all bonds ----
    double angle_sev = 0.0;
    for (const auto& hb : hbonds) {
        angle_sev = std::max(angle_sev, cerny_sev(hb.donor_angle,    "hbonds", "donor_angle"));
        angle_sev = std::max(angle_sev, cerny_sev(hb.acceptor_angle, "hbonds", "acceptor_angle"));
    }
    if (angle_sev > 0.0) issues["hbond_angles"] = angle_sev;

    // ---- H-bond count: graded (canonical - actual) / canonical ----
    {
        const double s = impl_->hbond_count_severity(bp, lw, num_base);
        if (s > 0.0) issues["incorrect_hbond_count"] = s;
    }

    // Build weighted penalty list in kIssueKeys order, then sort by descending penalty.
    std::vector<IssuePenalty> weighted;
    double penalty = 0.0;
    for (const char* k : kIssueKeys) {
        const auto it = issues.find(k);
        if (it == issues.end()) continue;
        const double w = impl_->pair_weights.at(k) * it->second;
        weighted.push_back({k, w});
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
    return out;
}

namespace {

// Map each res_id to its backbone predecessor's res_id: the residue
// immediately preceding it within the same chain, ordered by (seq, icode),
// when their integer sequence numbers differ by at most one (equal when only
// the insertion code advances, e.g. 20 -> 20A -> 21). Mirrors the Python
// features.build_predecessor_index. Negative seqs / insertion codes included.
std::unordered_map<std::string, std::string> build_predecessor_index(
    const std::unordered_map<std::string, ResidueTorsions>& tors) {
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
            if (prev_num == cur_num || prev_num == cur_num - 1)
                predecessor[std::get<2>(residues[i])] = std::get<2>(residues[i - 1]);
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

}  // namespace

StructureScore Scorer::score_structure(
    const std::vector<classification::ScoredCandidate>& selected, const core::Structure& structure,
    hbond::HBondFinder& finder, const hbond::HBondChemistry& chem,
    const core::BaseTyping& typing) const {
    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    const auto tors       = compute_all_torsions(structure);
    const auto pred_index = build_predecessor_index(tors);

    StructureScore ss;
    for (const auto& sc : selected) {
        if (!is_scorable(sc)) { ++ss.skipped_pairs; continue; }
        const auto r1 = res_by_id.find(sc.res_id1);
        const auto r2 = res_by_id.find(sc.res_id2);
        if (r1 == res_by_id.end() || r2 == res_by_id.end()) continue;
        ss.pair_scores.push_back(score_pair(sc, *r1->second, *r2->second, finder, chem, typing));
    }
    for (const auto& [rid, t] : tors) {
        (void)t;
        const auto suite = build_suite(rid, tors, pred_index);
        if (!suite) continue;
        const double suiteness = impl_->richardson->suiteness(*suite);
        ResidueScore rs;
        rs.res_id   = rid;
        rs.suiteness = suiteness;
        rs.score    = std::max(0.0, std::min(100.0, suiteness * 100.0));
        ss.residue_scores.push_back(rs);
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
    double overall;
    if (!ss.pair_scores.empty() && !ss.residue_scores.empty())
        overall = impl_->w_pairs * pm + impl_->w_residues * rm;
    else if (!ss.pair_scores.empty())
        overall = pm;
    else if (!ss.residue_scores.empty())
        overall = rm;
    else
        overall = 0.0;
    ss.overall = round2(overall);
    return ss;
}

}  // namespace pairfinder::scoring
