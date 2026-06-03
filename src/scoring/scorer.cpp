/**
 * @file scorer.cpp
 * @brief Empirical pair scorer (port of scorer.py + issues.py + features.py + canonical.py).
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

constexpr double kProscoPreferred = 25.0;
constexpr double kProscoOfConcern = 5.0;
constexpr double kPi = 3.14159265358979323846;

const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
const std::unordered_set<std::string> kSelfReciprocal = {
    "cWW", "tWW", "cHH", "tHH", "cSS", "tSS"};
const std::unordered_set<std::string> kSignFlip = {"shear", "stagger", "buckle"};
const std::array<const char*, 6> kIssueKeys = {
    "misaligned", "non_coplanar", "rotational_distortion",
    "bad_hbond_distance", "bad_hbond_angles", "incorrect_hbond_count"};

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

double severity_from_prosco(std::optional<double> p) {
    if (!p) return 0.0;
    if (*p >= kProscoPreferred) return 0.0;
    if (*p <= kProscoOfConcern) return 1.0;
    return (kProscoPreferred - *p) / (kProscoPreferred - kProscoOfConcern);
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
    json prosco;        // prosco_distributions
    json categorical;   // categorical_distributions
    std::unordered_map<std::string, double> pair_weights;
    double w_pairs = 0.5, w_residues = 0.5;
    std::unique_ptr<RichardsonClassifier> richardson;

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
    {
        std::ifstream in(weights_path);
        json d = json::parse(in);
        const json& raw = d["penalty_weights"];
        double total = 0.0;
        for (const char* k : kIssueKeys) total += raw.value(k, 0.0);
        for (const char* k : kIssueKeys)
            impl_->pair_weights[k] = total <= 0 ? 100.0 / kIssueKeys.size()
                                                : raw.value(k, 0.0) * 100.0 / total;
        if (d.contains("bucket_weights")) {
            impl_->w_pairs = d["bucket_weights"].value("W_PAIRS", 0.5);
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
            f.donor_angle = hbond_angle(*dres, hb.donor_atom, da->coords, aa->coords, chem);
            f.acceptor_angle = hbond_angle(*ares, hb.acceptor_atom, aa->coords, da->coords, chem);
        }
        hbonds.push_back(f);
        ++num_base;
    }

    auto sev_pair = [&](const char* param, double value) {
        return severity_from_prosco(impl_->prosco_lookup(value, "pairs", param, bp, lw));
    };
    auto sev_hb = [&](const char* param, std::optional<double> value) {
        return severity_from_prosco(impl_->prosco_lookup(value, "hbonds", param, bp, lw));
    };

    std::unordered_map<std::string, double> issues;
    double s;
    s = std::max(sev_pair("shear", params.shear), sev_pair("stretch", params.stretch));
    if (s > 0) issues["misaligned"] = s;
    s = std::max(sev_pair("stagger", params.stagger), sev_pair("buckle", params.buckle));
    if (s > 0) issues["non_coplanar"] = s;
    s = std::max(sev_pair("propeller", params.propeller), sev_pair("opening", params.opening));
    if (s > 0) issues["rotational_distortion"] = s;

    double bad_dist = 0.0, bad_ang = 0.0;
    for (const auto& hb : hbonds) {
        bad_dist = std::max(bad_dist, sev_hb("distance", hb.distance));
        bad_ang = std::max({bad_ang, sev_hb("donor_angle", hb.donor_angle),
                            sev_hb("acceptor_angle", hb.acceptor_angle)});
    }
    if (bad_dist > 0) issues["bad_hbond_distance"] = bad_dist;
    if (bad_ang > 0) issues["bad_hbond_angles"] = bad_ang;

    s = impl_->hbond_count_severity(bp, lw, num_base);
    if (s > 0) issues["incorrect_hbond_count"] = s;

    // Weighted penalties; issue list ordered by descending weight (scorer.py).
    // Iterate in canonical kIssueKeys order (= Python dict insertion order) then
    // stable-sort by descending weight so ties break exactly like Python.
    std::vector<std::pair<std::string, double>> weighted;
    double penalty = 0.0;
    for (const char* k : kIssueKeys) {
        const auto it = issues.find(k);
        if (it == issues.end()) continue;
        const double w = impl_->pair_weights.at(k) * it->second;
        weighted.emplace_back(k, w);
        penalty += w;
    }
    std::stable_sort(weighted.begin(), weighted.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
    const double score = std::max(0.0, std::min(100.0, 100.0 - penalty));

    PairScore out;
    out.res_id1 = pair.res_id1;
    out.res_id2 = pair.res_id2;
    out.bp_type = bp;
    out.lw_class = lw;
    out.score = score;
    out.penalty = penalty;
    for (const auto& [k, w] : weighted) {
        (void)w;
        out.issues.push_back(k);
    }
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
            const long cur_num = std::get<0>(residues[i]);
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
    const auto cit = tors.find(res_id);
    const auto prit = tors.find(pit->second);
    if (cit == tors.end() || prit == tors.end()) return std::nullopt;
    const ResidueTorsions& cur = cit->second;
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

    const auto tors = compute_all_torsions(structure);
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
        rs.res_id = rid;
        rs.suiteness = suiteness;
        rs.score = std::max(0.0, std::min(100.0, suiteness * 100.0));
        ss.residue_scores.push_back(rs);
    }

    auto mean = [](const auto& v, auto getter) {
        if (v.empty()) return 0.0;
        double s = 0.0;
        for (const auto& x : v) s += getter(x);
        return s / v.size();
    };
    const double pm = mean(ss.pair_scores, [](const PairScore& p) { return p.score; });
    const double rm = mean(ss.residue_scores, [](const ResidueScore& r) { return r.score; });
    ss.pairs_score = round2(pm);
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
