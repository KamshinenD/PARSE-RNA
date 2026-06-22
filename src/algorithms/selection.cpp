/**
 * @file selection.cpp
 * @brief Two-phase pair selection (port of finder/selection.py + helix.py).
 */
#include <pairfinder/algorithms/selection.hpp>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace pairfinder::algorithms::selection {

namespace {

using classification::ScoredCandidate;
using Slot = std::pair<std::string, char>;
using SlotSet = std::set<Slot>;
using Cand = ScoredCandidate;

const std::unordered_set<std::string> kCanonicalSeq = {"AU", "UA", "GC", "CG", "GU", "UG"};

std::string sequence(const Cand& c) { return c.res_name1 + c.res_name2; }

int face_priority(char f) { return f == 'W' ? 0 : (f == 'H' ? 1 : 2); }

std::vector<std::string> split_pipe(const std::string& lw) {
    std::vector<std::string> out;
    std::size_t start = 0, p;
    while ((p = lw.find('|', start)) != std::string::npos) {
        out.push_back(lw.substr(start, p - start));
        start = p + 1;
    }
    out.push_back(lw.substr(start));
    return out;
}

// ---- stacking-signature filter (validation/lw_class_filters.py) ----
bool is_stacking_signature(const Cand& c) {
    static const std::unordered_map<std::string, std::pair<double, double>> thr = {
        {"cWW", {7.5, 1.5}}, {"cWH", {7.5, 1.5}}, {"cHH", {7.5, 1.5}}, {"tWW", {7.5, 1.5}},
        {"tHH", {7.5, 1.5}}, {"cHS", {7.0, 1.5}}, {"tHS", {6.5, 1.5}}, {"tSS", {6.5, 1.5}}};
    const auto it = thr.find(c.lw_class);
    if (it == thr.end()) return false;
    return c.validation.dNN < it->second.first && c.validation.d_v > it->second.second;
}

// ---- slot helpers (selection.py) ----
bool same_chain_adjacent(const std::string& a, const std::string& b) {
    auto split = [](const std::string& r, std::string& chain, long& num) {
        const auto p2 = r.rfind('-');
        if (p2 == std::string::npos || p2 == 0) return false;
        const auto p1 = r.rfind('-', p2 - 1);
        if (p1 == std::string::npos) return false;
        chain = r.substr(0, p1);
        try { num = std::stol(r.substr(p2 + 1)); } catch (...) { return false; }
        return true;
    };
    std::string c1, c2;
    long n1 = 0, n2 = 0;
    if (!split(a, c1, n1) || !split(b, c2, n2)) return false;
    return c1 == c2 && std::abs(n1 - n2) == 1;
}

std::optional<SlotSet> candidate_slots(const Cand& c) {
    if (c.lw_class.empty()) return std::nullopt;
    SlotSet slots;
    for (const auto& part : split_pipe(c.lw_class)) {
        if (part.size() < 3) continue;
        slots.insert({c.res_id1, part[1]});
        slots.insert({c.res_id2, part[2]});
    }
    if (slots.empty()) return std::nullopt;
    return slots;
}

char face_from_lw(const std::string& lw, int which) {
    if (lw.size() >= 3) return which == 1 ? lw[1] : lw[2];
    return 'W';
}

// resolve_for_selection (MUTATES c.lw_class). Returns slots to claim or nullopt.
std::optional<SlotSet> resolve_for_selection(Cand& c, const SlotSet& used) {
    const std::string lw = c.lw_class;
    const std::string& r1 = c.res_id1;
    const std::string& r2 = c.res_id2;
    if (lw.find('|') == std::string::npos) {
        const char f1 = face_from_lw(lw, 1), f2 = face_from_lw(lw, 2);
        if (used.count({r1, f1}) || used.count({r2, f2})) return std::nullopt;
        return SlotSet{{r1, f1}, {r2, f2}};
    }
    std::vector<std::tuple<std::string, char, char>> avail;
    for (const auto& part : split_pipe(lw)) {
        const char f1 = part.size() >= 3 ? part[1] : 'W';
        const char f2 = part.size() >= 3 ? part[2] : 'W';
        if (!used.count({r1, f1}) && !used.count({r2, f2})) avail.emplace_back(part, f1, f2);
    }
    if (avail.empty()) return std::nullopt;
    if (avail.size() == 1) {
        c.lw_class = std::get<0>(avail[0]);  // resolve
        return SlotSet{{r1, std::get<1>(avail[0])}, {r2, std::get<2>(avail[0])}};
    }
    SlotSet slots;
    for (const auto& [part, f1, f2] : avail) {
        slots.insert({r1, f1});
        slots.insert({r2, f2});
    }
    return slots;
}

std::vector<std::pair<Slot, Slot>> possible_face_slots(const std::string& lw,
                                                       const std::string& r1,
                                                       const std::string& r2) {
    if (lw.empty()) return {{{r1, 'W'}, {r2, 'W'}}};
    std::vector<std::pair<Slot, Slot>> out;
    for (const auto& part : split_pipe(lw)) {
        const char f1 = part.size() >= 3 ? part[1] : 'W';
        const char f2 = part.size() >= 3 ? part[2] : 'W';
        out.push_back({{r1, f1}, {r2, f2}});
    }
    return out;
}

// ---- HelixDetector (helix.py) ----
bool pairs_consecutive(const Cand& p1, const Cand& p2, const RNAChains& ch) {
    const std::string& a1 = p1.res_id1;
    const std::string& a2 = p1.res_id2;
    const std::string& b1 = p2.res_id1;
    const std::string& b2 = p2.res_id2;
    const std::string& next_1a = ch.next_residue(a1);
    const std::string& prev_2a = ch.prev_residue(a2);
    if (!next_1a.empty() && next_1a == b1 && !prev_2a.empty() && prev_2a == b2) return true;
    if (!next_1a.empty() && next_1a == b2 && !prev_2a.empty() && prev_2a == b1) return true;
    const std::string& prev_1a = ch.prev_residue(a1);
    const std::string& next_2a = ch.next_residue(a2);
    if (!prev_1a.empty() && prev_1a == b1 && !next_2a.empty() && next_2a == b2) return true;
    if (!prev_1a.empty() && prev_1a == b2 && !next_2a.empty() && next_2a == b1) return true;
    if (!next_2a.empty() && next_2a == b1 && !prev_1a.empty() && prev_1a == b2) return true;
    if (!next_2a.empty() && next_2a == b2 && !prev_1a.empty() && prev_1a == b1) return true;
    return false;
}

struct Segment {
    std::vector<int> pair_indices;  // into helix_cands
    double total_score = 0.0;
    bool is_isolated() const { return pair_indices.size() == 1; }
};

std::vector<Segment> detect_segments(const std::vector<Cand*>& cands, const RNAChains& ch) {
    const int n = static_cast<int>(cands.size());
    std::vector<std::vector<int>> adj(n);
    // O(n) edge build: a pair is consecutive only with the candidate whose
    // residues are {next(a1),prev(a2)} or {prev(a1),next(a2)}. Index by residue
    // pair and look those two up, instead of testing all O(n^2) pairs. Identical
    // edge set -> identical components/segments (pairs_consecutive kept as spec).
    auto pkey = [](const std::string& x, const std::string& y) {
        return x < y ? x + '\x01' + y : y + '\x01' + x;
    };
    std::unordered_map<std::string, std::vector<int>> by_pair;
    by_pair.reserve(static_cast<std::size_t>(n) * 2);
    for (int i = 0; i < n; ++i)
        by_pair[pkey(cands[i]->res_id1, cands[i]->res_id2)].push_back(i);
    auto link = [&](int i, const std::string& r1, const std::string& r2) {
        if (r1.empty() || r2.empty()) return;
        const auto it = by_pair.find(pkey(r1, r2));
        if (it == by_pair.end()) return;
        for (int j : it->second)
            if (j != i) { adj[i].push_back(j); adj[j].push_back(i); }
    };
    for (int i = 0; i < n; ++i) {
        const std::string& a1 = cands[i]->res_id1;
        const std::string& a2 = cands[i]->res_id2;
        link(i, ch.next_residue(a1), ch.prev_residue(a2));
        link(i, ch.prev_residue(a1), ch.next_residue(a2));
    }
    std::vector<char> visited(n, 0);
    std::vector<Segment> segs;
    for (int s = 0; s < n; ++s) {
        if (visited[s]) continue;
        std::vector<int> comp;
        std::vector<int> q = {s};
        visited[s] = 1;
        std::size_t head = 0;
        while (head < q.size()) {
            const int node = q[head++];
            comp.push_back(node);
            for (int nb : adj[node])
                if (!visited[nb]) { visited[nb] = 1; q.push_back(nb); }
        }
        auto sort_key = [&](int idx) -> std::pair<int, int> {
            auto p = ch.position(cands[idx]->res_id1);
            if (!p) p = ch.position(cands[idx]->res_id2);
            return p ? *p : std::pair<int, int>{999, 999};
        };
        std::stable_sort(comp.begin(), comp.end(),
                         [&](int a, int b) { return sort_key(a) < sort_key(b); });
        Segment seg;
        seg.pair_indices = comp;
        for (int i : comp) seg.total_score += cands[i]->quality_score;
        segs.push_back(std::move(seg));
    }
    std::stable_sort(segs.begin(), segs.end(), [](const Segment& a, const Segment& b) {
        if (a.pair_indices.size() != b.pair_indices.size())
            return a.pair_indices.size() > b.pair_indices.size();
        return a.total_score > b.total_score;
    });
    return segs;
}

// ---- _CWWHelixPhase ----
constexpr double kMinHelixQuality = 0.45;
constexpr double kMinIsolatedQuality = 0.45;

std::vector<Cand*> cww_helix_phase(const std::vector<Cand*>& cww, const RNAChains& ch) {
    std::vector<Cand*> canonical;
    for (Cand* c : cww)
        if (!same_chain_adjacent(c->res_id1, c->res_id2) && c->has_strong_base_hbond)
            canonical.push_back(c);
    if (canonical.empty()) return {};
    std::vector<Cand*> helix_cands;
    for (Cand* c : canonical)
        if (c->quality_score >= kMinHelixQuality) helix_cands.push_back(c);

    const auto segs = detect_segments(helix_cands, ch);

    std::unordered_map<std::string, double> best_for_res;
    for (Cand* c : helix_cands) {
        for (const std::string& r : {c->res_id1, c->res_id2}) {
            auto it = best_for_res.find(r);
            if (it == best_for_res.end()) best_for_res[r] = c->quality_score;
            else it->second = std::max(it->second, c->quality_score);
        }
    }

    std::vector<Cand*> selected;
    std::unordered_set<std::string> used;
    std::unordered_set<int> helix_idx;
    // Phase 1: helix pairs.
    for (const auto& seg : segs) {
        if (seg.is_isolated()) continue;
        for (int idx : seg.pair_indices) {
            Cand* c = helix_cands[idx];
            if (used.count(c->res_id1) || used.count(c->res_id2)) continue;
            if (c->quality_score < kMinHelixQuality) continue;
            bool dominated = false;
            for (const std::string& r : {c->res_id1, c->res_id2})
                if (best_for_res[r] - c->quality_score > 0.15) { dominated = true; break; }
            if (dominated) continue;
            selected.push_back(c);
            used.insert(c->res_id1);
            used.insert(c->res_id2);
            helix_idx.insert(idx);
        }
    }
    // Phase 2: isolated pairs — extension-aware.
    // potential extensions: helix ends -> candidate extensions.
    std::unordered_set<const Cand*> selected_set(selected.begin(), selected.end());
    std::vector<Cand*> extensions;  // flat list (we only need the union for overlap check)
    // O(n) via residue-pair indexes (same neighbor-pair trick as detect_segments):
    // a pair is consecutive only with candidates whose residues are its two
    // chain-stepped neighbor-pairs, so look those up instead of scanning all.
    auto pkey2 = [](const std::string& x, const std::string& y) {
        return x < y ? x + '\x01' + y : y + '\x01' + x;
    };
    std::unordered_map<std::string, std::vector<Cand*>> sel_by_pair, can_by_pair;
    for (Cand* c : selected) sel_by_pair[pkey2(c->res_id1, c->res_id2)].push_back(c);
    for (Cand* c : canonical) can_by_pair[pkey2(c->res_id1, c->res_id2)].push_back(c);
    for (Cand* pr : selected) {
        const std::string nbr[2][2] = {
            {ch.next_residue(pr->res_id1), ch.prev_residue(pr->res_id2)},
            {ch.prev_residue(pr->res_id1), ch.next_residue(pr->res_id2)}};
        int consec = 0;
        for (const auto& s : nbr) {
            if (s[0].empty() || s[1].empty()) continue;
            const auto it = sel_by_pair.find(pkey2(s[0], s[1]));
            if (it != sel_by_pair.end())
                for (Cand* o : it->second)
                    if (o != pr) ++consec;
        }
        if (consec > 1) continue;
        for (const auto& s : nbr) {
            if (s[0].empty() || s[1].empty()) continue;
            const auto it = can_by_pair.find(pkey2(s[0], s[1]));
            if (it == can_by_pair.end()) continue;
            for (Cand* cand : it->second) {
                if (selected_set.count(cand)) continue;
                if (used.count(cand->res_id1) || used.count(cand->res_id2)) continue;
                if (cand->quality_score < kMinHelixQuality) continue;
                extensions.push_back(cand);
            }
        }
    }
    std::vector<int> isolated;
    for (const auto& seg : segs)
        if (seg.is_isolated())
            for (int idx : seg.pair_indices) isolated.push_back(idx);
    std::stable_sort(isolated.begin(), isolated.end(), [&](int a, int b) {
        return helix_cands[a]->quality_score > helix_cands[b]->quality_score;
    });
    for (int idx : isolated) {
        Cand* c = helix_cands[idx];
        if (used.count(c->res_id1) || used.count(c->res_id2)) continue;
        if (c->quality_score < kMinIsolatedQuality) continue;
        bool blocked = false;
        for (Cand* ext : extensions) {
            const bool overlaps = c->res_id1 == ext->res_id1 || c->res_id1 == ext->res_id2 ||
                                  c->res_id2 == ext->res_id1 || c->res_id2 == ext->res_id2;
            if (overlaps && c->quality_score - ext->quality_score <= 0.15) { blocked = true; break; }
        }
        if (blocked) continue;
        selected.push_back(c);
        used.insert(c->res_id1);
        used.insert(c->res_id2);
    }
    return selected;
}

// ---- shared-residue tie-break: H-bond strength first, then planarity ----
// When several candidates contest one residue, prefer the strongest interaction
// (most base-base H-bonds), breaking ties by coplanarity (plane_angle, then d_v).
// Coplanarity alone over-weights flatness vs bond strength; H-bonds-first wins.
using PlanarityKey = std::tuple<int, double, double>;

PlanarityKey planarity_rank_key(const Cand& c) {
    return {-c.num_hbonds, c.validation.plane_angle, c.validation.d_v};
}

std::optional<std::vector<int>> resolve_shared_residue_conflict(
    const std::vector<Cand*>& candidates, const std::vector<int>& component,
    const std::vector<std::vector<int>>& conflicts) {
    if (component.size() < 2) return std::nullopt;

    std::unordered_set<std::string> shared = {candidates[component[0]]->res_id1,
                                              candidates[component[0]]->res_id2};
    for (std::size_t i = 1; i < component.size(); ++i) {
        std::unordered_set<std::string> pair_res = {candidates[component[i]]->res_id1,
                                                    candidates[component[i]]->res_id2};
        std::unordered_set<std::string> inter;
        for (const auto& r : shared)
            if (pair_res.count(r)) inter.insert(r);
        shared = std::move(inter);
        if (shared.empty()) return std::nullopt;
    }

    for (std::size_t ii = 0; ii < component.size(); ++ii)
        for (std::size_t jj = ii + 1; jj < component.size(); ++jj) {
            const int i = component[ii], j = component[jj];
            bool found = false;
            for (int nb : conflicts[i])
                if (nb == j) {
                    found = true;
                    break;
                }
            if (!found) return std::nullopt;
        }

    int best = component[0];
    PlanarityKey best_key = planarity_rank_key(*candidates[best]);
    for (std::size_t i = 1; i < component.size(); ++i) {
        const PlanarityKey key = planarity_rank_key(*candidates[component[i]]);
        if (key < best_key) {
            best_key = key;
            best = component[i];
        }
    }
    return std::vector<int>{best};
}

// ---- GlobalOptimalStrategy ----
constexpr double kMinQualityThreshold = 0.5;
constexpr int kMaxBacktrack = 30;

bool class_valid(const Cand& c, const std::string& lw) {
    if (lw == "cWW") return c.quality_score >= kMinQualityThreshold && c.has_strong_base_hbond;
    const auto parts = split_pipe(lw);
    bool is_ss = false, is_hh = false;
    for (const auto& p : parts) {
        if (p.size() == 3 && p[1] == 'S' && p[2] == 'S') is_ss = true;
        if (p.size() == 3 && p[1] == 'H' && p[2] == 'H') is_hh = true;
    }
    const bool ss_ok = is_ss && c.num_hbonds >= 2 && c.template_rmsd && *c.template_rmsd <= 1.5;
    const bool hh_ok = is_hh && c.template_rmsd && c.quality_score > 0.0;
    return c.has_strong_base_hbond || ss_ok || hh_ok;
}

std::vector<Cand*> global_optimal(const std::vector<Cand*>& input, double min_score) {
    std::vector<Cand*> valid;
    for (Cand* c : input) {
        if (!c->validation.is_valid || c->quality_score < min_score) continue;
        if (same_chain_adjacent(c->res_id1, c->res_id2)) continue;
        if (class_valid(*c, c->lw_class)) valid.push_back(c);
        else if (!c->precorr_lw.empty() && class_valid(*c, c->precorr_lw)) valid.push_back(c);
    }
    const int n = static_cast<int>(valid.size());
    if (n == 0) return {};

    // Face-aware conflict graph.
    std::vector<std::vector<int>> conflicts(n);
    std::vector<SlotSet> slots(n);
    for (int i = 0; i < n; ++i)
        for (const auto& [s1, s2] : possible_face_slots(valid[i]->lw_class, valid[i]->res_id1,
                                                        valid[i]->res_id2)) {
            slots[i].insert(s1);
            slots[i].insert(s2);
        }
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            bool overlap = false;
            for (const auto& s : slots[i])
                if (slots[j].count(s)) { overlap = true; break; }
            if (overlap) { conflicts[i].push_back(j); conflicts[j].push_back(i); }
        }

    // Connected components (BFS).
    std::vector<char> visited(n, 0);
    std::vector<int> best_indices;
    for (int s = 0; s < n; ++s) {
        if (visited[s]) continue;
        std::vector<int> comp;
        std::vector<int> q = {s};
        visited[s] = 1;
        std::size_t head = 0;
        while (head < q.size()) {
            const int node = q[head++];
            comp.push_back(node);
            for (int nb : conflicts[node])
                if (!visited[nb]) { visited[nb] = 1; q.push_back(nb); }
        }
        if (comp.size() == 1) { best_indices.push_back(comp[0]); continue; }

        if (auto planarity_pick = resolve_shared_residue_conflict(valid, comp, conflicts)) {
            best_indices.push_back((*planarity_pick)[0]);
            continue;
        }

        // Local conflict graph.
        const int m = static_cast<int>(comp.size());
        std::unordered_map<int, int> g2l;
        for (int l = 0; l < m; ++l) g2l[comp[l]] = l;
        std::vector<std::vector<int>> lc(m);
        std::vector<double> scores(m);
        for (int l = 0; l < m; ++l) {
            scores[l] = valid[comp[l]]->quality_score;
            for (int gj : conflicts[comp[l]])
                if (g2l.count(gj)) lc[l].push_back(g2l[gj]);
        }
        std::vector<int> chosen;
        if (m <= kMaxBacktrack) {
            std::vector<double> suffix(m + 1, 0.0);
            for (int i = m - 1; i >= 0; --i) suffix[i] = suffix[i + 1] + scores[i];
            std::vector<int> best_sel, cur;
            double best_score = 0.0;
            std::vector<char> excluded(m, 0);
            std::function<void(int, double)> search = [&](int idx, double cur_score) {
                if (cur_score + suffix[idx] <= best_score) return;
                if (idx == m) {
                    if (cur_score > best_score) { best_sel = cur; best_score = cur_score; }
                    return;
                }
                search(idx + 1, cur_score);  // skip
                if (!excluded[idx]) {
                    std::vector<int> newly;
                    for (int j : lc[idx])
                        if (!excluded[j]) { excluded[j] = 1; newly.push_back(j); }
                    cur.push_back(idx);
                    search(idx + 1, cur_score + scores[idx]);
                    cur.pop_back();
                    for (int j : newly) excluded[j] = 0;
                }
            };
            search(0, 0.0);
            chosen = best_sel;
        } else {
            std::vector<int> order(m);
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(),
                             [&](int a, int b) { return scores[a] > scores[b]; });
            std::vector<char> excl(m, 0);
            for (int l : order) {
                if (excl[l]) continue;
                chosen.push_back(l);
                excl[l] = 1;
                for (int nb : lc[l]) excl[nb] = 1;
            }
        }
        for (int l : chosen) best_indices.push_back(comp[l]);
    }

    // Resolve ambiguity in score order, building used from scratch.
    std::vector<Cand*> selected;
    for (int i : best_indices) selected.push_back(valid[i]);
    std::stable_sort(selected.begin(), selected.end(),
                     [](const Cand* a, const Cand* b) { return a->quality_score > b->quality_score; });
    SlotSet used;
    for (Cand* c : selected) {
        auto s = resolve_for_selection(*c, used);
        if (s) used.insert(s->begin(), s->end());
    }
    return selected;
}

}  // namespace

std::vector<ScoredCandidate> select_pairs(std::vector<ScoredCandidate> candidates,
                                          const RNAChains& chains, double min_score) {
    std::vector<ScoredCandidate>& work = candidates;  // stable addresses
    std::vector<Cand*> valid;
    for (auto& c : work)
        if (c.validation.is_valid && !is_stacking_signature(c)) valid.push_back(&c);

    std::vector<Cand*> cww;
    for (Cand* c : valid)
        if (c->lw_class == "cWW" && kCanonicalSeq.count(sequence(*c))) cww.push_back(c);

    std::vector<Cand*> phase1 = cww_helix_phase(cww, chains);
    std::unordered_set<const Cand*> phase1_set(phase1.begin(), phase1.end());

    SlotSet used;
    for (Cand* c : phase1) {
        auto s = resolve_for_selection(*c, used);
        if (s) used.insert(s->begin(), s->end());
    }

    std::vector<Cand*> phase2_input;
    for (Cand* c : valid) {
        if (phase1_set.count(c)) continue;
        auto needed = candidate_slots(*c);
        if (!needed) continue;
        bool blocked = false;
        for (const auto& s : *needed)
            if (used.count(s)) { blocked = true; break; }
        if (blocked) continue;
        phase2_input.push_back(c);
    }

    std::vector<Cand*> phase2 = global_optimal(phase2_input, min_score);
    for (Cand* c : phase2) {
        auto s = resolve_for_selection(*c, used);
        if (s) used.insert(s->begin(), s->end());
    }

    std::vector<ScoredCandidate> out;
    out.reserve(phase1.size() + phase2.size());
    for (Cand* c : phase1) out.push_back(*c);
    for (Cand* c : phase2) out.push_back(*c);
    return out;
}

}  // namespace pairfinder::algorithms::selection
