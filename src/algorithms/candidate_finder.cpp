/**
 * @file candidate_finder.cpp
 * @brief Pair-candidate finding + geometric validation (port of PairCache).
 */
#include <pairfinder/algorithms/candidate_finder.hpp>

#include <cmath>
#include <unordered_map>

namespace pairfinder::algorithms {

namespace {

constexpr double kPi = 3.14159265358979323846;
using geometry::Vector3d;

double clamp_unit(double x) { return x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x); }

// Single-letter res_name from res_id (matches PairCache._extract_residue_name:
// strip a leading 'D' from 2-char D-prefixed DNA codes; keep everything else).
std::string extract_res_name(const std::string& res_id) {
    const auto p1 = res_id.find('-');
    if (p1 == std::string::npos) return res_id;
    const auto p2 = res_id.find('-', p1 + 1);
    const std::string name =
        res_id.substr(p1 + 1, (p2 == std::string::npos ? res_id.size() : p2) - p1 - 1);
    if (name.size() == 2 && name[0] == 'D') return name.substr(1);
    return name;
}

}  // namespace

ValidationResult validate_pair(const core::ReferenceFrame& f1,
                               const core::ReferenceFrame& f2,
                               const Vector3d& n1n9_1, const Vector3d& n1n9_2,
                               const ValidationThresholds& thr) {
    ValidationResult r;
    const Vector3d dorg_vec = f1.origin - f2.origin;
    r.dorg = dorg_vec.norm();

    r.dir_x = f1.x_axis().dot(f2.x_axis());
    r.dir_y = f1.y_axis().dot(f2.y_axis());
    r.dir_z = f1.z_axis().dot(f2.z_axis());

    // Mean helix axis (get_bp_zoave): add z-axes if aligned, subtract if anti.
    Vector3d zave = r.dir_z > 0 ? f1.z_axis() + f2.z_axis() : f2.z_axis() - f1.z_axis();
    const double zlen = zave.norm();
    zave = zlen > 1e-10 ? zave / zlen : f1.z_axis();
    r.d_v = std::abs(dorg_vec.dot(zave));

    r.plane_angle = std::acos(std::abs(clamp_unit(r.dir_z))) * 180.0 / kPi;
    r.dNN = (n1n9_1 - n1n9_2).norm();
    r.quality_score =
        r.dorg + thr.d_v_weight * r.d_v + r.plane_angle / thr.plane_angle_divisor;

    r.distance_check = r.dorg <= thr.max_dorg;
    r.d_v_check = r.d_v <= thr.max_d_v;
    r.plane_angle_check = r.plane_angle <= thr.max_plane_angle;
    r.dNN_check = r.dNN >= thr.min_dNN;
    r.is_valid =
        r.distance_check && r.d_v_check && r.plane_angle_check && r.dNN_check;
    return r;
}

std::vector<PairCandidate> find_candidates(
    const core::Structure& structure,
    const std::vector<std::pair<std::string, core::ReferenceFrame>>& frames,
    const core::BaseTyping& typing, double max_distance,
    const ValidationThresholds& thr) {
    // res_id -> Residue* for glycosidic-N lookup.
    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    // Resolve glycosidic-N coordinate once per framed residue (nullptr-safe).
    struct Node {
        const std::string* res_id;
        const core::ReferenceFrame* frame;
        Vector3d gly;
        bool has_gly;
    };
    std::vector<Node> nodes;
    nodes.reserve(frames.size());
    for (const auto& [res_id, frame] : frames) {
        Node n{&res_id, &frame, {}, false};
        const auto it = res_by_id.find(res_id);
        if (it != res_by_id.end()) {
            const std::string name = typing.glycosidic_n_name(it->second->base_type());
            if (!name.empty()) {
                if (const core::Atom* a = it->second->get_atom(name)) {
                    n.gly = a->coords;
                    n.has_gly = true;
                }
            }
        }
        nodes.push_back(n);
    }

    std::vector<PairCandidate> out;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        for (std::size_t j = i + 1; j < nodes.size(); ++j) {
            const double dorg = (nodes[i].frame->origin - nodes[j].frame->origin).norm();
            if (dorg > max_distance) continue;  // KDTree query_ball_point radius
            if (!nodes[i].has_gly || !nodes[j].has_gly) continue;
            PairCandidate c;
            c.res_id1 = *nodes[i].res_id;
            c.res_id2 = *nodes[j].res_id;
            c.res_name1 = extract_res_name(c.res_id1);
            c.res_name2 = extract_res_name(c.res_id2);
            c.validation = validate_pair(*nodes[i].frame, *nodes[j].frame, nodes[i].gly,
                                         nodes[j].gly, thr);
            out.push_back(std::move(c));
        }
    }
    return out;
}

}  // namespace pairfinder::algorithms
