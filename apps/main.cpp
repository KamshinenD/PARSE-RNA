/**
 * @file main.cpp
 * @brief pairfinder CLI entry point.
 *
 * Usage (as the port fills in):
 *     pairfinder <input.pdb> [--out scores.json]   find + classify + score
 *     pairfinder dump-parse <input.pdb>            TSV of parsed atoms (for the
 *                                                  differential check vs Python)
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <pairfinder/algorithms/base_frame_calculator.hpp>
#include <pairfinder/algorithms/candidate_finder.hpp>
#include <pairfinder/algorithms/classify_pair.hpp>
#include <pairfinder/algorithms/confidence.hpp>
#include <pairfinder/algorithms/edge_classifier.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/rigid_body.hpp>
#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>
#include <pairfinder/algorithms/selection.hpp>
#include <pairfinder/algorithms/template_aligner.hpp>
#include <pairfinder/scoring/scorer.hpp>
#include <pairfinder/scoring/torsions.hpp>
#include <pairfinder/scoring/richardson.hpp>
#include <array>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/io/pdb_parser.hpp>

namespace {
namespace resources = pairfinder::resources;
constexpr const char* kVersion = "0.1.0";

// Resource paths default to the bundled resources/ dir (resolved exe-relative);
// each is overridable by its env var for advanced/relocated setups.
std::string env_or(const char* var, const std::filesystem::path& fallback) {
    if (const char* env = std::getenv(var)) return env;
    return fallback.string();
}

// Standard base templates (Atomic_<X>.pdb / Atomic.<x>.pdb).
std::string templates_dir() {
    return env_or("PAIRFINDER_TEMPLATES_DIR", resources::templates_dir());
}

// modified_nucleotides.json (residue typing table).
std::string nucleotides_config() {
    return env_or("PAIRFINDER_NT_CONFIG", resources::config("modified_nucleotides.json"));
}

// Emit "<res_id>\t<rms>\tR00..R22\ttx\tty\ttz" for diffing frames vs ls_fitting.
int dump_frames(const std::string& path) {
    const auto structure = pairfinder::io::parse_pdb(path);
    const pairfinder::core::NucleotideRegistry registry(nucleotides_config());
    pairfinder::algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = pairfinder::algorithms::calculate_all_frames(structure, lib);
    for (const auto& [res_id, f] : frames) {
        std::printf("%s\t%.6f", res_id.c_str(), f.rms_fit);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) std::printf("\t%.6f", f.rotation[i][j]);
        std::printf("\t%.6f\t%.6f\t%.6f\n", f.origin.x, f.origin.y, f.origin.z);
    }
    return 0;
}

// ligand_hbond_db.json (H-bond capacity extensions).
std::string ligand_db() {
    return env_or("PAIRFINDER_LIGAND_DB", resources::config("ligand_hbond_db.json"));
}

// Template directories for LW template-RMSD lookup.
std::string idealized_dir() {
    return env_or("PAIRFINDER_IDEALIZED_DIR", resources::basepair_idealized());
}
std::string exemplar_dir() {
    return env_or("PAIRFINDER_EXEMPLAR_DIR", resources::basepair_exemplars());
}

// Reference-data dir for the empirical scorer.
std::string reference_dir() {
    return env_or("PAIRFINDER_REFERENCE_DIR", resources::reference_dir());
}

// Richardson suites table (backbone suiteness).
std::string richardson_suites() {
    return env_or("PAIRFINDER_RICHARDSON_SUITES", resources::file("richardson_suites.json"));
}

// Full pipeline + empirical pair scoring. Emits per scorable selected pair:
// "<res_a>\t<res_b>\t<score>\t<penalty>" sorted.
int dump_scores(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::classification::TemplateAligner aligner(idealized_dir(), exemplar_dir());
    algorithms::classification::HBondPatterns patterns(idealized_dir(), chem);
    const algorithms::classification::QualityScorer qscorer;
    auto scored = algorithms::classification::score_candidates(
        structure, frames, typing, chem, aligner, patterns, qscorer);
    const auto chains = algorithms::RNAChains::from_structure(structure);
    const auto selected = algorithms::selection::select_pairs(std::move(scored), chains);

    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    const std::string ref = reference_dir();
    scoring::Scorer scorer(ref + "/parameter_distributions.json",
                           ref + "/prosco_distributions.json",
                           ref + "/penalty_weights.json", richardson_suites());
    algorithms::hbond::HBondFinder finder(chem);

    std::vector<std::string> rows;
    for (const auto& sc : selected) {
        if (!scoring::Scorer::is_scorable(sc)) continue;
        const core::Residue* r1 = res_by_id.count(sc.res_id1) ? res_by_id[sc.res_id1] : nullptr;
        const core::Residue* r2 = res_by_id.count(sc.res_id2) ? res_by_id[sc.res_id2] : nullptr;
        if (!r1 || !r2) continue;
        const auto ps = scorer.score_pair(sc, *r1, *r2, finder, chem, typing);
        const std::string ki = sc.res_id1 < sc.res_id2 ? sc.res_id1 : sc.res_id2;
        const std::string kj = sc.res_id1 < sc.res_id2 ? sc.res_id2 : sc.res_id1;
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%.2f\t%.4f", ki.c_str(), kj.c_str(),
                      ps.score, ps.penalty);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Debug: classify a 7-angle suite and print its suiteness.
int classify_suite(const std::vector<std::string>& angles) {
    std::array<double, 7> s{};
    for (int i = 0; i < 7; ++i) s[i] = std::stod(angles[i]);
    pairfinder::scoring::RichardsonClassifier cl(richardson_suites());
    std::printf("%.6f\n", cl.suiteness(s));
    return 0;
}

// Debug: emit "<res_id>\t<alpha>..<zeta>" (suite torsions) for diffing.
int dump_torsions(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);
    const auto tors = scoring::compute_all_torsions(structure);
    auto fmt = [](const std::optional<double>& v) {
        char b[24];
        if (v) std::snprintf(b, sizeof(b), "%.3f", *v);
        else std::snprintf(b, sizeof(b), "None");
        return std::string(b);
    };
    std::vector<std::string> rows;
    for (const auto& [rid, t] : tors)
        rows.push_back(rid + "\t" + fmt(t.alpha) + "\t" + fmt(t.beta) + "\t" + fmt(t.gamma) +
                       "\t" + fmt(t.delta) + "\t" + fmt(t.epsilon) + "\t" + fmt(t.zeta));
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Whole-structure score: "OVERALL\t<overall>\t<pairs>\t<residues>" then
// per-residue "RES\t<res_id>\t<score>" sorted.
int dump_structure(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::classification::TemplateAligner aligner(idealized_dir(), exemplar_dir());
    algorithms::classification::HBondPatterns patterns(idealized_dir(), chem);
    const algorithms::classification::QualityScorer qscorer;
    auto scored = algorithms::classification::score_candidates(
        structure, frames, typing, chem, aligner, patterns, qscorer);
    const auto chains = algorithms::RNAChains::from_structure(structure);
    const auto selected = algorithms::selection::select_pairs(std::move(scored), chains);

    const std::string ref = reference_dir();
    scoring::Scorer scorer(ref + "/parameter_distributions.json",
                           ref + "/prosco_distributions.json",
                           ref + "/penalty_weights.json", richardson_suites());
    algorithms::hbond::HBondFinder finder(chem);
    const auto ss = scorer.score_structure(selected, structure, finder, chem, typing);

    std::printf("OVERALL\t%.2f\t%.2f\t%.2f\n", ss.overall, ss.pairs_score, ss.residues_score);
    std::vector<std::string> rows;
    for (const auto& rs : ss.residue_scores) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "RES\t%s\t%.2f", rs.res_id.c_str(), rs.score);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Full pipeline through selection: emit the SELECTED pairs as
// "<res_a>\t<res_b>\t<lw_class>" sorted, for diffing vs find_pairs().pairs.
int dump_pairs(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::classification::TemplateAligner aligner(idealized_dir(), exemplar_dir());
    algorithms::classification::HBondPatterns patterns(idealized_dir(), chem);
    const algorithms::classification::QualityScorer qscorer;

    auto scored = algorithms::classification::score_candidates(
        structure, frames, typing, chem, aligner, patterns, qscorer);
    const auto chains = algorithms::RNAChains::from_structure(structure);
    const auto selected = algorithms::selection::select_pairs(std::move(scored), chains);

    std::vector<std::string> rows;
    for (const auto& sc : selected) {
        const std::string ki = sc.res_id1 < sc.res_id2 ? sc.res_id1 : sc.res_id2;
        const std::string kj = sc.res_id1 < sc.res_id2 ? sc.res_id2 : sc.res_id1;
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s", ki.c_str(), kj.c_str(),
                      sc.lw_class.empty() ? "None" : sc.lw_class.c_str());
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Full scoring pass (classify + score + O2 correction) over all valid candidates.
// Emits "<res_a>\t<res_b>\t<pair_cat>\t<lw>\t<rmsd>\t<quality>\t<num_hb>\t<strong>\t
//        <display>\t<conf_csv>\t<is_amb>\t<precorr>" sorted by valid base-base candidates.
int dump_scored(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::classification::TemplateAligner aligner(idealized_dir(), exemplar_dir());
    algorithms::classification::HBondPatterns patterns(idealized_dir(), chem);
    const algorithms::classification::QualityScorer qscorer;

    const auto scored = algorithms::classification::score_candidates(
        structure, frames, typing, chem, aligner, patterns, qscorer);

    std::vector<std::string> rows;
    for (const auto& sc : scored) {
        if (!sc.validation.is_valid) continue;
        // Canonical key (sorted res_id) for order-independent comparison.
        const std::string ki = sc.res_id1 < sc.res_id2 ? sc.res_id1 : sc.res_id2;
        const std::string kj = sc.res_id1 < sc.res_id2 ? sc.res_id2 : sc.res_id1;
        std::string rmsd_s = "None";
        if (sc.template_rmsd) {
            char b[32];
            if (std::isinf(*sc.template_rmsd)) std::snprintf(b, sizeof(b), "inf");
            else std::snprintf(b, sizeof(b), "%.4f", *sc.template_rmsd);
            rmsd_s = b;
        }
        std::string conf_csv;
        for (std::size_t k = 0; k < sc.lw_class_confidence.size(); ++k) {
            char cb[16];
            std::snprintf(cb, sizeof(cb), "%.3f", sc.lw_class_confidence[k]);
            if (k) conf_csv += ",";
            conf_csv += cb;
        }
        char buf[640];
        std::snprintf(buf, sizeof(buf),
                      "%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%s\t%s\t%d\t%s", ki.c_str(),
                      kj.c_str(), sc.pair_category.c_str(),
                      sc.lw_class.empty() ? "None" : sc.lw_class.c_str(), rmsd_s.c_str(),
                      sc.quality_score, sc.num_hbonds, sc.has_strong_base_hbond ? 1 : 0,
                      sc.lw_class_display.empty() ? "None" : sc.lw_class_display.c_str(),
                      conf_csv.empty() ? "None" : conf_csv.c_str(), sc.is_ambiguous ? 1 : 0,
                      sc.precorr_lw.empty() ? "None" : sc.precorr_lw.c_str());
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Full classification (geometry + template RMSD + ambiguity) per valid candidate.
// Emits "<res_a>\t<res_b>\t<lw_class>\t<rmsd|None|inf>\t<swapped>" sorted.
int dump_classified(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    const auto cands = algorithms::find_candidates(structure, frames, typing);
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;
    std::unordered_map<std::string, const core::ReferenceFrame*> frame_by_id;
    for (const auto& [rid, f] : frames) frame_by_id[rid] = &f;

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::hbond::HBondFinder finder(chem);
    algorithms::classification::TemplateAligner aligner(idealized_dir(), exemplar_dir());
    algorithms::classification::HBondPatterns patterns(idealized_dir(), chem);
    const algorithms::classification::QualityScorer qscorer;

    std::vector<std::string> rows;
    for (const auto& c : cands) {
        if (!c.validation.is_valid) continue;
        std::string ra = c.res_id1, rb = c.res_id2;
        if (rb < ra) std::swap(ra, rb);
        const core::Residue* r1 = res_by_id[ra];
        const core::Residue* r2 = res_by_id[rb];
        const core::ReferenceFrame* f1 = frame_by_id[ra];
        const core::ReferenceFrame* f2 = frame_by_id[rb];
        if (!r1 || !r2 || !f1 || !f2) continue;
        const auto voters = finder.find_face_voters(*r1, *r2);
        const auto strict = finder.find_between(*r1, *r2);
        const auto res = algorithms::classification::classify_pair(*f1, *f2, *r1, *r2, voters,
                                                                   strict, aligner, typing);
        std::string rmsd_s = "None";
        if (res.template_rmsd) {
            const double v = *res.template_rmsd;
            char rb2[32];
            if (std::isinf(v)) std::snprintf(rb2, sizeof(rb2), "inf");
            else std::snprintf(rb2, sizeof(rb2), "%.4f", v);
            rmsd_s = rb2;
        }
        // Confidence (v3): res order is post-swap (lw_class[1] describes cr1).
        const core::Residue* cr1 = res.swapped ? r2 : r1;
        const core::Residue* cr2 = res.swapped ? r1 : r2;
        const std::string sequence = cr1->base_type() + cr2->base_type();
        const auto conf = algorithms::classification::compute_classification_confidence(
            *cr1, *cr2, sequence, res.lw_class, strict, aligner, patterns, typing);
        std::string conf_csv;
        for (std::size_t k = 0; k < conf.confidences.size(); ++k) {
            char cb[16];
            std::snprintf(cb, sizeof(cb), "%.3f", conf.confidences[k]);
            if (k) conf_csv += ",";
            conf_csv += cb;
        }
        const double quality = algorithms::classification::compute_selection_quality(
            c.validation, sequence, strict, res, qscorer, patterns, typing);
        char buf[512];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%.3f", ra.c_str(),
                      rb.c_str(), res.lw_class.c_str(), rmsd_s.c_str(), res.swapped ? 1 : 0,
                      conf.display_label.c_str(), conf_csv.c_str(), conf.is_ambiguous ? 1 : 0,
                      quality);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Classify each valid candidate (canonical res_id order) and emit
// "<res_a>\t<res_b>\t<lw_class>\t<swapped>" sorted, for diffing the LW classifier.
int dump_lwclass(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    const auto cands = algorithms::find_candidates(structure, frames, typing);

    // Switch residues to pipeline base types (PairCache: DC->C, DT->T) — frames
    // were already computed from the original names.
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;
    std::unordered_map<std::string, const core::ReferenceFrame*> frame_by_id;
    for (const auto& [rid, f] : frames) frame_by_id[rid] = &f;

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::hbond::HBondFinder finder(chem);

    std::vector<std::string> rows;
    for (const auto& c : cands) {
        if (!c.validation.is_valid) continue;
        // Canonical (res_a < res_b) so the classifier's input order matches the oracle.
        std::string ra = c.res_id1, rb = c.res_id2;
        if (rb < ra) std::swap(ra, rb);
        const core::Residue* r1 = res_by_id[ra];
        const core::Residue* r2 = res_by_id[rb];
        const core::ReferenceFrame* f1 = frame_by_id[ra];
        const core::ReferenceFrame* f2 = frame_by_id[rb];
        if (!r1 || !r2 || !f1 || !f2) continue;
        const auto voters = finder.find_face_voters(*r1, *r2);
        const auto res = algorithms::classification::classify_lw_class(*f1, *f2, *r1, *r2,
                                                                       voters, typing);
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%d", ra.c_str(), rb.c_str(),
                      res.lw_class.c_str(), res.swapped ? 1 : 0);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Emit each pair candidate (within 15A) as a sorted TSV row for diffing against
// the Python PairCache: "<res_i>\t<res_j>\t dorg d_v plane dNN dx dy dz qscore valid".
int dump_candidates(const std::string& path) {
    using namespace pairfinder;
    const auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    const auto cands = algorithms::find_candidates(structure, frames, typing);

    std::vector<std::string> rows;
    rows.reserve(cands.size());
    for (const auto& c : cands) {
        const std::string& a = c.res_id1;
        const std::string& b = c.res_id2;
        const std::string& ki = a < b ? a : b;
        const std::string& kj = a < b ? b : a;
        const auto& v = c.validation;
        char buf[640];
        std::snprintf(buf, sizeof(buf),
                      "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.5f\t%.5f\t%.5f\t%d",
                      ki.c_str(), kj.c_str(), v.dorg, v.d_v, v.plane_angle, v.dNN,
                      v.dir_x, v.dir_y, v.dir_z, v.quality_score,
                      v.is_valid ? 1 : 0);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Run find_between on every residue pair (legacy order) and emit each H-bond as
// "<res_i>\t<res_j>\t<donor_res>\t<donor_atom>\t<acc_res>\t<acc_atom>\t<dist>\t
//  <h_idx>\t<lp_idx>\t<align>" — sorted — for diffing against the Python finder.
int dump_hbonds(const std::string& path) {
    using namespace pairfinder;
    const auto structure = io::parse_pdb(path);
    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::hbond::HBondFinder finder(chem);

    std::vector<const core::Residue*> residues;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) residues.push_back(&res);
    // Canonical order (sorted res_id) so find_between's argument order matches
    // the Python differential oracle exactly (selection can depend on it).
    std::sort(residues.begin(), residues.end(),
              [](const core::Residue* a, const core::Residue* b) {
                  return a->res_id() < b->res_id();
              });

    std::vector<std::string> rows;
    for (std::size_t i = 0; i < residues.size(); ++i) {
        for (std::size_t j = i + 1; j < residues.size(); ++j) {
            const auto hbonds = finder.find_between(*residues[i], *residues[j]);
            for (const auto& hb : hbonds) {
                const std::string& a = residues[i]->res_id();
                const std::string& b = residues[j]->res_id();
                const std::string& key_i = a < b ? a : b;
                const std::string& key_j = a < b ? b : a;
                char buf[512];
                std::snprintf(buf, sizeof(buf),
                              "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%.4f",
                              key_i.c_str(), key_j.c_str(), hb.donor_res_id.c_str(),
                              hb.donor_atom.c_str(), hb.acceptor_res_id.c_str(),
                              hb.acceptor_atom.c_str(), hb.distance, hb.h_slot_idx,
                              hb.lp_slot_idx, hb.alignment_score);
                rows.emplace_back(buf);
            }
        }
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// ---- Step 9: full `find` pipeline -> JSON (parity with cli/find_pairs.py) ----
using json = nlohmann::json;

double round_to(double v, int n) {
    const double f = std::pow(10.0, n);
    return std::round(v * f) / f;
}

json vec3_json(const pairfinder::geometry::Vector3d& v) {
    return json::array({round_to(v.x, 4), round_to(v.y, 4), round_to(v.z, 4)});
}

// Filename stem with structure extensions stripped (.gz, then .cif/.pdb/.ent).
std::string clean_pdb_id(const std::string& path) {
    std::string s = std::filesystem::path(path).filename().string();
    auto ends_with = [&](const std::string& ext) {
        if (s.size() < ext.size()) return false;
        std::string tail = s.substr(s.size() - ext.size());
        for (char& c : tail) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return tail == ext;
    };
    if (ends_with(".gz")) s.erase(s.size() - 3);
    for (const char* ext : {".cif", ".pdb", ".ent"})
        if (ends_with(ext)) { s.erase(s.size() - 4); break; }
    return s;
}

json frame_json(const pairfinder::core::ReferenceFrame& f) {
    return json{{"origin", vec3_json(f.origin)},
                {"x_axis", vec3_json(f.x_axis())},
                {"y_axis", vec3_json(f.y_axis())},
                {"z_axis", vec3_json(f.z_axis())}};
}

// Run the full pipeline and emit the find_pairs.py JSON (default + --details).
int find_json(const std::string& path, bool score, bool details, const std::string& out_path) {
    using namespace pairfinder;
    // Per-stage wall-clock when PAIRFINDER_PROFILE is set (-> stderr).
    const bool prof = std::getenv("PAIRFINDER_PROFILE") != nullptr;
    auto t0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* name) {
        if (!prof) return;
        const auto now = std::chrono::steady_clock::now();
        std::fprintf(stderr, "[profile] %-22s %8.3f s\n", name,
                     std::chrono::duration<double>(now - t0).count());
        t0 = now;
    };

    auto structure = io::parse_pdb(path);
    lap("parse");
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    lap("load templates/registry");
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    lap("frames");
    const core::BaseTyping typing(ligand_db());
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);
    lap("typing + pipeline_view");

    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::classification::TemplateAligner aligner(idealized_dir(), exemplar_dir());
    algorithms::classification::HBondPatterns patterns(idealized_dir(), chem);
    const algorithms::classification::QualityScorer qscorer;
    lap("load chem/aligner/patterns");
    auto scored = algorithms::classification::score_candidates(
        structure, frames, typing, chem, aligner, patterns, qscorer);
    lap("score_candidates (find+classify)");

    int candidates_total = static_cast<int>(scored.size());
    int candidates_valid = 0;
    for (const auto& sc : scored)
        if (sc.validation.is_valid) ++candidates_valid;

    const auto chains = algorithms::RNAChains::from_structure(structure);
    lap("rna_chains");
    const auto selected = algorithms::selection::select_pairs(std::move(scored), chains);
    lap("selection");

    // Empirical scoring (per-pair score/issues + structure aggregates).
    scoring::StructureScore ss;
    std::unordered_map<std::string, const scoring::PairScore*> score_by_key;
    if (score) {
        const std::string ref = reference_dir();
        scoring::Scorer scorer(ref + "/parameter_distributions.json",
                               ref + "/prosco_distributions.json",
                               ref + "/penalty_weights.json", richardson_suites());
        algorithms::hbond::HBondFinder sfinder(chem);
        ss = scorer.score_structure(selected, structure, sfinder, chem, typing);
        for (const auto& ps : ss.pair_scores)
            score_by_key[ps.res_id1 + '\x01' + ps.res_id2] = &ps;
        lap("score_structure");
    }

    json out;
    out["pdb_id"] = clean_pdb_id(path);
    out["candidates_total"] = candidates_total;
    out["candidates_valid"] = candidates_valid;
    out["n_pairs"] = static_cast<int>(selected.size());

    json pairs = json::array();
    for (const auto& sc : selected) {
        json p;
        p["res_id1"] = sc.res_id1;
        p["res_id2"] = sc.res_id2;
        p["sequence"] = sc.res_name1 + sc.res_name2;
        p["pair_category"] = sc.pair_category;
        p["num_hbonds"] = sc.num_hbonds;
        p["hbond_categories"] = sc.hbond_categories;

        if (sc.lw_class_display.empty() && sc.lw_class.empty()) {
            p["classification"] = nullptr;
        } else {
            json conf = json::array();
            for (double c : sc.lw_class_confidence) conf.push_back(round_to(c, 2));
            p["classification"] = json{
                {"lw_class", sc.lw_class_display.empty() ? sc.lw_class : sc.lw_class_display},
                {"lw_class_confidence", conf},
                {"is_ambiguous", sc.is_ambiguous}};
        }

        const auto& v = sc.validation;
        p["geometry"] = json{{"dorg", round_to(v.dorg, 4)},
                             {"d_v", round_to(v.d_v, 4)},
                             {"plane_angle", round_to(v.plane_angle, 4)},
                             {"dNN", round_to(v.dNN, 4)}};

        if (sc.frame1 && sc.frame2) {
            const auto rb = algorithms::compute_rigid_body_parameters(*sc.frame1, *sc.frame2);
            p["rigid_body"] = json{{"shear", round_to(rb.shear, 4)},
                                   {"stretch", round_to(rb.stretch, 4)},
                                   {"stagger", round_to(rb.stagger, 4)},
                                   {"buckle", round_to(rb.buckle, 4)},
                                   {"propeller", round_to(rb.propeller, 4)},
                                   {"opening", round_to(rb.opening, 4)}};
            p["frames"] = json{{"res1", frame_json(*sc.frame1)},
                               {"res2", frame_json(*sc.frame2)}};
        }

        if (details) {
            if (!sc.template_rmsd) p["template_rmsd"] = nullptr;
            else if (std::isinf(*sc.template_rmsd)) p["template_rmsd"] = "Infinity";
            else p["template_rmsd"] = round_to(*sc.template_rmsd, 4);
            p["quality_score"] = round_to(sc.quality_score, 4);
        }

        if (score) {
            const auto it = score_by_key.find(sc.res_id1 + '\x01' + sc.res_id2);
            if (it != score_by_key.end()) {
                p["score"] = round_to(it->second->score, 2);
                p["issues"] = it->second->issues;
            } else {
                p["score"] = nullptr;
                p["issues"] = json::array();
            }
        }
        pairs.push_back(std::move(p));
    }
    out["pairs"] = std::move(pairs);

    if (score) {
        out["overall_score"] = round_to(ss.overall, 2);
        out["pairs_score"] = round_to(ss.pairs_score, 2);
        out["residues_score"] = round_to(ss.residues_score, 2);
        if (details) {
            json summary = json::object();
            for (const auto& ps : ss.pair_scores)
                for (const auto& k : ps.issues)
                    summary[k] = summary.value(k, 0) + 1;
            out["score_details"] = json{
                {"pairs_score", round_to(ss.pairs_score, 2)},
                {"residues_score", round_to(ss.residues_score, 2)},
                {"issue_summary", summary},
                {"n_scored_pairs", static_cast<int>(ss.pair_scores.size())},
                {"n_residues_scored", static_cast<int>(ss.residue_scores.size())},
                {"n_skipped_pairs", ss.skipped_pairs}};
        }
    }

    const std::string text = out.dump(2);
    if (out_path.empty()) {
        std::cout << text << '\n';
    } else {
        std::ofstream(out_path) << text << '\n';
    }
    lap("json build + output");
    return 0;
}

// Emit "<res_id>\t<atom>\t%.3f\t%.3f\t%.3f", sorted, for diffing against Python.
int dump_parse(const std::string& path) {
    const auto structure = pairfinder::io::parse_pdb(path);
    std::vector<std::string> rows;
    for (const auto& chain : structure.chains()) {
        for (const auto& res : chain.residues()) {
            for (const auto& atom : res.atoms()) {
                char buf[256];
                std::snprintf(buf, sizeof(buf), "%s\t%s\t%.3f\t%.3f\t%.3f",
                              res.res_id().c_str(), atom.name.c_str(),
                              atom.coords.x, atom.coords.y, atom.coords.z);
                rows.emplace_back(buf);
            }
        }
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// ---- PDB-ID auto-fetch (input acquisition; upstream of the pipeline) ---------

// Cache dir for downloaded structures: --cache-dir, else $PAIRFINDER_CACHE_DIR,
// else ~/.cache/pairfinder, else <tmp>/pairfinder. Persistent + per-user.
std::filesystem::path cache_dir(const std::string& override_dir) {
    if (!override_dir.empty()) return override_dir;
    if (const char* e = std::getenv("PAIRFINDER_CACHE_DIR")) return e;
    if (const char* h = std::getenv("HOME"))
        return std::filesystem::path(h) / ".cache" / "pairfinder";
    return std::filesystem::temp_directory_path() / "pairfinder";
}

// A bare PDB id: short, all-alphanumeric (e.g. 6V3A, 1GID). Distinguishes an id
// from a file path (which has an extension and/or a directory separator).
bool looks_like_pdb_id(const std::string& s) {
    if (s.empty() || s.size() > 8) return false;
    for (char c : s)
        if (!std::isalnum(static_cast<unsigned char>(c))) return false;
    return true;
}

// Ensure a cached structure file for <ID> exists, downloading from RCSB if
// needed. Tries formats in priority order (gzipped mmCIF, then plain mmCIF) so a
// transient/missing variant falls back to the next. Returns the path, or empty
// on failure / when downloads are disabled.
std::filesystem::path fetch_structure(const std::string& id_raw,
                                      const std::filesystem::path& cdir, bool allow_download) {
    std::string id = id_raw;
    for (char& c : id) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    std::error_code ec;
    std::filesystem::create_directories(cdir, ec);

    // mmCIF is the canonical archive format (exists for every released entry,
    // incl. large structures with no legacy .pdb); .cif.gz first, .cif fallback.
    const std::array<const char*, 2> suffixes = {".cif.gz", ".cif"};

    // Cache hit: return the first already-downloaded variant.
    for (const char* suf : suffixes) {
        const std::filesystem::path p = cdir / (id + suf);
        if (std::filesystem::exists(p)) return p;
    }
    if (!allow_download) return {};

    // Download: try each variant until one succeeds.
    for (const char* suf : suffixes) {
        const std::filesystem::path dest = cdir / (id + suf);
        const std::filesystem::path part = dest.string() + ".part";
        const std::string url = "https://files.rcsb.org/download/" + id + suf;
        std::fprintf(stderr, "Fetching %s -> %s\n", url.c_str(), dest.string().c_str());
        const std::string cmd = "curl -fsSL -o '" + part.string() + "' '" + url + "'";
        if (std::system(cmd.c_str()) == 0) {
            std::filesystem::rename(part, dest, ec);
            if (!ec) return dest;
        }
        std::filesystem::remove(part, ec);
    }
    return {};
}

// Resolve the CLI argument to an actual structure file: existing path used as-is;
// otherwise treated as a PDB id and fetched into the cache.
std::filesystem::path resolve_input(const std::string& arg,
                                    const std::filesystem::path& cdir, bool allow_download) {
    if (std::filesystem::exists(arg)) return arg;
    if (looks_like_pdb_id(arg)) {
        auto p = fetch_structure(arg, cdir, allow_download);
        if (!p.empty()) return p;
        throw std::runtime_error("could not obtain PDB id '" + arg +
                                 "' (download failed or --no-download set)");
    }
    throw std::runtime_error("input not found: " + arg +
                             " (give a structure file path or a PDB id)");
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);

    if (args.empty() || args[0] == "-h" || args[0] == "--help") {
        std::cout << "pairfinder " << kVersion << "\n"
                  << "usage:\n"
                  << "  pairfinder <pdb-id | structure-file> [options]\n"
                  << "      find + classify + score -> JSON (stdout or --out FILE)\n"
                  << "      a bare PDB id (e.g. 6V3A) is downloaded from RCSB and cached;\n"
                  << "      an existing .cif/.pdb(.gz) path is used directly.\n"
                  << "  options: --out FILE  --details  --no-score  --no-download  --cache-dir DIR\n"
                  << "  cache:   $PAIRFINDER_CACHE_DIR or ~/.cache/pairfinder\n"
                  << "  pairfinder dump-parse <input.pdb>   # parsed atoms (TSV)\n";
        return args.empty() ? 1 : 0;
    }
    if (args[0] == "--version") {
        std::cout << kVersion << "\n";
        return 0;
    }
    if (args[0] == "dump-parse") {
        if (args.size() < 2) {
            std::cerr << "dump-parse: missing <input.pdb>\n";
            return 1;
        }
        return dump_parse(args[1]);
    }
    if (args[0] == "dump-frames") {
        if (args.size() < 2) {
            std::cerr << "dump-frames: missing <input.pdb>\n";
            return 1;
        }
        return dump_frames(args[1]);
    }
    if (args[0] == "dump-hbonds") {
        if (args.size() < 2) {
            std::cerr << "dump-hbonds: missing <input.pdb>\n";
            return 1;
        }
        return dump_hbonds(args[1]);
    }
    if (args[0] == "dump-candidates") {
        if (args.size() < 2) {
            std::cerr << "dump-candidates: missing <input.pdb>\n";
            return 1;
        }
        return dump_candidates(args[1]);
    }
    if (args[0] == "dump-lwclass") {
        if (args.size() < 2) {
            std::cerr << "dump-lwclass: missing <input.pdb>\n";
            return 1;
        }
        return dump_lwclass(args[1]);
    }
    if (args[0] == "dump-classified") {
        if (args.size() < 2) {
            std::cerr << "dump-classified: missing <input.pdb>\n";
            return 1;
        }
        return dump_classified(args[1]);
    }
    if (args[0] == "dump-scored") {
        if (args.size() < 2) {
            std::cerr << "dump-scored: missing <input.pdb>\n";
            return 1;
        }
        return dump_scored(args[1]);
    }
    if (args[0] == "dump-pairs") {
        if (args.size() < 2) {
            std::cerr << "dump-pairs: missing <input.pdb>\n";
            return 1;
        }
        return dump_pairs(args[1]);
    }
    if (args[0] == "dump-scores") {
        if (args.size() < 2) {
            std::cerr << "dump-scores: missing <input.pdb>\n";
            return 1;
        }
        return dump_scores(args[1]);
    }
    if (args[0] == "dump-structure") {
        if (args.size() < 2) {
            std::cerr << "dump-structure: missing <input.pdb>\n";
            return 1;
        }
        return dump_structure(args[1]);
    }
    if (args[0] == "classify-suite") {
        if (args.size() < 8) { std::cerr << "classify-suite: need 7 angles\n"; return 1; }
        return classify_suite(std::vector<std::string>(args.begin() + 1, args.begin() + 8));
    }
    if (args[0] == "dump-torsions") {
        if (args.size() < 2) {
            std::cerr << "dump-torsions: missing <input.pdb>\n";
            return 1;
        }
        return dump_torsions(args[1]);
    }

    // Default: `pairfinder <pdb-id | structure-file> [options]` runs the full
    // pipeline and emits the find_pairs.py-parity JSON. A bare PDB id is fetched
    // from RCSB (cached) unless an existing file path is given.
    const std::string input = args[0];
    bool score = true, details = false, allow_download = true;
    std::string out_path, cache_override;
    for (std::size_t i = 1; i < args.size(); ++i) {
        if (args[i] == "--no-score") score = false;
        else if (args[i] == "--score") score = true;
        else if (args[i] == "--details") details = true;
        else if (args[i] == "--no-download") allow_download = false;
        else if (args[i] == "--out" && i + 1 < args.size()) out_path = args[++i];
        else if (args[i] == "--cache-dir" && i + 1 < args.size()) cache_override = args[++i];
        else { std::cerr << "unknown option: " << args[i] << "\n"; return 1; }
    }
    try {
        const std::filesystem::path file =
            resolve_input(input, cache_dir(cache_override), allow_download);
        return find_json(file.string(), score, details, out_path);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
