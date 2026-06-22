/**
 * @file main.cpp
 * @brief pairfinder CLI entry point.
 *
 * Production path: resolve the input (a 4-char PDB id fetched from RCSB, or a
 * local .cif/.pdb(.gz)), run the find + classify + score pipeline, and emit the
 * find_pairs.py-parity JSON. The hidden `dump-*` subcommands used by the Python
 * differential tests live in debug_commands.cpp; resource-path resolution lives
 * in cli_support.cpp.
 */
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/rigid_body.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/geometry/vector3d.hpp>
#include <pairfinder/pipeline.hpp>
#include <pairfinder/scoring/scorer.hpp>

#include "support/cli_support.hpp"
#include "differential/differential_commands.hpp"

namespace {
constexpr const char* kVersion = "0.1.0";

// Maps raw parameter keys to human-readable report group names.
// Mirrors Python REPORT_GROUPS in pair_finder.scoring.issues.
const std::unordered_map<std::string, std::string> kReportGroup = {
    {"shear",                 "Pair In-Plane Displacement"},
    {"stretch",               "Pair In-Plane Displacement"},
    {"stagger",               "Pair Non-Coplanarity"},
    {"buckle",                "Pair Non-Coplanarity"},
    {"propeller",             "Propeller Twist"},
    {"opening",               "Pair Opening"},
    {"distance",              "H-Bond Distance"},
    {"hbond_angles",          "H-Bond Angles"},
    {"incorrect_hbond_count", "Incorrect H-Bond Count"},
};

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
int find_json(const std::string& path, bool score, bool details, bool emit_classified,
              const std::string& out_path) {
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

    Pipeline pipeline(pfcli::cli_pipeline_config());
    lap("load resources");
    auto prep = pipeline.prepare(path);
    lap("prepare (parse + frames + pipeline_view)");
    auto scored = pipeline.classify(prep);
    lap("score_candidates (find+classify)");

    int candidates_total = static_cast<int>(scored.size());
    int candidates_valid = 0;
    for (const auto& sc : scored)
        if (sc.validation.is_valid) ++candidates_valid;

    // Optionally capture every valid classified candidate (pre-selection) so a
    // downstream comparison can measure identification recall — whether a pair was
    // found + classified at all, independent of helix-priority selection.
    json classified = json::array();
    if (emit_classified) {
        for (const auto& sc : scored) {
            if (!sc.validation.is_valid) continue;
            if (sc.lw_class.empty() && sc.lw_class_display.empty()) continue;
            classified.push_back(json{
                {"res_id1", sc.res_id1},
                {"res_id2", sc.res_id2},
                {"lw_class", sc.lw_class_display.empty() ? sc.lw_class : sc.lw_class_display},
                {"pair_category", sc.pair_category},
                {"is_ambiguous", sc.is_ambiguous}});
        }
    }

    const auto selected = pipeline.select(std::move(scored), prep.structure);
    lap("rna_chains + selection");

    // Empirical scoring (per-pair score/issues + structure aggregates).
    scoring::StructureScore ss;
    std::unordered_map<std::string, const scoring::PairScore*> score_by_key;
    if (score) {
        auto scorer = pipeline.make_scorer();
        algorithms::hbond::HBondFinder sfinder(pipeline.chem());
        ss = scorer.score_structure(selected, prep.structure, sfinder, pipeline.chem(),
                                    pipeline.typing());
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
            p["classification"] = json{
                {"lw_class", sc.lw_class_display.empty() ? sc.lw_class : sc.lw_class_display},
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
                // Issues as a list of key strings — mirrors Python's
                // [ip.issue for ip in ps.issues] in finder.py.
                json iss = json::array();
                // Per-issue Černý provenance (ProSco / Z' / tier) for downstream
                // highlighting; same order as "issues" (descending penalty).
                json details = json::array();
                for (const auto& ip : it->second->issues) {
                    iss.push_back(ip.issue);
                    details.push_back(json{
                        {"issue", ip.issue},
                        {"severity", round_to(ip.severity, 4)},
                        {"prosco", ip.prosco ? json(round_to(*ip.prosco, 4)) : json(nullptr)},
                        {"zprime", ip.zprime ? json(round_to(*ip.zprime, 4)) : json(nullptr)},
                        {"tier", ip.tier}});
                }
                p["issues"] = std::move(iss);
                p["issue_details"] = std::move(details);
            } else {
                p["score"] = nullptr;
                p["issues"] = json::array();
                p["issue_details"] = json::array();
            }
        }
        pairs.push_back(std::move(p));
    }
    out["pairs"] = std::move(pairs);
    if (emit_classified) out["classified"] = std::move(classified);

    if (score) {
        out["overall_score"] = round_to(ss.overall, 2);
        out["pairs_score"] = round_to(ss.pairs_score, 2);
        out["residues_score"] = round_to(ss.residues_score, 2);
        if (details) {
            // Group by human-readable category names — mirrors Python issue_summary().
            json summary = json::object();
            for (const auto& ps : ss.pair_scores) {
                for (const auto& ip : ps.issues) {
                    const auto git = kReportGroup.find(ip.issue);
                    const std::string& grp = (git != kReportGroup.end()) ? git->second : ip.issue;
                    summary[grp] = summary.value(grp, 0) + 1;
                }
            }
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

// Hidden differential-test subcommands keyed by name. Each takes a single
// <input.pdb> path; classify-suite (7 angles) is handled apart.
using DumpFn = int (*)(const std::string&);
const std::pair<const char*, DumpFn> kDumpCommands[] = {
    {"dump-parse", pfcli::dump_parse},
    {"dump-frames", pfcli::dump_frames},
    {"dump-hbonds", pfcli::dump_hbonds},
    {"dump-candidates", pfcli::dump_candidates},
    {"dump-lwclass", pfcli::dump_lwclass},
    {"dump-classified", pfcli::dump_classified},
    {"dump-scored", pfcli::dump_scored},
    {"dump-pairs", pfcli::dump_pairs},
    {"dump-scores", pfcli::dump_scores},
    {"dump-structure", pfcli::dump_structure},
    {"dump-torsions", pfcli::dump_torsions},
};

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

    // Hidden dump-* commands: <command> <input.pdb>.
    for (const auto& [name, fn] : kDumpCommands) {
        if (args[0] == name) {
            if (args.size() < 2) {
                std::cerr << args[0] << ": missing <input.pdb>\n";
                return 1;
            }
            return fn(args[1]);
        }
    }
    if (args[0] == "classify-suite") {
        if (args.size() < 8) { std::cerr << "classify-suite: need 7 angles\n"; return 1; }
        return pfcli::classify_suite(std::vector<std::string>(args.begin() + 1, args.begin() + 8));
    }

    // Default: `pairfinder <pdb-id | structure-file> [options]` runs the full
    // pipeline and emits the find_pairs.py-parity JSON. A bare PDB id is fetched
    // from RCSB (cached) unless an existing file path is given.
    const std::string input = args[0];
    bool score = true, details = false, allow_download = true, emit_classified = false;
    std::string out_path, cache_override;
    for (std::size_t i = 1; i < args.size(); ++i) {
        if (args[i] == "--no-score") score = false;
        else if (args[i] == "--score") score = true;
        else if (args[i] == "--details") details = true;
        else if (args[i] == "--classified") emit_classified = true;
        else if (args[i] == "--no-download") allow_download = false;
        else if (args[i] == "--out" && i + 1 < args.size()) out_path = args[++i];
        else if (args[i] == "--cache-dir" && i + 1 < args.size()) cache_override = args[++i];
        else { std::cerr << "unknown option: " << args[i] << "\n"; return 1; }
    }
    try {
        const std::filesystem::path file =
            resolve_input(input, cache_dir(cache_override), allow_download);
        return find_json(file.string(), score, details, emit_classified, out_path);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
