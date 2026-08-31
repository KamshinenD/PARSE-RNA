/**
 * @file clashes.hpp
 * @brief Heavy-atom steric clash detection (mirror of scoring/clashes.py).
 *
 * A pair of heavy atoms overlaps when it sits closer than its minimum
 * acceptable separation:
 *
 *     min_ok  = kHBondFloor  when both atoms are polar (N/O/S)
 *             = r_i + r_j    otherwise
 *     overlap = min_ok - d
 *
 * An overlap is only a CLASH once it exceeds kClashThreshold; below that it is
 * ordinary close packing and carries no penalty. Severity then ramps linearly
 * from the threshold to kClashSeverityCap:
 *
 *     severity = clamp((overlap - kClashThreshold)
 *                      / (kClashSeverityCap - kClashThreshold), 0, 1)
 *
 * so a contact that has only just become a clash costs almost nothing.
 * pair_clash() applies the threshold; worst_overlap() does NOT — it is the
 * feature-extraction entry point, and scoring its return value directly is a
 * bug (it charges score for ordinary van der Waals proximity).
 *
 * Polar pairs get the H-bond floor rather than the van der Waals sum because a
 * donor/acceptor pair is *supposed* to sit inside vdW contact; only below
 * kHBondFloor is the distance physically impossible. Normal-range H-bond
 * distances are already scored by the `distance` parameter.
 *
 * Hydrogens are ignored: modelled in only ~8% of PDB entries, so counting them
 * would make the result depend on whether H were deposited.
 *
 * Performance notes (this is the hot path — the Python version is the
 * reference, this one has to be fast):
 *   - element/radius/polarity are resolved once per residue into a compact
 *     POD array, not re-derived per atom pair;
 *   - the inner loop compares SQUARED distances and never calls sqrt;
 *   - whole-structure detection uses a uniform grid (O(N)), not all-pairs;
 *   - the covalent-linkage test is two O(1) atom lookups, and the bond-graph
 *     walk only runs for the rare pair that is actually linked.
 */
#ifndef PAIRFINDER_SCORING_CLASHES_HPP
#define PAIRFINDER_SCORING_CLASHES_HPP

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::scoring {

/// Overlap (A) at which a contact is called a clash.
constexpr double kClashThreshold = 0.4;
/// Overlap (A) at which severity saturates at 1.0. The ramp runs BETWEEN the
/// two, so severity is 0 at the threshold rather than jumping to 0.8.
///
/// Was 0.5, which both left only a 0.1 A ramp and (because the ramp was not
/// rebased) made the first real clash cost ~80% of the full penalty. 0.9
/// restores a 0.5 A ramp: above the 90th percentile of observed clash overlaps
/// (0.83 good / 0.86 poor), saturating only ~7-9% of clashes, and where the
/// good/poor log-odds stops improving (1.1074; 1.1283 at cap 1.2).
constexpr double kClashSeverityCap = 0.9;
/// Minimum physically possible donor-acceptor separation (A). Real H-bonds
/// bottom out near 2.4-2.6; below 2.2 the geometry is impossible.
constexpr double kHBondFloor = 2.2;
/// Covalent bond length ceiling for the inferred bond graph.
constexpr double kBondMax = 1.95;
/// Contacts within this many covalent bonds are normal geometry, not clashes.
constexpr int kBondDepth = 4;

/// One steric clash between two heavy atoms.
struct ClashContact {
    std::string res_id1, atom1;
    std::string res_id2, atom2;
    double distance = 0.0;
    double overlap = 0.0;
};

/// Severity in [0,1] for a steric overlap. NaN / no overlap -> 0.
double clash_severity(double overlap);
double clash_severity(std::optional<double> overlap);

/// True when `res_name` is a standard or modified nucleotide. Uses BaseTyping
/// (the C++ mirror of Python's _EXTENDED_BASE_MAP) so both sides agree on which
/// residues count as RNA — including modified bases such as PSU/5MC/H2U/YYG.
bool is_rna_residue_name(std::string_view res_name, const core::BaseTyping& typing);
/// True when the residue name embedded in `res_id` is a nucleotide.
bool is_rna_res_id(std::string_view res_id, const core::BaseTyping& typing);

/// True for a base atom (not sugar, not phosphate). Defined negatively: sugar
/// atoms carry a prime, phosphate atoms are P / OP* / O*P. Avoids an enumerated
/// list so modified-base substituents (2MG's CM2, wybutosine's tail) count as
/// base rather than being silently dropped.
bool is_base_atom(std::string_view atom_name);

/// Worst SCORED overlap between two residues (signed, A): base-involved
/// contacts only. Negative when cleanly separated; nullopt when either residue
/// has no heavy atoms. Covalently bonded atom pairs are excluded.
std::optional<double> worst_overlap(const core::Residue& res1, const core::Residue& res2);

/// Worst SCORED clash between two paired residues (>=1 base atom involved),
/// or nullopt when below `threshold`. The penalty is charged to the base pair,
/// so the base must be implicated: base-base, base-sugar and base-phosphate
/// all count.
std::optional<ClashContact> pair_clash(const core::Residue& res1, const core::Residue& res2,
                                       double threshold = kClashThreshold);

/// Worst backbone-only clash between two paired residues (no base atom
/// involved). Report-only: a sugar-sugar / phosphate-phosphate contact is a
/// backbone defect, not a base-pair one.
std::optional<ClashContact> pair_backbone_clash(const core::Residue& res1,
                                                const core::Residue& res2,
                                                double threshold = kClashThreshold);

/// All RNA-RNA heavy-atom clashes in a structure, worst first.
///
/// Restricted to nucleotide-nucleotide contacts by design: validated on 896
/// PDB-REDO pairs, RNA-RNA clashes fall 29.9% under re-refinement while protein
/// contacts rise 189%, because heavy-atom criteria penalise the hydrogen-aware
/// side-chain packing that re-refinement produces.
std::vector<ClashContact> structure_clashes(
    const std::unordered_map<std::string, core::Residue>& residues,
    const core::BaseTyping& typing, double threshold = kClashThreshold);

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_CLASHES_HPP
