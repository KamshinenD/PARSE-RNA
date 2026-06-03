/**
 * @file base_typing.hpp
 * @brief Modified-base -> parent normalization and glycosidic-N selection.
 *
 * Faithful port of the Python core ``normalize_base_type`` / ``is_purine`` /
 * ``is_pyrimidine`` / ``get_glycosidic_n`` (core/residue.py + constants.py).
 * Lookup order: standard {A,G,C,U,T,I,P} -> DNA prefix (DA->A) -> curated
 * MODIFIED_BASE_MAP -> extended map built from ligand_hbond_db.json (parent
 * field, with atom-signature inference for parent=null) -> unchanged.
 *
 * Distinct from core/nucleotide_registry (modified_nucleotides.json), which
 * drives standard frame-template selection; this one matches the *Python*
 * pipeline's typing used by candidate finding, H-bonds and classification.
 */
#ifndef PAIRFINDER_CORE_BASE_TYPING_HPP
#define PAIRFINDER_CORE_BASE_TYPING_HPP

#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>

namespace pairfinder::core {

class BaseTyping {
public:
    /// Builds the curated maps and extends from the ligand DB at @p ligand_db_path
    /// (silently skipped if unreadable, matching the Python import-time behavior).
    explicit BaseTyping(const std::filesystem::path& ligand_db_path);

    /// Normalize a residue name to its parent base (single letter for known
    /// nucleotides; original name for unrecognized ligands).
    std::string normalize(std::string_view residue_name) const;

    bool is_purine(std::string_view base_type) const;
    bool is_pyrimidine(std::string_view base_type) const;

    /// Glycosidic-N atom name (N9 purine / N1 pyrimidine / C5 pseudouridine),
    /// or "" if unknown. Mirrors GLYCOSIDIC_N with the normalize fallback.
    std::string glycosidic_n_name(std::string_view base_type) const;

private:
    std::unordered_map<std::string, std::string> extended_;  ///< from ligand DB
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_BASE_TYPING_HPP
