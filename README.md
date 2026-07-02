# pairfinder

A fast, self-contained C++ tool that reads an RNA structure file (mmCIF or PDB),
**finds** base pairs, **classifies** them (Leontis–Westhof), and **scores** their
geometric quality — emitting per-pair and whole-structure quality as JSON.

A full 70S ribosome is annotated and scored in **under a second**.

---

## 1. Requirements

You need a C++ toolchain, CMake, git, and an internet connection for the first
build. No databases or other tools to install.

| Requirement | Notes |
|---|---|
| **C++20 compiler** | Apple Clang ≥ 14, Clang ≥ 13, or GCC ≥ 10 |
| **CMake ≥ 3.28** | check with `cmake --version` |
| **git** | used to fetch a dependency during the first build |
| **curl** | used at runtime to download PDB entries from RCSB (preinstalled on macOS and most Linux) |
| **Internet** | first build downloads & builds GEMMI (the structure parser); the tool also downloads PDB entries at runtime unless you pass a local file |

All other dependencies are bundled in the repo (`nlohmann/json`, the reference
scoring data) or fetched automatically by CMake (GEMMI). There is nothing to
`apt install` / `brew install` beyond a compiler + CMake.

---

## 2. Build

```bash
# 1. clone
git clone https://github.com/KamshinenD/pair-finder-cpp.git
cd pair-finder-cpp

# 2. configure (Release = optimized)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# 3. build (uses all cores)
cmake --build build -j
```

The first configure/build takes a few minutes because it fetches and compiles
GEMMI; later builds are fast. When it finishes, the executable is at
`build/pairfinder`, and the reference scoring data is staged next to it
automatically — the tool finds it with no extra setup.

---

## 3. Run

The basic form is `pairfinder <input> [options]`, where `<input>` is either a
4-character **PDB ID** (downloaded from RCSB and cached) or a path to a local
`.cif` / `.pdb` (`.gz` accepted).

```bash
# by PDB ID (downloaded + cached under ~/.cache/pairfinder)
./build/pairfinder 1EHZ

# write the JSON to a file (see output rules below)
./build/pairfinder 1EHZ --out 1EHZ.json

# from a local file (no download)
./build/pairfinder path/to/structure.cif --out structure.json

# a full ribosome — still sub-second
./build/pairfinder 7K00 --out 7K00.json
```

### Where the output goes

| you run | JSON goes to |
|---|---|
| no `--out` | **stdout** (redirect with `> file.json` if you like) |
| `--out name.json` (bare filename) | **`pairs/name.json`** — collected in a `pairs/` folder (created automatically) |
| `--out path/to/name.json` (has a directory) | exactly that path (created if needed) |

Downloaded structures are cached in **`~/.cache/pairfinder/`** — nothing is
written into your working directory unless you ask for it.

### Options

| Option | Effect |
|---|---|
| `--out FILE` | write JSON to `FILE` (bare name → `pairs/`; default: stdout) |
| `--details` | include extra per-pair fields (template RMSD, score breakdown) |
| `--no-score` | find + classify only, skip quality scoring |
| `--no-download` | never fetch from RCSB; only use a local file |
| `--cache-dir DIR` | where downloaded structures are cached (default `$PAIRFINDER_CACHE_DIR` or `~/.cache/pairfinder`) |

---

## 4. Output

A JSON object per structure:

```jsonc
{
  "pdb_id": "1EHZ",
  "candidates_total": 1234,      // pairs considered
  "candidates_valid": 200,       // passed geometric validation
  "n_pairs": 31,                 // final selected base pairs
  "overall_score": 88.36,        // whole-RNA quality (0–100)
  "pairs_score": 98.46,          // base-pair geometry component
  "residues_score": 69.57,       // backbone (suiteness) component
  "pairs": [
    {
      "res_id1": "A-G-1", "res_id2": "A-C-72",
      "sequence": "GC",
      "pair_category": "base-base",
      "classification": { "lw_class": "cWW", "is_ambiguous": false },
      "geometry":   { "dorg": 2.40, "d_v": 0.11, "plane_angle": 24.08, "dNN": 8.42 },
      "rigid_body": { "shear": 2.22, "stretch": 0.89, "stagger": -0.11, "buckle": -12.78, "propeller": 20.40, "opening": 0.33 },
      "num_hbonds": 2,
      "score": 95.12              // per-pair quality (0–100)
    }
    // ...
  ]
}
```

- **`overall_score` = 0.65 · `pairs_score` + 0.35 · `residues_score`.**
- Per-pair `score` and per-structure scores are 0–100 (higher = better geometry).
- A `classification.lw_class` containing `|` (e.g. `cWH|cWW`) with
  `is_ambiguous: true` means the edge assignment is genuinely undetermined
  (the same cases DSSR marks with a `.`).

---

## 5. Optional: PyMOL viewer

If you have **PyMOL** installed, you can open a structure, score it, and
highlight the lowest-quality pairs in one command:

```bash
./integrations/pymol/parse-view 1GID        # a PDB id
./integrations/pymol/parse-view file.cif    # or a local file
```

This also writes the scored JSON to `pairs/1GID_pairs.json`. PyMOL is only
needed for this viewer — the core `pairfinder` tool does not require it.

Once PyMOL is open, an interactive worklist lets you step through the
lowest-quality pairs (`parse_next` / `parse_prev`), jump to any pair or residue
(`parse_goto`), inspect a clicked residue (`parse_info`), overlay the idealized
target geometry (`parse_ideal`), and dump the scored pairs (`parse_dump`). The
full command reference and keyboard navigation are documented in
[`integrations/pymol/README.md`](integrations/pymol/README.md).

---

## 6. Notes

- **Internet:** needed once (first build, to fetch GEMMI) and at runtime only
  when you give a bare PDB ID (it downloads via `curl` and caches the file).
  Local-file runs need no network — use `--no-download` to enforce that.
- **Cache:** downloaded structures live in `~/.cache/pairfinder` (override with
  `--cache-dir` or `$PAIRFINDER_CACHE_DIR`). Safe to delete anytime.
- **Reference data** (`resources/`) ships with the repo and is staged next to the
  binary at build time; the tool finds it automatically.
- **Performance:** ~4,000 base pairs/second; a 70S ribosome (~2,100 pairs) scores
  in ~0.5 s, the largest RNA structures (~16,000 nt) in ~1.5 s.

---

## 7. Project layout

```
include/pairfinder/   public headers (geometry, core, io, algorithms, scoring)
src/                  implementations
apps/main.cpp         the CLI
resources/            bundled reference data + templates (ship with the binary)
integrations/pymol/   optional PyMOL viewer (parse-view, parse_pymol.py)
CMakeLists.txt        build (fetches GEMMI, vendors nlohmann/json)
```

---

## 8. Troubleshooting

| Symptom | Fix |
|---|---|
| `CMake 3.28 or higher is required` | upgrade CMake (`cmake --version` to check) |
| First build fails fetching GEMMI | check internet / proxy; re-run `cmake -S . -B build` |
| `download failed` / `HTTP 404` at runtime | the PDB ID isn't on RCSB, or no network — pass a local `.cif`/`.pdb`, or check `curl` is installed |
| Want a clean rebuild | `rm -rf build` and repeat the Build steps |
