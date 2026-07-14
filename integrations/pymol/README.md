# PARSE → PyMOL

Highlight problematic base pairs in PyMOL during refinement, using the
`parse` binary. Flagged pairs are shown as sticks on a continuous
yellow→red heat-map by Černý severity; re-run after an edit and a fixed pair
disappears. The original object is never modified.

## Quick start (one command)

From the repo, after building `parse`:

```bash
./integrations/pymol/parse-view 1GID                 # a PDB id (fetched from RCSB)
./integrations/pymol/parse-view /path/to/mymodel.pdb # or a local file
./integrations/pymol/parse-view 1GID acceptable      # optional tier: review|acceptable|all
```

This opens PyMOL, loads the plugin, loads the structure, scores it, and
highlights the issues — all in one go. (Override the PyMOL or parse
executable with `PYMOL=...` / `PARSE_BINARY=...`; `PARSE_VIEW_DRYRUN=1`
prints the PyMOL commands without launching.)

To drive it interactively inside an already-open PyMOL session instead, use the
commands below.

## Setup

1. Build `parse` (see the top-level README).
2. Tell the plugin where the binary is — either:
   - export `PARSE_BINARY=/path/to/build/parse` before launching PyMOL, or
   - run `parse_set_binary /path/to/build/parse` inside PyMOL.
3. In PyMOL: `run /path/to/integrations/pymol/parse_pymol.py`

## Use

```
parse_score                 # score current object; worklist = Review tier (score<75)
parse_score 6DN1            # score object named 6DN1 (Review tier)
parse_score 6DN1, acceptable  # the minor-issue middle (75<=score<100), a polish pass
parse_score 6DN1, all       # every non-Preferred pair
parse_score 6DN1, all, 0.5  # ...and only issues with color severity >= 0.5 (less noise)
parse_clear                 # remove the overlay
```

Each `parse_score` run:
1. saves the object's **current** coordinates to a temp PDB,
2. runs `parse` on it (sub-second, even on a ribosome),
3. builds a non-destructive overlay object `parse_overlay` containing only the
   flagged pairs (sticks, colored by severity), plus a `parse_problems`
   selection you can `zoom parse_problems` / iterate over,
4. prints a summary to the PyMOL log, ranked by each pair's worst-parameter
   |Z′| (ProSco-gated).

**Which pairs are flagged** is gated by the whole-pair **score tier** (`tier`
arg): `review` (score < 75, the default), `acceptable` (75 ≤ score < 100), or
`all` (both). The score-tier decides *membership* — who a refiner reviews — while
the ProSco/Z′ severity still decides *color and rank* within the worklist. So the
default `parse_score` shows only the genuinely-poor pairs, and `parse_score obj,
acceptable` (or `all`) widens the net for a thorough polishing pass. The printed
header shows the structure's tier counts (e.g. `104 Preferred · 32 Acceptable ·
10 Review`) so you know how big each pass is before you run it.

A typical loop: `parse_score` (Review) → fix the flagged pairs → re-run; when the
Review worklist is empty, `parse_score obj, acceptable` to polish the middle.

## What gets colored

The whole structure is colored **gray** (unflagged pairs recede), and the flagged
pairs are drawn as colored sticks on top. **Hue encodes the *kind* of defect;
within the warm ramp, intensity encodes severity:**

- **gray** — no flagged issue; the neutral background.
- **yellow → red** — a geometric defect, colored by the continuous severity of the
  worst ProSco/Z′-scored parameter (shear, stretch, stagger, buckle, propeller,
  opening, H-bond distance/angles). A near-boundary "Allowed" pair (severity ≈ 1)
  shows up nearly as red as an "Of Concern" one; the discrete Černý tier (NAR 2026)
  and the actual Z′ appear in the printed summary.
- **blue** — a bonding defect: `incorrect_hbond_count` (a missing/extra canonical
  H-bond). It has no |Z′|, so it is kept **off** the warm ramp (a blue pair can be
  geometrically perfect) and is placed at the **bottom** of the worklist. A pair
  with both a geometry and a count issue is colored by its geometry (warm); the
  count still shows in its breakdown.

Pass `parse_score target, review, 0, 0, 0` (last arg `gray_background=0`) to keep
the object's existing colors instead of graying. Only colors change — coordinates
are never modified.

## Commands

| command | effect |
|---|---|
| `parse_score [obj], [tier], [min_severity], [zoom]` | score + highlight the flagged pairs; `tier` = `review` (default) / `acceptable` / `all` |
| `parse_list [n]` | print the ranked worklist (default 25) |
| `parse_dump [path]` | write the current scored pairs (incl. edits) to JSON — default `<obj>_pairs.json` |
| `parse_load <path>` | rebuild the worklist from a `parse_dump` JSON — **no binary, no re-scoring** (for sharing a session) |
| `parse_next` / `parse_prev` | step through the worklist one pair at a time |
| `parse_overview` | zoom back out to the whole-structure overview (undim + re-show all flagged pairs) |
| `parse_goto <rank \| chain+number>` | jump to worklist entry `12`, **or** inspect any residue by `A-169` / `A169` / `A-G-169` — including a clean pair not in the worklist, printing its score + what's out of range |
| `parse_info [selection]` | inspect a clicked/selected residue → jump the worklist to its pair |
| `parse_ideal [on\|off]` | overlay the **idealized** version of the current pair (green ghost = target geometry) so a refiner sees how to fix it; sticky toggle, default off |
| `parse_keys [on\|off]` | bind/unbind the navigation keys |
| `parse_clear` | remove every overlay + reset |
| `parse_set_binary <path>` | set the parse binary path |
| `parse_set_ideals <path>` | set the idealized-template dir (default `resources/basepair-idealized`) |

## Guided navigation (keys, on-screen status, click-to-inspect)

After `parse_score`, the integration becomes a guided triage tool — you don't have to
keep typing commands:

- **Keyboard:** **→** = next pair, **←** = previous, **Ctrl-I** = inspect the atom you
  last clicked. (Keys bind automatically; `parse_keys off` unbinds them.)
- **On-screen status:** a one-line HUD at the bottom of the viewport always shows where
  you are — `PARSE 3/45  A226·U249 (chain A)  cWW AU  propeller 6.2σ [Of Concern]` — so
  you never have to scroll the log.
- **Click-to-inspect:** click any flagged stick in 3D, then `parse_info` (or **Ctrl-I**).
  The worklist jumps to that pair and prints its full breakdown + Phenix selection.
  Click a non-flagged residue and it tells you it's Preferred/unscored.

## Worklist (fix them one at a time)

`parse_score` lights up *every* flagged pair at once — a health map, good for seeing
where the problem regions cluster, but overwhelming on a large structure. To actually
work through them, use the worklist:

```
parse_score target      # build the queue + global map
parse_next              # frame the worst pair: bright sticks, faded 8 A context,
                        #   H-bonds, a label, and its distortion printed to the log
parse_next              # ... the next one, and so on
parse_overview          # zoom back out to the whole structure (undim, all pairs shown)
parse_goto 1A-G-169     # jump straight to a specific pair
```

The queue is **chain-grouped, worst-first within each chain**, so stepping stays in
one region instead of teleporting across the structure. On each stop the rest of the
structure is dimmed (still visible for orientation), and the log prints the failing
parameters (e.g. `shear Z'=+11.3, opening Z'=+10.3` — which tells you *how* it's
distorted) plus a **Phenix-style selection string** for the pair.

The queue is ranked by each pair's **single most-deviant parameter's |Z′|** (the
ProSco-gated Černý severity, *uncapped* so a 6.6σ pair outranks a 5.1σ one) — not
by the 0–100 quality score, which saturates near 100 for most pairs and so makes a
poor rank key. Each entry shows **both**: the pair's score (e.g. `score 65/100`)
*and* its worst parameter (e.g. `propeller 6.2σ`). `parse_overview` additionally
prints the whole-structure score (overall / pairs / backbone).

**Visual priority:** the worklist is a *visual* triage tool, so `hbond_angles` —
a bond-geometry (chemistry) strain that doesn't make a pair *look* distorted — is
down-weighted (×0.3) in the rank and highlight color. A visibly twisted/splayed
pair therefore outranks one whose only large deviation is H-bond angle, and the
headline names the *visible* defect. The true σ is still reported, and H-bond
*distance* keeps full weight (a long/short bond shows as a visibly splayed pair).
This affects only the PyMOL triage order/color — never the 0–100 score or scorer.

Typical loop: `parse_next` to a problem → understand it → fix it in Phenix / edit the
coordinates → reload and re-run `parse_score` → the fixed pair drops off the queue.

## Show the ideal pair (what it should look like)

Knowing a pair is distorted is half the job; a refiner also needs to see the
*target*. `parse_ideal` overlays the **idealized** version of the current pair —
a perfectly planar, H-bond-optimized template for that Leontis–Westhof class and
sequence — as a thin **green ghost** with its ideal H-bonds:

```
parse_next
parse_ideal            # green ghost of the ideal pair appears, fitted on the real one
parse_ideal on         # ... and stays on for every parse_next / parse_goto
parse_ideal off        # hide it again
```

The ideal is **best-fit** onto the actual pair (least-squares over the shared base
atoms of *both* bases), so each base shows how far it sits from where it should be
— flatten/untwist the real (colored) sticks toward the green ghost. Templates live
in `resources/basepair-idealized/<class>/<seq>.pdb`; point elsewhere with
`parse_set_ideals <dir>` or the `PARSE_IDEALS_DIR` env var. It's a sticky toggle,
**off by default** (so the view stays clean until you ask for the target). If no
template exists for a pair's class+sequence, it says so and shows nothing.

## Share a scored session (colleague navigates with no binary)

A `.pse` saves the *picture* but not the plugin or the worklist, so a colleague
can't navigate from it alone. To let someone step through a scored structure
**without the `parse` binary and without re-scoring**, send three files:

1. `1y26.pse`        — the session (`save 1y26.pse` after scoring)
2. `1y26_pairs.json` — the worklist (`parse_dump`)
3. `parse_pymol.py`  — the plugin

They then run, in PyMOL:

```
load 1y26.pse           # the structure + highlights
run parse_pymol.py      # registers the parse_* commands
parse_load 1y26_pairs.json   # rebuilds the worklist from the JSON (no binary)
parse_next              # ... and they can step through exactly as you did
```

`parse_load` reconstructs the queue from the dumped JSON, so `parse_next` / `parse_prev`
/ `parse_list` / `parse_overview` / `parse_dump` all work. Only `parse_score` itself
(scoring a *new* model) needs the binary. Tip: put those three lines in a `view.pml`
so they only have to type `@view.pml`.
