# PARSE → PyMOL

Highlight problematic base pairs in PyMOL during refinement, using the
`pairfinder` binary. Flagged pairs are shown as sticks on a continuous
yellow→red heat-map by Černý severity; re-run after an edit and a fixed pair
disappears. The original object is never modified.

## Quick start (one command)

From the repo, after building `pairfinder`:

```bash
./integrations/pymol/parse-view 1GID                 # a PDB id (fetched from RCSB)
./integrations/pymol/parse-view /path/to/mymodel.pdb # or a local file
```

This opens PyMOL, loads the plugin, loads the structure, scores it, and
highlights the issues — all in one go. (Override the PyMOL or pairfinder
executable with `PYMOL=...` / `PAIRFINDER_BINARY=...`; `PARSE_VIEW_DRYRUN=1`
prints the PyMOL commands without launching.)

To drive it interactively inside an already-open PyMOL session instead, use the
commands below.

## Setup

1. Build `pairfinder` (see the top-level README).
2. Tell the plugin where the binary is — either:
   - export `PAIRFINDER_BINARY=/path/to/build/pairfinder` before launching PyMOL, or
   - run `parse_set_binary /path/to/build/pairfinder` inside PyMOL.
3. In PyMOL: `run /path/to/integrations/pymol/parse_pymol.py`

## Use

```
parse_score                 # score the current object, highlight problems
parse_score 6DN1            # score object named 6DN1
parse_score 6DN1, 0.5       # only flag issues with severity >= 0.5 (less noise)
parse_score 6DN1, 0, 1      # also zoom onto the flagged pairs
parse_clear                 # remove the overlay
```

Each `parse_score` run:
1. saves the object's **current** coordinates to a temp PDB,
2. runs `pairfinder` on it (sub-second, even on a ribosome),
3. builds a non-destructive overlay object `parse_overlay` containing only the
   flagged pairs (sticks, colored by severity), plus a `parse_problems`
   selection you can `zoom parse_problems` / iterate over,
4. prints a summary to the PyMOL log, ranked by each pair's worst-parameter
   |Z′| (ProSco-gated). The 0–100 quality score is never used or shown — pairs
   are ranked and reported purely by their ProSco/Z′ geometric deviation.

A typical loop: `parse_score` → inspect the red pairs → fix geometry / run a
refinement round → reload → `parse_score` again; the highlights rebuild and
fixed pairs drop out.

## What gets colored

The whole structure is colored **gray** (all Preferred geometry recedes), and the
flagged pairs are drawn as colored sticks on top. Coloring follows the three-tier
scheme (Černý et al., NAR 2026), but the **color is the continuous severity**, not
the discrete tier — so a near-boundary "Allowed" pair (severity ≈ 1) shows up
nearly as red as an "Of Concern" pair. The discrete tier and the actual Z′ appear
in the printed summary.

- **gray** — Preferred geometry (no issue); the neutral background.
- **yellow → red** — increasing severity of the worst ProSco/Z′-scored parameter
  (shear, stretch, stagger, buckle, propeller, opening, H-bond distance/angles).
- **not colored** — `incorrect_hbond_count` (a missing canonical H-bond) is not
  ProSco/Z′-scored, so it has no tier; such pairs appear in the summary only.

Pass `parse_score target, 0, 0, 0` (last arg `gray_background=0`) to keep the
object's existing colors instead of graying. Only colors change — coordinates are
never modified.

## Commands

| command | effect |
|---|---|
| `parse_score [obj], [min_severity], [zoom]` | score + highlight all flagged pairs |
| `parse_list [n]` | print the ranked worklist (default 25) |
| `parse_next` / `parse_prev` | step through the worklist one pair at a time |
| `parse_goto <rank \| res-id>` | jump to a worklist entry (e.g. `parse_goto 12` or `parse_goto 1A-G-169`) |
| `parse_info [selection]` | inspect a clicked/selected residue → jump the worklist to its pair |
| `parse_keys [on\|off]` | bind/unbind the navigation keys |
| `parse_clear` | remove every overlay + reset |
| `parse_set_binary <path>` | set the pairfinder binary path |

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
parse_goto 1A-G-169     # jump straight to a specific pair
```

The queue is **chain-grouped, worst-first within each chain**, so stepping stays in
one region instead of teleporting across the structure. On each stop the rest of the
structure is dimmed (still visible for orientation), and the log prints the failing
parameters (e.g. `shear Z'=+11.3, opening Z'=+10.3` — which tells you *how* it's
distorted) plus a **Phenix-style selection string** for the pair.

The queue is ranked by each pair's **single most-deviant parameter's |Z′|** (the
ProSco-gated Černý severity, *uncapped* so a 6.6σ pair outranks a 5.1σ one) — never
by the 0–100 quality score, which is intentionally absent from all output. Each
entry's headline is that worst parameter, e.g. `propeller 6.2σ`.

Typical loop: `parse_next` to a problem → understand it → fix it in Phenix / edit the
coordinates → reload and re-run `parse_score` → the fixed pair drops off the queue.
