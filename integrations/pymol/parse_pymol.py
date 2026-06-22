"""PARSE -> PyMOL: highlight problematic base pairs during refinement.

Load in PyMOL:    run /path/to/parse_pymol.py
Then:             parse_score            # score the current object, highlight problems
                  parse_score 6DN1, 0.3  # object 6DN1, only severity >= 0.3
                  parse_clear            # remove the highlights

What it does each run:
  1. saves the object's nucleotides to a temp PDB (chains renamed to unique
     single characters so big multi-chain structures -- e.g. a ribosome with
     chain "1A" -- survive the legacy-PDB round-trip; res_ids are mapped back
     to the real chain names afterwards),
  2. runs the `pairfinder` binary on it (sub-second, even on a ribosome),
  3. colors the whole structure gray (so Preferred geometry recedes), then builds
     an overlay object `parse_overlay` of just the flagged base pairs, shown as
     sticks on a continuous yellow->red heat-map by Cerny severity (mild
     Allowed -> yellow, Of-Concern -> red),
  4. prints a ranked summary (incl. H-bond-count issues, which are summary-only).

Re-run after an edit/refinement round and the overlay rebuilds from scratch, so a
pair that's been fixed simply disappears. Only colors change — coordinates are
never modified (pass gray_background=0 to keep the object's existing colors).

Binary: set $PAIRFINDER_BINARY, or `parse_set_binary /path/to/pairfinder`.
Tiers follow Cerny et al. (NAR 2026): Of Concern |Z'|>5, Allowed |Z'|<=5,
Preferred (unflagged). H-bond-count issues carry no ProSco/Z' tier and are
reported in the summary only, not colored.
"""
import json
import os
import shutil
import string
import subprocess
import tempfile

# ---------------------------------------------------------------------------
# Pure logic (importable / testable without PyMOL)
# ---------------------------------------------------------------------------

#: issues that are not ProSco/Z'-scored -> excluded from coloring (summary only)
SUMMARY_ONLY = {"incorrect_hbond_count"}


def parse_res_id(res_id):
    """'A-5MU-54' -> ('A', '5MU', '54'); handles modified names and negative seq.

    Format is chain-resname-resseq. Only resseq may carry a sign, and resname
    never contains '-', so the first '-' splits chain, and the first '-' of the
    remainder splits resname from resseq (which keeps any leading '-').
    """
    i1 = res_id.index("-")
    chain = res_id[:i1]
    rest = res_id[i1 + 1:]
    i2 = rest.index("-")
    return chain, rest[:i2], rest[i2 + 1:]


def residue_selection(res_id):
    """PyMOL selection string for one PARSE res_id (negative resi escaped)."""
    chain, _name, seq = parse_res_id(res_id)
    resi = ("\\" + seq) if seq.startswith("-") else seq
    chain_tok = f"chain {chain} and " if chain else ""
    return f"({chain_tok}resi {resi})"


def colored_issues(issue_details):
    """ProSco/Z'-scored issues only (drops summary-only ones like hbond count)."""
    return [d for d in (issue_details or []) if d.get("issue") not in SUMMARY_ONLY]


_ZPRIME_THRESHOLD = 5.0    # Cerny: |Z'| = 5 is the Allowed/Of-Concern boundary
_PROSCO_PREFERRED = 5.0    # Cerny: ProSco >= 5 is Preferred (not flagged)


def rank_severity(detail):
    """Unclamped Cerny severity of one parameter, for ranking: |Z'| / 5.

    Mirrors the scorer's `min(1, |Z'|/5)` (when ProSco < 5) but WITHOUT the cap,
    so pairs don't all saturate to 1.0 -- a 6.6-sigma deviation outranks a
    5.1-sigma one. ProSco gates whether the parameter counts at all; |Z'| is the
    continuous magnitude. Returns 0 for Preferred parameters (ProSco >= 5).
    """
    ps, zp = detail.get("prosco"), detail.get("zprime")
    if ps is None or ps >= _PROSCO_PREFERRED:
        return 0.0
    if zp is None:
        return 1.0   # outside the distribution, no Z' magnitude available
    return abs(zp) / _ZPRIME_THRESHOLD


def worst_parameter(issue_details):
    """The colored issue with the largest unclamped |Z'| severity, or None."""
    ci = colored_issues(issue_details)
    if not ci:
        return None
    return max(ci, key=rank_severity)


def severity_to_rgb(sev):
    """Continuous yellow(0)->red(1) heat-map; clamps to [0,1]."""
    s = 0.0 if sev is None else max(0.0, min(1.0, float(sev)))
    return [1.0, 1.0 - s, 0.0]


def pair_label(f):
    """Compact, readable name for a flagged pair, e.g. 'A226.U249 (chain A)'.

    The raw res_ids ('A-A-235' / 'A-U-239') run together when joined, so build
    'resnameSeq.resnameSeq' and append the chain (or both chains if they differ).
    """
    c1, n1, s1 = parse_res_id(f["res_id1"])
    c2, n2, s2 = parse_res_id(f["res_id2"])
    if c1 == c2:
        return f"{n1}{s1}·{n2}{s2} (chain {c1})"
    return f"{c1}/{n1}{s1}·{c2}/{n2}{s2}"


def severity_headline(f):
    """'propeller 6.2sigma' -- the worst parameter and its |Z'| magnitude."""
    issue = f.get("worst_issue") or "?"
    z = f.get("worst_z")
    return f"{issue} {abs(z):.1f}σ" if z is not None else f"{issue} (no Z')"


def flagged_pairs(data, min_severity=0.0):
    """Pairs with a colorable issue, ranked worst-first by unclamped |Z'|.

    Ranking is the single most-deviant parameter's |Z'| (ProSco-gated), never the
    0-100 quality score. Returns dicts:
    {res_id1, res_id2, sequence, lw_class, rank, sev, worst_issue, worst_z,
     tier, issue, details}.  `rank` sorts; `sev` (clamped 0..1) colors.
    """
    out = []
    for p in data.get("pairs", []):
        w = worst_parameter(p.get("issue_details"))
        if w is None:
            continue
        rank = rank_severity(w)
        if rank <= 0.0:
            continue
        clamped = min(1.0, rank)            # 0..1, only for the color ramp
        if clamped < min_severity:
            continue
        out.append({
            "res_id1": p["res_id1"], "res_id2": p["res_id2"],
            "sequence": p.get("sequence"),
            "lw_class": (p.get("classification") or {}).get("lw_class"),
            "rank": rank, "sev": clamped,
            "worst_issue": w.get("issue"), "worst_z": w.get("zprime"),
            "tier": w.get("tier", ""),
            "issue": w.get("issue"), "details": p.get("issue_details", []),
        })
    out.sort(key=lambda r: -r["rank"])
    return out


def hbond_count_pairs(data):
    """Pairs whose only/also-present issue is incorrect_hbond_count (summary)."""
    out = []
    for p in data.get("pairs", []):
        hb = [d for d in (p.get("issue_details") or [])
              if d.get("issue") == "incorrect_hbond_count"]
        if hb:
            out.append({"res_id1": p["res_id1"], "res_id2": p["res_id2"],
                        "sequence": p.get("sequence"),
                        "lw_class": (p.get("classification") or {}).get("lw_class")})
    return out


def format_summary(data, flagged, hbcount, obj):
    """Human-readable ranked summary string (ranked by ProSco/|Z'|, no score)."""
    n_oc = sum(1 for f in flagged if f["tier"] == "Of Concern")
    n_al = sum(1 for f in flagged if f["tier"] == "Allowed")
    lines = [f"PARSE: {obj}  {len(data.get('pairs', []))} pairs, "
             f"{len(flagged)} flagged ({n_oc} Of Concern, {n_al} Allowed)"]
    if flagged:
        lines.append("worst pairs (ranked by worst-parameter |Z'|):")
        for f in flagged[:25]:
            lines.append(f"  {pair_label(f)}  "
                         f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
                         f"{severity_headline(f)}  [{failing_params(f['details'])}]")
    # Only surface count issues for pairs that are NOT otherwise flagged/colored
    # (a colored pair already shows its count issue in its [...] breakdown above).
    flagged_keys = {(f["res_id1"], f["res_id2"]) for f in flagged}
    only_count = [h for h in hbcount if (h["res_id1"], h["res_id2"]) not in flagged_keys]
    if only_count:
        lines.append(f"H-bond-count only ({len(only_count)}, summary only, not colored):")
        for h in only_count[:25]:
            lines.append(f"  {pair_label(h)}  "
                         f"{h.get('lw_class') or '?'} {h.get('sequence') or ''}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# PyMOL integration
# ---------------------------------------------------------------------------

_BINARY = os.environ.get("PAIRFINDER_BINARY", "pairfinder")
_OVERLAY = "parse_overlay"
_PROBLEMS = "parse_problems"
_FOCUS = "parse_focus"        # worklist: the current pair (bright sticks)
_CONTEXT = "parse_context"    # worklist: residues within a few A (faded lines)
_HB = "parse_hb"              # worklist: polar contacts within the current pair
_TMP_OBJ = "_parse_tmp"

#: worklist state set by parse_score, walked by parse_next / parse_prev / parse_goto
_QUEUE = []          # flagged pairs, ordered chain-grouped then worst-first
_POS = -1            # index of the current pair in _QUEUE (-1 = not started)
_SCORED_OBJ = None   # object the queue was built against
_CONTEXT_RADIUS = 8.0
#: single-char labels used to rename nucleotide chains before a legacy-PDB save
#: (A-Z, a-z, 0-9 = 62). Legacy PDB has a 1-char chain column and a 99999-atom
#: serial limit, so a big multi-char-chain structure (e.g. a ribosome with chain
#: "1A") would otherwise be saved with chains truncated/collapsed -> the scored
#: res_ids no longer match the in-PyMOL object and nothing gets highlighted.
_CHAIN_POOL = string.ascii_uppercase + string.ascii_lowercase + string.digits


def _resolve_binary():
    if os.path.sep in _BINARY and os.path.exists(_BINARY):
        return _BINARY
    found = shutil.which(_BINARY)
    if found:
        return found
    raise RuntimeError(
        f"pairfinder binary not found ('{_BINARY}'). "
        f"Set it with: parse_set_binary /path/to/pairfinder")


def _run_pairfinder(pdb_path):
    binary = _resolve_binary()
    out = subprocess.run([binary, pdb_path, "--no-download"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"pairfinder failed:\n{out.stderr.strip()}")
    return json.loads(out.stdout)


def _backmap_chain(res_id, inv):
    """Rewrite a res_id's chain via the inverse remap (single-char -> real)."""
    chain, name, seq = parse_res_id(res_id)
    return f"{inv.get(chain, chain)}-{name}-{seq}"


def _score_object(cmd, obj):
    """Save `obj` to a temp PDB, score it, return the parsed JSON.

    The nucleotides are copied to a scratch object whose chains are renamed to
    unique single characters first, so the legacy-PDB round-trip can't truncate
    or collapse multi-character chains (e.g. a ribosome's "1A"/"1B"). The scored
    res_ids are then mapped back to the object's real chain names so the overlay
    selection matches the untouched original object. Protein/solvent is dropped
    (the scorer only scores base pairs); coordinates are never modified.

    Falls back to a plain whole-object save when there are no nucleotide chains
    or too many to fit the single-char pool.
    """
    # Broad nucleotide selection: PyMOL's nucleic polymer plus any residue that
    # carries a C1' sugar atom, so modified / HETATM bases are never dropped.
    nsel = f"(({obj}) and (polymer.nucleic or byres (({obj}) and name C1')))"
    chains = cmd.get_chains(nsel)

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp.close()
    inv = None
    try:
        if chains and len(chains) <= len(_CHAIN_POOL):
            fwd = {c: _CHAIN_POOL[i] for i, c in enumerate(chains)}
            inv = {v: k for k, v in fwd.items()}
            cmd.delete(_TMP_OBJ)
            cmd.create(_TMP_OBJ, nsel)
            try:
                for orig, single in fwd.items():
                    cmd.alter(f'{_TMP_OBJ} and chain "{orig}"', f"chain='{single}'")
                cmd.sort(_TMP_OBJ)
                cmd.save(tmp.name, _TMP_OBJ)
            finally:
                cmd.delete(_TMP_OBJ)
        else:
            # No nucleotides, or more chains than the single-char pool: best-effort
            # direct save (small single-char-chain structures already round-trip).
            cmd.save(tmp.name, obj)
        data = _run_pairfinder(tmp.name)
    finally:
        os.unlink(tmp.name)

    if inv is not None:
        for p in data.get("pairs", []):
            p["res_id1"] = _backmap_chain(p["res_id1"], inv)
            p["res_id2"] = _backmap_chain(p["res_id2"], inv)
    return data


def parse_set_binary(path):
    """Point the plugin at a specific pairfinder binary."""
    global _BINARY
    _BINARY = path
    print(f"PARSE: binary = {path}")


def parse_clear(_self=None):
    """Remove every PARSE overlay/selection; the original object is untouched."""
    from pymol import cmd
    global _QUEUE, _POS, _SCORED_OBJ
    if _SCORED_OBJ is not None:
        cmd.set("cartoon_transparency", 0.0, _SCORED_OBJ)
        _set_hud(cmd, "")
    for o in (_OVERLAY, _PROBLEMS, _FOCUS, _CONTEXT, _HB):
        cmd.delete(o)
    _QUEUE, _POS, _SCORED_OBJ = [], -1, None
    print("PARSE: cleared")


def parse_score(obj=None, min_severity=0.0, zoom=0, gray_background=1, _self=None):
    """Score the current (or named) object and highlight problematic base pairs.

    obj             object to score (default: first enabled object)
    min_severity    only flag issues with severity >= this (0..1; default 0 = all)
    zoom            1 to zoom onto the flagged pairs
    gray_background 1 to color the whole structure gray (so Preferred geometry is
                    gray and only the flagged pairs carry color); 0 to leave the
                    object's existing colors. Coordinates are never modified.
    """
    from pymol import cmd
    min_severity = float(min_severity)

    if obj is None:
        objs = cmd.get_object_list()
        if not objs:
            print("PARSE: no object loaded")
            return
        obj = objs[0]

    try:
        data = _score_object(cmd, obj)
    except Exception as e:        # noqa: BLE001 - surface any failure to the user
        print(f"PARSE: error: {e}")
        return

    flagged = flagged_pairs(data, min_severity)
    hbcount = hbond_count_pairs(data)
    print(format_summary(data, flagged, hbcount, obj))

    # Reset any previous worklist focus, then build a fresh queue for this run.
    global _QUEUE, _POS, _SCORED_OBJ
    cmd.delete(_FOCUS)
    cmd.delete(_CONTEXT)
    cmd.delete(_HB)
    cmd.set("cartoon_transparency", 0.0, obj)
    _QUEUE = build_queue(flagged)
    _POS = -1
    _SCORED_OBJ = obj

    # Gray the whole structure so Preferred geometry recedes and the flagged
    # pairs (colored sticks, built below) stand out. Color only, not coordinates.
    if int(gray_background):
        cmd.color("gray70", obj)

    # Rebuild a fresh overlay of just the flagged pairs.
    cmd.delete(_OVERLAY)
    cmd.delete(_PROBLEMS)
    if not flagged:
        return
    _build_overlay(cmd, obj, flagged)

    cmd.select(_PROBLEMS, _OVERLAY)
    cmd.deselect()
    if int(zoom):
        cmd.zoom(_OVERLAY)

    # Guided HUD: bind nav keys (once) and show a starting status line.
    _bind_nav(cmd)
    n_oc = sum(1 for x in flagged if x["tier"] == "Of Concern")
    _set_hud(cmd, f"PARSE  {len(flagged)} flagged ({n_oc} Of Concern)   "
                  f"press → to start the worklist   (Ctrl-I = inspect a clicked pair)")


def _build_overlay(cmd, obj, flagged):
    """Fast, non-destructive severity overlay.

    Instead of one giant residue selection (slow to parse) and a per-pair color
    loop (hundreds of PyMOL round-trips), we work on a single *copy*: stamp each
    flagged residue's clamped severity into the COPY's B-factor, trim the copy to
    the flagged residues, and color everything in one `spectrum` call. The B-factor
    is written only on the throwaway copy, so the original object's B-factors,
    colors and coordinates are never touched. ~150x faster on a ribosome.
    """
    # {(chain, resi): worst clamped severity} for every flagged residue.
    sevmap = {}
    for f in flagged:
        for rid in (f["res_id1"], f["res_id2"]):
            ch, _, seq = parse_res_id(rid)
            sevmap[(ch, seq)] = max(sevmap.get((ch, seq), 0.0), f["sev"])

    # Copy the nucleotides (simple selection), stamp severity into the copy's B,
    # then drop everything that isn't flagged (sentinel -1 -> removed).
    nsel = f"({obj}) and (polymer.nucleic or byres (({obj}) and name C1'))"
    cmd.create(_OVERLAY, nsel)
    cmd.alter(_OVERLAY, "b = sevmap.get((chain, resi), -1.0)", space={"sevmap": sevmap})
    cmd.remove(f"{_OVERLAY} and b < -0.5")
    cmd.hide("everything", _OVERLAY)
    cmd.show("sticks", _OVERLAY)
    cmd.set("stick_radius", 0.18, _OVERLAY)
    # One pass: clamped severity 0->1 maps yellow->red, matching severity_to_rgb.
    cmd.spectrum("b", "yellow red", _OVERLAY, minimum=0.0, maximum=1.0)


# ---------------------------------------------------------------------------
# Worklist: walk the flagged pairs one at a time (parse_next / prev / goto)
# ---------------------------------------------------------------------------

def build_queue(flagged):
    """Order flagged pairs chain-grouped, worst-first within each chain.

    `flagged` is already globally severity-descending, so the first time a chain
    appears is its worst pair -> chains end up ordered by their worst pair, and
    each chain's pairs stay worst-first. Walking the queue then keeps consecutive
    steps in the same region instead of teleporting across the structure.
    """
    by_chain, order = {}, []
    for f in flagged:
        ch = parse_res_id(f["res_id1"])[0]
        if ch not in by_chain:
            by_chain[ch] = []
            order.append(ch)
        by_chain[ch].append(f)
    queue = []
    for ch in order:
        queue.extend(by_chain[ch])
    return queue


def failing_params(details, max_n=4):
    """Short 'which way is it distorted' string from the worst colored issues.

    e.g. 'shear Z'=+11.3, opening Z'=+10.3, propeller Z'=+8.5'. Of-Concern first,
    then by |Z'|; summary-only issues (hbond count) are excluded.
    """
    ci = [d for d in (details or []) if d.get("issue") not in SUMMARY_ONLY]
    ci.sort(key=lambda d: (d.get("tier") != "Of Concern",
                           -abs(d.get("zprime") or 0.0)))
    parts = []
    for d in ci[:max_n]:
        z = d.get("zprime")
        parts.append(f"{d['issue']} Z'={z:+.1f}" if z is not None else d["issue"])
    return ", ".join(parts)


def phenix_selection(f):
    """Phenix-style atom-selection string for a flagged pair (for hand-off)."""
    c1, _, s1 = parse_res_id(f["res_id1"])
    c2, _, s2 = parse_res_id(f["res_id2"])
    return f"(chain {c1} and resseq {s1}) or (chain {c2} and resseq {s2})"


def _focus(cmd, idx, dim=1):
    """Frame queue[idx]: bright pair + faded local context, H-bonds, label."""
    global _POS
    _POS = idx
    f = _QUEUE[idx]
    obj = _SCORED_OBJ
    s1, s2 = residue_selection(f["res_id1"]), residue_selection(f["res_id2"])
    pair_sel = f"({obj}) and ({s1} or {s2})"

    cmd.delete(_FOCUS)
    cmd.delete(_CONTEXT)
    cmd.delete(_HB)

    # De-emphasise the global view so the current pair reads clearly.
    cmd.disable(_OVERLAY)
    if int(dim):
        cmd.set("cartoon_transparency", 0.85, obj)

    # The pair itself: bright severity-colored sticks.
    cmd.create(_FOCUS, pair_sel)
    cmd.hide("everything", _FOCUS)
    cmd.show("sticks", _FOCUS)
    cmd.set("stick_radius", 0.25, _FOCUS)
    cmd.set_color("parse_focus_col", severity_to_rgb(f["sev"]))
    cmd.color("parse_focus_col", _FOCUS)

    # Local context: nearby residues as thin gray lines (what's crowding it).
    cmd.create(_CONTEXT, f"byres (({obj}) within {_CONTEXT_RADIUS} of ({pair_sel}))"
                         f" and not ({s1} or {s2})")
    cmd.hide("everything", _CONTEXT)
    cmd.show("lines", _CONTEXT)
    cmd.color("gray50", _CONTEXT)

    # Polar contacts within the pair (a stand-in for its H-bonds).
    cmd.distance(_HB, _FOCUS, _FOCUS, mode=2)
    cmd.hide("labels", _HB)

    # One label on the pair: rank, pair, class, and worst-parameter severity.
    label = (f"#{idx + 1} {pair_label(f)} "
             f"{f.get('lw_class') or '?'} {f.get('sequence') or ''} "
             f"— {severity_headline(f)}")
    cmd.label(f"{_FOCUS} and {s1} and name C1'", repr(label))

    cmd.zoom(_FOCUS, 5)

    # Console: worst-parameter headline + full breakdown + Phenix selection.
    n = len(_QUEUE)
    print(f"PARSE [{idx + 1}/{n}] {pair_label(f)}  "
          f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
          f"worst: {severity_headline(f)}  ({f['tier']})")
    print(f"  distortion: {failing_params(f['details']) or '(none colored)'}")
    print(f"  phenix:     {phenix_selection(f)}")

    # On-screen HUD (viewport title): always-visible position + current pair.
    _set_hud(cmd, f"PARSE {idx + 1}/{n}  {pair_label(f)}  "
                  f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
                  f"{severity_headline(f)} [{f['tier']}]   "
                  f"(→ next  ← prev  Ctrl-I = inspect click)")


def _set_hud(cmd, text):
    """Show a one-line status in the viewport (object title), or clear it."""
    if _SCORED_OBJ is None:
        return
    try:
        cmd.set_title(_SCORED_OBJ, 1, text)
    except Exception:        # noqa: BLE001 - title is cosmetic, never fail on it
        pass


def _require_queue():
    if not _QUEUE:
        print("PARSE: no worklist - run parse_score first")
        return False
    return True


def parse_next(_self=None):
    """Step to the next (less severe) flagged pair in the worklist."""
    from pymol import cmd
    if not _require_queue():
        return
    if _POS >= len(_QUEUE) - 1:
        print(f"PARSE: at the end of the worklist ({len(_QUEUE)} pairs)")
        return
    _focus(cmd, _POS + 1)


def parse_prev(_self=None):
    """Step back to the previous (more severe) flagged pair."""
    from pymol import cmd
    if not _require_queue():
        return
    if _POS <= 0:
        print("PARSE: at the start of the worklist")
        return
    _focus(cmd, _POS - 1)


def parse_goto(which, _self=None):
    """Jump to a worklist entry by 1-based rank, or by a res-id substring."""
    from pymol import cmd
    if not _require_queue():
        return
    token = str(which).strip()
    if token.isdigit():
        idx = int(token) - 1
        if not 0 <= idx < len(_QUEUE):
            print(f"PARSE: out of range 1..{len(_QUEUE)}")
            return
        _focus(cmd, idx)
        return
    for i, f in enumerate(_QUEUE):
        if token in f["res_id1"] or token in f["res_id2"]:
            _focus(cmd, i)
            return
    print(f"PARSE: no worklist pair matching '{token}'")


def parse_list(n=25, _self=None):
    """Print the ranked worklist (rank, pair, class, worst |Z'|, distortion)."""
    if not _require_queue():
        return
    n = int(n)
    print(f"PARSE worklist: {len(_QUEUE)} pairs (chain-grouped, worst-first by |Z'|)")
    for i, f in enumerate(_QUEUE[:n]):
        mark = ">" if i == _POS else " "
        print(f" {mark}{i + 1:>4}  {pair_label(f)}  "
              f"{f.get('lw_class') or '?'} {f.get('sequence') or ''}  "
              f"{severity_headline(f)} {f['tier']}  [{failing_params(f['details'], 2)}]")
    if len(_QUEUE) > n:
        print(f"  ... {len(_QUEUE) - n} more (parse_list {len(_QUEUE)} for all)")


def parse_info(selection="pk1", _self=None):
    """Inspect a clicked/selected residue: jump the worklist to its pair.

    Click an atom in the viewport (that makes the `pk1` selection), then run
    `parse_info` (or press `i`). If the residue is a flagged pair, the worklist
    jumps there and prints its breakdown; otherwise it says so.
    """
    from pymol import cmd
    if not _require_queue():
        return
    picks = []
    try:
        cmd.iterate(selection, "picks.append((chain, resi))", space={"picks": picks})
    except Exception:        # noqa: BLE001 - empty/invalid selection
        picks = []
    if not picks:
        print("PARSE: nothing selected — click an atom first, or pass a selection")
        return
    chain, resi = picks[0]
    for i, f in enumerate(_QUEUE):
        for rid in (f["res_id1"], f["res_id2"]):
            c, _name, s = parse_res_id(rid)
            if c == chain and s == resi:
                _focus(cmd, i)
                return
    print(f"PARSE: {chain}/{resi} is not a flagged pair (Preferred or unscored)")


#: nav keybindings — PyMOL's set_key only takes special/modifier keys (not plain
#: letters), so arrows drive the worklist and Ctrl-I inspects the clicked atom.
_NAV_KEYS = (("right", lambda *a: parse_next()),
             ("left", lambda *a: parse_prev()),
             ("CTRL-I", lambda *a: parse_info("pk1")))


def _bind_nav(cmd):
    """Bind → next / ← prev / Ctrl-I inspect (idempotent, never raises)."""
    for key, fn in _NAV_KEYS:
        try:
            cmd.set_key(key, fn)
        except Exception:        # noqa: BLE001 - some keys unmappable on a platform
            pass


def parse_keys(state="on", _self=None):
    """Bind/unbind worklist keys: → next, ← prev, Ctrl-I inspect clicked atom."""
    from pymol import cmd
    if str(state).lower() in ("off", "0", "false"):
        for key, _fn in _NAV_KEYS:
            try:
                cmd.set_key(key, lambda *a: None)
            except Exception:    # noqa: BLE001
                pass
        print("PARSE: keyboard navigation off")
        return
    _bind_nav(cmd)
    print("PARSE keys: → next   ← prev   Ctrl-I = inspect clicked atom "
          "(or click + `parse_info`).  parse_keys off to unbind.")


# Register as PyMOL commands (no-op if imported outside PyMOL).
try:
    from pymol import cmd as _cmd
    _cmd.extend("parse_score", parse_score)
    _cmd.extend("parse_clear", parse_clear)
    _cmd.extend("parse_set_binary", parse_set_binary)
    _cmd.extend("parse_next", parse_next)
    _cmd.extend("parse_prev", parse_prev)
    _cmd.extend("parse_goto", parse_goto)
    _cmd.extend("parse_list", parse_list)
    _cmd.extend("parse_info", parse_info)
    _cmd.extend("parse_keys", parse_keys)
    print("PARSE plugin loaded: parse_score / parse_clear / parse_set_binary / "
          "parse_next / parse_prev / parse_goto / parse_list / parse_info / parse_keys")
except ImportError:
    pass
