# Shared session memory

Cross-session memory for every Claude session (and human) on this project.
The goal is that a cold session becomes productive by reading **two small
files**, not by re-reading the codebase.

## The three layers

| Layer | Brain analogue | Lives in | Lifetime |
|---|---|---|---|
| Short-term | Episodic — "what happened recently" | `memory/journal/` (one file per working session) | Last **2–3** sessions kept verbatim; older ones distilled then deleted |
| Mid-term | Working memory — "what is true right now" | `memory/STATE.md` (a single overwritten snapshot) | Always current; rewritten at session end |
| Long-term | Semantic — "what we know and why" | see the table below | Permanent, curated |

### Long-term homes in THIS repo

This repo predates the memory system and already has a strong convention:
**decisions live next to the thing they decide, with the measurement that
justified them.** Keep doing that; the memory layer points at it.

| What | Home | Notes |
|---|---|---|
| Traps, invariants, session rules | `/CLAUDE.md` | auto-loaded; keep it short |
| Decisions + the evidence behind them | `/README.md` | the de-facto devlog: measured behaviours, layout designs, vetting results, API traps live here as narrative sections |
| Run-parameter decisions | `config.toml` comments | e.g. the monthly-DNI rationale, the 1-worker licence cap — the comment IS the record |
| Per-run / per-experiment rationale | `scripts/*.sh` and `scripts/*.py` headers | e.g. run_full8.sh documents why that run differs from full7 |
| Regression truth | `tests/verify.py` stages | a decision enforced by a test needs no prose |
| Work queue | session task list + STATE.md "open" section | no tasks/ directory; loose ends must land in STATE.md or die |

## Session-start protocol (cold session, human or agent coordinator)

1. `CLAUDE.md` is auto-loaded — costs nothing.
2. Read `memory/STATE.md`. That is the resume point.
3. Read the newest file in `memory/journal/`. Now you know what happened
   last time and why the state is what it is.
4. **Do not** re-read README or source files wholesale. Follow pointers
   from STATE.md only as the task requires.

Executor subagents keep getting self-contained task specs as before —
they generally don't need this directory, and their prompts say what to read.

## Session-end protocol (whoever coordinated the session)

1. **Rewrite `STATE.md`** — overwrite, never append. It answers: what's
   running, is the suite green, what's queued, what decisions are open,
   what changed lately. Hard cap ~80 lines.
2. **Write one journal entry** `journal/YYYY-MM-DD-<slug>.md`, ≤40 lines:
   what was attempted, what landed (with commit ids), what failed and why,
   owner decisions made, loose ends. Compression is the point — a journal
   is a summary you'd give a colleague, not a transcript.
3. **Promote before you delete** ("consolidation"): when journals rotate
   out (keep the newest 2–3), anything still relevant must already live in
   its long-term home per the table above. If it has no home, it wasn't
   worth remembering.
4. Commit. Memory that isn't committed dies with the machine (a power
   event killed a 24 h sweep here on 2026-07-30; the code survived only
   because it was on disk). No remote exists yet — pushing is an open item.

## Rules that keep this cheap

- STATE.md ≤ ~80 lines, journals ≤ ~40. If it doesn't fit, it belongs in a
  long-term home with a pointer.
- No content lives ONLY here except journals-in-window and the snapshot.
  This directory is an index and a diary, not a second wiki — README,
  config comments and script headers stay the single source of truth.
- Never paste code, test bodies, or file contents into memory files; link
  paths instead.
- Never put run DATA here — `analysis_output/` is git-ignored (34 GB) and
  its numbers belong in README sections once vetted.
- The coordinator's private user-preference memory (how the owner likes to
  work) lives outside the repo and is not duplicated here.
