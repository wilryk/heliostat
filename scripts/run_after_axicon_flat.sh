#!/usr/bin/env bash
# Chain watcher: when axicon_flat finishes CLEANLY, run the two prime-focus
# sweeps back to back. Safe to run any time; safe to run twice (each sweep
# script takes its own lock and exits politely if it is held).
#
# It only chains on a clean finish: lock gone AND the log's final "done"
# line present. A stale lock (e.g. the launching shell was killed) makes it
# wait forever rather than start a second seat-holder -- in that case check
# for a live python first, remove the stale .axicon_flat.lock by hand, and
# relaunch run_axicon_flat.sh (the sweep resumes from what it stored).
#
# If a chained sweep fails, the chain STOPS. Never retry a seat failure --
# each attempt pops a modal H0038 dialog and can leak the seat (run_full8.sh
# documents the measured disaster).

set -u

REPO="C:/gitlab/heliostats"
cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] chain: $*"; }

say "waiting for axicon_flat to release the seat"
while [ -d "analysis_output/.axicon_flat.lock" ]; do sleep 120; done

if ! tail -5 "analysis_output/axicon_flat.log" 2>/dev/null | grep -q "] done"; then
    say "axicon_flat lock is gone but its log has no 'done' -- it did not"
    say "finish cleanly. NOT chaining; resume axicon_flat first."
    exit 1
fi

say "axicon_flat done -- starting prime_focus (focused)"
bash scripts/run_prime_focus.sh || { say "prime_focus failed; chain stopped"; exit 1; }

if ! tail -5 "analysis_output/prime_focus.log" 2>/dev/null | grep -q "] done"; then
    say "prime_focus did not report done; chain stopped"
    exit 1
fi

say "prime_focus done -- starting prime_focus_flat"
bash scripts/run_prime_focus_flat.sh || { say "prime_focus_flat failed; chain stopped"; exit 1; }
say "chain complete: axicon_flat, prime_focus, prime_focus_flat all done"
