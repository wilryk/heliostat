#!/usr/bin/env bash
# Wait for the full4 sweep to finish, then run full5 with the corrected axicon
# shape formula so the two can be compared.
#
# The two runs must differ ONLY in the shape coefficients, so every other
# parameter here is copied from full4's manifest: the same four dates, all 645
# heliostats, 120,000 rays, one worker.
#
# Waiting is necessary rather than polite: the HASP key has one reachable seat,
# and starting the second sweep while the first still holds it raises a modal
# licence dialog that blocks until someone clicks it.

set -u

REPO="C:/gitlab"
BASE="$REPO/analysis_output/full4"
OUT="analysis_output/full5"
LOG="$REPO/analysis_output/full5_launch.log"
TOTAL=44
STALL_LIMIT=$((45 * 60))   # give up rather than launch blind into a dead sweep

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# Only one copy of this may ever reach the launch step. Two watchers both firing
# on full4's completion started two sweeps 33 s apart, both writing to the same
# run directory: duplicated summary rows and a raw-ray index that disagreed with
# the ray file by exactly one heliostat. mkdir is atomic, so it is the lock.
LOCKDIR="$REPO/analysis_output/.full5.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another watcher holds $LOCKDIR -- exiting rather than starting a second sweep"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT   # rm -rf, not rmdir: the lock holds a pid file
echo "$$" > "$LOCKDIR/pid"

say "watching $BASE for $TOTAL timesteps  (lock held by pid $$)"
last_count=-1
last_change=$(date +%s)

while true; do
    count=$(ls "$BASE/flux" 2>/dev/null | grep -c '\.npy$' || echo 0)
    now=$(date +%s)

    if [ "$count" != "$last_count" ]; then
        say "full4 at $count/$TOTAL"
        last_count=$count
        last_change=$now
    fi

    if [ "$count" -ge "$TOTAL" ]; then
        say "full4 complete"
        break
    fi

    if [ $((now - last_change)) -gt "$STALL_LIMIT" ]; then
        say "ABORT: no new timestep in $((STALL_LIMIT / 60)) min, stuck at $count/$TOTAL."
        say "full4 may have died. Not launching full5 -- check and start it by hand."
        exit 1
    fi

    sleep 60
done

# The seat is not released the instant the last file lands; the process still
# has to tear its Quadoa session down. README: do not reopen in quick succession.
say "waiting 180 s for the licence seat to be released"
sleep 180

if pgrep -f "beamdown sweep.*$OUT" >/dev/null 2>&1; then
    say "ABORT: a sweep is already writing $OUT. Refusing to start a second one."
    exit 1
fi

say "starting full5 (corrected axicon shape)"
python -m beamdown sweep \
    --all-heliostats \
    --dates 2026-03-20 2026-06-21 2026-09-22 2026-12-21 \
    --rays 120000 \
    --workers 1 \
    --output "$OUT" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}

if [ "$status" -ne 0 ]; then
    say "full5 exited with status $status"
    exit "$status"
fi

say "full5 complete -- comparison:"
python -m beamdown compare analysis_output/full4 "$OUT" \
    --labels original corrected --attribute 2>&1 | tee -a "$LOG"
say "done"
