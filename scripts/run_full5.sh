#!/usr/bin/env bash
# Run full5 (corrected axicon shape) and compare it against full4.
#
# Separate from run_full5_after_full4.sh because full4 has already finished --
# there is nothing left to wait for. Same lock, for the same reason: two of
# these writing one run directory produced duplicated summary rows and a
# raw-ray index that disagreed with the ray file by one heliostat's worth.

set -u

REPO="C:/gitlab/heliostats"
OUT="analysis_output/full5"
LOG="$REPO/analysis_output/full5_launch.log"
LOCKDIR="$REPO/analysis_output/.full5.lock"

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another run holds $LOCKDIR -- exiting"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT   # rm -rf, not rmdir: the lock holds a pid file
echo "$$" > "$LOCKDIR/pid"

if pgrep -f "beamdown sweep.*$OUT" >/dev/null 2>&1; then
    say "ABORT: a sweep is already writing $OUT"
    exit 1
fi

say "starting full5 (corrected axicon shape), lock held by pid $$"
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
