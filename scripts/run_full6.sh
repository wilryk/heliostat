#!/usr/bin/env bash
# full6: the whole field with shading and blocking traced as real geometry.
#
# Same dates, rays and worker count as full4/full5 so the three are comparable;
# the only difference is --occluders, which swaps in the model carrying occluder
# slots and puts each heliostat's neighbours in its ray path.
#
# The lock is not politeness. Two sweeps writing one run directory produced
# duplicated summary rows and a raw-ray index that disagreed with the ray file,
# and the corruption was only caught later by a crash in `compare`.

set -u

REPO="C:/gitlab"
OUT="analysis_output/full6"
LOG="$REPO/analysis_output/full6.log"
LOCKDIR="$REPO/analysis_output/.full6.lock"

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another run holds $LOCKDIR -- exiting"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT   # rm -rf, not rmdir: the lock holds a pid file
echo "$$" > "$LOCKDIR/pid"

say "starting full6 (occluders traced), lock held by pid $$"
python -u -m beamdown sweep \
    --occluders \
    --all-heliostats \
    --dates 2026-03-20 2026-06-21 2026-09-22 2026-12-21 \
    --rays 120000 \
    --workers 1 \
    --output "$OUT" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}

if [ "$status" -ne 0 ]; then
    say "full6 exited with status $status"
    exit "$status"
fi

say "full6 complete"
python -u -m beamdown compare analysis_output/full5 "$OUT" \
    --labels scalar traced 2>&1 | tee -a "$LOG"
say "done"
