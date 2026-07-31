#!/usr/bin/env bash
# prime_focus_flat: run_prime_focus.sh with --flat-mirrors -- identical
# pointing, Zernike c3/c4/c5 forced to zero. See that script's header for
# every other decision; see run_axicon_flat.sh for the union-form scalars.

set -u

REPO="C:/gitlab/heliostats"
NAME="prime_focus_flat"
OUT="analysis_output/$NAME"
LOG="$REPO/analysis_output/$NAME.log"
LOCKDIR="$REPO/analysis_output/.$NAME.lock"

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another run holds $LOCKDIR -- exiting"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT   # rm -rf, not rmdir: the lock holds a pid file
echo "$$" > "$LOCKDIR/pid"

say "starting $NAME (7 distinct declinations, FLAT mirrors, scalar occlusion), pid $$"
python -u -m beamdown sweep \
    --secondary prime_focus \
    --focus-height-mm 47000 \
    --n-mirrors 1 \
    --model-file models/heliostat_field_prime_focus.optx \
    --flat-mirrors \
    --all-heliostats \
    --dates 2026-12-21 2026-01-21 2026-02-20 2026-03-20 2026-04-21 2026-05-21 2026-06-21 \
    --rays 120000 \
    --rays-per-trace 120000 \
    --workers 1 \
    --output "$OUT" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}

if [ "$status" -ne 0 ]; then
    say "$NAME exited with status $status"
    exit "$status"
fi

say "$NAME complete -- annual energy"
python -u scripts/report_energy.py --run "$OUT" 2>&1 | tee -a "$LOG"
say "done"
