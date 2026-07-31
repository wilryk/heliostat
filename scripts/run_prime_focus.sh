#!/usr/bin/env bash
# prime_focus: focused heliostats aimed at F1 = (0,0,47000) mm -- the receiver
# 20 m ABOVE the axicon vertex, mirroring the beam-down receiver 20 m below.
#
#   - Model: models/heliostat_field_prime_focus.optx (built by guarded text
#     surgery, VERIFIED 10/10 by scripts/verify_prime_focus_model.py on
#     2026-07-31: pf_height reaches the trace, sequence 3 keeps literal
#     ray-count semantics, centroids land on-axis).
#   - --n-mirrors 1: one reflection in the path (no secondary bounce), so
#     throughput = reflectivity**1. This is the layout's whole advantage.
#   - NO --occluders: scalar occlusion, UNION form (see run_axicon_flat.sh).
#   - 7 explicit dates = the 7 distinct declinations (NOT --suggest-dates,
#     which ADDS to config.toml's 12 -- see run_axicon_flat.sh's header).
#   - 120,000 rays, one traceRays call per heliostat.
#
# ONE worker, deliberately -- see run_full8.sh for the measured seat-leak
# disaster behind this. The lock is not politeness.

set -u

REPO="C:/gitlab/heliostats"
NAME="prime_focus"
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

say "starting $NAME (7 distinct declinations, focused mirrors, scalar occlusion), pid $$"
python -u -m beamdown sweep \
    --secondary prime_focus \
    --focus-height-mm 47000 \
    --n-mirrors 1 \
    --model-file models/heliostat_field_prime_focus.optx \
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
