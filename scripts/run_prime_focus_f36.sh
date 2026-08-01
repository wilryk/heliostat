#!/usr/bin/env bash
# prime_focus_f36: prime focus with F1 LOWERED to (0,0,36000) mm -- 9 m above
# the axicon tip -- so the aim geometry matches what the axicon field actually
# does: its aim rays cross the axis between tip+2.6 m and tip+10.0 m (far
# heliostat; field mean +7.8 m; scripts/aoi_stats.py). Matching this keeps
# blocking/shadowing comparable to the axicon runs and holds mirror diameter
# and receiver position; the mirrors are slightly faster F/# than at 47 m.
#
# 4-day sample (both solstices + two mid declinations), not the full 7: this is
# the comparison baseline for the fixed-mirror-figure scenarios, all of which
# trace the SAME 4 days:
#   run_prime_focus_f36.sh          adaptive figure (re-solved every instant)
#   run_prime_focus_f36_sphere.sh   fixed spherical (correct RoC, no astig)
#   run_prime_focus_f36_meancos.sh  fixed astig, year mean weighted by cos(AOI)
#   run_prime_focus_f36_median.sh   fixed astig, year median
# Expect ~6 h each at ~615 ms/heliostat. Everything else follows
# run_prime_focus.sh (same model: pf_height is a parameter, the write moves
# the detector to 36000 with it -- that is the whole point of the parameter).
#
# ONE worker, deliberately -- see run_full8.sh. The lock is not politeness.

set -u

REPO="C:/gitlab/heliostats"
NAME="prime_focus_f36"
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

say "starting $NAME (4-day sample, F1 36000, adaptive figure), pid $$"
python -u -m beamdown sweep \
    --secondary prime_focus \
    --focus-height-mm 36000 \
    --n-mirrors 1 \
    --model-file models/heliostat_field_prime_focus.optx \
    --all-heliostats \
    --dates 2026-12-21 2026-02-20 2026-04-21 2026-06-21 \
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
