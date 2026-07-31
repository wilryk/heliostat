#!/usr/bin/env bash
# axicon_flat: first of the 5 comparison sweeps for the geometry paper.
#
#   - FLAT heliostats (--flat-mirrors: Zernike c3/c4/c5 forced to zero,
#     pointing identical). Axicon secondary = config.toml default.
#   - NO --occluders: neighbour shading/blocking applied as read-time
#     scalars in the UNION form (owner decision 2026-07-31; the product
#     form double-charged overlapping losses, 0.338% low annually; union
#     is 0.114%). Scalars stay invalid for low-sun instantaneous power
#     and per-mirror claims -- annual comparisons only.
#   - 7 declination-spaced dates, not 12 monthly ones: 5 of the monthly
#     dates are declination duplicates (full8 measured the pairs equal to
#     <=0.04%), so the efficiency surface is identical either way.
#   - 120,000 rays (single-heliostat SNR decision), one traceRays call
#     per heliostat (--rays-per-trace equal to --rays; chunking measured
#     as pure loss at this scale).
#
# ONE worker, deliberately. Measured: asking for 2 got "asked for 2, 1 got a
# license seat", and the failed second seat then leaked the first one, so the
# NEXT single-worker run found 0 seats and sat waiting 10 minutes. Each attempt
# also pops a modal H0038 dialog a human has to dismiss. Do not raise this.
#
# The lock is not politeness. Two sweeps writing one run directory produced
# duplicated summary rows and a raw-ray index that disagreed with the ray file.

set -u

REPO="C:/gitlab/heliostats"
NAME="axicon_flat"
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

say "starting $NAME (7 declination dates, flat mirrors, scalar occlusion), pid $$"
python -u -m beamdown sweep \
    --flat-mirrors \
    --all-heliostats \
    --suggest-dates 7 \
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
