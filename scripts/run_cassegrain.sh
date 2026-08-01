#!/usr/bin/env bash
# cassegrain: focused heliostats aimed at F1 = (0,0,38986) mm -- the hyperboloid's
# far (virtual) focus. The hyperboloid relays the bundle to the receiver at its
# near focus, z = 7000 mm (the beam-down receiver, same as the axicon's).
#
#   - Design (scripts/design_cassegrain.py, rim height chosen = the axicon's own
#     rim): rim r 15,000 mm at z 32,460; F1 z 38,986 (aperture-fill closure);
#     vertex z 29,589; conic K -5.8789; |R_vertex| 32,181.3 mm.
#   - Model: models/heliostat_field_cassegrain.optx -- built BY HAND in Quadoa
#     (the hyperboloid's conic constants are literals in the file; the Python
#     side never sees them, see beamdown/secondary/cassegrain.py). The model
#     MUST end sequence 3 on the receiver and MUST accept the two single_params
#     every session writes: sec_height (27000) and rec_offset (-20000), which
#     must leave the receiver at z 7000. Writing a param a model lacks is
#     SILENTLY ignored -- wrong wiring shows up only as plausible-wrong spot
#     sizes. Before the FIRST sweep, verify the model with an adaptation of
#     scripts/verify_prime_focus_model.py (needs a licence seat; the centroid
#     prediction (0,0) still holds -- a hyperboloid is stigmatic between foci).
#   - --n-mirrors 2: heliostat + hyperboloid, so throughput = reflectivity**2,
#     like the axicon (applied at READ time; the manifest records the truth).
#   - NO --occluders: scalar occlusion, UNION form (see run_axicon_flat.sh).
#   - 7 explicit dates = the 7 distinct declinations (NOT --suggest-dates,
#     which ADDS to config.toml's 12 -- see run_axicon_flat.sh's header).
#   - 120,000 rays, one traceRays call per heliostat.
#
# ONE worker, deliberately -- see run_full8.sh for the measured seat-leak
# disaster behind this. The lock is not politeness.

set -u

REPO="C:/gitlab/heliostats"
NAME="cassegrain"
OUT="analysis_output/$NAME"
LOG="$REPO/analysis_output/$NAME.log"
LOCKDIR="$REPO/analysis_output/.$NAME.lock"
MODEL="models/heliostat_field_cassegrain.optx"

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

if [ ! -f "$MODEL" ]; then
    say "$MODEL does not exist -- it is built BY HAND in Quadoa (design"
    say "numbers in this script's header / scripts/design_cassegrain.py)."
    say "Build it, verify it, then rerun. NOT starting."
    exit 1
fi

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another run holds $LOCKDIR -- exiting"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT   # rm -rf, not rmdir: the lock holds a pid file
echo "$$" > "$LOCKDIR/pid"

say "starting $NAME (7 distinct declinations, focused mirrors, scalar occlusion), pid $$"
python -u -m beamdown sweep \
    --secondary cassegrain \
    --focus-height-mm 38986 \
    --rim-height-mm 32460 \
    --n-mirrors 2 \
    --model-file "$MODEL" \
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
