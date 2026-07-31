#!/usr/bin/env bash
# full8: the annual-energy sweep. 12 monthly dates on the corrected time grid,
# every occlusion traced as real geometry (neighbours + axicon).
#
# Three things make this different from full7:
#
#   1. 12 dates instead of 4, so the efficiency surface spans 7 distinct
#      declinations instead of 3 and far less of the year is extrapolated
#      (full7 left 14.4% of daylight hours outside the traced hull).
#   2. The time grid samples uniformly across the true daylight window rather
#      than snapping inward to whole hours. full7's grid started 1.5 h after
#      sunrise on some dates and silently dropped ~20% of each day from the
#      energy integral.
#   3. It ASKED for 60,000 rays rather than 120,000. It did not get them, and
#      the speed claim this point used to make was never measured. Both are
#      corrected in full at the FOOT OF THIS FILE -- see "CORRECTION". In
#      short: the --rays override never reached the workers, so this run traced
#      config.toml's 120,000; the manifest still records 60,000; and nothing in
#      this repository measures what a ray costs. full7 also traced 120,000.
##
# ONE worker, deliberately. Measured: asking for 2 got "asked for 2, 1 got a
# license seat", and the failed second seat then leaked the first one, so the
# NEXT single-worker run found 0 seats and sat waiting 10 minutes. Each attempt
# also pops a modal H0038 dialog a human has to dismiss. Do not raise this.
#
# Note the grid shares no timesteps with full7 (06:42 vs 08:00), so `compare`
# cannot be run between the two -- it matches on timestep keys. The check that
# matters here is the energy cross-check, which report_energy.py prints.
#
# The lock is not politeness. Two sweeps writing one run directory produced
# duplicated summary rows and a raw-ray index that disagreed with the ray file.

set -u

REPO="C:/gitlab/heliostats"
OUT="analysis_output/full8"
LOG="$REPO/analysis_output/full8.log"
LOCKDIR="$REPO/analysis_output/.full8.lock"

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another run holds $LOCKDIR -- exiting"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT   # rm -rf, not rmdir: the lock holds a pid file
echo "$$" > "$LOCKDIR/pid"

say "starting full8 (12 dates, corrected grid, occluders + axicon), pid $$"
python -u -m beamdown sweep \
    --occluders \
    --all-heliostats \
    --rays 60000 \
    --workers 1 \
    --output "$OUT" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}

if [ "$status" -ne 0 ]; then
    say "full8 exited with status $status"
    exit "$status"
fi

say "full8 complete -- annual energy"
python -u scripts/report_energy.py --run "$OUT" 2>&1 | tee -a "$LOG"
say "done"

# ---------------------------------------------------------------------------
# CORRECTION -- appended after this run was launched, deliberately at the end of
# the file because bash was still executing it and shifting any earlier byte
# would have made it resume mid-command.
#
# WHAT WENT WRONG. This run was launched with "--rays 60000". At that time the
# override never reached the sweep workers: a worker process calls load_config
# itself, from disk (sweep._init_worker), so a value applied only to the
# driver's copy of the config reached the report and not the trace. full8
# therefore traced config.toml's rays_per_heliostat = 120,000 per heliostat,
# exactly like full7. Every stored run in analysis_output/ is in the same
# position: whatever any of them asked for, they all traced 120,000.
#
# THE BUG IS FIXED. config.apply_overrides is now replayed inside every worker
# and the whole override set is recorded in the run manifest, so "--rays 60000"
# launched today genuinely traces 60,000. beamdown/config.py, beamdown/cli.py
# and beamdown/sweep.py carry the details.
#
# READING full8's NUMBERS. Its manifest.json records rays_per_heliostat = 60000,
# because the manifest was written from the driver's config; the trace emitted
# 120,000. Watts per ray is source_power_w / rays_emitted, so anything that
# scales stored counts by the MANIFEST's ray budget reads 2x HIGH until that
# field is patched to 120000. Ray counts and every ratio built from them
# (transmission, eta_shade, eta_block) are unaffected -- only absolute power is.
# The Monte-Carlo argument is unaffected too, and cuts the same way at either
# budget: 0.0067/sqrt(645) = 0.03% at 120,000 rays, 0.04% at 60,000, both far
# under the 0.46% trapezoid residual.
#
# WHAT IS STILL UNKNOWN: the cost. "Halving the rays halves a 26 h run" was an
# assumption, not a measurement, and could not have been a measurement, because
# no run has ever traced anything other than 120,000 rays. There is no evidence
# here in either direction -- rays may dominate a trace's time or be a minor
# part of it, and the fixed per-heliostat and per-traceRays-call costs are
# equally unmeasured. Do not repeat the old claim, and do not replace it with
# the opposite one.
#
# HOW TO FIND OUT: scripts/probe_ray_cost.py. It times the fixed per-heliostat
# cost, the fixed per-traceRays-call cost and the marginal per-ray cost
# separately on one session, and prints the ray budget below which halving the
# rays saves less than 10%. It must NOT be run while a sweep holds the licence
# seat, and it refuses to start while any analysis_output/.*.lock exists.
# ---------------------------------------------------------------------------
