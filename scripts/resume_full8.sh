#!/usr/bin/env bash
# Resume full8 after the power event killed it at timestep 93/161 (exit 127 --
# a native-level death, consistent with the USB HASP key dropping out from
# under the live session; python and the CLI were verified fine afterwards).
#
# Deliberately NOT scripts/run_full8.sh: that file was patched byte-neutrally
# while bash was executing it, and while its error path demonstrably still
# worked, there is no reason to re-enter a file with that history.
#
# Three flags differ from the original launch, all deliberate:
#
#   --rays 120000        What the first 93 timesteps ACTUALLY traced. The
#                        original launch asked for 60,000 but the override
#                        never reached the workers, so config.toml's 120,000
#                        was used. Asking for it explicitly keeps the run
#                        homogeneous AND makes the manifest true: write_manifest
#                        runs at every sweep start (sweep.py:298), so this
#                        resume rewrites rays_per_heliostat = 120000 and clears
#                        the 2x error that count-scaled readers would otherwise
#                        carry.
#
#   --rays-per-trace 120000   One traceRays call instead of two. MEASURED by
#                        scripts/probe_ray_cost.py on this machine: 120,000 rays
#                        cost 646 ms in 1 call, 698 ms in 2, 776 ms in 4. The
#                        split is pure loss here -- same rays, same source, same
#                        statistics, ~7% less wall clock. Chunking exists for
#                        memory and 120,000 rays is nowhere near needing it.
#
#   --workers 1          Unchanged, and not negotiable: one reachable HASP seat,
#                        and a failed second request leaks the first.
#
# Resume is on by default, so the 93 stored timesteps are skipped and only the
# remaining 68 are traced. At the measured 698 ms/heliostat that is ~8.5 h;
# expect a little more, since the per-timestep binning and writing sit outside
# the traced figure.

set -u

REPO="C:/gitlab"
OUT="analysis_output/full8"
LOG="$REPO/analysis_output/full8.log"
LOCKDIR="$REPO/analysis_output/.full8.lock"

cd "$REPO" || exit 1
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    say "another run holds $LOCKDIR -- exiting"
    exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT
echo "$$" > "$LOCKDIR/pid"

say "resuming full8 from 93/161 (120,000 rays in 1 call), pid $$"
python -u -m beamdown sweep \
    --occluders \
    --all-heliostats \
    --rays 120000 \
    --rays-per-trace 120000 \
    --workers 1 \
    --output "$OUT" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}

if [ "$status" -ne 0 ]; then
    say "full8 resume exited with status $status"
    exit "$status"
fi

say "full8 complete -- annual energy"
python -u scripts/report_energy.py --run "$OUT" 2>&1 | tee -a "$LOG"
say "done"
