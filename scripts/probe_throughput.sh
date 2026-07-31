#!/usr/bin/env bash
# How fast can we actually trace, and does a second Quadoa session help?
#
# The 12-date grid is ~104,000 traces. At full7's measured 918 ms/heliostat that
# is 26 hours, so the two levers -- worker count and ray count -- are worth
# measuring rather than guessing. This machine has 4 physical cores and a single
# traceRays call already runs at ~4x parallelism, so a second session may only
# contend; the only honest way to know is to time it.
#
# One date (13 timesteps) over the downselect, occluders on, so the numbers
# carry over to the real run.

set -u
cd "C:/gitlab/heliostats" || exit 1
LOG="analysis_output/probe_throughput.log"
: > "$LOG"

run() {
    tag="$1"; w="$2"; r="$3"
    echo "=== $tag : workers=$w rays=$r ===" | tee -a "$LOG"
    rm -rf "analysis_output/probe_$tag"
    t0=$(date +%s)
    python -u -m beamdown sweep --occluders --dates 2026-03-20 \
        --rays "$r" --workers "$w" --no-resume \
        --output "analysis_output/probe_$tag" >> "$LOG" 2>&1
    status=$?
    t1=$(date +%s)
    echo "WALL $tag = $((t1 - t0)) s   (exit $status)" | tee -a "$LOG"
}

run w1r120 1 120000
run w2r120 2 120000
run w1r60  1 60000
run w2r60  2 60000
echo "ALL DONE" | tee -a "$LOG"
