#!/usr/bin/env bash
# Weekend chain with the owner's Sunday rule built in: NO analysis runs on
# Sunday, midnight to midnight, local time (religious observance -- a hard
# constraint, not a preference).
#
# Queue, in order (each ~6-9 h, each takes the ONE licence seat):
#   run_cassegrain.sh          settled design, focused      (~9 h)
#   run_cassegrain_flat.sh     same, flat mirrors           (~9 h)
#   run_prime_focus_f36.sh     F1 36 m, adaptive figure     (~6 h)
#   run_prime_focus_f36_sphere.sh / _meancos.sh / _median.sh (fixed figures)
#
# Sunday handling:
#   - Before each launch: if it is Sunday (or Saturday past 23:50), wait
#     until Monday 00:05.
#   - While a run is live: at Saturday 23:55 it is STOPPED (process tree
#     killed, lock removed). Sweeps RESUME from stored steps -- this is the
#     documented stale-lock recovery path, losing at most one in-flight
#     step -- and the same script is relaunched after Sunday ends.
#   - A run that fails for any OTHER reason stops the chain. In particular
#     a seat failure is never retried (H0038 modal; see run_full8.sh).
#
# Safe to run twice: each run script takes its own lock and exits politely
# if it is held; the chain treats that polite exit as "seat busy" and stops
# rather than queueing behind an unknown seat holder.

set -u

REPO="C:/gitlab/heliostats"
cd "$REPO" || exit 1
LOG="$REPO/analysis_output/weekend_chain.log"
say() { echo "[$(date '+%a %H:%M:%S')] chain: $*" | tee -a "$LOG"; }

RUNS=(cassegrain cassegrain_flat prime_focus_f36 prime_focus_f36_sphere
      prime_focus_f36_meancos prime_focus_f36_median)

sunday_now() {
    # Sunday all day, or Saturday from 23:50 (no point starting anything).
    local dow hm
    dow=$(date +%u); hm=$(( 10#$(date +%H%M) ))
    [ "$dow" = 7 ] || { [ "$dow" = 6 ] && [ "$hm" -ge 2350 ]; }
}

wait_out_sunday() {
    if sunday_now; then
        say "Sunday (or too close to it) -- waiting until Monday 00:05"
        while sunday_now || { [ "$(date +%u)" = 1 ] && [ $(( 10#$(date +%H%M) )) -lt 5 ]; }; do
            sleep 300
        done
        say "Monday -- resuming"
    fi
}

run_one() {
    local name="$1" script="scripts/run_$1.sh" lockdir="analysis_output/.$1.lock"
    while true; do
        wait_out_sunday
        say "launching $script"
        bash "$script" >> "$LOG" 2>&1 &
        local pid=$!
        local stopped=0
        while kill -0 "$pid" 2>/dev/null; do
            if sunday_now; then
                say "Saturday 23:55 -- stopping $name for Sunday (resumes Monday)"
                taskkill //PID "$pid" //T //F >> "$LOG" 2>&1
                sleep 10
                rm -rf "$lockdir"
                stopped=1
                break
            fi
            sleep 60
        done
        wait "$pid" 2>/dev/null
        if [ "$stopped" = 1 ]; then
            continue                       # wait out Sunday, relaunch, resume
        fi
        if tail -5 "analysis_output/$name.log" 2>/dev/null | grep -q "] done"; then
            say "$name finished cleanly"
            return 0
        fi
        if tail -20 "analysis_output/$name.log" 2>/dev/null | grep -q "another run holds"; then
            say "$name found the seat busy -- NOT queueing blindly; chain stopped"
        else
            say "$name did not finish cleanly -- chain stopped (never retry a seat failure)"
        fi
        return 1
    done
}

say "weekend chain started; queue: ${RUNS[*]}"
for name in "${RUNS[@]}"; do
    if tail -5 "analysis_output/$name.log" 2>/dev/null | grep -q "] done"; then
        say "$name already done -- skipping"
        continue
    fi
    run_one "$name" || exit 1
done
say "chain complete: all queued runs done"
