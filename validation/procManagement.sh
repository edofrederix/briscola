#!/bin/bash

# Define an array to keep track of the number of procs used per job

declare -A PROCS_PER_JOB

# Total number of procs used

PROCS_COUNT=0

# Function to update the resource usages. Free up procs in the total proc count
# when a pid is no longer alive.

check_jobs() {
    for PID in "${!PROCS_PER_JOB[@]}"; do
        if ! kill -0 "$PID" 2>/dev/null; then
            ((PROCS_COUNT -= PROCS_PER_JOB[$PID]))
            unset PROCS_PER_JOB[$PID]
        fi
    done
}

# Function to wait for resources to come available. First arg: procs needed,
# second arg: total available procs.

wait_for_procs() {
    while ((PROCS_COUNT + NPROC > $2)); do
        echo "Procs running = $PROCS_COUNT, Procs needed = $1"
        check_jobs
        sleep 1
    done
}

# Function to store the resource usage for a task. First arg: number of procs
# used, second arg: pid

store_procs() {
    PROCS_PER_JOB[$2]=$1
    ((PROCS_COUNT += NPROC))
}
