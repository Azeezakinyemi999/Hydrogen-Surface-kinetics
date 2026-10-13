#!/usr/bin/env bash
# Run the analyses for one or more run folders, one after another.
#
#   bash run.sh runs/ni_first_try                         # everything listed in its settings.py
#   bash run.sh runs/try_2/Cu runs/try_2/Fe runs/try_2/Ni  # several runs, in this order
#   bash run.sh runs/try_*/Ni/settings.py                 # a folder or its settings.py; globs work
#   bash run.sh runs/ni_first_try --only curves,sweep     # a subset
#   bash run.sh runs/ni_first_try --workers 4 --recompute
#   bash run.sh runs/try_2/* check                        # only check the settings (fast)
#   bash run.sh runs/ni_first_try status                  # what is cached / saved
#
# The runs come first, then check or status (default: run), then the options, which apply to
# every run. Every run is looked up before the first starts. A failed run does not stop the
# rest; the exit status is that of the last failed run (0 if none). Ctrl-C stops them all.
# Each run's output goes to the screen and to <run folder>/run.log (appended).
# Environment variables: PYTHON (default: the rmg4.0 conda env's python).
set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,16p' "$0"; exit 1
fi

SETTINGS=()
while [[ $# -gt 0 && "$1" != -* && "$1" != check && "$1" != status ]]; do
    ARG="${1%/}"; shift
    if [[ -d "$ARG" ]]; then
        [[ -f "$ARG/settings.py" ]] || { echo "no settings.py in $ARG (an RMG input.py is not the package settings file)" >&2; exit 1; }
        ARG="$ARG/settings.py"
    fi
    [[ -f "$ARG" ]] || { echo "no run folder or settings file at $ARG" >&2; exit 1; }
    SETTINGS+=("$ARG")
done
[[ ${#SETTINGS[@]} -gt 0 ]] || { echo "give the run folders first: bash run.sh runs/<name> $1" >&2; exit 1; }

COMMAND=run
if [[ $# -gt 0 && ( "$1" == check || "$1" == status ) ]]; then
    COMMAND="$1"; shift
fi

PYTHON="${PYTHON:-$HOME/anaconda3/envs/rmg4.0/bin/python}"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python)"

# One thread per worker process: the parallelism comes from the worker pool.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export MPLBACKEND=Agg PYTHONUNBUFFERED=1

# Keep the Mac awake for long runs. -i alone is not enough: maintenance sleep pauses the
# workers and the job time limits then fire falsely. -dims blocks display, idle, disk and
# system sleep (on AC power) until the run ends.
KEEP_AWAKE=()
if [[ "$COMMAND" == run ]] && command -v caffeinate >/dev/null; then
    KEEP_AWAKE=(caffeinate -dims)
fi

# Without this, a Ctrl-C that Python catches would only end the current run and start the next.
LOG=/dev/null
trap 'echo "=== $(date "+%Y-%m-%d %H:%M:%S")  interrupted" | tee -a "$LOG"; exit 130' INT

FAILED=()
LAST_BAD=0
for SET in "${SETTINGS[@]}"; do
    LOG="$(dirname "$SET")/run.log"
    {
        echo "=== $(date '+%Y-%m-%d %H:%M:%S')  $COMMAND $SET $*"
        echo "    python: $PYTHON"
    } | tee -a "$LOG"

    set +e
    ${KEEP_AWAKE[@]+"${KEEP_AWAKE[@]}"} "$PYTHON" -m hydrogen_surface_kinetics "$COMMAND" "$SET" "$@" 2>&1 | tee -a "$LOG"
    STATUS=${PIPESTATUS[0]}
    set -e
    echo "=== $(date '+%Y-%m-%d %H:%M:%S')  exit $STATUS" | tee -a "$LOG"
    if [[ $STATUS -ne 0 ]]; then
        FAILED+=("$SET (exit $STATUS)")
        LAST_BAD=$STATUS
    fi
done

if [[ ${#SETTINGS[@]} -gt 1 ]]; then
    echo "=== $(( ${#SETTINGS[@]} - ${#FAILED[@]} )) of ${#SETTINGS[@]} runs OK"
    for F in ${FAILED[@]+"${FAILED[@]}"}; do
        echo "    failed: $F"
    done
fi
exit "$LAST_BAD"
