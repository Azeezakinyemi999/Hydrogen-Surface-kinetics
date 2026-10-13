#!/usr/bin/env bash
# Run the analyses for one run folder.
#
#   bash run.sh runs/ni_first_try                       # everything listed in its settings.py
#   bash run.sh runs/ni_first_try --only curves,sweep    # a subset
#   bash run.sh runs/ni_first_try --workers 4 --recompute
#   bash run.sh runs/ni_first_try check                  # only check the settings (fast)
#   bash run.sh runs/ni_first_try status                 # what is cached / saved
#
# The output goes to the screen and to <run folder>/run.log (appended).
# Environment variables: PYTHON (default: the rmg4.0 conda env's python).
set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,10p' "$0"; exit 1
fi
RUN_DIR="${1%/}"; shift
SETTINGS="$RUN_DIR/settings.py"
[[ -f "$SETTINGS" ]] || { echo "no settings.py in $RUN_DIR (an RMG input.py is not the package settings file)" >&2; exit 1; }

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

LOG="$RUN_DIR/run.log"
{
    echo "=== $(date '+%Y-%m-%d %H:%M:%S')  $COMMAND $SETTINGS $*"
    echo "    python: $PYTHON"
} | tee -a "$LOG"

set +e
${KEEP_AWAKE[@]+"${KEEP_AWAKE[@]}"} "$PYTHON" -m hydrogen_surface_kinetics "$COMMAND" "$SETTINGS" "$@" 2>&1 | tee -a "$LOG"
STATUS=${PIPESTATUS[0]}
set -e
echo "=== $(date '+%Y-%m-%d %H:%M:%S')  exit $STATUS" | tee -a "$LOG"
exit "$STATUS"
