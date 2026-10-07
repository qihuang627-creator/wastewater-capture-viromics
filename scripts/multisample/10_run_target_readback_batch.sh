#!/usr/bin/env bash
set -euo pipefail

RUNLIST="${1:-config/multisample/run_lists/all_runs.txt}"

while read -r RUN
do
    [[ -z "${RUN}" ]] && continue

    echo
    echo "##########################################"
    echo "START TARGET READ-BACK ${RUN}"
    date
    echo "##########################################"

    ./scripts/multisample/09_target_readback_one.sh "${RUN}"

    echo
    echo "FINISHED TARGET READ-BACK ${RUN}"
    date

done < "${RUNLIST}"

echo
echo "ALL MULTISAMPLE TARGET READ-BACK COMPLETE"
