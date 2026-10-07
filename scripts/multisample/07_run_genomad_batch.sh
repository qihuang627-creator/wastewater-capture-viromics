#!/usr/bin/env bash
set -euo pipefail

RUNLIST="${1:-config/multisample/run_lists/all_runs.txt}"

while read -r RUN
do
    [[ -z "${RUN}" ]] && continue

    echo
    echo "##########################################"
    echo "START geNomad ${RUN}"
    date
    echo "##########################################"

    ./scripts/multisample/06_genomad_one.sh "${RUN}"

    echo
    echo "FINISHED geNomad ${RUN}"
    date

done < "${RUNLIST}"

echo
echo "ALL MULTISAMPLE geNomad RUNS COMPLETE"
