#!/usr/bin/env bash
set -euo pipefail

RUNLIST="${1:?Usage: $0 runlist}"

while read -r RUN
do
    [[ -z "${RUN}" ]] && continue

    echo
    echo "##########################################"
    echo "START ${RUN}"
    date
    echo "##########################################"

    ./scripts/multisample/02_preprocess_one.sh "${RUN}"

    echo
    echo "FINISHED ${RUN}"
    date

done < "${RUNLIST}"

echo
echo "ALL MULTISAMPLE PREPROCESSING COMPLETE"
