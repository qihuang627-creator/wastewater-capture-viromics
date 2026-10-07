#!/usr/bin/env bash
set -euo pipefail

RUNLIST="${1:?Usage: $0 runlist}"

while read -r RUN
do
    [[ -z "${RUN}" ]] && continue

    echo
    echo "##########################################"
    echo "START ASSEMBLY ${RUN}"
    date
    echo "##########################################"

    ./scripts/multisample/04_assembly_one.sh "${RUN}"

    echo "FINISHED ASSEMBLY ${RUN}"
    date

done < "${RUNLIST}"

echo
echo "ALL MULTISAMPLE ASSEMBLIES COMPLETE"
