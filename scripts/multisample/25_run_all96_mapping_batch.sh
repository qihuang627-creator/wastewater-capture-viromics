#!/usr/bin/env bash
set -euo pipefail

SAMPLES=(
ERR14788990
ERR14788988
ERR14788991
ERR14788992
ERR14788995
ERR14789044
ERR14788864
ERR14788954
)

for RUN in "${SAMPLES[@]}"
do
    echo
    echo "########################################"
    echo "ALL96 competitive mapping: ${RUN}"
    echo "########################################"

    ./scripts/multisample/24_all96_mapping_one.sh "${RUN}"
done
