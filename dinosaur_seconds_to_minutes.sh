#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Usage: ${0##*/} input.tsv output.tsv"
    exit 1
fi

input=$1
output=$2

 grep -v "^#" ${input} | awk -v OFS='\t' '{if (NR > 1) $2=$2/60; print}' > ${output}
