#!/bin/bash
## Basic script to sort a vcf, remove the chr prefix, compress and tabix index
# Params: input output is_pre_zipped
# Author: Job van Schipstal
if [ "$3" == "True" ]; then
  bgzip -dc "$1" | sed "s/^chr//" | awk '$1 ~ /^#/ {print $0;next}
  {print $0 | "sort -k1,1V -k2,2n"}' | bgzip -c > "$2"
else
  cat "$1" | sed "s/^chr//" | awk '$1 ~ /^#/ {print $0;next}
  {print $0 | "sort -k1,1V -k2,2n"}' | bgzip -c > "$2"
fi
tabix -p vcf -f "$2"
