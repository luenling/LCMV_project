#!/bin/bash

shopt -s extglob

for i in *full.idxstats; do
  j=$(basename $i _full.idxstats)
  SMP=${j##BSF_*([0-9])_}
  BSF=${j%%_S*([0-9])}
  if [ ! -e full_readcounts_${BSF}.tab ]; then
    echo -e "genome\nmouse\nlcmv\nunmapped" > full_readcounts_${BSF}.tab
  fi
  paste <(cat full_readcounts_${BSF}.tab)  <(awk -v SMP="$SMP" 'BEGIN{mouse=0; lcmv=0; unmapped=0} /^[^LS\*]/ {mouse+=$3} /^[LS]/ {lcmv+=$3} /^\*/ {unmapped+=$4} END{print SMP; print mouse; print lcmv; print unmapped}' $i)  > int_file
  mv -f int_file full_readcounts_${BSF}.tab
done
