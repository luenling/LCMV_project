#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# gets strand specific coverages
#--------------
#
source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"
mkdir -p Coverages

while read i ; do
  SM=$($SAMTOOLS view -H $i | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g")
  $SAMTOOLS view -q 1 -h -F 256 -f 2 -G 96  $i  | $SAMTOOLS view -b -G 144 -  | $SAMTOOLS depth -m 100000000 -a - > Coverages/${SM}_plus.depth
  $SAMTOOLS view -q 1 -h -F 256 -f 2 -G 80  $i  | $SAMTOOLS view -b -G 160 -  | $SAMTOOLS depth -m 100000000 -a - > Coverages/${SM}_minus.depth
done < $1

for i in *_plus.depth; do
  A=$(basename $i _plus.depth);
  echo $A | cat - <(cut -f 3 $i) | paste  bsf_${RUN_ID}_all_plus.depths - > temp;
  mv -f bsf_${RUN_ID}_all_plus.depths;
done

for i in *_minus.depth; do
  A=$(basename $i _minus.depth);
  echo $A | cat - <(cut -f 3 $i) | paste  bsf_${RUN_ID}_all_minus.depths - > temp;
  mv -f bsf_${RUN_ID}_all_minus.depths;
done

exit

Rscript $(dirname $BASH_SOURCE)"/coverage_diagrams.R" bsf_${RUN_ID}_all_plus.depths bsf_${RUN_ID}_all_minus.depths
