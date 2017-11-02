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
  SM=$($SAMTOOLS view -H BSF_0355_S01_real_viterbi_IDQS_bqsr.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g")
  $SAMTOOLS view -q 1 -h -F 256 -f 2 -G 96  $i  | $SAMTOOLS view -b -G 144 -  | $SAMTOOLS depth -m 100000000 -a - > Coverages/${SM}_plus.depth
  $SAMTOOLS view -q 1 -h -F 256 -f 2 -G 80  $i  | $SAMTOOLS view -b -G 160 -  | $SAMTOOLS depth -m 100000000 -a - > Coverages/${SM}_minus.depth

done < $1
