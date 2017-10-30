#!/bin/bash
#----------
# author: Lukas Endler
# created: 9/25/2017, 2:34:27 PM
# Time-stamp: 9/25/2017, 2:34:37 PM lukasendler>
# takes a bam file, reheaders it to only the two viral segments in short format
#--------------
source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"

FN=`basename $1 .bam`

LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

echo reheader $FN at `date` >> $LOGFILE
echo $SAMTOOLS reheader \<\($SAMTOOLS view -H $1 \| grep -Ev \''^@SQ.*SN:[^g]'\' \| sed \''s/gi|86440167|gb|DQ361066.1|/L/; s/gi|116563461|gb|DQ361065.2|/S/'\'\) $1 \> ${FN}_rh.bam >> $LOGFILE
$SAMTOOLS reheader <($SAMTOOLS view -H $1 | grep -Ev '^@SQ.*SN:[^g]' | sed 's/gi|86440167|gb|DQ361066.1|/L/; s/gi|116563461|gb|DQ361065.2|/S/') $1 > ${FN}_rh.bam 2>> $ERRORLOG
ES=$?
echo finished reheader at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

echo rm -f $1 >> $LOGFILE
rm -f $1  2>> $ERRORLOG

$SAMTOOLS index ${FN}_rh.bam
echo flagstat ${FN}_rh.bam >> $LOGFILE
$SAMTOOLS flagstat ${FN}_rh.bam >> $LOGFILE
echo idxstats ${FN}_rh.bam >> $LOGFILE
$SAMTOOLS idxstats ${FN}_rh.bam >> $LOGFILE
sam-stats ${FN}_rh.bam > ${FN}_rh.stats


#echo creating coverages at `date` >> $LOGFILE
#if [ ! -d Coverages ]; then mkdir Coverages ; fi
