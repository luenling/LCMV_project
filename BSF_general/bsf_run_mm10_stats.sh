#!/bin/bash
#----------
# author: Lukas Endler
# date: 2017-09-24T21:35:14.964+02:00
# Time-stamp: 9/25/2017, 11:22:27 AM Lukas Endler
# takes one fastq file and an outdir name (def. MAPPING), derives second filename runs bwa with 12 threads and outputs a bam file called like the fastq run and sample ID.
# call with bash command blub_1.fq outdir >> logfile.log 2>> log.error.log
# sort the file afterwards
#--------------


source $(dirname $BASH_SOURCE)"/bsf_params.sh"

OUTDIR='MAPPING/'
if [ $2 ]; then
    OUTDIR=$2'/'
    if [ ! -d $2 ]; then mkdir $2 ; fi
fi

STATS=${OUTDIR}/Stats
if [ ! -d $STATS ]; then mkdir $STATS ; fi
LOGFILE=${OUTDIR}/$OUT.log
ERRORFILE=${OUTDIR}/$OUT.err.log

FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# second read - should be kinda save
# R2=`echo $1 | sed 's/_1\.f/_2.f/'`
# # create read group string out of filename
# RG=$(echo $FN | perl -ne 'chomp; $a = $_; if ($a=~/BSF_(\d+)_([^_]+)_(\d+)_S_(\d+)_/) { $SM=sprintf("S%02d",$4);print "\@RG\\tID:$a\\tLB:BSF_$1_$SM\\tSM:$SM\\tPL:illumina\\tPU:$2.$3"} else {print""}')
# [ $RG == "" ] && exit 1 "$FN not of the format BSF_##__HL##_#_S_#_"
# # new filename from now on run and samplename
FN=$(echo $FN | perl -ne 'chomp; $a = $_; if ($a=~/BSF_(\d+)_([^_]+)_(\d+)_S_(\d+)_/) { $SM=sprintf("S%02d",$4);print "BSF_$1_$SM"} else {print"$a"}')
FN=${OUTDIR}/$FN
# LOGFILE=${FN}.log
# ERRORLOG=${FN}.err.log
# echo "start bwa mem  at" `date` >> $LOGFILE
# echo $BWA mem -R $RG -k 17 -r 1.25 -M -t 17 $REFGENOME_FULL $1 $R2  2\>\> $ERRORLOG  \| $SAMTOOLS view -Shb - \| $SAMTOOLS sort -T ${FN}_temp - \> $FN"_sorted.bam"  >> $LOGFILE
#
# $BWA mem -R $RG  -k 17 -r 1.25 -M -t 17 $REFGENOME_FULL $1 $R2 2>> $ERRORLOG | $SAMTOOLS view -Shb - |  $SAMTOOLS sort -T ${FN}_temp - > $FN"_sorted.bam"
# ES=$?
# echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES
# $SAMTOOLS index ${FN}_sorted.bam
BN=`basename $FN`
{
#echo flagstat >> $LOGFILE
$SAMTOOLS flagstat ${FN}_sorted.bam > $STATS/$BN"_full.flagstat"
#echo idxstats >> $LOGFILE
$SAMTOOLS idxstats ${FN}_sorted.bam > $STATS/$BN"_full.idxstats"
# #echo stats >> $LOGFILE
$SAMTOOLS stats ${FN}_sorted.bam > $STATS/$BN"_full.stats"
} &
# # extract only the viral genomes
# $SAMTOOLS view -bh -f 2 -F 256 ${FN}_sorted.bam 'L' 'S' > ${FN}_sorted_viral.bam
# $SAMTOOLS index ${FN}_sorted_viral.bam

# reheader

FN=${FN}_sorted_viral

# LOGFILE=${FN}.log
# ERRORLOG=${FN}.err.log

# echo reheader $FN at `date` >> $LOGFILE
# echo $SAMTOOLS reheader \<\($SAMTOOLS view -H ${FN}.bam \| grep -Ev \''^@SQ.*SN:[^LS]'\' \) ${FN}.bam \> ${FN}_rh.bam >> $LOGFILE
# $SAMTOOLS reheader <($SAMTOOLS view -H ${FN}.bam | grep -Ev '^@SQ.*SN:[^LS]') ${FN}.bam > ${FN}_rh.bam 2>> $ERRORLOG
# ES=$?
# echo finished reheader at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES
#
# #echo rm -f ${FN}.bam >> $LOGFILE
# #rm -f ${FN}.bam  2>> $ERRORLOG
#
# $SAMTOOLS index ${FN}_rh.bam
BN=`basename $FN`
 {

$SAMTOOLS flagstat ${FN}_rh.bam > $STATS/${BN}_rh.flagstat
$SAMTOOLS idxstats ${FN}_rh.bam > $STATS/${BN}_rh.idxstats
$SAMTOOLS stats ${FN}_rh.bam > $STATS/${BN}_rh.stats
#
} &
