#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <9/28/2017, 9:35:20 PM lukasendler>
# date: 20.9.2015 at 12:23
# takes a file with list of bam files and vcf from samtools for depths and calls variants with lofreq2
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
SCRIPTS=${BASEDIR}/LCMV_project/BSF_0355/
REFGENOME=$BASEDIR/References/viruses_short.fasta
FN=`basename $1 .list`
FN=${FN/_noprime*/}
FN=${FN/_sorted*/}
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

# do call lofreq
while read file; do
    bash $SCRIPTS/call_lofreq.sh $file
done < $1

bash $SCRIPTS/combine_lofreq_vars.sh

bash $SCRIPTS/call_lofreq_withbed.sh ./all_samp_bcf.bed $2
