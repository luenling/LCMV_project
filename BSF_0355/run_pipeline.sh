#!/bin/bash
#----------
# author: Lukas Endler
# created: 9/25/2017, 2:34:27 PM
# Time-stamp: 9/25/2017, 2:34:37 PM lukasendler>
# takes a bam file, reheaders it to only the two viral segments in short format
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
SCRIPTS=$BASEDIR/LCMV_project/BSF_0355/
LOGFILE=pipeline.log

# go through all bam files to create fastqs
for i in BSF_0355*_S_[1-9]*.bam ; do
  echo at `date`:  >> $LOGFILE
  echo bash ${SCRIPTS}/bsf0355_preprocess.sh $i FASTQ >> $LOGFILE
  bash ${SCRIPTS}/bsf0355_preprocess.sh $i FASTQ
  ES=$?
  echo finished at `date` with exit state $ES >> $LOGFILE
  [ $ES -eq 0 ] || exit $ES
done
# map to mouse and extract viruses

for i in FASTQ/BSF*_npa_1.fq.gz ; do
  echo at `date`:
  echo bash ${SCRIPTS}/bsf0355_run_mm10.sh $i MAPPING >> $LOGFILE
  bash ${SCRIPTS}/bsf0355_run_mm10.sh.sh $i MAPPING
  ES=$?
  echo finished at `date` with exit state $ES >> $LOGFILE
  [ $ES -eq 0 ] || exit $ES
done
