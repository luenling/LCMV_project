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
RUNBASE=$BASEDIR/Run_0355
RUN_ID=0355

# go through all bam files to create fastqs
# for i in BSF_0355*_S_[1-9]*.bam ; do
#   echo at `date`:  >> $LOGFILE
#   echo bash ${SCRIPTS}/bsf0355_preprocess.sh $i FASTQ >> $LOGFILE
#   bash ${SCRIPTS}/bsf0355_preprocess.sh $i FASTQ
#   ES=$?
#   echo finished at `date` with exit state $ES >> $LOGFILE
#   [ $ES -eq 0 ] || exit $ES
# done
# map to mouse and extract viruses

# for i in FASTQ/BSF*_npa_1.fq.gz ; do
#   echo at `date`: >> $LOGFILE
#   echo bash ${SCRIPTS}/bsf0355_run_mm10.sh $i MAPPING >> $LOGFILE
#   bash ${SCRIPTS}/bsf0355_run_mm10.sh $i MAPPING
#   ES=$?
#   echo finished at `date` with exit state $ES >> $LOGFILE
#   [ $ES -eq 0 ] || exit $ES
# done

# get coverages
echo at `date`: >> $LOGFILE
echo get coverage data >> $LOGFILE
mkdir -P MAPPING/Coverages
cd MAPPING/Coverages
ES=$?
echo finished at `date` with exit state $ES >> $RUNBASE/$LOGFILE
[ $ES -eq 0 ] || exit $ES
ls ../*viral_rh.bam | xargs -L 1 -P 5 -I{} bash -c "genomeCoverageBed -d -ibam '{}'  > `basename '{}'`.coverage"
ES=$?
echo finished at `date` with exit state $ES >> $RUNBASE/$LOGFILE
[ $ES -eq 0 ] || exit $ES

for i in *viral_rh.bam.coverage; do
  A=${i/_sorted*};
  A=${A#*${RUN_ID}_};
  echo $A | cat - <(cut -f 3 $i) | paste  bsf_${RUN_ID}_all.coverages - > temp;
  mv -f temp bsf_${RUN_ID}_all.coverages;
done

# create coverage diagrams
Rscript $SCRIPTS/coverage_diagrams.R bsf_${RUN_ID}_all.coverages
ES=$?
echo finished at `date` with exit state $ES >> $RUNBASE/$LOGFILE
[ $ES -eq 0 ] || exit $ES
cd $RUNBASE



exit 0



ls ../*viral_rh.bam | xargs -L 1 -P 5 -I{} bash -c "genomeCoverageBed -d -ibam '{}'  > '{}'.coverage"

echo -e "CHROM\tBPS" | cat - <(cut -f 1,2 S_10_S19422_sorted_viral_rh.bam.coverage) > bsf_0277_all.coverages

for i in *viral_rh.bam.coverage; do A=${i/_S19*}; echo $A | cat - <(cut -f 3 $i) | paste  bsf_0277_all.coverages - > temp; mv -f temp bsf_0277_all.coverages; done

java -Xmx10g -jar ~/LCMV_project/Tools/GenomeAnalysisTK.jar -T DepthOfCoverage  -nt 10 -R ~/LCMV_project/References/viruses.fasta -o allbams.depth -I bam.list  --countType COUNT_FRAGMENTS --omitIntervalStatistics --omitPerSampleStats

~/LCMV_project/Run_0277/Mappings/Maxcov10K  > for i in ../*_rh.bam; do python ../../../Scripts/downsample_to_cov.py -b $i -c 10000; samtools index `basename $i .bam`_max_cov_10000.bam; done

nohup bash -c ' for i in ../*_max_cov_10000.bam; do bash \
/Volumes/Temp/Lukas/LCMV_project/Scripts/call_lofreq.sh $i; done; \
bash /Volumes/Temp/Lukas/LCMV_project/Scripts/combine_lofreq_vars.sh;\
bash /Volumes/Temp/Lukas/LCMV_project/Scripts/call_lofreq_withbed.sh all_samp_bcf.bed'
