#!/bin/bash
#----------
# author: Lukas Endler
# created: 9/25/2017, 2:34:27 PM
# Time-stamp: 9/25/2017, 2:34:37 PM lukasendler>
# takes a bam file, reheaders it to only the two viral segments in short format
#--------------

source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"

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

# # get coverages
# echo at `date`: >> $LOGFILE
# echo get coverage data >> $LOGFILE
# mkdir -p MAPPING/Coverages
# cd MAPPING/Coverages
# ES=$?
# echo finished at `date` with exit state $ES >> $RUNBASE/$LOGFILE
# [ $ES -eq 0 ] || exit $ES
# ls ../*viral_rh.bam | xargs -L 1 -P 5 -I{} bash -c "genomeCoverageBed -d -ibam '{}'  > `basename '{}'`.coverage"
# ES=$?
# echo finished at `date` with exit state $ES >> $RUNBASE/$LOGFILE
# [ $ES -eq 0 ] || exit $ES
#
# for i in *viral_rh.bam.coverage; do
#   A=${i/_sorted*};
#   A=${A#*${RUN_ID}_};
#   echo $A | cat - <(cut -f 3 $i) | paste  bsf_${RUN_ID}_all.coverages - > temp;
#   mv -f temp bsf_${RUN_ID}_all.coverages;
# done
#
# # create coverage diagrams
# Rscript $SCRIPTS/coverage_diagrams.R bsf_${RUN_ID}_all.coverages
# ES=$?
# echo finished at `date` with exit state $ES >> $RUNBASE/$LOGFILE
# [ $ES -eq 0 ] || exit $ES
# cd $RUNBASE

#Downsample to max cov 50K,not really necessary I guess
# mkdir -p $RUNBASE/MAPPINGS/MaxCov50K
# for i in ../*S3[1789]*_rh.bam; do
#   python ${SCRIPTS}/downsample_to_cov.py -b $i -c 50000;
#   samtools index `basename $i .bam`_max_cov_50K.bam;
# done

#ls $RUNBASE/MAPPING/*_rh.bam > bam_list.txt
#echo $SAMTOOLS mpileup -f $REFGENOME -q 30 -Q 30 -I -v -d 100000000 -u -k DP,AD,ADF,ADR,SP -o  all_samps_samtools.vcf -b  bam_list.txt >> $LOGFILE
#$SAMTOOLS mpileup -f $REFGENOME -q 30 -Q 30 -I -v -d 100000000 -u -k DP,AD,ADF,ADR,SP -o  all_samps_samtools.vcf -b  bam_list.txt
#bgzip all_samps_samtools.vcf
#tabix -p vcf all_samps_samtools.vcf.gz

# Do lofreq 2 calling
# mkdir -p $RUNBASE/LOFREQ2
# cd LOFREQ2
# bash $SCRIPTS/do_full_lofreq.sh ../bam_list.txt ../all_samps_samtools.vcf.gz


# echo mkdir -p $RUNBASE/BQSR >> $LOGFILE
# mkdir -p $RUNBASE/BQSR
# echo cd $RUNBASE/BQSR >> $LOGFILE
# cd $RUNBASE/BQSR
# echo  ls $RUNBASE/LOFREQ2/*.bam \> bam.list >> $LOGFILE
# ls $RUNBASE/LOFREQ2/*.bam > bam.list
# echo vcf2bed \< $RUNBASE/LOFREQ2/lofreq2_all_samp_bed_norm_0.05.vcf \| cut -f 1-5 - \> lofreq2_all_samp_bed_norm_0.05_5col.bed $LOGFILE
# vcf2bed < $RUNBASE/LOFREQ2/lofreq2_all_samp_bed_norm_0.05.vcf | cut -f 1-5 - > lofreq2_all_samp_bed_norm_0.05_5col.bed
#
# echo java -Xmx20G -jar $GATK -T BaseRecalibrator  -R $REFGENOME -knownSites lofreq2_all_samp_bed_norm_0.05_5col.bed -o recal_afs_0.005.tab -I bam.list >> $LOGFILE
# java -Xmx20G -jar $GATK -T BaseRecalibrator  -R $REFGENOME -knownSites lofreq2_all_samp_bed_norm_0.05_5col.bed -o recal_afs_0.005.tab -I bam.list 2>> bqsr.err.log
# echo java -Xmx20G -jar $GATK -T BaseRecalibrator  -R $REFGENOME -knownSites lofreq2_all_samp_bed_norm_0.05_5col.bed -BQSR recal_afs_0.005.tab -I bam.list -o recal_afs_0.005_secondpass.table >> $LOGFILE
# java -Xmx20G -jar $GATK -T BaseRecalibrator  -R $REFGENOME -knownSites lofreq2_all_samp_bed_norm_0.05_5col.bed -BQSR recal_afs_0.005.tab -I bam.list -o recal_afs_0.005_secondpass.table  2>> bqsr.err.log
# echo java -jar $GATK -T AnalyzeCovariates -R $REFGENOME  -before recal_afs_0.005.tab -after recal_afs_0.005_secondpass.table -plots BQSR.pdf >> $LOGFILE
# java -jar $GATK -T AnalyzeCovariates -R $REFGENOME  -before recal_afs_0.005.tab -after recal_afs_0.005_secondpass.table -plots BQSR.pdf  2>> bqsr.err.log
# while read p || [[ -n $p ]]; do
#   FN=`basename $p .bam`
#   echo java -jar $GATK -T PrintReads  -R $REFGENOME -dt NONE -BQSR recal_afs_0.005.tab -I $p -o ${FN}_bqsr.bam >> $LOGFILE
#   java -jar $GATK -T PrintReads  -R $REFGENOME -dt NONE -BQSR recal_afs_0.005.tab -I $p -o ${FN}_bqsr.bam 2>> bqsr.err.log
# done < bam.list

cd $RUNBASE/BQSR
ls "*.bam" > bam_list.txt
bash $SCRIPTS/do_full_lofreq.sh ../bam_list.txt ../all_samps_samtools.vcf.gz

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

 ~/LCMV_project/Tools/VarDictJava/VarDict_java -G ../../References/viruses_short.fasta -f 0.01 -b ../BQSR/BSF_0355_S10_real_viterbi_IDQS_bqsr.bam -c 1 -S 2 -E 3 -z 1 viruses.bed | teststrandbias.R | var2vcf_valid.pl -N S10 -E -f 0.01 > test_S10.vcf

java -jar $GATK -T AnalyzeCovariates -R $REFGENOME  -BQSR recal_afs_0.005.tab -plots BQSR.pdf
