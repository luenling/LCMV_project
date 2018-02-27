#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of vcf files (bgzipped and indexed) with one sample each and runs SNPgenie for each
#--------------

source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"
LOGFILE=`pwd`/"bsf_"${RUN_ID}"_snpgenie.log"

SNPGENIE="perl $BASEDIR/Tools/snpgenie/snpgenie.pl --minfreq 0.0001 --vcfformat 2 "
declare -A REFS_rc
declare -A GTFS_rc
declare -A REFS
declare -A GTFS

REFS[L]=$BASEDIR/References/L.fasta
GTFS[L]=$BASEDIR/References/L.gtf
REFS_rc[L]=$BASEDIR/References/L_revcom.fasta
GTFS_rc[L]=$BASEDIR/References/L_revcom.gtf
REFS[S]=$BASEDIR/References/S.fasta
GTFS[S]=$BASEDIR/References/S.gtf
REFS_rc[S]=$BASEDIR/References/S_revcom.fasta
GTFS_rc[S]=$BASEDIR/References/S_revcom.gtf

while read fn; do
  SAMP=`bcftools view $fn | grep "#CHR" | cut -f 10`
  # create folder
  mkdir ${SAMP}_snpgenie
  cd ${SAMP}_snpgenie
  for CHR in L S; do
    # go through each chromosome
    bcftools view -i 'TYPE="snp"' $fn $CHR | bcftools norm -m +snps - > ${SAMP}_${CHR}.vcf
    python ${BASEDIR}/LCMV_project/rev_comp_vcf.py --in ${SAMP}_${CHR}.vcf > ${SAMP}_${CHR}_rc.vcf
    # forward
    $SNPGENIE --fastafile ${REFS[$CHR]} --gtffile ${GTFS[$CHR]} --minfreq 0.001 --snpreport ${SAMP}_${CHR}.vcf --vcfformat 2
    mv SNPGenie_Results ${SAMP}_${CHR}_fw
    # reverse
    $SNPGENIE --fastafile ${REFS_rc[$CHR]} --gtffile ${GTFS_rc[$CHR]} --minfreq 0.001 --snpreport ${SAMP}_${CHR}_rc.vcf --vcfformat 2
    mv SNPGenie_Results ${SAMP}_${CHR}_rc
    # combine results

    # awk 'NR == 1; NR > 1 {print $0 | "sort -n"}'

  done
  cd ..
done < $1


exit 0

 perl ../../../Tools/CHASeq/vcf2revcom.pl ../../../References/L.fasta ../../../References/L.gtf  S48_samp_L.lofreq.bed.vcf
"
