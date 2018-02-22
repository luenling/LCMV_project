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


while read fn; do
  SAMP = `bcftools view $fn | grep "#CHR" | cut -f 10`
  # create folder
  mkdir ${SAMP}_snpgenie
  cd ${SAMP}_snpgenie
  for CHR in L S; do
    # go through each chromosome
    bcftools view -i 'TYPE="snp"' ../Run_0355/BQSR/S07_samp.lofreq.bed.vcf.gz L | bcftools norm -m +snps - > ${SAMP}_${CHR}.vcf
    # forward
    $SNPGENIE --fastafile /Volumes/Temp/Lukas/LCMV/References/L.fasta --gtffile /Volumes/Temp/Lukas/LCMV/References/L.gtf

    # reverse

    # combine results

  done

done < $1

 perl ../../../Tools/CHASeq/vcf2revcom.pl ../../../References/L.fasta ../../../References/L.gtf  S48_samp_L.lofreq.bed.vcf
"
