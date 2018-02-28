#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of vcf files (bgzipped and indexed) with one sample each and runs SNPgenie for each
#--------------

source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"
LOGFILE=`pwd`/"bsf_"${RUN_ID}"_snpgenie.log"


set -x
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
    # add sample and chrom fields and sort for site, codon and sliding window site_results
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/site_results.txt |  awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_site_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/codon_results.txt |  awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_codon_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/sliding_window_length9_results.txt |  awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_sliding_window_length9_results.txt
    # reverse
    $SNPGENIE --fastafile ${REFS_rc[$CHR]} --gtffile ${GTFS_rc[$CHR]} --minfreq 0.001 --snpreport ${SAMP}_${CHR}_rc.vcf --vcfformat 2
    mv SNPGenie_Results ${SAMP}_${CHR}_rc
    # reverse site, codon results and sliding windows and add sample and chromosome info
    python ${BASEDIR}/LCMV_project/rev_comp_vcf.py -s $CHR --in ${SAMP}_${CHR}_rc/site_results.txt | awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 && ! /noncoding/ ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' |  awk ' NR > 1 {print $0 | "sort -k3,3n"}' >> ${SAMP}_site_results.txt
    python ${BASEDIR}/LCMV_project/rev_comp_vcf.py -s $CHR --in ${SAMP}_${CHR}_rc/codon_results.txt | awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 && ! /noncoding/ ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' |  awk ' NR > 1 {print $0 | "sort -k3,3n"}' >> ${SAMP}_codon_results.txt
    python ${BASEDIR}/LCMV_project/rev_comp_vcf.py -s $CHR --in ${SAMP}_${CHR}_rc/sliding_window_length9_results.txt | awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 && ! /noncoding/ ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' |  awk ' NR > 1 {print $0 | "sort -k3,3n"}' >> ${SAMP}_sliding_window_length9_results.txt
    # combine results
    cat ${SAMP}_site_results.txt >> ../all_samps_site_results.txt
    cat ${SAMP}_codon_results.txt >> ../all_samps_codon_results.txt
    cat ${SAMP}_sliding_window_length9_results.txt >> ../all_samps_sliding_window_length9_results.txt
    # awk 'NR == 1; NR > 1 {print $0 | "sort -n"}'
  done
  cd ..
done < $1
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_site_results.txt > all_samps_site_results_sorted.txt
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_codon_results.txt > all_samps_codon_results_sorted.txt
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_sliding_window_length9_results.txt > all_samps_sliding_window_length9_results_sorted.txt

exit 0

 perl ../../../Tools/CHASeq/vcf2revcom.pl ../../../References/L.fasta ../../../References/L.gtf  S48_samp_L.lofreq.bed.vcf
"
