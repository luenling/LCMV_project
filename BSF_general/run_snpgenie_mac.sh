#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of vcf files (bgzipped and indexed with full path) with one sample each and runs SNPgenie for each
#--------------
shopt -s extglob



#source $(dirname $BASH_SOURCE)"/bsf_params.sh"
LOGFILE=`pwd`/$( basename $1 )"_snpgenie.log"
BN=${1//.*}
SLIDEWIN=15
SNPGENIEDIR="${HOME}/Tools/snpgenie"
LCMVPROJ="${HOME}/Data/LCMV_project/"
REFDIR="${HOME}/Data/P_2023_LCEV_LCMV-Evolution_LE/Resources/References/"
#set -x
SNPGENIE="perl ${SNPGENIEDIR}/snpgenie.pl --minfreq 0.01 --vcfformat 2 --slidingwindow $SLIDEWIN "
declare -A REFS_rc
declare -A GTFS_rc
declare -A REFS
declare -A GTFS

REFS[L]=$REFDIR/L.fasta
GTFS[L]=$REFDIR/L.gtf
REFS_rc[L]=$REFDIR/L_revcom.fasta
GTFS_rc[L]=$REFDIR/L_revcom.gtf
REFS[S]=$REFDIR/S.fasta
GTFS[S]=$REFDIR/S.gtf
REFS_rc[S]=$REFDIR/S_revcom.fasta
GTFS_rc[S]=$REFDIR/S_revcom.gtf

while read fn; do
  #get the first sample field as sample name
  SAMP=`bcftools view $fn | grep "#CHR" | cut -f 10`
  # create folder
  mkdir ${SAMP}_snpgenie
  cd ${SAMP}_snpgenie
  for CHR in L S; do
    # go through each chromosome
    bcftools view -i 'TYPE="snp"' $fn $CHR | bcftools norm -m +snps - > ${SAMP}_${CHR}.vcf
    python ${LCMVPROJ}/rev_comp_vcf.py --in ${SAMP}_${CHR}.vcf > ${SAMP}_${CHR}_rc.vcf
    # forward
    $SNPGENIE --fastafile ${REFS[$CHR]} --gtffile ${GTFS[$CHR]} --minfreq 0.01 --snpreport ${SAMP}_${CHR}.vcf --vcfformat 2
    mv SNPGenie_Results ${SAMP}_${CHR}_fw
    # add sample and chrom fields and sort for site, codon and sliding window site_results
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; 
    (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/site_results.txt |  \
    awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_site_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; 
    (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/codon_results.txt |  \
    awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_codon_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; 
    (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/sliding_window_length${SLIDEWIN}_results.txt |  \
    awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_sliding_window_length${SLIDEWIN}_results.txt
    awk -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {$1="Sample"; print $0 }; (NR > 1 ) {$1=samp; print $0 } ' ${SAMP}_${CHR}_fw/product_results.txt  >> ../all_samps_product_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {print "Sample","CHR",$0 }; 
    (NR > 1 ) {print samp, chrom, $0 } ' ${SAMP}_${CHR}_fw/population_summary.txt  >> ../all_samps_population_summary.txt
   # reverse
    $SNPGENIE --fastafile ${REFS_rc[$CHR]} --gtffile ${GTFS_rc[$CHR]} --minfreq 0.001 --snpreport ${SAMP}_${CHR}_rc.vcf --vcfformat 2
    mv SNPGenie_Results ${SAMP}_${CHR}_rc
    # reverse site, only coding, codon results and sliding windows and add sample and chromosome info
    python ${LCMVPROJ}/rev_comp_vcf.py -s $CHR --in ${SAMP}_${CHR}_rc/site_results.txt | \
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; 
    print "Sample","CHR",$0 }; (NR > 1 && ! /noncoding/ ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' |  \
    awk ' NR > 1 {print $0 | "sort -k3,3n"}' >> ${SAMP}_site_results.txt
    python ${LCMVPROJ}/rev_comp_vcf.py -s $CHR --in ${SAMP}_${CHR}_rc/codon_results.txt | \
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site;
    print "Sample","CHR",$0 }; (NR > 1 && ! /noncoding/ ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' |  \
    awk ' NR > 1 {print $0 | "sort -k3,3n"}' >> ${SAMP}_codon_results.txt
    python ${LCMVPROJ}/rev_comp_vcf.py -s $CHR --in ${SAMP}_${CHR}_rc/sliding_window_length${SLIDEWIN}_results.txt | \
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 };
    (NR > 1 && ! /noncoding/ ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' | \
    awk ' NR > 1 {print $0 | "sort -k3,3n"}' >> ${SAMP}_sliding_window_length${SLIDEWIN}_results.txt
    awk -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {$1="Sample"; print $0 }; 
    (NR > 1 ) {$1=samp; print $0 } ' ${SAMP}_${CHR}_rc/product_results.txt  >> ../all_samps_product_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {print "Sample","CHR",$0 }; 
    (NR > 1 ) {print samp, chrom, $0 } ' ${SAMP}_${CHR}_rc/population_summary.txt  >> ../all_samps_population_summary.txt
    # combine results
    cat ${SAMP}_site_results.txt >> ../all_samps_site_results.txt
    cat ${SAMP}_codon_results.txt >> ../all_samps_codon_results.txt
    cat ${SAMP}_sliding_window_length${SLIDEWIN}_results.txt >> ../all_samps_sliding_window_length${SLIDEWIN}_results.txt
    # awk 'NR == 1; NR > 1 {print $0 | "sort -n"}'
    rm -f *.vcf *.txt
  done
  cd ..
  zip -r ${SAMP}_snpgenie.zip ${SAMP}_snpgenie
  rm -rf ${SAMP}_snpgenie
done < $1
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_site_results.txt > ${BN}_site_results_sorted.txt
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_codon_results.txt > ${BN}_codon_results_sorted.txt
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_sliding_window_length${SLIDEWIN}_results.txt > ${BN}_sliding_window_length${SLIDEWIN}_results_sorted.txt
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort"}'  all_samps_product_results.txt > ${BN}_product_results_sorted.txt
awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort"}'  all_samps_population_summary.txt > ${BN}_population_summary_sorted.txt

rm -f all_samps_*_results.txt all_samps_population_summary.txt
# collect all product results
#find . -name "product_results.txt" -exec cat {} +  | sort -u | sed 's/file/Sample/g; s/_[LS].*.vcf//g' > all_samps_product_results.txt

exit 0

 perl ../../../Tools/CHASeq/vcf2revcom.pl ../../../References/L.fasta ../../../References/L.gtf  S48_samp_L.lofreq.bed.vcf
"
