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
#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of vcf files (bgzipped and indexed) with one sample each and runs SNPgenie for each
#--------------

source $(dirname $BASH_SOURCE)"/bsf_params.sh"
LOGFILE=`pwd`/"bsf_"${RUN_ID}"_snpgenie.log"


set -x
SNPGENIE="perl $BASEDIR/Tools/snpgenie/snpgenie.pl --minfreq 0.0001 --vcfformat 2 --slidingwindow 9 "

REFS=$BASEDIR/References/L.fasta
GTFS=$BASEDIR/References/L.gtf


while read fn; do
  SAMP=`bcftools view $fn | grep "#CHR" | cut -f 10`
  # create folder
  mkdir ${SAMP}_snpgenie
  cd ${SAMP}_snpgenie
    bcftools view -i 'TYPE="snp"' | bcftools norm -m +snps - > ${SAMP}.vcf
    $SNPGENIE --fastafile ${REFS[$CHR]} --gtffile ${GTFS[$CHR]} --minfreq 0.001 --snpreport ${SAMP}_${CHR}.vcf --vcfformat 2
    mv SNPGenie_Results ${SAMP}_${CHR}_fw
    # add sample and chrom fields and sort for site, codon and sliding window site_results
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/site_results.txt |  awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_site_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/codon_results.txt |  awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_codon_results.txt
    awk -v chrom="$CHR" -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {site=$3; $3=$1; $1=site; print "Sample","CHR",$0 }; (NR > 1 ) {site=$3; $3=$1; $1=site; print samp,chrom,$0 } ' ${SAMP}_${CHR}_fw/sliding_window_length9_results.txt |  awk ' NR == 1; NR > 1 {print $0 | "sort -k3,3n"}' > ${SAMP}_sliding_window_length9_results.txt
    awk -v samp="$SAMP" -v OFS="\t" -F"\t" '(NR == 1) {$1="Sample"; print $0 }; (NR > 1 ) {$1=samp; print $0 } ' ${SAMP}_${CHR}_fw/product_results.txt  >> ../all_samps_product_results.txt
    # combine results
     cat ${SAMP}_site_results.txt >> ../all_samps_site_results.txt
     cat ${SAMP}_codon_results.txt >> ../all_samps_codon_results.txt
     cat ${SAMP}_sliding_window_length9_results.txt >> ../all_samps_sliding_window_length9_results.txt
     # awk 'NR == 1; NR > 1 {print $0 | "sort -n"}'
   cd ..
 done < $1
 awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_site_results.txt > all_samps_site_results_sorted.txt
 awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_codon_results.txt > all_samps_codon_results_sorted.txt
 awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort -u -k1,2 -k3,3n"}'  all_samps_sliding_window_length9_results.txt > all_samps_sliding_window_length9_results_sorted.txt

 # collect all product results
 #find . -name "product_results.txt" -exec cat {} +  | sort -u | sed 's/file/Sample/g; s/_[LS].*.vcf//g' > all_samps_product_results.txt
 awk 'NR == 1 ;NR > 1 && !/^Sample\t/ {print $0 | "sort"}'  all_samps_product_results.txt > all_samps_product_results_sorted.txt

 exit 0
