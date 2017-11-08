#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of bam files and runs breseq over them to predict
#--------------

source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"
LOGFILE=`pwd`/"bsf_"${RUN_ID}"_breseq.log"

while read file; do
  SM=$($SAMTOOLS view -H $file | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g")
  BN="bsf_"${RUN_ID}"_"$SM
  ERRLOG=$BN"_breseq.err.log"
  mkdir $BN
  cd $BN
  echo java -jar $PICARD SamToFastq I=$file F=${SM}_1.fastq.gz F2=${SM}_2.fastq.gz FU=/dev/null 2\>\> $ERRLOG >> $LOGFILE
  java -jar $PICARD SamToFastq I=$file F=${SM}_1.fastq.gz F2=${SM}_2.fastq.gz FU=/dev/null 2>> $ERRLOG
  echo start breseq at `date` >> $LOGFILE
  echo $BRESEQ -n $SM -j 50 -p --polymorphism-frequency-cutoff 0.01 --polymorphism-minimum-coverage-each-strand 2 --require-match-fraction 0.85 --junction-minimum-side-match 10 --junction-indel-split-length 3 -r $REFGENOME -r ${REFGENOME/%\.fasta/\.gff3} ${SM}_1.fastq.gz ${SM}_2.fastq.gz >> $LOGFILE
  $BRESEQ -n $SM -j 50 -p --polymorphism-frequency-cutoff 0.01 --polymorphism-minimum-coverage-each-strand 2 --require-match-fraction 0.85 --junction-minimum-side-match 10 --junction-indel-split-length 3 -r $REFGENOME -r ${REFGENOME/%\.fasta/\.gff3} ${SM}_1.fastq.gz ${SM}_2.fastq.gz 2>> $ERRLOG
  ES=$?
  echo finished breseq at `date` with exit state $ES >> $LOGFILE
  echo $GDTOOLS gd2vcf -r $REFGENOME -r ${REFGENOME/%\.fasta/\.gff3} -o ${BN}.vcf output/output.gd >> $LOGFILE
  $GDTOOLS gd2vcf -r $REFGENOME -r ${REFGENOME/%\.fasta/\.gff3} -o ${BN}.vcf output/output.gd
  awk -v OFS="\t" '/^\#\#[^I]/ {print} /^\#\#INFO/ {print $0; sub("INFO","FORMAT",$0); print $0; } /\#CH/ {print $0,"FORMAT","'$SM'"} !/^\#/ {form=$8;  gsub(/=[^A-Z]+/,":",form); gsub(/;/,":",form); sub(/:$/,"",form); sub(/INDEL:/,"",form);  samp=$8; gsub(/[A-Z4]+=/,"",samp); gsub(/;/,":",samp); sub(/INDEL:/,"",samp); print $0,form,samp}' ${BN}.vcf | $BCFTOOLS norm -f $REFGENOME - 2>> $ERRLOG | bgzip -c > ${BN}_samp_breseq.vcf.gz
  tabix -p vcf ${BN}_samp_breseq.vcf.gz
  cd ..
done < $1

$BCFTOOLS merge -i "DP:sum,AF:max,AD:sum"  -m none -O v "bsf_"*/*_samp_breseq.vcf.gz > all_samp_breseq.vcf

DPS=$2

if [[ -e $DPS ]];
then
  mv all_samp_breseq.vcf all_samp_bed_tmp.vcf
  bgzip all_samp_tmp.vcf
  tabix -f -p vcf all_samp_tmp.vcf.gz
  $BCFTOOLS annotate -a $DPS -c "+FORMAT/DP" all_samp_tmp.vcf.gz > all_samp_breseq.vcf
fi
