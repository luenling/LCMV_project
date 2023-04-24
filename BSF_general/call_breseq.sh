#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of bam files and runs breseq over them to predict
#--------------

source $(dirname $BASH_SOURCE)"/bsf_params.sh"
LOGFILE=`pwd`/"bsf_"${RUN_ID}"_breseq.log"

while read file; do
  SM=$($SAMTOOLS view -H $file | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g")
  BN="bsf_"${RUN_ID}"_"$SM
  ERRLOG=$BN"_breseq.err.log"
  mkdir -p $BN
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
  awk -v OFS="\t" '/^\#\#[^I]/ {print} /^\#\#INFO/ {sub("AD,Number=1","AD,Number=A",$0); print $0; sub("INFO","FORMAT",$0); print $0; } /\#CH/ {print "##FORMAT=<ID=PQ,Number=1,Type=Float,Description=\"Breseq variant quality score (log10(pvariant/Pn)-log10(total length or references))\">"; print $0,"FORMAT","'$SM'"} !/^\#/ {form=$8;  gsub(/=[^A-Z]+/,":",form); gsub(/;/,":",form); sub(/:$/,"",form); sub(/INDEL:/,"",form);  samp=$8; gsub(/[A-Z4]+=/,"",samp); gsub(/;/,":",samp); sub(/INDEL:/,"",samp); print $0,form":PQ",samp":"$6}' ${BN}.vcf | sed '/ID\=[AD][DP]/ s/Float/Integer/g' | $BCFTOOLS norm -f $REFGENOME - 2>> $ERRLOG | bgzip -c > ${BN}_samp_breseq.vcf.gz
  rm -f ${BN}_samp_breseq.vcf.gz.tbi
  tabix -p vcf ${BN}_samp_breseq.vcf.gz
  cd ..
done < $1


gdtools COMPARE -r $REFGENOME -r ${REFGENOME/%\.fasta/\.gff3} -o all_samp_breseq_comp.html --repeat-header 25 `ls bsf_${RUN_ID}_S*/output/output.gd`

gdtools COMPARE -r $REFGENOME -r ${REFGENOME/%\.fasta/\.gff3} -o all_samp_breseq_comp.tsv -f TSV `ls bsf_${RUN_ID}_S*/output/output.gd`


$BCFTOOLS merge -i "DP:sum,AF:max,AD:sum"  -m none -O v "bsf_"*/*_samp_breseq.vcf.gz > all_samp_breseq.vcf

DPS=$2

if [[ -e $DPS ]];
then
  mv all_samp_breseq.vcf all_samp_tmp.vcf
  bgzip all_samp_tmp.vcf
  tabix -f -p vcf all_samp_tmp.vcf.gz
  $BCFTOOLS annotate -a $DPS -c "+FORMAT/DP" all_samp_tmp.vcf.gz > all_samp_breseq.vcf
fi

for i in 0.1 0.05 0.001 ; do
  LFVCF=all_samp_breseq_${i}.vcf
  $BCFTOOLS view -i "AF>$i" all_samp_breseq.vcf | $BCFTOOLS norm -f $REFGENOME -m+any - > $LFVCF
  FN=`basename $LFVCF .vcf`
  LOGFILE=${FN}.log
  ERRORLOG=${FN}.err.log
  # run snpeff
  echo starting snpeff at `date` >> $LOGFILE
  echo  $SNPEFF -c $SNPEFF_CONF -dataDir $SNPEFF_DATA lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html $LFVCF \> ${FN}_snpeff.vcf >> $LOGFILE
  $SNPEFF  -c $SNPEFF_CONF -dataDir $SNPEFF_DATA lcmv -no-intergenic  -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html $LFVCF > ${FN}_snpeff.vcf
  ES=$?
  echo finished snpeff at `date` with exit state $ES >> $LOGFILE
  [ $ES -eq 0 ] || exit $ES

  # extract annotations and combine them with AFs and statistics
  $SNPSIFT extractFields ${FN}_snpeff.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" | sed ' s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' > ${FN}_snpeff.tab
  $BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF\t%DP\t%AD\t%PQ]\n"  $LFVCF > ${FN}.stats.tab
  $BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"  $LFVCF > ${FN}.afs.tab
  $BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n"  $LFVCF > ${FN}.dp.tab
  $BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%PQ]\n"  $LFVCF > ${FN}.pq.tab
  paste ${FN}.afs.tab <(cut -f5- ${FN}.dp.tab) <(cut -f5- ${FN}.pq.tab) >  ${FN}.afs.dp.pq.tab
  python ${BASEDIR}/Scripts/get_positions_from_sync.py -a ${FN}.afs.dp.pq.tab -b ${FN}_snpeff.tab --both | cat <(head -1  ${FN}.afs.dp.pq.tab )  -  >  ${FN}.afs.anno.tab ;
  #
  # sed ' s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.afs.anno.tab >  ${FN}.afs.anno_alt.tab ;
  # sed 's/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.stats.tab >  ${FN}.stats_alt.tab ;
done
