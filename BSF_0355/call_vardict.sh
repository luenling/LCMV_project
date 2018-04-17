#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: 11/6/2017, 1:05:34 PM
# date: 11/6/2017, 1:04:06 PM
# takes a list of bam files and a vcf file from samtools and calls variants with VarDict and annotates missing depths from the samtools vcf
#--------------

source $(dirname $BASH_SOURCE)"/bsf_0355_params.sh"
LOGFILE=${RUN_ID}"_vardict.log"

while read file; do
  SM=$($SAMTOOLS view -H $file | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g")
  BN=${RUN_ID}"_"$SM
  ERRLOG=$BN"_vardict.err.log"
  echo start vardict at `date` >> $LOGFILE
  echo $VARDICT/VarDict_java -G $REFGENOME -th 20 -f 0.001 -m 10 -P 10 -B 3 -b $file -c 1 -S 2 -E 3 -z 1 $BASEDIR/References/viruses_short.bed \| $VARDICT/VarDict/teststrandbias.R \| $VARDICT/VarDict/var2vcf_valid.pl -A -N $SM -E -f 0.001 \| sed 's/VCFv4.1/VCFv4.2/; s/\([AV][FD]\),Number=1/\1,Number=A/; s/AD,Number=\./AD,Number=R/;' \| bgzip -c  \> ${BN}_vardict_filt.vcf.gz >> $LOGFILE
  $VARDICT/VarDict_java -G $REFGENOME -th 20 -f 0.001 -m 10 -P 10 -B 3 -b $file -c 1 -S 2 -E 3 -z 1 $BASEDIR/References/viruses_short.bed | $VARDICT/VarDict/teststrandbias.R | $VARDICT/VarDict/var2vcf_valid.pl  -A -N $SM -E -f 0.001 |  sed 's/VCFv4.1/VCFv4.2/; s/\([AV][FD]\),Number=1/\1,Number=A/; s/AD,Number=\./AD,Number=R/;'  | bgzip -c  > ${BN}_vardict_filt.vcf.gz
  ES=$?
  echo finished vardict at `date` with exit state $ES >> $LOGFILE
  echo tabix -p vcf ${BN}_vardict_filt.vcf.gz >> $LOGFILE
  tabix -p vcf ${BN}_vardict_filt.vcf.gz
  echo $BCFTOOLS norm -f $REFGENOME -m+any ${BN}_vardict_filt.vcf.gz \| $BCFTOOLS view --apply-filters "PASS" -O z - \>  ${BN}_vardict_filt_norm.vcf.gz >> $LOGFILE
  $BCFTOOLS norm -f $REFGENOME -m+any ${BN}_vardict_filt.vcf.gz | $BCFTOOLS view --apply-filters "PASS" -O z - >  ${BN}_vardict_filt_norm.vcf.gz
  tabix -p vcf ${BN}_vardict_filt_norm.vcf.gz
done < $1

echo $BCFTOOLS merge -i "SAMPLE:join,DP:sum,AF:max,SBF:min,MQ:avg,QUAL:avg,PMEAN:avg,PSTD:min,BIAS:join,REFBIAS:join,VARBIAS:join"  -m none -O v \*_vardict_filt_norm.vcf.gz \> all_samps_vardict_filt_norm.vcf >> $LOGFILE
$BCFTOOLS merge -i "SAMPLE:join,DP:sum,AF:max,SBF:min,MQ:avg,QUAL:avg,PMEAN:avg,PSTD:min,BIAS:join,REFBIAS:join,VARBIAS:join"  -m none -O v *_vardict_filt_norm.vcf.gz > all_samps_vardict_filt_norm.vcf

DPS=$2

if [[ -e $DPS ]];
then
  mv all_samps_vardict_filt_norm.vcf all_samps_vardict_filt_tmp.vcf
  bgzip -f all_samps_vardict_filt_tmp.vcf
  tabix -f -p vcf all_samps_vardict_filt_tmp.vcf.gz
  echo $BCFTOOLS annotate --collapse all -a $DPS -c \"+FORMAT/DP,FORMAT/AD2:=FORMAT/AD\" all_samps_vardict_filt_tmp.vcf.gz \| sed \''/ID=AD2/ {s/\".*\"/\"Alleleic Depths as from samtools mpileup\"/;s/Number\=R/Number\=\./}' \' \> all_samps_vardict_filt_norm.vcf >> $LOGFILE
  $BCFTOOLS annotate --collapse all -a $DPS -c "+FORMAT/DP,FORMAT/AD2:=FORMAT/AD" all_samps_vardict_filt_tmp.vcf.gz | sed '/ID=AD2/ {s/\".*\"/\"Alleleic Depths as from samtools mpileup\"/;s/Number\=R/Number\=\./}' > all_samps_vardict_filt_norm.vcf
fi

for i in 0.1 0.05 0.01 0.001 ; do
  LFVCF=all_samps_vardict_filt_norm_${i}.vcf
  $BCFTOOLS view -i "AF>$i" all_samps_vardict_filt_norm.vcf | $BCFTOOLS norm -f $REFGENOME -m+any - > $LFVCF
  FN=`basename $LFVCF .vcf`
  LOGFILE=${FN}.log
  ERRORLOG=${FN}.err.log
  echo startingrunning snpeff at `date` >> $LOGFILE
  echo  $SNPEFF -c $SNPEFF_CONF -dataDir $SNPEFF_DATA lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html $LFVCF \> ${FN}_snpeff.vcf >> $LOGFILE
  $SNPEFF  -c $SNPEFF_CONF -dataDir $SNPEFF_DATA lcmv -no-intergenic  -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html $LFVCF > ${FN}_snpeff.vcf

  ES=$?
  echo finished snpeff at `date` with exit state $ES >> $LOGFILE
  [ $ES -eq 0 ] || exit $ES

  $SNPSIFT extractFields ${FN}_snpeff.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" > ${FN}_snpeff.tab
  #$BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%SB\t%DP\t%AF\t%PQ]\n"  $LFVCF > ${FN}.stats.tab
  $BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"  $LFVCF > ${FN}.afs.tab
  $BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n"  $LFVCF > ${FN}.dp.tab
  #$BCFTOOLS query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%PQ]\n"  $LFVCF > ${FN}.pq.tab
  paste ${FN}.afs.tab <(cut -f5- ${FN}.dp.tab) >  ${FN}.afs.dp.tab
  python ${BASEDIR}/Scripts/get_positions_from_sync.py -a ${FN}.afs.dp.tab -b ${FN}_snpeff.tab --both | cat <(head -1  ${FN}.afs.dp.tab )  -  >  ${FN}.afs.anno.tab ;
  sed ' s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.afs.anno.tab >  ${FN}.afs.anno_alt.tab ;

  #sed 's/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.stats.tab >  ${FN}.stats_alt.tab ;

done


exit
