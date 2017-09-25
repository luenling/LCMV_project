#!/bin/bash
# author Lukas Endler
# created: 2017-09-24T13:54:09.196+02:00
# Changed: 9/25/2017, 11:21:52 AM
# takes a bam file of reads with the BSF format BSF_##__HL37TBBXX_6_S_44_0_S28607.bam and trimmes and clips to get 2 fq.gz files
# turn on extended globbing behavior
shopt -s extglob
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_Mus_musculus.GRCm38.fa.gz
PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local//Cellar/bwa/0.7.15/bin/bwa
BBDUK=$BASEDIR/Tools/bbmap/bbduk.sh
BBMERGE=$BASEDIR/Tools/bbmap/bbmerge.sh
ADAPTS=$BASEDIR/Tools/bbmap/resources/adapters.fa
PRIMERS=$BASEDIR/References/primers.fna
#REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
INFILE=$1
FN=`basename $1 .bam`
FN=`basename $FN .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
#INF2=${INF1/_1.f/_2.f}
OUT=$FN
#ls ~/LCMV_project/Run_0277/biomedical-sequencing.at/dna/BSF_0277_HFNMYBBXX_5_samples_41fc5015a64547a9a8193a1d6de705d2/*.bam | xargs -P 5 -L 1 ../../Scripts/bsf0277_bbduk_trim_cutadapts.sh  > run.log &

LOGFILE=$OUT.log
ERRORFILE=$OUT.err.log
#TRIMLOG=$OUT.trim.log
OUTDIR='FASTQ/'
if [ $2 ]; then
    OUTDIR=$2'/'
    if [ ! -d $2 ]; then mkdir $2 ; fi
fi

OUT=${OUTDIR}${OUT}

# echo converting bam to fastq at `date` >> ${LOGFILE}
# echo java -Xmx4g -jar ${PICARD} SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=$INFILE FASTQ=${OUTDIR}${OUT}_1.fastq SECOND_END_FASTQ=${OUTDIR}${OUT}_2.fastq \>\> ${LOGFILE} 2\> ${ERRORFILE} >> ${LOGFILE}
# java -Xmx4g -jar ${PICARD} SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=$INFILE FASTQ=${OUTDIR}${OUT}_1.fastq SECOND_END_FASTQ=${OUTDIR}${OUT}_2.fastq  >> ${LOGFILE} 2>> ${ERRORFILE}
# ES=$?
# echo finished at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES
# echo gzipping fastq ${OUTDIR}${OUT}_1.fastq at `date` >> ${LOGFILE}
# gzip -fq ${OUTDIR}${OUT}_1.fastq  >> ${LOGFILE} 2>> ${ERRORFILE}
# ES=$?
# echo finished at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES
# echo gzipping fastq ${OUTDIR}${OUT}_2.fastq at `date` >> ${LOGFILE}
# gzip -fq ${OUTDIR}${OUT}_2.fastq >> ${LOGFILE} 2>> ${ERRORFILE}
# ES=$?
# echo finished at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES

# INF1=${OUTDIR}${OUT}_1.fastq.gz
# INF2=${OUTDIR}${OUT}_2.fastq.gz



echo trimming and removing adaptors at `date` >> ${LOGFILE}
echo $SAMTOOLS fastq $1 | $BBDUK -Xmx4g in=stdin.fq out1=${OUT}_npa_1.fq.gz out2=${OUT}_npa_2.fq.gz qtrim=r minlen=75 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t \>\> ${LOGFILE} 2\> ${ERRORFILE} >> ${LOGFILE}
$SAMTOOLS fastq $1 | $BBDUK -Xmx4g in=stdin.fq out1=${OUT}_npa_1.fq.gz out2=${OUT}_npa_2.fq.gz qtrim=r minlen=55 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t >> ${LOGFILE} 2> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
if [ ! -d ${OUTDIR}/fastqc ]; then mkdir ${OUTDIR}/fastqc ; fi
fastqc -q -o ${OUTDIR}/fastqc ${OUT}_npa_1.fq.gz ${OUT}_npa_2.fq.gz &


#
#
# if [ ! -d Merged ]; then mkdir Merged ; fi
# echo merging overlapping reads at `date` >> ${LOGFILE}
# echo $BBMERGE -Xmx4g  in1=${OUT}_npa_1.fq.gz in2=${OUT}_npa_2.fq.gz out=Merged/${OUT}_m.fq.gz outu1=Merged/${OUT}_m_1.fq.gz outu2=Merged/${OUT}_m_2.fq.gz strict=t qtrim2=r trimq=10 ouq=t ihist=Merged/${OUT}_m.hist 2\> Merged/${OUT}_m.log >>  $LOGFILE
# $BBMERGE -Xmx4g  in1=${OUT}_npa_1.fq.gz in2=${OUT}_npa_2.fq.gz out=Merged/${OUT}_m.fq.gz outu1=Merged/${OUT}_m_1.fq.gz outu2=Merged/${OUT}_m_2.fq.gz strict=t qtrim2=r trimq=10 ouq=t ihist=Merged/${OUT}_m.hist 2> Merged/${OUT}_m.log
# ES=$?
# echo finished at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES
