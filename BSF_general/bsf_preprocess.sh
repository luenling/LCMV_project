#!/bin/bash
# author Lukas Endler
# created: 2017-09-24T13:54:09.196+02:00
# Changed: 9/25/2017, 11:21:52 AM
# takes a bam file of reads with the BSF format name and the name of an outputdir (default: FASTQ) BSF_##__HL37TBBXX_6_S_44_0_S28607.bam and trimmes and clips to get 2 fq.gz files
# turn on extended globbing behavior
shopt -s extglob

source $(dirname $BASH_SOURCE)"/bsf_params.sh"

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

OUTDIR='FASTQ/'
if [ $2 ]; then
    OUTDIR=$2'/'
    if [ ! -d $2 ]; then mkdir $2 ; fi
fi
LOGFILE=${OUTDIR}/$OUT.log
ERRORFILE=${OUTDIR}/$OUT.err.log

OUT=${OUTDIR}${OUT}
STATDIR=${OUTDIR}/Stats
if [ ! -d $STATDIR ]; then mkdir $STATDIR ; fi

echo trimming and removing adaptors at `date` >> ${LOGFILE}
echo $SAMTOOLS fastq $1 | $BBDUK -Xmx4g in=stdin.fq out1=${OUT}_npa_1.fq.gz out2=${OUT}_npa_2.fq.gz qtrim=r minlen=75 ktrim=r ref=\"$ADAPTS\" int=t k=21 mink=11 trimq=15 overwrite=t tbo=t tpe=t lhist=$STATDIR/${FN}.bbuk.lhist stats=$STATDIR/${FN}.bbuk.stats \>\> ${LOGFILE} 2\> ${ERRORFILE} >> ${LOGFILE}
$SAMTOOLS fastq $1 | $BBDUK -Xmx4g in=stdin.fq out1=${OUT}_npa_1.fq.gz out2=${OUT}_npa_2.fq.gz qtrim=r minlen=55 ktrim=r ref=\"$ADAPTS\" int=t k=21 mink=11 trimq=15 overwrite=t tbo=t tpe=t lhist=$STATDIR/${FN}.bbuk.lhist stats=$STATDIR/${FN}.bbuk.stats >> ${LOGFILE} 2> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
if [ ! -d ${OUTDIR}/fastqc ]; then mkdir ${OUTDIR}/fastqc ; fi
fastqc -q -o ${OUTDIR}/fastqc ${OUT}_npa_1.fq.gz ${OUT}_npa_2.fq.gz &
