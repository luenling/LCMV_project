#!/bin/bash
# author Lukas Endler
# created: 10/30/2017, 9:17:18 AM
# Changed: Time-stamp: <10/30/2017, 9:21:31 AM>
# Variables for scripts, looks for different users (vetgrid01 or vetlinux01)
# set runid with export RUN_ID=0355
# set runid with export RUN_ID=
shopt -s extglob

if [ -z ${RUNBASE+z} ] ; then RUNBASE=`pwd` ; fi
if [ -z ${RUN_ID+z} ] ; then RUN_ID=${RUNBASE//*_/} ; fi

BASEDIR=$(dirname $BASH_SOURCE)
BASEDIR=${BASEDIR%\/LCMV_project\/BSF_general}


SCRIPTS=${BASEDIR}/LCMV_project/BSF_general/
REFGENOME=$BASEDIR/References/viruses_short.fasta
REFGENOME_FULL=$BASEDIR/References/viruses_short_Mus_musculus.GRCm38.dna.toplevel.fa.gz
PRIMERS=$BASEDIR/References/primers.fna

# setup everything for vetgrid01

if [ $USER == "vetgrid01" ]; then
  # general tools
  PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
  SAMTOOLS=/usr/local/bin/samtools
  BWA=/usr/local//Cellar/bwa/0.7.15/bin/bwa
  BCFTOOLS=/usr/local/bin/bcftools
  # Lofreq executable
  LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq
  #varscan2
  VARSCAN=$BASEDIR/Tools/varscan/VarScan.v2.4.0.jar
  #VarDict java base dir
  VARDICT=$BASEDIR/Tools/VarDictJava/
  # BRESEQ path
  BRESEQ=/usr/local/Cellar/breseq/0.31.1/bin/breseq
  GDTOOLS=/usr/local/Cellar/breseq/0.31.1/bin/gdtools

  # setup variables for SNPEFF
  SNPEFF_DATA="/usr/local/Cellar/snpeff/4.2/share/snpeff/data"
  SNPEFF_CONF="/usr/local/Cellar/snpeff/4.2/share/snpeff/snpEff.config"
  SNPEFF="snpEff"
  SNPSIFT="SnpSift"
  # setup GATK path
  GATK=$BASEDIR/Tools/GenomeAnalysisTK_3.5.0.jar
  # setup path to BBDUK
  BBDUK=$BASEDIR/Tools/bbmap/bbduk.sh
  BBMERGE=$BASEDIR/Tools/bbmap/bbmerge.sh
  ADAPTS=$BASEDIR/Tools/bbmap/resources/adapters.fa
  # VPHASER2
  VPHASER2=$BASEDIR/Tools/viral-ngs/intrahost_alt.py
  #BAMTOOLS23=$BASEDIR/Tools/viral-ngs/tools/tools/binaries/V-Phaser-2.0/bamtools-2.3.0/lib
  CONDAPATH=

elif [ $USER == "vetlinux01" ] ; then
  # general tools
  PICARD=/home/vetlinux01/.linuxbrew/Cellar/picard-tools/2.12.1/share/java/picard.jar
  SAMTOOLS=/home/vetlinux01/.linuxbrew/bin/samtools
  BWA=/home/vetlinux01/.linuxbrew/bin/bwa
  BCFTOOLS=/home/vetlinux01/.linuxbrew/bin/bcftools
  # Lofreq executable
  LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq
  #varscan2
  VARSCAN=$BASEDIR/Tools/VarScan.v2.4.2.jar
  #VarDict java base dir
  VARDICT=$BASEDIR/Tools/VarDictJava
  # BRESEQ path
  BRESEQ=/home/vetlinux01/.linuxbrew/Cellar/breseq/0.31.1/bin/breseq
  GDTOOLS=/home/vetlinux01/.linuxbrew/Cellar/breseq/0.31.1/bin/gdtools
  # setup variables for SNPEFF
  SNPEFF_DATA="/home/vetlinux01/.linuxbrew/Cellar/snpeff/4.3p/share/snpeff/data"
  SNPEFF_CONF="/home/vetlinux01/.linuxbrew/Cellar/snpeff/4.3p/share/snpeff/snpEff.config"
  SNPEFF="snpEff"
  SNPSIFT="SnpSift"
  # setup GATK path
  GATK=/home/vetlinux01/Tools/GenomeAnalysisTK.jar
  # setup path to BBDUK
  BBDUK=/home/vetlinux01/Tools/bbmap/bbduk.sh
  BBMERGE=/home/vetlinux01/Tools/bbmap/bbmerge.sh
  BBREPAIR=/home/vetlinux01/Tools/bbmap/repair.sh
  ADAPTS=/home/vetlinux01/Tools/bbmap/resources/adapters.fa
  # VPHASER2
  VPHASER2=$BASEDIR/Tools/viral-ngs/intrahost_alt.py
  BAMTOOLS23=$BASEDIR/Tools/viral-ngs/tools/tools/binaries/V-Phaser-2.0/bamtools-2.3.0/lib
  CONDAPATH=/Volumes/Temp/Lukas/miniconda/bin

fi

# if snpeff database not set up, do that
if [ ! -d  "$SNPEFF_DATA/lcmv" ] ; then
  SNP_LOG="snpeff.log"
  echo creating new lcmv db entry at date >> $SNP_LOG
  echo mkdir -p $SNPEFF_DATA/lcmv >> $SNP_LOG
  mkdir -p $SNPEFF_DATA/lcmv
  echo cp $REFGENOME $SNPEFF_DATA/lcmv/sequences.fa >> $SNP_LOG
  cp $REFGENOME $SNPEFF_DATA/lcmv/sequences.fa
  echo cp ${REFGENOME/%\.fasta/\.gff} $SNPEFF_DATA/lcmv/genes.gff >> $SNP_LOG
  cp ${REFGENOME/%\.fasta/\.gff} $SNPEFF_DATA/lcmv/genes.gff
  echo \echo '-e \"# LCMV\nlcmv.genome : lcmv\n\' \>\> $SNPEFF_CONF >> $SNP_LOG
  echo -e "# LCMV\nlcmv.genome : lcmv\n" >> $SNPEFF_CONF
  echo $SNPEFF build  -c $SNPEFF_CONF -gff3 -dataDir $SNPEFF_DATA >> $SNP_LOG
  $SNPEFF build  -c $SNPEFF_CONF -gff3 -dataDir $SNPEFF_DATA lcmv 2>> $SNP_LOG
  ES=$?
  echo finished creating snpeff db at `date` with exit state $ES >> $SNP_LOG
  [ $ES -eq 0 ] || exit $ES
fi
