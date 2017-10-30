#!/bin/bash
# author Lukas Endler
# created: 10/30/2017, 9:17:18 AM
# Changed: Time-stamp: <10/30/2017, 9:21:31 AM>
# Variables for scripts, looks for different users (vetgrid01 or vetlinux01)
shopt -s extglob

BASEDIR=$(dirname $BASH_SOURCE)
BASEDIR=${BASEDIR%\/LCMV_project\/BSF_0355}

RUN_ID=0355
SCRIPTS=${BASEDIR}/LCMV_project/BSF_0355/
REFGENOME=$BASEDIR/References/viruses_short.fasta
REFGENOME_FULL=$BASEDIR/References/viruses_Mus_musculus.GRCm38.fa.gz
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
elif [ $USER == "vetlinux01" ] ; then

fi


exit
